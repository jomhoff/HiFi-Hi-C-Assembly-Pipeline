#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * PacBio HiFi + Hi-C chromosome-scale assembly:
 * hifiasm -> draft QC -> BWA/samtools -> YaHS -> final QC -> Juicebox .hic
 */

params.reads = null
params.hic_r1 = null
params.hic_r2 = null
params.sample = 'hifi_assembly'
params.outdir = 'results'
params.hifiasm_prefix = "${params.sample}.asm"
params.hifiasm_extra = '-l 2'
params.yahs_extra = ''
params.yahs_mapq = 10
params.lineage = 'sauropsida'
params.odb = 'odb12'
params.compleasm_library = null
params.compleasm_extra = ''
params.download_lineage = false
params.run_compleasm_on = 'total,scaffolded'
params.juicer_tools_jar = null
params.juicer_tools_url = 'https://s3.amazonaws.com/hicfiles.tc4ga.com/public/juicer/juicer_tools_1.22.01.jar'
params.juicer_java_memory = '48g'
params.make_hic = true

workflow {
    if (!params.reads) {
        error "Missing required parameter: --reads '/path/to/*.fastq.gz'"
    }

    useHiC = params.hic_r1 && params.hic_r2
    if ((params.hic_r1 && !params.hic_r2) || (!params.hic_r1 && params.hic_r2)) {
        error "Hi-C mode needs both --hic_r1 and --hic_r2"
    }

    runCompleasmOn = params.run_compleasm_on
        .split(',')
        .collect { it.trim() }
        .findAll { it }
    if (!runCompleasmOn.every { it in ['hap1', 'hap2', 'total', 'scaffolded', 'all'] }) {
        error "--run_compleasm_on must contain hap1, hap2, total, scaffolded, or all"
    }
    if ('all' in runCompleasmOn) {
        runCompleasmOn = useHiC
            ? ['hap1', 'hap2', 'total', 'scaffolded']
            : ['hap1', 'hap2', 'total']
    }
    if (!useHiC && 'scaffolded' in runCompleasmOn) {
        runCompleasmOn = runCompleasmOn.findAll { it != 'scaffolded' }
    }

    downloadLineage = params.download_lineage.toString().toBoolean()
    makeHic = params.make_hic.toString().toBoolean()

    reads_ch = Channel.fromPath(params.reads, checkIfExists: true).collect()
    hic_r1_ch = useHiC
        ? Channel.fromPath(params.hic_r1, checkIfExists: true).collect()
        : Channel.value([])
    hic_r2_ch = useHiC
        ? Channel.fromPath(params.hic_r2, checkIfExists: true).collect()
        : Channel.value([])

    HIFIASM(reads_ch, hic_r1_ch, hic_r2_ch)
    draft_gfas_ch = HIFIASM.out.hap1_gfa
        .mix(HIFIASM.out.hap2_gfa)
        .mix(HIFIASM.out.total_gfa)
    CHECK_HIFIASM_OUTPUTS(draft_gfas_ch)
    GFA_TO_FASTA(draft_gfas_ch)
    CHECK_FASTA(GFA_TO_FASTA.out.fasta)

    all_fastas_ch = GFA_TO_FASTA.out.fasta
    hic_bam_checks_ch = Channel.value([])
    yahs_checks_ch = Channel.value([])
    contact_map_checks_ch = Channel.value([])

    if (useHiC) {
        total_fasta_ch = GFA_TO_FASTA.out.fasta.filter { asm_type, fasta ->
            asm_type == 'total'
        }

        HIC_ALIGN(total_fasta_ch, hic_r1_ch, hic_r2_ch)
        PREPARE_HIC_BAM(HIC_ALIGN.out.raw_bam)
        CHECK_HIC_BAM(PREPARE_HIC_BAM.out.final_bam)
        hic_bam_checks_ch = CHECK_HIC_BAM.out.report.collect()

        YAHS(total_fasta_ch, PREPARE_HIC_BAM.out.final_bam)
        CHECK_YAHS(YAHS.out.scaffold)
        yahs_checks_ch = CHECK_YAHS.out.report.collect()

        scaffold_fasta_ch = YAHS.out.scaffold.map { asm_type, fasta, agp, contig_fai ->
            tuple(asm_type, fasta)
        }
        all_fastas_ch = GFA_TO_FASTA.out.fasta.mix(scaffold_fasta_ch)

        if (makeHic) {
            JUICER_PRE(YAHS.out.scaffold, PREPARE_HIC_BAM.out.final_bam)

            if (params.juicer_tools_jar) {
                juicer_jar_ch = Channel.fromPath(params.juicer_tools_jar, checkIfExists: true)
            } else {
                DOWNLOAD_JUICER_TOOLS()
                juicer_jar_ch = DOWNLOAD_JUICER_TOOLS.out.jar
            }

            BUILD_HIC_MAP(JUICER_PRE.out.pre, juicer_jar_ch)
            CHECK_CONTACT_MAP(BUILD_HIC_MAP.out.hic)
            contact_map_checks_ch = CHECK_CONTACT_MAP.out.report.collect()
        }
    }

    ASSEMBLY_STATS(all_fastas_ch)
    CHECK_ASSEMBLY_STATS(ASSEMBLY_STATS.out.stats)

    compleasm_input_ch = all_fastas_ch.filter { asm_type, fasta ->
        runCompleasmOn.contains(asm_type)
    }
    if (downloadLineage) {
        DOWNLOAD_COMPLEASM_LINEAGE()
        COMPLEASM(compleasm_input_ch, DOWNLOAD_COMPLEASM_LINEAGE.out.done)
    } else {
        COMPLEASM(compleasm_input_ch, Channel.value('lineage_download_not_requested'))
    }
    CHECK_COMPLEASM(COMPLEASM.out.summary)

    REPORT(
        CHECK_HIFIASM_OUTPUTS.out.report.collect(),
        CHECK_FASTA.out.report.collect(),
        ASSEMBLY_STATS.out.stats.map { asm_type, stats_json -> stats_json }.collect(),
        CHECK_COMPLEASM.out.report.collect(),
        hic_bam_checks_ch,
        yahs_checks_ch,
        contact_map_checks_ch
    )
}

process HIFIASM {
    tag "${params.sample}"
    label 'hifiasm'
    publishDir "${params.outdir}/01_hifiasm", mode: 'copy'

    input:
    path reads
    path hic_r1
    path hic_r2

    output:
    tuple val('hap1'), path("${params.sample}.hap1.p_ctg.gfa"), emit: hap1_gfa
    tuple val('hap2'), path("${params.sample}.hap2.p_ctg.gfa"), emit: hap2_gfa
    tuple val('total'), path("${params.sample}.total.p_ctg.gfa"), emit: total_gfa

    script:
    def hicPrep = hic_r1
        ? """
          cat ${hic_r1.collect { "\"${it}\"" }.join(' ')} > hifiasm_hic_R1.fastq.gz
          cat ${hic_r2.collect { "\"${it}\"" }.join(' ')} > hifiasm_hic_R2.fastq.gz
          """
        : ''
    def hicArgs = hic_r1
        ? '--h1 hifiasm_hic_R1.fastq.gz --h2 hifiasm_hic_R2.fastq.gz'
        : ''
    def sourceStem = hic_r1 ? 'hic' : 'bp'
    """
    ${hicPrep}
    hifiasm -o ${params.hifiasm_prefix} ${params.hifiasm_extra} -t ${task.cpus} \
        ${hicArgs} ${reads.collect { "\"${it}\"" }.join(' ')}

    for gfa in \
        "${params.hifiasm_prefix}.${sourceStem}.hap1.p_ctg.gfa" \
        "${params.hifiasm_prefix}.${sourceStem}.hap2.p_ctg.gfa" \
        "${params.hifiasm_prefix}.${sourceStem}.p_ctg.gfa"
    do
        test -s "\$gfa"
    done

    cp "${params.hifiasm_prefix}.${sourceStem}.hap1.p_ctg.gfa" "${params.sample}.hap1.p_ctg.gfa"
    cp "${params.hifiasm_prefix}.${sourceStem}.hap2.p_ctg.gfa" "${params.sample}.hap2.p_ctg.gfa"
    cp "${params.hifiasm_prefix}.${sourceStem}.p_ctg.gfa" "${params.sample}.total.p_ctg.gfa"
    """
}

process CHECK_HIFIASM_OUTPUTS {
    tag "$asm_type"
    label 'small'
    publishDir "${params.outdir}/checks/01_hifiasm", mode: 'copy'

    input:
    tuple val(asm_type), path(gfa)

    output:
    path "${params.sample}.${asm_type}.hifiasm_check.tsv", emit: report

    script:
    """
    awk -v sample="${params.sample}" -v asm_type="${asm_type}" '
        BEGIN { OFS="\\t"; segments=0; links=0; bases=0 }
        /^S/ { segments++; bases += length(\$3) }
        /^L/ { links++ }
        END {
            if (segments == 0 || bases == 0) exit 1
            print "sample","assembly","gfa","segments","links","segment_bases","status"
            print sample,asm_type,FILENAME,segments,links,bases,"PASS"
        }
    ' "${gfa}" > "${params.sample}.${asm_type}.hifiasm_check.tsv"
    """
}

process GFA_TO_FASTA {
    tag "$asm_type"
    label 'small'
    publishDir "${params.outdir}/02_draft_fastas", mode: 'copy'

    input:
    tuple val(asm_type), path(gfa)

    output:
    tuple val(asm_type), path("${params.sample}.${asm_type}.ctg.fa"), emit: fasta

    script:
    """
    awk '/^S/{print ">"\$2; print \$3}' "${gfa}" > "${params.sample}.${asm_type}.ctg.fa"
    test -s "${params.sample}.${asm_type}.ctg.fa"
    """
}

process CHECK_FASTA {
    tag "$asm_type"
    label 'small'
    publishDir "${params.outdir}/checks/02_draft_fastas", mode: 'copy'

    input:
    tuple val(asm_type), path(fasta)

    output:
    path "${params.sample}.${asm_type}.fasta_check.tsv", emit: report

    script:
    """
    awk -v sample="${params.sample}" -v asm_type="${asm_type}" '
        BEGIN { OFS="\\t"; seqs=0; bases=0 }
        /^>/ { seqs++; next }
        { bases += length(\$0) }
        END {
            if (seqs == 0 || bases == 0) exit 1
            print "sample","assembly","fasta","sequences","bases","status"
            print sample,asm_type,FILENAME,seqs,bases,"PASS"
        }
    ' "${fasta}" > "${params.sample}.${asm_type}.fasta_check.tsv"
    """
}

process HIC_ALIGN {
    tag "${params.sample}"
    label 'hic_mapping'
    publishDir "${params.outdir}/03_hic_alignment", mode: 'copy'

    input:
    tuple val(asm_type), path(fasta)
    path hic_r1
    path hic_r2

    output:
    tuple val(asm_type), path("${params.sample}.hic.raw.bam"), emit: raw_bam

    script:
    """
    cat ${hic_r1.collect { "\"${it}\"" }.join(' ')} > hic_R1.fastq.gz
    cat ${hic_r2.collect { "\"${it}\"" }.join(' ')} > hic_R2.fastq.gz

    bwa index "${fasta}"
    bwa mem -t ${task.cpus} -5SP "${fasta}" hic_R1.fastq.gz hic_R2.fastq.gz \
        | samtools view -@ ${task.cpus} -b -o "${params.sample}.hic.raw.bam" -

    samtools quickcheck -v "${params.sample}.hic.raw.bam"
    """
}

process PREPARE_HIC_BAM {
    tag "${params.sample}"
    label 'hic_mapping'
    publishDir "${params.outdir}/03_hic_alignment", mode: 'copy'

    input:
    tuple val(asm_type), path(raw_bam)

    output:
    tuple val(asm_type), path("${params.sample}.hic.final.qname.bam"), emit: final_bam
    path "${params.sample}.hic.markdup.metrics.txt", emit: metrics

    script:
    """
    samtools sort -n -@ ${task.cpus} -o hic.qnamesorted.bam "${raw_bam}"
    samtools fixmate -m -@ ${task.cpus} hic.qnamesorted.bam hic.fixmate.bam
    samtools sort -@ ${task.cpus} -o hic.fixmate.sorted.bam hic.fixmate.bam
    samtools index -@ ${task.cpus} hic.fixmate.sorted.bam
    samtools markdup -@ ${task.cpus} -r -s \
        -f "${params.sample}.hic.markdup.metrics.txt" \
        hic.fixmate.sorted.bam hic.dedup.bam
    samtools sort -n -@ ${task.cpus} \
        -o "${params.sample}.hic.final.qname.bam" hic.dedup.bam

    samtools quickcheck -v "${params.sample}.hic.final.qname.bam"
    """
}

process CHECK_HIC_BAM {
    tag "${params.sample}"
    label 'small'
    publishDir "${params.outdir}/checks/03_hic_alignment", mode: 'copy'

    input:
    tuple val(asm_type), path(bam)

    output:
    path "${params.sample}.hic_bam_check.tsv", emit: report

    script:
    """
    samtools quickcheck -v "${bam}"
    total=\$(samtools view -c "${bam}")
    mapped=\$(samtools view -c -F 4 "${bam}")
    test "\$total" -gt 0
    test "\$mapped" -gt 0
    printf 'sample\\tassembly\\tbam\\ttotal_records\\tmapped_records\\tstatus\\n' \
        > "${params.sample}.hic_bam_check.tsv"
    printf '${params.sample}\\t${asm_type}\\t%s\\t%s\\t%s\\tPASS\\n' \
        "${bam}" "\$total" "\$mapped" >> "${params.sample}.hic_bam_check.tsv"
    """
}

process YAHS {
    tag "${params.sample}"
    label 'yahs'
    publishDir "${params.outdir}/04_yahs", mode: 'copy'

    input:
    tuple val(asm_type), path(fasta)
    tuple val(bam_asm_type), path(bam)

    output:
    tuple val('scaffolded'),
        path("${params.sample}.scaffolds_final.fa"),
        path("${params.sample}.scaffolds_final.agp"),
        path("${params.sample}.total.ctg.fa.fai"),
        emit: scaffold
    path "${params.sample}.yahs.bin", optional: true, emit: bin

    script:
    """
    cp "${fasta}" "${params.sample}.total.ctg.fa"
    samtools faidx "${params.sample}.total.ctg.fa"

    yahs "${params.sample}.total.ctg.fa" "${bam}" \
        -o "${params.sample}.yahs" -q ${params.yahs_mapq} ${params.yahs_extra}

    cp "${params.sample}.yahs_scaffolds_final.fa" "${params.sample}.scaffolds_final.fa"
    cp "${params.sample}.yahs_scaffolds_final.agp" "${params.sample}.scaffolds_final.agp"
    test -s "${params.sample}.scaffolds_final.fa"
    test -s "${params.sample}.scaffolds_final.agp"
    """
}

process CHECK_YAHS {
    tag "${params.sample}"
    label 'small'
    publishDir "${params.outdir}/checks/04_yahs", mode: 'copy'

    input:
    tuple val(asm_type), path(fasta), path(agp), path(contig_fai)

    output:
    path "${params.sample}.yahs_check.tsv", emit: report

    script:
    """
    samtools faidx "${fasta}"
    scaffolds=\$(cut -f1 "${fasta}.fai" | wc -l | awk '{print \$1}')
    bases=\$(awk '{sum += \$2} END {print sum+0}' "${fasta}.fai")
    components=\$(awk '\$5 == "W" {n++} END {print n+0}' "${agp}")
    gaps=\$(awk '\$5 == "N" || \$5 == "U" {n++} END {print n+0}' "${agp}")
    test "\$scaffolds" -gt 0
    test "\$bases" -gt 0
    test "\$components" -gt 0
    printf 'sample\\tassembly\\tfasta\\tagp\\tscaffolds\\tbases\\tcomponents\\tgaps\\tstatus\\n' \
        > "${params.sample}.yahs_check.tsv"
    printf '${params.sample}\\t${asm_type}\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\tPASS\\n' \
        "${fasta}" "${agp}" "\$scaffolds" "\$bases" "\$components" "\$gaps" \
        >> "${params.sample}.yahs_check.tsv"
    """
}

process ASSEMBLY_STATS {
    tag "$asm_type"
    label 'stats'
    publishDir "${params.outdir}/05_assembly_stats", mode: 'copy'

    input:
    tuple val(asm_type), path(fasta)

    output:
    tuple val(asm_type), path("${params.sample}.${asm_type}.assembly_stats.json"), emit: stats

    script:
    """
    assemblystats.py "${fasta}" > "${params.sample}.${asm_type}.assembly_stats.json"
    test -s "${params.sample}.${asm_type}.assembly_stats.json"
    """
}

process CHECK_ASSEMBLY_STATS {
    tag "$asm_type"
    label 'small'
    publishDir "${params.outdir}/checks/05_assembly_stats", mode: 'copy'

    input:
    tuple val(asm_type), path(stats_json)

    output:
    path "${params.sample}.${asm_type}.assembly_stats_check.tsv", emit: report

    script:
    """
    python3 ${projectDir}/bin/check_assembly_stats.py \
        --sample "${params.sample}" \
        --assembly "${asm_type}" \
        --stats "${stats_json}" \
        > "${params.sample}.${asm_type}.assembly_stats_check.tsv"
    """
}

process DOWNLOAD_COMPLEASM_LINEAGE {
    tag "${params.lineage}_${params.odb}"
    label 'compleasm'
    publishDir "${params.outdir}/00_compleasm_lineages", mode: 'copy'

    output:
    path "download_${params.lineage}_${params.odb}.done", emit: done

    script:
    def libraryPath = params.compleasm_library ?:
        "${launchDir}/${params.outdir}/00_compleasm_lineages/library"
    def libraryArg = "-L \"${libraryPath}\""
    """
    compleasm download ${params.lineage} --odb ${params.odb} ${libraryArg}
    touch "download_${params.lineage}_${params.odb}.done"
    """
}

process COMPLEASM {
    tag "$asm_type"
    label 'compleasm'
    publishDir "${params.outdir}/06_compleasm", mode: 'copy'

    input:
    tuple val(asm_type), path(fasta)
    val lineage_ready

    output:
    tuple val(asm_type), path("${params.sample}.${asm_type}.compleasm"), emit: outdir
    tuple val(asm_type), path("${params.sample}.${asm_type}.compleasm/summary.txt"), emit: summary

    script:
    def libraryArg = params.compleasm_library
        ? "-L \"${params.compleasm_library}\""
        : (params.download_lineage.toString().toBoolean()
            ? "-L \"${launchDir}/${params.outdir}/00_compleasm_lineages/library\""
            : '')
    """
    compleasm run -a "${fasta}" -o "${params.sample}.${asm_type}.compleasm" \
        -l "${params.lineage}" --odb "${params.odb}" -t ${task.cpus} \
        ${libraryArg} ${params.compleasm_extra}
    test -s "${params.sample}.${asm_type}.compleasm/summary.txt"
    """
}

process CHECK_COMPLEASM {
    tag "$asm_type"
    label 'small'
    publishDir "${params.outdir}/checks/06_compleasm", mode: 'copy'

    input:
    tuple val(asm_type), path(summary)

    output:
    path "${params.sample}.${asm_type}.compleasm_check.tsv", emit: report

    script:
    """
    python3 ${projectDir}/bin/check_compleasm_summary.py \
        --sample "${params.sample}" \
        --assembly "${asm_type}" \
        --summary "${summary}" \
        > "${params.sample}.${asm_type}.compleasm_check.tsv"
    """
}

process JUICER_PRE {
    tag "${params.sample}"
    label 'contact_map'
    publishDir "${params.outdir}/07_contact_map", mode: 'copy'

    input:
    tuple val(asm_type), path(fasta), path(agp), path(contig_fai)
    tuple val(bam_asm_type), path(bam)

    output:
    tuple val(asm_type),
        path("${params.sample}.alignments_sorted.txt"),
        path("${params.sample}.scaffolds_final.chrom.sizes"),
        emit: pre

    script:
    """
    juicer pre "${bam}" "${agp}" "${contig_fai}" \
        | sort -k2,2d -k6,6d -T . --parallel=${task.cpus} -S${task.memory.toMega()}M \
        | awk 'NF' > "${params.sample}.alignments_sorted.txt.part"
    mv "${params.sample}.alignments_sorted.txt.part" "${params.sample}.alignments_sorted.txt"

    samtools faidx "${fasta}"
    cut -f1,2 "${fasta}.fai" > "${params.sample}.scaffolds_final.chrom.sizes"
    test -s "${params.sample}.alignments_sorted.txt"
    test -s "${params.sample}.scaffolds_final.chrom.sizes"
    """
}

process DOWNLOAD_JUICER_TOOLS {
    tag 'juicer_tools_1.22.01'
    label 'download'
    publishDir "${params.outdir}/00_tools", mode: 'copy'

    output:
    path "juicer_tools_1.22.01.jar", emit: jar

    script:
    """
    curl -L --fail --retry 3 \
        -o juicer_tools_1.22.01.jar "${params.juicer_tools_url}"
    test -s juicer_tools_1.22.01.jar
    """
}

process BUILD_HIC_MAP {
    tag "${params.sample}"
    label 'contact_map'
    publishDir "${params.outdir}/07_contact_map", mode: 'copy'

    input:
    tuple val(asm_type), path(alignments), path(chrom_sizes)
    path juicer_jar

    output:
    tuple val(asm_type), path("${params.sample}.scaffolded.hic"), emit: hic

    script:
    """
    java -Xmx${params.juicer_java_memory} -jar "${juicer_jar}" pre \
        "${alignments}" "${params.sample}.scaffolded.hic.part" "${chrom_sizes}"
    mv "${params.sample}.scaffolded.hic.part" "${params.sample}.scaffolded.hic"
    test -s "${params.sample}.scaffolded.hic"
    """
}

process CHECK_CONTACT_MAP {
    tag "$asm_type"
    label 'small'
    publishDir "${params.outdir}/checks/07_contact_map", mode: 'copy'

    input:
    tuple val(asm_type), path(hic)

    output:
    path "${params.sample}.${asm_type}.contact_map_check.tsv", emit: report

    script:
    """
    python3 ${projectDir}/bin/check_contact_map.py \
        --sample "${params.sample}" \
        --assembly "${asm_type}" \
        --map "${hic}" \
        > "${params.sample}.${asm_type}.contact_map_check.tsv"
    """
}

process REPORT {
    label 'small'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path hifiasm_checks
    path fasta_checks
    path stats_json
    path compleasm_checks
    path hic_bam_checks
    path yahs_checks
    path contact_map_checks

    output:
    path "${params.sample}.pipeline_report.md"

    script:
    """
    python3 ${projectDir}/bin/build_report.py \
        --sample "${params.sample}" \
        --hifiasm-checks ${hifiasm_checks.join(' ')} \
        --fasta-checks ${fasta_checks.join(' ')} \
        --assembly-stats ${stats_json.join(' ')} \
        --compleasm-checks ${compleasm_checks.join(' ')} \
        --hic-bam-checks ${hic_bam_checks.join(' ')} \
        --yahs-checks ${yahs_checks.join(' ')} \
        --contact-map-checks ${contact_map_checks.join(' ')} \
        > "${params.sample}.pipeline_report.md"
    """
}
