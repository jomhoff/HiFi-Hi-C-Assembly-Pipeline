# HiFi + Hi-C Chromosome Assembly Pipeline

An automated Nextflow workflow that turns PacBio HiFi and Hi-C reads into a
**final Hi-C-scaffolded, chromosome-scale genome assembly**.

## Primary Result: Final Scaffolded Genome

The main deliverable is the YaHS-scaffolded genome, not only the initial
hifiasm contigs:

```text
04_yahs/<sample>.scaffolds_final.fa
04_yahs/<sample>.scaffolds_final.agp
07_contact_map/<sample>.scaffolded.hic
```

- `scaffolds_final.fa` is the final genome assembly after Hi-C scaffolding.
- `scaffolds_final.agp` records contig placement, order, orientation, and gaps.
- `scaffolded.hic` is the matching contact map for inspection in Juicebox.

The pipeline validates the scaffolded FASTA and AGP, recalculates assembly
statistics, and runs Compleasm on the final genome so the report directly
compares draft-contig and post-scaffolding quality.

This pipeline automates the complete workflow in the original *Plestiodon
fasciatus* assembly notes:

1. assemble PacBio HiFi reads with `hifiasm -l 2`
2. create haplotype 1, haplotype 2, and total primary-contig FASTAs
3. validate every GFA and FASTA
4. calculate `assemblystats.py` metrics for all draft assemblies
5. run Compleasm instead of BUSCO
6. align Hi-C reads to the total assembly with `bwa mem -5SP`
7. name-sort, fix mates, coordinate-sort, index, remove duplicates, and
   name-sort the Hi-C BAM for YaHS
8. scaffold the assembly with YaHS
9. validate the final FASTA and AGP and calculate post-scaffolding statistics
10. run Compleasm on the final scaffolded genome
11. run YaHS `juicer pre` and Juicer Tools to create a Juicebox-compatible
    `.hic` contact map
12. write one Markdown QC report

Nextflow runs independent work concurrently. Draft FASTA validation,
assembly-statistics jobs, and selected Compleasm jobs fan out by assembly.
Hi-C processing starts as soon as the total FASTA is available. Final
statistics and final Compleasm start as soon as YaHS finishes.

## Full Run

```bash
nextflow run . \
  -profile conda \
  --reads '/path/to/hifi_reads/*.fastq.gz' \
  --hic_r1 '/path/to/hic/*_R1*.fastq.gz' \
  --hic_r2 '/path/to/hic/*_R2*.fastq.gz' \
  --sample pfas \
  --outdir results_pfas \
  --lineage sauropsida \
  --odb odb12 \
  --run_compleasm_on total,scaffolded
```

With Hi-C reads supplied, this is a chromosome-scale assembly workflow:
hifiasm performs Hi-C-aware phasing, YaHS scaffolds the total primary-contig
assembly, and the final `.hic` map is built in scaffold coordinates.

## Offline Cluster Run

Compleasm databases and Juicer Tools should be prepared on a login node when
compute nodes cannot access the internet.

Download the Compleasm lineage:

```bash
compleasm download sauropsida --odb odb12 -L /path/to/compleasm_lineages
```

Download `juicer_tools_1.22.01.jar` from the official Juicer download page,
then provide both local paths:

```bash
nextflow run . \
  -profile conda,slurm \
  --reads '/path/to/hifi_reads/*.fastq.gz' \
  --hic_r1 '/path/to/hic/*_R1*.fastq.gz' \
  --hic_r2 '/path/to/hic/*_R2*.fastq.gz' \
  --sample pfas \
  --compleasm_library /path/to/compleasm_lineages \
  --juicer_tools_jar /path/to/juicer_tools_1.22.01.jar
```

When `--juicer_tools_jar` is omitted, the pipeline downloads version 1.22.01
to `<outdir>/00_tools/`. This requires internet access in that task.

## HiFi-Only Run

Hi-C is optional. Without `--hic_r1` and `--hic_r2`, the workflow stops after
the hifiasm draft statistics and Compleasm checks:

```bash
nextflow run . \
  -profile conda \
  --reads '/path/to/hifi_reads/*.fastq.gz' \
  --sample pfas \
  --run_compleasm_on total
```

## Key Parameters

| Parameter | Default | Meaning |
| --- | --- | --- |
| `--reads` | required | HiFi FASTQ/FASTQ.GZ path or glob |
| `--hic_r1` | unset | Hi-C R1 FASTQ/FASTQ.GZ path or glob |
| `--hic_r2` | unset | Hi-C R2 FASTQ/FASTQ.GZ path or glob |
| `--sample` | `hifi_assembly` | output filename prefix |
| `--outdir` | `results` | output directory |
| `--hifiasm_extra` | `-l 2` | additional hifiasm arguments |
| `--yahs_mapq` | `10` | YaHS minimum mapping quality |
| `--yahs_extra` | unset | additional YaHS arguments |
| `--lineage` | `sauropsida` | Compleasm lineage |
| `--odb` | `odb12` | Compleasm OrthoDB release |
| `--compleasm_library` | unset | local Compleasm database library |
| `--download_lineage` | `false` | download the lineage in the workflow |
| `--run_compleasm_on` | `total,scaffolded` | any of `hap1,hap2,total,scaffolded,all` |
| `--juicer_tools_jar` | unset | local Juicer Tools jar |
| `--juicer_tools_url` | official 1.22.01 jar | download URL used when no jar is supplied |
| `--juicer_java_memory` | `48g` | Java heap for Juicer Tools |
| `--make_hic` | `true` | create the final Juicebox `.hic` map |

## Outputs

| Directory | Contents |
| --- | --- |
| `01_hifiasm/` | hap1, hap2, and total primary-contig GFAs |
| `02_draft_fastas/` | draft hap1, hap2, and total FASTAs |
| `03_hic_alignment/` | raw BAM, deduplicated name-sorted BAM, markdup metrics |
| `04_yahs/` | final scaffold FASTA, final AGP, and YaHS binary links when available |
| `05_assembly_stats/` | draft and scaffolded assembly-statistics JSON |
| `06_compleasm/` | selected draft and scaffolded Compleasm results |
| `07_contact_map/` | sorted Juicer input, chromosome sizes, and final `.hic` map |
| `checks/` | machine-readable PASS/FAIL checks for each stage |
| `<sample>.pipeline_report.md` | combined assembly QC report |

The `.hic` file can be opened directly in Juicebox.

## SLURM

The `slurm` profile uses the resource labels in `nextflow.config`:

- `hifiasm`: 32 CPUs, 50 GB, 30 hours
- `hic_mapping`: 25 CPUs, 50 GB, 100 hours
- `yahs`: 1 CPU, 10 GB, 10 hours
- `compleasm`: 16 CPUs, 64 GB, 48 hours
- `contact_map`: 8 CPUs, 50 GB, 48 hours

Adjust these values for the genome size and cluster limits.

## Software

The Conda environment includes hifiasm, Compleasm, BWA, samtools, YaHS,
OpenJDK, and the Compleasm dependencies. Juicer Tools is supplied with
`--juicer_tools_jar` or downloaded by the pipeline.

The environment pins Python to 3.9 because Compleasm depends on `sepp` and
`dendropy` versions that do not currently solve cleanly with newer Python
releases on `linux-64`.

It also pins OpenJDK to 17 because Nextflow requires Java 8 or later but does
not support newer Java releases such as OpenJDK 25.

The Conda profile is the supported complete environment. The generic
Docker/Singularity profiles remain suitable only when their container is
replaced with an image containing the full tool set.

## References

- [hifiasm documentation](https://hifiasm.readthedocs.io/)
- [Compleasm](https://github.com/huangnengCSU/compleasm)
- [YaHS](https://github.com/c-zhou/yahs)
- [Juicer Tools downloads](https://github.com/aidenlab/juicer/wiki/Download)
