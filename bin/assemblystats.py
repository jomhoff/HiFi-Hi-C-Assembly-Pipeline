#!/usr/bin/env python3

import json
import sys
from itertools import groupby


def fasta_iter(fasta_file):
    with open(fasta_file) as fh:
        fa_iter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
        for header in fa_iter:
            header = next(header)[1:].strip()
            seq = "".join(s.upper().strip() for s in next(fa_iter))
            yield header, seq


def read_genome(fasta_file):
    gc = 0
    total_len = 0
    contig_lens = []
    scaffold_lens = []

    for _, seq in fasta_iter(fasta_file):
        scaffold_lens.append(len(seq))
        contig_list = seq.split("NN") if "NN" in seq else [seq]
        for contig in contig_list:
            if contig:
                gc += contig.count("G") + contig.count("C")
                total_len += len(contig)
                contig_lens.append(len(contig))

    if not contig_lens or not scaffold_lens or total_len == 0:
        raise ValueError(f"No sequence bases found in {fasta_file}")

    gc_cont = (gc / total_len) * 100
    return contig_lens, scaffold_lens, gc_cont


def calculate_stats(seq_lens, gc_cont):
    stats = {}
    sorted_lens = sorted(seq_lens, reverse=True)
    sequence_count = len(sorted_lens)

    stats["sequence_count"] = sequence_count
    stats["gc_content"] = float(gc_cont)
    stats["longest"] = int(sorted_lens[0])
    stats["shortest"] = int(sorted_lens[-1])
    mid = sequence_count // 2
    if sequence_count % 2:
        stats["median"] = float(sorted_lens[mid])
    else:
        stats["median"] = float((sorted_lens[mid - 1] + sorted_lens[mid]) / 2)
    stats["total_bps"] = int(sum(sorted_lens))
    stats["mean"] = float(stats["total_bps"] / sequence_count)

    csum = []
    running_total = 0
    for length in sorted_lens:
        running_total += length
        csum.append(running_total)

    for level in [10, 20, 30, 40, 50]:
        nx = int(stats["total_bps"] * (level / 100))
        hit_index = next(index for index, value in enumerate(csum) if value >= nx)
        stats[f"L{level}"] = hit_index
        stats[f"N{level}"] = int(sorted_lens[hit_index])

    return stats


def main():
    if len(sys.argv) != 2:
        raise SystemExit("Usage: assemblystats.py assembly.fa")

    infilename = sys.argv[1]
    contig_lens, scaffold_lens, gc_cont = read_genome(infilename)
    stat_output = {
        "Contig Stats": calculate_stats(contig_lens, gc_cont),
        "Scaffold Stats": calculate_stats(scaffold_lens, gc_cont),
    }
    print(json.dumps(stat_output, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
