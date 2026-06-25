#!/usr/bin/env python3

import argparse
import json
from pathlib import Path


def table_from_tsv(path):
    lines = [line.rstrip("\n").split("\t") for line in Path(path).read_text().splitlines()]
    if len(lines) < 2:
        return ""
    header, rows = lines[0], lines[1:]
    out = ["| " + " | ".join(header) + " |", "| " + " | ".join(["---"] * len(header)) + " |"]
    out.extend("| " + " | ".join(row) + " |" for row in rows)
    return "\n".join(out)


def stats_summary(paths):
    rows = []
    for path in paths:
        asm_type = Path(path).name.split(".")[-3]
        with open(path) as handle:
            stats = json.load(handle)["Scaffold Stats"]
        rows.append(
            [
                asm_type,
                str(stats["sequence_count"]),
                str(stats["total_bps"]),
                str(stats["N50"]),
                str(stats["L50"]),
                str(stats["longest"]),
                f"{float(stats['gc_content']):.2f}",
            ]
        )
    header = ["assembly", "sequences", "total_bps", "N50", "L50", "longest", "gc_pct"]
    out = ["| " + " | ".join(header) + " |", "| " + " | ".join(["---"] * len(header)) + " |"]
    out.extend("| " + " | ".join(row) + " |" for row in sorted(rows))
    return "\n".join(out)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--hifiasm-checks", nargs="*", default=[])
    parser.add_argument("--fasta-checks", nargs="*", default=[])
    parser.add_argument("--assembly-stats", nargs="*", default=[])
    parser.add_argument("--compleasm-checks", nargs="*", default=[])
    parser.add_argument("--hic-bam-checks", nargs="*", default=[])
    parser.add_argument("--yahs-checks", nargs="*", default=[])
    parser.add_argument("--contact-map-checks", nargs="*", default=[])
    args = parser.parse_args()

    print(f"# {args.sample} Assembly Pipeline Report")
    print()
    print("## Hifiasm GFA Checks")
    for path in sorted(args.hifiasm_checks):
        print()
        print(table_from_tsv(path))

    print()
    print("## FASTA Checks")
    for path in sorted(args.fasta_checks):
        print()
        print(table_from_tsv(path))

    print()
    print("## Hi-C Alignment and Duplicate Removal")
    if args.hic_bam_checks:
        for path in sorted(args.hic_bam_checks):
            print()
            print(table_from_tsv(path))
    else:
        print()
        print("No Hi-C alignment was run.")

    print()
    print("## YaHS Scaffolding")
    if args.yahs_checks:
        for path in sorted(args.yahs_checks):
            print()
            print(table_from_tsv(path))
    else:
        print()
        print("No YaHS scaffolding was run.")

    print()
    print("## Assembly Stats")
    print()
    print(stats_summary(args.assembly_stats))

    print()
    print("## Compleasm")
    if args.compleasm_checks:
        for path in sorted(args.compleasm_checks):
            print()
            print(table_from_tsv(path))
    else:
        print()
        print("No Compleasm summaries were produced.")

    print()
    print("## Hi-C Contact Map")
    if args.contact_map_checks:
        for path in sorted(args.contact_map_checks):
            print()
            print(table_from_tsv(path))
    else:
        print()
        print("No Hi-C contact map was produced.")


if __name__ == "__main__":
    main()
