#!/usr/bin/env python3

import argparse
import json
import sys


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--stats", required=True)
    args = parser.parse_args()

    with open(args.stats) as handle:
        data = json.load(handle)

    scaffold = data.get("Scaffold Stats", {})
    required = ["sequence_count", "total_bps", "N50", "L50", "gc_content", "longest"]
    missing = [key for key in required if key not in scaffold]
    if missing:
        raise SystemExit(f"Missing scaffold stats in {args.stats}: {', '.join(missing)}")
    if int(scaffold["sequence_count"]) <= 0 or int(scaffold["total_bps"]) <= 0:
        raise SystemExit(f"Invalid zero-sized assembly stats in {args.stats}")

    sys.stdout.write(
        "\t".join(
            [
                "sample",
                "assembly",
                "sequences",
                "total_bps",
                "N50",
                "L50",
                "longest",
                "gc_content",
                "status",
            ]
        )
        + "\n"
    )
    sys.stdout.write(
        "\t".join(
            [
                args.sample,
                args.assembly,
                str(scaffold["sequence_count"]),
                str(scaffold["total_bps"]),
                str(scaffold["N50"]),
                str(scaffold["L50"]),
                str(scaffold["longest"]),
                f"{float(scaffold['gc_content']):.4f}",
                "PASS",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
