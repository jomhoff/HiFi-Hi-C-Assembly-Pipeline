#!/usr/bin/env python3

import argparse
import re
import sys


RESULT_RE = re.compile(
    r"C:(?P<C>[0-9.]+)%.*?S:(?P<S>[0-9.]+)%.*?D:(?P<D>[0-9.]+)%.*?"
    r"F:(?P<F>[0-9.]+)%.*?M:(?P<M>[0-9.]+)%.*?n:(?P<n>[0-9]+)"
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--summary", required=True)
    args = parser.parse_args()

    with open(args.summary) as handle:
        text = handle.read()

    match = RESULT_RE.search(text.replace("\n", " "))
    if not match:
        raise SystemExit(f"Could not parse Compleasm summary: {args.summary}")

    fields = match.groupdict()
    sys.stdout.write(
        "\t".join(
            [
                "sample",
                "assembly",
                "complete_pct",
                "single_copy_pct",
                "duplicated_pct",
                "fragmented_pct",
                "missing_pct",
                "orthologs",
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
                fields["C"],
                fields["S"],
                fields["D"],
                fields["F"],
                fields["M"],
                fields["n"],
                "PASS",
            ]
        )
        + "\n"
    )


if __name__ == "__main__":
    main()
