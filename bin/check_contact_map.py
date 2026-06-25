#!/usr/bin/env python3

import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--assembly", required=True)
    parser.add_argument("--map", required=True)
    args = parser.parse_args()

    map_path = Path(args.map)
    size = map_path.stat().st_size
    if size <= 0:
        raise SystemExit(f"Empty contact map file: {args.map}")

    print("\t".join(["sample", "assembly", "contact_map", "format", "bytes", "status"]))
    print(
        "\t".join(
            [
                args.sample,
                args.assembly,
                str(map_path),
                map_path.suffix.lstrip("."),
                str(size),
                "PASS",
            ]
        )
    )


if __name__ == "__main__":
    main()
