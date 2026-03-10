#!/usr/bin/env python3
from __future__ import annotations

import argparse
import gzip
from pathlib import Path


def open_text(path: Path, mode: str):
    if path.suffix == ".gz":
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return path.open(mode, encoding="utf-8", newline="")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Remove malformed VCF records with unexpected column counts."
    )
    parser.add_argument("--in-vcf", required=True)
    parser.add_argument("--out-vcf", required=True)
    parser.add_argument("--max-log", type=int, default=5)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    in_path = Path(args.in_vcf)
    out_path = Path(args.out_vcf)

    if not in_path.exists():
        raise SystemExit(f"Input VCF not found: {in_path}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    expected_cols = None
    total_records = 0
    kept_records = 0
    dropped_records = 0
    logged = 0

    with open_text(in_path, "rt") as src, open_text(out_path, "wt") as dst:
        for line_no, raw in enumerate(src, start=1):
            if raw.startswith("##"):
                dst.write(raw)
                continue

            if raw.startswith("#CHROM"):
                header = raw.rstrip("\r\n").split("\t")
                expected_cols = len(header)
                dst.write(raw)
                continue

            if expected_cols is None:
                raise SystemExit("VCF header (#CHROM) not found before variant lines.")

            total_records += 1
            row = raw.rstrip("\r\n").split("\t")
            if len(row) != expected_cols:
                dropped_records += 1
                if logged < args.max_log:
                    print(
                        f"Dropping malformed line {line_no}: "
                        f"expected {expected_cols} columns, found {len(row)}"
                    )
                    logged += 1
                continue

            dst.write(raw)
            kept_records += 1

    print(f"input={in_path}")
    print(f"output={out_path}")
    print(f"variants_total={total_records}")
    print(f"variants_kept={kept_records}")
    print(f"variants_dropped={dropped_records}")


if __name__ == "__main__":
    main()
