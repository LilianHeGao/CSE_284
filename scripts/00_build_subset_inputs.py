#!/usr/bin/env python3
import argparse
import csv
from collections import Counter
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Build phenotype and keep files for a configurable accession subset."
    )
    parser.add_argument("--values-csv", required=True)
    parser.add_argument("--subset-name", default="all")
    parser.add_argument("--countries", default="")
    parser.add_argument("--pheno-name", default="PHENO")
    parser.add_argument("--out-pheno", required=True)
    parser.add_argument("--out-keep", required=True)
    parser.add_argument("--out-samples", required=True)
    return parser.parse_args()


def parse_countries(raw: str) -> set[str]:
    return {country.strip() for country in raw.split(",") if country.strip()}


def sort_key(row: dict[str, str]) -> tuple[int, object]:
    iid = row["IID"]
    if iid.isdigit():
        return (0, int(iid))
    return (1, iid)


def write_tsv(path: Path, header: list[str], rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        writer.writerows(rows)


def main() -> None:
    args = parse_args()
    subset_name = args.subset_name.strip() or "all"
    allowed_countries = parse_countries(args.countries)

    if subset_name != "all" and not allowed_countries:
        raise SystemExit("--countries is required when --subset-name is not 'all'")

    samples: dict[str, dict[str, object]] = {}
    duplicates = 0

    with Path(args.values_csv).open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            accession_id = row["accession_id"].strip()
            country = row["accession_country"].strip()
            phenotype_value = row["phenotype_value"].strip()

            if not accession_id or not phenotype_value:
                continue

            if subset_name != "all" and country not in allowed_countries:
                continue

            try:
                phenotype_float = float(phenotype_value)
            except ValueError as exc:
                raise SystemExit(
                    f"Non-numeric phenotype value for accession {accession_id}: {phenotype_value}"
                ) from exc

            if accession_id not in samples:
                samples[accession_id] = {
                    "FID": accession_id,
                    "IID": accession_id,
                    "COUNTRY": country,
                    "values": [phenotype_float],
                }
                continue

            duplicates += 1
            sample = samples[accession_id]
            if sample["COUNTRY"] != country:
                raise SystemExit(
                    f"Accession {accession_id} has conflicting countries: "
                    f"{sample['COUNTRY']} vs {country}"
                )
            sample["values"].append(phenotype_float)

    ordered = []
    for sample in samples.values():
        values = sample.pop("values")
        mean_value = sum(values) / len(values)
        sample[args.pheno_name] = f"{mean_value:g}"
        ordered.append(sample)

    ordered = sorted(ordered, key=sort_key)
    if not ordered:
        raise SystemExit("No samples matched the configured subset.")

    pheno_rows = [[row["FID"], row["IID"], row[args.pheno_name]] for row in ordered]
    keep_rows = [[row["FID"], row["IID"]] for row in ordered]
    sample_rows = [[row["FID"], row["IID"], row["COUNTRY"], row[args.pheno_name]] for row in ordered]

    write_tsv(Path(args.out_pheno), ["FID", "IID", args.pheno_name], pheno_rows)

    out_keep = Path(args.out_keep)
    out_keep.parent.mkdir(parents=True, exist_ok=True)
    with out_keep.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerows(keep_rows)

    write_tsv(Path(args.out_samples), ["FID", "IID", "COUNTRY", args.pheno_name], sample_rows)

    country_counts = Counter(row["COUNTRY"] for row in ordered)
    print(f"subset={subset_name}")
    print(f"samples={len(ordered)}")
    print(f"duplicate_rows_aggregated={duplicates}")
    for country, count in sorted(country_counts.items()):
        print(f"{country}\t{count}")


if __name__ == "__main__":
    main()
