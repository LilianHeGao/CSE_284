#!/usr/bin/env bash
set -euo pipefail

CFG=${1:-config/project.env}
source "$CFG"

if [[ -f "${CFG}.local" ]]; then
  source "${CFG}.local"
fi

mkdir -p "$RAW_DIR" "$PROC_DIR" "$RESULTS_DIR"

RUN_LABEL=${RUN_LABEL:-${SUBSET_NAME:-all}}
PHENO_SOURCE_CSV=${PHENO_SOURCE_CSV:-"$RAW_DIR/values.csv"}
PHENO_FILE=${PHENO_FILE:-"$PROC_DIR/${RUN_LABEL}.pheno.txt"}
KEEP_IDS_FILE=${KEEP_IDS_FILE:-"$PROC_DIR/${RUN_LABEL}.keep_ids.txt"}
PHENO_SAMPLES_FILE=${PHENO_SAMPLES_FILE:-"$PROC_DIR/${RUN_LABEL}.samples.tsv"}
PFILE_PREFIX=${PFILE_PREFIX:-"$PROC_DIR/chr${CHR}.${RUN_LABEL}"}
BFILE_PREFIX=${BFILE_PREFIX:-"$PROC_DIR/chr${CHR}.${RUN_LABEL}.qc"}
CLEAN_VCF_GZ=${CLEAN_VCF_GZ:-"$RAW_DIR/1001G.Chr${CHR}_clean.vcf.gz"}

echo "[1/3] Building phenotype and keep files for subset '$RUN_LABEL'"
python scripts/00_build_subset_inputs.py \
  --values-csv "$PHENO_SOURCE_CSV" \
  --subset-name "${SUBSET_NAME:-all}" \
  --countries "${SUBSET_COUNTRIES:-}" \
  --pheno-name "$PHENO_NAME" \
  --out-pheno "$PHENO_FILE" \
  --out-keep "$KEEP_IDS_FILE" \
  --out-samples "$PHENO_SAMPLES_FILE"

echo "[2/3] Normalizing VCF input"
VCF_FOR_PLINK="$CLEAN_VCF_GZ"
if command -v "$BCFTOOLS" >/dev/null 2>&1; then
  "$BCFTOOLS" view -v snps -m2 -M2 "$VCF_GZ" -Oz -o "$CLEAN_VCF_GZ"
else
  echo "bcftools not found; sanitizing VCF with Python fallback."
  python scripts/sanitize_vcf.py --in-vcf "$VCF_GZ" --out-vcf "$CLEAN_VCF_GZ"
fi

echo "[3/3] Converting VCF to PLINK and applying QC"
"$PLINK2" \
  --vcf "$VCF_FOR_PLINK" \
  --keep "$KEEP_IDS_FILE" \
  --chr-set 5 \
  --allow-extra-chr \
  --set-all-var-ids @:#:\$r:\$a \
  --make-pgen \
  --out "$PFILE_PREFIX"

"$PLINK2" \
  --pfile "$PFILE_PREFIX" \
  --maf "$MAF" \
  --geno "$GENO_MISSING" \
  --make-bed \
  --out "$BFILE_PREFIX"
