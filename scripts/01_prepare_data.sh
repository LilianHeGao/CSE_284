#!/usr/bin/env bash
set -euo pipefail

CFG=${1:-config/project.env}
source "$CFG"

mkdir -p "$RAW_DIR" "$PROC_DIR" "$RESULTS_DIR"

echo "[1/3] Subsetting chromosome ${CHR} from VCF"
"$BCFTOOLS" view -r "$CHR" -Oz -o "$RAW_DIR/chr${CHR}.vcf.gz" "$VCF_GZ"
"$BCFTOOLS" index -f "$RAW_DIR/chr${CHR}.vcf.gz"

echo "[2/3] Converting VCF to PLINK2 pfile"
"$PLINK2" \
  --vcf "$RAW_DIR/chr${CHR}.vcf.gz" \
  --set-all-var-ids @:#:\$r:\$a \
  --make-pgen \
  --out "$PROC_DIR/chr${CHR}"

echo "[3/3] Applying basic QC and exporting bed/bim/fam"
"$PLINK2" \
  --pfile "$PROC_DIR/chr${CHR}" \
  --maf "$MAF" \
  --geno "$GENO_MISSING" \
  --make-bed \
  --out "$PROC_DIR/chr${CHR}.qc"

echo "Done: $PROC_DIR/chr${CHR}.qc.{bed,bim,fam}"
