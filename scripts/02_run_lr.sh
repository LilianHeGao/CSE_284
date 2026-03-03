#!/usr/bin/env bash
set -euo pipefail

CFG=${1:-config/project.env}
source "$CFG"

mkdir -p "$RESULTS_DIR/lr" "$RESULTS_DIR/lr_pcs"

BFILE="$PROC_DIR/chr${CHR}.qc"

echo "[1/3] Naive LR GWAS (no structure correction)"
"$PLINK2" \
  --bfile "$BFILE" \
  --pheno "$PHENO_FILE" \
  --pheno-name "$PHENO_NAME" \
  --glm hide-covar \
  --out "$RESULTS_DIR/lr/lr"

echo "[2/3] PCA for LR+PC covariates"
"$PLINK2" \
  --bfile "$BFILE" \
  --pca "$NUM_PCS" \
  --out "$RESULTS_DIR/lr_pcs/pca"

awk 'BEGIN{OFS="\t"} NR==1{$1="IID"; print; next} {print $2, $3, $4, $5, $6, $7}' \
  "$RESULTS_DIR/lr_pcs/pca.eigenvec" > "$RESULTS_DIR/lr_pcs/covars.tsv"

echo "[3/3] LR+PC GWAS"
"$PLINK2" \
  --bfile "$BFILE" \
  --pheno "$PHENO_FILE" \
  --pheno-name "$PHENO_NAME" \
  --covar "$RESULTS_DIR/lr_pcs/covars.tsv" \
  --glm hide-covar \
  --out "$RESULTS_DIR/lr_pcs/lr_pcs"

echo "Done: LR and LR+PC results written to $RESULTS_DIR"
