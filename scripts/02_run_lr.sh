#!/usr/bin/env bash
set -euo pipefail

CFG=${1:-config/project.env}
source "$CFG"

if [[ -f "${CFG}.local" ]]; then
  source "${CFG}.local"
fi

mkdir -p "$RESULTS_DIR/lr" "$RESULTS_DIR/lr_pcs"

RUN_LABEL=${RUN_LABEL:-${SUBSET_NAME:-all}}
BFILE=${BFILE_PREFIX:-"$PROC_DIR/chr${CHR}.${RUN_LABEL}.qc"}
LR_PREFIX="$RESULTS_DIR/lr/${RUN_LABEL}_lr"
PCA_PREFIX="$RESULTS_DIR/lr_pcs/${RUN_LABEL}_pca"
PCA_FREQ_PREFIX="$RESULTS_DIR/lr_pcs/${RUN_LABEL}_pca_freq"
LR_PCS_PREFIX="$RESULTS_DIR/lr_pcs/${RUN_LABEL}_lr_pcs"

echo "[1/3] Naive LR GWAS (no structure correction)"
"$PLINK2" \
  --bfile "$BFILE" \
  --pheno "$PHENO_FILE" \
  --pheno-name "$PHENO_NAME" \
  --glm hide-covar allow-no-covars \
  --out "$LR_PREFIX"

echo "[2/3] PCA for LR+PC covariates"
"$PLINK2" \
  --bfile "$BFILE" \
  --freq \
  --out "$PCA_FREQ_PREFIX"

"$PLINK2" \
  --bfile "$BFILE" \
  --pca "$NUM_PCS" \
  --read-freq "$PCA_FREQ_PREFIX.afreq" \
  --out "$PCA_PREFIX"

echo "[3/3] LR+PC GWAS"
"$PLINK2" \
  --bfile "$BFILE" \
  --pheno "$PHENO_FILE" \
  --pheno-name "$PHENO_NAME" \
  --covar "$PCA_PREFIX.eigenvec" \
  --glm hide-covar \
  --out "$LR_PCS_PREFIX"

echo "Done: LR and LR+PC results written to $RESULTS_DIR"
