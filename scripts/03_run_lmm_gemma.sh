#!/usr/bin/env bash
set -euo pipefail

CFG=${1:-config/project.env}
source "$CFG"

mkdir -p "$RESULTS_DIR/lmm"
BFILE="$PROC_DIR/chr${CHR}.qc"

echo "[1/3] Building kinship matrix"
"$GEMMA" -bfile "$BFILE" -gk 1 -o kinship
mv output/kinship.cXX.txt "$RESULTS_DIR/lmm/kinship.cXX.txt"

echo "[2/3] Building PCs for fixed-effect sensitivity (optional for LMM)"
"$PLINK2" --bfile "$BFILE" --pca "$NUM_PCS" --out "$RESULTS_DIR/lmm/pca"
awk 'BEGIN{OFS="\t"} {for(i=3;i<=NF;i++) printf "%s%s", $i, (i==NF?"\n":"\t")}' \
  "$RESULTS_DIR/lmm/pca.eigenvec" > "$RESULTS_DIR/lmm/lmm_covars.txt"

echo "[3/3] LMM association test"
"$GEMMA" \
  -bfile "$BFILE" \
  -k "$RESULTS_DIR/lmm/kinship.cXX.txt" \
  -lmm 4 \
  -c "$RESULTS_DIR/lmm/lmm_covars.txt" \
  -o lmm
mv output/lmm.assoc.txt "$RESULTS_DIR/lmm/lmm.assoc.txt"

echo "Done: $RESULTS_DIR/lmm/lmm.assoc.txt"
