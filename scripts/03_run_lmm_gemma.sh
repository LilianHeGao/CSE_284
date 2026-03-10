#!/usr/bin/env bash
set -euo pipefail

CFG=${1:-config/project.env}
source "$CFG"

if [[ -f "${CFG}.local" ]]; then
  source "${CFG}.local"
fi

mkdir -p "$RESULTS_DIR/lmm"

RUN_LABEL=${RUN_LABEL:-${SUBSET_NAME:-all}}
BFILE=${BFILE_PREFIX:-"$PROC_DIR/chr${CHR}.${RUN_LABEL}.qc"}
KINSHIP_PREFIX="${RUN_LABEL}_kinship"
PCA_PREFIX="$RESULTS_DIR/lmm/${RUN_LABEL}_pca"
LMM_COVARS="$RESULTS_DIR/lmm/${RUN_LABEL}_lmm_covars.txt"
LMM_PREFIX="${RUN_LABEL}_lmm"
KINSHIP_OUT="$RESULTS_DIR/lmm/${KINSHIP_PREFIX}.cXX.txt"
LMM_OUT="$RESULTS_DIR/lmm/${LMM_PREFIX}.assoc.txt"

echo "[1/3] Building kinship matrix"
"$GEMMA" -bfile "$BFILE" -gk 1 -o "$KINSHIP_PREFIX"
mv "output/${KINSHIP_PREFIX}.cXX.txt" "$KINSHIP_OUT"

echo "[2/3] Building PCs for fixed-effect sensitivity (optional for LMM)"
"$PLINK2" --bfile "$BFILE" --pca "$NUM_PCS" --out "$PCA_PREFIX"
awk 'NR==1 {next} {for(i=3;i<=NF;i++) printf "%s%s", $i, (i==NF?"\n":"\t")}' \
  "$PCA_PREFIX.eigenvec" > "$LMM_COVARS"

echo "[3/3] LMM association test"
"$GEMMA" \
  -bfile "$BFILE" \
  -k "$KINSHIP_OUT" \
  -lmm 4 \
  -c "$LMM_COVARS" \
  -o "$LMM_PREFIX"
mv "output/${LMM_PREFIX}.assoc.txt" "$LMM_OUT"

echo "Done: $LMM_OUT"
