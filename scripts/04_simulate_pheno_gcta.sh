#!/usr/bin/env bash
set -euo pipefail

CFG=${1:-config/project.env}
source "$CFG"

if [[ -f "${CFG}.local" ]]; then
  source "${CFG}.local"
fi

mkdir -p "$RESULTS_DIR/sim"

RUN_LABEL=${RUN_LABEL:-${SUBSET_NAME:-all}}
BFILE=${BFILE_PREFIX:-"$PROC_DIR/chr${CHR}.${RUN_LABEL}.qc"}

# Create a random list of causal SNPs from BIM (example: 50 causal SNPs)
shuf "$BFILE.bim" | head -n 50 | awk '{print $2}' > "$RESULTS_DIR/sim/causal.snplist"

echo "Simulating quantitative trait with h2=0.5 (single replicate)"
"$GCTA" \
  --bfile "$BFILE" \
  --simu-qt \
  --simu-causal-loci "$RESULTS_DIR/sim/causal.snplist" \
  --simu-hsq 0.5 \
  --simu-rep 1 \
  --out "$RESULTS_DIR/sim/sim"

echo "GCTA outputs written to $RESULTS_DIR/sim"
