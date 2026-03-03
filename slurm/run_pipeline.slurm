#!/usr/bin/env bash
#SBATCH --job-name=cse284_gwas
#SBATCH --partition=hotel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=08:00:00
#SBATCH --output=results/pipeline_%j.out
#SBATCH --error=results/pipeline_%j.err

set -euo pipefail

bash scripts/00_check_tools.sh
bash scripts/01_prepare_data.sh config/project.env
bash scripts/02_run_lr.sh config/project.env
bash scripts/03_run_lmm_gemma.sh config/project.env

python scripts/05_evaluate.py \
  --lr results/lr/lr.PHENO1.glm.linear \
  --lr-pcs results/lr_pcs/lr_pcs.PHENO1.glm.linear \
  --lmm results/lmm/lmm.assoc.txt \
  --out-prefix results/plots/benchmark
