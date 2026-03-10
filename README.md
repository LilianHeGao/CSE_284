# Plant GWAS Benchmark: Standard Linear Regression vs. Linear Mixed Models

This repository implements a benchmark comparing SNP-wise GWAS using:
- Naive linear regression (`LR`) with no population-structure correction
- Linear regression + genotype PCs as fixed covariates (`LR+PCs`)
- Linear mixed model (`LMM`) with kinship random effect

Dataset target: *Arabidopsis thaliana* 1001 Genomes (initial focus: chromosome 4).

## Why this project
Population structure and relatedness can inflate GWAS p-values. We compare calibration, signal consistency, and runtime between LR and LMM approaches on real plant genotype structure.

## Tools (addressing proposal feedback)
- `bcftools`: chromosome subsetting and variant filtering in VCF
- `plink2`: genotype QC, PCA, and SNP-wise linear regression GWAS
- `GEMMA`: kinship matrix + LMM association testing
- `gcta64` (optional simulation backend): quantitative phenotype simulation from real genotypes
- `python` (`pandas`, `numpy`, `matplotlib`, `scipy`): QQ/Manhattan plots, lambda GC, top-hit overlap, runtime summary

## Repository layout
- `scripts/01_prepare_data.sh`: subset chr4 + convert to PLINK + light QC
- `scripts/02_run_lr.sh`: naive LR and LR+PC GWAS (PLINK2)
- `scripts/03_run_lmm_gemma.sh`: kinship and LMM GWAS (GEMMA)
- `scripts/04_simulate_pheno_gcta.sh`: simulate quantitative phenotype (GCTA template)
- `scripts/05_evaluate.py`: QQ/Manhattan, lambda GC, top-hit overlap
- `slurm/run_pipeline.slurm`: HPC wrapper for end-to-end run
- `config/project.env.example`: central paths and parameters

## Quick start
1. Copy config and edit paths:
```bash
cp config/project.env.example config/project.env
```
2. Run stepwise:
```bash
bash scripts/01_prepare_data.sh config/project.env
bash scripts/02_run_lr.sh config/project.env
bash scripts/03_run_lmm_gemma.sh config/project.env
python scripts/05_evaluate.py \
  --lr results/lr/lr.PHENO1.glm.linear \
  --lr-pcs results/lr_pcs/lr_pcs.PHENO1.glm.linear \
  --lmm results/lmm/lmm.assoc.txt \
  --out-prefix results/plots/benchmark
```
3. Optional simulation:
```bash
bash scripts/04_simulate_pheno_gcta.sh config/project.env
```

## Windows cmd.exe
If you are running from a Windows `cmd.exe` or Conda prompt, use the Python runner instead of `bash`:
```bat
python scripts\run_pipeline.py --config config\project.env check-tools
python scripts\run_pipeline.py --config config\project.env prepare-data
python scripts\run_pipeline.py --config config\project.env run-lr
python scripts\run_pipeline.py --config config\project.env run-lmm
python scripts\run_pipeline.py --config config\project.env evaluate
```
`plink2` and `python` are required. `bcftools` (VCF normalization) and `gemma` (LMM) are optional; if `gemma` is missing, run `evaluate` after `run-lr` to compare LR vs. LR+PCs only.
When `bcftools` is missing, the pipeline automatically runs `scripts/sanitize_vcf.py` to drop malformed VCF lines before PLINK import.

## Results so far
- Pipeline scaffold and reproducible scripts are now in place.
- Methods and toolchain are fully specified (no longer vague).
- Full benchmark outputs are pending first data run on chr4 subset.

## Evaluation outputs planned
- QQ plots for LR, LR+PCs, and LMM
- Manhattan plots for each method
- Genomic inflation factor (`lambda_GC`) per method
- Top-hit overlap and effect/p-value concordance tables
- Runtime comparison table

## Remaining work
1. Finalize exact accession subset and phenotype file harmonization.
2. Run full chr4 benchmark and collect first-pass plots/tables.
3. Tune variant QC thresholds (MAF/missingness) and assess sensitivity.
4. Validate simulation settings (heritability, number of causal SNPs, effect-size distribution).
5. Decide whether to include a second trait/chromosome for robustness.

## Challenges
- Best practice for choosing number of PCs in LR+PC baseline.
- Fair runtime comparison settings between PLINK2 and GEMMA.
- Recommended simulation design for realistic LD-aware causal architectures.

## Reference
Alonso-Blanco, Carlos, et al. "1,135 genomes reveal the global pattern of polymorphism in Arabidopsis thaliana." *Cell* 166.2 (2016): 481-491.
