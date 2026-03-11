# Plant GWAS Benchmark Report

## Overview

This report summarizes the completed benchmark in `/tscc/projects/ps-gaultonlab/junxif/CSE_284`, comparing three association strategies on *Arabidopsis thaliana* chromosome 4 genotype data:

1. Naive linear regression (`LR`)
2. Linear regression with genotype principal components as fixed covariates (`LR+PCs`)
3. Linear mixed model with a kinship matrix (`LMM`)

The main goal of the project was to evaluate how progressively stronger correction for relatedness and population structure changes GWAS calibration and agreement of top association signals.

## Data Source and Provenance

The project used a chromosome 4 subset derived from the 1001 Genomes *A. thaliana* resource. The local input file was:

- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/data/36a6dac7ec9b5de7f6bee054f3a59ae9_fullgenome.vcf`

During preprocessing, this file was found to contain a valid VCF prefix followed by embedded HTML from a remote `500 Internal Server Error` response. The valid portion ended at line 115898, after which the file was no longer parseable as VCF. To complete the analysis reproducibly, the intact prefix was preserved as:

- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/data/raw/chr4_clean.vcf`

This salvage step was necessary because standard tools (`bcftools`, `plink2`) could not process the corrupted tail of the original download. The cleaned file still contained only chromosome 4 records and was sufficient for the benchmark run.

## Phenotype Strategy

No curated real phenotype file was present in the project at run time, so the benchmark was completed using a simulated quantitative phenotype generated from the filtered genotype matrix with GCTA. This means the current results measure method behavior in a controlled simulation setting rather than recovery on a real biological trait.

Simulation settings:

- Tool: `gcta64`
- Genotype basis: post-QC chromosome 4 PLINK bed/bim/fam
- Number of causal SNPs: 50
- Heritability (`h2`): 0.5
- Replicates: 1

Output phenotype files:

- PLINK-ready phenotype: `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/sim/sim_pheno.tsv`
- GEMMA-ready phenotype: `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/sim/sim_gemma.txt`

## Genotype Processing Methods

Genotype preparation was performed in three stages:

1. The cleaned chromosome 4 VCF was compressed and indexed with `bcftools`.
2. The VCF was converted to PLINK2 pgen/pvar/psam format.
3. Basic QC was applied in PLINK2.

QC thresholds:

- Minor allele frequency (`MAF`) >= 0.05
- Variant missingness (`GENO_MISSING`) <= 0.05
- Restriction to biallelic variants (`--max-alleles 2`)

Observed filtering outcome:

- Samples retained: 1135
- Variants before QC import: 115882
- Variants remaining after restricting to biallelic sites and QC: 461

The large reduction in site count indicates that the chosen QC thresholds were stringent relative to the local chromosome 4 subset. This is not inherently wrong, but it means the final benchmark was run on a compact marker set.

## Association Methods

### 1. Naive Linear Regression

The first baseline used PLINK2 linear regression with no population-structure covariates:

- Tool: `plink2`
- Model: additive single-SNP linear regression
- Phenotype: simulated quantitative trait

This method intentionally omits structure correction and serves as a calibration baseline.

### 2. Linear Regression with Principal Components

The second baseline used the same PLINK2 regression framework, but included the top 5 genotype principal components as fixed-effect covariates:

- Tool: `plink2`
- PCA count: 5
- Covariate source: PLINK2 `--pca`
- Model: additive single-SNP linear regression with PC adjustment

This method is intended to reduce confounding from broad ancestry or structure while remaining simpler than a mixed model.

### 3. Linear Mixed Model

The third method used GEMMA:

- Tool: `gemma`
- Kinship: centered relatedness matrix (`-gk 1`)
- Association test: LMM (`-lmm 4`)
- Covariates: top 5 PCs
- Phenotype: simulated quantitative trait

This approach models polygenic background and relatedness through the kinship matrix, making it the most statistically conservative method of the three.

## Evaluation Metrics

The evaluation script compared outputs from all three methods using:

- QQ plot
- Manhattan plot
- Genomic inflation factor (`lambda_GC`)
- Overlap of top 100 associated SNPs between methods

Output files:

- QQ plot: `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/benchmark.qq.png`
- Manhattan plot: `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/benchmark.manhattan.png`
- Inflation summary: `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/benchmark.lambda_gc.tsv`
- Top-hit overlap: `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/benchmark.top_overlap.tsv`

Each method produced 461 tested variants, matching the final QC-filtered marker count.

## Observed Results

### Genomic Inflation

The most important result in this benchmark is the genomic inflation pattern:

- `LR`: 19.2501
- `LR+PCs`: 11.6591
- `LMM`: 0.9209

Interpretation:

- The naive linear regression is extremely inflated, indicating severe confounding or structure-induced signal distortion.
- Adding 5 PCs reduces the inflation materially, but the result remains heavily inflated.
- The mixed model is close to the ideal range around 1.0 and is therefore substantially better calibrated than either linear-regression baseline.

This is the central statistical result of the run. In this dataset, broad PC adjustment helps, but it does not solve the confounding nearly as effectively as the kinship-based mixed model.

### Top-Hit Overlap

Overlap among the top 100 SNPs:

- `LMM` vs `LR`: 30
- `LMM` vs `LR+PCs`: 49
- `LR` vs `LR+PCs`: 32

Interpretation:

- `LR+PCs` is more similar to `LMM` than naive `LR` is, which is the expected direction.
- Even so, the overlap between `LR+PCs` and `LMM` is only 49 out of 100, indicating that the choice of association model meaningfully changes which SNPs rank as top hits.
- The relatively low `LR` vs `LMM` overlap is consistent with the strong inflation seen in naive regression.

### Mixed-Model Variance Explained

GEMMA reported:

- `pve estimate = 0.509036`
- `se(pve) = 0.0527948`

This is very close to the simulated trait heritability of 0.5, which is a useful internal consistency check. It suggests that the simulation and mixed-model fitting are aligned as expected.

### GEMMA Eigenvalue Warning

GEMMA reported:

- `Matrix G has 449 eigenvalues close to zero`

Given that only 461 variants remained after QC, this is not surprising. The kinship matrix is being estimated from a relatively small marker set, so near-singularity is expected. The model still completed successfully, but this warning reinforces that the run should be interpreted as a proof-of-pipeline benchmark rather than a fully powered GWAS.

## Interpretation

The observed results support the standard expectation for structured genotype data:

- Uncorrected linear regression is badly miscalibrated.
- Principal component adjustment improves calibration but remains insufficient.
- A kinship-based linear mixed model is much better calibrated and more defensible for inference.

In this run, the mixed model not only reduced inflation dramatically, but also produced a materially different ranking of top signals. This indicates that the confounding present in the genotype structure is strong enough to alter apparent discoveries under simpler methods.

## Practical Limitations

Several limitations should be kept in mind:

1. The phenotype was simulated, not biological.
2. Only a chromosome 4 subset was analyzed.
3. The usable VCF had to be salvaged from a corrupted download.
4. Only 461 variants survived QC, which is a small marker panel for kinship estimation and GWAS benchmarking.
5. The project currently uses a single simulation replicate, so the reported result pattern could vary under repeated simulation.

These limitations do not invalidate the result, but they define its scope. The current benchmark is best interpreted as a validated working pipeline plus a first-pass demonstration of calibration differences across methods.

## Recommended Next Steps

The strongest follow-up improvements would be:

1. Replace the simulated phenotype with a real accession-matched trait.
2. Re-download or regenerate a clean chromosome 4 VCF from a reliable source.
3. Relax or reassess QC thresholds if the goal is denser marker coverage.
4. Run multiple simulation replicates and summarize mean and variance of `lambda_GC`, top-hit overlap, and runtime.
5. Add a second chromosome or second trait to test robustness.

## Conclusion

This project now runs successfully from local input through analysis and figure generation. The main observed result is clear: for this chromosome 4 benchmark, naive linear regression is severely inflated, PC-adjusted regression improves but remains miscalibrated, and the mixed model provides near-ideal genomic control and the most credible association profile.

That makes the current run a successful proof of concept for the original project question: population structure and relatedness matter substantially, and a mixed-model GWAS is much better behaved than simpler regression baselines on these data.
