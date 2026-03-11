# GWAS Benchmark Report on Arabidopsis Chromosome 4

## Abstract

This project benchmarks three genome-wide association strategies on *Arabidopsis thaliana* chromosome 4 genotype data from the 1001 Genomes resource: naive linear regression (`LR`), linear regression with genotype principal components (`LR+PCs`), and a kinship-based linear mixed model (`LMM`). Two phenotype settings were analyzed: a simulated quantitative trait generated from the observed genotype matrix, and a real flowering-time phenotype (`FT16`) provided in `data/values.csv`. The benchmark was designed to evaluate calibration, agreement of top-ranked associations, and the practical effect of population-structure correction. Across both phenotype settings, naive regression was strongly inflated, PC adjustment improved but did not eliminate inflation, and the mixed model remained close to the expected calibration range. The main conclusion is consistent across both analyses: mixed-model correction is materially more reliable than naive or PC-only regression on this structured genotype set.

## 1. Introduction

Genome-wide association studies are sensitive to population structure, relatedness, and cryptic confounding. In structured panels such as *Arabidopsis thaliana* accessions, these effects can induce inflated association statistics and unstable rankings of top loci when simple regression methods are used. The purpose of this project was to implement a reproducible benchmark comparing:

1. single-marker linear regression without correction,
2. linear regression with genotype principal components as fixed covariates, and
3. a linear mixed model with a kinship matrix.

The specific project objective was not only to run these methods, but to compare their empirical behavior on a shared genotype set under two phenotype settings:

- a simulated trait with known polygenic signal,
- a real phenotype derived from the `FT16` measurements.

## 2. Data Sources

### 2.1 Genotype Data

The genotype input originated from a chromosome 4 subset derived from the 1001 Genomes *A. thaliana* release. The local file initially used was:

- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/data/36a6dac7ec9b5de7f6bee054f3a59ae9_fullgenome.vcf`

During preprocessing, this file was found to contain valid VCF records followed by embedded HTML from a remote `500 Internal Server Error` response. The valid portion extended through line 115898. To preserve reproducibility, the intact prefix was extracted and saved as:

- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/data/raw/chr4_clean.vcf`

This cleaned chromosome 4 VCF was used for all downstream analysis.

### 2.2 Real Phenotype Data

The real phenotype input was:

- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/data/values.csv`

This file contained a single phenotype, `FT16`, keyed by `accession_id`. Direct overlap with genotype sample IDs was sufficient for phenotype alignment even though the local `accessions.csv` file was incomplete and could not support a full Assembly_ID to Accession_ID remapping.

Aligned phenotype file:

- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/data/processed/FT16_pheno.tsv`

Alignment summary:

- total genotype samples in QC set: 1135
- matched FT16 phenotype samples: 970
- unmatched genotype samples: 165
- duplicate accession IDs averaged: 1

### 2.3 Simulated Phenotype Data

Because the initial project state lacked a curated phenotype file, a simulated trait was also generated to support method validation:

- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/sim/sim_pheno.tsv`

This simulation remains a useful benchmark setting because it tests whether the implemented pipeline behaves as expected under a phenotype generated directly from the observed genotype matrix.

## 3. Methods

### 3.1 Genotype Preprocessing

Genotype preprocessing was performed using `bcftools` and `plink2`:

1. chromosome 4 VCF compression and indexing,
2. conversion from VCF to PLINK2 pgen/pvar/psam,
3. variant-level QC and export to bed/bim/fam.

QC thresholds:

- minor allele frequency (`MAF`) >= 0.05
- variant missingness <= 0.05
- restriction to biallelic variants (`--max-alleles 2`)

Post-QC genotype summary:

- samples retained in the full genotype set: 1135
- variants retained after QC: 461

These 461 variants constituted the common starting point for all association analyses. For the real-phenotype `LMM`, the genotype matrix was additionally subset to the 970 phenotyped samples before kinship estimation and model fitting.

### 3.2 Simulated Phenotype Generation

The simulated phenotype was generated with GCTA:

- tool: `gcta64`
- genotype basis: post-QC chromosome 4 bed/bim/fam
- number of causal SNPs: 50
- simulated heritability (`h2`): 0.5
- replicates: 1

This produced a quantitative phenotype for all 1135 individuals. The simulation was intended to provide a controlled condition where method calibration could be compared independently of phenotype missingness.

### 3.3 Real Phenotype Preparation

The real `FT16` phenotype file was generated by:

1. selecting rows with `phenotype_name = FT16`,
2. matching `accession_id` directly to genotype `IID`,
3. averaging duplicate phenotype entries for the same accession,
4. writing a PLINK/GEMMA-ready phenotype file in genotype sample order,
5. assigning `-9` to unmatched samples for the PLINK analyses,
6. explicitly subsetting the LMM analysis to phenotyped samples only.

The explicit subset step was necessary because GEMMA did not drop the placeholder `-9` values in the same way PLINK did; without this correction, the mixed-model run would have analyzed all 1135 samples rather than the intended 970 phenotyped individuals.

### 3.4 Association Methods

Three methods were benchmarked.

#### 3.4.1 Naive Linear Regression

Naive association testing was performed with PLINK2 additive linear regression and no structure covariates:

- tool: `plink2`
- model: additive single-SNP regression
- phenotype: simulated trait or `FT16`

This serves as the uncorrected baseline.

#### 3.4.2 Linear Regression with Principal Components

The second analysis used PLINK2 with the top 5 genotype PCs as fixed covariates:

- tool: `plink2`
- covariates: 5 genotype PCs
- model: additive regression with PC adjustment

This method aims to reduce broad ancestry or structure confounding while remaining computationally simple.

#### 3.4.3 Linear Mixed Model

The third analysis used GEMMA:

- tool: `gemma`
- kinship matrix: centered relatedness (`-gk 1`)
- association test: `-lmm 4`
- covariates: 5 genotype PCs

For the simulated phenotype, the LMM used the full 1135-sample genotype set. For the real `FT16` phenotype, the genotype data were first subset to the 970 phenotyped samples and then passed to GEMMA.

### 3.5 Evaluation Metrics

Method outputs were compared with:

- QQ plots
- Manhattan plots
- genomic inflation factor (`lambda_GC`)
- overlap among the top 100 associated SNPs

Evaluation output files:

- simulated trait:
  - `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/benchmark.qq.png`
  - `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/benchmark.manhattan.png`
  - `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/benchmark.lambda_gc.tsv`
  - `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/benchmark.top_overlap.tsv`
- FT16:
  - `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/ft16_benchmark.qq.png`
  - `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/ft16_benchmark.manhattan.png`
  - `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/ft16_benchmark.lambda_gc.tsv`
  - `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/ft16_benchmark.top_overlap.tsv`

## 4. Results

### 4.1 Simulated Phenotype Benchmark

Sample and variant counts:

- LR tested variants: 461
- LR+PCs tested variants: 461
- LMM tested variants: 461

Genomic inflation:

- `LR`: 19.2501
- `LR+PCs`: 11.6591
- `LMM`: 0.9209

Top-100 hit overlap:

- `LMM` vs `LR`: 30
- `LMM` vs `LR+PCs`: 49
- `LR` vs `LR+PCs`: 32

GEMMA variance explained estimate:

- `pve estimate = 0.509036`
- `se(pve) = 0.0527948`

Interpretation:

The simulated trait benchmark showed extreme inflation under naive regression, substantial but incomplete improvement with PC adjustment, and near-ideal calibration under the mixed model. The mixed-model `pve` estimate closely matched the simulation setting of `h2 = 0.5`, which supports the internal consistency of the simulation and analysis pipeline.

### 4.2 Real FT16 Phenotype Benchmark

Sample and variant counts:

- PLINK analyses used 970 phenotyped samples
- LR tested variants: 461
- LR+PCs tested variants: 461
- LMM tested variants: 412

The lower variant count in GEMMA reflects additional internal filtering on the phenotyped subset.

Genomic inflation:

- `LR`: 10.2228
- `LR+PCs`: 5.2349
- `LMM`: 0.9089

Top-100 hit overlap:

- `LMM` vs `LR`: 26
- `LMM` vs `LR+PCs`: 36
- `LR` vs `LR+PCs`: 41

GEMMA variance explained estimate:

- `pve estimate = 0.645208`
- `se(pve) = 0.0490408`

Interpretation:

The real phenotype produced the same qualitative pattern as the simulated benchmark. Inflation remained large in naive regression, PC adjustment reduced but did not remove inflation, and the mixed model remained near the desired calibration range. The real phenotype showed less inflation than the simulated trait in the regression baselines, but those baselines still remained far from well calibrated.

### 4.3 Direct Comparison Across Phenotype Settings

The main calibration metrics are summarized below.

| Phenotype | LR lambda_GC | LR+PCs lambda_GC | LMM lambda_GC |
| --- | ---: | ---: | ---: |
| Simulated | 19.2501 | 11.6591 | 0.9209 |
| FT16 | 10.2228 | 5.2349 | 0.9089 |

Top-hit agreement:

| Phenotype | LMM vs LR | LMM vs LR+PCs | LR vs LR+PCs |
| --- | ---: | ---: | ---: |
| Simulated | 30 | 49 | 32 |
| FT16 | 26 | 36 | 41 |

These results show that the precise degree of inflation depends on the phenotype, but the ranking of methods is stable across both settings:

`LMM` performs best, `LR+PCs` is intermediate, and `LR` performs worst.

## 5. Discussion

### 5.1 Main Scientific Conclusion

Across both the simulated phenotype and the real FT16 phenotype, the benchmark supports the same substantive conclusion: population structure correction matters substantially in this accession panel, and a kinship-based mixed model is much better behaved than either naive regression or regression with PCs alone.

The consistency of this result across two phenotype settings is important. The simulated benchmark confirms the expected statistical behavior in a controlled setting, while the real phenotype benchmark shows that the same calibration pattern persists in practice.

### 5.2 Interpretation of PC Adjustment

The PC-adjusted regression clearly improved over naive regression in both runs. However, the residual inflation was still substantial:

- simulated: 11.66
- FT16: 5.23

This indicates that broad principal components capture some structure, but not enough to fully model relatedness and fine-scale confounding in this dataset. In other words, PC adjustment helps, but it is not an adequate substitute for a mixed model here.

### 5.3 Interpretation of the Mixed Model

In both runs, the mixed model achieved `lambda_GC` values near 1:

- simulated: 0.9209
- FT16: 0.9089

This is the strongest evidence in the project that the LMM is the most statistically defensible method among the three. The mixed model also produced top-hit rankings that differed materially from the regression baselines, which suggests that method choice affects biological conclusions and is not merely changing p-value scale.

### 5.4 Effect of Phenotype Setting

The real phenotype run was less inflated overall than the simulated run in the two regression baselines. Several factors could explain this:

- different signal architecture,
- different sample count due to phenotype missingness,
- different effective variant subset in the LMM branch,
- phenotype-specific coupling to population structure.

Even so, the relative ordering of methods was stable, which increases confidence in the general conclusion.

## 6. Limitations

This project has several important limitations.

### 6.1 Corrupted Original Genotype Download

The original local VCF was not fully valid and had to be salvaged by truncating at the end of the parseable region. The benchmark therefore depends on a cleaned local derivative rather than a pristine upstream file.

### 6.2 Chromosome-Restricted Analysis

Only chromosome 4 was analyzed. This is acceptable for a focused benchmark, but it limits the generality of variant counts, kinship estimation, and inferred association behavior.

### 6.3 Small Post-QC Marker Set

Only 461 variants remained after QC. This is a small marker panel for both GWAS and kinship estimation, and it likely contributes to the repeated GEMMA warnings about many near-zero eigenvalues.

### 6.4 Incomplete Metadata Mapping

The local `accessions.csv` file contained only a small subset of accession mappings. As a result, the FT16 rerun used only direct `accession_id` overlap rather than a broader Assembly_ID remapping strategy. Fortunately, direct overlap still yielded 970 usable samples.

### 6.5 Single Simulation Replicate

The simulated benchmark used a single replicate. That is sufficient for a first-pass project report, but not for estimating variability across repeated phenotype realizations.

## 7. Recommended Next Steps

The most valuable project extensions would be:

1. obtain a clean full mapping file for Assembly_ID to Accession_ID,
2. rerun FT16 after full phenotype harmonization to maximize sample count,
3. assess sensitivity to QC thresholds,
4. run multiple simulated replicates,
5. extend the benchmark beyond chromosome 4,
6. compare runtime and memory across methods explicitly.

## 8. Conclusion

This project now provides a working, end-to-end GWAS benchmark in two complementary settings: a controlled simulated phenotype and a real FT16 phenotype. Both analyses lead to the same main result. Naive regression is strongly inflated, PC-adjusted regression improves but remains miscalibrated, and the mixed model is consistently the best-calibrated and most defensible approach on this structured *Arabidopsis* genotype panel.

That consistency across both phenotype settings is the strongest outcome of the project. It turns the repository from a scaffold into a completed comparative benchmark with a clear statistical conclusion.
