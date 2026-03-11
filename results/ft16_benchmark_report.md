# FT16 Benchmark Report

## Summary

This report documents the rerun of the chromosome 4 GWAS benchmark using the real phenotype table in `data/values.csv` instead of the earlier simulated trait. The phenotype in `values.csv` is `FT16`, and it was aligned directly to genotype sample IDs using the `accession_id` column.

## Phenotype Alignment

Input phenotype file:

- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/data/values.csv`

Prepared phenotype file:

- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/data/processed/FT16_pheno.tsv`

Alignment details:

- Genotype samples in QC set: 1135
- Samples with matched FT16 values: 970
- Samples without matched FT16 values: 165
- Duplicate accession IDs averaged: 1 (`8424`)

The local `accessions.csv` file was too incomplete to perform a broader Assembly_ID to Accession_ID remapping, so this rerun used direct `accession_id` matches only. That still yielded a large usable cohort of 970 samples.

## Genotype Processing

The genotype input and QC procedure were unchanged from the earlier run:

- chromosome 4 VCF cleaned from the valid prefix of the downloaded file
- PLINK2 conversion
- MAF >= 0.05
- missingness <= 0.05
- biallelic restriction

Final post-QC genotype set:

- 1135 individuals in the full genotype matrix
- 461 variants after QC

For LMM analysis, the genotype matrix was explicitly subset to the 970 phenotyped individuals before kinship estimation and association testing.

## Methods

Three methods were run:

1. Naive linear regression with PLINK2
2. Linear regression with 5 genotype PCs as fixed covariates
3. GEMMA linear mixed model with kinship and 5 PCs

Outputs:

- LR: `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/lr/lr.FT16.glm.linear`
- LR+PCs: `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/lr_pcs/lr_pcs.FT16.glm.linear`
- LMM: `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/lmm/lmm.assoc.txt`

Evaluation outputs:

- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/ft16_benchmark.qq.png`
- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/ft16_benchmark.manhattan.png`
- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/ft16_benchmark.lambda_gc.tsv`
- `/tscc/projects/ps-gaultonlab/junxif/CSE_284/results/plots/ft16_benchmark.top_overlap.tsv`

## Observed Results

### Sample and Variant Counts

- LR tested variants: 461
- LR+PCs tested variants: 461
- LMM tested variants: 412

GEMMA analyzed 412 variants on the phenotyped subset, which is expected because its internal filtering on the subset can differ slightly from PLINK2's plain association pass.

### Genomic Inflation

Observed `lambda_GC`:

- `LR`: 10.2228
- `LR+PCs`: 5.2349
- `LMM`: 0.9089

Interpretation:

- Naive regression remains strongly inflated.
- PC adjustment reduces inflation substantially, but the result is still clearly miscalibrated.
- The mixed model remains close to the ideal range around 1 and is the best-calibrated method in this rerun as well.

### Top-Hit Overlap

Top-100 SNP overlap:

- `LMM` vs `LR`: 26
- `LMM` vs `LR+PCs`: 36
- `LR` vs `LR+PCs`: 41

This indicates that model choice still changes the ranking of leading SNPs meaningfully under the real FT16 phenotype.

### LMM Variance Explained

GEMMA reported:

- `pve estimate = 0.645208`
- `se(pve) = 0.0490408`

This suggests a substantial fraction of FT16 variance is captured by the chromosome 4 genotype subset plus the kinship structure induced from it, though this should be interpreted cautiously because only one chromosome and a small post-QC marker set were analyzed.

## Comparison to the Earlier Simulated Run

Relative to the earlier simulated-trait benchmark:

- inflation decreased in the LR baselines
- LMM remained the best-calibrated method
- the conclusion did not change: mixed-model correction is still much more reliable than naive or PC-only regression

Updated lambda comparison:

- Simulated trait: `LR 19.25`, `LR+PCs 11.66`, `LMM 0.92`
- FT16 trait: `LR 10.22`, `LR+PCs 5.23`, `LMM 0.91`

The real phenotype appears less extreme than the simulated trait in terms of inflation under the linear-regression baselines, but those baselines are still far from acceptable calibration.

## Conclusion

The end-to-end rerun with the real FT16 phenotype completed successfully. The main conclusion is unchanged from the simulated benchmark: population structure correction matters substantially, principal components help but do not fully solve the problem, and the kinship-based mixed model is the most defensible of the three methods on this dataset.
