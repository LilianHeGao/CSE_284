#!/usr/bin/env bash
set -euo pipefail

CFG=${1:-config/project.env}
source "$CFG"

mkdir -p "$RAW_DIR" "$PROC_DIR" "$RESULTS_DIR"

tail -n +2 $PROC_DIR/pheno.txt > "$PROC_DIR/pheno.keep_ids.txt"

"$BCFTOOLS" view -v snps -m2 -M2 "$VCF_GZ" -Oz -o "$RAW_DIR/1001G.Chr4_clean.vcf.gz"

# convert VCF to pfile and filerinbg accession to only those with phenotypes
"$PLINK2" \
  --vcf "$VCF_GZ" \
  --keep "$PROC_DIR/pheno.keep_ids.txt" \
  --chr-set 5 \
  --allow-extra-chr \
  --set-all-var-ids @:#:\$r:\$a \
  --make-pgen \
  --out "$PROC_DIR/chr${CHR}"

# QC
"$PLINK2" \
  --pfile "$PROC_DIR/chr${CHR}" \
  --maf "$MAF" \
  --geno "$GENO_MISSING" \
  --make-bed \
  --out "$PROC_DIR/chr${CHR}.qc"
