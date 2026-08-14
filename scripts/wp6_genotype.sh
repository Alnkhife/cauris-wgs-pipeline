#!/usr/bin/env bash
# WP6: genotype context samples at the study SNP sites and build combined phylogeny.
set -o pipefail

ENV=$CONDA_PREFIX/bin
ROOT=$HOME/Projects/cauris-wgs-analysis
REF=$ROOT/data/raw/reference/B8441_refseq.fna
BAM=$ROOT/analysis/07_context/bam
SITES=$ROOT/analysis/07_context/sites.tsv
GENO=$ROOT/analysis/07_context/genotypes
OUT=$ROOT/analysis/07_context
mkdir -p "$GENO"

# per-context-sample: mpileup at the study sites -> haploid genotype VCF
for bam in "$BAM"/*.bam; do
  base=$(basename "$bam" .bam)
  OUTV=$GENO/$base.vcf
  [ -f "$OUTV" ] && continue
  echo "genotyping $base"
  "$ENV/bcftools" mpileup -f "$REF" -Q 20 -q 20 -R "$SITES" -d 500 "$bam" 2>/dev/null \
    | "$ENV/bcftools" call -m --ploidy 1 -v -o "$OUTV" 2>/dev/null
done
echo "genotyping done: $(ls "$GENO"/*.vcf 2>/dev/null | wc -l) samples"

# combine into the study alignment
"$ENV/python" "$ROOT/scripts/wp6_build_combined.py" "$OUT"
