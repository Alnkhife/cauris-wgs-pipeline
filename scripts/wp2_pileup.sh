#!/usr/bin/env bash
# WP2 stage B: independent read-level variant confirmation at resistance loci.
# bcftools mpileup across each gene's CDS (+100 bp flanks), call variants,
# then annotate amino-acid changes with python (stage C in wp2_annotate.py).
set -euo pipefail

ENV=$CONDA_PREFIX/bin
ROOT=$HOME/Projects/cauris-wgs-analysis
REF=$ROOT/data/raw/reference/B8441_refseq.fna
BAMDIR=$ROOT/analysis/01_alignment/bam
OUT=$ROOT/analysis/02_variant_confirm
mkdir -p "$OUT/vcfs"

# gene -> chrom:start-end (CDS span)
REGIONS=(
  "ERG11:NC_140807.1:1474622-1476396"
  "CDR1:NC_140808.1:348828-353554"
  "FCY1:NC_140806.1:1098491-1099191"
  "FKS1:NC_140808.1:2088086-2093952"
)

BAMS=()
for s in Sample01 Sample03 Sample04 Sample05 Sample06 Sample07 Sample08 Sample09 Sample10 Sample11 Sample12 Sample13 Sample14 Sample15 Sample16 Sample17 Sample18 Sample19 Sample20; do
  BAMS+=("$BAMDIR/$s.sorted.md.bam")
done

for entry in "${REGIONS[@]}"; do
  gene="${entry%%:*}"
  region="${entry#*:}"
  chrom="${region%%:*}"
  span="${region#*:}"
  "$ENV/bcftools" mpileup -f "$REF" -r "$chrom:$span" -Q 20 -q 20 --annotate FORMAT/AD,FORMAT/DP -d 500 \
      "${BAMS[@]}" 2>"$OUT/vcfs/$gene.mpileup.log" \
    | "$ENV/bcftools" call -m --ploidy 1 -v -o "$OUT/vcfs/$gene.raw.vcf" 2>>"$OUT/vcfs/$gene.mpileup.log" || true
  "$ENV/bcftools" view -h "$OUT/vcfs/$gene.raw.vcf" > /dev/null 2>&1 || echo "WARN: $gene vcf empty"
  echo "$gene: $(grep -vc '^#' "$OUT/vcfs/$gene.raw.vcf" 2>/dev/null || echo 0) variant records"
done
echo "WP2 pileup done -> $OUT/vcfs/"
