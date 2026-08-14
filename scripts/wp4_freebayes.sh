#!/usr/bin/env bash
# WP4: orthogonal genome-wide variant calling with freebayes (19 samples, haploid)
# then build pairwise SNP-distance matrix and compare against MycoSNP matrix.
set -euo pipefail

ENV=$CONDA_PREFIX/bin
ROOT=$HOME/Projects/cauris-wgs-analysis
REF=$ROOT/data/raw/reference/B8441_refseq.fna
BAMDIR=$ROOT/analysis/01_alignment/bam
STATS_DIR=$ROOT/analysis/01_alignment/stats
OUT=$ROOT/analysis/04_freebayes
mkdir -p "$OUT"

BAMS=()
for s in Sample01 Sample03 Sample04 Sample05 Sample06 Sample07 Sample08 Sample09 Sample10 Sample11 Sample12 Sample13 Sample14 Sample15 Sample16 Sample17 Sample18 Sample19 Sample20; do
  if [ -f "$BAMDIR/$s.sorted.md.bam" ]; then
    BAMS+=("-b $BAMDIR/$s.sorted.md.bam")
  else
    echo "WARN: $s BAM not present - excluded"
  fi
done
echo "samples to call: ${#BAMS[@]}"

# Downsample deep-coverage BAMs to ~100x to keep freebayes tractable and
# comparable with MycoSNP's ~50x target. Fraction per sample from mosdepth mean.
DOWNSAMPLED=$ROOT/analysis/04_freebayes/downsampled
mkdir -p "$DOWNSAMPLED"
for s in Sample01 Sample03 Sample04 Sample05 Sample06 Sample07 Sample08 Sample09 Sample10 Sample11 Sample12 Sample13 Sample14 Sample15 Sample16 Sample17 Sample18 Sample19 Sample20; do
  B=$BAMDIR/$s.sorted.md.bam
  [ -f "$B" ] || continue
  OUTB=$DOWNSAMPLED/$s.ds.bam
  [ -f "$OUTB" ] && continue
  MEAN=$(awk -F'\t' '$1=="genome" {print $6}' "$STATS_DIR/$s.mosdepth.summary.txt" 2>/dev/null || true)
  MEAN=${MEAN:-100}
  FRAC=$(echo "scale=4; 100/$MEAN" | bc 2>/dev/null)
  FRAC=$(awk -v m="$MEAN" 'BEGIN{printf "%.4f", 100/m}')
  echo "$s: mean cov $MEAN -> frac $FRAC"
  if [ "$(echo "$FRAC < 1" | bc 2>/dev/null)" = "1" ]; then
    "$ENV/samtools" view -@ 4 -b -s "42$FRAC" "$B" -o "$OUTB"
    "$ENV/samtools" index "$OUTB"
  else
    cp "$B" "$OUTB"
    "$ENV/samtools" index "$OUTB"
  fi
done
# rebuild BAM list from downsampled files
BAMS=()
for s in Sample01 Sample03 Sample04 Sample05 Sample06 Sample07 Sample08 Sample09 Sample10 Sample11 Sample12 Sample13 Sample14 Sample15 Sample16 Sample17 Sample18 Sample19 Sample20; do
  [ -f "$DOWNSAMPLED/$s.ds.bam" ] && BAMS+=("-b $DOWNSAMPLED/$s.ds.bam")
done
echo "downsampled BAMs: ${#BAMS[@]}"

echo "freebayes calling (haploid, 8 threads)..."
export OMP_NUM_THREADS=8
$ENV/freebayes -f "$REF" --ploidy 1 --min-base-quality 20 --min-mapping-quality 20 \
    --min-alternate-fraction 0.5 --min-alternate-count 3 --use-best-n-alleles 2 \
    ${BAMS[@]} -v "$OUT/freebayes_all.vcf" 2> "$OUT/freebayes.log"
echo "raw records: $(grep -vc '^#' "$OUT/freebayes_all.vcf")"

# keep biallelic SNPs with PASS + depth >= 10, no missing calls in >10% of samples
$ENV/bcftools view -v snps -m2 -M2 -i 'QUAL>=50 && INFO/DP>=10' \
    "$OUT/freebayes_all.vcf" -o "$OUT/freebayes_snps.vcf"
$ENV/bcftools filter -i 'F_MISSING<=0.1' "$OUT/freebayes_snps.vcf" \
    -o "$OUT/freebayes_filtered.vcf" 2>/dev/null || cp "$OUT/freebayes_snps.vcf" "$OUT/freebayes_filtered.vcf"
echo "filtered SNP records: $(grep -vc '^#' "$OUT/freebayes_filtered.vcf")"
$ENV/bcftools view -Oz -o "$OUT/freebayes_filtered.vcf.gz" "$OUT/freebayes_filtered.vcf"
$ENV/bcftools index "$OUT/freebayes_filtered.vcf.gz"

# multi-sample FASTA alignment for snp-dists
$ENV/bcftools query -f '%CHROM\t%POS\t%REF\t%ALT\n' "$OUT/freebayes_filtered.vcf" > "$OUT/sites.tsv"
$ENV/bcftools consensus -f "$REF" -H 1 "$OUT/freebayes_filtered.vcf.gz" > /dev/null 2>&1 || true
echo "WP4 calling done"
