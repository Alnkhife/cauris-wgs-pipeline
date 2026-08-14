#!/usr/bin/env bash
# WP6: align Saudi context genomes to B8441 and genotype them at the
# 1,301 study SNP sites. Then build the combined phylogeny.
set -o pipefail

ENV=$CONDA_PREFIX/bin
ROOT=$HOME/Projects/cauris-wgs-analysis
FASTQ=$ROOT/data/raw/context_fastq
REF=$ROOT/data/raw/reference/B8441_refseq.fna
SITES_VCF=$ROOT/analysis/04_freebayes/freebayes_filtered.vcf
BAM=$ROOT/analysis/07_context/bam
SITES=$ROOT/analysis/07_context/sites.tsv
mkdir -p "$BAM"

# site list: chrom<TAB>pos from the study SNP set
$ENV/bcftools query -f '%CHROM\t%POS\n' "$SITES_VCF" > "$SITES"
echo "sites: $(wc -l < "$SITES")"

for f in "$FASTQ"/*_1.fastq.gz; do
  base=$(basename "$f" _1.fastq.gz)
  r1="$f"; r2="$FASTQ/${base}_2.fastq.gz"
  [ -f "$r2" ] || { echo "MISSING $r2"; continue; }
  OUTB=$BAM/$base.bam
  if [ -f "$OUTB" ] && "$ENV/samtools" quickcheck "$OUTB" 2>/dev/null; then
    echo "SKIP $base"
    continue
  fi
  echo "align $base"
  "$ENV/bwa" mem -t 4 -R "@RG\tID:$base\tSM:$base\tPL:ILLUMINA" "$REF" "$r1" "$r2" 2>/dev/null \
    | "$ENV/samtools" sort -@ 4 -o "$OUTB" - && \
    "$ENV/samtools" index "$OUTB"
done
echo "alignment done: $(ls "$BAM"/*.bam 2>/dev/null | wc -l) BAMs"
