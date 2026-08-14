#!/usr/bin/env bash
# Per-gene mean coverage for all samples using samtools bedcov (mosdepth was broken).
set -euo pipefail
ENV=$CONDA_PREFIX/bin
ROOT=$HOME/Projects/cauris-wgs-analysis
GENES=$ROOT/data/raw/reference/genes.bed
BAMDIR=$ROOT/analysis/01_alignment/bam
OUT=$ROOT/analysis/01_alignment/stats
mkdir -p "$OUT"

for s in Sample01 Sample03 Sample04 Sample05 Sample06 Sample07 Sample08 Sample09 Sample10 Sample11 Sample12 Sample13 Sample14 Sample15 Sample16 Sample17 Sample18 Sample19 Sample20; do
  OUTB=$OUT/$s.per-gene.tsv
  [ -f "$OUTB" ] && continue
  echo -e "chrom\tstart\tend\tgene\ttotal_depth\tlength\tmean_cov" > "$OUTB"
  "$ENV/samtools" bedcov "$GENES" "$BAMDIR/$s.sorted.md.bam" 2>/dev/null \
    | awk -F'\t' '{len=$3-$2; printf "%s\t%s\t%s\t%s\t%s\t%s\t%.2f\n", $1,$2,$3,$4,$5,len,$5/len}' \
    >> "$OUTB"
  echo "$s done ($(wc -l < "$OUTB") genes)"
done
echo "per-gene coverage complete"
