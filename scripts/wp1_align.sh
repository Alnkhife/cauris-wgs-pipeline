#!/usr/bin/env bash
# WP1: independent alignment of 19 C. auris samples to B8441 with bwa
# Produces per-sample: sorted, markdup BAM + mosdepth coverage + flagstat.
set -euo pipefail

ENV=$CONDA_PREFIX/bin
ROOT=$HOME/Projects/cauris-wgs-analysis
FASTQ=$ROOT/data/raw/fastq
REF=$ROOT/data/raw/reference
BAM=$ROOT/analysis/01_alignment/bam
STATS=$ROOT/analysis/01_alignment/stats
mkdir -p "$BAM" "$STATS"

THREADS=6
SAMPLES=(01 03 04 05 06 07 08 09 10 11 12 13 14 15 16 17 18 19 20)

log() { echo "[$(date +%H:%M:%S)] $*"; }

for s in "${SAMPLES[@]}"; do
    R1=$FASTQ/Candida_Auris_Sample${s}_R1.fastq.gz
    R2=$FASTQ/Candida_Auris_Sample${s}_R2.fastq.gz
    OUT=$BAM/Sample${s}
    if [ -f "$OUT.sorted.md.bam" ] && "$ENV/samtools" quickcheck "$OUT.sorted.md.bam" 2>/dev/null \
        && [ "$("$ENV/samtools" view -c "$OUT.sorted.md.bam" 2>/dev/null || echo 0)" -gt 100000 ]; then
        log "SKIP Sample${s} (done)"
        continue
    fi
    # remove any partial/corrupt outputs for this sample
    rm -f "$OUT.name.bam" "$OUT.name.fx.bam" "$OUT.sorted.bam" "$OUT.sorted.bam.bai" \
          "$OUT.sorted.md.bam" "$OUT.sorted.md.bam.bai" \
          "$OUT.sorted.fx.bam" "$OUT.sorted.fx.bam.bai" "$OUT.sorted.fx.sorted.bam"
    if [ ! -f "$R1" ] || [ ! -f "$R2" ]; then
        log "WAIT Sample${s} (fastq not downloaded yet)"
        continue
    fi
    log "align Sample${s}"
    "$ENV/bwa" mem -t "$THREADS" -R "@RG\tID:Sample${s}\tSM:Sample${s}\tPL:ILLUMINA" \
        "$REF/B8441_refseq.fna" "$R1" "$R2" 2>"$STATS/Sample${s}.bwa.log" \
    | "$ENV/samtools" sort -n -@ "$THREADS" -o "$OUT.name.bam" - && \
    "$ENV/samtools" fixmate -m -@ "$THREADS" "$OUT.name.bam" "$OUT.name.fx.bam" && \
    "$ENV/samtools" sort -@ "$THREADS" -o "$OUT.sorted.bam" "$OUT.name.fx.bam" && \
    "$ENV/samtools" index "$OUT.sorted.bam" && \
    "$ENV/samtools" markdup -r -@ "$THREADS" "$OUT.sorted.bam" "$OUT.sorted.md.bam" 2>"$STATS/Sample${s}.markdup.log" && \
    "$ENV/samtools" index "$OUT.sorted.md.bam" && \
    rm -f "$OUT.name.bam" "$OUT.name.fx.bam" "$OUT.sorted.bam" "$OUT.sorted.bam.bai"
    log "Sample${s} done"
done

# coverage + flagstat for all
for s in "${SAMPLES[@]}"; do
    OUT=$BAM/Sample${s}
    if [ ! -f "$STATS/Sample${s}.mosdepth.regions.bed.gz" ]; then
        "$ENV/mosdepth" -t 4 --by 1000 "$STATS/Sample${s}.mosdepth" "$OUT.sorted.md.bam" >/dev/null 2>&1 || true
    fi
    if [ ! -f "$STATS/Sample${s}.mosdepth.per-gene.bed.gz" ]; then
        "$ENV/mosdepth" -t 4 --by "$ROOT/data/raw/reference/genes.bed" \
            "$STATS/Sample${s}.mosdepth.per-gene" "$OUT.sorted.md.bam" >/dev/null 2>&1 || true
    fi
    if [ ! -f "$STATS/Sample${s}.flagstat.txt" ]; then
        "$ENV/samtools" flagstat "$OUT.sorted.md.bam" > "$STATS/Sample${s}.flagstat.txt"
    fi
done
log "WP1 alignment complete"
