#!/usr/bin/env bash
# WP6 runner: waits for downloads, then chains align -> genotype -> combined tree.
set -o pipefail

ROOT=$HOME/Projects/cauris-wgs-analysis
FASTQ=$ROOT/data/raw/context_fastq
PROG=$ROOT/scripts/progress.sh
LOG=$ROOT/logs/wp6_run.log

notify() { "$PROG" notify "$1" "$2"; }
log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

# wait for all 50 files (5 min checks), notifying every 10
last=0
while true; do
  n=$(ls "$FASTQ"/*.fastq.gz 2>/dev/null | wc -l | tr -d ' ')
  if [ "$n" -ge 50 ]; then
    log "all 50 context FASTQs present"
    break
  fi
  if [ $((n / 10)) -gt $((last / 10)) ]; then
    notify "WP6 downloads: $n/50 files" "$(( $(du -sm "$FASTQ" | cut -f1) / 1024 )) GB of ~19 GB"
  fi
  last=$n
  sleep 300
done
"$PROG" stage wp6 "downloads done"
notify "WP6 downloads complete (50/50)" "Starting alignment of 27 Saudi genomes"

# 1. align
bash "$ROOT/scripts/wp6_align.sh" >> "$LOG" 2>&1
nb=$(ls "$ROOT/analysis/07_context/bam"/*.bam 2>/dev/null | wc -l | tr -d ' ')
"$PROG" stage wp6 "aligned $nb/27"
notify "WP6: $nb/27 aligned" "Starting site genotyping"

# 2. genotype at study sites
bash "$ROOT/scripts/wp6_genotype.sh" >> "$LOG" 2>&1
ng=$(ls "$ROOT/analysis/07_context/genotypes"/*.vcf 2>/dev/null | wc -l | tr -d ' ')
"$PROG" stage wp6 "genotyped $ng/27"
notify "WP6: $ng/27 genotyped" "Building combined phylogeny"

# 3. combined tree
bash "$ROOT/scripts/wp6_build_combined.py" >> "$LOG" 2>&1 || \
  "$ROOT/miniforge3/envs/cauris/bin/python" "$ROOT/scripts/wp6_build_combined.py" >> "$LOG" 2>&1
if [ -f "$ROOT/analysis/07_context/tree_combined.contree" ]; then
  "$PROG" done wp6
  notify "WP6 COMPLETE" "Combined 46-isolate tree ready (19 study + 27 Saudi)"
else
  notify "WP6 tree FAILED" "Check $LOG"
fi
