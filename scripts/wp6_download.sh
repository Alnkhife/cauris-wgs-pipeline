#!/usr/bin/env bash
# WP6: download Saudi C. auris context genomes (23 ENA + 4 NCBI SRR).
set -o pipefail

ROOT=$HOME/Projects/cauris-wgs-analysis
OUT=$ROOT/data/raw/context_fastq
LOG=$ROOT/logs/wp6_download.log
mkdir -p "$OUT"
: > "$LOG"

log() { echo "[$(date +%H:%M:%S)] $*" | tee -a "$LOG"; }

# --- ENA runs: (accession) ---
ENA=(ERR14122511 ERR14122512 ERR14122513 ERR14122514 ERR14122515 ERR14122516 ERR14122517 ERR14122518 ERR14122519 ERR14122520 ERR14122521 ERR14122522 ERR14122523 ERR14122524 ERR14122525 ERR14122526 ERR14122527 ERR14122528 ERR14122529 ERR14122530 ERR14122531 ERR14122532 ERR14122533)

while IFS=$'\t' read -r acc ftp _; do
  [ -z "$acc" ] && continue
  IFS=';' read -r -a files <<< "$ftp"
  for f in "${files[@]}"; do
    name=$(basename "$f")
    dest="$OUT/$name"
    if [ -f "$dest" ] && gzip -t "$dest" 2>/dev/null; then
      log "SKIP $name"
      continue
    fi
    log "DL $name"
    curl -sL --retry 3 -o "$dest" "https://$f" || { log "FAIL $name"; rm -f "$dest"; }
    gzip -t "$dest" 2>/dev/null && log "OK  $name" || { log "CORRUPT $name"; rm -f "$dest"; }
  done
done < "$ROOT/data/raw/ena_manifest.tsv"

# --- NCBI SRR runs via sra-tools ---
SRR=(SRR10461149 SRR10461151 SRR10461152 SRR10461153)
for acc in "${SRR[@]}"; do
  r1="$OUT/${acc}_1.fastq.gz"; r2="$OUT/${acc}_2.fastq.gz"
  if [ -f "$r1" ] && [ -f "$r2" ] && gzip -t "$r1" 2>/dev/null && gzip -t "$r2" 2>/dev/null; then
    log "SKIP $acc"
    continue
  fi
  log "SRA $acc (fasterq-dump)"
  mkdir -p /tmp/sra_$acc && cd /tmp/sra_$acc
  docker run --rm -v "$PWD":/data quay.io/biocontainers/sra-tools:2.11.0--pl5262h314213e_0 \
    fasterq-dump "$acc" --split-files --outdir /data -e 4 2>>"$LOG" || { log "FAIL $acc"; continue; }
  gzip -c /tmp/sra_$acc/${acc}_1.fastq > "$r1" 2>/dev/null
  gzip -c /tmp/sra_$acc/${acc}_2.fastq > "$r2" 2>/dev/null
  rm -rf /tmp/sra_$acc
  gzip -t "$r1" 2>/dev/null && gzip -t "$r2" 2>/dev/null && log "OK  $acc" || log "CORRUPT $acc"
done

log "WP6 download complete: $(ls "$OUT" | wc -l) files"
