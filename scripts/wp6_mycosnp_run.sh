#!/usr/bin/env bash
# Run the official CDC MycoSNP-NF pipeline locally (docker) on the 19 C. auris samples.
# Produces the reference MycoSNP/GATK call set for the WP4 cross-validation.
set -euo pipefail

ROOT=$HOME/Projects/cauris-wgs-analysis
PIPE=$ROOT/analysis/mycosnp-nf
OUT=$ROOT/analysis/06_mycosnp
mkdir -p "$OUT"

export NXF_OPTS="-Xms1g -Xmx6g"

cd "$PIPE"
nextflow run . \
    -profile docker \
    --input "$ROOT/data/raw/samplesheet.csv" \
    --fasta "$ROOT/data/raw/reference/B8441_refseq.fna" \
    --coverage 50 \
    --outdir "$OUT" \
    -with-report "$ROOT/logs/mycosnp_report.html" \
    -with-trace "$ROOT/logs/mycosnp_trace.txt" \
    -resume
echo "MycoSNP-NF run finished"
