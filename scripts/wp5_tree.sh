#!/usr/bin/env bash
# WP5: rebuilt maximum-likelihood phylogeny (IQ-TREE2) + SNP-distance matrix
# Runs on both the MycoSNP core-SNP alignment (if available) and the freebayes SNP set.
set -euo pipefail

ENV=$CONDA_PREFIX/bin
ROOT=$HOME/Projects/cauris-wgs-analysis
ALN_DIR=$ROOT/analysis/05_phylogeny
mkdir -p "$ALN_DIR"

# --- 1. freebayes-derived SNP alignment + snp-dists ---
$ENV/python $ROOT/scripts/wp5_build_aln.py \
    $ROOT/analysis/04_freebayes/freebayes_filtered.vcf \
    $ROOT/data/raw/reference/B8441_refseq.fna \
    "$ALN_DIR/freebayes_snps.fasta"

$ENV/snp-dists "$ALN_DIR/freebayes_snps.fasta" > "$ALN_DIR/freebayes_snpdists.tsv"
echo "freebayes snp-dists matrix written"

# --- 2. IQ-TREE2 on freebayes SNP alignment (ModelFinder + 1000 UFBoot) ---
$ENV/iqtree2 -s "$ALN_DIR/freebayes_snps.fasta" \
    -m MFP -bb 1000 -nt AUTO --prefix "$ALN_DIR/tree_freebayes" > "$ALN_DIR/iqtree_freebayes.log" 2>&1
echo "freebayes tree done"

# --- 3. MycoSNP alignment (from Drive results) if present ---
MSNP="$ROOT/data/raw/mycosnp_combined/phylogeny"
if [ -d "$MSNP" ]; then
    FA=$(find "$MSNP" -name "*.sth" -o -name "*.fasta" -o -name "*.fa" | head -1)
    if [ -n "$FA" ]; then
        cp "$FA" "$ALN_DIR/mycosnp_alignment.sth"
        $ENV/iqtree2 -s "$ALN_DIR/mycosnp_alignment.sth" \
            -m MFP -bb 1000 -nt AUTO --prefix "$ALN_DIR/tree_mycosnp" > "$ALN_DIR/iqtree_mycosnp.log" 2>&1
        echo "mycosnp tree done (alignment: $FA)"
    else
        echo "no mycosnp alignment fasta found in $MSNP"
    fi
fi
echo "WP5 done"
