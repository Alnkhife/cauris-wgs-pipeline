#!/usr/bin/env bash
# iTOL batch upload + export using the account API key.
set -o pipefail

ROOT=$HOME/Projects/cauris-wgs-analysis
OUT=$ROOT/output
KEY=xgFHdqQjVFuF6l0AgE0ZFg
ITOL=https://itol.embl.de

# --- clade colour assignments (matches thesis figures) ---
CLADE1="Sample01 Sample04 Sample05 Sample14 Sample15 Sample16 Sample19 Sample20"
CLADE2="Sample06 Sample07 Sample09 Sample10 Sample13"
CLADE3="Sample03 Sample08 Sample11 Sample12 Sample17 Sample18"
CONTEXT="ERR14122511 ERR14122512 ERR14122513 ERR14122514 ERR14122515 ERR14122516 ERR14122517 ERR14122518 ERR14122519 ERR14122520 ERR14122521 ERR14122522 ERR14122523 ERR14122524 ERR14122525 ERR14122526 ERR14122527 ERR14122528 ERR14122529 ERR14122530 ERR14122531 ERR14122532 ERR14122533 SRR10461149 SRR10461151 SRR10461152 SRR10461153"

make_colors() { # make_colors <file> [with-context]
  local f=$1; shift
  { echo "TREE_COLORS"; echo "SEPARATOR TAB"; echo "DATA"; }
  for s in $CLADE1; do echo -e "$s\tlabel\t#1a237e"; echo -e "$s\tbranch\t#1a237e"; done
  for s in $CLADE2; do echo -e "$s\tlabel\t#00897b"; echo -e "$s\tbranch\t#00897b"; done
  for s in $CLADE3; do echo -e "$s\tlabel\t#d62728"; echo -e "$s\tbranch\t#d62728"; done
  if [ "$1" = "with-context" ]; then
    for s in $CONTEXT; do echo -e "$s\tlabel\t#9e9e9e"; echo -e "$s\tbranch\t#d9d9d9"; done
  fi
} > "$f"

make_colors /tmp/colors19.txt
make_colors /tmp/colors46.txt with-context

for spec in "tree19|analysis/05_phylogeny/tree_freebayes.contree|colors19.txt|Figure_4_phylogeny_iTOL" \
            "tree46|analysis/07_context/tree_combined.contree|colors46.txt|Figure_WP6_combined_tree_iTOL"; do
  name=${spec%%|*}; rest=${spec#*|}
  tree=${rest%%|*}; rest=${rest#*|}
  ann=${rest%%|*}; outbase=${rest#*|}
  echo "=== $name ==="
  RESP=$(curl -s -F "APIkey=$KEY" -F "treeFile=@$ROOT/$tree" -F "annotationFile=@/tmp/$ann" "$ITOL/batch_upload.cgi")
  echo "$RESP"
  TID=$(echo "$RESP" | grep -oE '^SUCCESS: Tree [0-9]+' | tail -1 | grep -oE '[0-9]+')
  if [ -n "$TID" ]; then
    curl -s "$ITOL/export.cgi?tree=$TID&format=svg" -o "$OUT/$outbase.svg"
    curl -s "$ITOL/export.cgi?tree=$TID&format=png" -o "$OUT/$outbase.png"
    for fmt in svg png; do
      f="$OUT/$outbase.$fmt"
      sz=$(stat -f%z "$f" 2>/dev/null || echo 0)
      echo "  $fmt: $sz bytes"
      if [ "$sz" -lt 5000 ]; then head -c 200 "$f"; echo; fi
    done
  else
    echo "  !! no tree ID from batch upload"
  fi
done
