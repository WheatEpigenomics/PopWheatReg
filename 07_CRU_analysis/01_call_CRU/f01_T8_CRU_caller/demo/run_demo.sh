#!/usr/bin/env bash
# Toy run: 20-Mb one-chromosome genome, 508 ACRs, 1,467 agLinks, no R needed.
#
#   bash demo/run_demo.sh [outdir]
#
# Runs:
#   1) 01_run_cru.py   : build CRUs on the toy chromosome
#   2) 10_score_cru.py : score them with the frozen T7-A anchors
#
# Output defaults to $TMPDIR/cru_demo because the caller refuses to write
# inside its own source directory.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="${PYTHON:-python3}"
OUT="${1:-${TMPDIR:-/tmp}/cru_demo}"

cd "$ROOT"
"$PY" "$HERE/make_demo_data.py"
rm -rf "$OUT"

echo "== calling =="
"$PY" -B ./01_run_cru.py \
  --aglink "$HERE/data/aglink.tsv" \
  --gene-bed "$HERE/data/genes.bed" \
  --acr-bed "$HERE/data/acrs.bed" \
  --genome-sizes "$HERE/data/genome.sizes" \
  --config ./f01_rules/CRU_rules.json \
  --chroms 1A --workers 1 --control none \
  --mode all --plot none --block none \
  --outdir "$OUT/call"

echo
echo "== scoring =="
"$PY" -B ./10_score_cru.py \
  --cru-bed "$OUT/call/f04_merged/CRU.bed" \
  --aglink "$HERE/data/aglink.tsv" \
  --acr-bed "$HERE/data/acrs.bed" \
  --gene-bed "$HERE/data/genes.bed" \
  --genome-sizes "$HERE/data/genome.sizes" \
  --coacr "$HERE/data/coacr.tsv" \
  --outdir "$OUT/score"

echo
echo "CRUs   : $OUT/call/f04_merged/CRU.bed"
echo "scores : $OUT/score/04_CRU_scores_simple.tsv"
