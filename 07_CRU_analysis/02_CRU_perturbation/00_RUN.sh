#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-demo}"
fail() { echo "Error: $*" >&2; exit 2; }
[[ $# -le 2 ]] || fail 'Usage: bash 00_RUN.sh [demo|project|full] [OUTPUT_DIR]'
case "$MODE" in
  demo) SOURCE="$HERE/demo"; OUT="${2:-$HERE/output_demo}" ;;
  project) SOURCE="${PROJECT_ROOT:-}"; [[ -n "$SOURCE" ]] || fail 'Set PROJECT_ROOT to a prepared CRU project containing data/ and docs/, or run bash 00_RUN.sh demo.'; OUT="${2:-$HERE/output}" ;;
  full) RAW="${INPUT_ROOT:-}"; [[ -n "$RAW" ]] || fail 'Set INPUT_ROOT to the raw 01_data directory, or run bash 00_RUN.sh demo.'; OUT="${2:-$HERE/output_full}" ;;
  -h|--help) echo 'Usage: bash 00_RUN.sh [demo|project|full] [OUTPUT_DIR]; project requires PROJECT_ROOT, full requires INPUT_ROOT.'; exit 0 ;;
  *) fail 'Unknown mode. Use demo, project, or full.' ;;
esac
if [[ "$MODE" == full ]]; then
  for file in 14_CRU/csv2_T8_CRU.bed 04_acr/05_acr_with_name_fil110.bed 04_acr/07_200_spike_ACR_fil110_signal.txt 03_exp/05_202hex_csV2_rna_spike_ave_exp_sum1_clean.txt 01_ref/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3; do
    [[ -f "$RAW/$file" ]] || fail "Missing $RAW/$file. Set INPUT_ROOT to the raw 01_data directory."
  done
  SOURCE="$RAW"
else
  for file in data/accessions.tsv data/sentinel_selection.tsv data/sentinel_accessibility.tsv.gz data/distance_matched_pairs.tsv.gz data/feature_effects.tsv.gz data/feature_effects_high.tsv.gz data/sentinel_high_groups.tsv docs/analysis_parameters.tsv docs/high_contrast_parameters.tsv docs/raw_input_manifest.tsv; do
    [[ -f "$SOURCE/$file" ]] || fail "Missing $SOURCE/$file. Set PROJECT_ROOT to the prepared CRU project."
  done
  BED="$SOURCE/data/cru_input.bed"
  if [[ ! -f "$BED" && -n "${INPUT_ROOT:-}" ]]; then BED="$INPUT_ROOT/14_CRU/csv2_T8_CRU.bed"; fi
  [[ -f "$BED" ]] || fail 'Missing data/cru_input.bed. Add the CRU BED to PROJECT_ROOT/data/ or set INPUT_ROOT to the raw data directory.'
fi
OUT="$(realpath -m -- "$OUT")"; SOURCE="$(realpath -m -- "$SOURCE")"
[[ "$OUT" != "$SOURCE" && "$SOURCE" != "$OUT/"* && "$OUT" != "$SOURCE/data" && "$OUT" != "$SOURCE/docs" ]] || fail 'Output must be separate from the input data directory.'
RSCRIPT="${RSCRIPT_BIN:-Rscript}"
command -v "$RSCRIPT" >/dev/null || fail 'Rscript not found. Install R or set RSCRIPT_BIN.'
if [[ "$MODE" == full ]]; then
  "${PYTHON_BIN:-python3}" "$HERE/01_prepare_cru.py" --input-root "$RAW" --outdir "$OUT"
else
  mkdir -p -- "$OUT/data" "$OUT/docs"
  cp -- "$SOURCE"/data/* "$OUT/data/"
  cp -- "$BED" "$OUT/data/cru_input.bed"
  cp -- "$SOURCE"/docs/{analysis_parameters,high_contrast_parameters,raw_input_manifest}.tsv "$OUT/docs/"
fi
export XDG_CACHE_HOME="$OUT/.cache"
mkdir -p -- "$XDG_CACHE_HOME"
"$RSCRIPT" "$HERE/02_analyze_cru.R" all "$OUT"
echo "Output: $OUT"
