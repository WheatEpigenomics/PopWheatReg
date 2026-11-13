#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-demo}"
fail() { echo "Error: $*" >&2; exit 2; }
[[ $# -le 2 ]] || fail 'Usage: bash 00_RUN.sh [demo|project|full] [OUTPUT_DIR]'
case "$MODE" in
  demo) SOURCE="$HERE/demo"; OUT="${2:-$HERE/output_demo}" ;;
  project|full) SOURCE="${PROJECT_ROOT:-}"; [[ -n "$SOURCE" ]] || fail 'Set PROJECT_ROOT to a prepared project containing data/ and docs/, or run bash 00_RUN.sh demo.'; OUT="${2:-$HERE/output}" ;;
  -h|--help) echo 'Usage: bash 00_RUN.sh [demo|project|full] [OUTPUT_DIR]; full requires PROJECT_ROOT.'; exit 0 ;;
  *) fail 'Unknown mode. Use demo, project, or full.' ;;
esac
for file in data/analysis_blocks.tsv.gz data/SNP_density_2Mb.tsv.gz docs/analysis_settings.tsv docs/data_dictionary.tsv docs/SNP_covariate_provenance.tsv; do
  [[ -f "$SOURCE/$file" ]] || fail "Missing $SOURCE/$file. Set PROJECT_ROOT to the prepared data directory."
done
OUT="$(realpath -m -- "$OUT")"; SOURCE="$(realpath -m -- "$SOURCE")"
[[ "$OUT" != "$SOURCE" && "$SOURCE" != "$OUT/"* && "$OUT" != "$SOURCE/data" && "$OUT" != "$SOURCE/docs" ]] || fail 'Output must be separate from the input project.'
RSCRIPT="${RSCRIPT_BIN:-Rscript}"
command -v "$RSCRIPT" >/dev/null || fail 'Rscript not found. Install R or set RSCRIPT_BIN.'
mkdir -p -- "$OUT/data" "$OUT/docs" "$OUT/.cache"
cp -- "$SOURCE"/data/{analysis_blocks,SNP_density_2Mb}.tsv.gz "$OUT/data/"
cp -- "$SOURCE"/docs/{analysis_settings,data_dictionary,SNP_covariate_provenance}.tsv "$OUT/docs/"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
export XDG_CACHE_HOME="$OUT/.cache"
"$RSCRIPT" "$HERE/01_fit_models.R" "$OUT"
"$RSCRIPT" "$HERE/02_plot_validate.R" all "$OUT"
echo "Output: $OUT"
