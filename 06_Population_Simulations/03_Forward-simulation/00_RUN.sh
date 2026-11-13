#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-demo}"
fail() { echo "Error: $*" >&2; exit 2; }
[[ $# -le 2 ]] || fail 'Usage: bash 00_RUN.sh [demo|project|full] [OUTPUT_DIR]'
case "$MODE" in
  demo) OUT="${2:-$HERE/output_demo}"; EXTRA=(--demo) ;;
  project|full) [[ -n "${RAW_DIR:-}" ]] || fail 'Set RAW_DIR to the complete raw input directory, or run bash 00_RUN.sh demo.'; OUT="${2:-$HERE/output}"; EXTRA=() ;;
  -h|--help) echo 'Usage: bash 00_RUN.sh [demo|project|full] [OUTPUT_DIR]; full requires RAW_DIR.'; exit 0 ;;
  *) fail 'Unknown mode. Use demo, project, or full.' ;;
esac
if [[ "$MODE" != demo ]]; then
  for file in 05_acr_with_name_fil110.bed 06_aglink_sig_only.tsv haplotype_block.bed Ta_CS_v2.1.gff3; do
    [[ -f "$RAW_DIR/$file" ]] || fail "Missing $RAW_DIR/$file. Set RAW_DIR to the complete raw input directory."
  done
  [[ -d "$RAW_DIR/vcf_file" ]] || fail "Missing $RAW_DIR/vcf_file. Set RAW_DIR to the complete raw input directory."
fi
OUT="$(realpath -m -- "$OUT")"
PYTHON="${PYTHON_BIN:-python3}"
"$PYTHON" -c 'import numpy, matplotlib' >/dev/null 2>&1 || fail 'Python with NumPy and Matplotlib is required. Install these packages or set PYTHON_BIN.'
export MPLBACKEND=Agg OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
"$PYTHON" "$HERE/01_forward_simulation.py" --config "$HERE/config.json" --output "$OUT" "${EXTRA[@]}" full
"$PYTHON" "$HERE/02_plot_results.py" --config "$OUT/config.json" --input "$OUT" --output "$OUT"
