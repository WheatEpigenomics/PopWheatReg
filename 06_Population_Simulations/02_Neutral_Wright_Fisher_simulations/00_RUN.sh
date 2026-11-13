#!/usr/bin/env bash
set -euo pipefail
HERE="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
MODE="${1:-demo}"
fail() { echo "Error: $*" >&2; exit 2; }
[[ $# -le 2 ]] || fail 'Usage: bash 00_RUN.sh [demo|project|full] [OUTPUT_DIR]'
case "$MODE" in
  demo) SOURCE="$HERE/demo"; OUT="${2:-$HERE/output_demo}"; EXTRA=(--demo) ;;
  project|full) SOURCE="${INPUT_ROOT:-${SOURCE_RUN:-}}"; [[ -n "$SOURCE" ]] || fail 'Set INPUT_ROOT to the prepared neutral data directory (see settings.ini), or run bash 00_RUN.sh demo.'; OUT="${2:-$HERE/output}"; EXTRA=() ;;
  -h|--help) echo 'Usage: bash 00_RUN.sh [demo|project|full] [OUTPUT_DIR]; full requires INPUT_ROOT (SOURCE_RUN is also accepted).'; exit 0 ;;
  *) fail 'Unknown mode. Use demo, project, or full.' ;;
esac
[[ -d "$SOURCE" ]] || fail "Input directory not found: $SOURCE. Set INPUT_ROOT to your prepared data directory."
PYTHON="${PYTHON_BIN:-python3}"
"$PYTHON" -c 'import numpy' >/dev/null 2>&1 || fail 'Python with NumPy is required. Install numpy or set PYTHON_BIN.'
export OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
exec "$PYTHON" "$HERE/01_run_neutral.py" --input "$SOURCE" --output "$OUT" --jobs "${JOBS:-1}" "${EXTRA[@]}"
