#!/usr/bin/env bash
# End-to-end demo on synthetic data (no real data, no network).
#
#   bash demo/run_demo.sh
#
# Runs:
#   1) global_pc        : y ~ G + 3 PCs, all pairs, no cis filtering
#   2) subgroup_omnibus : per-subgroup y ~ G + exact omnibus, safe prefilter ON
#   3) merge            : assemble the final subgroup parquet
#
# Set PYTHON=/path/to/python to use a specific interpreter.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="${PYTHON:-python}"
OUT="$HERE/out"

cd "$HERE"
"$PY" make_demo_data.py

rm -rf "$OUT"
mkdir -p "$OUT"

echo
echo "=============================================================="
echo "1) global_pc   y ~ G + PC1..PC3   (all SNP x all ACR, no cis filter)"
echo "=============================================================="
"$PY" "$ROOT/00_TensorQTL_add_subgroup_asso.py" \
  --mode global_pc \
  --genotype      data/genotype.parquet \
  --phenotype     data/phenotype.tsv \
  --covariates    data/covariates.tsv \
  --maf_threshold 0.05 \
  --pval_threshold 1e-3 \
  --batch_size 500 \
  --out "$OUT/global_pc.tsv"

echo
echo "=============================================================="
echo "2) subgroup_omnibus   per-subgroup y ~ G + exact omnibus"
echo "   (safe GPU prefilter ON by default)"
echo "=============================================================="
"$PY" "$ROOT/00_TensorQTL_add_subgroup_asso.py" \
  --mode subgroup_omnibus \
  --genotype      data/genotype.parquet \
  --phenotype     data/phenotype.tsv \
  --subgroups     data/subgroups.tsv \
  --maf_threshold 0.05 \
  --pval_threshold 1e-3 \
  --batch_size 500 \
  --pheno_chunk 200 \
  --out "$OUT/sgom.parquet"

echo
echo "=============================================================="
echo "3) merge subgroup parts -> final parquet"
echo "=============================================================="
"$PY" "$ROOT/03_merge_subgroup_parts.py" --out "$OUT/sgom.parquet"

echo
echo "=============================================================="
echo "results"
echo "=============================================================="
"$PY" - "$OUT" <<'PY'
import sys
from pathlib import Path

import pandas as pd

out = Path(sys.argv[1])
glob = pd.read_csv(out / "global_pc.tsv", sep="\t")
sgom = pd.read_parquet(out / "sgom.parquet")

print(f"global_pc: {len(glob):,} pairs; columns = {list(glob.columns)}")
print(glob.head(5).to_string(index=False))

print()
print(f"subgroup_omnibus: {len(sgom):,} pairs; columns = {list(sgom.columns)}")
show = [c for c in ("variant_id", "phenotype_id", "P_W", "W", "K_est",
                    "P_min", "strongest_subgroup") if c in sgom.columns]
print(sgom[show].head(5).to_string(index=False))
PY

echo
echo "Done. Outputs are in $OUT"
