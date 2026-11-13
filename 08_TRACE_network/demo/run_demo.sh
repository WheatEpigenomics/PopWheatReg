#!/usr/bin/env bash
# End-to-end demo on a toy dataset: two TF genes, two target genes, six ACRs
# and two SNPs, small enough to run in a few seconds with no external data.
#
#   bash demo/run_demo.sh [outdir]
#
# Same seven blocks, same variable names and same output layout as
# call_trace_main.sh; only the paths differ. Output defaults to demo/out.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(dirname "$HERE")"
PY="${PYTHON_BIN:-python3}"
CODE="$ROOT/TRACE_v2_caller"

eqtl="$HERE/data/eqtl_fil_beta.bed"
caqtl="$HERE/data/caqtl_source_beta.tsv.gz"
aglink_raw="$HERE/data/aglink_merge_distance.txt"
coacr_raw="$HERE/data/coacces_acr.tsv"
binding_raw="$HERE/data/tf_acr_pcc.tsv"
gene_bed="$HERE/data/genes.bed"
tf_list="$HERE/data/TF_subfamilies.txt"
all_ann="$HERE/data/snpeff_allANN.tsv"

out="${1:-$HERE/out}"
edge=${out}/f01_edge
rm -rf "${out}"
mkdir -p "${edge}" "${out}/f02_local" "${out}/f03_tf" "${out}/f04_tmp"

echo "== 02 standard edge tables =="
"$PY" -B "$CODE/01_prepare_edge_tables.py" edges \
  --kind eqtl --input "${eqtl}" --gene-bed "${gene_bed}" \
  --out "${edge}/01_eqtl_edge.tsv" --qc "${edge}/01_eqtl_edge.qc.tsv"
"$PY" -B "$CODE/01_prepare_edge_tables.py" edges \
  --kind caqtl --input "${caqtl}" --gene-bed "${gene_bed}" \
  --out "${edge}/02_caqtl_edge.tsv" --qc "${edge}/02_caqtl_edge.qc.tsv"
"$PY" -B "$CODE/01_prepare_edge_tables.py" aglink \
  --input "${aglink_raw}" \
  --out "${edge}/03_aglink_edge.tsv" --qc "${edge}/03_aglink_edge.qc.tsv"
"$PY" -B "$CODE/01_prepare_edge_tables.py" coacr \
  --input "${coacr_raw}" \
  --out "${edge}/04_coacr_edge.tsv" --qc "${edge}/04_coacr_edge.qc.tsv"
"$PY" -B "$CODE/01_prepare_edge_tables.py" binding \
  --input "${binding_raw}" --min-abs-pcc 0.35 \
  --out "${edge}/05_binding_edge.tsv" --qc "${edge}/05_binding_edge.qc.tsv"
"$PY" -B "$CODE/01_prepare_edge_tables.py" normalize-ann \
  --input "${all_ann}" --out "${edge}/08_snpeff_allANN_canonical.tsv" \
  --qc "${edge}/08_snpeff_allANN.qc.tsv"

echo "== 03 local units =="
"$PY" -B "$CODE/02_build_local_units.py" \
  --edges "${edge}/01_eqtl_edge.tsv" "${edge}/02_caqtl_edge.tsv" \
          "${edge}/03_aglink_edge.tsv" "${edge}/04_coacr_edge.tsv" \
  --gene-bed "${gene_bed}" --span 500000 \
  --out "${out}/f02_local/01_local_units.tsv" \
  --qc "${out}/f02_local/01_local_units.qc.tsv"

echo "== 04 SNP-TF edges =="
"$PY" -B "$CODE/03_build_snp_tf_edges.py" \
  --local-units "${out}/f02_local/01_local_units.tsv" \
  --snpeff "${edge}/08_snpeff_allANN_canonical.tsv" \
  --tf-list "${tf_list}" \
  --out "${out}/f03_tf/01_TF_SNP.edge.tsv" \
  --qc "${out}/f03_tf/02_TF_SNP_QC.tsv"

echo "== 05 C-agLINK network =="
"$PY" -B "$CODE/05_build_caglink_network.py" \
  --eqtl "${edge}/01_eqtl_edge.tsv" \
  --caqtl "${edge}/02_caqtl_edge.tsv" \
  --aglink "${edge}/03_aglink_edge.tsv" \
  --coacr "${edge}/04_coacr_edge.tsv" \
  --binding "${edge}/05_binding_edge.tsv" \
  --tf-snp "${out}/f03_tf/01_TF_SNP.edge.tsv" \
  --gene-bed "${gene_bed}" \
  --out-prefix "${out}/f03_tf/03_TF_network" \
  --checkpoint-db "${out}/f04_tmp/network.duckdb" --tmp-dir "${out}/f04_tmp/network_spill" \
  --ram-limit 4GB --threads 1 --batch-snps 100 --evidence-part-net-ids 1000

echo "== 06 evidence to text =="
"$PY" -B "$CODE/06_finalize_evidence.py" \
  --prefix "${out}/f03_tf/03_TF_network" \
  --out "${out}/f03_tf/04_TF_network.evidence.tsv.gz" \
  --qc "${out}/f03_tf/05_evidence_text_QC.tsv" --delete-parts

echo "== 07 score =="
"$PY" -B "$CODE/07_score_network.py" \
  --input "${out}/f03_tf/03_TF_network.network.tsv" \
  --output "${out}/f03_tf/06_TF_network_main_scored_simple.tsv" --ram-limit 4GB

echo
echo "network: ${out}/f03_tf/06_TF_network_main_scored_simple.tsv"
column -t -s $'\t' "${out}/f03_tf/06_TF_network_main_scored_simple.tsv" 2>/dev/null \
  || cat "${out}/f03_tf/06_TF_network_main_scored_simple.tsv"
