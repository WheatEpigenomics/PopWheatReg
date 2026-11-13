# Runtime: the bundled demo finishes in seconds; a full run is dominated by the network build and the evidence conversion.
##*## TRACE V2-2 (C-agLINK) regulatory network: copy one block at a time and run it here.
##*## Needs python3 with duckdb + pandas. Inputs are prepared by the user; only the
##*## standard formats in README.md are required. No result files are shipped.


################# 01 inputs #################
## Raw network inputs. See README.md for the exact columns of each file.
eqtl=./data/eqtl_fil_beta.bed          # 9 cols, no header
caqtl=./data/caqtl_source_beta.tsv.gz  # header, source-annotated
aglink_raw=./data/aglink_merge_distance.txt   # gene, ACR, value, p/NA, aglink|dist
coacr_raw=./data/coacces_acr.tsv       # id1, id2, correlation
binding_raw=./data/tf_acr_pcc.tsv      # TF, motif, ACR, PCC, effect
gene_bed=./data/genes.bed              # BED6
tf_list=./data/TF_subfamilies.txt      # TSV with a TF column
all_ann=./data/snpeff_allANN.tsv       # SnpEff ANN table

out=./out_net
edge=${out}/f01_edge
mkdir -p ${edge} ${out}/f02_local ${out}/f03_tf ${out}/f04_tmp


################# 02 standard edge tables #################
## All six network inputs are built here, so nothing has to be pre-made.
## The converters never overwrite an existing output; delete it or use a new name.

## eQTL: chr/start/end/SNP/gene/P/cis-tra/beta/se -> SNP-Gene edge
python3 -B ./TRACE_v2_caller/01_prepare_edge_tables.py edges \
  --kind eqtl --input ${eqtl} --gene-bed ${gene_bed} \
  --out ${edge}/01_eqtl_edge.tsv --qc ${edge}/01_eqtl_edge.qc.tsv

## caQTL: needs the source-annotated table, not the 9-column BED
python3 -B ./TRACE_v2_caller/01_prepare_edge_tables.py edges \
  --kind caqtl --input ${caqtl} --gene-bed ${gene_bed} \
  --out ${edge}/02_caqtl_edge.tsv --qc ${edge}/02_caqtl_edge.qc.tsv

## agLINK: type=aglink/calink -> agLINK (score=PCC); type=dist -> distance (score=NA)
python3 -B ./TRACE_v2_caller/01_prepare_edge_tables.py aglink \
  --input ${aglink_raw} \
  --out ${edge}/03_aglink_edge.tsv --qc ${edge}/03_aglink_edge.qc.tsv

## coACR: ACR-ACR co-accessibility, header line is optional
python3 -B ./TRACE_v2_caller/01_prepare_edge_tables.py coacr \
  --input ${coacr_raw} \
  --out ${edge}/04_coacr_edge.tsv --qc ${edge}/04_coacr_edge.qc.tsv

## TF binding: signed PCC kept, |PCC| < 0.35 and Un dropped by default
python3 -B ./TRACE_v2_caller/01_prepare_edge_tables.py binding \
  --input ${binding_raw} --min-abs-pcc 0.35 \
  --out ${edge}/05_binding_edge.tsv --qc ${edge}/05_binding_edge.qc.tsv

## SNP annotation: append canonical_gene_id, used by the SNP-TF step
python3 -B ./TRACE_v2_caller/01_prepare_edge_tables.py normalize-ann \
  --input ${all_ann} --out ${edge}/08_snpeff_allANN_canonical.tsv \
  --qc ${edge}/08_snpeff_allANN.qc.tsv


################# 03 Local units (L1/L2) #################
## Reads the two QTL edge tables plus the agLINK and coACR edge tables.
## Whole span <= 500 kb: SNP position / ACR floor midpoint / gene strand-aware TSS.
python3 -u ./TRACE_v2_caller/02_build_local_units.py \
  --edges ${edge}/01_eqtl_edge.tsv ${edge}/02_caqtl_edge.tsv \
          ${edge}/03_aglink_edge.tsv ${edge}/04_coacr_edge.tsv \
  --gene-bed ${gene_bed} --span 500000 \
  --out ${out}/f02_local/01_local_units.tsv \
  --qc ${out}/f02_local/01_local_units.qc.tsv \
  --checkpoint-db ${out}/f04_tmp/local.duckdb --tmp-dir ${out}/f04_tmp/local_spill \
  --threads 4

## Only L1/L2 units, no binding: binding never enters the Local span calculation.


################# 04 SNP-TF edges #################
## Two sources only: Local units and non-intergenic proximal SnpEff annotation.
python3 -B ./TRACE_v2_caller/03_build_snp_tf_edges.py \
  --local-units ${out}/f02_local/01_local_units.tsv \
  --snpeff ${edge}/08_snpeff_allANN_canonical.tsv \
  --tf-list ${tf_list} \
  --out ${out}/f03_tf/01_TF_SNP.edge.tsv \
  --qc ${out}/f03_tf/02_TF_SNP_QC.tsv


################# 05 TF network #################
## Routes 1/2 are unchanged; Route 3 keeps a SNP-TF pair when it comes from a
## Local unit, or when it is proximal AND within 5 kb AND supported by a
## same-SNP cis-eQTL or a cis-caQTL -> ACR -> direct agLINK chain (<= 500 kb).
python3 -u ./TRACE_v2_caller/05_build_caglink_network.py \
  --eqtl ${edge}/01_eqtl_edge.tsv \
  --caqtl ${edge}/02_caqtl_edge.tsv \
  --aglink ${edge}/03_aglink_edge.tsv \
  --coacr ${edge}/04_coacr_edge.tsv \
  --binding ${edge}/05_binding_edge.tsv \
  --tf-snp ${out}/f03_tf/01_TF_SNP.edge.tsv \
  --gene-bed ${gene_bed} \
  --out-prefix ${out}/f03_tf/03_TF_network \
  --same-chr-trans-min-distance 500000 --local-span 500000 \
  --checkpoint-db ${out}/f04_tmp/network.duckdb --tmp-dir ${out}/f04_tmp/network_spill \
  --ram-limit 256GB --threads 4 --batch-snps 1000 --evidence-part-net-ids 200000

## ${out}/f03_tf/03_TF_network.network.tsv        main table, TFN000000001 ids
## ${out}/f03_tf/03_TF_network.evidence.part*.parquet
## ${out}/f03_tf/03_TF_network.filter_qc.tsv      per-rule losses
## ${out}/f03_tf/01_Caglink_SNP_TF_entrances.tsv / 02_Caglink_entrance_QC.tsv


################# 06 trans evidence format #################
python3 -u ./TRACE_v2_caller/06_finalize_evidence.py \
  --prefix ${out}/f03_tf/03_TF_network \
  --out ${out}/f03_tf/04_TF_network.evidence.tsv.gz \
  --qc ${out}/f03_tf/05_evidence_text_QC.tsv \
  --batch-rows 256 --compresslevel 1 --delete-parts


################# 07 score #################
## score = number of distinct Route labels behind one TF-gene relation.
## sorted by score DESC, n_support DESC; net_id is kept unchanged.
python3 -u ./TRACE_v2_caller/07_score_network.py \
  --input ${out}/f03_tf/03_TF_network.network.tsv \
  --output ${out}/f03_tf/06_TF_network_main_scored_simple.tsv \
  --ram-limit 64GB


################# 08 demo #################
## Toy dataset: one example file per input type, same seven blocks, output to demo/out.
bash ./demo/run_demo.sh

## To run the blocks above on the demo data instead, point the block 01 variables
## at demo/data/ and keep everything else unchanged:
#  eqtl=./demo/data/eqtl_fil_beta.bed
#  caqtl=./demo/data/caqtl_source_beta.tsv.gz
#  aglink_raw=./demo/data/aglink_merge_distance.txt
#  coacr_raw=./demo/data/coacces_acr.tsv
#  binding_raw=./demo/data/tf_acr_pcc.tsv
#  gene_bed=./demo/data/genes.bed
#  tf_list=./demo/data/TF_subfamilies.txt
#  all_ann=./demo/data/snpeff_allANN.tsv
#  out=./demo/out_net
