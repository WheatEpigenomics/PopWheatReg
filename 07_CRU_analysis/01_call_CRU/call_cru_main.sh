# Runtime: the demo finishes in under a minute; full calling parallelises over chromosomes.
##*## CRU construction and scoring: copy one block at a time and run it here.
##*## Needs python3 with numpy + scipy. Blocks 05/06 and the figures need R with
##*## ggplot2 and a cairo device; calling itself never depends on R.
##*## Inputs are the same files used by the ACR-gene link and co-accessibility steps.


################# 01 inputs #################
aglink=./data/aglink.tsv            # gene, ACR, PCC, p            (no header)
gene_bed=./data/genes.bed           # BED6: chr start end gene_id score strand
acr_bed=./data/acrs.bed             # BED4: chr start end acr_id
genome_sizes=./data/genome.sizes    # chr <TAB> size
block_bed=./data/blocks.bed         # BED3 haplotype blocks, for block test and figures
coacr=./data/coacr.tsv              # id1, id2, correlation, figures only

cru_out=./out_cru
score_out=./out_score


################# 02 call CRUs #################
## calling only: no figures, no block test. Add --chroms 1A to try one chromosome.
python3 -B ./f01_T8_CRU_caller/01_run_cru.py \
  --aglink ${aglink} --gene-bed ${gene_bed} --acr-bed ${acr_bed} \
  --genome-sizes ${genome_sizes} --config ./f01_T8_CRU_caller/f01_rules/CRU_rules.json \
  --chroms all --workers 21 --control none \
  --mode all --plot none --block none \
  --outdir ${cru_out} --resume

## CRUs:        ${cru_out}/f04_merged/CRU.bed
## per chr:     ${cru_out}/f02_chromosomes/<chr>/f02_calling/<chr>.CRU.bed
## diagnostics: ${cru_out}/f02_chromosomes/<chr>/f03_diagnostics/09_summary.tsv
## rerun after a crash: same command, finished chromosomes are skipped by MD5


################# 03 haplotype-block test #################
## 1000 permutations per chromosome; independent of calling, reads the merged BED
python3 -B ./f01_T8_CRU_caller/01_run_cru.py \
  --aglink ${aglink} --gene-bed ${gene_bed} --acr-bed ${acr_bed} \
  --genome-sizes ${genome_sizes} --config ./f01_T8_CRU_caller/f01_rules/CRU_rules.json \
  --block-bed ${block_bed} --chroms all --workers 21 --control none \
  --mode block --plot none --block all --block-perm-n 1000 --local-window-bp 10000000 \
  --outdir ${cru_out} --resume

## ${cru_out}/f06_block_enrichment/f01_statistics/02_enrichment_summary.tsv


################# 04 region figures (optional) #################
## 5-Mb windows per chromosome. A failed figure is reported as PLOT FAILED and
## never blocks calling or merging; redraw later with the same command.
python3 -B ./f01_T8_CRU_caller/01_run_cru.py \
  --aglink ${aglink} --gene-bed ${gene_bed} --acr-bed ${acr_bed} \
  --genome-sizes ${genome_sizes} --config ./f01_T8_CRU_caller/f01_rules/CRU_rules.json \
  --block-bed ${block_bed} --coacr ${coacr} --chroms all --workers 21 --control none \
  --mode plot --plot all --block none --window-bp 5000000 \
  --rscript "$(command -v Rscript)" --font-family 'Liberation Sans' --dpi 120 \
  --outdir ${cru_out} --resume

## ${cru_out}/f02_chromosomes/<chr>/f05_figures/


################# 05 score the fixed CRUs #################
## Reuses the frozen f25 formula with the T7-A normalization anchors; the CRU
## BED is never modified. Weights: agLINK density .40, co-accessibility density .25,
## agLINK independence .20, co-accessibility independence .15.
python3 -B ./f01_T8_CRU_caller/10_score_cru.py \
  --cru-bed ${cru_out}/f04_merged/CRU.bed \
  --aglink ${aglink} --acr-bed ${acr_bed} --gene-bed ${gene_bed} \
  --genome-sizes ${genome_sizes} --coacr ${coacr} \
  --outdir ${score_out}

## all metrics:  ${score_out}/f03_scores/f01_scores/01_CRU_regulatory_scores.tsv
## BED6, score:  ${score_out}/02_CRU_scored.bed          (col 5 = round(score*10), 0-1000)
## simple table: ${score_out}/04_CRU_scores_simple.tsv   (chr start end cru_id score)


################# 06 score figure (optional) #################
Rscript --vanilla ./f01_T8_CRU_caller/11_plot_score.R \
  --scores ${score_out}/f03_scores/f01_scores/01_CRU_regulatory_scores.tsv \
  --outdir ${score_out}/f05_figures --font-family 'Liberation Sans'


################# 07 demo #################
## 20-Mb toy genome, 508 ACRs, 1,467 agLinks; calling and scoring take seconds.
sh ./f01_T8_CRU_caller/demo/run_demo.sh

