#### cis-eQTL/caQTL analysis
plink_prefix_path=input_snp
pheno_bed=input.bed
prefix=output

python3 -m tensorqtl \
  ${plink_prefix_path} ${pheno_bed} ${prefix} \
  --mode cis_nominal \
  --window 5000000 \
  --maf_threshold 0.02 \
  --output_dir ./output/cis


#### trans-eQTL/caQTL analysis
mkdir -p ./output/trans
plink_prefix_path=input_snp
pheno_bed=input.bed
prefix=output

python3 -m tensorqtl \
  ${plink_prefix_path} ${pheno_bed} ${prefix} \
  --mode trans \
  --pval_threshold 1e-6 \
  --maf_threshold 0.02 \
  --batch_size 1024 \
  --output_dir ./output/trans
