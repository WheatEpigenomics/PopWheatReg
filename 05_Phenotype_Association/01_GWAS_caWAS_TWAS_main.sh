#### GWAS
#calcu kinship
input_bed="input"
out="out"
emmax-kin-intel64 $input_bed -v -d 10 -o $out.BN.kinf

#gwas
input_tped="input"
input_kin="input.BN.kinf"
input_pheno="pheno.txt"
gwas_out="gwas_out"
emmax-intel64 -t $input_tped -o $gwas_out -p $input_pheno -k $input_kin


#### caWAS
PHENO="pheno.txt"
ACR="acr_sig.matrix"
outpath="out"

python2 PrediXcan.py \
  --assoc \
  --pheno "$PHENO" \
  --mpheno 1 \
  --missing-phenotype NA \
  --pred_exp "$ACR" \
  --output_dir "$outpath"


#### TWAS
PHENO="pheno.txt"
EXP="exp.matrix"
outpath="out"

python2 PrediXcan.py \
  --assoc \
  --pheno "$PHENO" \
  --mpheno 1 \
  --missing-phenotype NA \
  --pred_exp "$EXP" \
  --output_dir "$outpath"