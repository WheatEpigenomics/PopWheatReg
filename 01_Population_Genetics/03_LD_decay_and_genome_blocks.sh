# Runtime: PLINK r2 on thinned SNPs and block detection per subgenome; production scale is in README.md.
#### Genome-wide haplotype blocks from LD decay
# manual run; append new species/dataset sections below
# requires: bcftools, plink 1.9, bedtools, python (+ ld_decay.py next to this file)


################# wheat csv2 ##################
VCF=./input.vcf.gz              # biallelic SNPs, MAF >= 0.02, missing <= 0.3
OUT=./wheat_csv2_block
LD_DECAY=./ld_decay.py          # bins LD by distance, reports mean r2
THREADS=16
SEED=20240623
BLOCK_TAG=wheat_csv2

SPLIT=${OUT}/f01_split_subgenome
LD=${OUT}/f02_sub_calcu_ld
BLK=${OUT}/f03_sub_block
mkdir -p ${SPLIT} ${LD} ${BLK}


### step1: split A/B/D subgenomes, then convert to PLINK bed
bcftools view -t 1A,2A,3A,4A,5A,6A,7A -Oz -o ${SPLIT}/chrA.vcf.gz ${VCF} &
bcftools view -t 1B,2B,3B,4B,5B,6B,7B -Oz -o ${SPLIT}/chrB.vcf.gz ${VCF} &
bcftools view -t 1D,2D,3D,4D,5D,6D,7D -Oz -o ${SPLIT}/chrD.vcf.gz ${VCF} &
wait
bcftools index -f ${SPLIT}/chrA.vcf.gz &
bcftools index -f ${SPLIT}/chrB.vcf.gz &
bcftools index -f ${SPLIT}/chrD.vcf.gz &
wait

plink --vcf ${SPLIT}/chrA.vcf.gz --double-id --make-bed --allow-extra-chr --out ${SPLIT}/chrA &
plink --vcf ${SPLIT}/chrB.vcf.gz --double-id --make-bed --allow-extra-chr --out ${SPLIT}/chrB &
plink --vcf ${SPLIT}/chrD.vcf.gz --double-id --make-bed --allow-extra-chr --out ${SPLIT}/chrD &
wait

## subsample SNPs for LD decay (full data too slow); fixed seed for reproducibility
plink --bfile ${SPLIT}/chrA --thin 0.05 --seed ${SEED} --allow-extra-chr --make-bed --out ${SPLIT}/chrA.thin &
plink --bfile ${SPLIT}/chrB --thin 0.05 --seed ${SEED} --allow-extra-chr --make-bed --out ${SPLIT}/chrB.thin &
plink --bfile ${SPLIT}/chrD --thin 0.20 --seed ${SEED} --allow-extra-chr --make-bed --out ${SPLIT}/chrD.thin &
wait


### step2: pairwise r2 within subgenome (A/B <= 5 Mb, D <= 30 Mb)
plink --bfile ${SPLIT}/chrA.thin --r2 gz --ld-window-kb 5000 --ld-window 99999 --ld-window-r2 0 --allow-extra-chr --threads ${THREADS} --out ${LD}/chrA.ld &
plink --bfile ${SPLIT}/chrB.thin --r2 gz --ld-window-kb 5000 --ld-window 99999 --ld-window-r2 0 --allow-extra-chr --threads ${THREADS} --out ${LD}/chrB.ld &
plink --bfile ${SPLIT}/chrD.thin --r2 gz --ld-window-kb 30000 --ld-window 99999 --ld-window-r2 0 --allow-extra-chr --threads ${THREADS} --out ${LD}/chrD.ld &
wait


### step3: LD decay; report the distance where mean r2 drops below 0.2
## bin 10 kb for A/B, 20 kb for D
python ${LD_DECAY} ${LD}/chrA.ld.ld.gz 10000 0.2 && mv ld_decay.tsv ${LD}/chrA_ld_decay.tsv && mv ld_decay.png ${LD}/chrA_ld_decay.png
python ${LD_DECAY} ${LD}/chrB.ld.ld.gz 10000 0.2 && mv ld_decay.tsv ${LD}/chrB_ld_decay.tsv && mv ld_decay.png ${LD}/chrB_ld_decay.png
python ${LD_DECAY} ${LD}/chrD.ld.ld.gz 20000 0.2 && mv ld_decay.tsv ${LD}/chrD_ld_decay.tsv && mv ld_decay.png ${LD}/chrD_ld_decay.png

## --blocks-max-kb from the decay distance: A/B ~5-8x, D ~2-3x
## wheat csv2 final: A=9000, B=10000, D=15000 (kb)


### step4: PLINK haplotype blocks
## V1 parameters; stricter settings (strong-lowci 0.68 etc.) gave too many small blocks
plink --bfile ${SPLIT}/chrA --allow-extra-chr \
  --blocks no-pheno-req no-small-max-span --blocks-min-maf 0.02 --blocks-max-kb 9000 \
  --blocks-strong-lowci 0.62 --blocks-strong-highci 0.92 --blocks-recomb-highci 0.92 \
  --blocks-inform-frac 0.80 --threads ${THREADS} --out ${BLK}/chrA.blocks &
plink --bfile ${SPLIT}/chrB --allow-extra-chr \
  --blocks no-pheno-req no-small-max-span --blocks-min-maf 0.02 --blocks-max-kb 10000 \
  --blocks-strong-lowci 0.62 --blocks-strong-highci 0.92 --blocks-recomb-highci 0.92 \
  --blocks-inform-frac 0.80 --threads ${THREADS} --out ${BLK}/chrB.blocks &
plink --bfile ${SPLIT}/chrD --allow-extra-chr \
  --blocks no-pheno-req no-small-max-span --blocks-min-maf 0.02 --blocks-max-kb 15000 \
  --blocks-strong-lowci 0.62 --blocks-strong-highci 0.92 --blocks-recomb-highci 0.92 \
  --blocks-inform-frac 0.80 --threads ${THREADS} --out ${BLK}/chrD.blocks &
wait


## merge blocks, add block IDs, sort
cat ${BLK}/chrA.blocks.blocks.det ${BLK}/chrB.blocks.blocks.det ${BLK}/chrD.blocks.blocks.det \
  | grep -v "CHR" \
  | sort -k1,1 -k2,2 \
  | sed 's/ /\t/g' | sed 's/\t\{1,\}/\t/g' | sed 's/^[^\t]*\t//' \
  > ${BLK}/merge.block

bedtools sort -i ${BLK}/merge.block \
  | awk -v tag="${BLOCK_TAG}" '{printf("%s\t%s\t%s\t%s_block%04d:%s_%s_%s\n", $1,$2,$3,tag,NR,$1,$2,$3)}' \
  > ${OUT}/${BLOCK_TAG}.blocks.bed

awk '{print $3-$2}' ${OUT}/${BLOCK_TAG}.blocks.bed > ${OUT}/${BLOCK_TAG}.blocks.length.txt
