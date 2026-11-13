# Runtime: the per-gene Pearson scan plus the 10,000-permutation null dominates; demo runs in seconds, production scale is in README.md.
#### agLink: gene-ACR association (manual run, 4-script version)
# Copy the blocks below and run them step by step.
# Requires: python3 (numpy/scipy/pandas), bedtools, shuf.
# 04_co_access_analysis.py needs CUDA; the SMR block needs SMR v1.4.1 + PLINK bed.


################# demo test #################
mkdir -p demo_out
python 01_calcu_acr_gene_pcc.py demo/aglink/demo_gene_real_matrix.txt > demo_out/demo.pcc.txt
paste demo/aglink/demo_acr_good.list demo_out/demo.pcc.txt \
  | awk '{print "demo_gene",$1,$2}' > demo_out/demo.real_cor.txt

# permutation control from the same toy matrix (100 permutations x 4 ACRs)
python - demo/aglink/demo_gene_real_matrix.txt demo_out/demo.fake_cor.txt <<'PY'
import sys
import numpy as np
m = np.loadtxt(sys.argv[1])
gene, acrs = m[0], m[1:]
rng = np.random.default_rng(1)
vals = []
for _ in range(100):
    for a in acrs:
        r = np.corrcoef(gene, rng.permutation(a))[0, 1]
        if np.isfinite(r):
            vals.append(r)
with open(sys.argv[2], "w") as f:
    f.write("demo_gene " + " ".join(f"{v:.4f}" for v in vals) + "\n")
PY

sh 02_fill_column_for_ztest.sh demo_out/demo.fake_cor.txt demo_out/demo.fake_cor_fill.txt 401 0
python 03_aglink_ztest.py demo_out/demo.fake_cor_fill.txt demo_out/demo.real_cor.txt demo_out/demo.ztest.txt
python 04_co_access_analysis.py \
  --matrix demo/coacr/demo_acr_matrix.txt \
  --bed demo/coacr/demo_acr.bed \
  --window 5000 --method pearson --cutoff 0.6 \
  --gpu-mem auto --device cuda:0 \
  --output demo_out/demo.coaccess.tsv


################# agLink #################
## inputs
expre_file=./gene_expression_matrix.txt     # genes x samples, space-delimited, no header
acr_file=./acr_signal_matrix.txt            # ACRs x samples, first col = ACR id, space-delimited, no header
gene_bed=./genes.bed                        # BED6: chr start end gene_id score strand
acr_bed=./acrs.bed                          # BED4: chr start end acr_id
outdir=./aglink_out

n_fake=10000
window=1000000
pcc_cut=0.4
p_cut=0.05
seed=20240623

mkdir -p ${outdir}/step1_gene_acr
mkdir -p ${outdir}/step2_matrix/real ${outdir}/step2_matrix/fake
mkdir -p ${outdir}/step3_person/real ${outdir}/step3_person/fake
mkdir -p ${outdir}/step4_z ${outdir}/step5_sig

##*## Step1: candidate ACRs within +/-1 Mb of each gene TSS
awk -v w=${window} 'BEGIN{OFS="\t"} {
  if ($6=="-") { s=($3-w)<1?1:($3-w); e=$3 }
  else         { s=($2-w)<1?1:($2-w); e=$2+w }
  print $1,s,e,$4,$5,$6
}' ${gene_bed} > ${outdir}/step1_gene_acr/tss_1m.bed
cut -f4 ${outdir}/step1_gene_acr/tss_1m.bed > ${outdir}/step1_gene_acr/gene.list

for i in $(cat ${outdir}/step1_gene_acr/gene.list)
do
  grep -m1 -w $i ${outdir}/step1_gene_acr/tss_1m.bed \
    | bedtools intersect -a - -b ${acr_bed} -wa -wb \
    | cut -f10 | sort -u > ${outdir}/step1_gene_acr/${i}.acr.list
done

##*## Step2/3: real and fake matrices + PCC
rm -f ${outdir}/step3_person/real/all_real_cor.txt
rm -f ${outdir}/step3_person/fake/all_fake_cor.txt

for i in $(cat ${outdir}/step1_gene_acr/gene.list)
do
  chr=$(grep -m1 -w $i ${outdir}/step1_gene_acr/tss_1m.bed | cut -f1)
  acr_list=${outdir}/step1_gene_acr/${i}.acr.list
  [ -s $acr_list ] || continue

  # real matrix: gene expression + candidate ACRs
  grep -m1 -w $i ${expre_file} | cut -d' ' -f2- > ${outdir}/step2_matrix/real/${i}.real_matrix.txt
  for x in $(cat $acr_list)
  do
    grep -m1 -w $x ${acr_file} | cut -d' ' -f2- >> ${outdir}/step2_matrix/real/${i}.real_matrix.txt
  done
  python 01_calcu_acr_gene_pcc.py ${outdir}/step2_matrix/real/${i}.real_matrix.txt > ${outdir}/step3_person/real/${i}.pcc.txt
  paste $acr_list ${outdir}/step3_person/real/${i}.pcc.txt \
    | awk -v g=$i '{print g"\t"$1"\t"$2}' >> ${outdir}/step3_person/real/all_real_cor.txt

  # fake matrix: same gene + random ACRs from other chromosomes
  cp ${outdir}/step2_matrix/real/${i}.real_matrix.txt ${outdir}/step2_matrix/fake/${i}.fake_matrix.txt
  grep -v -w $chr ${acr_file} | cut -d' ' -f2- \
    | shuf -n ${n_fake} --random-source=<(yes ${seed}) >> ${outdir}/step2_matrix/fake/${i}.fake_matrix.txt
  fk_cor=$(python 01_calcu_acr_gene_pcc.py ${outdir}/step2_matrix/fake/${i}.fake_matrix.txt | grep -v nan | xargs)
  printf "$i $fk_cor\n" >> ${outdir}/step3_person/fake/all_fake_cor.txt

  rm -f ${outdir}/step2_matrix/real/${i}.real_matrix.txt ${outdir}/step2_matrix/fake/${i}.fake_matrix.txt
done

##*## Step4: pad control rows to n_fake values
sh 02_fill_column_for_ztest.sh \
  ${outdir}/step3_person/fake/all_fake_cor.txt \
  ${outdir}/step3_person/fake/all_fake_cor_fill.txt \
  $((n_fake+1)) 0

##*## Step5: z-test (real PCC vs control PCC)
python 03_aglink_ztest.py \
  ${outdir}/step3_person/fake/all_fake_cor_fill.txt \
  ${outdir}/step3_person/real/all_real_cor.txt \
  ${outdir}/step4_z/aglink_ztest.txt

##*## Step6: significant agLinks
awk -F'\t' -v c=${pcc_cut} -v p=${p_cut} 'NR==1 || (($6>c)||($6<-c)) && ($7<p)' \
  ${outdir}/step4_z/aglink_ztest.txt > ${outdir}/step5_sig/aglink_significant.txt



################# aglink SMR/HEIDI #################
## cis_all legacy output (10 columns, no header); SMR v1.4.1; final HEIDI nsnp >= 3
eqtl_cis=./eqtl_cis_legacy.tsv.gz
caqtl_cis=./caqtl_cis_legacy.tsv.gz
bim=./csv2_snp_203hex_maf002_miss03.bim
bfile=./csv2_snp_203hex_maf002_miss03
gene_bed=./Ta_CS_v2.1.gene.bed
aglink=./aglink_out/step5_sig/aglink_significant.txt

mkdir -p ./smr_out
cd ./smr_out

##*## convert to matrix-eQTL (gene=$1 SNP=$2 p=$8 beta=$9 se=$10)
gunzip -c ${eqtl_cis} | LC_ALL=C awk '
  BEGIN{OFS="\t"; print "SNP","gene","beta","t-stat","p-value","FDR"}
  NF>=10 { beta=$9+0; se=$10+0; p=$8;
           if (se>0 && p!="NA" && p!="nan" && p!="NaN" && p>0) print $2,$1,beta,beta/se,p,1 }' \
  > eqtl.matrixeqtl.txt
gunzip -c ${caqtl_cis} | LC_ALL=C awk '
  BEGIN{OFS="\t"; print "SNP","gene","beta","t-stat","p-value","FDR"}
  NF>=10 { beta=$9+0; se=$10+0; p=$8;
           if (se>0 && p!="NA" && p!="nan" && p!="NaN" && p>0) print $2,$1,beta,beta/se,p,1 }' \
  > caqtl.matrixeqtl.txt

## drop t-stat = 0
awk 'NR==1 || $4!=0' eqtl.matrixeqtl.txt > eqtl.matrixeqtl.nonzero.txt
awk 'NR==1 || $4!=0' caqtl.matrixeqtl.txt > caqtl.matrixeqtl.nonzero.txt

##*## BESD
smr --eqtl-summary eqtl.matrixeqtl.nonzero.txt  --matrix-eqtl-format --make-besd --out eqtl_input
smr --eqtl-summary caqtl.matrixeqtl.nonzero.txt --matrix-eqtl-format --make-besd --out caqtl_input

##*## fill esi (SNP position/alleles from BIM)
for PREFIX in eqtl_input caqtl_input
do
  cp ${PREFIX}.esi ${PREFIX}.esi.raw
  awk 'BEGIN{OFS="\t"}
       NR==FNR { split($2,a,"-"); info[$2]=a[1] OFS $2 OFS 0 OFS $4 OFS $5 OFS $6 OFS "NA"; next }
       { if ($2 in info) print info[$2]; else print "NA",$2,0,"NA","NA","NA","NA" }' \
    ${bim} ${PREFIX}.esi.raw > ${PREFIX}.esi
done

##*## fill epi: eQTL gene TSS (strand-aware)
cp eqtl_input.epi eqtl_input.epi.raw
awk 'BEGIN{OFS="\t"}
     NR==FNR { s=($6=="-")?$3:$2+1; info[$4]=$1 OFS $4 OFS 0 OFS s OFS $4 OFS $6; next }
     { print (($2 in info)?info[$2]:"NA" OFS $2 OFS 0 OFS "NA" OFS $2 OFS "NA") }' \
  ${gene_bed} eqtl_input.epi.raw > eqtl_input.epi

##*## fill epi: caQTL ACR midpoint (id encodes chr_start_end)
cp caqtl_input.epi caqtl_input.epi.raw
awk 'BEGIN{OFS="\t"}
     { split($2,a,":"); split(a[2],b,"_");
       if (b[1]!="" && b[2]!="" && b[3]!="") print b[1],$2,0,int((b[2]+b[3])/2),$2,"NA";
       else print "NA",$2,0,"NA",$2,"NA" }' \
  caqtl_input.epi.raw > caqtl_input.epi

##*## agLink gene / ACR lists
cut -f1 ${aglink} | tail -n +2 | sort -u > aglink_gene.list
cut -f2 ${aglink} | tail -n +2 | sort -u > aglink_acr.list

##*## SMR: caQTL as exposure, eQTL as outcome
smr --bfile ${bfile} \
  --beqtl-summary caqtl_input \
  --beqtl-summary eqtl_input \
  --extract-exposure-probe aglink_acr.list \
  --extract-outcome-probe aglink_gene.list \
  --cis-wind 1000 --thread-num 48 \
  --out aglink_caqtl_eqtl_smr

##*## final SMR/HEIDI pass: p_SMR<0.05, p_HEIDI>0.05, nsnp_HEIDI>=3
## .smr columns: 23=p_SMR, 24=p_HEIDI, 25=nsnp_HEIDI
awk -F'\t' -v OFS='\t' 'NR==1 || ($23<0.05 && $24>0.05 && $25>=3)' \
  aglink_caqtl_eqtl_smr.smr > aglink_smr_heidi_pass.tsv



############### coacr #################
## ACR-ACR co-accessibility (needs CUDA)
coacr_acr_file=./acr_signal_matrix.txt      # ACRs x samples, first col = ACR id
coacr_bed=./acrs.bed                        # BED4: chr start end acr_id
coacr_out=./coacr_out
coacr_window=1000000
coacr_cutoff=0.6

mkdir -p ${coacr_out}

## 04 reads a header row: ACR_ID + sample columns
awk 'NR==1{printf "ACR_ID"; for(i=2;i<=NF;i++) printf " S%d", i-1; print ""} {print}' \
  ${coacr_acr_file} > ${coacr_out}/acr_matrix.txt

python 04_co_access_analysis.py \
  --matrix ${coacr_out}/acr_matrix.txt \
  --bed ${coacr_bed} \
  --window ${coacr_window} --method pearson --cutoff ${coacr_cutoff} \
  --gpu-mem auto --device cuda:0 \
  --output ${coacr_out}/coaccess.tsv
