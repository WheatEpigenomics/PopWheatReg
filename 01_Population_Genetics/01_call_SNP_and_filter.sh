# Runtime: dominated by HaplotypeCaller and the 100-600 GB joint-genotyping steps; production scale and memory are in README.md.
#### GATK calling (GVCF -> joint genotyping)
REF=./ref.fa
JAVA_MEM=40G
THREADS=16
OUT=./snp_filter

mkdir -p tmp gvcf
gatk --java-options "-Xmx30G -Djava.io.tmpdir=./tmp" HaplotypeCaller \
  -R ${REF} -I ${name}.bam -O gvcf/${name}.g.vcf \
  -ERC GVCF -stand-call-conf 30

find gvcf -name "*.g.vcf" > gvcf.list
gatk --java-options -Xmx600G CombineGVCFs -R ${REF} -V gvcf.list -O combined.g.vcf
gatk --java-options -Xmx100G GenotypeGVCFs -R ${REF} -V combined.g.vcf -O raw.vcf


#### SNP filtering (same flow/tools as the production csv2 run)
mkdir -p ${OUT}/01_chr_raw ${OUT}/02_mark ${OUT}/03_pass ${OUT}/04_biallelic ${OUT}/05_mask ${OUT}/06_final

for chr in $(bcftools query -f '%CHROM\n' raw.vcf | sort -u)
do
  ## split by chromosome
  bcftools view -r ${chr} -Ov -o ${OUT}/01_chr_raw/${chr}.vcf raw.vcf

  ## hard filter
  gatk --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=./tmp" VariantFiltration \
    -R ${REF} -V ${OUT}/01_chr_raw/${chr}.vcf -O ${OUT}/02_mark/${chr}.vcf \
    --filter-expression 'QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0' \
    --filter-name lowQualFilter --missing-values-evaluate-as-failing \
    --create-output-variant-index false

  ## keep PASS variants only
  gatk --java-options "-Xmx${JAVA_MEM} -Djava.io.tmpdir=./tmp" SelectVariants \
    -R ${REF} -V ${OUT}/02_mark/${chr}.vcf --exclude-filtered true \
    -O ${OUT}/03_pass/${chr}.vcf --create-output-variant-index false

  ## biallelic SNPs
  bcftools view -m2 -M2 -v snps -Oz -o ${OUT}/04_biallelic/${chr}.vcf.gz ${OUT}/03_pass/${chr}.vcf
  bcftools index -c -f ${OUT}/04_biallelic/${chr}.vcf.gz

  ## genotype-level masking (DP<3, DP>1000, heterozygous) + tags
  bcftools +setGT ${OUT}/04_biallelic/${chr}.vcf.gz -Ou -- -t q -n . -i 'FMT/DP<3' \
    | bcftools +setGT -Ou -- -t q -n . -i 'FMT/DP>1000' \
    | bcftools +setGT -Ou -- -t q -n . -i 'FMT/GT ~ "0[|/]1"' \
    | bcftools +setGT -Ou -- -t q -n . -i 'FMT/GT ~ "1[|/]0"' \
    | bcftools +fill-tags -Oz -o ${OUT}/05_mask/${chr}.vcf.gz -- -t AC,AN,AF,MAF,F_MISSING
  bcftools index -c -f ${OUT}/05_mask/${chr}.vcf.gz

  ## SNP ID + recalculated tags, then final MAF / missingness filter
  bcftools annotate --set-id '%CHROM-%POS' -Ou ${OUT}/05_mask/${chr}.vcf.gz \
    | bcftools +fill-tags -Oz -o ${OUT}/06_final/${chr}.tagged.vcf.gz -- -t AC,AN,AF,MAF,F_MISSING
  bcftools index -c -f ${OUT}/06_final/${chr}.tagged.vcf.gz

  ## optional: sample subset (production splits 210all / 203hex at this point)
  # bcftools view -S samples.txt -Oz -o ${OUT}/06_final/${chr}.subset.vcf.gz ${OUT}/06_final/${chr}.tagged.vcf.gz

  bcftools view -i "AN>0 && AC>0 && MAF>=0.02 && F_MISSING<=0.3" \
    -Oz -o ${OUT}/06_final/${chr}.maf002.miss03.vcf.gz ${OUT}/06_final/${chr}.tagged.vcf.gz
  bcftools index -c -f ${OUT}/06_final/${chr}.maf002.miss03.vcf.gz
done

## concatenate chromosomes
ls ${OUT}/06_final/*.maf002.miss03.vcf.gz > ${OUT}/concat.list
bcftools concat -f ${OUT}/concat.list -Oz -o ${OUT}/final.maf002.miss03.vcf.gz
bcftools index -c -f ${OUT}/final.maf002.miss03.vcf.gz
