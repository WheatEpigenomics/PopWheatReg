#### gatk call gvcf
# GATK HaplotypeCaller call snp
mkdir ./tmp
gatk --java-options "-Xmx200G -Djava.io.tmpdir=./tmp" HaplotypeCaller -R Ta_MC.fa -I ${name}.clean.bam -O ${name}.gvcf -ERC GVCF -stand-call-conf 30

#### Merge gvcf and convert it to vcf
#write all gvcf to gvcf.list
path=
find $path -name *gvcf > gvcf.list
#combine gvcf
gatk --java-options -Xmx600G CombineGVCFs -R Ta_MC.fa -V gvcf.list -O hebing.gvcf 
#detect SNP (trans gvcf to vcf)
gatk --java-options -Xmx600G GenotypeGVCFs -R Ta_MC.fa -V hebing.gvcf -O hebing.vcf

#### vcf filter
#extract SNP
gatk --java-options -Xmx600G SelectVariants -R Ta_MC.fa -O hebing.SNPs.vcf --variant hebing.vcf --select-type-to-include SNP --create-output-variant-index true
#hard filter
gatk VariantFiltration -V hebing.SNPs.vcf -O hebing.SNPs.mark.vcf --filter-expression 'QD < 2.0 || FS > 60.0 || SOR > 3.0 || MQ < 40.0' --filter-name lowQualFilter --missing-values-evaluate-as-failing
#Retaining biallelic SNPs
vcftools --vcf hebing.SNPs.mark.vcf --min-alleles 2 --max-alleles 2 --recode --recode-INFO-all --out hebing.SNPs.biall
#Quality filtering based on sequencing depth and heterozygosity
grep -v "lowQualFilter"  hebing.SNPs.biall.recode.vcf | bcftools filter -S . -e 'FORMAT/DP < 3' | bcftools filter -S . -e 'FORMAT/DP > 1000' | bcftools +setGT -- -t q -i 'FMT/GT ~ "0[|/]1"' -n . | bcftools +setGT -- -t q -i 'FMT/GT ~ "1[|/]0"' -n . > snp.fil.vcf
