#### Population structure
#calculation
for i in {1..10}
do 
 admixture --cv input.bed $i -j24 | tee log${i}.out
done 
#find minimum CV
grep -h 'CV'  log*.out


#### Phylogenetic tree
VCF2Dis -i input.vcf -o output.mat
fastme -i 210snp_miss02_maf02.mat  -o 210snp_miss02_maf02_OLSME.nwk -m O

