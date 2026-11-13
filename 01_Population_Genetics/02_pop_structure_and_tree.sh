# Runtime: ADMIXTURE for K=1..10 plus tree building; production scale is in README.md.
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
fastme -i output.mat  -o output.nwk -m O

