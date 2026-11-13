# Population genetics

Variant calling, population structure, phylogenetic tree, LD decay and haplotype
blocks. Three independent scripts; run them in that order.

## `01_call_SNP_and_filter.sh`

GATK GVCF calling followed by the same hard-filter flow used for the production
csv2 callset.

```
HaplotypeCaller per sample (GVCF, stand-call-conf 30)
  -> CombineGVCFs -> GenotypeGVCFs
  -> per chromosome:
       VariantFiltration   QD < 2 || FS > 60 || SOR > 3 || MQ < 40
       SelectVariants      keep PASS only
       bcftools view       biallelic SNPs
       bcftools +setGT     DP < 3, DP > 1000 and heterozygous calls -> missing
       bcftools annotate   ID = CHROM-POS, tags recomputed
       bcftools view       AN > 0 && AC > 0 && MAF >= 0.02 && F_MISSING <= 0.3
  -> bcftools concat
```

Inputs: per-sample BAM files and the reference FASTA. Output: one VCF per
chromosome plus the concatenated `final.maf002.miss03.vcf.gz`.

Sample subsetting (210 all / 203 hexaploid) happens at the tagged step; the
commented line shows where.

## `02_pop_structure_and_tree.sh`

ADMIXTURE for K = 1..10 with cross-validation, then the minimum-CV K is chosen
from `log*.out`. The tree is built from the VCF with VCF2Dis (distance matrix)
and FastME (NJ, `-m O`).

## `03_LD_decay_and_genome_blocks.sh`

LD decay and PLINK haplotype blocks per subgenome. The wheat csv2 section splits
A/B/D, thins SNPs with a fixed seed, computes pairwise r2 within subgenome
(A/B within 5 Mb, D within 30 Mb), bins the decay curve with `ld_decay.py`, and
runs PLINK `--blocks` with the V1 parameters.

`ld_decay.py` sits next to the script; it prints the distance where mean r2 drops
below 0.2, which sets `--blocks-max-kb`.

Production parameters for wheat csv2: thin 0.05 / 0.05 / 0.20, LD windows
5 / 5 / 30 Mb, decay bins 10 / 10 / 20 kb, blocks-max-kb 9000 / 10000 / 15000,
strong-lowci 0.62, strong-highci 0.92, recomb-highci 0.92, inform-frac 0.80.

## Requirements and runtime

GATK 4.1, bcftools, PLINK 1.9, ADMIXTURE, VCF2Dis, FastME, Python 3 with numpy and
matplotlib.

No demo is shipped: a meaningful run needs these tools and full-size inputs.
Production scale is 210 accessions over 21 chromosomes. Memory is the limiting
factor. The calling step needs 30 GB for HaplotypeCaller,
600 GB heap for CombineGVCFs and 100 GB for GenotypeGVCFs; the per-chromosome
filter steps use 40 GB.
