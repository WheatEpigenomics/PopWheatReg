# Phenotype association

GWAS with EMMAX, and caWAS / TWAS with PrediXcan. One script, three independent
blocks.

## `01_GWAS_caWAS_TWAS_main.sh`

| Block | Tool | Command |
|---|---|---|
| GWAS | EMMAX | `emmax-kin-intel64` for the kinship matrix, then `emmax-intel64` per phenotype |
| caWAS | PrediXcan | `PrediXcan.py --assoc --pred_exp <ACR matrix>` |
| TWAS | PrediXcan | `PrediXcan.py --assoc --pred_exp <expression matrix>` |

Inputs: PLINK tped/tfam genotype files, a phenotype table, a kinship matrix
(`*.BN.kinf`), and the predicted ACR accessibility and gene expression matrices.

PrediXcan is Python 2 in the production environment.

## Demo

No demo is shipped: a run needs EMMAX and PrediXcan plus full-size genotype and
phenotype data.

## Requirements and runtime

EMMAX, PrediXcan (Python 2).

A demo is not shipped, since a run needs EMMAX and PrediXcan plus full-size
genotype and phenotype data. Production scale is 203 accessions; cost grows
linearly with the number of phenotypes.
