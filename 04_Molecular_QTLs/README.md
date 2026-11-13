# Molecular QTLs

caQTL and eQTL association with TensorQTL, all pair, no cis filtering.

## Layout

```
04_Molecular_QTLs/
├── 01_caQTL_eQTL_main.sh          the three association runs plus the merge
└── TensorQTL_subgroup_omnibus/    all-pair association code and its own README
```

## Runs in `01_caQTL_eQTL_main.sh`

| Run | Model | Notes |
|---|---|---|
| eQTL `global_pc` | `y ~ G + PC1..PC3` | every SNP x gene pair |
| caQTL `global_pc` | `y ~ G + PC1..PC3` | every SNP x ACR pair |
| caQTL `subgroup_omnibus` | per subgroup `y ~ G`, no PCs, then an exact omnibus across subgroups | same INT caQTL matrix |

The third run writes one parquet per subgroup and needs
`TensorQTL_subgroup_omnibus/03_merge_subgroup_parts.py` afterwards.

Phenotypes are inverse-normal transformed beforehand: rows = phenotypes, columns =
samples. Genotypes are plink bed/bim/fam.

## `subgroup_omnibus` in one paragraph

Each subgroup is fitted separately with its own sample size and its own MAF and
imputation rules; a subgroup that is not estimable is skipped. Per-subgroup
`P_k` is turned into `Z_k`, `W = sum Z_k^2` is compared with a chi-square on the
estimated number of subgroups, and the minimum subgroup P is reported as well.
A safe GPU prefilter skips pairs that cannot reach either threshold without
changing the statistics.

## Demo

`TensorQTL_subgroup_omnibus/demo/` builds a small synthetic dataset and runs
`global_pc`, `subgroup_omnibus` and the merge in a few seconds.

## Requirements and runtime

TensorQTL 1.0.10 with PyTorch, pandas, numpy and scipy; a GPU speeds up the
all-pair scan but the code also runs on CPU.

Production scale is 202 accessions for eQTL and 200 for caQTL over 21 chromosomes.
The demo runs both all-pair scans on a synthetic dataset in under a minute on CPU;
a full run scales with the number of SNP x phenotype pairs.
