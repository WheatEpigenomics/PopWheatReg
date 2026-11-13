#!/usr/bin/env bash
# Runtime: all-pair TensorQTL scans, GPU strongly recommended; production scale is in README.md.
# caQTL / eQTL association runs
#
# Run this script from the directory that contains it and the
# TensorQTL_subgroup_omnibus/ folder.
#
# All phenotype matrices are inverse-normal transformed (INT) beforehand:
#   rows = phenotypes, columns = samples, values already INT.
#
#   eQTL : global_pc         y ~ G + PC1..PC3                 once
#   caQTL: global_pc         y ~ G + PC1..PC3                 once
#   caQTL: subgroup_omnibus  per-subgroup y ~ G, no PC         once
#           (same INT caQTL phenotype matrix as the caQTL global run)

PY=python
CODE=./TensorQTL_subgroup_omnibus

# ---- inputs ----
# INT-transformed phenotype matrices: phenotypes x samples.
# Accepted suffixes: .txt, .tsv, .parquet (format inferred from suffix).
PHENO_EQTL_INT=/path/to/eqtl_phenotype_INT.txt
PHENO_CAQTL_INT=/path/to/caqtl_phenotype_INT.txt

GENO_EQTL=/path/to/eqtl_genotype
GENO_CAQTL=/path/to/caqtl_genotype

PC3=/path/to/pc1_3.tsv
SUBGROUPS=/path/to/subgroups.tsv

OUT=./output
mkdir -p "${OUT}"


# ---- eQTL: global + 3 PCs  ----
# Input: INT-transformed eQTL phenotype matrix + PC1..PC3 covariates.
"${PY}" "${CODE}/00_TensorQTL_add_subgroup_asso.py" \
  --mode global_pc \
  --genotype  "${GENO_EQTL}" \
  --phenotype "${PHENO_EQTL_INT}" \
  --covariates "${PC3}" \
  --maf_threshold 0.05 \
  --pval_threshold 1e-5 \
  --batch_size 20000 \
  --out "${OUT}/eqtl.global_pc.tsv"


# ---- caQTL: global + 3 PCs ----
# Input: INT-transformed caQTL phenotype matrix + PC1..PC3 covariates.
"${PY}" "${CODE}/00_TensorQTL_add_subgroup_asso.py" \
  --mode global_pc \
  --genotype  "${GENO_CAQTL}" \
  --phenotype "${PHENO_CAQTL_INT}" \
  --covariates "${PC3}" \
  --maf_threshold 0.05 \
  --pval_threshold 1e-5 \
  --batch_size 20000 \
  --out "${OUT}/caqtl.global_pc.tsv"


# ---- caQTL: subgroup_omnibus (no PC) ----
# Input: the SAME INT-transformed caQTL phenotype matrix used above.
# Model: within each subgroup y ~ G; no covariates are passed.
"${PY}" "${CODE}/00_TensorQTL_add_subgroup_asso.py" \
  --mode subgroup_omnibus \
  --genotype  "${GENO_CAQTL}" \
  --phenotype "${PHENO_CAQTL_INT}" \
  --subgroups "${SUBGROUPS}" \
  --maf_threshold 0.05 \
  --pval_threshold 1e-5 \
  --batch_size 20000 \
  --pheno_chunk 2000 \
  --out "${OUT}/caqtl.sgom.parquet"

"${PY}" "${CODE}/03_merge_subgroup_parts.py" --out "${OUT}/caqtl.sgom.parquet"
