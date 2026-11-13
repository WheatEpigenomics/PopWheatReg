# TensorQTL subgroup-omnibus workflow

A small, standalone layer on top of **TensorQTL 1.0.10** for caQTL / eQTL
association. It adds two association modes that TensorQTL does not provide
(`global_pc`, `subgroup_omnibus`), plus two convenience modes (`cis_all`,
`selected`) and a merge helper.

It never modifies TensorQTL and never imports the non-upstream
`tensorqtl_NOfilcis.py`; it only calls `tensorqtl.core` / `tensorqtl.genotypeio`
building blocks (`Residualizer`, `calculate_corr`, `filter_maf`, `impute_mean`,
`get_t_pval`, `GenotypeGeneratorTrans`, `load_genotypes`).

All-pair modes (`global_pc`, `subgroup_omnibus`, `selected`) deliberately do **no
cis filtering**: every tested SNP × phenotype pair is a candidate. cis/trans
annotation is left to downstream scripts.

---

## Requirements

- Python ≥ 3.8
- `tensorqtl==1.0.10` (the association math is tied to this release)
- `torch` (CPU is fine; GPU is used when available), `scipy`, `numpy`, `pandas`, `pyarrow`

```bash
pip install tensorqtl==1.0.10
```

Optional GPU: any CUDA build of torch. If `torch.cuda.is_available()` is False the
scripts run on CPU automatically.

---

## Files

| File | Role |
|---|---|
| `00_TensorQTL_add_subgroup_asso.py` | single entry point; dispatches the four modes |
| `01_cis_all_mode_run.py` | implementation of `--mode cis_all` (loaded internally) |
| `02_selected_mode_run.py` | implementation of `--mode selected` (loaded internally) |
| `03_merge_subgroup_parts.py` | merges the streaming `--mode subgroup_omnibus` parts |
| `demo/make_demo_data.py` | generates a tiny synthetic dataset |
| `demo/run_demo.sh` | end-to-end demo (global_pc + subgroup_omnibus + merge) |

`00`, `01`, `02`, `03` **must stay in the same directory**: `00` loads `01`/`02`
by path, and `02` loads `01` by path.

---

## Quick start (synthetic demo, ~10 s)

```bash
bash demo/run_demo.sh
```

The demo writes `demo/data/` (genotypes, INT phenotypes, 3 PCs, subgroup labels)
and `demo/out/` (one `global_pc` table and one `subgroup_omnibus` parquet), then
prints the first rows. Nothing outside `demo/` is touched.

---

## Input formats

### Genotypes (`--genotype`)
- PLINK 1 prefix (`*.bed/*.bim/*.fam`)
- PLINK 2 prefix (`*.pgen/*.psam/*.pvar`)
- parquet matrix, `variants × samples` (index = variant id)
- gzipped TSV matrix, `variants × samples`

Missing genotypes must be `-9`; they are mean-imputed **within the analysis
samples** (within each subgroup in `subgroup_omnibus`).

### Phenotypes (`--phenotype`)
`phenotypes × samples`, index = phenotype id, columns = sample id:
- BED (`.bed`, `.bed.gz`) — for cis modes, or a plain phenotype BED matrix
- parquet
- TSV (`index_col=0`)

Values are expected to be already inverse-normal transformed (INT).

### Covariates (`--covariates`)
`covariates × samples` (this file is transposed internally to `samples ×
covariates`). Optional for `global_pc`; **omit it** for `subgroup_omnibus`
(which is a no-PC, within-group model).

### Subgroups (`--subgroups`)
Headerless two-column TSV: `sample_id <TAB> subgroup`. Every analysis sample
must appear exactly once. Subgroups must be disjoint; each needs ≥ 3 samples.

### cis positions (`--phenotype_positions`)
Only for `cis_all` / `selected` with matrix genotypes:
- eQTL: BED6 with strand (TSS = `start+1` on `+`, `end` on `-`)
- caQTL: BED4 ACR intervals (center = `floor((start+end)/2)+1`); alternatively
  ACR IDs that already encode `chr_start_end`

Coordinates are **BED 0-based half-open**; SNP positions are 1-based.

For matrix genotypes in cis/selected modes, `--variant_metadata` is required:
`variant_id, chrom, pos, effect_allele, other_allele[, ref, alt]`.

---

## Modes and templates

All examples assume `PY=python` and the scripts in the current directory.

### 1. `global_pc` — all-pair `y ~ G + covariates`

```bash
python 00_TensorQTL_add_subgroup_asso.py \
  --mode global_pc \
  --genotype   GENO \
  --phenotype  PHENO_INT.tsv.gz \
  --covariates PC1_3.tsv \
  --maf_threshold 0.05 --pval_threshold 1e-5 \
  --batch_size 20000 \
  --out OUT/global_pc.tsv
```

- Omit `--covariates` to run the no-PC model `y ~ G` (`dof = n - 2`).
- Output columns: `variant_id, phenotype_id, pval, beta, beta_se, af`.
- Only pairs with `pval < pval_threshold` are emitted.
- Identical to TensorQTL trans (non-interaction) **without** cis filtering.

### 2. `subgroup_omnibus` — per-subgroup association + exact omnibus

```bash
python 00_TensorQTL_add_subgroup_asso.py \
  --mode subgroup_omnibus \
  --genotype   GENO \
  --phenotype  PHENO_INT.tsv.gz \
  --subgroups  SUBGROUPS.tsv \
  --maf_threshold 0.05 \
  --pval_threshold 1e-5 \
  --batch_size 20000 --pheno_chunk 2000 \
  --out OUT/sgom.parquet

# REQUIRED second step: merge the verified streaming parts
python 03_merge_subgroup_parts.py --out OUT/sgom.parquet
```

- `--out` must end in `.parquet`; its parent directory must already exist.
  The run writes `<out>.tmp_parts/` and the merge step assembles + deletes them.
- Separate thresholds: `--omnibus_p_threshold` and `--pmin_threshold`
  (both default to `--pval_threshold`).
- `--pval_threshold 1` = output-all mode: every `K_est ≥ 1` pair is emitted
  using the exact path, ignoring the PW/Pmin gates and the prefilter.
- `--no_prefilter` disables the safe GPU prefilter (slower, numerically identical).
- `--resume` continues an interrupted run from its tile markers.
- Output columns: `variant_id, phenotype_id, P_W, W, K_est, P_min,
  strongest_subgroup`, then `{group}_beta/_se/_p/_maf` per subgroup.

### 3. `cis_all` — nominal cis scan (all computable cis pairs)

```bash
# eQTL (gene TSS)
python 00_TensorQTL_add_subgroup_asso.py \
  --mode cis_all --phenotype_kind eqtl \
  --genotype GENO --phenotype PHENO_INT.tsv.gz \
  --phenotype_positions GENES_BED6.bed \
  --window 1000000 \
  --out OUT/cis_eqtl          # must be an EXISTING EMPTY directory

# caQTL (ACR midpoint)
python 00_TensorQTL_add_subgroup_asso.py \
  --mode cis_all --phenotype_kind caqtl \
  --genotype GENO --phenotype ACR_INT.tsv.gz \
  --phenotype_positions ACR_BED4.bed \
  --window 1000000 \
  --out OUT/cis_caqtl
```

- Writes one `NN_<chrom>.cis_all.tsv.gz` per chromosome, a phenotype QC table,
  a legacy 10-column nominal table, and a completion manifest.
- `--pval_threshold` is ignored here (all computable pairs are exported).

### 4. `selected` — all-pair association for selected genes/ACRs

```bash
python 00_TensorQTL_add_subgroup_asso.py \
  --mode selected --method global \
  --genotype GENO --phenotype PHENO_INT.tsv.gz \
  --covariates PC1_3.tsv \
  --gene-list GENES.txt \
  --pval_threshold 1 \
  --out OUT/batch_global          # new batch DIRECTORY

python 00_TensorQTL_add_subgroup_asso.py \
  --mode selected --method pw \
  --genotype GENO --phenotype PHENO_INT.tsv.gz \
  --subgroups SUBGROUPS.tsv \
  --acr-list ACRS.txt \
  --pval_threshold 1 \
  --out OUT/batch_pw
```

- `--method` ∈ `global | pw | sg1 | sg2 | sg3 | sg5`.
- `--pval_threshold` **must be exactly 1** (this mode exports all associations).
- `global` accepts `--covariates`; `pw`/`sg*` are within-group no-PC models and
  must omit covariates; `pw` requires all four groups `sg1,sg2,sg3,sg5`.
- Selectors: `--acr-id/--acr-list` or `--gene-id/--gene-list`.
- `--out` is a new/empty batch directory with `f01_input … f05_logs/`.
- Optional Manhattan plots: `--plot`, `--plot-only`, `--plot-region`.

---

## Method summary

**`global_pc`** (identical to TensorQTL trans, no cis filter): residualize `G` and
`Y` on the covariates, compute the Pearson correlation `r`, then

```
t      = r * sqrt(dof / (1 - r²))
beta   = r * sqrt(var(Y_res) / var(G_res))
beta_se= beta / t
pval   = 2 * t.sf(|t|, dof),    dof = n - 2 - n_covariates
```

**`subgroup_omnibus`** — per subgroup `k` with its own `n_k`, `dof_k = n_k - 2`,
`MAF_k`, imputation and estimability mask:

```
t_k   = r_k * sqrt(dof_k / (1 - r_k²))
P_k   = betainc(dof_k/2, 0.5, dof_k/(dof_k + t_k²))     # == 2 * t.sf(|t_k|, dof_k)
Z_k   = -ndtri(P_k / 2)                                 # == norm.isf(P_k/2), tail-safe
W     = Σ_k Z_k²                (over estimable subgroups only)
P_W   = gammaincc(K_est/2, W/2)                         # == chi2.sf(W, df=K_est)
P_min = min_k P_k ;  strongest_subgroup = argmin_k P_k
PASS  = (P_W < omnibus_p_threshold) OR (P_min < pmin_threshold)
```

`K_est` is computed **per variant × phenotype pair** (1..K). A subgroup is not
estimable if the variant is monomorphic there, its MAF fails, `gvar` is
non-finite, or the statistic is `NaN`.

**Numerical edge cases.** A perfect fit `|r_k| = 1` is a valid extreme result
(`SE = 0`, `P_k = 0`, `Z_k = +∞`, `W = +∞`, `P_W = 0` → PASS), and so is a
finite-`t` `P_k` underflow to 0. `|r_k| > 1 + 1e-6` and `NaN` are invalid. A
near-perfect but strictly `|r| < 1` correlation is **not** snapped to 1; it is
computed as measured and only counted (`near_perfect_pairs` in the log).

**Safe GPU prefilter (default ON).** Because the exact t-CDF runs on the CPU, the
implementation first selects candidates on the GPU with a provably necessary
condition:

```
zgate_W   = min over m = 1..Kmax of sqrt(chi2.isf(omnibus_p_threshold, m) / m)
zgate_min = norm.isf(pmin_threshold / 2)
zgate     = min(zgate_W, zgate_min)
r_gate_k  = (t.isf(p_gate/2, dof_k)) / sqrt(dof_k + t.isf(p_gate/2, dof_k)²),
            p_gate = 2 * norm.sf(zgate)
candidate |= estimable_k & (|r_k| >= r_gate_k)
```

The exact transform is then applied only to candidate pairs. This is a necessary
condition for either PASS branch, so **ON and OFF produce identical results**;
`--no_prefilter` is the reference path.

---

## Notes and caveats

- **No cis filtering** in `global_pc` / `subgroup_omnibus` / `selected`.
- **Disjoint subgroups.** `Σ Z_k² ~ χ²(K)` assumes independent subgroups.
  Overlap is warned about but not forbidden.
- **Exact t-CDF is CPU-only.** `torch` 2.x has no `betainc` and
  `StudentT.cdf` is not implemented, so the exact transform uses
  `scipy.special` (`betainc`, `ndtri`, `gammaincc`) in float64. GPU is still used
  for the correlation matrix multiplication.
- **`subgroup_omnibus` is streaming.** A run writes `<out>.tmp_parts/`; the final
  parquet appears only after `03_merge_subgroup_parts.py`.
- **Existing outputs are protected.** Modes refuse to overwrite existing results
  (`subgroup_omnibus` needs a new `--out`, `cis_all` an empty directory,
  `selected` a new batch directory).
- **Memory knobs:** `--batch_size` (variants per batch) and `--pheno_chunk`
  (phenotypes per tile).

---

## License / citation

The association math is derived from TensorQTL (BSD 3-Clause). Please cite the
TensorQTL paper together with the study that uses this workflow. Add a
`LICENSE` file (BSD-3-Clause is the natural choice) before publishing the
repository.
