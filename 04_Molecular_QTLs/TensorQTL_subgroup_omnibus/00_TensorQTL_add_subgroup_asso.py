#!/usr/bin/env python3
"""
Standalone all-pair caQTL association built on the *installed* TensorQTL 1.0.10.

================================================================================
IMPORTANT
================================================================================
This script does NOT modify site-packages, does NOT copy the TensorQTL package,
and does NOT import the local non-upstream file `tensorqtl_NOfilcis.py`.
It only imports low-level building blocks from the installed TensorQTL:

    tensorqtl.core.Residualizer
    tensorqtl.core.calculate_corr
    tensorqtl.core.filter_maf
    tensorqtl.core.impute_mean
    tensorqtl.core.get_t_pval
    tensorqtl.core.SimpleLogger
    tensorqtl.genotypeio.GenotypeGeneratorTrans
    tensorqtl.genotypeio.load_genotypes

The association math (correlation -> t -> beta -> SE -> P) is identical to
`tensorqtl.trans.map_trans` (non-interaction branch). See README.md for the method
description, input/output schemas and usage templates.

Association modes
---------
  global_pc:
      y_int ~ G + covariates  (identical to TensorQTL trans non-interaction)
      all SNP x all ACR pairs, NO cis filtering.
      outputs (within pval_threshold): variant_id, phenotype_id, pval, beta, beta_se, af

  subgroup_omnibus:
      For each subgroup k (disjoint sample strata), independently:
          y_int ~ G          (no PC in this stage)
      using each subgroup's own n_k, dof_k = n_k - 2, MAF_k, per-subgroup missing
      imputation, and a per-variant "estimable" mask.

      Per estimable subgroup: t_k, beta_k, SE_k, P_k.
      P_min = min_k(P_k), strongest_subgroup = argmin.
      Exact omnibus:
          P_k = 2 * t.sf(|t_k|, dof_k)          (== betainc(dof_k/2, 0.5, dof_k/(dof_k+t_k^2)))
          Z_k = norm.isf(P_k / 2)               (== -ndtri(P_k/2); direction not used)
          W   = sum_k Z_k^2
          P_W = chi2.sf(W, df=K_est)            (== gammaincc(K_est/2, W/2))
      K_est is computed per variant x phenotype pair (number of estimable subgroups).

      PASS = (P_W < omnibus_p_threshold) OR (P_min < pmin_threshold)
      (both thresholds default to --pval_threshold).
      For pval_threshold < 1, only PASS pairs are emitted. At
      pval_threshold == 1, output-all mode uses the exact path and emits every
      K_est >= 1 pair with every observed subgroup statistic.

      SAFE GPU PRE-FILTER (default ON; disable with --no_prefilter):
      the exact finite-df transform (betainc/ndtri/gammaincc) is expensive, so
      candidates are first selected on the GPU with a provably necessary condition
      on |r_k| (see _prefilter_zgate). Any pair that can pass either branch is kept;
      the exact statistics and the final PASS set are IDENTICAL with and without the
      pre-filter. --no_prefilter runs the exact transform on every pair and is the
      numerical golden reference.

See README.md for input/output schemas and examples.

  cis_all:
      Independent nominal cis scan around strand-aware gene TSS or ACR midpoint.
      Outputs all computable pairs, including nonsignificant and zero correlation.
      See README.md for this mode's template; it does not change the two existing
      implementations.

  selected:
      Select gene/ACR IDs and run global, PW or one subgroup. All estimable
      associations, text and Manhattan inputs; plotting is opt-in. See README.md.
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import sys
import time
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

import numpy as np
import pandas as pd
import torch
import scipy.special as sp
import scipy.stats as stats

from tensorqtl.core import (
    Residualizer,
    calculate_corr,
    filter_maf,
    impute_mean,
    get_t_pval,
    read_phenotype_bed,
    SimpleLogger,
)
from tensorqtl.genotypeio import GenotypeGeneratorTrans, load_genotypes



def _preview(values, limit=5):
    values = list(values)
    shown = ", ".join(map(str, values[:limit]))
    return shown + (f", ... ({len(values)} total)" if len(values) > limit else "")


def _check_unique(index, name):
    index = pd.Index(index)
    duplicates = index[index.duplicated(keep=False)].unique().tolist()
    if duplicates:
        raise ValueError(f"Duplicate {name}: {len(duplicates)} duplicate IDs; "
                         f"examples: {_preview(duplicates)}")


def validate_parameters(pval_threshold, maf_threshold, batch_size,
                        pheno_chunk=None, omnibus_p_threshold=None, pmin_threshold=None):
    if not (0 < pval_threshold <= 1):
        raise ValueError(f"pval_threshold must satisfy 0 < pval_threshold <= 1; got {pval_threshold!r}.")
    for name, value in (("omnibus_p_threshold", omnibus_p_threshold),
                        ("pmin_threshold", pmin_threshold)):
        if value is not None and not (0 < value < 1):
            raise ValueError(f"{name} must satisfy 0 < {name} < 1; got {value!r}.")
    if not (0 <= maf_threshold < 0.5):
        raise ValueError(f"maf_threshold must satisfy 0 <= maf_threshold < 0.5; got {maf_threshold!r}.")
    if batch_size <= 0:
        raise ValueError(f"batch_size must be > 0; got {batch_size!r}.")
    if pheno_chunk is not None and pheno_chunk <= 0:
        raise ValueError(f"pheno_chunk must be > 0; got {pheno_chunk!r}.")


def validate_matrices(genotype_df, phenotype_df, covariates_df=None):
    _check_unique(genotype_df.columns, "genotype sample IDs")
    _check_unique(phenotype_df.columns, "phenotype sample IDs")
    _check_unique(phenotype_df.index, "phenotype IDs")
    if covariates_df is not None:
        _check_unique(covariates_df.index, "covariate sample IDs")
    for name, frame in (("Phenotype", phenotype_df), ("Covariate", covariates_df)):
        if frame is not None:
            try:
                finite = np.isfinite(frame.to_numpy(dtype=np.float64))
            except (TypeError, ValueError) as exc:
                raise ValueError(f"{name} matrix must contain numeric, finite values.") from exc
            if not finite.all():
                n_bad = int((~finite).sum())
                raise ValueError(f"{name} matrix contains {n_bad} NaN or infinite values; "
                                 "check the input file before running.")


def validate_subgroups(subgroup_s, analysis_samples):
    _check_unique(subgroup_s.index, "subgroup sample IDs")
    analysis_samples = pd.Index(analysis_samples)
    missing = analysis_samples.difference(subgroup_s.index, sort=False)
    if len(missing):
        raise ValueError(f"Missing subgroup assignment for {len(missing)} analysis samples; "
                         f"examples: {_preview(missing)}")
    assigned = subgroup_s.reindex(analysis_samples)
    invalid = assigned.isna() | assigned.astype(str).str.strip().eq("")
    if invalid.any():
        examples = [f"{sample}={label!r}" for sample, label in assigned[invalid].iloc[:5].items()]
        raise ValueError(f"Empty subgroup label for {int(invalid.sum())} analysis samples; "
                         f"examples: {_preview(examples)}")
    assigned = assigned.astype(str)
    labels = sorted(assigned.unique())
    counts = {lab: int((assigned == lab).sum()) for lab in labels}
    for lab, n in counts.items():
        if n - 2 <= 0:  # current subgroup model has q_k = 0
            raise ValueError(f"{lab} has {n} analysis samples; dof_k = n_k - 2 - q_k "
                             f"= {n - 2} (q_k=0) must be > 0.")
    return counts


def validate_output_parent(path):
    parent = os.path.dirname(os.path.abspath(path))
    if not os.path.isdir(parent):
        raise ValueError(f"Output parent directory does not exist: {parent}")


def check_raw_sample_header(path, offset, name):
    """Catch duplicate TSV/BED sample names before pandas mangles them."""
    opener = gzip.open if path.endswith(".gz") else open
    with opener(path, "rt") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
    if len(header) <= offset:
        raise ValueError(f"{name} file has no sample columns: {path}")
    _check_unique(header[offset:], f"{name} sample IDs")


def plink_sample_ids(path):
    """Read PLINK sample IDs before loading genotypes, preserving sample subsetting."""
    if all(os.path.exists(f"{path}.{ext}") for ext in ("bed", "bim", "fam")):
        ids = pd.read_csv(f"{path}.fam", sep=r"\s+", header=None, dtype=str, usecols=[1]).iloc[:, 0]
        return pd.Index(ids, name="sample_id")
    if all(os.path.exists(f"{path}.{ext}") for ext in ("pgen", "psam", "pvar")):
        table = pd.read_csv(f"{path}.psam", sep=r"\s+", dtype=str)
        column = next((c for c in ("IID", "#IID") if c in table.columns), None)
        if column is None:
            raise ValueError(f"PLINK2 sample file has no IID/#IID column: {path}.psam")
        return pd.Index(table[column], name="sample_id")
    return None


def validate_alignment(genotype_df, phenotype_df, covariates_df=None):
    common = phenotype_df.columns.intersection(genotype_df.columns, sort=False)
    if len(common) == 0:
        raise ValueError("Genotype and phenotype have 0 common sample IDs; check the input files.")
    missing_g = phenotype_df.columns.difference(genotype_df.columns, sort=False)
    if len(missing_g):
        raise ValueError(f"{len(missing_g)} phenotype samples are missing from the genotype matrix; "
                         f"examples: {_preview(missing_g)}")
    if covariates_df is not None:
        missing_c = phenotype_df.columns.difference(covariates_df.index, sort=False)
        if len(missing_c):
            raise ValueError(f"{len(missing_c)} analysis samples are missing from covariates; "
                             f"examples: {_preview(missing_c)}")
    return common


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------
def _get_device():
    return torch.device("cuda" if torch.cuda.is_available() else "cpu")


def load_phenotypes(path):
    """Load phenotype matrix (phenotypes x samples).

    Supported: BED (.bed / .bed.gz, positions dropped), parquet, tab-delimited
    TSV with index_col=0.  Index = phenotype_id, columns = sample_id.
    """
    if path.lower().endswith((".bed", ".bed.gz")):
        check_raw_sample_header(path, 4, "phenotype")
        phenotype_df, _ = read_phenotype_bed(path)
    elif path.endswith(".parquet"):
        phenotype_df = pd.read_parquet(path)
    else:
        check_raw_sample_header(path, 1, "phenotype")
        phenotype_df = pd.read_csv(path, sep="\t", index_col=0)
    return phenotype_df


def load_covariates(path):
    """Load covariates. File layout follows TensorQTL CLI: covariates x samples.

    Returns a DataFrame with index = sample_id and columns = covariate names
    (i.e. transposed relative to the input file).
    """
    check_raw_sample_header(path, 1, "covariate")
    covariates_df = pd.read_csv(path, sep="\t", index_col=0).T
    return covariates_df


def load_subgroups(path):
    """Load subgroup mapping: headerless TSV with two columns: sample_id, subgroup.

    Returns pd.Series with index = sample_id, values = subgroup label.
    """
    df = pd.read_csv(path, sep="\t", header=None, comment="#",
                     dtype=str, keep_default_na=False)
    if df.shape[1] < 2:
        raise ValueError(
            f"Subgroup file must have >= 2 columns (sample_id, subgroup): {path}")
    sample_ids = df[0].astype(str)
    labels = df[1].astype(str)
    repeated = sample_ids[sample_ids.duplicated(keep=False)].unique().tolist()
    if repeated:
        conflicts = []
        for sample in repeated:
            found = labels[sample_ids == sample].unique().tolist()
            if len(found) > 1:
                conflicts.append(f"{sample}: {', '.join(found)}")
        detail = f"; conflicting labels: {_preview(conflicts)}" if conflicts else ""
        raise ValueError(f"Duplicate subgroup sample IDs: {len(repeated)} duplicate samples; "
                         f"examples: {_preview(repeated)}{detail}")
    return pd.Series(labels.values, index=sample_ids.values, name="subgroup")


# ---------------------------------------------------------------------------
# MODE 1: global_pc
# ---------------------------------------------------------------------------
def run_global_pc(genotype_df, phenotype_df, covariates_df=None,
                  pval_threshold=1e-5, maf_threshold=0.05, batch_size=20000,
                  pheno_chunk=None, logger=None, verbose=True):
    """All-pair association, y ~ g + covariates. No cis filtering.

    Mirrors `tensorqtl.trans.map_trans` non-interaction sparse branch exactly
    (same residualization, same correlation shortcut, same dof, same formulas),
    except that (a) it never applies filter_cis and (b) it renames columns
    b -> beta, b_se -> beta_se and drops r2. Optional phenotype chunks
    bound device memory without changing the association model.

    Returns DataFrame columns: variant_id, phenotype_id, pval, beta, beta_se, af
    """
    validate_parameters(pval_threshold, maf_threshold, batch_size, pheno_chunk)
    validate_matrices(genotype_df, phenotype_df, covariates_df)
    validate_alignment(genotype_df, phenotype_df, covariates_df)
    n_samples = phenotype_df.shape[1]
    dof = n_samples - 2 - (0 if covariates_df is None else covariates_df.shape[1])
    if dof <= 0:
        raise ValueError(f"global_pc needs positive dof: n_samples={n_samples}, "
                         f"n_covariates={0 if covariates_df is None else covariates_df.shape[1]}, "
                         f"dof={dof}.")
    device = _get_device()
    if logger is None:
        logger = SimpleLogger(verbose=verbose)

    if covariates_df is not None:
        if not phenotype_df.columns.equals(covariates_df.index):
            raise ValueError("Covariate sample order must match phenotype sample order; "
                             "reindex covariates to the analysis samples before running.")
        residualizer = Residualizer(
            torch.tensor(covariates_df.values, dtype=torch.float32).to(device))
    else:
        residualizer = None

    logger.write("all-pair association (global_pc)")
    logger.write(f"  * {n_samples} samples")
    logger.write(f"  * {phenotype_df.shape[0]} phenotypes")
    logger.write(f"  * {0 if covariates_df is None else covariates_df.shape[1]} covariates")
    logger.write(f"  * dof = {dof}")
    logger.write(f"  * pval_threshold = {pval_threshold:.3g}, maf_threshold = {maf_threshold:.3g}")

    n_phenotypes = phenotype_df.shape[0]
    if pheno_chunk is None:
        chunks = [(0, n_phenotypes)]
        phenotypes_t = torch.tensor(phenotype_df.values, dtype=torch.float32).to(device)
        phenotype_values = None
    else:
        chunks = [(i, min(i + pheno_chunk, n_phenotypes))
                  for i in range(0, n_phenotypes, pheno_chunk)]
        phenotypes_t = None
        phenotype_values = phenotype_df.values
        logger.write(f"  * phenotype chunk size = {pheno_chunk} ({len(chunks)} chunks)")
    genotype_ix = np.array([genotype_df.columns.tolist().index(i) for i in phenotype_df.columns])
    genotype_ix_t = torch.from_numpy(genotype_ix).to(device)

    # pval threshold -> correlation threshold (identical to trans.py L98-L100)
    tstat_threshold = -stats.t.ppf(pval_threshold / 2, dof)
    r_threshold = tstat_threshold / np.sqrt(dof + tstat_threshold ** 2)

    ggt = GenotypeGeneratorTrans(genotype_df, batch_size=batch_size)
    res = []
    n_variants_maf_pass = 0
    start_time = time.time()
    for k, (genotypes, variant_ids) in enumerate(ggt.generate_data(verbose=verbose), 1):
        genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
        genotypes_t = genotypes_t[:, genotype_ix_t]
        impute_mean(genotypes_t)
        genotypes_t, variant_ids, af_t = filter_maf(genotypes_t, variant_ids, maf_threshold)
        n_variants_maf_pass += genotypes_t.shape[0]
        if genotypes_t.shape[0] == 0:
            del genotypes_t
            continue

        for c0, c1 in chunks:
            if pheno_chunk is None:
                phenotype_chunk_t = phenotypes_t
            else:
                phenotype_chunk_t = torch.tensor(
                    phenotype_values[c0:c1], dtype=torch.float32).to(device)
            r_t, genotype_var_t, phenotype_var_t = calculate_corr(
                genotypes_t, phenotype_chunk_t, residualizer=residualizer, return_var=True)
            if pheno_chunk is not None:
                del phenotype_chunk_t

            m = r_t.abs() >= r_threshold
            ix_t = m.nonzero(as_tuple=False)
            if ix_t.shape[0] == 0:
                del r_t, m, ix_t, genotype_var_t, phenotype_var_t
                continue
            ix = ix_t.cpu().numpy()

            r_t = r_t.masked_select(m).type(torch.float64)
            r2_t = r_t.pow(2)
            tstat_t = r_t * torch.sqrt(dof / (1 - r2_t))
            std_ratio_t = torch.sqrt(phenotype_var_t[ix_t[:, 1]] / genotype_var_t[ix_t[:, 0]])
            b_t = r_t * std_ratio_t
            b_se_t = (b_t / tstat_t).type(torch.float32)

            res.append(np.c_[
                variant_ids[ix[:, 0]], phenotype_df.index[c0:c1][ix[:, 1]],
                tstat_t.cpu(), b_t.cpu(), b_se_t.cpu(),
                af_t[ix_t[:, 0]].cpu()
            ])
            del r_t, m, ix_t, ix, genotype_var_t, phenotype_var_t
        del genotypes_t

    logger.write(f"  * variants passing MAF >= {maf_threshold:.3g} = {n_variants_maf_pass:,} "
                 f"/ {genotype_df.shape[0]:,}")
    logger.write(f"    elapsed time: {(time.time() - start_time) / 60:.2f} min")
    del phenotypes_t, phenotype_values, residualizer

    cols = ["variant_id", "phenotype_id", "pval", "beta", "beta_se", "af"]
    if len(res) == 0:
        return pd.DataFrame(columns=cols)
    res = np.concatenate(res)
    res[:, 2] = get_t_pval(res[:, 2].astype(np.float64), dof, log=False)
    out = pd.DataFrame(res, columns=cols)
    out["pval"] = out["pval"].astype(np.float64)
    out["beta"] = out["beta"].astype(np.float32)
    out["beta_se"] = out["beta_se"].astype(np.float32)
    out["af"] = out["af"].astype(np.float32)
    logger.write("done.")
    return out


# ---------------------------------------------------------------------------
# MODE 2: subgroup_omnibus (exact finite-df statistics + safe GPU prefilter)
# ---------------------------------------------------------------------------
def _grouped_chi2_sf(W, K_est):
    """P_W = chi2.sf(W, K_est), grouped by the (few, integer) K_est values.

    Grouping lets gammaincc receive a *scalar* shape argument (at most Kmax calls)
    instead of an element-wise shape array over the whole tile. Numerically this is
    the same function: chi2.sf(W, k) == gammaincc(k/2, W/2).
    """
    P_W = np.full(W.shape, np.nan, dtype=np.float64)
    for k in np.unique(K_est):
        kk = int(k)
        if kk <= 0:
            continue
        m = (K_est == kk)
        with np.errstate(divide="ignore", invalid="ignore"):
            P_W[m] = sp.gammaincc(0.5 * float(kk), 0.5 * W[m])
    return P_W


# A perfect fit (|r| == 1) gives SE = 0 and is a legitimate extreme result, not an
# error. float32 correlation can also overshoot 1 by a rounding ulp, so any |r|
# within R_TOL ABOVE 1 is treated as the same perfect fit; |r| > 1 + R_TOL is
# clearly illegal and is marked invalid.
#
# R_TOL is used ONLY for that |r| > 1 tolerance. A near-perfect but strictly
# |r| < 1 correlation is NOT snapped to a perfect fit (no artificial SE = 0): it is
# computed normally and merely COUNTED (R_NEAR_PERFECT) so those few pairs can be
# re-run in float64 later if that is ever needed.
R_TOL = 1e-6
R_NEAR_PERFECT = 1.0 - 1e-6


def _exact_stats_block(r_src, est_src, dof_k, pc_broadcast=None):
    """Shared exact transform + accumulation (identical for both paths).

    r_src[k]   : float32 correlation; shape (ng, Pc) when pc_broadcast=Pc, else (n,)
    est_src[k] : bool estimability; shape (ng,) when pc_broadcast=Pc, else (n,)
    dof_k      : (K,) float64 degrees of freedom

    Iteration order over subgroups is fixed (k = 0..K-1), so P_min / strongest tie
    handling is deterministic and identical for both paths.

    Numerical policy:
      * Estimability is decided FIRST by `est` (per-variant: MAF passes, genotype
        variance > 0, finite gvar/maf). Any pair that reaches this point has finite,
        positive gvar/pvar, so beta is always well defined.
      * P_k = betainc(dof_k/2, 0.5, dof_k/(dof_k + t_k^2))   (== 2 * t.sf(|t|, dof_k))
      * Z_k = -ndtri(P_k / 2)                                (tail-safe; NOT ndtri(1 - P/2))
      * A perfect fit |r_k| == 1 is a VALID extreme result, not an error: 1 - r^2 is
        clamped at 0, so t_k = +/-inf, P_k = 0, Z_k = +inf, W = +inf -> P_W = 0 ->
        PASS, with beta kept and SE = 0. Numerical underflow of P_k to 0 at a finite
        t is kept the same way. The strongest possible signal is never deleted.
      * Only NaN in r/t/P/Z (monomorphic / all-missing) or a clearly illegal
        |r_k| > 1 + R_TOL is invalid (excluded from K_est, counted in n_invalid_stats).
      * |r_k| < 1 but very close to 1 is left EXACTLY as computed -- it is never
        snapped to a perfect fit. Such near-perfect effects are counted (n_near_perfect)
        so those few pairs can be re-checked in float64 later if needed. Perfect fits
        (|r_k| >= 1) are counted separately (n_perfect_fit).

    Returns W, K_est, P_min, strongest, P_list, valid_list,
            n_invalid_stats, n_perfect_fit, n_near_perfect
    where the last two are per-PAIR counts (pairs with >= 1 such subgroup effect).
    """
    K = len(r_src)
    n = int(np.asarray(r_src[0]).size)
    W = np.zeros(n, dtype=np.float64)
    K_est = np.zeros(n, dtype=np.int32)
    P_min = np.full(n, np.inf, dtype=np.float64)
    strongest = np.full(n, -1, dtype=np.int32)
    P_list = [None] * K
    valid_list = [None] * K
    n_invalid_stats = 0
    perfect_any = None      # per-pair: any subgroup is a perfect fit (|r| >= 1)
    near_any = None         # per-pair: any subgroup is |r| < 1 but near-perfect

    for k in range(K):
        r = np.asarray(r_src[k], dtype=np.float64).reshape(-1)
        est = np.asarray(est_src[k]).reshape(-1)
        if pc_broadcast is not None:
            est = np.repeat(est, pc_broadcast)
        dof = dof_k[k]
        r2 = r * r
        # clamp 1 - r^2 at 0: |r| >= 1 is a perfect fit -> t = +/-inf -> P = 0 -> Z = +inf
        denom = np.maximum(1.0 - r2, 0.0)
        with np.errstate(divide="ignore", invalid="ignore"):
            t = r * np.sqrt(dof / denom)
            t2 = t * t
            P = sp.betainc(dof / 2.0, 0.5, dof / (dof + t2))
            Z = -sp.ndtri(P / 2.0)
        # |r| <= 1: normal.  |r| in (1, 1+R_TOL]: float32 overshoot of a perfect fit.
        # |r| > 1+R_TOL: clearly illegal -> invalid.  +inf Z (P underflow or a perfect
        # fit) is KEPT: W = +inf -> P_W = 0 -> PASS.
        abs_r = np.abs(r)
        illegal = abs_r > (1.0 + R_TOL)
        bad = np.isnan(t) | np.isnan(P) | np.isnan(Z) | illegal
        valid = est & ~bad
        n_invalid_stats += int(np.count_nonzero(est & bad))
        # telemetry only (never changes the statistics): perfect fits, and
        # near-perfect effects that are strictly |r| < 1 (float64-recompute candidates)
        perfect_k = valid & (abs_r >= 1.0)
        near_k = valid & (abs_r > R_NEAR_PERFECT) & (abs_r < 1.0)
        if perfect_any is None:
            perfect_any, near_any = perfect_k, near_k
        else:
            perfect_any |= perfect_k
            near_any |= near_k
        with np.errstate(invalid="ignore"):
            W += np.where(valid, Z * Z, 0.0)
        K_est += valid.astype(np.int32)
        P_masked = np.where(valid, P, np.inf)
        new_min = P_masked < P_min
        P_min = np.where(new_min, P_masked, P_min)
        strongest = np.where(new_min, k, strongest)
        P_list[k] = P
        valid_list[k] = valid
        del r, est, t, t2, P, Z, denom, abs_r, illegal, bad, valid, perfect_k, near_k, P_masked, new_min
    n_perfect_fit = int(np.count_nonzero(perfect_any)) if perfect_any is not None else 0
    n_near_perfect = int(np.count_nonzero(near_any)) if near_any is not None else 0
    return (W, K_est, P_min, strongest, P_list, valid_list,
            n_invalid_stats, n_perfect_fit, n_near_perfect)


def _prefilter_zgate(omnibus_p_threshold, pmin_threshold, Kmax):
    """SAFE combined z-gate for PASS = (P_W < omnibus_thr) OR (P_min < pmin_thr).

    omnibus branch: if W > Wcrit_m for m valid subgroups then, since
        W = sum(Z_k^2) <= m * max(Z_k^2),
    some valid subgroup has |Z_k| >= sqrt(Wcrit_m / m) >= min_m sqrt(Wcrit_m/m).
    pmin branch: P_min < pmin_thr  =>  some valid |Z_k| > norm.isf(pmin_thr/2).

    The combined gate is the more permissive (smaller) of the two, so no final
    PASS pair can be excluded. Wcrit_m = chi2.isf(omnibus_thr, df=m) = chi2.ppf(1-thr, m).
    """
    zgate_W = min(np.sqrt(stats.chi2.isf(omnibus_p_threshold, df=m) / float(m))
                  for m in range(1, Kmax + 1))
    zgate_min = stats.norm.isf(pmin_threshold / 2.0)
    return float(min(zgate_W, zgate_min))


def _subgroup_r_gate(zgate, dof_k):
    """TensorQTL correlation gate: |r| at which |Z(r)| == zgate, for this dof."""
    p_gate = 2.0 * stats.norm.sf(zgate)
    t_gate = stats.t.isf(p_gate / 2.0, dof_k)
    return float(t_gate / np.sqrt(dof_k + t_gate ** 2))


def _merge_pass_stats(stats, P_W, P_min, pass_pos, omnibus_p_threshold, pmin_threshold):
    stats["pass_pairs"] += int(pass_pos.size)
    if pass_pos.size == 0:
        return
    pw = P_W[pass_pos] < omnibus_p_threshold
    pm = P_min[pass_pos] < pmin_threshold
    stats["pw_only"] += int(np.count_nonzero(pw & ~pm))
    stats["pmin_only"] += int(np.count_nonzero(~pw & pm))
    stats["both"] += int(np.count_nonzero(pw & pm))


def _pack_table(pass_pos, ci, cj, variant_ids, phenotype_ids, P_W, W, K_est,
                P_min, strongest, P_list, r_c, gvar_c, pvar_c, maf_c,
                valid_list, sub_names, dof_k, schema, observed_maf_c=None):
    """Assemble only passing pairs as Arrow columns; no per-pair Python objects."""
    if observed_maf_c is None:  # preserve existing direct helper callers
        observed_maf_c = maf_c
    a = np.asarray(pass_pos, dtype=np.int64)
    ii = np.asarray(ci)[a]
    jj = np.asarray(cj)[a]
    name_array = np.asarray(sub_names, dtype=object)
    strongest_ix = np.asarray(strongest)[a]
    strongest_name = np.empty(len(a), dtype=object)
    strongest_name[:] = None
    good_strongest = strongest_ix >= 0
    strongest_name[good_strongest] = name_array[strongest_ix[good_strongest]]
    cols = {
        "variant_id": np.asarray(variant_ids, dtype=object)[ii],
        "phenotype_id": np.asarray(phenotype_ids, dtype=object)[jj],
        "P_W": np.asarray(P_W[a], dtype=np.float64),
        "W": np.asarray(W[a], dtype=np.float64),
        "K_est": np.asarray(K_est[a], dtype=np.int32),
        "P_min": np.asarray(P_min[a], dtype=np.float64),
        "strongest_subgroup": strongest_name,
    }
    for k, lab in enumerate(sub_names):
        valid = np.asarray(valid_list[k])[a] & np.isfinite(np.asarray(r_c[k])[a])
        r = np.asarray(r_c[k], dtype=np.float64)[a]
        gv = np.asarray(gvar_c[k], dtype=np.float64)[a]
        pv = np.asarray(pvar_c[k], dtype=np.float64)[a]
        observed_mf = np.asarray(observed_maf_c[k], dtype=np.float64)[a]
        with np.errstate(divide="ignore", invalid="ignore"):
            beta = r * np.sqrt(pv / gv)
            se = np.sqrt(np.maximum(1.0 - r * r, 0.0) * pv / (dof_k[k] * gv))
        # Association statistics retain their existing values and estimability mask.
        cols[f"{lab}_beta"] = pa.array(np.where(valid, beta, np.nan).astype(np.float32),
                                       mask=~valid, type=pa.float32())
        cols[f"{lab}_se"] = pa.array(np.where(valid, se, np.nan).astype(np.float32),
                                     mask=~valid, type=pa.float32())
        cols[f"{lab}_p"] = pa.array(np.where(valid, np.asarray(P_list[k])[a], np.nan),
                                    mask=~valid, type=pa.float64())
        # Output MAF always comes from called genotypes before imputation.
        # The separate post-imputation MAF remains only in the original estimability gate.
        display_maf = observed_mf.astype(np.float32)
        cols[f"{lab}_maf"] = pa.array(display_maf, mask=~np.isfinite(display_maf),
                                      type=pa.float32())
    return pa.Table.from_pydict(cols, schema=schema)


def _t_upper_gate(r, valid, dof, tcrit):
    """GPU per-subgroup upper bound. float64 and outward tolerance keep boundaries safe."""
    rd = r.to(torch.float64)
    abs_r = rd.abs()
    usable = valid[:, None] & torch.isfinite(rd) & (abs_r <= 1.0 + R_TOL)
    r2 = rd.square()
    with torch.no_grad():
        t2 = torch.where(r2 >= 1.0, torch.full_like(r2, float("inf")),
                         dof * r2 / (1.0 - r2))
        t2 = torch.where(usable, t2, torch.zeros_like(t2))
        pmin = usable & (torch.sqrt(t2) >= tcrit * (1.0 - 1e-10))
    return t2, usable.to(torch.uint8), pmin

def _tile_exact(sub_g, sub_maf, sub_observed_maf, pheno_chunk_t, sub_idx_t, dof_k, maf_threshold,
                omnibus_p_threshold, pmin_threshold,
                variant_ids_batch, pheno_ids_chunk, sub_names, stats, schema, output_all=False):
    """Exact golden path (--no_prefilter): exact transform on EVERY pair."""
    K = len(sub_g)
    ng = sub_g[0].shape[0]
    Pc = pheno_chunk_t.shape[0]
    stats["total_pairs"] += ng * Pc
    stats["exact_pairs"] += ng * Pc

    t0 = time.time()
    r_full = [None] * K
    gvar_full = [None] * K
    pvar_full = [None] * K
    maf_full = [None] * K
    observed_maf_full = [None] * K
    est_full = [None] * K
    for k in range(K):
        r, gvar, pvar = calculate_corr(sub_g[k], pheno_chunk_t[:, sub_idx_t[k]],
                                       residualizer=None, return_var=True)
        r_full[k] = r.cpu().numpy()                              # float32 (ng, Pc)
        gvar_full[k] = gvar.cpu().numpy().astype(np.float64)     # (ng,)
        pvar_full[k] = pvar.cpu().numpy().astype(np.float64)     # (Pc,)
        maf_full[k] = sub_maf[k].cpu().numpy().astype(np.float64)
        observed_maf_full[k] = sub_observed_maf[k].cpu().numpy()
        est_full[k] = ((maf_full[k] >= maf_threshold) & (gvar_full[k] > 0.0)
                       & np.isfinite(gvar_full[k]) & np.isfinite(maf_full[k]))
        del r, gvar, pvar
    stats["t_gpu"] += time.time() - t0

    t0 = time.time()
    W, K_est, P_min, strongest, P_list, valid_list, n_invalid, n_perfect, n_near = _exact_stats_block(
        r_full, est_full, dof_k, pc_broadcast=Pc)
    P_W = _grouped_chi2_sf(W, K_est)
    stats["invalid_stats"] += n_invalid
    stats["perfect_fit_pairs"] += n_perfect
    stats["near_perfect_pairs"] += n_near
    stats["t_cpu"] += time.time() - t0

    PASS = (K_est >= 1) if output_all else ((P_W < omnibus_p_threshold) | (P_min < pmin_threshold))
    pass_pos = np.flatnonzero(PASS)
    if output_all:
        stats["pass_pairs"] += int(pass_pos.size)
    else:
        _merge_pass_stats(stats, P_W, P_min, pass_pos, omnibus_p_threshold, pmin_threshold)
    if pass_pos.size == 0:
        return pa.Table.from_pylist([], schema=schema)
    t0 = time.time()
    table = _pack_table(pass_pos, np.repeat(np.arange(ng), Pc),
                        np.tile(np.arange(Pc), ng), variant_ids_batch,
                        pheno_ids_chunk, P_W, W, K_est, P_min, strongest,
                        P_list, [x.reshape(-1) for x in r_full],
                        [np.repeat(x, Pc) for x in gvar_full],
                        [np.tile(x, ng) for x in pvar_full],
                        [np.repeat(x, Pc) for x in maf_full],
                        valid_list, sub_names, dof_k, schema,
                        observed_maf_c=[np.repeat(x, Pc) for x in observed_maf_full])
    stats["t_assembly"] += time.time() - t0
    return table


def _tile_prefilter(sub_g, sub_maf, sub_observed_maf, pheno_chunk_t, sub_idx_t, dof_k, maf_threshold,
                    omnibus_p_threshold, pmin_threshold, r_gate_k, wcrit_k, tcrit_k,
                    variant_ids_batch, pheno_ids_chunk, sub_names, stats, schema, prefilter_mode):
    """SAFE GPU prefilter path (two-pass tile).

    PASS A: per-subgroup calculate_corr -> candidate mask built on GPU only, no
            CPU special function is evaluated.
    PASS B: recompute per-subgroup calculate_corr (one dense matrix at a time) and
            gather ONLY the candidate elements, then run the exact transform on the
            candidate subset. No dense per-subgroup t/beta/SE matrix is retained.
    """
    K = len(sub_g)
    ng = sub_g[0].shape[0]
    Pc = pheno_chunk_t.shape[0]
    stats["total_pairs"] += ng * Pc

    # ---------------- PASS A: GPU candidate mask (cheap, safe superset) ----------------
    t0 = time.time()
    candidate = None
    T_upper = K_upper = possible_pmin = None
    for k in range(K):
        r, gvar, pvar = calculate_corr(sub_g[k], pheno_chunk_t[:, sub_idx_t[k]],
                                       residualizer=None, return_var=True)
        maf = sub_maf[k]
        est = (maf >= maf_threshold) & (gvar > 0.0) & torch.isfinite(gvar) & torch.isfinite(maf)
        if prefilter_mode == "legacy":
            gate = torch.as_tensor(float(r_gate_k[k]), dtype=r.dtype, device=r.device)
            hit = est[:, None] & (r.abs() >= gate)
            candidate = hit if candidate is None else (candidate | hit)
            del hit, gate
        else:
            t2, valid, hit = _t_upper_gate(r, est, float(dof_k[k]), float(tcrit_k[k]))
            if T_upper is None:
                T_upper, K_upper, possible_pmin = t2, valid, hit
            else:
                T_upper += t2
                K_upper += valid
                possible_pmin |= hit
            del t2, valid, hit
        del r, gvar, pvar, maf, est
    if prefilter_mode == "t_upper":
        crit = torch.as_tensor(wcrit_k, dtype=torch.float64, device=T_upper.device)
        candidate = ((K_upper > 0) & (T_upper >= crit[K_upper.long()] * (1.0 - 1e-10))) | possible_pmin
        del T_upper, K_upper, possible_pmin, crit
    stats["t_gpu"] += time.time() - t0

    ncand = int(candidate.sum().item())
    stats["cand_pairs"] += ncand
    if ncand == 0:
        del candidate
        return pa.Table.from_pylist([], schema=schema)
    ix = candidate.nonzero(as_tuple=False)          # (ncand, 2) on device
    del candidate
    ci = ix[:, 0]
    cj = ix[:, 1]

    # ---------------- PASS B: gather candidate elements, then exact on CPU ----------------
    t0 = time.time()
    r_c = [None] * K
    gvar_c = [None] * K
    pvar_c = [None] * K
    maf_c = [None] * K
    observed_maf_c = [None] * K
    est_c = [None] * K
    for k in range(K):
        r, gvar, pvar = calculate_corr(sub_g[k], pheno_chunk_t[:, sub_idx_t[k]],
                                       residualizer=None, return_var=True)
        rc = r[ci, cj].cpu().numpy()                            # float32 (ncand,)
        gv = gvar[ci].cpu().numpy().astype(np.float64)
        pv = pvar[cj].cpu().numpy().astype(np.float64)
        mf = sub_maf[k][ci].cpu().numpy().astype(np.float64)
        observed_maf_c[k] = sub_observed_maf[k][ci].cpu().numpy()
        r_c[k] = rc; gvar_c[k] = gv; pvar_c[k] = pv; maf_c[k] = mf
        est_c[k] = (mf >= maf_threshold) & (gv > 0.0) & np.isfinite(gv) & np.isfinite(mf)
        del r, gvar, pvar
    stats["t_gpu"] += time.time() - t0

    t0 = time.time()
    stats["exact_pairs"] += ncand
    W, K_est, P_min, strongest, P_list, valid_list, n_invalid, n_perfect, n_near = _exact_stats_block(
        r_c, est_c, dof_k, pc_broadcast=None)
    P_W = _grouped_chi2_sf(W, K_est)
    stats["invalid_stats"] += n_invalid
    stats["perfect_fit_pairs"] += n_perfect
    stats["near_perfect_pairs"] += n_near
    stats["t_cpu"] += time.time() - t0

    PASS = (P_W < omnibus_p_threshold) | (P_min < pmin_threshold)
    pass_pos = np.flatnonzero(PASS)
    _merge_pass_stats(stats, P_W, P_min, pass_pos, omnibus_p_threshold, pmin_threshold)
    if pass_pos.size == 0:
        return pa.Table.from_pylist([], schema=schema)
    t0 = time.time()
    ci_np = ci.cpu().numpy(); cj_np = cj.cpu().numpy()
    table = _pack_table(pass_pos, ci_np, cj_np, variant_ids_batch, pheno_ids_chunk,
                        P_W, W, K_est, P_min, strongest, P_list, r_c, gvar_c,
                        pvar_c, maf_c, valid_list, sub_names, dof_k, schema,
                        observed_maf_c=observed_maf_c)
    stats["t_assembly"] += time.time() - t0
    return table


def _subgroup_parquet_schema(labels):
    fields = [pa.field("variant_id", pa.string()),
              pa.field("phenotype_id", pa.string()),
              pa.field("P_W", pa.float64()), pa.field("W", pa.float64()),
              pa.field("K_est", pa.int32()), pa.field("P_min", pa.float64()),
              pa.field("strongest_subgroup", pa.string())]
    for lab in labels:
        for suffix, dtype in (("beta", pa.float32()), ("se", pa.float32()),
                              ("p", pa.float64()), ("maf", pa.float32())):
            fields.append(pa.field(f"{lab}_{suffix}", dtype))
    return pa.schema(fields)


def _id_digest(values):
    h = hashlib.sha256()
    for value in values:
        h.update(str(value).encode("utf-8"))
        h.update(b"\0")
    return h.hexdigest()


def _source_fingerprints(input_paths):
    records = {}
    for key, path in (input_paths or {}).items():
        p = Path(path).resolve()
        st = p.stat()
        records[key] = dict(path=str(p), size=st.st_size, mtime_ns=st.st_mtime_ns)
    return records


def run_subgroup_omnibus(genotype_df, phenotype_df, subgroup_s,
                         pval_threshold=1e-5, maf_threshold=0.05, batch_size=20000,
                         pheno_chunk=None, omnibus_p_threshold=None, pmin_threshold=None,
                         prefilter=True, stats_out=None, logger=None, verbose=True,
                         parts_dir=None, max_rows_per_part=100000,
                         prefilter_mode="t_upper", resume=False, input_paths=None):
    """Exact subgroup association with safe GPU gates and tile-level Parquet checkpoints."""
    validate_parameters(pval_threshold, maf_threshold, batch_size, pheno_chunk,
                        omnibus_p_threshold, pmin_threshold)
    if prefilter_mode not in ("legacy", "t_upper"):
        raise ValueError(f"prefilter_mode must be legacy or t_upper; got {prefilter_mode!r}.")
    if resume and parts_dir is None:
        raise ValueError("--resume requires a subgroup parts directory.")
    validate_matrices(genotype_df, phenotype_df)
    validate_alignment(genotype_df, phenotype_df)
    counts = validate_subgroups(subgroup_s, phenotype_df.columns)
    device = _get_device()
    if logger is None:
        logger = SimpleLogger(verbose=verbose)
    if omnibus_p_threshold is None:
        omnibus_p_threshold = pval_threshold
    if pmin_threshold is None:
        pmin_threshold = pval_threshold
    output_all = pval_threshold == 1.0
    if output_all:
        prefilter = False
    labels = list(counts)
    K = len(labels)
    assigned = subgroup_s.reindex(phenotype_df.columns).astype(str).to_numpy()
    sub_idx = [np.flatnonzero(assigned == lab).astype(np.int64) for lab in labels]
    n_k = np.array([len(ix) for ix in sub_idx], dtype=np.int64)
    if np.any(n_k <= 2):
        bad = ", ".join(f"{lab}={n}" for lab, n in zip(labels, n_k) if n <= 2)
        raise ValueError(f"Each subgroup needs n > 2 samples (dof = n - 2 > 0). Problematic: {bad}")
    dof_k = (n_k - 2).astype(np.float64)
    logger.write("all-pair association (subgroup_omnibus)")
    logger.write(f"  * subgroups: {K} ({', '.join(labels)}), Kmax={K}")
    for lab, n, d in zip(labels, n_k, dof_k):
        logger.write(f"      {lab}: n={n}, dof={d:.0f}")
    logger.write(f"  * maf_threshold = {maf_threshold:.3g}")
    logger.write(f"  * omnibus_p_threshold = {omnibus_p_threshold:.3g}, "
                 f"pmin_threshold = {pmin_threshold:.3g}")
    if output_all:
        logger.write("  * OUTPUT-ALL: exact statistics for every K_est>=1 pair")
    if prefilter:
        logger.write(f"  * PRE-FILTER: {prefilter_mode}")
        if prefilter_mode == "legacy":
            zgate = _prefilter_zgate(omnibus_p_threshold, pmin_threshold, K)
            r_gate_k = np.array([_subgroup_r_gate(zgate, d) for d in dof_k])
            logger.write(f"      zgate={zgate:.6g}; r gates: " +
                         ", ".join(f"{lab}={g:.6g}" for lab, g in zip(labels, r_gate_k)))
            wcrit_k = tcrit_k = None
        else:
            r_gate_k = None
            wcrit_k = np.r_[np.inf, [stats.chi2.isf(omnibus_p_threshold, df=k)
                                      for k in range(1, K + 1)]]
            tcrit_k = stats.t.isf(pmin_threshold / 2.0, dof_k)
            logger.write("      safe bound: W_exact <= sum(t_k^2); direct t gate for P_min")
    else:
        r_gate_k = wcrit_k = tcrit_k = None
        logger.write("  * PRE-FILTER: DISABLED (exact reference)")

    phenotypes_t = torch.tensor(phenotype_df.values, dtype=torch.float32).to(device)
    genotype_ix = np.array([genotype_df.columns.tolist().index(i) for i in phenotype_df.columns])
    genotype_ix_t = torch.from_numpy(genotype_ix).to(device)
    sub_idx_t = [torch.from_numpy(ix).to(device) for ix in sub_idx]
    P = phenotype_df.shape[0]
    chunks = [(0, P)] if pheno_chunk is None else [
        (i, min(i + pheno_chunk, P)) for i in range(0, P, pheno_chunk)]
    schema = _subgroup_parquet_schema(labels)
    stat_keys = ("total_pairs", "cand_pairs", "exact_pairs", "pass_pairs", "pw_only",
                 "pmin_only", "both", "invalid_stats", "perfect_fit_pairs",
                 "near_perfect_pairs", "t_gpu", "t_cpu", "t_assembly", "t_write")
    numeric = {key: 0.0 if key.startswith("t_") else 0 for key in stat_keys}
    all_tables = [] if parts_dir is None else None
    parts_path = Path(parts_dir) if parts_dir is not None else None
    if parts_path is not None:
        signature = dict(format_version=2, mode="subgroup_omnibus",
                         code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                         variant_count=len(genotype_df), phenotype_count=P,
                         genotype_ids_sha256=_id_digest(genotype_df.index),
                         phenotype_ids_sha256=_id_digest(phenotype_df.index),
                         samples_sha256=_id_digest(phenotype_df.columns),
                         subgroup_sha256=_id_digest(assigned), labels=labels,
                         pval_threshold=pval_threshold,
                         omnibus_p_threshold=omnibus_p_threshold,
                         pmin_threshold=pmin_threshold, maf_threshold=maf_threshold,
                         batch_size=batch_size, pheno_chunk=pheno_chunk,
                         prefilter_mode=prefilter_mode if prefilter else "off",
                         input_sources=_source_fingerprints(input_paths))
        run_marker = parts_path / "RUN.json"
        if resume:
            if not run_marker.is_file():
                raise ValueError(f"Cannot resume: run manifest missing: {run_marker}")
            stored = json.loads(run_marker.read_text())
            if stored != signature:
                raise ValueError("Cannot resume: inputs, tile layout, or parameters differ from RUN.json.")
            if (parts_path / "COMPLETE.json").exists():
                raise ValueError("This subgroup run is already complete; run 03_merge_subgroup_parts.py.")
            logger.write(f"  * RESUME: validating completed tiles in {parts_path}")
        else:
            if parts_path.exists():
                raise ValueError(f"Temporary parts directory already exists: {parts_path}. "
                                 "Use --resume for this run or choose a new --out path.")
            parts_path.mkdir()
            run_marker.write_text(json.dumps(signature, indent=2) + "\n")
    ggt = GenotypeGeneratorTrans(genotype_df, batch_size=batch_size)
    start_time = time.time()
    resumed_tiles = 0
    tile_entries = []
    for batch_id, (genotypes, variant_ids_batch) in enumerate(ggt.generate_data(verbose=verbose), 1):
        # Skip genotype transfer and subgroup imputation when every phenotype tile is complete.
        if parts_path is not None and resume:
            markers = [parts_path / f"tile-g{batch_id:06d}-p{j:06d}.json"
                       for j in range(1, len(chunks) + 1)]
            if all(path.is_file() for path in markers):
                for path in markers:
                    entry = json.loads(path.read_text())
                    _verify_subgroup_tile(parts_path, entry, schema)
                    for key in stat_keys:
                        numeric[key] += entry["stats"][key]
                    tile_entries.append(entry)
                    resumed_tiles += 1
                continue
        vid = list(variant_ids_batch)
        genotypes_t = torch.tensor(genotypes, dtype=torch.float).to(device)
        genotypes_t = genotypes_t[:, genotype_ix_t]
        sub_g = [None] * K
        sub_maf = [None] * K
        sub_observed_maf = [None] * K
        for kk in range(K):
            g = genotypes_t[:, sub_idx_t[kk]]
            # Output-only variant MAF from real calls, before TensorQTL imputation.
            called = (g != -9) & torch.isfinite(g)
            n_called = called.sum(1)
            allele_sum = torch.where(called, g, torch.zeros_like(g)).sum(1)
            observed_af = allele_sum / (2.0 * n_called)
            observed_maf = torch.minimum(observed_af, 1.0 - observed_af)
            sub_observed_maf[kk] = torch.where(
                n_called > 0, observed_maf, torch.full_like(observed_maf, float("nan")))
            # Preserve the original imputation and MAF used for estimability.
            impute_mean(g)
            af = g.sum(1) / (2.0 * float(n_k[kk]))
            sub_g[kk] = g
            sub_maf[kk] = torch.where(af > 0.5, 1 - af, af)
            del called, n_called, allele_sum, observed_af, observed_maf
        for chunk_id, (c0, c1) in enumerate(chunks, 1):
            tile_id = f"tile-g{batch_id:06d}-p{chunk_id:06d}"
            if parts_path is not None:
                tile_path = parts_path / f"{tile_id}.parquet"
                tile_marker = parts_path / f"{tile_id}.json"
                if resume and tile_marker.is_file():
                    entry = json.loads(tile_marker.read_text())
                    _verify_subgroup_tile(parts_path, entry, schema)
                    for key in stat_keys:
                        numeric[key] += entry["stats"][key]
                    tile_entries.append(entry)
                    resumed_tiles += 1
                    continue
                # A parquet without its atomic completion marker is not a complete tile.
                if tile_path.exists():
                    tile_path.unlink()
                partial = parts_path / f"{tile_id}.parquet.partial"
                if partial.exists():
                    partial.unlink()
            before = numeric.copy()
            tile_args = (sub_g, sub_maf, sub_observed_maf, phenotypes_t[c0:c1], sub_idx_t, dof_k,
                         maf_threshold, omnibus_p_threshold, pmin_threshold)
            pid = list(phenotype_df.index[c0:c1])
            if prefilter:
                table = _tile_prefilter(*tile_args, r_gate_k, wcrit_k, tcrit_k,
                                        vid, pid, labels, numeric, schema, prefilter_mode)
            else:
                table = _tile_exact(*tile_args, vid, pid, labels, numeric, schema,
                                    output_all=output_all)
            if parts_path is None:
                all_tables.append(table)
            else:
                t0 = time.time()
                pq.write_table(table, partial, compression="snappy")
                os.replace(partial, tile_path)
                numeric["t_write"] += time.time() - t0
                delta = {key: numeric[key] - before[key] for key in stat_keys}
                entry = dict(tile_id=tile_id, batch_id=batch_id, chunk_id=chunk_id,
                             name=tile_path.name, marker=tile_marker.name,
                             rows=table.num_rows, stats=delta)
                mark_partial = parts_path / f"{tile_id}.json.partial"
                mark_partial.write_text(json.dumps(entry) + "\n")
                os.replace(mark_partial, tile_marker)
                tile_entries.append(entry)
            del table
        del genotypes_t, sub_g, sub_maf, sub_observed_maf
    total_wall = time.time() - start_time
    del phenotypes_t
    total = numeric["total_pairs"]
    logger.write("  * QC / performance:")
    logger.write(f"      total_pairs_examined      = {total:,}")
    if prefilter:
        frac = (numeric["cand_pairs"] / total) if total else float("nan")
        logger.write(f"      prefilter_candidate_pairs = {numeric['cand_pairs']:,}")
        logger.write(f"      prefilter_fraction        = {frac:.6g}")
    logger.write(f"      exact_pairs_processed     = {numeric['exact_pairs']:,}")
    logger.write(f"      final_pass_pairs          = {numeric['pass_pairs']:,}")
    for label, key in (("PW_only", "pw_only"), ("Pmin_only", "pmin_only"),
                       ("both", "both"), ("invalid_statistics", "invalid_stats"),
                       ("perfect_fit_pairs", "perfect_fit_pairs"),
                       ("near_perfect_pairs", "near_perfect_pairs")):
        logger.write(f"      {label:25s}= {numeric[key]:,}")
    for label, key in (("GPU association/prefilter", "t_gpu"),
                       ("CPU exact transform", "t_cpu"),
                       ("row assembly", "t_assembly"), ("parquet write", "t_write")):
        logger.write(f"      wall time {label:25s}= {numeric[key]:.3f} s")
    logger.write(f"      wall time this invocation      = {total_wall:.3f} s")
    if stats_out is not None:
        stats_out.update(numeric)
        stats_out.update(total_wall=total_wall, prefilter=bool(prefilter),
                         prefilter_mode=prefilter_mode if prefilter else "off",
                         output_all=output_all, Kmax=K, labels=labels,
                         resumed_tiles=resumed_tiles)
    if parts_path is not None:
        expected_tiles = ((len(genotype_df) + batch_size - 1) // batch_size) * len(chunks)
        if len(tile_entries) != expected_tiles:
            raise RuntimeError(f"Completed tile count {len(tile_entries)} != expected {expected_tiles}.")
        if sum(x["rows"] for x in tile_entries) != numeric["pass_pairs"]:
            raise RuntimeError("Parquet tile row counts differ from final_pass_pairs.")
        manifest = dict(format_version=2, mode="subgroup_omnibus",
                        parts_dir=str(parts_path.resolve()), tiles=tile_entries,
                        expected_tiles=expected_tiles, n_rows=numeric["pass_pairs"],
                        labels=labels, stats=numeric)
        marker_tmp = parts_path / "COMPLETE.json.partial"
        marker_tmp.write_text(json.dumps(manifest) + "\n")
        os.replace(marker_tmp, parts_path / "COMPLETE.json")
        logger.write(f"  * completed {len(tile_entries):,} tiles ({resumed_tiles:,} resumed)")
        logger.write("done; merge parts to create the final parquet.")
        return None
    if not all_tables:
        return pa.Table.from_pylist([], schema=schema).to_pandas()
    return pa.concat_tables(all_tables).to_pandas()


def _verify_subgroup_tile(parts_path, entry, schema):
    tile_id = entry.get("tile_id")
    expect = f"tile-g{entry.get('batch_id'):06d}-p{entry.get('chunk_id'):06d}"
    if tile_id != expect or entry.get("name") != expect + ".parquet" or entry.get("marker") != expect + ".json":
        raise ValueError(f"Invalid tile completion marker: {entry}")
    file = parts_path / entry["name"]
    if not file.is_file():
        raise ValueError(f"Completed tile parquet is missing: {file}")
    reader = pq.ParquetFile(file)
    try:
        if reader.metadata.num_rows != entry["rows"] or not reader.schema_arrow.equals(schema):
            raise ValueError(f"Completed tile row count/schema mismatch: {file}")
    finally:
        reader.close()


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def write_output(df, path):
    if path.endswith(".parquet"):
        df.to_parquet(path, index=False)
    else:
        compression = "gzip" if path.endswith(".gz") else None
        df.to_csv(path, sep="\t", index=False, compression=compression, float_format="%.6g")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main():
    p = argparse.ArgumentParser(
        description="Standalone all-pair caQTL association on installed TensorQTL "
                    "(modes: global_pc, subgroup_omnibus, cis_all, selected).")
    p.add_argument("--mode", required=True, choices=["global_pc", "subgroup_omnibus", "cis_all", "selected"])
    p.add_argument("--genotype",
                   help="Genotypes: PLINK prefix (bed/bim/fam or pgen/psam/pvar), parquet, or tsv.gz")
    p.add_argument("--phenotype",
                   help="INT phenotype matrix (phenotypes x samples): BED, parquet, or TSV")
    p.add_argument("--covariates", default=None,
                   help="Covariates (covariates x samples). Optional for global_pc: "
                        "if omitted the model is y ~ G (no PCs, dof = n - 2).")
    p.add_argument("--subgroups", default=None,
                   help="Subgroup file: headerless TSV (sample_id \\t subgroup). "
                        "Observed analysis labels define the groups; each needs >=3 samples.")
    p.add_argument("--maf_threshold", type=float, default=0.05)
    p.add_argument("--pval_threshold", type=float, default=1e-5,
                   help="0 < p <= 1. At p=1, subgroup_omnibus outputs every K_est>=1 pair "
                        "using the exact path, ignoring prefilter and PW/Pmin filters. "
                        "Otherwise defaults both subgroup thresholds when omitted.")
    p.add_argument("--omnibus_p_threshold", type=float, default=None,
                   help="PASS threshold for P_W (subgroup_omnibus). Default: --pval_threshold")
    p.add_argument("--pmin_threshold", type=float, default=None,
                   help="PASS threshold for P_min (subgroup_omnibus). Default: --pval_threshold")
    p.add_argument("--no_prefilter", action="store_true",
                   help="subgroup_omnibus: disable the safe GPU prefilter and run the exact "
                        "transform on every pair (numerical golden path). Default: prefilter ON.")
    p.add_argument("--prefilter_mode", choices=["legacy", "t_upper"], default="t_upper",
                   help="subgroup_omnibus safe GPU gate (default: t_upper).")
    p.add_argument("--resume", action="store_true",
                   help="Resume a matching subgroup run from complete tile markers.")
    p.add_argument("--batch_size", type=int, default=20000)
    p.add_argument("--pheno_chunk", type=int, default=None,
                   help="Phenotypes per global/subgroup tile (default: all at once); cis_all processes one phenotype at a time")
    p.add_argument("--out", required=True,
                   help="Final output path. subgroup_omnibus requires .parquet and writes "
                        "temporary parts first; global_pc also supports TSV/CSV. "
                        "cis_all requires an existing empty output directory; selected uses a new batch directory.")
    p.add_argument("--logfile", default=None, help="Optional log file.")
    p.add_argument("--window", type=int, default=1000000,
                   help="cis_all only: inclusive distance around TSS/ACR midpoint (default: 1000000 bp)")
    p.add_argument("--phenotype_kind", "--phenotype-kind", choices=["eqtl", "caqtl"], default=None,
                   help="cis_all: choose TSS/ACR midpoint; selected: optional ID kind, otherwise inferred from selectors")
    p.add_argument("--phenotype_positions", "--phenotype-positions", default=None,
                   help="cis_all only: BED6 gene annotation or BED4 ACR positions; BED coordinates are 0-based half-open")
    p.add_argument("--variant_metadata", "--variant-metadata", default=None,
                   help="cis_all/selected: explicit 1-based variant position/counted-allele TSV for matrix genotypes")
    p.add_argument("--method", choices=["global", "pw", "sg1", "sg2", "sg3", "sg5"],
                   help="selected only: one association route")
    p.add_argument("--acr-id", action="append", default=[], help="selected: repeatable short or full ACR ID")
    p.add_argument("--acr-list", help="selected: one ACR ID per line")
    p.add_argument("--gene-id", action="append", default=[], help="selected: repeatable exact gene ID")
    p.add_argument("--gene-list", help="selected: one exact gene ID per line")
    plots = p.add_mutually_exclusive_group()
    plots.add_argument("--plot", action="store_true", help="selected: also draw genome Manhattan PDF/PNG")
    plots.add_argument("--plot-only", action="store_true", help="selected: draw existing results without loading genotype/phenotype")
    p.add_argument("--plot-region", action="append", default=[], metavar="CHR:START-END",
                   help="selected plotting: optional local plots, 1-based inclusive bp; repeatable")
    args = p.parse_args()

    if not (args.mode == "selected" and args.plot_only):
        missing = ["--"+name for name in ("genotype", "phenotype") if getattr(args, name) is None]
        if missing:
            p.error("the following arguments are required: " + ", ".join(missing))
    try:
        if args.mode == "selected":
            import importlib.util
            spec = importlib.util.spec_from_file_location("tensorqtl_selected", Path(__file__).with_name("02_selected_mode_run.py"))
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module.run(args, sys.modules[__name__])
            return
        if args.method or args.acr_id or args.acr_list or args.gene_id or args.gene_list or args.plot or args.plot_only or args.plot_region:
            raise ValueError("ID selectors, --method and plotting options require --mode selected.")
        if args.mode == "cis_all":
            import importlib.util
            spec = importlib.util.spec_from_file_location("tensorqtl_cis_all", Path(__file__).with_name("01_cis_all_mode_run.py"))
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)
            module.run(args, sys.modules[__name__])
            return
        validate_parameters(args.pval_threshold, args.maf_threshold, args.batch_size,
                            args.pheno_chunk, args.omnibus_p_threshold, args.pmin_threshold)
        validate_output_parent(args.out)
        if args.mode == "subgroup_omnibus":
            if args.subgroups is None:
                raise ValueError("--subgroups is required for mode 'subgroup_omnibus'.")
            if not args.out.endswith(".parquet"):
                raise ValueError("subgroup_omnibus streaming output requires --out ending in .parquet.")
            if os.path.exists(args.out):
                raise ValueError(f"Output already exists: {args.out}. Choose a new --out path.")
            parts_dir = f"{args.out}.tmp_parts"
            if os.path.exists(parts_dir) and not args.resume:
                raise ValueError(f"Temporary parts directory already exists: {parts_dir}. "
                                 "Use --resume for this run or a new --out path.")

        logger = SimpleLogger(logfile=args.logfile, verbose=True)
        logger.write(f"device: {_get_device()}")

        # Check raw input IDs before any reindexing or sample selection.
        phenotype_df = load_phenotypes(args.phenotype)
        _check_unique(phenotype_df.columns, "phenotype sample IDs")
        _check_unique(phenotype_df.index, "phenotype IDs")
        covariates_df = load_covariates(args.covariates) if args.covariates else None
        if covariates_df is not None:
            _check_unique(covariates_df.index, "covariate sample IDs")
        subgroup_s = load_subgroups(args.subgroups) if args.subgroups else None

        if args.genotype.endswith(".gz"):
            check_raw_sample_header(args.genotype, 4 if args.genotype.endswith(".bed.gz") else 1,
                                    "genotype")
        genotype_ids = plink_sample_ids(args.genotype)
        if genotype_ids is not None:
            _check_unique(genotype_ids, "genotype sample IDs")
            selected = phenotype_df.columns.intersection(genotype_ids, sort=False)
            if len(selected) == 0:
                raise ValueError("Genotype and phenotype have 0 common sample IDs; check the input files.")
            select_samples = selected.tolist()
            n_genotype_samples = len(genotype_ids)
        else:
            select_samples = phenotype_df.columns.tolist()
            n_genotype_samples = None

        genotype_df, _ = load_genotypes(args.genotype, select_samples=select_samples)
        if n_genotype_samples is None:
            n_genotype_samples = genotype_df.shape[1]
        validate_matrices(genotype_df, phenotype_df, covariates_df)
        common = phenotype_df.columns.intersection(genotype_df.columns, sort=False)
        if len(common) == 0:
            raise ValueError("Genotype and phenotype have 0 common sample IDs; check the input files.")

        logger.write(f"  * n_genotype_samples = {n_genotype_samples}")
        logger.write(f"  * n_phenotype_samples = {phenotype_df.shape[1]}")
        logger.write(f"  * n_common_samples = {len(common)}")
        phenotype_df = phenotype_df.loc[:, common]
        if covariates_df is not None:
            missing_c = common.difference(covariates_df.index, sort=False)
            if len(missing_c):
                raise ValueError(f"{len(missing_c)} analysis samples are missing from covariates; "
                                 f"examples: {_preview(missing_c)}")
            covariates_df = covariates_df.loc[common]
            logger.write(f"  * n_covariates = {covariates_df.shape[1]}")
        if args.mode == "subgroup_omnibus":
            counts = validate_subgroups(subgroup_s, common)
            for lab, n in counts.items():
                logger.write(f"  * {lab} n = {n}")

        validate_alignment(genotype_df, phenotype_df, covariates_df)
        logger.write(f"  * {genotype_df.shape[0]} variants x {genotype_df.shape[1]} loaded genotype samples")
        logger.write(f"  * {phenotype_df.shape[0]} phenotypes x {phenotype_df.shape[1]} analysis samples")

        if args.mode == "global_pc":
            out = run_global_pc(genotype_df, phenotype_df, covariates_df=covariates_df,
                                pval_threshold=args.pval_threshold, maf_threshold=args.maf_threshold,
                                batch_size=args.batch_size, pheno_chunk=args.pheno_chunk,
                                logger=logger)
        else:
            out = run_subgroup_omnibus(genotype_df, phenotype_df, subgroup_s,
                                       pval_threshold=args.pval_threshold, maf_threshold=args.maf_threshold,
                                       batch_size=args.batch_size, pheno_chunk=args.pheno_chunk,
                                       omnibus_p_threshold=args.omnibus_p_threshold,
                                       pmin_threshold=args.pmin_threshold,
                                       prefilter=(not args.no_prefilter), logger=logger,
                                       parts_dir=parts_dir, prefilter_mode=args.prefilter_mode,
                                       resume=args.resume,
                                       input_paths={"phenotype": args.phenotype,
                                                    "subgroups": args.subgroups,
                                                    "genotype_pgen": args.genotype + ".pgen"} if
                                       os.path.exists(args.genotype + ".pgen") else
                                       {"phenotype": args.phenotype,
                                        "subgroups": args.subgroups,
                                        "genotype": args.genotype})
            logger.write("  * merge command: python "
                         f"{Path(__file__).with_name('03_merge_subgroup_parts.py')} "
                         f"--out {args.out}")
            return

        logger.write(f"  * writing {out.shape[0]} pairs to {args.out}")
        write_output(out, args.out)
        logger.write("done.")
    except (ValueError, FileNotFoundError) as exc:
        p.exit(2, f"ERROR: {exc}\n")


if __name__ == "__main__":
    main()
