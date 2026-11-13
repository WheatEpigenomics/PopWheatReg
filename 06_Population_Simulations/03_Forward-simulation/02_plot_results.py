#!/usr/bin/env python3
"""Redraw the original forward-simulation figures from existing merged results."""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
from typing import Sequence, Tuple

import numpy as np
import matplotlib

matplotlib.use("Agg")


MODELS = ["background_only", "random_aglink_selection", "top1_aglink_selection"]


MODEL_LABELS = {
    "background_only": "Background only",
    "random_aglink_selection": "Random agLink selection",
    "top1_aglink_selection": "PCC-selected agLink selection",
}


def model_label(cfg: dict, model: str) -> str:
    if model == "top1_aglink_selection":
        return str(cfg.get("selection_label", MODEL_LABELS[model]))
    return MODEL_LABELS[model]


def mkdir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def mean_ci(values: Sequence[float]) -> Tuple[float, float, float, int]:
    a = np.asarray(values, dtype=float)
    a = a[np.isfinite(a)]
    n = len(a)
    if n == 0:
        return np.nan, np.nan, np.nan, 0
    m = float(np.mean(a))
    if n == 1:
        return m, np.nan, np.nan, 1
    se = float(np.std(a, ddof=1) / math.sqrt(n))
    return m, m - 1.96 * se, m + 1.96 * se, n


def plot(cfg: dict, chrom: str):
    import matplotlib.pyplot as plt
    base = Path(cfg["output_root"]) / chrom
    z = np.load(base / "03_merged" / "merged.npz")
    checkpoints = np.array(cfg["checkpoints"], dtype=int)
    out = mkdir(Path(cfg["plot_output_root"]) / chrom / "05_plots")

    # Fig1b dynamic PCC
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    for mi, model in enumerate(MODELS):
        means, los, his = [], [], []
        for ti in range(len(checkpoints)):
            m, lo, hi, _ = mean_ci(z["dyn_pcc"][:, mi, ti])
            means.append(m); los.append(lo); his.append(hi)
        ax.plot(checkpoints, means, marker="o", label=model_label(cfg, model))
        ax.fill_between(checkpoints, los, his, alpha=0.15)
    ax.axhline(0, ls=":", lw=1)
    ax.set_xlabel("Generation")
    ax.set_ylabel("PCC with real $r^2$")
    ax.set_title(f"{chrom}: simulation–real LD similarity\n(changing complete-pair set; diagnostic)")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout()
    fig.savefig(out / "Fig1b_dynamic_PCC.pdf")
    plt.close(fig)

    # Fig1c fixed common mean r2
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    for mi, model in enumerate(MODELS):
        means, los, his = [], [], []
        for ti in range(len(checkpoints)):
            m, lo, hi, _ = mean_ci(z["fixed_mean"][:, mi, ti])
            means.append(m); los.append(lo); his.append(hi)
        ax.plot(checkpoints, means, marker="o", label=model_label(cfg, model))
        ax.fill_between(checkpoints, los, his, alpha=0.12)
    real_ref = float(np.nanmean(z["fixed_real_mean"]))
    ax.axhline(real_ref, ls="--", lw=1.5, color="0.25", label=f"Real population, same fixed links = {real_ref:.3f}")
    ax.axhline(0, ls=":", lw=1, color="0.7")
    ax.set_xlabel("Generation")
    ax.set_ylabel("Mean $r^2$ on fixed common links")
    ax.set_title(f"{chrom}: mean simulated LD vs real LD\n(fixed common links, G0–G{int(checkpoints[-1])})")
    ax.legend(frameon=False, fontsize=8)
    fig.tight_layout()
    fig.savefig(out / "Fig1c_fixed_mean_r2_vs_real.pdf")
    plt.close(fig)

    # Diagnostics: heterozygosity and burden
    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    for mi, model in enumerate(MODELS):
        ax.plot(checkpoints, np.nanmean(z["heterozygosity"][:, mi, :], axis=0), marker="o", label=model_label(cfg, model))
    ax.set_xlabel("Generation"); ax.set_ylabel("Heterozygosity")
    ax.set_title(f"{chrom}: heterozygosity")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(out / "Fig2_heterozygosity.pdf"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.0, 5.2))
    for mi, model in enumerate(MODELS):
        if model == "background_only":
            continue
        ax.plot(checkpoints, np.nanmean(z["burden"][:, mi, :], axis=0), marker="o", color=f"C{mi}", label=model_label(cfg, model))
    ax.set_xlabel("Generation"); ax.set_ylabel("Mean target disruption")
    ax.set_title(f"{chrom}: selection-target disruption burden")
    ax.legend(frameon=False, fontsize=9)
    fig.tight_layout(); fig.savefig(out / "Fig3_target_disruption.pdf"); plt.close(fig)
    print(f"[{chrom}] plots written")


def combine(cfg: dict):
    import matplotlib.pyplot as plt
    root = Path(cfg["output_root"])
    chroms = cfg["chromosomes"]
    out = mkdir(Path(cfg["plot_output_root"]) / "combined")
    # Two chromosomes per row, each with Fig1b and Fig1c. Layout expands to all
    # configured chromosomes (21 wheat chromosomes in production).
    nrows = int(math.ceil(len(chroms) / 2.0))
    fig, axes = plt.subplots(nrows, 4, figsize=(18, 3.8 * nrows), squeeze=False)
    endpoints = []
    for ci, chrom in enumerate(chroms):
        row = ci // 2
        col0 = (ci % 2) * 2
        axb, axc = axes[row, col0], axes[row, col0 + 1]
        z = np.load(root / chrom / "03_merged" / "merged.npz")
        cps = np.array(cfg["checkpoints"], dtype=int)
        for mi, model in enumerate(MODELS):
            means = []; los = []; his = []
            for ti in range(len(cps)):
                m, lo, hi, _ = mean_ci(z["dyn_pcc"][:, mi, ti])
                means.append(m); los.append(lo); his.append(hi)
            axb.plot(cps, means, marker="o", ms=3, label=model_label(cfg, model))
            axb.fill_between(cps, los, his, alpha=0.12)
            fmeans = []; flos = []; fhis = []
            for ti in range(len(cps)):
                m, lo, hi, _ = mean_ci(z["fixed_mean"][:, mi, ti])
                fmeans.append(m); flos.append(lo); fhis.append(hi)
            axc.plot(cps, fmeans, marker="o", ms=3, label=model_label(cfg, model))
            axc.fill_between(cps, flos, fhis, alpha=0.10)
            endpoints.append({
                "chrom": chrom, "model": model, "generation": int(cps[-1]),
                "dynamic_pcc_mean": means[-1], "fixed_mean_r2": fmeans[-1],
                "real_fixed_mean_r2": float(np.nanmean(z["fixed_real_mean"])),
            })
        axb.axhline(0, ls=":", lw=0.8, color="0.5")
        axb.set_title(f"{chrom} — Fig1b dynamic PCC")
        axb.set_xlabel("Generation"); axb.set_ylabel("PCC with real $r^2$")
        real_ref = float(np.nanmean(z["fixed_real_mean"]))
        axc.axhline(real_ref, ls="--", lw=1.2, color="0.25")
        axc.set_title(f"{chrom} — Fig1c fixed mean $r^2$")
        axc.set_xlabel("Generation"); axc.set_ylabel("Mean $r^2$")
    for ci in range(len(chroms), nrows * 2):
        row = ci // 2
        col0 = (ci % 2) * 2
        axes[row, col0].axis("off")
        axes[row, col0 + 1].axis("off")
    handles, labels = axes[0, 0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 0.995))
    fig.suptitle(f"{len(chroms)}-chromosome forward simulation: {cfg.get('analysis_name', 'PCC-selected agLinks')}", y=1.02, fontsize=16)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    combined_pdf = out / "Combined_Fig1bc_all_chromosomes.pdf"
    fig.savefig(combined_pdf, bbox_inches="tight")
    plt.close(fig)

    with open(out / "all_chromosome_endpoint_summary.tsv", "w", newline="") as f:
        fields = list(endpoints[0].keys())
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader(); w.writerows(endpoints)

    # Equal-weight mean across all configured chromosomes at each generation.
    # Replicates are averaged within chromosome first, so chromosomes with more
    # complete links do not receive extra weight in this summary.
    cps = np.array(cfg["checkpoints"], dtype=int)
    chrom_dyn = np.full((len(chroms), len(MODELS), len(cps)), np.nan)
    chrom_fixed = np.full_like(chrom_dyn, np.nan)
    chrom_real = np.full(len(chroms), np.nan)
    for ci, chrom in enumerate(chroms):
        z = np.load(root / chrom / "03_merged" / "merged.npz")
        chrom_dyn[ci] = np.nanmean(z["dyn_pcc"], axis=0)
        chrom_fixed[ci] = np.nanmean(z["fixed_mean"], axis=0)
        chrom_real[ci] = float(np.nanmean(z["fixed_real_mean"]))

    fig, (axb, axc) = plt.subplots(1, 2, figsize=(13.5, 5.2))
    rows = []
    real_mean, real_lo, real_hi, real_n = mean_ci(chrom_real)
    for mi, model in enumerate(MODELS):
        dmean, dlo, dhi = [], [], []
        fmean, flo, fhi = [], [], []
        for ti, gen in enumerate(cps):
            dm, dl, dh, dn = mean_ci(chrom_dyn[:, mi, ti])
            fm, fl, fh, fn = mean_ci(chrom_fixed[:, mi, ti])
            dmean.append(dm); dlo.append(dl); dhi.append(dh)
            fmean.append(fm); flo.append(fl); fhi.append(fh)
            rows.append({
                "generation": int(gen), "model": model,
                "dynamic_pcc_chrom_mean": dm,
                "dynamic_pcc_ci_low": dl, "dynamic_pcc_ci_high": dh,
                "dynamic_pcc_n_chromosomes": dn,
                "fixed_mean_r2_chrom_mean": fm,
                "fixed_mean_r2_ci_low": fl, "fixed_mean_r2_ci_high": fh,
                "fixed_mean_r2_n_chromosomes": fn,
                "real_fixed_mean_r2_chrom_mean": real_mean,
                "real_fixed_mean_r2_n_chromosomes": real_n,
            })
        axb.plot(cps, dmean, marker="o", ms=4, label=model_label(cfg, model))
        axb.fill_between(cps, dlo, dhi, alpha=0.15)
        axc.plot(cps, fmean, marker="o", ms=4, label=model_label(cfg, model))
        axc.fill_between(cps, flo, fhi, alpha=0.12)
    axb.axhline(0, ls=":", lw=1, color="0.5")
    axb.set_title(f"Fig1b: mean across {len(chroms)} chromosomes")
    axb.set_xlabel("Generation"); axb.set_ylabel("Mean PCC with real $r^2$")
    axc.axhline(real_mean, ls="--", lw=1.5, color="0.25",
                label=f"Mean real population $r^2$ = {real_mean:.3f}")
    axc.fill_between(cps, real_lo, real_hi, color="0.25", alpha=0.08)
    axc.set_title(f"Fig1c: mean across {len(chroms)} chromosomes")
    axc.set_xlabel("Generation"); axc.set_ylabel("Mean $r^2$ on fixed common links")
    handles, labels = axb.get_legend_handles_labels()
    handles2, labels2 = axc.get_legend_handles_labels()
    fig.legend(handles + handles2[-1:], labels + labels2[-1:], loc="upper center",
               ncol=4, frameon=False, bbox_to_anchor=(0.5, 1.02))
    fig.suptitle(f"Forward simulation ({cfg.get('analysis_name', 'PCC-selected agLinks')}): chromosome mean by generation", y=1.08, fontsize=15)
    fig.tight_layout()
    mean_pdf = out / "Mean_Fig1bc_all_chromosomes_by_generation.pdf"
    fig.savefig(mean_pdf, bbox_inches="tight")
    plt.close(fig)

    mean_tsv = out / "all_chromosome_generation_mean.tsv"
    with open(mean_tsv, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)
    print(f"Combined figure: {combined_pdf}")
    print(f"Chromosome-mean figure: {mean_pdf}")


def validate_results(cfg: dict, chrom: str):
    """Check saved metadata before assigning labels or generations to results."""
    base = Path(cfg["output_root"]) / chrom
    with (base / "01_prepare" / "prepare_qc.json").open() as f:
        qc = json.load(f)
    for key in ("selection_scope", "top_pcc_mode", "per_block_top_fraction"):
        if qc.get(key) != cfg[key]:
            raise ValueError(f"{chrom}: {key} in saved results does not match config")
    if qc.get("per_block_top_count") is not None:
        raise ValueError(f"{chrom}: saved results use a fixed link count, not top 1%")
    with (base / "03_merged" / "metrics_summary.tsv").open() as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    expected = {(chrom, int(g), model) for g in cfg["checkpoints"] for model in MODELS}
    actual = [(r["chrom"], int(r["generation"]), r["model"]) for r in rows]
    if len(actual) != len(expected) or set(actual) != expected:
        raise ValueError(f"{chrom}: saved generations/models do not match config")
    shape = (int(cfg["replicates"]), len(MODELS), len(cfg["checkpoints"]))
    with np.load(base / "03_merged" / "merged.npz", allow_pickle=False) as z:
        for name in ("dyn_pcc", "fixed_mean", "heterozygosity", "burden"):
            if z[name].shape != shape:
                raise ValueError(f"{chrom}: {name} has shape {z[name].shape}; expected {shape}")
        if z["fixed_real_mean"].shape != (shape[0],):
            raise ValueError(f"{chrom}: fixed_real_mean does not match replicate count")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--input", type=Path, help="Existing category result directory")
    parser.add_argument("--output", type=Path, help="Destination for plots and combined tables; defaults to input")
    parser.add_argument("--chrom", help="Redraw one chromosome only, without combined figures")
    args = parser.parse_args()
    config_path = args.config.resolve()
    with config_path.open() as f:
        cfg = json.load(f)
    root = (config_path.parent / Path(cfg["output_root"]).expanduser()).resolve()
    if os.environ.get("OUTPUT_DIR"):
        root = Path(os.environ["OUTPUT_DIR"]).expanduser().resolve() / root.name
    if args.input is not None:
        root = args.input.expanduser().resolve()
    cfg["output_root"] = str(root)
    cfg["plot_output_root"] = str(args.output.expanduser().resolve() if args.output is not None else root)
    if (cfg.get("selection_scope") != "acr_block" or cfg.get("per_block_top_count") is not None
            or cfg.get("per_block_top_fraction") != 0.01):
        parser.error("only per-block top 1% selection is supported")
    if not cfg["chromosomes"] or len(set(cfg["chromosomes"])) != len(cfg["chromosomes"]):
        parser.error("config chromosomes must be non-empty and unique")
    cps = cfg["checkpoints"]
    if (not cps or any(not isinstance(g, int) or isinstance(g, bool) for g in cps)
            or cps != sorted(set(cps)) or cps[0] != 0 or cps[-1] > cfg["generations"]):
        parser.error("checkpoints must be increasing unique integers from 0 through generations")
    if args.chrom is not None and args.chrom not in cfg["chromosomes"]:
        parser.error(f"chromosome {args.chrom} is not in config chromosomes")
    chroms = [args.chrom] if args.chrom is not None else cfg["chromosomes"]
    for chrom in chroms:
        merged = root / chrom / "03_merged" / "merged.npz"
        if not merged.is_file():
            parser.error(f"merged result not found: {merged}")
        try:
            validate_results(cfg, chrom)
        except (OSError, ValueError, KeyError) as exc:
            parser.error(f"result validation failed: {exc}")
    for chrom in chroms:
        plot(cfg, chrom)
    if args.chrom is None:
        combine(cfg)


if __name__ == "__main__":
    main()
