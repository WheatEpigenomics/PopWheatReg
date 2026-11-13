#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ACR–ACR co-accessibility (block + window) — GPU version

Features
- Per-chromosome processing, sliding window (±window around each ACR center)
- Pearson/Spearman, row-wise standardization, dot-product correlation
- Streaming write with |r| >= cutoff (reduces disk I/O)
- VRAM-aware batching: auto estimate pairs-per-chunk from --gpu-mem (or auto)
- De-duplicate: only compute pairs where i < j
"""

import argparse
import math
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

# ---------------------------
# Utilities
# ---------------------------

def parse_gpu_mem(mem_str):
    """Return user limit in bytes or None if 'auto'."""
    s = str(mem_str).strip().lower()
    if s in ("auto", "", "none"):
        return None
    if s.endswith("gb"):
        return int(float(s[:-2]) * (1024**3))
    if s.endswith("mb"):
        return int(float(s[:-2]) * (1024**2))
    # plain number -> treat as GB for convenience
    try:
        return int(float(s) * (1024**3))
    except:
        raise ValueError(f"Unrecognized --gpu-mem: {mem_str}")

def get_gpu_budget_bytes(device, user_limit_bytes=None, safety_ratio=0.7):
    """Get safe VRAM budget (bytes) using free*ratio, capped by user limit if provided."""
    if device.type != "cuda":
        raise RuntimeError("CUDA device not available; please set --device correctly or install CUDA.")
    free_b, total_b = torch.cuda.mem_get_info(device)
    usable = int(free_b * safety_ratio)
    if user_limit_bytes is not None:
        usable = min(usable, user_limit_bytes)
    return max(usable, 16 * 1024 * 1024)  # at least 16MB

def estimate_pairs_per_chunk(n_samples, vram_budget_bytes, safety_ratio=0.6, hard_cap=5_000_000):
    """
    Estimate how many PAIRS we can process in one chunk.

    Memory per pair ≈ (2 * n_samples floats + 1 float result).
    float32 = 4 bytes.
    Use extra safety_ratio because we also hold Z_chunk etc.
    """
    budget = int(vram_budget_bytes * safety_ratio)
    bytes_per_pair = (2 * n_samples + 1) * 4
    P = max(1, budget // bytes_per_pair)
    return min(P, hard_cap)

def rowwise_rank_torch(x):
    """
    Row-wise rank (0..S-1). Tie-handling: ordinal ranks (fast).
    For stricter Spearman with average ranks on ties, replace with a slower routine if needed.
    """
    return x.argsort(dim=1).argsort(dim=1).float()

def normalize_rows(x, eps=1e-8):
    mu = x.mean(dim=1, keepdim=True)
    sd = x.std(dim=1, keepdim=True)
    return (x - mu) / (sd + eps), sd.squeeze(1)

# ---------------------------
# Main
# ---------------------------

def main():
    ap = argparse.ArgumentParser(description="GPU ACR–ACR co-accessibility (block + window)")
    ap.add_argument("--matrix", required=True,
                    help="ACR matrix (rows=ACR IDs, cols=samples). Whitespace- or tab-separated.")
    ap.add_argument("--bed", required=True,
                    help="BED-like: chr start end id (tab-separated). IDs must match matrix row index.")
    ap.add_argument("--window", type=int, default=1_000_000,
                    help="Window size in bp (default: 1,000,000).")
    ap.add_argument("--method", choices=["pearson","spearman"], default="pearson",
                    help="Correlation method.")
    ap.add_argument("--cutoff", type=float, default=0.6,
                    help="Absolute correlation cutoff to write.")
    ap.add_argument("--min-std", type=float, default=1e-8,
                    help="Drop rows whose per-row std < min-std to avoid NaNs (after ranking if spearman).")
    ap.add_argument("--gpu-mem", default="auto",
                    help='Target VRAM (e.g. "20GB", "20480MB", number=GB). Default "auto" uses ~70%% of free.')
    ap.add_argument("--device", default="cuda:0", help='CUDA device, e.g. "cuda:0".')
    ap.add_argument("--max-pairs-per-chunk", type=int, default=0,
                    help="Optional hard cap for pairs per GPU chunk (override auto). 0 = auto.")
    ap.add_argument("--output", required=True, help="Output TSV.")
    args = ap.parse_args()

    # Device
    device = torch.device(args.device if torch.cuda.is_available() else "cpu")
    if device.type != "cuda":
        raise RuntimeError("No CUDA device available, but GPU version requested. Please use a CUDA device or CPU script.")

    # Read matrix (float32) — support space or tab
    print(f"[INFO] Reading matrix: {args.matrix}")
    df = pd.read_csv(args.matrix, sep=r"\s+", index_col=0)
    df = df.astype(np.float32)
    n_rows, n_samples = df.shape
    if n_samples == 0:
        raise ValueError("Matrix has 0 columns. Check delimiter/header. Use whitespace or tab with header row.")
    print(f"[INFO] Matrix loaded: {n_rows:,} ACRs x {n_samples} samples")

    # ID -> row index
    id2row = {rid: i for i, rid in enumerate(df.index.tolist())}

    # Read BED
    print(f"[INFO] Reading BED: {args.bed}")
    bed = pd.read_csv(args.bed, sep="\t", header=None, names=["chr","start","end","id"])
    # Keep only IDs that exist in matrix
    bed = bed[bed["id"].isin(id2row)].copy()
    if bed.empty:
        raise ValueError("No BED ids matched matrix row IDs. Check BED columns order (chr start end id).")
    # center position
    bed["center"] = ((bed["start"] + bed["end"]) // 2).astype(np.int64)

    # VRAM-based batching
    user_limit_bytes = parse_gpu_mem(args.gpu_mem)
    vram_budget_bytes = get_gpu_budget_bytes(device, user_limit_bytes=user_limit_bytes, safety_ratio=0.7)
    auto_pairs = estimate_pairs_per_chunk(n_samples, vram_budget_bytes, safety_ratio=0.6, hard_cap=5_000_000)
    if args.max_pairs_per_chunk and args.max_pairs_per_chunk > 0:
        pairs_per_chunk = min(auto_pairs, args.max_pairs_per_chunk)
    else:
        pairs_per_chunk = auto_pairs
    print(f"[INFO] VRAM budget ≈ {vram_budget_bytes/1024**3:.2f} GB; pairs-per-chunk ≈ {pairs_per_chunk:,}")

    # Open output
    out = open(args.output, "w")
    out.write("id1\tid2\tcorrelation\n")

    # Process per chromosome
    for chrom, group in bed.groupby("chr", sort=False):
        grp = group.sort_values("center").reset_index(drop=True)
        ids = grp["id"].tolist()
        centers = grp["center"].to_numpy()

        # Map to matrix row indices
        row_idx = np.fromiter((id2row[i] for i in ids), dtype=np.int64)
        n_chr = len(row_idx)
        if n_chr < 2:
            continue

        print(f"[INFO] Chromosome {chrom}: {n_chr:,} ACRs")

        # Pull this chromosome block to GPU
        Xc_cpu = df.iloc[row_idx].to_numpy(dtype=np.float32, copy=False)
        Xc = torch.from_numpy(Xc_cpu).to(device, non_blocking=True)

        # Spearman: rank per row first
        if args.method == "spearman":
            Xc = rowwise_rank_torch(Xc)

        # Drop near-constant rows (std < min-std)
        Zc, stds = normalize_rows(Xc)
        valid_mask = stds >= args.min_std
        if valid_mask.sum().item() < 2:
            print(f"[WARN] Chromosome {chrom}: <2 valid rows after std filter; skipping.")
            del Xc, Zc
            torch.cuda.empty_cache()
            continue

        if valid_mask.sum().item() < len(valid_mask):
            # Filter IDs/centers accordingly
            valid_idx = valid_mask.nonzero(as_tuple=False).squeeze(1).cpu().numpy()
            Zc = Zc[valid_idx]
            ids = [ids[i] for i in valid_idx]
            centers = centers[valid_idx]
            n_chr = len(ids)
            print(f"[INFO] Chromosome {chrom}: {len(valid_idx)} rows kept after std filter")

        S = Zc.shape[1]
        assert S == n_samples, "Sample dimension mismatch"

        # Sliding window with i<j; batch pairs up to pairs_per_chunk
        r = 0
        i = 0
        # Progress bar over anchors (i)
        pbar = tqdm(total=n_chr, desc=f"{chrom}", unit="anchor", leave=False)
        buf_lines = []
        while i < n_chr:
            # expand right pointer
            if r < i + 1:
                r = i + 1
            while r < n_chr and (centers[r] - centers[i]) <= args.window:
                r += 1
            j_start, j_end = i + 1, r  # [j_start, j_end)

            # process pairs (i, j_start..j_end-1) in chunks
            if j_end > j_start:
                A = Zc[i:i+1]  # 1 x S
                j = j_start
                while j < j_end:
                    k = min(pairs_per_chunk, j_end - j)
                    B = Zc[j:j+k]  # k x S
                    # corr = (A*B).sum(1)/S
                    rvec = torch.sum(A * B, dim=1) / float(S)
                    if args.cutoff > 0.0:
                        mask = torch.abs(rvec) >= args.cutoff
                    else:
                        mask = torch.ones_like(rvec, dtype=torch.bool)

                    if mask.any():
                        sel = torch.nonzero(mask, as_tuple=False).squeeze(1).tolist()
                        if sel:
                            # collect lines, then write in one go
                            for off in sel:
                                corr_val = rvec[off].item()
                                buf_lines.append(f"{ids[i]}\t{ids[j+off]}\t{corr_val:.6f}\n")
                            if len(buf_lines) >= 100000:
                                out.writelines(buf_lines)
                                buf_lines.clear()
                    j += k

            i += 1
            pbar.update(1)

        if buf_lines:
            out.writelines(buf_lines)
            buf_lines.clear()
        pbar.close()

        # free GPU memory
        del Xc, Zc
        torch.cuda.empty_cache()

    out.close()
    print(f"[INFO] Done. Results written to {args.output}")

if __name__ == "__main__":
    torch.set_grad_enabled(False)
    main()
