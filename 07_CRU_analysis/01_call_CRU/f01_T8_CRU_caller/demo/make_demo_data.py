#!/usr/bin/env python3
"""Make a toy one-chromosome dataset so the caller can be tried in a few seconds.

The genome is 20 Mb with three gene-dense islands; ACRs are 500 bp wide and
every agLink runs from one ACR to a gene within 150 kb. Numbers are arbitrary:
the demo only shows the input formats and that the caller runs end to end.
"""
import argparse
from pathlib import Path

import numpy as np

ap = argparse.ArgumentParser(description=__doc__)
ap.add_argument("--outdir", type=Path, default=Path(__file__).resolve().parent / "data")
ap.add_argument("--chrom", default="1A")
ap.add_argument("--size", type=int, default=20_000_000)
ap.add_argument("--islands", default="3000000-8000000,11000000-15000000,16500000-19500000")
ap.add_argument("--seed", type=int, default=20240623)
a = ap.parse_args()

rng = np.random.default_rng(a.seed)
a.outdir.mkdir(parents=True, exist_ok=True)
islands = [tuple(int(x) for x in part.split("-")) for part in a.islands.split(",")]

# ACRs: one every 25 kb inside the islands, sparse outside
acrs = []
for lo, hi in islands:
    for start in range(lo, hi, 25_000):
        acrs.append((start + int(rng.integers(0, 5_000)), 500))
for start in range(500_000, a.size - 500_000, 250_000):
    if not any(lo <= start < hi for lo, hi in islands):
        acrs.append((start + int(rng.integers(0, 5_000)), 500))
acrs.sort()
acr_id = [f"acr{i:07d}:{a.chrom}_{s}_{s + w}" for i, (s, w) in enumerate(acrs, 1)]

# genes: two per 100 kb inside islands, one per 400 kb outside
genes = []
for lo, hi in islands:
    for start in range(lo, hi, 50_000):
        genes.append(start + int(rng.integers(0, 20_000)))
for start in range(400_000, a.size - 400_000, 400_000):
    if not any(lo <= start < hi for lo, hi in islands):
        genes.append(start)
genes.sort()
gene_id = [f"TraesCS1A03G{i:07d}" for i in range(len(genes))]

with (a.outdir / "genome.sizes").open("w") as f:
    f.write(f"{a.chrom}\t{a.size}\n")

with (a.outdir / "genes.bed").open("w") as f:
    for gid, tss in zip(gene_id, genes):
        strand = "+" if rng.random() < .5 else "-"
        start = tss if strand == "+" else max(0, tss - 3_000)
        f.write(f"{a.chrom}\t{start}\t{start + 3_000}\t{gid}\t.\t{strand}\n")

with (a.outdir / "acrs.bed").open("w") as f:
    for aid, (start, width) in zip(acr_id, acrs):
        f.write(f"{a.chrom}\t{start}\t{start + width}\t{aid}\n")

# agLinks: link each ACR to up to three genes within 150 kb
rows = []
gene_starts = np.asarray(genes)
for aid, (start, width) in zip(acr_id, acrs):
    mid = start + width // 2
    lo = np.searchsorted(gene_starts, mid - 150_000)
    hi = np.searchsorted(gene_starts, mid + 150_000)
    if hi <= lo:
        continue
    for j in rng.choice(np.arange(lo, hi), size=min(3, hi - lo), replace=False):
        rows.append(f"{gene_id[j]}\t{aid}\t{rng.uniform(.3, .95):.4f}\t{rng.uniform(1e-12, 1e-4):.3e}\n")
with (a.outdir / "aglink.tsv").open("w") as f:
    f.writelines(rows)

# haplotype blocks and a small co-accessibility table (plotting only)
with (a.outdir / "blocks.bed").open("w") as f:
    for start in range(0, a.size, 2_000_000):
        f.write(f"{a.chrom}\t{start}\t{min(a.size, start + 2_000_000)}\n")

with (a.outdir / "coacr.tsv").open("w") as f:
    f.write("id1\tid2\tcorrelation\n")
    for i in range(0, min(400, len(acr_id) - 2), 2):
        f.write(f"{acr_id[i]}\t{acr_id[i + 2]}\t{rng.uniform(.6, .95):.4f}\n")

print(f"wrote {a.outdir}: {len(acr_id)} ACRs, {len(gene_id)} genes, {len(rows)} agLinks")
