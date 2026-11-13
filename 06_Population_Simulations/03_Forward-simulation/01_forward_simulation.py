#!/usr/bin/env python3
from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
import os
import random
import re
import shutil
from collections import defaultdict
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np


MODELS = ["background_only", "random_aglink_selection", "top1_aglink_selection"]


def selected_count_for_block(cfg: dict, eligible_n: int) -> int:
    """Keep the top 1% of eligible links in each non-empty ACR block."""
    if eligible_n == 0:
        return 0
    return max(1, int(math.ceil(eligible_n * float(cfg["per_block_top_fraction"]))))


def select_within_acr_blocks(cfg: dict, eligible: Sequence[PairRec]) -> Tuple[List[PairRec], Dict[str, int], Dict[str, int]]:
    """Rank and select independently inside each source ACR haplotype block."""
    mode = cfg.get("top_pcc_mode", "raw")
    groups: Dict[str, List[PairRec]] = defaultdict(list)
    for rec in eligible:
        groups[rec.acr_block].append(rec)
    selected: List[PairRec] = []
    eligible_counts: Dict[str, int] = {}
    selected_counts: Dict[str, int] = {}
    for block_id in sorted(groups):
        rows = groups[block_id]
        if mode == "abs":
            rows.sort(key=lambda r: (abs(r.pcc), r.pcc), reverse=True)
        else:
            rows.sort(key=lambda r: r.pcc, reverse=True)
        k = selected_count_for_block(cfg, len(rows))
        eligible_counts[block_id] = len(rows)
        selected_counts[block_id] = k
        selected.extend(rows[:k])
    # Stable global ordering is only for deterministic output identifiers.
    if mode == "abs":
        selected.sort(key=lambda r: (abs(r.pcc), r.pcc, r.acr_block), reverse=True)
    else:
        selected.sort(key=lambda r: (r.pcc, r.acr_block), reverse=True)
    return selected, eligible_counts, selected_counts


def norm_chrom(s: str) -> str:
    s = s.strip()
    if s.lower().startswith("chr"):
        s = s[3:]
    return s


def mkdir(p: Path) -> Path:
    p.mkdir(parents=True, exist_ok=True)
    return p


def load_json(path: Path) -> dict:
    with open(path) as f:
        return json.load(f)


def dump_json(obj, path: Path):
    mkdir(path.parent)
    with open(path, "w") as f:
        json.dump(obj, f, indent=2, sort_keys=True)


def stable_int_seed(*parts, modulo=2**32 - 1) -> int:
    h = hashlib.sha256("|".join(map(str, parts)).encode()).hexdigest()
    return int(h[:16], 16) % modulo


def overlap(a0: int, a1: int, b0: int, b1: int) -> bool:
    return a0 < b1 and b0 < a1


def pearson_safe(x, y) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    if m.sum() < 3:
        return float("nan")
    xx, yy = x[m], y[m]
    if np.nanstd(xx) <= 0 or np.nanstd(yy) <= 0:
        return float("nan")
    return float(np.corrcoef(xx, yy)[0, 1])


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


@dataclass
class ACR:
    chrom: str
    start: int
    end: int
    acr_id: str

    @property
    def center(self) -> int:
        return (self.start + self.end) // 2


@dataclass
class Gene:
    chrom: str
    start: int  # 0-based half-open
    end: int
    strand: str
    gene_id: str

    @property
    def center(self) -> int:
        return (self.start + self.end) // 2

    @property
    def tss(self) -> int:
        return self.start if self.strand == "+" else self.end - 1

    def promoter(self, upstream_bp: int) -> Tuple[int, int]:
        if self.strand == "+":
            return max(0, self.start - upstream_bp), self.start
        return self.end, self.end + upstream_bp


@dataclass
class Block:
    chrom: str
    start: int
    end: int
    block_id: str


@dataclass
class AgLink:
    gene_id: str
    acr_id: str
    pcc: float
    pvalue: float


@dataclass
class PairRec:
    pair_id: str
    gene_id: str
    acr_id: str
    pcc: float
    acr_center: int
    gene_center: int
    distance_bp: int
    acr_block: str
    gene_block: str
    acr_snp_pos: int
    gene_snp_pos: int
    target_parity: int
    real_r2: float = float("nan")


class BlockIndex:
    def __init__(self, blocks: Sequence[Block]):
        self.blocks = sorted(blocks, key=lambda b: (b.start, b.end))
        self.starts = [b.start for b in self.blocks]

    def containing(self, pos: int) -> Optional[Block]:
        i = bisect.bisect_right(self.starts, pos) - 1
        if i >= 0:
            b = self.blocks[i]
            if b.start <= pos < b.end:
                return b
        return None


def parse_attrs(s: str) -> Dict[str, str]:
    out = {}
    for item in s.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
    return out


def gene_aliases(gid: str) -> List[str]:
    vals = {gid}
    if gid.startswith("gene:"):
        vals.add(gid[5:])
    else:
        vals.add("gene:" + gid)
    # Also support accidental versioned IDs while preserving full ID.
    vals.add(re.sub(r"\.\d+$", "", gid))
    return [x for x in vals if x]


def read_genes(gff: Path, wanted: set[str]) -> Dict[str, Dict[str, Gene]]:
    genes = {c: {} for c in wanted}
    with open(gff, errors="replace") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != "gene":
                continue
            c = norm_chrom(p[0])
            if c not in wanted:
                continue
            try:
                start = int(p[3]) - 1
                end = int(p[4])
            except ValueError:
                continue
            attrs = parse_attrs(p[8])
            gid = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name")
            if not gid:
                continue
            g = Gene(c, start, end, p[6], gid.replace("gene:", "", 1) if gid.startswith("gene:") else gid)
            for alias in gene_aliases(gid):
                genes[c][alias] = g
            for k in ("Name", "gene_id"):
                if attrs.get(k):
                    for alias in gene_aliases(attrs[k]):
                        genes[c][alias] = g
    return genes


def read_acrs(path: Path, wanted: set[str]) -> Dict[str, Dict[str, ACR]]:
    out = {c: {} for c in wanted}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 4:
                continue
            c = norm_chrom(p[0])
            if c not in wanted:
                continue
            try:
                a = ACR(c, int(p[1]), int(p[2]), p[3])
            except ValueError:
                continue
            out[c][a.acr_id] = a
    return out


def read_blocks(path: Path, wanted: set[str]) -> Dict[str, List[Block]]:
    out = {c: [] for c in wanted}
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            c = norm_chrom(p[0])
            if c not in wanted:
                continue
            try:
                start, end = int(p[1]), int(p[2])
            except ValueError:
                continue
            bid = p[3] if len(p) > 3 else f"{c}:{start}-{end}"
            out[c].append(Block(c, start, end, bid))
    return out


def read_aglinks(path: Path, acr_by_chr: Dict[str, Dict[str, ACR]], genes_by_chr: Dict[str, Dict[str, Gene]], wanted: set[str]) -> Dict[str, List[AgLink]]:
    out = {c: [] for c in wanted}
    # Determine chromosome by resolving ACR first, then gene.
    acr_to_chr = {}
    for c, d in acr_by_chr.items():
        acr_to_chr.update({k: c for k in d})
    with open(path) as f:
        for line in f:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 3:
                continue
            try:
                pcc = float(p[2])
                pv = float(p[3]) if len(p) > 3 else float("nan")
            except ValueError:
                continue  # header
            gid, aid = p[0], p[1]
            c = acr_to_chr.get(aid)
            if c not in wanted:
                continue
            g = genes_by_chr[c].get(gid) or genes_by_chr[c].get(gid.replace("gene:", "", 1))
            if g is None:
                # Keep unresolved out of pool; QC counts resolution later.
                continue
            out[c].append(AgLink(g.gene_id, aid, pcc, pv))
    return out


def find_vcf(vcf_dir: Path, chrom: str) -> Path:
    pats = [
        f"*{chrom}*.vcf.gz",
        f"*{chrom}*.vcf",
    ]
    candidates = []
    for pat in pats:
        candidates.extend(vcf_dir.glob(pat))
    # Avoid accidental index files and prefer exact-looking names.
    candidates = sorted({p for p in candidates if not str(p).endswith((".tbi", ".csi"))})
    if not candidates:
        raise FileNotFoundError(f"No VCF found for {chrom} under {vcf_dir}")
    scored = []
    for p in candidates:
        name = p.name
        score = 0
        if f".{chrom}." in name:
            score += 5
        if "SNP" in name:
            score += 2
        if "hex203" in name:
            score += 2
        if "maf002" in name:
            score += 1
        if "miss02" in name:
            score += 1
        scored.append((score, len(name), name, p))
    scored.sort(key=lambda x: (-x[0], x[1], x[2]))
    return scored[0][3]


def open_text_auto(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return open(path)


def scan_vcf_positions(path: Path, chrom: str) -> np.ndarray:
    vals = []
    with open_text_auto(path) as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            p = line.split("\t", 3)
            if len(p) < 2:
                continue
            if norm_chrom(p[0]) != chrom:
                continue
            try:
                vals.append(int(p[1]))
            except ValueError:
                pass
    if not vals:
        raise RuntimeError(f"No SNP positions found in {path} for {chrom}")
    return np.array(sorted(set(vals)), dtype=np.int64)


def parse_dosage(fmt: str, samples: Sequence[str]) -> np.ndarray:
    keys = fmt.split(":")
    idx_gt = keys.index("GT") if "GT" in keys else None
    idx_pgt = keys.index("PGT") if "PGT" in keys else None
    out = np.full(len(samples), np.nan, dtype=float)
    for i, s in enumerate(samples):
        vals = s.split(":")
        token = None
        if idx_gt is not None and idx_gt < len(vals):
            g = vals[idx_gt]
            if g not in (".", "./.", ".|."):
                token = g
        if token is None and idx_pgt is not None and idx_pgt < len(vals):
            g = vals[idx_pgt]
            if g not in (".", "./.", ".|."):
                token = g
        if token is None:
            continue
        a = re.split(r"[|/]", token)
        if len(a) != 2 or any(x == "." for x in a):
            continue
        try:
            aa = [int(x) for x in a]
        except ValueError:
            continue
        if any(x not in (0, 1) for x in aa):
            continue
        out[i] = aa[0] + aa[1]
    return out


def load_vcf_genotypes(path: Path, chrom: str, positions: set[int]) -> Tuple[Dict[int, np.ndarray], List[str]]:
    geno: Dict[int, np.ndarray] = {}
    samples: List[str] = []
    with open_text_auto(path) as f:
        for line in f:
            if line.startswith("#CHROM"):
                p = line.rstrip("\n").split("\t")
                samples = p[9:]
                continue
            if not line or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 10 or norm_chrom(p[0]) != chrom:
                continue
            pos = int(p[1])
            if pos not in positions:
                continue
            # Only biallelic records are supported by this forward model.
            if "," in p[4]:
                continue
            geno[pos] = parse_dosage(p[8], p[9:])
            if len(geno) == len(positions):
                break
    return geno, samples


def nearest_snp_in_block(pos: int, block: Block, snps: np.ndarray, window: int) -> Optional[int]:
    # BED block [start,end) corresponds to 1-based VCF positions start+1 ... end (inclusive).
    # searchsorted upper bound is exclusive, hence end+1.
    lo = max(block.start + 1, pos - window)
    hi = min(block.end + 1, pos + window + 1)
    i = int(np.searchsorted(snps, lo, side="left"))
    j = int(np.searchsorted(snps, hi, side="left"))
    if i >= j:
        return None
    arr = snps[i:j]
    k = int(np.argmin(np.abs(arr - pos)))
    return int(arr[k])


def endpoint_mapping(center0: int, block_index: BlockIndex, snps: np.ndarray, window: int):
    b = block_index.containing(center0)
    if b is None:
        return None
    # convert BED center0 to approximately matching 1-based VCF position
    q = center0 + 1
    s = nearest_snp_in_block(q, b, snps, window)
    if s is None:
        return None
    return b, s


def real_r2(g1: np.ndarray, g2: np.ndarray) -> float:
    m = np.isfinite(g1) & np.isfinite(g2)
    if m.sum() < 5:
        return float("nan")
    x, y = g1[m], g2[m]
    if np.std(x) <= 0 or np.std(y) <= 0:
        return float("nan")
    r = np.corrcoef(x, y)[0, 1]
    return float(r * r)


def af_from_dosage(g: np.ndarray) -> float:
    m = np.isfinite(g)
    if m.sum() == 0:
        return float("nan")
    return float(np.mean(g[m]) / 2.0)


def promoter_overlap(acr: ACR, gene: Gene, upstream: int) -> bool:
    p0, p1 = gene.promoter(upstream)
    return overlap(acr.start, acr.end, p0, p1)


def pair_to_dict(p: PairRec) -> dict:
    d = asdict(p)
    return d


def write_pairs(path: Path, rows: Sequence[PairRec]):
    mkdir(path.parent)
    fields = list(pair_to_dict(rows[0]).keys()) if rows else [
        "pair_id", "gene_id", "acr_id", "pcc", "acr_center", "gene_center", "distance_bp",
        "acr_block", "gene_block", "acr_snp_pos", "gene_snp_pos", "target_parity", "real_r2"
    ]
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(pair_to_dict(r))


def prepare_common(cfg: dict):
    raw = Path(cfg["raw_dir"])
    wanted = set(cfg["chromosomes"])
    fcfg = cfg["files"]
    genes = read_genes(raw / fcfg["gff"], wanted)
    acrs = read_acrs(raw / fcfg["acrs"], wanted)
    blocks = read_blocks(raw / fcfg["blocks"], wanted)
    aglinks = read_aglinks(raw / fcfg["aglinks"], acrs, genes, wanted)
    return genes, acrs, blocks, aglinks


def make_pair_record(a: ACR, g: Gene, pcc: float, bidx: BlockIndex, snps: np.ndarray, window: int, pair_id: str, seed: int) -> Optional[PairRec]:
    ma = endpoint_mapping(a.center, bidx, snps, window)
    mg = endpoint_mapping(g.center, bidx, snps, window)
    if ma is None or mg is None:
        return None
    ba, sa = ma
    bg, sg = mg
    if sa == sg:
        return None
    parity = stable_int_seed(seed, pair_id, a.acr_id, g.gene_id) % 2
    return PairRec(
        pair_id=pair_id,
        gene_id=g.gene_id,
        acr_id=a.acr_id,
        pcc=float(pcc),
        acr_center=a.center,
        gene_center=g.center,
        distance_bp=abs(a.center - g.center),
        acr_block=ba.block_id,
        gene_block=bg.block_id,
        acr_snp_pos=sa,
        gene_snp_pos=sg,
        target_parity=int(parity),
    )


def choose_top1_and_random(cfg: dict, chrom: str, genes: Dict[str, Gene], acrs: Dict[str, ACR], blocks: List[Block], aglinks: List[AgLink], vcf: Path):
    upstream = int(cfg["promoter_upstream_bp"])
    window = int(cfg["representative_snp_window_bp"])
    seed = int(cfg["seed"])
    snps = scan_vcf_positions(vcf, chrom)
    bidx = BlockIndex(blocks)

    known_pairs = {(x.acr_id, x.gene_id) for x in aglinks}
    n_input = len(aglinks)
    n_promoter_masked = 0
    n_missing = 0
    n_mapping_fail = 0
    n_same_snp = 0
    eligible: List[PairRec] = []
    seen_rep_pairs = set()

    for i, al in enumerate(aglinks):
        a = acrs.get(al.acr_id)
        g = genes.get(al.gene_id)
        if a is None or g is None:
            n_missing += 1
            continue
        if promoter_overlap(a, g, upstream):
            n_promoter_masked += 1
            continue
        rec = make_pair_record(a, g, al.pcc, bidx, snps, window, f"true_{i:08d}", seed)
        if rec is None:
            # Distinguish same SNP is not essential for QC here.
            n_mapping_fail += 1
            continue
        key = tuple(sorted((rec.acr_snp_pos, rec.gene_snp_pos)))
        if key in seen_rep_pairs:
            continue
        seen_rep_pairs.add(key)
        eligible.append(rec)

    if not eligible:
        raise RuntimeError(f"{chrom}: no promoter-masked strictly mappable agLinks")

    mode = cfg.get("top_pcc_mode", "raw")
    top, eligible_by_block, selected_by_block = select_within_acr_blocks(cfg, eligible)
    k = len(top)
    for j, r in enumerate(top, 1):
        r.pair_id = f"BLOCKSEL_{chrom}_{j:05d}"
        r.target_parity = stable_int_seed(seed, chrom, "TOP", r.acr_id, r.gene_id) % 2

    sel_acrs = {r.acr_id for r in top}
    sel_genes = {r.gene_id for r in top}
    sel_markers = {p for r in top for p in (r.acr_snp_pos, r.gene_snp_pos)}
    sel_rep_pairs = {tuple(sorted((r.acr_snp_pos, r.gene_snp_pos))) for r in top}

    # Sorted unique genes by center for distance-targeted random generation.
    gene_unique = {}
    for g in genes.values():
        gene_unique[g.gene_id] = g
    gene_list = sorted(gene_unique.values(), key=lambda g: g.center)
    gene_centers = [g.center for g in gene_list]
    acr_list = list(acrs.values())
    rng = random.Random(stable_int_seed(seed, chrom, "random_pairs"))

    random_rows: List[PairRec] = []
    used_random_pairs = set()
    used_random_acrs = set()
    used_random_genes = set()
    used_random_markers = set()
    tolerances = list(map(float, cfg["random_distance_relative_tolerances"]))
    min_abs = int(cfg["random_distance_min_abs_bp"])

    # cache endpoint mappings to avoid repeated block/SNP lookups
    amap_cache = {}
    gmap_cache = {}

    def amap(a: ACR):
        if a.acr_id not in amap_cache:
            amap_cache[a.acr_id] = endpoint_mapping(a.center, bidx, snps, window)
        return amap_cache[a.acr_id]

    def gmap(g: Gene):
        if g.gene_id not in gmap_cache:
            gmap_cache[g.gene_id] = endpoint_mapping(g.center, bidx, snps, window)
        return gmap_cache[g.gene_id]

    for idx, target in enumerate(top, 1):
        target_d = max(1, target.distance_bp)
        found = None
        for tol in tolerances:
            delta = max(min_abs, int(target_d * tol))
            # Increasing attempts at broader tolerance.
            for _ in range(12000):
                a = rng.choice(acr_list)
                if cfg.get("random_disallow_shared_acr", True) and (a.acr_id in sel_acrs or a.acr_id in used_random_acrs):
                    continue
                ma = amap(a)
                if ma is None:
                    continue
                ba, sa = ma
                if cfg.get("random_disallow_shared_marker", True) and (sa in sel_markers or sa in used_random_markers):
                    continue
                sign = -1 if rng.random() < 0.5 else 1
                target_center = a.center + sign * target_d
                lo_c, hi_c = target_center - delta, target_center + delta
                lo = bisect.bisect_left(gene_centers, lo_c)
                hi = bisect.bisect_right(gene_centers, hi_c)
                if lo >= hi:
                    continue
                g = gene_list[rng.randrange(lo, hi)]
                if cfg.get("random_disallow_shared_gene", True) and (g.gene_id in sel_genes or g.gene_id in used_random_genes):
                    continue
                if (a.acr_id, g.gene_id) in known_pairs:
                    continue
                if promoter_overlap(a, g, upstream):
                    continue
                mg = gmap(g)
                if mg is None:
                    continue
                bg, sg = mg
                if sa == sg:
                    continue
                if cfg.get("random_disallow_shared_marker", True) and (sg in sel_markers or sg in used_random_markers):
                    continue
                rp = tuple(sorted((sa, sg)))
                if rp in sel_rep_pairs or rp in used_random_pairs:
                    continue
                d = abs(a.center - g.center)
                if abs(d - target_d) > delta:
                    continue
                rec = PairRec(
                    pair_id=f"RAND_{chrom}_{idx:05d}", gene_id=g.gene_id, acr_id=a.acr_id, pcc=float("nan"),
                    acr_center=a.center, gene_center=g.center, distance_bp=d,
                    acr_block=ba.block_id, gene_block=bg.block_id,
                    acr_snp_pos=sa, gene_snp_pos=sg,
                    target_parity=int(stable_int_seed(seed, chrom, "RAND", a.acr_id, g.gene_id) % 2),
                )
                found = rec
                break
            if found is not None:
                break
        if found is None:
            raise RuntimeError(
                f"{chrom}: could not generate random pair {idx}/{k} with distance similar to {target_d} bp. "
                "Review random overlap constraints or distance tolerances in the selected configuration."
            )
        random_rows.append(found)
        used_random_pairs.add(tuple(sorted((found.acr_snp_pos, found.gene_snp_pos))))
        used_random_acrs.add(found.acr_id)
        used_random_genes.add(found.gene_id)
        used_random_markers.update((found.acr_snp_pos, found.gene_snp_pos))

    all_positions = {p for r in top + random_rows for p in (r.acr_snp_pos, r.gene_snp_pos)}
    geno, samples = load_vcf_genotypes(vcf, chrom, all_positions)
    missing_g = sorted(all_positions - set(geno))
    if missing_g:
        raise RuntimeError(f"{chrom}: {len(missing_g)} selected representative SNPs could not be genotype-loaded: {missing_g[:10]}")

    # Add real r2 to both true and random pairs; random real r2 is diagnostic only.
    for r in top + random_rows:
        r.real_r2 = real_r2(geno[r.acr_snp_pos], geno[r.gene_snp_pos])

    marker_positions = np.array(sorted(all_positions), dtype=np.int64)
    af = np.array([af_from_dosage(geno[int(p)]) for p in marker_positions], dtype=float)
    bad_af = ~np.isfinite(af) | (af <= 0.0) | (af >= 1.0)
    if bad_af.any():
        # MAF-filtered VCF should make this rare. Fail rather than silently alter the target set.
        badp = marker_positions[bad_af]
        raise RuntimeError(f"{chrom}: monomorphic/unusable selected marker(s): {badp[:10].tolist()}")

    pos_to_idx = {int(p): i for i, p in enumerate(marker_positions)}
    true_a = np.array([pos_to_idx[r.acr_snp_pos] for r in top], dtype=np.int32)
    true_b = np.array([pos_to_idx[r.gene_snp_pos] for r in top], dtype=np.int32)
    true_par = np.array([r.target_parity for r in top], dtype=np.int8)
    true_real = np.array([r.real_r2 for r in top], dtype=float)
    rand_a = np.array([pos_to_idx[r.acr_snp_pos] for r in random_rows], dtype=np.int32)
    rand_b = np.array([pos_to_idx[r.gene_snp_pos] for r in random_rows], dtype=np.int32)
    rand_par = np.array([r.target_parity for r in random_rows], dtype=np.int8)
    rand_real = np.array([r.real_r2 for r in random_rows], dtype=float)

    dist_ratio = np.array([r.distance_bp / max(1, t.distance_bp) for r, t in zip(random_rows, top)], dtype=float)
    qc = {
        "version": "v17_raw_standalone",
        "chrom": chrom,
        "vcf": str(vcf),
        "n_real_samples": len(samples),
        "n_input_resolved_aglinks": n_input,
        "n_promoter_masked": n_promoter_masked,
        "n_missing_endpoint": n_missing,
        "n_strict_mapping_fail": n_mapping_fail,
        "n_eligible_after_mask_and_mapping_and_rep_pair_dedup": len(eligible),
        "selection_scope": "acr_block",
        "per_block_top_fraction": float(cfg.get("per_block_top_fraction", 0.0)),
        "per_block_top_count": cfg.get("per_block_top_count"),
        "selection_label": cfg.get("selection_label", "PCC-selected agLink selection"),
        "n_eligible_acr_blocks": len(eligible_by_block),
        "eligible_count_by_acr_block": eligible_by_block,
        "selected_count_by_acr_block": selected_by_block,
        "top_pcc_mode": mode,
        "n_selected_top1": len(top),
        "n_random": len(random_rows),
        "n_unique_markers": len(marker_positions),
        "selected_pcc_min": float(min(r.pcc for r in top)),
        "selected_pcc_max": float(max(r.pcc for r in top)),
        "random_distance_ratio_median": float(np.median(dist_ratio)),
        "random_distance_ratio_min": float(np.min(dist_ratio)),
        "random_distance_ratio_max": float(np.max(dist_ratio)),
        "promoter_overlap_selected": int(sum(promoter_overlap(acrs[r.acr_id], genes[r.gene_id], upstream) for r in top)),
        "promoter_overlap_random": int(sum(promoter_overlap(acrs[r.acr_id], genes[r.gene_id], upstream) for r in random_rows)),
        "random_known_sig_overlap": int(sum((r.acr_id, r.gene_id) in known_pairs for r in random_rows)),
        "random_shared_acr_with_top": int(sum(r.acr_id in sel_acrs for r in random_rows)),
        "random_shared_gene_with_top": int(sum(r.gene_id in sel_genes for r in random_rows)),
        "random_shared_marker_with_top": int(sum((r.acr_snp_pos in sel_markers) or (r.gene_snp_pos in sel_markers) for r in random_rows)),
    }
    arrays = dict(
        marker_positions=marker_positions,
        af=af,
        true_a=true_a, true_b=true_b, true_par=true_par, true_real_r2=true_real,
        rand_a=rand_a, rand_b=rand_b, rand_par=rand_par, rand_real_r2=rand_real,
    )
    return top, random_rows, qc, arrays


def prepare_chrom(cfg: dict, chrom: str, genes_all=None, acrs_all=None, blocks_all=None, aglinks_all=None, force=False):
    if chrom not in cfg["chromosomes"]:
        raise ValueError(f"chromosome {chrom} is not in config chromosomes")
    raw = Path(cfg["raw_dir"])
    outroot = Path(cfg["output_root"]) / chrom
    prep = outroot / "01_prepare"
    if force and prep.exists():
        shutil.rmtree(prep)
    mkdir(prep)
    if genes_all is None:
        genes_all, acrs_all, blocks_all, aglinks_all = prepare_common(cfg)
    genes = genes_all[chrom]
    acrs = acrs_all[chrom]
    blocks = blocks_all[chrom]
    aglinks = aglinks_all[chrom]
    vcf = find_vcf(raw / cfg["files"]["vcf_dir"], chrom)
    top, rnd, qc, arrays = choose_top1_and_random(cfg, chrom, genes, acrs, blocks, aglinks, vcf)
    write_pairs(prep / "selected_top1_aglinks.tsv", top)
    write_pairs(prep / "random_aglinks.tsv", rnd)
    np.savez_compressed(prep / "prepared.npz", **arrays)
    dump_json(qc, prep / "prepare_qc.json")
    print(f"[{chrom}] prepared: eligible={qc['n_eligible_after_mask_and_mapping_and_rep_pair_dedup']}, blocks={qc['n_eligible_acr_blocks']}, selected={len(top)}, random={len(rnd)}, markers={len(arrays['marker_positions'])}")


def prepare_all(cfg: dict, force=False):
    genes, acrs, blocks, aglinks = prepare_common(cfg)
    for chrom in cfg["chromosomes"]:
        prepare_chrom(cfg, chrom, genes, acrs, blocks, aglinks, force=force)


def recomb_probs(positions: np.ndarray, rate: float) -> np.ndarray:
    if len(positions) < 2:
        return np.zeros(0, dtype=float)
    d = np.diff(positions.astype(float))
    return 0.5 * (1.0 - np.exp(-2.0 * rate * d))


def founder_population(N: int, af: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    M = len(af)
    pop = np.zeros((N, 2, M), dtype=np.uint8)
    for j, p in enumerate(af):
        n_alt = int(np.rint(2 * N * float(p)))
        n_alt = min(max(n_alt, 1), 2 * N - 1)
        idx = rng.choice(2 * N, size=n_alt, replace=False)
        # Build a contiguous 2N-copy vector explicitly.  Assigning through
        # pop[:,:,j].reshape(-1) is unsafe because that slice is strided and
        # reshape may return a copy rather than a writable view.
        copies = np.zeros(2 * N, dtype=np.uint8)
        copies[idx] = 1
        pop[:, :, j] = copies.reshape(N, 2)
    return pop


def make_gametes(pop: np.ndarray, parents: np.ndarray, rp: np.ndarray, rng: np.random.Generator, chunk: int) -> np.ndarray:
    C = len(parents)
    M = pop.shape[2]
    out = np.empty((C, M), dtype=np.uint8)
    arange_M = np.arange(M)
    for s in range(0, C, chunk):
        e = min(C, s + chunk)
        n = e - s
        state = np.empty((n, M), dtype=np.int8)
        state[:, 0] = rng.integers(0, 2, size=n, dtype=np.int8)
        if M > 1:
            switches = rng.random((n, M - 1)) < rp[None, :]
            # cumulative XOR from starting homolog
            state[:, 1:] = state[:, [0]] ^ (np.cumsum(switches, axis=1) % 2).astype(np.int8)
        pp = parents[s:e]
        # advanced indexing: rows n x M
        out[s:e] = pop[pp[:, None], state, arange_M[None, :]]
    return out


def disruption(pop: np.ndarray, a: np.ndarray, b: np.ndarray, parity: np.ndarray) -> np.ndarray:
    if len(a) == 0:
        return np.zeros(pop.shape[0], dtype=float)
    bad0 = (np.bitwise_xor(pop[:, 0, a], pop[:, 0, b]) != parity[None, :])
    bad1 = (np.bitwise_xor(pop[:, 1, a], pop[:, 1, b]) != parity[None, :])
    q = 0.5 * (bad0.astype(float) + bad1.astype(float))
    return q.mean(axis=1)


def next_generation(pop: np.ndarray, rp: np.ndarray, selfing: float, cand_mult: float, S: float,
                    a: np.ndarray, b: np.ndarray, parity: np.ndarray, rng: np.random.Generator, chunk: int) -> Tuple[np.ndarray, float]:
    N = pop.shape[0]
    C = max(N, int(math.ceil(N * cand_mult)))
    p1 = rng.integers(0, N, size=C)
    is_self = rng.random(C) < selfing
    p2 = p1.copy()
    out_idx = np.where(~is_self)[0]
    if len(out_idx):
        vals = rng.integers(0, N - 1, size=len(out_idx))
        p2vals = vals + (vals >= p1[out_idx])
        p2[out_idx] = p2vals
    g1 = make_gametes(pop, p1, rp, rng, chunk)
    g2 = make_gametes(pop, p2, rp, rng, chunk)
    cand = np.stack([g1, g2], axis=1)
    if S > 0 and len(a) > 0:
        D = disruption(cand, a, b, parity)
        fit = np.exp(-S * D)
        sf = fit.sum()
        probs = fit / sf if np.isfinite(sf) and sf > 0 else None
        chosen = rng.choice(C, size=N, replace=False, p=probs)
        mean_burden = float(np.mean(D[chosen]))
    else:
        chosen = rng.choice(C, size=N, replace=False)
        mean_burden = 0.0
    return cand[chosen].copy(), mean_burden


def pair_r2_from_pop(pop: np.ndarray, a: np.ndarray, b: np.ndarray, sample_idx: np.ndarray) -> np.ndarray:
    sub = pop[sample_idx]
    dosage = sub.sum(axis=1).astype(float)  # n x M
    out = np.full(len(a), np.nan, dtype=float)
    for i, (xj, yj) in enumerate(zip(a, b)):
        x, y = dosage[:, int(xj)], dosage[:, int(yj)]
        sx, sy = np.std(x), np.std(y)
        if sx <= 0 or sy <= 0:
            continue
        r = np.corrcoef(x, y)[0, 1]
        out[i] = r * r
    return out


def heterozygosity(pop: np.ndarray) -> float:
    return float(np.mean(pop[:, 0, :] != pop[:, 1, :]))


def simulate_rep(cfg: dict, chrom: str, prep: dict, rep: int):
    N = int(cfg["population_size"])
    sample_n = min(int(cfg["ld_sample_size"]), N)
    checkpoints = np.array(cfg["checkpoints"], dtype=int)
    G = int(cfg["generations"])
    if checkpoints[0] != 0 or checkpoints[-1] > G:
        raise ValueError("checkpoints must start at 0 and not exceed generations")
    positions = prep["marker_positions"]
    af = prep["af"]
    rp = recomb_probs(positions, float(cfg["recomb_rate_per_bp"]))
    base_seed = stable_int_seed(cfg["seed"], chrom, "rep", rep)
    frng = np.random.default_rng(base_seed)
    founder = founder_population(N, af, frng)
    g0_sample = np.random.default_rng(stable_int_seed(base_seed, "g0sample")).choice(N, size=sample_n, replace=False)

    K = len(prep["true_a"])
    T = len(checkpoints)
    r2_true = np.full((len(MODELS), T, K), np.nan, dtype=np.float32)
    het = np.full((len(MODELS), T), np.nan, dtype=np.float32)
    burden = np.full((len(MODELS), T), np.nan, dtype=np.float32)

    for mi, model in enumerate(MODELS):
        pop = founder.copy()
        rng = np.random.default_rng(stable_int_seed(base_seed, model))
        if model == "top1_aglink_selection":
            sa, sb, sp = prep["true_a"], prep["true_b"], prep["true_par"]
            S = float(cfg["selection_strength"])
        elif model == "random_aglink_selection":
            sa, sb, sp = prep["rand_a"], prep["rand_b"], prep["rand_par"]
            S = float(cfg["selection_strength"])
        else:
            sa = sb = sp = np.array([], dtype=np.int32)
            S = 0.0

        cp_to_idx = {int(g): i for i, g in enumerate(checkpoints)}
        ti = cp_to_idx[0]
        r2_true[mi, ti] = pair_r2_from_pop(pop, prep["true_a"], prep["true_b"], g0_sample)
        het[mi, ti] = heterozygosity(pop)
        burden[mi, ti] = float(np.mean(disruption(pop, sa, sb, sp))) if len(sa) else 0.0

        for gen in range(1, G + 1):
            pop, _ = next_generation(
                pop, rp, float(cfg["selfing_rate"]), float(cfg["candidate_multiplier"]), S,
                sa, sb, sp, rng, int(cfg.get("gamete_chunk_size", 512))
            )
            if gen in cp_to_idx:
                ti = cp_to_idx[gen]
                sidx = rng.choice(N, size=sample_n, replace=False)
                r2_true[mi, ti] = pair_r2_from_pop(pop, prep["true_a"], prep["true_b"], sidx)
                het[mi, ti] = heterozygosity(pop)
                burden[mi, ti] = float(np.mean(disruption(pop, sa, sb, sp))) if len(sa) else 0.0
    return r2_true, het, burden


def simulate(cfg: dict, chrom: str, rep_start: int, rep_end: int):
    prep_path = Path(cfg["output_root"]) / chrom / "01_prepare" / "prepared.npz"
    if not prep_path.exists():
        raise FileNotFoundError(f"Run prepare first: {prep_path}")
    prep_npz = np.load(prep_path)
    prep = {k: prep_npz[k] for k in prep_npz.files}
    R = int(cfg["replicates"])
    start = max(0, int(rep_start))
    end = min(R - 1, int(rep_end))
    if start > end:
        print(f"[{chrom}] no replicates in requested range {rep_start}-{rep_end}; configured R={R}")
        return
    outdir = mkdir(Path(cfg["output_root"]) / chrom / "02_sim")
    reps = []
    all_r2, all_het, all_burden = [], [], []
    for rep in range(start, end + 1):
        print(f"[{chrom}] replicate {rep+1}/{R}", flush=True)
        r2, h, b = simulate_rep(cfg, chrom, prep, rep)
        reps.append(rep)
        all_r2.append(r2)
        all_het.append(h)
        all_burden.append(b)
    path = outdir / f"shard_{start:03d}_{end:03d}.npz"
    np.savez_compressed(path,
                        reps=np.array(reps, dtype=np.int32),
                        r2_true=np.stack(all_r2),
                        heterozygosity=np.stack(all_het),
                        burden=np.stack(all_burden))
    print(f"[{chrom}] wrote {path}")


def merge(cfg: dict, chrom: str):
    base = Path(cfg["output_root"]) / chrom
    files = sorted((base / "02_sim").glob("shard_*.npz"))
    if not files:
        raise FileNotFoundError(f"No shards for {chrom}")
    seen = {}
    for f in files:
        z = np.load(f)
        for i, rep in enumerate(z["reps"].tolist()):
            if rep in seen:
                raise RuntimeError(f"Duplicate replicate {rep} in shards for {chrom}")
            seen[rep] = (z["r2_true"][i], z["heterozygosity"][i], z["burden"][i])
    R = int(cfg["replicates"])
    missing = [r for r in range(R) if r not in seen]
    if missing:
        raise RuntimeError(f"{chrom}: missing replicates: {missing}")
    r2 = np.stack([seen[r][0] for r in range(R)])
    het = np.stack([seen[r][1] for r in range(R)])
    burden = np.stack([seen[r][2] for r in range(R)])
    prep = np.load(base / "01_prepare" / "prepared.npz")
    real = prep["true_real_r2"].astype(float)
    checkpoints = np.array(cfg["checkpoints"], dtype=int)

    # Metrics per replicate/model/generation.
    dyn_pcc = np.full((R, len(MODELS), len(checkpoints)), np.nan)
    dyn_mean = np.full_like(dyn_pcc, np.nan)
    ncomplete = np.zeros_like(dyn_pcc, dtype=int)
    fixed_mean = np.full_like(dyn_pcc, np.nan)
    fixed_pcc = np.full_like(dyn_pcc, np.nan)
    fixed_real_mean = np.full(R, np.nan)
    fixed_n = np.zeros(R, dtype=int)
    for rr in range(R):
        # Same common true-link set for all models/times within a replicate.
        common = np.isfinite(real)
        common &= np.all(np.isfinite(r2[rr]), axis=(0, 1))
        fixed_n[rr] = int(common.sum())
        if common.any():
            fixed_real_mean[rr] = float(np.mean(real[common]))
        for mi in range(len(MODELS)):
            for ti in range(len(checkpoints)):
                x = r2[rr, mi, ti].astype(float)
                m = np.isfinite(x) & np.isfinite(real)
                ncomplete[rr, mi, ti] = int(m.sum())
                if m.any():
                    dyn_mean[rr, mi, ti] = float(np.mean(x[m]))
                dyn_pcc[rr, mi, ti] = pearson_safe(x[m], real[m]) if m.sum() >= 3 else np.nan
                if common.any():
                    fixed_mean[rr, mi, ti] = float(np.mean(x[common]))
                    fixed_pcc[rr, mi, ti] = pearson_safe(x[common], real[common]) if common.sum() >= 3 else np.nan

    out = mkdir(base / "03_merged")
    np.savez_compressed(out / "merged.npz", r2_true=r2, heterozygosity=het, burden=burden,
                        dyn_pcc=dyn_pcc, dyn_mean=dyn_mean, ncomplete=ncomplete,
                        fixed_mean=fixed_mean, fixed_pcc=fixed_pcc,
                        fixed_real_mean=fixed_real_mean, fixed_n=fixed_n)

    fields = ["chrom", "generation", "model",
              "dynamic_pcc_mean", "dynamic_pcc_ci_low", "dynamic_pcc_ci_high", "dynamic_pcc_n",
              "dynamic_mean_r2_mean", "dynamic_mean_r2_ci_low", "dynamic_mean_r2_ci_high",
              "fixed_mean_r2_mean", "fixed_mean_r2_ci_low", "fixed_mean_r2_ci_high",
              "fixed_pcc_mean", "fixed_pcc_ci_low", "fixed_pcc_ci_high",
              "real_fixed_mean_r2", "fixed_link_count_median",
              "heterozygosity_mean", "target_burden_mean"]
    with open(out / "metrics_summary.tsv", "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=fields)
        w.writeheader()
        for ti, gen in enumerate(checkpoints):
            for mi, model in enumerate(MODELS):
                dpm, dplo, dphi, dpn = mean_ci(dyn_pcc[:, mi, ti])
                dmm, dmlo, dmhi, _ = mean_ci(dyn_mean[:, mi, ti])
                fmm, fmlo, fmhi, _ = mean_ci(fixed_mean[:, mi, ti])
                fpm, fplo, fphi, _ = mean_ci(fixed_pcc[:, mi, ti])
                hm, _, _, _ = mean_ci(het[:, mi, ti])
                bm, _, _, _ = mean_ci(burden[:, mi, ti])
                w.writerow({
                    "chrom": chrom, "generation": int(gen), "model": model,
                    "dynamic_pcc_mean": dpm, "dynamic_pcc_ci_low": dplo, "dynamic_pcc_ci_high": dphi, "dynamic_pcc_n": dpn,
                    "dynamic_mean_r2_mean": dmm, "dynamic_mean_r2_ci_low": dmlo, "dynamic_mean_r2_ci_high": dmhi,
                    "fixed_mean_r2_mean": fmm, "fixed_mean_r2_ci_low": fmlo, "fixed_mean_r2_ci_high": fmhi,
                    "fixed_pcc_mean": fpm, "fixed_pcc_ci_low": fplo, "fixed_pcc_ci_high": fphi,
                    "real_fixed_mean_r2": float(np.nanmean(fixed_real_mean)),
                    "fixed_link_count_median": float(np.median(fixed_n)),
                    "heterozygosity_mean": hm, "target_burden_mean": bm,
                })
    print(f"[{chrom}] merged {R} replicates")


def read_tsv(path: Path) -> List[dict]:
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def audit(cfg: dict, chrom: str):
    base = Path(cfg["output_root"]) / chrom
    qc = load_json(base / "01_prepare" / "prepare_qc.json")
    z = np.load(base / "03_merged" / "merged.npz")
    problems = []
    eligible_by_block = qc.get("eligible_count_by_acr_block", {})
    selected_by_block = qc.get("selected_count_by_acr_block", {})
    expected_by_block = {b: selected_count_for_block(cfg, int(n)) for b, n in eligible_by_block.items()}
    expected_top = int(sum(expected_by_block.values()))
    if qc["n_selected_top1"] != expected_top:
        problems.append(f"selected count {qc['n_selected_top1']} != per-block expected {expected_top}")
    if selected_by_block != expected_by_block:
        bad = [b for b in expected_by_block if int(selected_by_block.get(b, -1)) != int(expected_by_block[b])]
        problems.append(f"per-block selected counts mismatch in {len(bad)} blocks: {bad[:10]}")
    if qc["n_random"] != qc["n_selected_top1"]:
        problems.append("random count != top1 count")
    for k in ["promoter_overlap_selected", "promoter_overlap_random", "random_known_sig_overlap"]:
        if qc.get(k, 0) != 0:
            problems.append(f"{k}={qc[k]}")
    if cfg.get("random_disallow_shared_acr", True) and qc.get("random_shared_acr_with_top", 0):
        problems.append("random shares ACR with top1")
    if cfg.get("random_disallow_shared_gene", True) and qc.get("random_shared_gene_with_top", 0):
        problems.append("random shares gene with top1")
    if cfg.get("random_disallow_shared_marker", True) and qc.get("random_shared_marker_with_top", 0):
        problems.append("random shares representative marker with top1")
    r2 = z["r2_true"]
    # G0 must be exactly equal because all models start from identical founder and sample.
    max_g0 = float(np.nanmax(np.abs(r2[:, :, 0, :] - r2[:, [0], 0, :])))
    if max_g0 > 0:
        problems.append(f"G0 identity mismatch max={max_g0}")
    if r2.shape[0] != int(cfg["replicates"]):
        problems.append(f"replicate count {r2.shape[0]} != configured {cfg['replicates']}")
    if r2.shape[1] != 3:
        problems.append("not exactly 3 models")
    report = {
        "version": "v17_raw_standalone",
        "chrom": chrom,
        "status": "PASS" if not problems else "FAIL",
        "problems": problems,
        "n_selected_top1": qc["n_selected_top1"],
        "n_random": qc["n_random"],
        "n_unique_markers": qc["n_unique_markers"],
        "max_g0_identity_error": max_g0,
        "fixed_link_count_min": int(np.min(z["fixed_n"])),
        "fixed_link_count_median": float(np.median(z["fixed_n"])),
        "fixed_link_count_max": int(np.max(z["fixed_n"])),
        "prepare_qc": qc,
    }
    dump_json(report, base / "04_audit" / "audit_report.json")
    print(f"[{chrom}] audit {report['status']}")
    if problems:
        for p in problems:
            print("  -", p)
        raise RuntimeError(f"{chrom} audit failed")


def check(cfg: dict):
    raw = Path(cfg["raw_dir"])
    problems = []
    for key in ["acrs", "aglinks", "blocks", "gff", "vcf_dir"]:
        p = raw / cfg["files"][key]
        if not p.exists():
            problems.append(str(p))
    if problems:
        raise FileNotFoundError("Missing raw input(s):\n" + "\n".join(problems))
    for c in cfg["chromosomes"]:
        print(f"{c}\t{find_vcf(raw / cfg['files']['vcf_dir'], c)}")
    print("Raw input check PASS")


def load_config(path: Path, demo: bool = False) -> dict:
    cfg = load_json(path)
    if demo:
        cfg.update(raw_dir=str(path.parent / 'demo'), output_root='output_demo',
                   chromosomes=['1A'], population_size=203, ld_sample_size=203,
                   replicates=3, generations=2, checkpoints=[0,1,2])
    elif not cfg.get('raw_dir') and not os.environ.get('RAW_DIR'):
        raise ValueError('Set RAW_DIR to the complete raw input directory, or run bash 00_RUN.sh demo.')
    for key in ("raw_dir", "output_root"):
        value = Path(cfg[key]).expanduser()
        cfg[key] = str((path.parent / value).resolve())
    if os.environ.get("RAW_DIR") and not demo:
        cfg["raw_dir"] = str(Path(os.environ["RAW_DIR"]).expanduser().resolve())
    if os.environ.get("OUTPUT_DIR"):
        cfg["output_root"] = str(Path(os.environ["OUTPUT_DIR"]).expanduser().resolve() / Path(cfg["output_root"]).name)
    if cfg.get("selection_scope") != "acr_block":
        raise ValueError("selection_scope must be acr_block")
    if cfg.get("top_pcc_mode", "raw") not in ("raw", "abs"):
        raise ValueError("top_pcc_mode must be raw or abs")
    if not cfg["chromosomes"] or len(set(cfg["chromosomes"])) != len(cfg["chromosomes"]):
        raise ValueError("chromosomes must be non-empty and unique")
    for key, minimum in (("population_size", 3), ("ld_sample_size", 3),
                         ("replicates", 1), ("generations", 0), ("gamete_chunk_size", 1)):
        value = cfg[key]
        if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
            raise ValueError(f"{key} must be an integer >= {minimum}")
    cps = cfg["checkpoints"]
    if (not cps or any(not isinstance(g, int) or isinstance(g, bool) for g in cps)
            or cps != sorted(set(cps)) or cps[0] != 0 or cps[-1] > cfg["generations"]):
        raise ValueError("checkpoints must be increasing unique integers from 0 through generations")
    if not 0 <= float(cfg["selfing_rate"]) <= 1:
        raise ValueError("selfing_rate must be between 0 and 1")
    multiplier = float(cfg["candidate_multiplier"])
    if not math.isfinite(multiplier) or multiplier < 1:
        raise ValueError("candidate_multiplier must be finite and >= 1")
    for key in ("selection_strength", "recomb_rate_per_bp", "promoter_upstream_bp",
                "representative_snp_window_bp", "random_distance_min_abs_bp"):
        if not math.isfinite(float(cfg[key])) or float(cfg[key]) < 0:
            raise ValueError(f"{key} must be finite and non-negative")
    tolerances = cfg["random_distance_relative_tolerances"]
    if not tolerances or any(not math.isfinite(float(t)) or float(t) < 0 for t in tolerances):
        raise ValueError("random distance tolerances must be finite and non-negative")
    if cfg.get("per_block_top_count") is not None or cfg.get("per_block_top_fraction") != 0.01:
        raise ValueError("Only per-block top 1% selection is supported (per_block_top_fraction=0.01)")
    return cfg


def full(cfg: dict):
    """Run the original analysis stages in a fresh output directory."""
    root = Path(cfg["output_root"])
    if root.exists() and (not root.is_dir() or any(root.iterdir())):
        raise FileExistsError(f"Output is not empty: {root}. Set OUTPUT_DIR to a new directory.")
    check(cfg)
    dump_json(cfg, root / "config.json")
    prepare_all(cfg)
    for chrom in cfg["chromosomes"]:
        simulate(cfg, chrom, 0, int(cfg["replicates"]) - 1)
        merge(cfg, chrom)
        audit(cfg, chrom)


def main():
    p = argparse.ArgumentParser(description="Forward agLink simulation with per-ACR-block PCC top 1% selection")
    p.add_argument("--config", type=Path, required=True)
    p.add_argument("--demo", action="store_true", help="Use the bundled 1A subset with 3 replicates, N=203, G0-G2")
    p.add_argument("--output", type=Path, help="Override the result directory")
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("check")
    pa = sub.add_parser("prepare-all"); pa.add_argument("--force", action="store_true")
    pc = sub.add_parser("prepare"); pc.add_argument("--chrom", required=True); pc.add_argument("--force", action="store_true")
    ps = sub.add_parser("simulate"); ps.add_argument("--chrom", required=True); ps.add_argument("--rep-start", type=int, default=0); ps.add_argument("--rep-end", type=int, default=999999)
    pm = sub.add_parser("merge"); pm.add_argument("--chrom", required=True)
    pu = sub.add_parser("audit"); pu.add_argument("--chrom", required=True)
    sub.add_parser("full")
    args = p.parse_args()
    try:
        cfg = load_config(args.config.resolve(), demo=args.demo)
    except (OSError, ValueError) as error:
        p.error(str(error))
    if args.output is not None:
        cfg['output_root'] = str(args.output.expanduser().resolve())
    if hasattr(args, "chrom") and args.chrom not in cfg["chromosomes"]:
        p.error(f"chromosome {args.chrom} is not in config chromosomes")
    if args.cmd == "check": check(cfg)
    elif args.cmd == "prepare-all": prepare_all(cfg, force=args.force)
    elif args.cmd == "prepare": prepare_chrom(cfg, args.chrom, force=args.force)
    elif args.cmd == "simulate": simulate(cfg, args.chrom, args.rep_start, args.rep_end)
    elif args.cmd == "merge": merge(cfg, args.chrom)
    elif args.cmd == "audit": audit(cfg, args.chrom)
    elif args.cmd == "full": full(cfg)


if __name__ == "__main__":
    try:
        main()
    except (FileNotFoundError, FileExistsError) as error:
        raise SystemExit(f'Error: {error}\nCheck RAW_DIR and the output directory; use bash 00_RUN.sh demo for bundled inputs.')
