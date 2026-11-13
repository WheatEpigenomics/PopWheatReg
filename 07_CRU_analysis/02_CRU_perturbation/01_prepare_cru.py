#!/usr/bin/env python3
"""CRU perturbation: prepare sentinel-low/high effects from raw ATAC/RNA matrices.

Consolidates analysis_helpers.py, prepare_data.py and prepare_high_groups.py
from csv2_T8_CRU_20261109. Scientific rules and random seeds are preserved.
Requires Python >=3.8 and NumPy. See README.md and --help for input formats.
"""
import argparse
import bisect
import csv
import gzip
import hashlib
import os
import platform
import random
import shutil
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np


# ========================================================================
# Coordinates, sentinel membership and distance matching
# ========================================================================

def log(msg):
    sys.stderr.write(f"[CRU-natural-perturbation] {msg}\n")
    sys.stderr.flush()

def norm_chrom(chrom):
    chrom = chrom.strip()
    if chrom.startswith("Chr"):
        chrom = chrom[3:]
    return chrom

def parse_attrs(attr):
    out = {}
    for item in attr.split(";"):
        if item and "=" in item:
            key, val = item.split("=", 1)
            out[key] = val
    return out

def read_matrix_header(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
    if len(header) < 2:
        raise ValueError(f"Matrix has fewer than two columns: {path}")
    return header

def read_id_set(path):
    ids = set()
    with open(path) as fh:
        next(fh, None)
        for line in fh:
            if line.strip():
                ids.add(line.split("\t", 1)[0])
    return ids

def read_crus(path):
    crus = {}
    by_chrom = defaultdict(list)
    with open(path) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                raise ValueError(f"Malformed CRU BED line {line_no}")
            chrom = norm_chrom(parts[0])
            start = int(parts[1])
            end = int(parts[2])
            cru_id = parts[3]
            if end <= start:
                raise ValueError(f"Invalid CRU coordinates at line {line_no}")
            row = {
                "cru_id": cru_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "length_bp": end - start,
            }
            if cru_id in crus:
                raise ValueError("Duplicate CRU identifier: " + cru_id)
            crus[cru_id] = row
            by_chrom[chrom].append((start, end, cru_id))
    index = {}
    for chrom, rows in by_chrom.items():
        rows.sort()
        index[chrom] = {"starts": [r[0] for r in rows], "rows": rows}
    if not crus:
        raise ValueError("CRU BED contains no intervals")
    for chrom, entry in index.items():
        rows = entry["rows"]
        if any(left[1] > right[0] for left, right in zip(rows, rows[1:])):
            raise ValueError("Overlapping CRU intervals on chromosome " + chrom)
    return crus, index

def find_cru(index, chrom, pos0):
    chrom = norm_chrom(chrom)
    if chrom not in index:
        return None
    starts = index[chrom]["starts"]
    rows = index[chrom]["rows"]
    j = bisect.bisect_right(starts, pos0) - 1
    if j >= 0 and pos0 < rows[j][1]:
        return rows[j][2]
    return None

def read_acr_coordinates(acr_bed, cru_index):
    acr_to_cru = {}
    acr_info = {}
    acrs_by_cru = defaultdict(list)
    acrs_by_chrom = defaultdict(list)
    total = 0
    in_cru = 0
    with open(acr_bed) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                raise ValueError(f"Malformed ACR BED line {line_no}")
            chrom = norm_chrom(parts[0])
            start = int(parts[1])
            end = int(parts[2])
            acr_id = parts[3]
            midpoint = (start + end) // 2
            cru_id = find_cru(cru_index, chrom, midpoint)
            acrs_by_chrom[chrom].append((midpoint, acr_id, cru_id))
            total += 1
            if cru_id is not None:
                acr_to_cru[acr_id] = cru_id
                acr_info[acr_id] = {
                    "acr_id": acr_id,
                    "chrom": chrom,
                    "start": start,
                    "end": end,
                    "midpoint": midpoint,
                    "cru_id": cru_id,
                }
                acrs_by_cru[cru_id].append(acr_id)
                in_cru += 1
    for chrom, rows in list(acrs_by_chrom.items()):
        rows.sort()
        acrs_by_chrom[chrom] = {"mids": [r[0] for r in rows], "rows": rows}
    return acr_to_cru, acr_info, acrs_by_cru, acrs_by_chrom, total, in_cru

def read_gene_coordinates(gff3, cru_index, expressed_gene_ids):
    gene_info = {}
    genes_by_cru = defaultdict(list)
    genes_by_chrom = defaultdict(list)
    total_gene_features = 0
    expressed_with_coords = 0
    expressed_in_cru = 0
    with open(gff3) as fh:
        for line_no, line in enumerate(fh, 1):
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[2] != "gene":
                continue
            total_gene_features += 1
            attrs = parse_attrs(parts[8])
            gene_id = attrs.get("ID") or attrs.get("Name")
            if not gene_id or gene_id not in expressed_gene_ids:
                continue
            chrom = norm_chrom(parts[0])
            start = int(parts[3]) - 1
            end = int(parts[4])
            strand = parts[6]
            if strand == "+":
                tss = start
            elif strand == "-":
                tss = end - 1
            else:
                tss = (start + end) // 2
            cru_id = find_cru(cru_index, chrom, tss)
            gene_info[gene_id] = {
                "gene_id": gene_id,
                "chrom": chrom,
                "start": start,
                "end": end,
                "tss": tss,
                "strand": strand,
                "cru_id": cru_id,
            }
            genes_by_chrom[chrom].append((tss, gene_id, cru_id))
            expressed_with_coords += 1
            if cru_id is not None:
                genes_by_cru[cru_id].append(gene_id)
                expressed_in_cru += 1
    for chrom, rows in list(genes_by_chrom.items()):
        rows.sort()
        genes_by_chrom[chrom] = {"mids": [r[0] for r in rows], "rows": rows}
    return (
        gene_info,
        genes_by_cru,
        genes_by_chrom,
        total_gene_features,
        expressed_with_coords,
        expressed_in_cru,
    )

def nearest_outside_record(index, target_pos, exclude_cru, tolerance_bp, exclude_id=None):
    if not index:
        return None
    mids = index["mids"]
    rows = index["rows"]
    j = bisect.bisect_left(mids, target_pos)
    best = None
    for radius in range(0, 200):
        candidates = []
        left = j - radius - 1
        right = j + radius
        if left >= 0:
            candidates.append(rows[left])
        if right < len(rows):
            candidates.append(rows[right])
        if not candidates and left < 0 and right >= len(rows):
            break
        for midpoint, rec_id, cru_id in candidates:
            if rec_id == exclude_id:
                continue
            if cru_id == exclude_cru:
                continue
            dist = abs(midpoint - target_pos)
            if dist <= tolerance_bp:
                if best is None or dist < best[0]:
                    best = (dist, midpoint, rec_id, cru_id)
        if best is not None:
            return best
    return None

def make_distance_matched_pairs(
    sentinels,
    eligible_crus,
    acrs_by_cru,
    acr_info,
    acrs_by_chrom,
    genes_by_cru,
    gene_info,
    genes_by_chrom,
    rng,
    max_acr_controls_per_cru,
    max_gene_controls_per_cru,
    min_tolerance_bp,
    tolerance_fraction,
):
    acr_pairs = []
    gene_pairs = []
    for cru_id in sorted(eligible_crus):
        sentinel = sentinels[cru_id]
        chrom = sentinel["sentinel_chrom"]
        sentinel_mid = sentinel["sentinel_midpoint"]
        other_acrs = [a for a in acrs_by_cru[cru_id] if a != sentinel["sentinel_acr"]]
        if len(other_acrs) > max_acr_controls_per_cru:
            other_acrs = rng.sample(other_acrs, max_acr_controls_per_cru)
        for acr_id in other_acrs:
            d = abs(acr_info[acr_id]["midpoint"] - sentinel_mid)
            tolerance = max(min_tolerance_bp, int(d * tolerance_fraction))
            side = -1 if rng.random() < 0.5 else 1
            targets = [sentinel_mid + side * d, sentinel_mid - side * d]
            match = None
            for target in targets:
                match = nearest_outside_record(
                    acrs_by_chrom.get(chrom, []),
                    target,
                    cru_id,
                    tolerance,
                    exclude_id=sentinel["sentinel_acr"],
                )
                if match is not None:
                    break
            if match is None:
                continue
            acr_pairs.append(
                {
                    "cru_id": cru_id,
                    "sentinel_acr": sentinel["sentinel_acr"],
                    "observed_acr": acr_id,
                    "control_acr": match[2],
                    "observed_distance_bp": d,
                    "control_distance_error_bp": match[0],
                    "control_cru_id": match[3] or "outside_any_cru",
                }
            )

        target_genes = genes_by_cru.get(cru_id, [])
        if len(target_genes) > max_gene_controls_per_cru:
            target_genes = rng.sample(target_genes, max_gene_controls_per_cru)
        for gene_id in target_genes:
            d = abs(gene_info[gene_id]["tss"] - sentinel_mid)
            tolerance = max(min_tolerance_bp, int(d * tolerance_fraction))
            side = -1 if rng.random() < 0.5 else 1
            targets = [sentinel_mid + side * d, sentinel_mid - side * d]
            match = None
            for target in targets:
                match = nearest_outside_record(
                    genes_by_chrom.get(chrom, []),
                    target,
                    cru_id,
                    tolerance,
                    exclude_id=gene_id,
                )
                if match is not None:
                    break
            if match is None:
                continue
            gene_pairs.append(
                {
                    "cru_id": cru_id,
                    "sentinel_acr": sentinel["sentinel_acr"],
                    "observed_gene": gene_id,
                    "control_gene": match[2],
                    "observed_distance_bp": d,
                    "control_distance_error_bp": match[0],
                    "control_cru_id": match[3] or "outside_any_cru",
                }
            )
    return acr_pairs, gene_pairs


# ========================================================================
# Sentinel selection and low-versus-upper-half effects
# ========================================================================

def write_table(path, fields, rows):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def read_selected_matrix(path, wanted, columns):
    values = np.empty((len(wanted), len(columns)), dtype=np.float64)
    names = []
    with open(path) as handle:
        next(handle)
        for line_no, line in enumerate(handle, 1):
            feature, _, raw = line.rstrip("\n").partition("\t")
            if feature not in wanted:
                continue
            row = np.fromstring(raw.replace("NA", "nan"), sep="\t")
            if len(row) <= max(columns):
                raise ValueError("Malformed matrix row: " + feature)
            row = row[columns]
            if np.any(row[np.isfinite(row)] < 0) or np.any(np.isinf(row)):
                raise ValueError("Negative or infinite measurement: " + feature)
            values[len(names)] = row
            names.append(feature)
            if len(names) % 100000 == 0:
                log("Loaded matrix rows: " + str(len(names)))
    if len(names) != len(set(names)):
        raise ValueError("Duplicate feature identifiers in " + str(path))
    return names, values[:len(names)], {name: i for i, name in enumerate(names)}


def select_sentinels(names, values, to_cru, info, min_group=15):
    by_cru = defaultdict(list)
    for i, name in enumerate(names):
        by_cru[to_cru[name]].append(i)
    sentinels = {}
    for cru, indices in by_cru.items():
        block = values[indices]
        count = np.isfinite(block).sum(axis=1)
        ordered = np.sort(block, axis=1)
        quantiles = []
        for q in (.20, .50, .80):
            location = q * np.maximum(count - 1, 0)
            lo = np.floor(location).astype(int)
            hi = np.ceil(location).astype(int)
            fraction = location - lo
            rows = np.arange(len(indices))
            quantiles.append(ordered[rows, lo] * (1-fraction) + ordered[rows, hi] * fraction)
        low, normal, high = quantiles
        zeros = (block <= 0).sum(axis=1) / np.maximum(count, 1)
        nlow = (block <= low[:, None]).sum(axis=1)
        nnormal = (block >= normal[:, None]).sum(axis=1)
        valid = ((count >= 2*min_group) & (zeros <= .80) & (low < normal)
                 & (nlow >= min_group) & (nnormal >= min_group))
        if not np.any(valid):
            continue
        # argmax returns the first matrix row in a tie, as in the historical scan.
        j = int(np.argmax(np.where(valid, high-low, -np.inf)))
        row_index = indices[j]
        name = names[row_index]
        raw = values[row_index]
        ilow = np.flatnonzero(raw <= low[j]).tolist()
        inormal = np.flatnonzero(raw >= normal[j]).tolist()
        delta = float(np.mean(np.log2(raw[ilow]+1))-np.mean(np.log2(raw[inormal]+1)))
        if not delta < 0:
            raise ValueError("Nonnegative sentinel contrast: " + cru)
        sentinels[cru] = dict(cru_id=cru, sentinel_acr=name,
            sentinel_chrom=info[name]["chrom"], sentinel_start=info[name]["start"],
            sentinel_end=info[name]["end"], sentinel_midpoint=info[name]["midpoint"],
            low_threshold=float(low[j]), normal_threshold=float(normal[j]),
            zero_fraction=float(zeros[j]), dynamic_range_q80_minus_q20=float(high[j]-low[j]),
            anchor_delta_log2_low_vs_normal=delta, low_idx=ilow, normal_idx=inormal,
            n_low=len(ilow), n_normal=len(inormal), matrix_row=row_index)
    return sentinels


def effects(matrix, indices, sentinel):
    block = matrix[indices]
    low = np.log2(block[:, sentinel["low_idx"]] + 1.0)
    normal = np.log2(block[:, sentinel["normal_idx"]] + 1.0)
    nl = np.isfinite(low).sum(axis=1)
    nn = np.isfinite(normal).sum(axis=1)
    ml = np.nansum(low, axis=1)/np.maximum(nl, 1)
    mn = np.nansum(normal, axis=1)/np.maximum(nn, 1)
    return ml, mn, ml-mn, nl, nn


def prepare_low(args):
    inp, out = args.input_root.resolve(), args.outdir.resolve()
    paths = {
        "cru": args.cru_bed or inp/"14_CRU/csv2_T8_CRU.bed",
        "acr_bed": args.acr_bed or inp/"04_acr/05_acr_with_name_fil110.bed",
        "acr_signal": args.acr_matrix or inp/"04_acr/07_200_spike_ACR_fil110_signal.txt",
        "expression": args.expression_matrix or inp/"03_exp/05_202hex_csV2_rna_spike_ave_exp_sum1_clean.txt",
        "gff": args.gff or inp/"01_ref/iwgsc_refseqv2.1_gene_annotation_200916/iwgsc_refseqv2.1_annotation_200916_HC.gff3",
    }
    paths = {key: value.resolve() for key, value in paths.items()}
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    ah = read_matrix_header(paths["acr_signal"])[1:]
    gh = read_matrix_header(paths["expression"])[1:]
    if len(set(ah)) != len(ah) or len(set(gh)) != len(gh):
        raise ValueError("Duplicate accession identifiers")
    samples = [name for name in ah if name in set(gh)]
    ai = [ah.index(s) for s in samples]
    gi = [gh.index(s) for s in samples]
    if len(samples) < 30:
        raise ValueError("At least 30 shared accessions are needed; normally >=75 to obtain a 15-sample 20% group")
    log("Shared accessions: " + str(len(samples)))
    write_table(out/"data/accessions.tsv", ["accession", "ATAC_column_1based", "RNA_column_1based"],
        [dict(accession=s, ATAC_column_1based=ah.index(s)+2, RNA_column_1based=gh.index(s)+2) for s in samples])
    crus, cru_index = read_crus(paths["cru"])
    a2c, ainfo, abycru, aby_chr, total_acr, inside_acr = read_acr_coordinates(paths["acr_bed"], cru_index)
    expressed = read_id_set(paths["expression"])
    ginfo, gbycru, gby_chr, ngene, ngcoords, nginside = read_gene_coordinates(paths["gff"], cru_index, expressed)
    log("Loading inside-CRU accessibility: " + str(inside_acr) + " rows")
    anames, amatrix, aindex = read_selected_matrix(paths["acr_signal"], a2c, ai)
    if set(anames) != set(a2c):
        raise ValueError("ACR matrix/coordinate identifiers disagree")
    log("Selecting the largest Q80-Q20 sentinel per CRU")
    sentinels = select_sentinels(anames, amatrix, a2c, ainfo)
    eligible = {c for c in sentinels if len(abycru[c])-1 >= 3 and len(gbycru[c]) >= 1}
    log("Sentinels: %s; eligible CRUs: %s" % (len(sentinels), len(eligible)))
    if not eligible:
        raise ValueError("No eligible CRUs: need a valid sentinel, >=3 other ACRs and >=1 expressed gene")
    rng = random.Random(20260812)
    apairs, gpairs = make_distance_matched_pairs(sentinels, eligible, abycru, ainfo, aby_chr,
        gbycru, ginfo, gby_chr, rng, 25, 50, 10000, .20)
    ac = {p["cru_id"] for p in apairs}
    gc = {p["cru_id"] for p in gpairs}
    log("Pairs: ACR=%s, gene=%s; CRUs with BOTH control types=%s/%s" %
        (len(apairs), len(gpairs), len(ac & gc), len(eligible)))
    if not ac.intersection(gc):
        raise ValueError("No CRU has both ACR and gene distance-matched controls")
    sentinel_rows = []
    for c in sorted(crus):
        row = dict(crus[c])
        row["n_acrs_in_cru"] = len(abycru[c])
        row["n_genes_in_cru"] = len(gbycru[c])
        row["eligible"] = int(c in eligible)
        row["has_sentinel"] = int(c in sentinels)
        row["exclusion_reason"] = "NA" if c in eligible else (
            "no_valid_sentinel" if c not in sentinels else "insufficient_inside_ACRs_or_genes")
        if c in sentinels:
            row.update({k: v for k, v in sentinels[c].items() if k not in
                        ("matrix_row", "low_idx", "normal_idx", "cru_id")})
        sentinel_rows.append(row)
    sf = list(sentinel_rows[0])
    for row in sentinel_rows:
        for key in row:
            if key not in sf:
                sf.append(key)
    for row in sentinel_rows:
        for key in sf:
            row.setdefault(key, "NA")
    write_table(out/"data/sentinel_selection.tsv", sf, sentinel_rows)
    write_table(out/"data/sentinel_accessibility.tsv.gz", ["cru_id"]+samples,
        [dict(cru_id=c, **{s: format(amatrix[sentinels[c]["matrix_row"], i], ".17g")
                            for i, s in enumerate(samples)}) for c in sorted(sentinels)])

    # Add coordinates and true distance errors, without changing pair-selection RNG.
    pairs = []
    for kind, rows, positions in (("ACR", apairs, aby_chr), ("gene", gpairs, gby_chr)):
        needed = {p["control_acr" if kind == "ACR" else "control_gene"] for p in rows}
        pos = {feature: coordinate for index in positions.values() for coordinate, feature, _ in index["rows"] if feature in needed}
        for p in rows:
            c = p["cru_id"]
            observed = p["observed_acr" if kind == "ACR" else "observed_gene"]
            control = p["control_acr" if kind == "ACR" else "control_gene"]
            distance = abs(pos[control]-sentinels[c]["sentinel_midpoint"])
            pairs.append(dict(cru_id=c, feature_type=kind, observed_feature=observed, control_feature=control,
                control_cru_id=p["control_cru_id"], observed_distance_bp=p["observed_distance_bp"],
                control_distance_bp=distance, control_distance_error_bp=p["control_distance_error_bp"],
                absolute_distance_error_bp=abs(distance-p["observed_distance_bp"]),
                tolerance_bp=max(10000, int(p["observed_distance_bp"]*.20))))
    write_table(out/"data/distance_matched_pairs.tsv.gz", list(pairs[0]), pairs)
    outside_needed = {p["control_acr"] for p in apairs if p["control_acr"] not in aindex}
    log("Loading matched ACR controls outside original inside-CRU matrix: " + str(len(outside_needed)))
    cnames, cmatrix, cindex = read_selected_matrix(paths["acr_signal"], outside_needed, ai)
    if set(cnames) != outside_needed:
        raise ValueError("Missing ACR controls in the signal matrix")
    gnames, gmatrix, gindex = read_selected_matrix(paths["expression"], ginfo, gi)
    if set(gnames) != set(ginfo):
        raise ValueError("Missing gene expression rows")
    controls = defaultdict(set)
    for p in pairs:
        controls[(p["cru_id"],p["feature_type"])].add(p["control_feature"])
    effect_path = out/"data/feature_effects.tsv.gz"
    effect_fields = ["cru_id","feature_type","role","feature_id","n_low","n_normal",
                     "mean_low_log2","mean_normal_log2","delta_log2"]
    counts = defaultdict(int)
    with gzip.open(effect_path, "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(effect_fields)
        for no, c in enumerate(sorted(eligible), 1):
            s = sentinels[c]
            for kind, role, ids in (
                ("ACR","inside",[a for a in abycru[c] if a != s["sentinel_acr"]]),
                ("gene","inside",gbycru[c]),
                ("ACR","control",sorted(controls[(c,"ACR")])),
                ("gene","control",sorted(controls[(c,"gene")]))):
                batches = []
                if kind == "gene":
                    batches.append((ids, gmatrix, gindex))
                else:
                    batches.append(([a for a in ids if a in aindex], amatrix, aindex))
                    batches.append(([a for a in ids if a in cindex], cmatrix, cindex))
                for names, matrix, index in batches:
                    if not names:
                        continue
                    ml, mn, delta, nl, nn = effects(matrix, [index[a] for a in names], s)
                    for j, name in enumerate(names):
                        if nl[j] < 15 or nn[j] < 15:
                            continue
                        writer.writerow([c,kind,role,name,int(nl[j]),int(nn[j]),
                            format(ml[j],".17g"),format(mn[j],".17g"),format(delta[j],".17g")])
                        counts[(kind,role)] += 1
            if no % 1000 == 0:
                log("CRU effect summaries written: " + str(no))
    log("Hashing read-only raw inputs for provenance")
    manifests = []
    for label, path in paths.items():
        digest = hashlib.sha256()
        with open(path,"rb") as handle:
            for chunk in iter(lambda: handle.read(4*1024*1024), b""):
                digest.update(chunk)
        manifests.append(dict(input=label, relative_path=os.path.relpath(path, inp),
            bytes=path.stat().st_size, sha256=digest.hexdigest()))
    write_table(out/"docs/raw_input_manifest.tsv", list(manifests[0]), manifests)
    metadata = dict(n_ATAC_accessions=len(ah), n_RNA_accessions=len(gh), n_shared_accessions=len(samples),
        n_input_CRUs=len(crus), n_sentinels=len(sentinels), n_eligible_CRUs=len(eligible),
        n_ACR_pairs=len(apairs), n_gene_pairs=len(gpairs), n_CRUs_both_control_types=len(ac & gc),
        low_quantile=.20, normal_quantile=.50, min_group_size=15, min_other_ACRs=3, min_genes=1,
        max_zero_fraction=.80, max_ACR_pairs_per_CRU=25, max_gene_pairs_per_CRU=50,
        minimum_tolerance_bp=10000, tolerance_fraction=.20, control_search_radius_records=200,
        matching_seed=20260812, permutation_seed=20261021, pooled_permutation_seed=20261022, permutations=10000,
        python_version=platform.python_version(), numpy_version=np.__version__)
    write_table(out/"docs/analysis_parameters.tsv", ["parameter","value"],
        [dict(parameter=k,value=v) for k,v in metadata.items()])
    bed_copy = out/"data/cru_input.bed"
    if paths["cru"] != bed_copy.resolve():
        shutil.copyfile(paths["cru"], bed_copy)
    log("Finished raw-data recomputation; feature counts: " + str(dict(counts)))


# ========================================================================
# High-versus-lower-half effects using the same sentinels and pairs
# ========================================================================

def read(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt", newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def prepare_high(args):
    root = args.outdir.resolve()
    # Verify that the raw matrices are exactly those used for the low contrast.
    inputs = {}
    for row in read(root/"docs/raw_input_manifest.tsv"):
        if row["input"] not in ("acr_signal", "expression"):
            continue
        path = args.input_root/row["relative_path"]
        digest = hashlib.sha256()
        with open(path, "rb") as handle:
            for chunk in iter(lambda: handle.read(4*1024*1024), b""):
                digest.update(chunk)
        if path.stat().st_size != int(row["bytes"]) or digest.hexdigest() != row["sha256"]:
            raise ValueError("Raw input differs from the frozen low analysis: " + str(path))
        inputs[row["input"]] = path
    samples = [r["accession"] for r in read(root/"data/accessions.tsv")]
    selected = {r["cru_id"]: r for r in read(root/"data/sentinel_selection.tsv")}
    groups = {}; group_rows = []
    for row in read(root/"data/sentinel_accessibility.tsv.gz"):
        cru = row["cru_id"]
        raw = np.array([float(row[s]) for s in samples])
        ordered = np.sort(raw[np.isfinite(raw)])
        def quantile(q):
            position = q*(len(ordered)-1)
            lo, hi = int(np.floor(position)), int(np.ceil(position))
            f = position-lo
            return float(ordered[lo]*(1-f)+ordered[hi]*f)
        high, normal = quantile(.8), quantile(.5)
        hi = np.flatnonzero(raw >= high)
        ni = np.flatnonzero(raw <= normal)
        valid = high > normal and len(hi)>=15 and len(ni)>=15
        delta = float(np.mean(np.log2(raw[hi]+1))-np.mean(np.log2(raw[ni]+1))) if valid else float("nan")
        if valid and not delta>0:
            raise ValueError("Sentinel-high effect is not positive: " + cru)
        groups[cru] = dict(high_idx=hi, normal_idx=ni, valid=valid)
        group_rows.append(dict(cru_id=cru, sentinel_acr=selected[cru]["sentinel_acr"],
            base_eligible=selected[cru]["eligible"], high_threshold=format(high,".17g"),
            normal_threshold=format(normal,".17g"), n_high=len(hi), n_normal=len(ni),
            high_groups_valid=int(valid), sentinel_delta_high_vs_normal=format(delta,".17g"),
            exclusion_reason="NA" if valid else "overlapping_thresholds_or_small_groups"))
    write_table(root/"data/sentinel_high_groups.tsv", list(group_rows[0]), group_rows)
    pairs = list(read(root/"data/distance_matched_pairs.tsv.gz"))
    features = defaultdict(set)
    for p in pairs:
        cru = p["cru_id"]
        if not groups[cru]["valid"]:
            continue
        features[cru,p["feature_type"],"inside"].add(p["observed_feature"])
        features[cru,p["feature_type"],"control"].add(p["control_feature"])
    matrices = {}
    for kind, label in (("ACR","acr_signal"),("gene","expression")):
        wanted = set().union(*(v for (c,k,r),v in features.items() if k==kind))
        header = read_matrix_header(inputs[label])[1:]
        columns = [header.index(s) for s in samples]
        log("Loading raw HIGH-contrast " + kind + " features: " + str(len(wanted)))
        names, matrix, index = read_selected_matrix(inputs[label], wanted, columns)
        if set(names)!=wanted:
            raise ValueError("Missing raw feature measurements for " + kind)
        matrices[kind] = matrix, index
    fields = ["cru_id","feature_type","role","feature_id","n_high","n_normal",
              "mean_high_log2","mean_normal_log2","delta_log2"]
    count = 0
    with gzip.open(root/"data/feature_effects_high.tsv.gz", "wt", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        for (cru,kind,role), ids in sorted(features.items()):
            matrix, index = matrices[kind]
            names = sorted(ids)
            block = matrix[[index[name] for name in names]]
            hi = np.log2(block[:,groups[cru]["high_idx"]]+1)
            no = np.log2(block[:,groups[cru]["normal_idx"]]+1)
            nh, nn = np.isfinite(hi).sum(axis=1), np.isfinite(no).sum(axis=1)
            mh = np.nansum(hi,axis=1)/np.maximum(nh,1)
            mn = np.nansum(no,axis=1)/np.maximum(nn,1)
            for j, name in enumerate(names):
                if nh[j]<15 or nn[j]<15:
                    continue
                writer.writerow([cru,kind,role,name,int(nh[j]),int(nn[j]),
                    format(mh[j],".17g"),format(mn[j],".17g"),format(mh[j]-mn[j],".17g")])
                count += 1
    write_table(root/"docs/high_contrast_parameters.tsv", ["parameter","value"],
        [dict(parameter=k,value=v) for k,v in dict(
            sentinel_selection="same original maximum Q80-Q20 sentinel", high_quantile=.8,
            normal_quantile=.5, high_rule=">=Q80", normal_rule="<=Q50", min_group_size=15,
            matched_features="same frozen genomic pairs as low contrast", feature_rows=count,
            high_groups_valid=sum(g["valid"] for g in groups.values()),
            contrast="mean(log2(value+1),high)-mean(log2(value+1),lower_half)",
            numpy_version=np.__version__).items()])
    log("Finished true high-vs-lower-half computation: " + str(count) + " feature effects")

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-root", type=Path, default=Path("."),
                        help="Raw input root with the original 01_data layout (default: current directory)")
    parser.add_argument("--outdir", type=Path, required=True, help="Output directory, separate from raw inputs")
    parser.add_argument("--cru-bed", type=Path, help="Override CRU BED4 path")
    parser.add_argument("--acr-bed", type=Path, help="Override ACR BED4 path")
    parser.add_argument("--acr-matrix", type=Path, help="Override tab-separated ACR accessibility matrix")
    parser.add_argument("--expression-matrix", type=Path, help="Override tab-separated gene expression matrix")
    parser.add_argument("--gff", type=Path, help="Override gene annotation GFF3 path")
    args = parser.parse_args()
    try:
        prepare_low(args)
        prepare_high(args)
    except FileNotFoundError as error:
        parser.error(f'Missing raw input: {error}. Set INPUT_ROOT (or --input-root) to the raw 01_data directory; use bash 00_RUN.sh demo for bundled data.')


if __name__ == "__main__":
    main()
