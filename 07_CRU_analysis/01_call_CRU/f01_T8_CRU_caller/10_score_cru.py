#!/usr/bin/env python3
"""Score fixed CRU intervals with the frozen f25 formula and fixed T7-A anchors.

The scoring formula, component weights, evidence flags and ranking are the ones
shipped in `f03_scoring/02_f25_scorer.py`, executed unmodified. Only the
normalization callback is replaced, so that every metric is rescaled with the
T7-A log1p q05/q95 anchors in `f03_scoring/01_norm_anchors.tsv` instead of being
re-estimated on the cohort being scored. CRU coordinates are never touched.

Two networks are scored: agLINK (ACR-gene) and co-accessibility (ACR-ACR), each
as a density component and an independence component.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
from datetime import datetime
from pathlib import Path
from types import SimpleNamespace as NS
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

METRICS = ("ag_per_mb", "ag_per_node_degree", "co_per_mb", "co_per_node_degree")
WEIGHTS = {"ag_density_component": .40, "co_density_component": .25,
           "ag_independence_component": .20, "co_independence_component": .15}
MIN_EVIDENCE, MIN_INTERNAL = 5, 3
COACR_MIN_ABS_PCC, COACR_MAX_DISTANCE_BP = .6, 1_000_000


def log(*items):
    print(datetime.now().isoformat(timespec="seconds"), *items, flush=True)


def md5(path):
    h = hashlib.md5()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    with Path(path).open() as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write(path, rows, fields=None):
    rows = list(rows)
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        Path(path).write_text("")
        return
    names = fields or list(dict.fromkeys(k for row in rows for k in row))
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, names, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def coords(row):
    return row["chrom"], int(row["start"]), int(row["end"])


def load_scorer():
    """Load the shipped, unedited f25 scorer by path."""
    path = HERE / "f03_scoring/02_f25_scorer.py"
    spec = importlib.util.spec_from_file_location("f25_scorer", str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_anchors(path):
    rows = read(path)
    anchors = {row["metric"]: row for row in rows}
    if set(anchors) != set(METRICS) or len(rows) != len(METRICS):
        raise ValueError("Anchor table must hold exactly the four density metrics")
    for row in rows:
        if row["transform"] != "log1p":
            raise ValueError("Anchors must be on the log1p scale")
    return anchors


class FrozenNormalize:
    """Replace the scorer's own quantiles with the fixed T7-A anchors."""

    def __init__(self, anchors):
        self.anchors = anchors
        self.calls = []

    def __call__(self, values, low, high):
        metric = METRICS[len(self.calls)]
        anchor = self.anchors[metric]
        if not (np.isclose(low, float(anchor["q_low"])) and np.isclose(high, float(anchor["q_high"]))):
            raise ValueError("Scorer requested unexpected quantiles for " + metric)
        z = np.asarray(values, float)
        ok = np.isfinite(z)
        out = np.full(len(z), np.nan)
        a, b = float(anchor["lower"]), float(anchor["upper"])
        if b > a:
            out[ok] = np.clip((z[ok] - a) / (b - a), 0, 1)
        else:
            out[ok] = np.where(z[ok] > 0, .5, 0)
        self.calls.append(dict(metric=metric, target_population_n=len(z),
                               target_finite_n=int(ok.sum())))
        return out, a, b


def build_aglink_bedpe(aglink, gene_bed, acr_bed, sizes, chroms, outdir):
    """Rebuild the annotated, ID-bearing agLINK table with the frozen loader."""
    from f02_modules import m01_density_signal as base
    inputs = NS(aglink=aglink, gene_bed=gene_bed, acr_bed=acr_bed,
                chroms=list(chroms), bedpe_reference=None)
    base.load_count_only_edges(inputs, {c: sizes[c] for c in chroms}, outdir)
    return outdir / "01_count_only_aglink.annotated.bedpe"


def build_coacr_cache(coacr, acr_bed, chroms, outpath):
    """Keep same-chromosome pairs with |PCC| >= .6 and <= 1 Mb, strongest per pair."""
    meta = {}
    with Path(acr_bed).open() as handle:
        for line in handle:
            row = line.split()
            if len(row) < 4 or row[0] not in chroms:
                continue
            meta[row[3]] = (row[0], (int(row[1]) + int(row[2])) // 2)
    best, audit = {}, dict(input_rows=0, kept=0)
    with Path(coacr).open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"id1", "id2", "correlation"} <= set(reader.fieldnames or []):
            raise ValueError("co-accessibility input needs id1, id2, correlation columns")
        for row in reader:
            audit["input_rows"] += 1
            a, b = row["id1"], row["id2"]
            if a == b or a not in meta or b not in meta:
                continue
            chrom, p = meta[a]
            other, q = meta[b]
            value = float(row["correlation"])
            if chrom != other or abs(value) < COACR_MIN_ABS_PCC or abs(p - q) > COACR_MAX_DISTANCE_BP:
                continue
            key = tuple(sorted((a, b)))
            if key not in best or abs(value) > abs(best[key][3]):
                best[key] = (chrom, p, q, value)
    rows = []
    for (a, b), (chrom, p, q, value) in best.items():
        id1, id2, p1, p2 = (a, b, p, q) if p <= q else (b, a, q, p)
        rows.append(dict(chrom=chrom, id1=id1, id2=id2, p1=p1, p2=p2, correlation=value))
    rows.sort(key=lambda r: (r["chrom"], r["p1"], r["p2"]))
    audit["kept"] = len(rows)
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(outpath, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, ["chrom", "id1", "id2", "p1", "p2", "correlation"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    return outpath, audit


def describe(label, rows, scope="Genomewide"):
    values = np.array([float(r["regulatory_score"]) for r in rows])
    finite = values[np.isfinite(values)]
    result = dict(dataset=label, scope=scope, cru_n=len(rows), finite_score_n=len(finite),
                  evaluable_cru_n=sum(r["evidence_status"] == "evaluable" for r in rows),
                  descriptive_score_n=sum(r["score_status"] == "descriptive" for r in rows))
    result["evaluable_cru_fraction"] = result["evaluable_cru_n"] / len(rows)
    result["mean_score"] = float(np.mean(finite))
    for name, q in (("min", 0), ("p05", .05), ("p25", .25), ("median", .5),
                    ("p75", .75), ("p95", .95), ("max", 1)):
        result[name + "_score"] = float(np.quantile(finite, q))
    return result


def compare_tables(original, replay, output):
    old, new = read(original), read(replay)
    a, b = {coords(r): r for r in old}, {coords(r): r for r in new}
    if not (len(a) == len(old) and len(b) == len(new) and a.keys() == b.keys()):
        raise ValueError("Reference replay covers different CRUs; cannot validate the scorer")
    if set(old[0]) != set(new[0]):
        raise ValueError("Reference replay wrote different columns")
    checks = []
    for field in old[0]:
        errors, maximum = 0, 0.
        for key in a:
            x, y = a[key][field], b[key][field]
            try:
                u, v = float(x), float(y)
                good = bool(np.isclose(u, v, rtol=1e-12, atol=1e-10, equal_nan=True))
                if np.isfinite(u) and np.isfinite(v):
                    maximum = max(maximum, abs(u - v))
            except ValueError:
                good = x == y
            errors += int(not good)
        checks.append(dict(field=field, compared_cru_n=len(a), mismatch_n=errors,
                           max_absolute_numeric_difference=maximum))
    write(output, checks)
    bad = [r for r in checks if r["mismatch_n"]]
    if bad:
        raise ValueError("Reference replay differs; refusing to score: " + str(bad[:3]))
    return len(a), max(r["max_absolute_numeric_difference"] for r in checks)


def settings_file(directory):
    path = directory / "f01_inputs/00_scoring_settings.json"
    write_json(path, dict(weights=WEIGHTS, min_evidence_edges=MIN_EVIDENCE,
                          min_internal_edges=MIN_INTERNAL,
                          normalization=dict(transform="log1p", quantiles=[.05, .95],
                                             population="fixed T7-A anchor cohort"),
                          formula="density = mean(z(per_mb), z(per_node_degree)); "
                                  "independence = (internal+1)/(touching+2)"))
    return path


def write_json(path, payload):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def main(args):
    root = args.outdir.resolve()
    root.mkdir(parents=True, exist_ok=True)
    sizes = {}
    with args.genome_sizes.open() as handle:
        for line in handle:
            if line.strip():
                chrom, size = line.split()[:2]
                sizes[chrom] = int(size)
    seen = []
    for line in args.cru_bed.read_text().splitlines():
        if line.strip() and line.split()[0] in sizes and line.split()[0] not in seen:
            seen.append(line.split()[0])
    chroms = args.chroms or seen
    if not chroms:
        raise ValueError("No chromosome in the CRU BED matches the genome sizes file")
    log("chromosomes:", ",".join(chroms))

    inputs = root / "f01_inputs"
    inputs.mkdir(exist_ok=True)
    log("rebuilding the annotated agLINK table")
    bedpe = build_aglink_bedpe(args.aglink, args.gene_bed, args.acr_bed, sizes, chroms, inputs)
    if args.coacr is not None:
        log("rebuilding the co-accessibility cache")
        cache, audit = build_coacr_cache(args.coacr, args.acr_bed, set(chroms),
                                         inputs / "02_coacr_cache.tsv.gz")
        log("co-accessibility pairs: %d kept of %d" % (audit["kept"], audit["input_rows"]))
    else:
        cache = inputs / "02_coacr_cache.tsv.gz"
        with gzip.open(cache, "wt", newline="") as handle:
            handle.write("chrom\tid1\tid2\tp1\tp2\tcorrelation\n")
        audit = dict(input_rows=0, kept=0)
        print("warning: no --coacr given, the co-accessibility network will be empty", flush=True)

    anchors = load_anchors(args.anchors)
    scorer = load_scorer()
    if scorer.COMPONENTS != WEIGHTS:
        raise ValueError("Frozen scorer weights changed")

    def options(bed, out, label):
        return NS(outdir=out, cru_bed=bed, scoring_settings=settings_file(root),
                  min_evidence_edges=MIN_EVIDENCE, min_internal_edges=MIN_INTERNAL,
                  norm_low=.05, norm_high=.95, genome_sizes=args.genome_sizes,
                  chroms=chroms, acr_bed=args.acr_bed, gene_bed=args.gene_bed,
                  aglink_bedpe=bedpe, coacr_cache=cache, label=label)

    difference = None
    if args.reference_bed is not None and args.reference_scores is not None:
        log("stage 1: replay the unedited scorer on the reference cohort")
        replay_dir = root / "f02_reference_replay"
        scorer.main(options(args.reference_bed, replay_dir, "reference_replay"))
        expected = read(args.reference_scores)
        if len(expected) and "regulatory_score" in expected[0]:
            difference = compare_tables(args.reference_scores,
                                        replay_dir / "f01_scores/01_CRU_regulatory_scores.tsv",
                                        replay_dir / "01_replay_checks.tsv")[1]
            log("PASS replay: %d reference CRUs agree; max numeric difference %s"
                % (len(expected), difference))
        recovered = {r["metric"]: r for r in read(replay_dir / "f01_scores/02_normalization_parameters.tsv")}
        for metric in METRICS:
            for key in ("lower", "upper"):
                if not np.isclose(float(recovered[metric][key]), float(anchors[metric][key]),
                                  rtol=1e-12, atol=1e-12):
                    raise ValueError("Anchor is not reproduced by the frozen scorer: " + metric)

    log("stage 2: score the fixed CRUs with the frozen T7-A anchors")
    target = root / "f03_scores"
    frozen_normalize = FrozenNormalize(anchors)
    scorer.normalize = frozen_normalize
    scorer.main(options(args.cru_bed, target, "fixed_CRU_f25"))
    if [c["metric"] for c in frozen_normalize.calls] != list(METRICS):
        raise ValueError("Frozen normalization was not used for all four metrics")

    rows = read(target / "f01_scores/01_CRU_regulatory_scores.tsv")
    if args.expected_cru_n and len(rows) != args.expected_cru_n:
        raise ValueError("Expected %d scored CRUs, got %d" % (args.expected_cru_n, len(rows)))
    compact = []
    with (root / "02_CRU_scored.bed").open("w") as handle:
        for row in rows:
            score = float(row["regulatory_score"])
            if not (0 <= score <= 100 + 1e-9):
                raise ValueError("Score out of range: " + str(score))
            entry = {k: row[k] for k in ("cru_id", "chrom", "start", "end", "length_bp", "gene_n", "acr_n")}
            for name in ("ag", "co"):
                for key in ("internal_n", "touching_n", "per_mb", "per_node_degree"):
                    entry[name + "_" + key] = row[name + "_" + key]
                entry[name + "_independence_raw"] = row[name + "_independence_component"]
                for key in ("density_component", "independence_component"):
                    entry[name + "_" + key + "_score_0_100"] = 100 * float(row[name + "_" + key])
            entry.update({k: row[k] for k in ("available_weight_fraction", "regulatory_score",
                                              "rank", "evidence_status", "score_status", "evidence_notes")})
            compact.append(entry)
            handle.write("\t".join([row["chrom"], row["start"], row["end"], row["cru_id"],
                                    str(int(round(score * 10))), "."]) + "\n")
    write(root / "01_CRU_scores.tsv", compact)
    simple = []
    for row in rows:
        entry = {k: row[k] for k in ("chrom", "start", "end", "cru_id")}
        entry["score"] = row["regulatory_score"]
        simple.append(entry)
    write(root / "04_CRU_scores_simple.tsv", simple)

    summaries = [describe("scored_cohort", rows)]
    summary_dir = root / "f04_summary"
    write(summary_dir / "01_score_summary.tsv", summaries)
    write(summary_dir / "02_score_by_chromosome.tsv",
          [describe("scored_cohort", [r for r in rows if r["chrom"] == c], c) for c in chroms])
    counts, bins = np.histogram([float(r["regulatory_score"]) for r in rows], bins=np.arange(0, 105, 5))
    write(summary_dir / "03_score_distribution.tsv",
          [dict(score_bin_start=int(a), score_bin_end=int(b), cru_n=int(n), fraction=float(n) / len(rows))
           for a, b, n in zip(bins, bins[1:], counts)])
    write(root / "05_input_MD5.tsv",
          [dict(input=k, path=str(p), md5=md5(p)) for k, p in
           [("cru_bed", args.cru_bed), ("aglink", args.aglink), ("coacr", args.coacr),
            ("gene_bed", args.gene_bed), ("acr_bed", args.acr_bed), ("genome_sizes", args.genome_sizes)]])
    log("COMPLETE %d scored CRUs" % len(rows))
    for r in summaries:
        log(r)
    return rows


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--cru-bed", type=Path, required=True, help="merged CRU.bed from 01_run_cru.py")
    p.add_argument("--aglink", type=Path, required=True, help="agLINK table: gene, ACR, PCC, p")
    p.add_argument("--acr-bed", type=Path, required=True)
    p.add_argument("--gene-bed", type=Path, required=True)
    p.add_argument("--genome-sizes", type=Path, required=True)
    p.add_argument("--coacr", type=Path, default=None,
                   help="co-accessibility table: id1, id2, correlation (optional)")
    p.add_argument("--outdir", type=Path, required=True)
    p.add_argument("--anchors", type=Path, default=HERE / "f03_scoring/01_norm_anchors.tsv")
    p.add_argument("--reference-bed", type=Path, default=None,
                   help="optional: original CRU BED replayed before scoring, for validation")
    p.add_argument("--reference-scores", type=Path, default=None,
                   help="optional: score table the replay must reproduce")
    p.add_argument("--chroms", nargs="+", default=None)
    p.add_argument("--expected-cru-n", type=int, default=0)
    return p


if __name__ == "__main__":
    main(parser().parse_args())
