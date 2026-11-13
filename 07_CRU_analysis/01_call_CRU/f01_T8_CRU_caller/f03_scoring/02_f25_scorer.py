#!/usr/bin/env python3
"""Recompute the frozen V11/f25 CRU regulatory score for an arbitrary CRU BED.

The formula, normalization and weights are copied from
f15_fixed_V11_scores/01_freeze_and_score.py (== f25_fixed_V11_300kb_regulatory_score
== f17_V11_pipeline/06_score_V11.py). Only the driving loop changed: the scored
cohort is the BED you pass, so T7-A CRUs can be scored without touching any
coordinate and without any score-based filtering.

Nothing here is allowed to modify the CRU BED: it is hashed before and after.
"""
import argparse
import csv
import gzip
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

COMPONENTS = {'ag_density_component': .40, 'co_density_component': .25,
              'ag_independence_component': .20, 'co_independence_component': .15}


def md5(path):
    h = hashlib.md5()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        yield from csv.DictReader(f, delimiter='\t')


def write(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(dict.fromkeys(k for r in rows for k in r))
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt', newline='') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def read_bed(path, chroms):
    result = {c: [] for c in chroms}
    with Path(path).open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            r = line.split()
            c, a, b, k = r[:4]
            if c not in result:
                continue
            a, b = int(a), int(b)
            if not 0 <= a < b:
                raise ValueError('Invalid CRU interval: ' + line.strip())
            result[c].append({'chrom': c, 'start': a, 'end': b, 'cru_id': k})
    for c in chroms:
        result[c].sort(key=lambda r: r['start'])
        assert all(0 <= r['start'] < r['end'] for r in result[c])
        assert all(a['end'] <= b['start'] for a, b in zip(result[c], result[c][1:])), \
            'CRUs overlap on ' + c
    return result


def metadata(path, chroms, gene=False):
    result, positions = {}, {c: [] for c in chroms}
    with Path(path).open() as f:
        for line in f:
            r = line.split()
            if not r or r[0] not in positions:
                continue
            c, a, b, k = r[:4]
            a, b = int(a), int(b)
            if gene and (len(r) < 6 or r[5] not in ['+', '-']):
                raise ValueError('Reliable gene strand is required')
            p = (b - 1 if r[5] == '-' else a) if gene else (a + b) // 2
            if k in result:
                assert result[k] == (c, p), 'Conflicting annotation ID: ' + k
                continue
            result[k] = (c, p)
            positions[c].append(p)
    return result, {c: np.sort(v) for c, v in positions.items()}


def assign(points, crus):
    if not len(crus):
        return np.full(len(points), -1, dtype=int)
    starts = np.array([r['start'] for r in crus])
    ends = np.array([r['end'] for r in crus])
    idx = np.searchsorted(starts, points, side='right') - 1
    return np.where((idx >= 0) & (np.asarray(points) < ends[np.maximum(idx, 0)]), idx, -1)


def endpoint_counts(left, right, n):
    internal = (left == right) & (left >= 0)
    other = (left >= 0) & (right >= 0) & (left != right)
    count = lambda x: np.bincount(x, minlength=n)
    I = count(left[internal])
    X = count(np.r_[left[other], right[other]])
    O = count(np.r_[left[(left >= 0) & (right < 0)], right[(right >= 0) & (left < 0)]])
    return I, X, O


def normalize(values, low, high):
    z = np.asarray(values, float)
    ok = np.isfinite(z)
    out = np.full(len(z), np.nan)
    if not np.any(ok):
        return out, np.nan, np.nan
    a, b = np.quantile(z[ok], [low, high])
    if b > a:
        out[ok] = np.clip((z[ok] - a) / (b - a), 0, 1)
    else:
        out[ok] = np.where(z[ok] > 0, .5, 0)
    return out, float(a), float(b)


def main(args):
    root = args.outdir.resolve()
    bed = args.cru_bed.resolve()
    before = md5(bed)
    cfg = json.loads(args.scoring_settings.read_text()) if args.scoring_settings else {}
    weights = dict(COMPONENTS) if not cfg.get('weights') else dict(cfg['weights'])
    min_evidence = int(cfg.get('min_evidence_edges', args.min_evidence_edges))
    min_internal = int(cfg.get('min_internal_edges', args.min_internal_edges))
    if abs(sum(weights.values()) - 1) > 1e-9:
        raise ValueError('Component weights must sum to 1')
    if not 0 <= args.norm_low < args.norm_high <= 1:
        raise ValueError('Invalid normalization quantiles')
    sizes = {}
    with args.genome_sizes.open() as f:
        for line in f:
            if line.strip():
                c, n = line.split()[:2]
                sizes[c] = int(n)
    chroms = list(sizes) if args.chroms == ['all'] else args.chroms
    if len(set(chroms)) != len(chroms) or set(chroms) - set(sizes):
        raise ValueError('Invalid chromosome list')
    crus = read_bed(bed, chroms)
    cohort_n = sum(len(v) for v in crus.values())
    if not cohort_n:
        raise ValueError('No CRUs to score')
    acr, ax = metadata(args.acr_bed, chroms)
    gene, gx = metadata(args.gene_bed, chroms, gene=True)
    networks = {name: {c: [] for c in chroms} for name in ['ag', 'co']}
    seen = set()
    dup = 0
    with args.aglink_bedpe.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            c = r['chrom']
            if c not in networks['ag']:
                continue
            key = (r['acr_id'], r['gene_id'])
            if key in seen:
                dup += 1
                continue
            seen.add(key)
            p, q = int(r['p1']), int(r['p2'])
            assert acr[key[0]] == (c, p) and gene[key[1]] == (c, q), 'Anchor does not match annotation'
            networks['ag'][c].append((p, q))
    seen = set()
    for r in read(args.coacr_cache):
        c = r['chrom']
        if c not in networks['co']:
            continue
        p, q = int(r['p1']), int(r['p2'])
        u, v = r['id1'], r['id2']
        cor = float(r['correlation'])
        assert c in chroms and acr[u] == (c, p) and acr[v] == (c, q)
        assert abs(cor) >= .6 and abs(p - q) <= 1000000 and u != v
        key = tuple(sorted((u, v)))
        if key in seen:
            raise ValueError('coACR cache is not unique')
        seen.add(key)
        networks['co'][c].append((p, q))
    rows, qc = [], []
    for c in chroms:
        rr = crus[c]
        n = len(rr)
        if not n:
            continue
        start = np.array([r['start'] for r in rr])
        end = np.array([r['end'] for r in rr])
        na = np.searchsorted(ax[c], end) - np.searchsorted(ax[c], start)
        ng = np.searchsorted(gx[c], end) - np.searchsorted(gx[c], start)
        stats = {}
        for name in ['ag', 'co']:
            ee = networks[name][c]
            if ee:
                p = np.array([e[0] for e in ee])
                q = np.array([e[1] for e in ee])
                l, r = assign(p, rr), assign(q, rr)
                I, X, O = endpoint_counts(l, r, n)
            else:
                I = X = O = np.zeros(n, dtype=int)
                l = q = np.array([], dtype=int)
            stats[name] = (I, X, O)
            qc.append({'chrom': c, 'network': name, 'unique_edges': len(ee),
                       'internal_edges': int(I.sum()),
                       'between_two_CRUs_edges': int(X.sum() // 2),
                       'one_CRU_one_unassigned_edges': int(O.sum()),
                       'both_unassigned_edges': int(np.count_nonzero((l < 0) & (q < 0)))
                       if ee else 0})
        for i, cru in enumerate(rr):
            length = cru['end'] - cru['start']
            row = {**cru, 'length_bp': length, 'acr_n': int(na[i]), 'gene_n': int(ng[i])}
            for name in ['ag', 'co']:
                I, X, O = (int(v[i]) for v in stats[name])
                T = I + X + O
                row.update({name + '_internal_n': I, name + '_other_CRU_n': X,
                            name + '_outside_CRU_n': O, name + '_touching_n': T,
                            name + '_internal_fraction': I / T if T else np.nan,
                            name + '_cross_fraction': (X + O) / T if T else np.nan,
                            name + '_per_mb': I * 1e6 / length,
                            name + '_independence_component': (I + 1) / (T + 2) if T else np.nan})
                nodes = int(na[i] + ng[i]) if name == 'ag' else int(na[i])
                valid = (na[i] > 0 and ng[i] > 0) if name == 'ag' else na[i] >= 2
                row[name + '_per_node_degree'] = 2 * I / nodes if valid else np.nan
                if not valid:
                    row[name + '_per_mb'] = np.nan
            flags = []
            if na[i] < 2:
                flags.append('few_ACRs')
            if ng[i] < 1:
                flags.append('no_gene_TSS')
            for name in ['ag', 'co']:
                if row[name + '_touching_n'] < min_evidence:
                    flags.append(name + '_sparse_incident_edges')
                if row[name + '_internal_n'] < min_internal:
                    flags.append(name + '_sparse_internal_edges')
            row['evidence_status'] = 'evaluable' if not flags else 'insufficient_evidence'
            row['evidence_notes'] = ';'.join(flags)
            rows.append(row)
    normalization = []
    for name in ['ag', 'co']:
        scaled = []
        for quantity in ['per_mb', 'per_node_degree']:
            key = name + '_' + quantity
            z = np.log1p([r[key] for r in rows])
            v, a, b = normalize(z, args.norm_low, args.norm_high)
            scaled.append(v)
            for r, x, log, norm in zip(rows, [r[key] for r in rows], z, v):
                r[key + '_log1p'] = float(log)
                r[key + '_normalized'] = float(norm)
            normalization.append({'metric': key, 'transform': 'log1p', 'q_low': args.norm_low,
                                  'q_high': args.norm_high, 'lower': a, 'upper': b,
                                  'population_n': len(rows),
                                  'finite_n': int(np.isfinite(z).sum())})
        for j, r in enumerate(rows):
            r[name + '_density_component'] = float((scaled[0][j] + scaled[1][j]) / 2)
    for r in rows:
        available = [k for k in weights if np.isfinite(r[k])]
        coverage = sum(weights[k] for k in available)
        r['available_weight_fraction'] = coverage
        r['regulatory_score'] = 100 * sum(weights[k] * r[k] for k in available) / coverage if coverage else np.nan
        r['score_status'] = ('descriptive' if r['evidence_status'] == 'evaluable' and coverage == 1
                             else 'provisional_insufficient_evidence')
        for k, w in weights.items():
            r[k + '_effective_weight'] = w / coverage if k in available and coverage else 0
    values = np.array([r['regulatory_score'] for r in rows])
    ok = np.isfinite(values)
    ranks = np.full(len(rows), np.nan)
    ranks[ok] = rankdata(-values[ok], method='min')
    for r, k in zip(rows, ranks):
        r['rank'] = float(k)
    assert np.all((values[ok] >= 0) & (values[ok] <= 100 + 1e-9))
    root.mkdir(parents=True, exist_ok=True)
    for d in ['f01_scores', 'f02_summary', 'f03_formula']:
        (root / d).mkdir(exist_ok=True)
    write(root / 'f01_scores/01_CRU_regulatory_scores.tsv', rows)
    write(root / 'f01_scores/02_normalization_parameters.tsv', normalization)
    write(root / 'f01_scores/03_endpoint_assignment_QC.tsv', qc)
    write(root / 'f01_scores/04_CRU_regulatory_scores_simple.tsv',
          ({'chrom': r['chrom'], 'start': r['start'], 'end': r['end'], 'cru_id': r['cru_id'],
            'length_bp': r['length_bp'], 'regulatory_score': r['regulatory_score'],
            'score_status': r['score_status'], 'rank': r['rank']} for r in rows))
    summary = []
    for scope in chroms + ['Tested_chromosomes']:
        rr = [r for r in rows if scope == 'Tested_chromosomes' or r['chrom'] == scope]
        if not rr:
            continue
        v = np.array([r['regulatory_score'] for r in rr], float)
        summary.append({'scope': scope, 'cru_n': len(rr),
                        'finite_score_n': int(np.isfinite(v).sum()),
                        'median_score': float(np.nanmedian(v)) if np.isfinite(v).any() else np.nan,
                        'min_score': float(np.nanmin(v)) if np.isfinite(v).any() else np.nan,
                        'max_score': float(np.nanmax(v)) if np.isfinite(v).any() else np.nan,
                        'insufficient_evidence_n': sum(r['evidence_status'] != 'evaluable' for r in rr)})
    write(root / 'f02_summary/05_score_summary.tsv', summary)
    params = {'label': args.label, 'cru_bed': str(bed), 'cru_bed_md5': before,
              'cohort_n': cohort_n, 'chroms': chroms,
              'aglink_bedpe': str(args.aglink_bedpe), 'coacr_cache': str(args.coacr_cache),
              'gene_bed': str(args.gene_bed), 'acr_bed': str(args.acr_bed),
              'genome_sizes': str(args.genome_sizes),
              'weights': weights, 'normalization': {'transform': 'log1p',
                                                    'quantiles': [args.norm_low, args.norm_high],
                                                    'population': 'all scored CRUs of this BED'},
              'minimum': {'min_evidence_edges': min_evidence, 'min_internal_edges': min_internal},
              'formula': 'per network X in {ag,co}: X_density = mean(z(X_per_mb), z(X_per_node_degree)); '
                         'X_independence = (internal+1)/(touching+2); '
                         'score = 100*sum(w*X)/sum(w over available)',
              'inputs_md5': {str(p): md5(p) for p in
                             [args.aglink_bedpe, args.coacr_cache, args.gene_bed, args.acr_bed,
                              args.genome_sizes]}}
    (root / '00_scoring_parameters.json').write_text(json.dumps(params, indent=2) + '\n')
    (root / 'f03_formula/FORMULA.md').write_text(FORMULA.format(**{
        'label': args.label, 'cohort_n': cohort_n, 'low': args.norm_low, 'high': args.norm_high,
        'min_ev': min_evidence, 'min_in': min_internal,
        'w_ag': weights['ag_density_component'], 'w_co': weights['co_density_component'],
        'w_agi': weights['ag_independence_component'],
        'w_ci': weights['co_independence_component'],
        'cru_bed': bed, 'bed_md5': before}))
    if md5(bed) != before:
        raise ValueError('CRU BED changed during scoring')
    print('scored', cohort_n, 'CRUs from', bed)
    for r in summary:
        if r['scope'] == 'Tested_chromosomes':
            print('SCORED', r, flush=True)


FORMULA = """# {label} CRU regulatory score (recomputed)

CRU BED: `{cru_bed}` (md5 `{bed_md5}`), cohort n = {cohort_n}.
The BED is hashed before and after scoring and is never modified; no CRU is
dropped, moved or filtered by score.

## 1. Raw endpoint metrics (per CRU)

For each network `X` in {{ag, co}} (ag = agLINK, co = coACR):

* Every edge has two genomic endpoints; each endpoint is *assigned* to the CRU
  whose `[start, end)` contains it, or to `none`.
* `I` = edges with both endpoints in this CRU (`internal`)
* `X` = endpoint slots of edges whose two endpoints sit in two different CRUs
* `O` = endpoint slots of edges with exactly one endpoint assigned
* `T = I + X + O`  (`touching`)

```
X_per_mb              = I * 1e6 / length_bp
X_internal_fraction   = I / T
X_cross_fraction      = (X + O) / T
X_independence        = (I + 1) / (T + 2)          # Laplace-smoothed
X_per_node_degree     = 2 * I / nodes              # nodes = acr_n+gene_n (ag), acr_n (co)
```

`nodes` is only used when the CRU is meaningful for that network: `ag` requires
`acr_n > 0 and gene_n > 0`; `co` requires `acr_n >= 2`. Otherwise
`X_per_node_degree` and `X_per_mb` are undefined (NaN).

Also recorded per CRU: `acr_n`, `gene_n` (annotation midpoints/TSS inside the CRU).

## 2. Normalization (cohort-relative, computed once)

For every metric above (`ag_per_mb`, `ag_per_node_degree`, `co_per_mb`,
`co_per_node_degree`):

```
z      = log1p(metric)
lower  = quantile(z, {low})        # over the finite values of this cohort
upper  = quantile(z, {high})
norm   = clip((z - lower) / (upper - lower), 0, 1)     # if upper > lower
norm   = 0.5 where z > 0 else 0                        # if upper == lower
```

The population is exactly the CRU set of the supplied BED, so scores are
cohort-relative and must not be compared across different cohorts.

## 3. Components

```
ag_density_component        = (norm(ag_per_mb) + norm(ag_per_node_degree)) / 2
co_density_component        = (norm(co_per_mb) + norm(co_per_node_degree)) / 2
ag_independence_component   = (ag_I + 1) / (ag_T + 2)
co_independence_component   = (co_I + 1) / (co_T + 2)
```

If either term of a density component is NaN the component is NaN.

## 4. Final score

```
weights = {{ag_density: {w_ag}, co_density: {w_co},
           ag_independence: {w_agi}, co_independence: {w_ci}}}     # sum = 1
available = components that are finite
coverage  = sum of weights of available components
regulatory_score = 100 * sum(w_k * comp_k for available k) / coverage
```

`available_weight_fraction` = coverage; `*_effective_weight` = `w_k / coverage`;
`score_status = descriptive` only when `evidence_status == evaluable` and
`coverage == 1`, otherwise `provisional_insufficient_evidence`.
`rank` = 1 for the highest score (`rankdata(-score, method='min')`).

## 5. Evidence flags (informational only, no filtering)

A CRU is `insufficient_evidence` if any of: `acr_n < 2` (`few_ACRs`),
`gene_n < 1` (`no_gene_TSS`), `X_touching < {min_ev}` (`X_sparse_incident_edges`),
`X_internal < {min_in}` (`X_sparse_internal_edges`) for X in {{ag, co}}.
Such CRUs are still scored with the renormalized weights above.
"""


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cru-bed', type=Path, required=True)
    p.add_argument('--aglink-bedpe', type=Path, required=True,
                   help='calling/f01_shared_input/01_count_only_aglink.annotated.bedpe')
    p.add_argument('--coacr-cache', type=Path, required=True,
                   help='coACR input dir 07_selected_coacr_pairs.tsv.gz')
    p.add_argument('--gene-bed', type=Path, required=True)
    p.add_argument('--acr-bed', type=Path, required=True)
    p.add_argument('--genome-sizes', type=Path, required=True)
    p.add_argument('--outdir', type=Path, required=True)
    p.add_argument('--scoring-settings', type=Path, default=None,
                   help='JSON with weights/min_evidence_edges/min_internal_edges (frozen defaults if omitted)')
    p.add_argument('--min-evidence-edges', type=int, default=5)
    p.add_argument('--min-internal-edges', type=int, default=3)
    p.add_argument('--norm-low', type=float, default=.05)
    p.add_argument('--norm-high', type=float, default=.95)
    p.add_argument('--label', default='CRU')
    p.add_argument('--chroms', nargs='+', default=['all'])
    main(p.parse_args())
