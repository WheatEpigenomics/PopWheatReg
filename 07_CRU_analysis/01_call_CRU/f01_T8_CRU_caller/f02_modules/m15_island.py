# m15_island: T7-A island rescue
# Frozen rule module, copied unchanged from the published calling chain
# (source: f21_T7A_pipeline/01_call_T7A.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
from ._load import load as _load
"""Genome-wide T7-A/T7-B low-background island rescue, packaged like f17 V11.

Semantics are byte-for-byte the f20_T7_island_rescue/01_test_T7_islands.py
algorithm: identical generate/shape/decide/summarise helpers, identical
settings, identical T6 baseline handling. Only the driver changed:

  * one chromosome per worker (multiprocessing, default 7 workers),
  * per-chromosome checkpoints with MD5 inventory so --resume is exact,
  * one merged result set instead of a single-process selected-chromosome run,
  * every chromosome's candidate decisions, density and agLINK tracks are
    written so 03_plot_T7_windows.R can render the 5-Mb atlas.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
import multiprocessing
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np

HERE = Path(__file__).resolve().parent
from . import m14_support_recall as S
P, STEP = S.P, S.STEP
METHODS = ['T6', 'T7_A', 'T7_B']
CTX = {}


def md5(path):
    h = hashlib.md5()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


class CoIndex:
    def __init__(self, rows):
        unique = {}
        for r in rows:
            key = tuple(sorted((r['id1'], r['id2'])))
            unique[key] = r
        self.rows = list(unique.values())
        self.p1 = np.array([r['p1'] for r in self.rows], dtype=np.int64)
        self.p2 = np.array([r['p2'] for r in self.rows], dtype=np.int64)

    def support(self, start, end, ag_acrs):
        inside = (self.p1 >= start) & (self.p1 < end) & (self.p2 >= start) & (self.p2 < end)
        rows = [self.rows[i] for i in np.flatnonzero(inside)]
        parent = {}

        def find(x):
            parent.setdefault(x, x)
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x

        for r in rows:
            parent[find(r['id1'])] = find(r['id2'])
        nodes, counts = {}, {}
        for r in rows:
            k = find(r['id1'])
            nodes.setdefault(k, set()).update([r['id1'], r['id2']])
            counts[k] = counts.get(k, 0) + 1
        connected = [(counts[k], len(nodes[k])) for k in counts if nodes[k] & ag_acrs]
        positions = {r[k] for r in rows for k in ['p1', 'p2']}
        return dict(coacr_internal_n=len(rows),
                    coacr_acr_n=len({r[k] for r in rows for k in ['id1', 'id2']}),
                    coacr_linked_component_edges=max((v[0] for v in connected), default=0),
                    coacr_linked_component_nodes=max((v[1] for v in connected), default=0),
                    coacr_short500_n=sum(abs(r['p1'] - r['p2']) <= 500000 for r in rows),
                    coacr_endpoint_extent_bp=max(positions) - min(positions) if positions else 0)


def load_coacr(path, chroms):
    rows = {c: [] for c in chroms}
    with gzip.open(path, 'rt') as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['chrom'] not in rows:
                continue
            r['p1'], r['p2'] = int(r['p1']), int(r['p2'])
            r['correlation'] = float(r['correlation'])
            if abs(r['correlation']) < .6 - 1e-12 or abs(r['p1'] - r['p2']) > 1000000:
                raise ValueError('coACR cache does not follow the frozen selection')
            rows[r['chrom']].append(r)
    return {c: CoIndex(rr) for c, rr in rows.items()}


def generate(raw, smooth, bg, size, occupied, zero_gap, opt):
    # Numerical noise floor only; no chromosome-wide positive-density quantile.
    active = (raw > opt['raw_noise_floor']) & (smooth[0] > opt['fine_noise_floor'])
    active &= ~occupied & ~zero_gap
    active = P.close_supported(active, ~occupied & ~zero_gap, opt['close_hole_bp'])
    return [(a, b) for a, b in P.runs(active)]


def shape(a, b, raw, smooth, bg, size, opt):
    row = P.describe(raw, smooth, bg, np.zeros(len(raw)), a, b, size)
    flank_n = opt['flank_bp'] // STEP
    left = np.arange(max(0, a - flank_n), a)
    right = np.arange(b, min(len(raw), b + flank_n))
    reference = np.maximum(smooth[0], smooth[1])
    level = float(np.quantile(reference[a:b], .5))
    threshold = max(opt['flank_absolute_near_zero'], opt['flank_relative_max'] * level)
    row.update(island_level=level, flank_low_threshold=threshold,
               sustained_fraction=float(np.mean((raw[a:b] > opt['raw_noise_floor']) &
                                                (smooth[1, a:b] > bg[a:b] +
                                                 opt['background_excess_sigma'] * np.sqrt(bg[a:b] + 1)))),
               fine_peak_excess=float(np.max(smooth[0, a:b] - bg[a:b])))
    for side, ii in [('left', left), ('right', right)]:
        row[side + '_flank_bp'] = len(ii) * STEP
        row[side + '_flank_median'] = float(np.median(reference[ii])) if len(ii) else np.nan
        row[side + '_flank_p75'] = float(np.quantile(reference[ii], .75)) if len(ii) else np.nan
        row[side + '_low_fraction'] = float(np.mean(reference[ii] <= threshold)) if len(ii) else 0.
        row[side + '_contrast_regularized'] = (level + .1) / (row[side + '_flank_median'] + .1) if len(ii) else 0.
    edge_n = min(opt['edge_transition_bp'] // STEP, b - a)
    left_inner, right_inner = float(np.mean(smooth[0, a:a + edge_n])), float(np.mean(smooth[0, b - edge_n:b]))
    left_outer = float(np.mean(smooth[0, max(0, a - edge_n):a])) if a else np.nan
    right_outer = float(np.mean(smooth[0, b:min(len(raw), b + edge_n)])) if b < len(raw) else np.nan
    row.update(left_edge_change=left_inner - left_outer, right_edge_change=right_inner - right_outer)
    return row


def decide(row, opt):
    shared = []
    for key, requirement in [('length_bp', 'min_length_bp'), ('sustained_fraction', 'min_sustained_fraction'),
                             ('internal_span_covered_fraction', 'min_internal_span_fraction')]:
        if row[key] < opt[requirement]:
            shared.append(key)
    if row['max_zero_run_bp'] > opt['max_zero_run_bp']:
        shared.append('wide_zero_gap')
    if row['protected_zero_overlap_bp']:
        shared.append('protected_hard_gap')
    for side in ['left', 'right']:
        if row[side + '_flank_bp'] < opt['min_flank_bp']:
            shared.append(side + '_flank_unobservable')
        if row[side + '_low_fraction'] < opt['min_flank_low_fraction']:
            shared.append(side + '_flank_not_low')
        if row[side + '_flank_p75'] > row['flank_low_threshold']:
            shared.append(side + '_flank_p75')
        if row[side + '_contrast_regularized'] < opt['min_island_flank_contrast']:
            shared.append(side + '_insufficient_contrast')
        if not np.isfinite(row[side + '_edge_change']) or row[side + '_edge_change'] <= opt['min_edge_change']:
            shared.append(side + '_no_density_transition')
    a = list(shared)
    if row['independent_pair_n'] < 2:
        a.append('minimum_two_internal_pairs')
    if row['largest_component_edges'] < 2 and row['endpoint_100kb_max_unique_nodes'] < 3:
        a.append('no_internal_network_or_endpoint_cluster')
    if row['local_1mb_internal_n'] < 1:
        a.append('no_local_internal_aglink')
    b = list(shared)
    if row['independent_pair_n'] != 1:
        b.append('one_pair_rescue_requires_exactly_one_internal_pair')
    if row['single_aglink_span_bp'] > opt['single_aglink_max_bp']:
        b.append('isolated_long_aglink')
    for key, requirement in [('coacr_internal_n', 'min_coacr_pairs'), ('coacr_acr_n', 'min_coacr_acrs'),
                             ('coacr_linked_component_edges', 'min_coacr_pairs'),
                             ('coacr_linked_component_nodes', 'min_coacr_acrs'),
                             ('coacr_short500_n', 'min_short_coacr'),
                             ('coacr_endpoint_extent_bp', 'min_coacr_extent_bp')]:
        if row[key] < opt[requirement]:
            b.append(key)
    return a, b


def summarise(method, c, rows, baseline, size):
    values = np.array([r['end'] - r['start'] for r in rows], dtype=np.int64)
    return dict(method=method, chrom=c, chromosome_bp=size, cru_n=len(rows), t6_cru_n=len(baseline),
                added_cru_n=len(rows) - len(baseline), covered_bp=int(values.sum()),
                coverage_fraction=float(values.sum() / size),
                median_bp=float(np.median(values)) if len(values) else np.nan,
                mean_bp=float(np.mean(values)) if len(values) else np.nan,
                p25_bp=float(np.quantile(values, .25)) if len(values) else np.nan,
                p75_bp=float(np.quantile(values, .75)) if len(values) else np.nan,
                p90_bp=float(np.quantile(values, .9)) if len(values) else np.nan,
                p95_bp=float(np.quantile(values, .95)) if len(values) else np.nan,
                max_bp=int(values.max()) if len(values) else 0,
                length_ge5mb_n=int((values >= 5000000).sum()),
                length_ge10mb_n=int((values >= 10000000).sum()))


def worker(chrom):
    root, source, t6 = CTX['root'], CTX['source'], CTX['t6_root']
    d = root / 'chromosomes' / chrom
    d.mkdir(parents=True, exist_ok=True)
    done = d / 'DONE.json'
    if done.exists() and CTX['resume']:
        old = json.loads(done.read_text())
        if old.get('signature') == CTX['signature'] and all(
                Path(p).is_file() and md5(p) == v for p, v in old['artifacts'].items()):
            return chrom, 'reused'
    opt, size, es = CTX['opt'], CTX['sizes'][chrom], CTX['edges'][chrom]
    index = S.EdgeIndex(es)
    baseline = sorted(P.read_bed(t6 / 'f03_tests/T6' / f'{chrom}.CRU.bed'), key=lambda r: r['start'])
    assert P.disjoint(baseline)
    with np.load(source / 'f02_inputs' / f'{chrom}.profiles.npz') as z:
        raw, smooth, bg = z['full_raw'], z['full_smooth'], z['full_background']
    protected = np.zeros(len(raw), bool)
    for a, b in P.runs(raw <= 1e-9):
        if min(size, b * STEP) - a * STEP >= 100000:
            protected[a:b] = True
    gaps = S.table(source / 'chromosomes' / chrom / 'refined_hard_gaps.tsv')
    for g in gaps:
        g['start'], g['end'] = int(g['start']), int(g['end'])
    gap_mask = P.mask_intervals(gaps, size)
    occupied = P.mask_intervals(baseline, size)
    additions = {'T7_A': [], 'T7_B': []}
    decisions, added = [], []
    added_n = 0
    for a, b in generate(raw, smooth, bg, size, occupied, protected, opt):
        row = dict(chrom=chrom, **shape(a, b, raw, smooth, bg, size, opt))
        row.update(index.describe(row['start'], row['end']))
        for prefix in ['raw', 'unique']:
            denom = sum(row[prefix + '_' + k + '_n'] for k in ['both_internal', 'single_internal', 'span_only'])
            row[prefix + '_region_edge_n'] = denom
            for k in ['both_internal', 'single_internal', 'span_only']:
                row[prefix + '_' + k + '_fraction'] = row[prefix + '_' + k + '_n'] / denom if denom else np.nan
        internal = [e for e in index.unique if row['start'] <= e['left'] and e['right'] < row['end']]
        row.update(single_aglink_span_bp=internal[0]['right'] - internal[0]['left'] if len(internal) == 1 else 0,
                   protected_zero_overlap_bp=int(protected[a:b].sum() * STEP),
                   original_hard_gap_fraction=float(gap_mask[a:b].mean()),
                   annotated_gene_n=int(np.searchsorted(CTX['gp'][chrom], row['end']) -
                                        np.searchsorted(CTX['gp'][chrom], row['start'])),
                   annotated_acr_n=int(np.searchsorted(CTX['ap'][chrom], row['end']) -
                                       np.searchsorted(CTX['ap'][chrom], row['start'])))
        row.update(CTX['co'][chrom].support(row['start'], row['end'], {e['acr_id'] for e in internal}))
        fail_a, fail_b = decide(row, opt)
        row.update(accepted_A=int(not fail_a), accepted_B=int(not fail_a or not fail_b),
                   reject_A=','.join(fail_a) or 'accepted_independent_aglink_island',
                   reject_B=','.join(fail_b) if fail_a and fail_b else
                   ('accepted_A_inherited' if not fail_a else 'accepted_one_aglink_with_coacr'),
                   method='T7_A' if not fail_a else 'T7_B', accepted=int(not fail_a or not fail_b))
        row['reject_reason'] = row['reject_A'] if not fail_a else row['reject_B']
        decisions.append(row)
        if not fail_a or not fail_b:
            added_n += 1
            bedrow = dict(chrom=chrom, start=row['start'], end=row['end'],
                          cru_id=f'{chrom}_T7_island_{added_n}')
            additions['T7_B'].append(bedrow)
            if not fail_a:
                additions['T7_A'].append(bedrow.copy())
            added.append(dict(**row, cru_id=bedrow['cru_id'],
                              acceptance_route='A_independent_aglink' if not fail_a else 'B_one_aglink_coacr'))
    summaries, invariants, support_rows = [], [], []
    allbeds = {}
    for method in METHODS:
        rows = sorted(baseline + additions.get(method, []), key=lambda r: r['start'])
        assert P.disjoint(rows) and all(r in rows for r in baseline), (chrom, method, 'inherited CRU lost')
        new_mask = P.mask_intervals(additions.get(method, []), size)
        assert not (new_mask & protected).any(), (chrom, method, 'protected zero gap filled')
        assert not (new_mask & occupied).any(), (chrom, method, 'T6 overlap')
        allbeds[method] = rows
        P.save_bed(root / 'f03_tests' / method / f'{chrom}.CRU.bed', rows)
        P.save_bed(root / 'f03_tests' / method / f'{chrom}.added.bed', additions.get(method, []))
        summaries.append(summarise(method, chrom, rows, baseline, size))
        for r in rows:
            support_rows.append(dict(method=method, **r, **index.describe(r['start'], r['end'])))
        invariants.append(dict(method=method, chrom=chrom,
                               missing_t6_records=sum(r not in rows for r in baseline),
                               overlap_n=sum(x['end'] > y['start'] for x, y in zip(rows, rows[1:])),
                               protected_hard_gap_filled_bp=int((new_mask & protected).sum() * STEP),
                               added_overlap_t6_bp=int((new_mask & occupied).sum() * STEP),
                               bed_n=len(rows), reported_n=summaries[-1]['cru_n']))
    c = chrom
    density = [dict(chrom=c, position=i * STEP + min(STEP, size - i * STEP) / 2, raw=float(raw[i]),
                    density25=float(smooth[0, i]), density105=float(smooth[1, i]),
                    density225=float(smooth[2, i]), background=float(bg[i])) for i in range(len(raw))]
    P.write(root / 'f06_plot_data' / f'{c}.density.tsv.gz', density)
    P.write(root / 'f06_plot_data' / f'{c}.aglink.tsv.gz',
            [dict(chrom=c, left=e['left'], right=e['right'], acr_id=e['acr_id'], gene_id=e['gene_id']) for e in es],
            ['chrom', 'left', 'right', 'acr_id', 'gene_id'])
    P.write(root / 'f04_summary/decisions' / f'{c}.candidate_decisions.tsv.gz', decisions)
    P.write(root / 'f04_summary/summaries' / f'{c}.summary.tsv', summaries)
    P.write(root / 'f04_summary/invariants' / f'{c}.invariants.tsv', invariants)
    P.write(root / 'f04_summary/support' / f'{c}.CRU_support.tsv.gz', support_rows)
    if added:
        P.write(root / 'f04_summary/added' / f'{c}.added_support.tsv', added)
    artifacts = [root / 'f03_tests' / m / f'{c}.CRU.bed' for m in METHODS] + \
                [root / 'f03_tests' / m / f'{c}.added.bed' for m in METHODS] + \
                [root / 'f06_plot_data' / f'{c}.density.tsv.gz',
                 root / 'f06_plot_data' / f'{c}.aglink.tsv.gz',
                 root / 'f04_summary/decisions' / f'{c}.candidate_decisions.tsv.gz',
                 root / 'f04_summary/summaries' / f'{c}.summary.tsv',
                 root / 'f04_summary/invariants' / f'{c}.invariants.tsv',
                 root / 'f04_summary/support' / f'{c}.CRU_support.tsv.gz']
    if added:
        artifacts.append(root / 'f04_summary/added' / f'{c}.added_support.tsv')
    temp = Path(str(done) + '.tmp')
    temp.write_text(json.dumps({'signature': CTX['signature'],
                                'artifacts': {str(p): md5(p) for p in artifacts}}, indent=2) + '\n')
    temp.replace(done)
    print('T7 DONE', c, [(m, len(allbeds[m])) for m in METHODS], flush=True)
    return chrom, 'called'


def main(args):
    root = args.outdir.resolve()
    source = args.source_dir.resolve()
    t6 = args.t6_dir.resolve()
    opt = json.loads(args.settings.read_text())
    sizes = {s.split()[0]: int(s.split()[1]) for s in args.genome_sizes.read_text().splitlines() if s.strip()}
    chroms = args.chroms
    if len(set(chroms)) != len(chroms):
        raise ValueError('Duplicate chromosomes')
    for c in chroms:
        if c not in sizes:
            raise ValueError('Unknown chromosome: ' + c)
        if not (t6 / f'f03_tests/T6/{c}.CRU.bed').exists():
            raise ValueError('Saved T6 missing for ' + c + '; never substitute V11 (build T6 first)')
    if not opt.get('min_length_bp'):
        raise ValueError('Settings file is not a T7 island settings file')
    original = source / 'f03_results/f06_V11_300kb/02_CRU.bed'
    co_path = args.coacr_dir / '07_selected_coacr_pairs.tsv.gz'
    co_meta = json.loads((args.coacr_dir / '03_CACHE.DONE.json').read_text())
    if P.md5(co_path) != co_meta['cache_md5']:
        raise ValueError('Frozen coACR cache changed')
    validation = json.loads((source.parent / 'f04_block_enrichment/01_parameters.json').read_text())
    for path in [args.genome_sizes, args.block_bed, args.block_helper]:
        if validation['inputs'].get(str(path)) != P.md5(path):
            raise ValueError('Block validation input changed: ' + str(path))
    if (args.perm_n, args.seed, args.local_window_bp) != (validation['perm_n'], validation['seed'],
                                                          validation['local_window_bp']):
        raise ValueError('Keep original local-null settings (perm-n/seed/local-window-bp)')
    files = [original, source / 'f01_shared_input/01_count_only_aglink.annotated.bedpe', co_path,
             args.genome_sizes, args.block_bed, args.block_helper, args.gene_bed, args.acr_bed, args.settings,
             t6 / '01_parameters.json', t6 / 'f04_summary/01_chromosome_comparison.tsv',
             Path(S.__file__), Path(P.__file__)]
    for c in chroms:
        files += [t6 / f'f03_tests/T6/{c}.CRU.bed', source / f'f02_inputs/{c}.profiles.npz',
                  source / f'chromosomes/{c}/refined_hard_gaps.tsv']
    files += [p for p in HERE.iterdir() if p.is_file() and p.suffix in ('.py', '.R', '.json', '.md')]
    hashes = {str(p): P.md5(p) for p in files}
    saved_hashes = {r['path']: r['after_md5'] for r in S.table(t6 / '02_input_MD5.tsv')}
    for p in [source / 'f01_shared_input/01_count_only_aglink.annotated.bedpe', args.gene_bed, args.acr_bed]:
        if saved_hashes.get(str(p)) != hashes[str(p)]:
            raise ValueError('Data differ from fixed T6: ' + str(p))
    for c in chroms:
        p = source / f'f02_inputs/{c}.profiles.npz'
        if saved_hashes.get(str(p)) != hashes[str(p)]:
            raise ValueError('Density differs from fixed T6: ' + c)
    signature = {'inputs': hashes, 'chroms': chroms, 'settings': opt}
    if root.exists() and any(root.iterdir()) and not (root / '00_T7_signature.json').exists():
        raise ValueError('Unknown nonempty output directory')
    root.mkdir(parents=True, exist_ok=True)
    marker = root / '00_T7_signature.json'
    if marker.exists():
        if not args.resume or json.loads(marker.read_text()) != signature:
            raise ValueError('T7 inputs/settings differ; use a new directory')
    else:
        marker.write_text(json.dumps(signature, indent=2) + '\n')
    for d in ['chromosomes', 'f01_scripts', 'f03_tests', 'f04_summary', 'f06_plot_data',
              'f03_tests/T6', 'f03_tests/T7_A', 'f03_tests/T7_B',
              'f04_summary/decisions', 'f04_summary/summaries', 'f04_summary/invariants',
              'f04_summary/support', 'f04_summary/added']:
        (root / d).mkdir(parents=True, exist_ok=True)
    for p in HERE.iterdir():
        if p.is_file():
            shutil.copy2(p, root / 'f01_scripts' / p.name)
    (root / '01_parameters.json').write_text(json.dumps(
        dict(source_dir=str(source), t6_dir=str(t6), block_bed=str(args.block_bed),
             options={k: (str(v) if isinstance(v, Path) else v) for k, v in vars(args).items()},
             settings=opt, scope='Selected_chromosomes'), default=str, indent=2) + '\n')
    # T6 BEDs are copied byte-for-byte; T7 never edits the inherited baseline.
    for c in chroms:
        src = t6 / f'f03_tests/T6/{c}.CRU.bed'
        dst = root / f'f03_tests/T6/{c}.CRU.bed'
        shutil.copy2(src, dst)
        if md5(src) != md5(dst):
            raise ValueError('T6 copy differs: ' + c)
    edges = S.load_edges(source / 'f01_shared_input/01_count_only_aglink.annotated.bedpe', chroms)
    co = load_coacr(co_path, chroms)
    gene_meta = _load('m01_density_signal').metadata(args.gene_bed, chroms)
    acr_meta = _load('m01_density_signal').metadata(args.acr_bed, chroms)
    gp, ap = {}, {}
    for c in chroms:
        gp[c] = np.sort([b - 1 if strand == '-' else a for cc, a, b, strand in gene_meta.values() if cc == c])
        ap[c] = np.sort([(a + b) // 2 for cc, a, b, strand in acr_meta.values() if cc == c])
    CTX.update(root=root, source=source, t6_root=t6, sizes=sizes, edges=edges, co=co,
               gp=gp, ap=ap, opt=opt, resume=args.resume, signature=signature)
    with multiprocessing.get_context('fork').Pool(min(args.workers, len(chroms))) as pool:
        for chrom, state in pool.imap_unordered(worker, chroms):
            print('T7', state, chrom, flush=True)
    # ---- merge per-chromosome products ----------------------------------
    stats, invariants, added, support = [], [], [], []
    for c in chroms:
        stats += S.table(root / f'f04_summary/summaries/{c}.summary.tsv')
        invariants += S.table(root / f'f04_summary/invariants/{c}.invariants.tsv')
        support += S.table(root / f'f04_summary/support/{c}.CRU_support.tsv.gz')
        p = root / f'f04_summary/added/{c}.added_support.tsv'
        if p.exists():
            added += S.table(p)
    allbeds = {}
    for method in METHODS:
        rows = [r for c in chroms for r in P.read_bed(root / f'f03_tests/{method}/{c}.CRU.bed')]
        allbeds[method] = rows
        total = summarise(method, 'Selected_chromosomes', rows, allbeds['T6'], sum(sizes[c] for c in chroms))
        stats.append(total)
        P.save_bed(root / f'f03_tests/{method}/01_all_chromosomes.CRU.bed', rows)
    P.write(root / 'f04_summary/01_chromosome_comparison.tsv', stats)
    P.write(root / 'f04_summary/02_candidate_decisions.tsv',
            [r for c in chroms for r in S.table(root / f'f04_summary/decisions/{c}.candidate_decisions.tsv.gz')])
    P.write(root / 'f04_summary/03_added_CRU_support.tsv', added,
            None if added else ['chrom', 'start', 'end', 'method', 'accepted', 'reject_reason'])
    P.write(root / 'f04_summary/04_invariants.tsv', invariants)
    distributions = []
    for method in ['T7_A', 'T7_B']:
        rr = [r for r in added if method == 'T7_B' or int(r['accepted_A'])]
        for key in ['length_bp', 'annotated_gene_n', 'annotated_acr_n', 'independent_pair_n',
                    'raw_both_internal_n', 'raw_single_internal_n', 'raw_span_only_n',
                    'unique_both_internal_n', 'unique_single_internal_n', 'unique_span_only_n']:
            vv = [float(r[key]) for r in rr]
            distributions.append(dict(method=method, field=key, n=len(vv),
                                      median=float(np.median(vv)) if vv else np.nan,
                                      mean=float(np.mean(vv)) if vv else np.nan,
                                      min=min(vv) if vv else np.nan, max=max(vv) if vv else np.nan))
    P.write(root / 'f04_summary/05_support_distributions.tsv', distributions)
    P.write(root / 'f04_summary/06_all_CRU_aglink_support.tsv', support)
    report = ['# T7 low-background island rescue (genome-wide package)', '',
              'Scope: ' + ', '.join(chroms) + ' (%d chromosomes).' % len(chroms),
              'Saved T6 is read directly; no existing CRU is reconstructed or moved.',
              'Protected hard gaps are the unchanged T6 raw-zero >=100-kb regions.',
              'T7-B includes T7-A plus exactly-one-agLINK/coACR-supported additions; coACR never sets boundaries.',
              '', '| Method | CRUs | Added | Coverage | Median kb | >=5 Mb | >=10 Mb |',
              '| --- | ---: | ---: | ---: | ---: | ---: | ---: |']
    for r in stats:
        if r['chrom'] == 'Selected_chromosomes':
            report.append('| {method} | {cru_n} | {added_cru_n} | {pc:.2f}% | {med:.1f} | {ge5} | {ge10} |'.format(
                **r, pc=float(r['coverage_fraction']) * 100, med=float(r['median_bp']) / 1000,
                ge5=r['length_ge5mb_n'], ge10=r['length_ge10mb_n']))
    report += ['', '## Per-chromosome CRU counts', '', '| Chr | T6 | T7-A | T7-B | T7-A added | T7-B added |',
               '| --- | ---: | ---: | ---: | ---: | ---: |']
    for c in chroms:
        rr = {r['method']: r for r in stats if r['chrom'] == c}
        report.append('| {} | {} | {} | {} | {} | {} |'.format(
            c, rr['T6']['cru_n'], rr['T7_A']['cru_n'], rr['T7_B']['cru_n'],
            rr['T7_A']['added_cru_n'], rr['T7_B']['added_cru_n']))
    report += ['', 'New islands (T7-B route, of which T7-A is a subset): %d.' % len(added),
               'No final-version replacement or parameter selection based on block enrichment was performed.']
    (root / '04_REPORT.md').write_text('\n'.join(report) + '\n')
    if any(P.md5(p) != h for p, h in hashes.items()):
        raise ValueError('Input changed during run')
    P.write(root / '02_input_MD5.tsv',
            [dict(path=p, before_md5=h, after_md5=P.md5(p)) for p, h in hashes.items()])
    (root / '03_PREPARATION.DONE').write_text('T7 completed; fixed T6 unchanged.\n')
    for r in stats:
        if r['chrom'] == 'Selected_chromosomes':
            print('SUMMARY', r, flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for key in ['source-dir', 't6-dir', 'coacr-dir', 'gene-bed', 'acr-bed', 'genome-sizes', 'block-bed',
                'block-helper', 'settings', 'outdir']:
        p.add_argument('--' + key, type=Path, required=True)
    p.add_argument('--chroms', nargs='+', required=True)
    p.add_argument('--workers', type=int, default=7)
    p.add_argument('--perm-n', type=int, default=1000)
    p.add_argument('--seed', type=int, default=20260805)
    p.add_argument('--local-window-bp', type=int, default=10000000)
    p.add_argument('--resume', action='store_true')
    args = p.parse_args()
    if args.workers < 1:
        raise SystemExit('--workers must be positive')
    main(args)
