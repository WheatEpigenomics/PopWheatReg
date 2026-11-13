# m13_low_density_recall: T3 low-density recall
# Frozen rule module, copied unchanged from the published calling chain
# (source: f18_low_density_recall/01_low_density_recall.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
from ._load import load as _load
"""Read frozen V11 signals; diagnose exclusions and append isolated recall tests.

No existing CRU is edited. Blocks are read only after all candidate decisions.
"""
import argparse
import csv
import gzip
import hashlib
import importlib.util
import json
from pathlib import Path
import random
import runpy
import shutil
from types import SimpleNamespace

import numpy as np
from scipy.ndimage import uniform_filter1d


HERE = Path(__file__).resolve().parent
STEP = 5000


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def md5(path):
    h = hashlib.md5()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def write(path, rows, fields=None):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = fields or list(dict.fromkeys(k for r in rows for k in r)) or ['chrom']
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'wt') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n')
        w.writeheader()
        w.writerows(rows)


def read_bed(path):
    rows = []
    with Path(path).open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            c, a, b, ident = line.split()[:4]
            a, b = int(a), int(b)
            if not 0 <= a < b:
                raise ValueError('Invalid BED coordinates')
            rows.append(dict(chrom=c, start=a, end=b, cru_id=ident))
    return rows


def save_bed(path, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as f:
        for r in rows:
            f.write('{chrom}\t{start}\t{end}\t{cru_id}\n'.format(**r))


def runs(mask):
    d = np.diff(np.r_[False, mask, False].astype(int))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)))


def mask_intervals(rows, size):
    mask = np.zeros((size + STEP - 1) // STEP, dtype=bool)
    for r in rows:
        if not 0 <= r['start'] < r['end'] <= size:
            raise ValueError('Interval outside chromosome bounds')
        # V11 cuts can be bin centers. Protect every intersected bin; never
        # round the original BED coordinates or recall its partial edge bins.
        mask[r['start']//STEP:(r['end']+STEP-1)//STEP] = True
    return mask


def close_supported(mask, available, max_hole_bp):
    out = mask.copy()
    for a, b in runs(~mask):
        if a > 0 and b < len(mask) and (b-a)*STEP <= max_hole_bp and available[a:b].all():
            out[a:b] = True
    return out & available


def support(edges, start, end):
    inside = [e for e in edges if start <= e['left'] and e['right'] < end]
    unique = {(e['acr_id'], e['gene_id']): e for e in inside}
    ii = list(unique.values())
    acr = {e['acr_id'] for e in ii}
    gene = {e['gene_id'] for e in ii}
    positions = {e['p1'] for e in ii} | {e['p2'] for e in ii}
    spans = sorted((e['left'], e['right']+1) for e in ii)
    union = 0
    lo = hi = None
    for a, b in spans:
        if hi is None or a > hi:
            if hi is not None:
                union += hi-lo
            lo, hi = a, b
        else:
            hi = max(hi, b)
    if hi is not None:
        union += hi-lo
    return dict(internal_aglink_n=len(inside), independent_pair_n=len(ii), acr_n=len(acr),
                gene_n=len(gene), acr_position_n=len({e['p1'] for e in ii}),
                tss_position_n=len({e['p2'] for e in ii}),
                endpoint_extent_bp=max(positions)-min(positions) if positions else 0,
                internal_span_covered_fraction=union/(end-start),
                independent_pairs_per_100kb=len(ii)*100000/(end-start),
                spanning_entire_region_n=sum(e['left'] < start and e['right'] >= end for e in edges),
                touching_aglink_n=sum(start <= e['p1'] < end or start <= e['p2'] < end for e in edges))


def describe(raw, smooth, bg, endpoint, a, b, size):
    end = min(size, b*STEP)
    lengths = np.minimum(STEP, size-np.arange(a, b)*STEP)
    zero_runs = runs(raw[a:b] <= 1e-9)
    return dict(start=a*STEP, end=end, length_bp=end-a*STEP,
                raw_mean=float(np.average(raw[a:b], weights=lengths)),
                raw_median=float(np.median(raw[a:b])), raw_p90=float(np.quantile(raw[a:b], .9)),
                density25_median=float(np.median(smooth[0, a:b])),
                density105_median=float(np.median(smooth[1, a:b])),
                density225_median=float(np.median(smooth[2, a:b])),
                background_median=float(np.median(bg[a:b])),
                background_near_zero_fraction=float(np.mean(bg[a:b] <= .1)),
                positive_fraction=float(np.average(raw[a:b] > 1e-9, weights=lengths)),
                max_zero_run_bp=max((int(lengths[x:y].sum()) for x, y in zero_runs), default=0),
                endpoint225_median=float(np.median(endpoint[a:b])))


def proposal_mask(raw, smooth, bg, available, opt, low_background=False):
    # Background is regularized as bg + sigma*sqrt(bg+1), never divided by zero.
    threshold = np.maximum(opt['absolute_floor'], bg + opt['background_sigma']*np.sqrt(bg+1))
    good = (raw >= opt['absolute_floor']-1e-9) & (np.sum(smooth >= threshold-1e-9, axis=0) >= 2)
    if low_background:
        good &= bg <= opt['near_zero_background']
    return close_supported(good & available, available, opt['close_hole_bp'])


def reasons_for(row, opt, supported_fraction):
    fail = []
    for name, key in [('length_bp', 'min_territory_bp'), ('independent_pair_n', 'min_pairs'),
                      ('acr_n', 'min_acrs'), ('gene_n', 'min_genes'),
                      ('acr_position_n', 'min_acrs'), ('tss_position_n', 'min_genes'),
                      ('endpoint_extent_bp', 'min_endpoint_extent_bp'),
                      ('internal_span_covered_fraction', 'min_internal_span_fraction'),
                      ('independent_pairs_per_100kb', 'min_pairs_per_100kb')]:
        if row[name] < opt[key]:
            fail.append(name)
    if supported_fraction < opt['min_supported_fraction']:
        fail.append('sustained_density_support')
    if row['max_zero_run_bp'] > opt['max_zero_run_bp']:
        fail.append('wide_zero_gap')
    return fail


def disjoint(rows):
    rows = sorted(rows, key=lambda r: (r['start'], r['end']))
    return all(a['end'] <= b['start'] for a, b in zip(rows, rows[1:]))


def boundary_groups(old, added, territories, size):
    def points(rows):
        return {p for r in rows for p in [r['start'], r['end']] if 0 < p < size}
    allrows = old + added
    starts, ends = {r['start'] for r in allrows}, {r['end'] for r in allrows}
    total = points(allrows)
    internal = starts & ends & total
    added_all = points(added)-points(old)
    new_outer = (points(territories)-points(old)) & added_all
    return dict(outer=total-internal, internal=internal, all_final=total,
                added_outer=new_outer, added_internal=added_all-new_outer, added_all=added_all)


def validate_blocks(args, summary_groups, sizes, out):
    helper = runpy.run_path(str(args.block_helper))
    blocks = helper['load_block_bed'](str(args.block_bed), args.chroms)
    index = helper['block_test_index'](blocks, args.chroms)
    def nearest(pp):
        values = []
        for c in args.chroms:
            pos = np.array([p['pos'] for p in pp if p['chrom'] == c], dtype=np.int64)
            bb = np.array(index[c]['boundaries'], dtype=np.int64)
            if not len(pos):
                continue
            if not len(bb):
                raise ValueError('Missing blocks: '+c)
            j = np.searchsorted(bb, pos)
            dd = np.minimum(abs(pos-bb[np.maximum(0, j-1)]), abs(pos-bb[np.minimum(len(bb)-1, j)]))
            values.extend(zip([c]*len(pos), pos.tolist(), dd.tolist()))
        return values
    summaries, detail, replicates = [], [], []
    for (method, group), bychrom in summary_groups.items():
        for scope in args.chroms+['Selected_chromosomes']:
            pp = [dict(chrom=c, pos=int(p)) for c in args.chroms for p in sorted(bychrom.get(c, []))
                  if scope == 'Selected_chromosomes' or scope == c]
            n = len(pp)
            dd = nearest(pp)
            if scope == 'Selected_chromosomes':
                detail.extend(dict(method=method, boundary_group=group, chrom=c, position=p, distance_bp=d) for c, p, d in dd)
            observed = np.array([sum(d <= k for _, _, d in dd) for k in args.block_cutoffs])
            null = np.zeros((args.perm_n, len(args.block_cutoffs)), dtype=int)
            rng = random.Random(args.seed)
            for rep in range(args.perm_n if n else 0):
                shuffled = helper['local_random_points'](pp, sizes, args.local_window_bp, rng)
                assert len(shuffled) == n, 'Local null must preserve query count'
                distances = nearest(shuffled)
                for j, k in enumerate(args.block_cutoffs):
                    nn = sum(d <= k for _, _, d in distances)
                    null[rep, j] = nn
                    replicates.append(dict(method=method, boundary_group=group, scope=scope, replicate=rep+1,
                                           cutoff_bp=k, boundary_n=n, random_match_n=nn, random_fraction=nn/n))
            for j, k in enumerate(args.block_cutoffs):
                mean = float(null[:, j].mean()) if n else np.nan
                low, high = np.quantile(null[:, j], [.025, .975]) if n else (np.nan, np.nan)
                summaries.append(dict(method=method, boundary_group=group, scope=scope, cutoff_bp=k, boundary_n=n,
                    observed_n=int(observed[j]), observed_fraction=observed[j]/n if n else np.nan,
                    random_match_n_mean=mean, random_match_n_p025=low, random_match_n_p975=high,
                    random_fraction_mean=mean/n if n else np.nan, random_fraction_p025=low/n if n else np.nan,
                    random_fraction_p975=high/n if n else np.nan, enrichment_fold=observed[j]/mean if mean > 0 else np.nan,
                    empirical_p=(1+np.count_nonzero(null[:, j] >= observed[j]))/(args.perm_n+1) if n else np.nan,
                    permutation_n=args.perm_n))
    write(out/'01_block_enrichment.tsv', summaries)
    write(out/'02_random_replicates.tsv.gz', replicates)
    write(out/'03_observed_distances.tsv', detail)


def main(args):
    if args.perm_n < 1 or args.min_daughter_bp != 300000 or args.window_bp < STEP:
        raise ValueError('Positive permutation count and frozen 300-kb daughter limit required')
    source = args.source_dir.resolve()
    rules_path = source/'f04_statistics/07_effective_rules.json'
    rules = json.loads(rules_path.read_text())
    config = json.loads(args.settings.read_text())
    for stage in ['T1', 'T2', 'T3_rescue']:
        opt = config[stage]
        if opt['min_territory_bp'] < 300000 or opt['absolute_floor'] <= 0 or opt['min_pairs'] < 4:
            raise ValueError('Invalid recall length/density/independent-pair protection')
        if min(opt['min_acrs'], opt['min_genes']) < 2 or opt['max_zero_run_bp'] > 10000:
            raise ValueError('Insufficient node diversity / wide zero gap bridging requested')
    saved_signature = json.loads((source/'00_call_signature.json').read_text())
    for path, expected in saved_signature['inputs'].items():
        if path.endswith('.py') and md5(path) != expected:
            raise ValueError('Frozen calling source changed: '+path)
    engine = _load('m11_engine')
    M, H, C, D = engine.M, engine.H, engine.H.C, engine.H.D
    old_path = source/'f03_results/f06_V11_300kb/02_CRU.bed'
    old_all = read_bed(old_path)
    if len(old_all) != args.expected_cru_n:
        raise ValueError('Frozen population differs from expected CRU count')
    sizes = {}
    with args.genome_sizes.open() as f:
        for line in f:
            if line.strip():
                c, n = line.split()[:2]
                sizes[c] = int(n)
    if len(args.chroms) != len(set(args.chroms)) or set(args.chroms)-set(sizes):
        raise ValueError('Invalid chromosomes')
    inputs = [old_path, rules_path, args.settings, args.gene_bed, args.acr_bed, args.genome_sizes,
              args.aglink, args.block_bed, args.block_helper, source/'f01_shared_input/01_count_only_aglink.annotated.bedpe']
    inputs += [source/f'f02_inputs/{c}.profiles.npz' for c in args.chroms]
    inputs += [source/f'chromosomes/{c}/refined_hard_gaps.tsv' for c in args.chroms]
    inputs += list(HERE.glob('*.py'))+list(HERE.glob('*.R'))+list(HERE.glob('*.md'))
    hashes = {str(p): md5(p) for p in inputs}
    validation_path = source.parent/'f04_block_enrichment/01_parameters.json'
    if not validation_path.is_file():
        raise ValueError('Original block-validation parameter file is required')
    validation = json.loads(validation_path.read_text())
    for p in [args.block_bed, args.block_helper, args.genome_sizes]:
        if validation['inputs'].get(str(p)) != md5(p):
            raise ValueError('Block-validation input differs from original V11: '+str(p))
    if not args.skip_block and (args.perm_n, args.seed, args.local_window_bp) != (
            validation['perm_n'], validation['seed'], validation['local_window_bp']):
        raise ValueError('Keep the original permutation count, seed, and local window')
    hashes[str(validation_path)] = md5(validation_path)
    # Confirm supplied raw input is the same one used by the frozen run.
    for supplied in [args.aglink, args.gene_bed, args.acr_bed, args.genome_sizes]:
        wanted = [v for p, v in saved_signature['inputs'].items() if Path(p).resolve() == supplied.resolve()]
        if len(wanted) != 1 or wanted[0] != hashes[str(supplied)]:
            raise ValueError('Input does not match original V11: '+str(supplied))
    root = args.outdir.resolve()
    if root == source or source in root.parents:
        raise ValueError('Use an independent sibling result directory')
    if root.exists() and any(p.name != 'f01_scripts' for p in root.iterdir()):
        previous = root/'01_parameters.json'
        if not args.resume or not previous.is_file():
            raise ValueError('Output exists; use --resume only for a preparation interrupted before test BEDs')
        old_settings = json.loads(previous.read_text())
        current = json.loads(json.dumps(vars(args), default=str))
        if any(old_settings.get(k) != v for k, v in current.items() if k != 'resume'):
            raise ValueError('Interrupted run arguments differ; use a fresh numbered directory')
        if old_settings.get('settings_content') != config or old_settings.get('frozen_rules') != rules:
            raise ValueError('Interrupted run rules changed; use a fresh numbered directory')
        if (root/'03_PREPARATION.DONE').exists() or any((root/'f03_tests').rglob('*.bed')):
            raise ValueError('Existing test BEDs/results will not be overwritten; use a fresh numbered directory')
    for name in ['f01_scripts', 'f02_diagnostics', 'f03_tests', 'f04_summary', 'f05_block_validation', 'f06_plot_data', 'f07_figures']:
        (root/name).mkdir(parents=True, exist_ok=True)
    for p in HERE.iterdir():
        if p.is_file():
            shutil.copy2(p, root/'f01_scripts'/p.name)
    (root/'01_parameters.json').write_text(json.dumps({**vars(args), 'settings_content': config,
        'population': 'only requested chromosomes; original global CRUs kept intact',
        'frozen_rules': rules}, default=str, indent=2)+'\n')
    (root/'f01_scripts/05_frozen_calling_rules.json').write_text(json.dumps(rules, indent=2)+'\n')
    edges = D.load_edges(source/'f01_shared_input/01_count_only_aglink.annotated.bedpe', args.chroms)
    stats, candidates, bin_diagnostics, exclusion_rows, signal_islands, all_added, windows = [], [], [], [], [], [], []
    block_groups = {}
    for chrom in args.chroms:
        print('Preparing diagnostic/sensitivity:', chrom, flush=True)
        size, es = sizes[chrom], edges[chrom]
        old = sorted([r for r in old_all if r['chrom'] == chrom], key=lambda r: r['start'])
        if not disjoint(old):
            raise ValueError('Frozen CRUs overlap')
        with np.load(source/f'f02_inputs/{chrom}.profiles.npz') as z:
            raw, smooth, bg = z['full_raw'].copy(), z['full_smooth'].copy(), z['full_background'].copy()
        n = len(raw)
        if n != (size+STEP-1)//STEP or smooth.shape != (3, n):
            raise ValueError('Unexpected saved profile shape')
        endpoints = np.zeros(n)
        for e in es:
            endpoints[e['p1']//STEP] += 1
            endpoints[e['p2']//STEP] += 1
        ep = uniform_filter1d(endpoints, 45, mode='constant')*45
        sig = D.territory_signals(raw, smooth, ep, SimpleNamespace(**rules['primary']))
        if not np.allclose(bg, sig['background'], rtol=1e-9, atol=1e-9):
            raise ValueError('Saved/recomputed local background mismatch')
        hs, gapmask, gaps, cores = C.identify_gaps(chrom, size, sig, es, SimpleNamespace(**rules['hard_gap']), D)
        _, actions = C.base_components(chrom, size, gapmask, sig, es, SimpleNamespace(**rules['hard_gap']), D)
        refined, events = H.E.refine_edges(chrom, sig, es, gaps, cores, SimpleNamespace(**rules['A_edge']))
        with (source/f'chromosomes/{chrom}/refined_hard_gaps.tsv').open() as f:
            saved = [(int(r['start']), int(r['end'])) for r in csv.DictReader(f, delimiter='\t')]
        if saved != [(r['start'], r['end']) for r in refined]:
            raise ValueError('Hard-gap replay differs from saved final gaps: '+chrom)
        diag = root/'f02_diagnostics'/chrom
        write(diag/'01_hard_gap_cores.tsv', cores)
        write(diag/'02_original_hard_gaps.tsv', gaps)
        write(diag/'03_base_component_decisions.tsv', actions)
        write(diag/'04_A_edge_decisions.tsv', events)
        write(diag/'05_verified_refined_gaps.tsv', refined)
        occupied = mask_intervals(old, size)
        available = ~occupied
        refined_mask = mask_intervals(refined, size)
        failed_length = np.zeros(n, bool)
        failed_links = np.zeros(n, bool)
        for row in actions:
            a, b = row['start']//STEP, (row['end']+STEP-1)//STEP
            if row['end']-row['start'] < rules['hard_gap']['min_cru_bp']:
                failed_length[a:b] = True
            if row['internal_n'] < rules['hard_gap']['min_internal']:
                failed_links[a:b] = True
        flags = dict(hard_absolute=hs['hard_absolute'], hard_relative=hs['hard_relative'],
                     original_hard_gap=gapmask, refined_hard_gap=refined_mask,
                     base_too_short=failed_length, base_insufficient_internal=failed_links,
                     legacy_seed_fail=~sig['seed'], raw_zero=raw <= 1e-9)
        for name, mask in flags.items():
            lengths = np.minimum(STEP, size-np.arange(n)*STEP)
            bin_diagnostics.append(dict(chrom=chrom, criterion=name, criterion_active_for_territory=int(name != 'legacy_seed_fail'),
                total_flagged_bp=int(lengths[mask].sum()), flagged_uncovered_bp=int(lengths[mask & available].sum()),
                flagged_uncovered_positive_bp=int(lengths[mask & available & (raw > 1e-9)].sum())))
        for a, b in runs(available):
            row = dict(chrom=chrom, **describe(raw, smooth, bg, ep, a, b, size))
            row.update(support(es, row['start'], row['end']))
            if row['positive_fraction'] == 0:
                cause = 'no_span_signal'
            else:
                cause = ','.join(name for name in ['hard_absolute', 'hard_relative', 'refined_hard_gap', 'base_too_short', 'base_insufficient_internal'] if flags[name][a:b].any()) or 'zero_trim_or_unassigned_margin'
            row['original_exclusion_causes'] = cause
            row.update({name+'_bin_fraction': float(mask[a:b].mean()) for name, mask in flags.items()})
            exclusion_rows.append(row)
        for a, b in runs(available & (raw > 1e-9)):
            row = dict(chrom=chrom, **describe(raw, smooth, bg, ep, a, b, size))
            row.update(support(es, row['start'], row['end']))
            row.update({name+'_bin_fraction': float(mask[a:b].mean()) for name, mask in flags.items()})
            row['signal_description'] = ('span_only_without_internal_pair' if not row['independent_pair_n']
                                         else 'positive_span_with_internal_links')
            signal_islands.append(row)
        selected, territory_by_method, cut_by_method = {}, {'T0': []}, {'T0': []}
        for method in ['T1', 'T2', 'T3']:
            stages = ['T1'] if method == 'T1' else ['T2'] if method == 'T2' else ['T2', 'T3_rescue']
            territories, cuts, additions = [], [], []
            accepted_mask = np.zeros(n, bool)
            for stage in stages:
                opt = config[stage]
                remaining = available & ~accepted_mask
                proposed = proposal_mask(raw, smooth, bg, remaining, opt, stage == 'T3_rescue')
                for a, b in runs(proposed):
                    row = dict(method=method, stage=stage, chrom=chrom, **describe(raw, smooth, bg, ep, a, b, size))
                    row.update(support(es, row['start'], row['end']))
                    thresh = np.maximum(opt['absolute_floor'], bg[a:b]+opt['background_sigma']*np.sqrt(bg[a:b]+1))
                    supported = (raw[a:b] >= opt['absolute_floor']-1e-9) & (np.sum(smooth[:, a:b] >= thresh-1e-9, axis=0) >= 2)
                    row['supported_bin_fraction'] = float(supported.mean())
                    fail = reasons_for(row, opt, row['supported_bin_fraction'])
                    row.update(accepted=int(not fail), reject_reason=','.join(fail) or 'supported_new_territory',
                        original_hard_gap_fraction=float(refined_mask[a:b].mean()),
                        old_hard_absolute_fraction=float(hs['hard_absolute'][a:b].mean()),
                        old_hard_relative_fraction=float(hs['hard_relative'][a:b].mean()))
                    candidates.append(row)
                    if fail:
                        continue
                    parent = dict(chrom=chrom, start=row['start'], end=row['end'],
                                  cru_id=f'{chrom}_{method}_recalled_territory_{len(territories)+1}')
                    territories.append(parent)
                    accepted_mask[a:b] = True
                    if args.split_new == 'frozen':
                        acr, genes = H.annotation_positions(chrom, args.gene_bed, args.acr_bed)
                        local = [e for e in es if e['left'] < parent['end'] and e['right'] >= parent['start']]
                        rr, bb, cc, _ = M.solve(parent, sig, local, size, SimpleNamespace(**rules['primary']),
                            M.options(rules['shape_settings'], args.min_daughter_bp), acr, genes, method, rules['max_rounds'])
                        write(diag/f'new_{method}_{len(territories):04d}_split_candidates.tsv', cc)
                        cuts.extend(bb)
                    else:
                        rr = [parent.copy()]
                    rr.sort(key=lambda r: r['start'])
                    assert rr[0]['start'] == parent['start'] and rr[-1]['end'] == parent['end']
                    assert all(x['end'] == y['start'] for x, y in zip(rr, rr[1:]))
                    assert sum(r['end']-r['start'] for r in rr) == parent['end']-parent['start']
                    for r in rr:
                        r.update(cru_id=f'{chrom}_{method}_new_CRU_{len(additions)+1}', territory_id=parent['cru_id'])
                        additions.append(r)
            assert disjoint(old+additions)
            selected[method] = additions
            territory_by_method[method], cut_by_method[method] = territories, cuts
            all_added.extend(dict(method=method, **r) for r in additions)
        selected['T0'] = []
        for method in ['T0', 'T1', 'T2', 'T3']:
            additions = selected[method]
            combined = sorted(old+additions, key=lambda r: r['start'])
            assert all(r in combined for r in old)
            target = root/'f03_tests'/method
            save_bed(target/f'{chrom}.CRU.bed', combined)
            save_bed(target/f'{chrom}.added_CRU.bed', additions)
            save_bed(target/f'{chrom}.added_territory.bed', territory_by_method[method])
            lengths = [r['end']-r['start'] for r in combined]
            support_new = sum(support(es, r['start'], r['end'])['independent_pair_n'] for r in territory_by_method[method])
            stats.append(dict(method=method, chrom=chrom, chromosome_bp=size, original_cru_n=len(old),
                cru_n=len(combined), added_cru_n=len(additions), new_territory_n=len(territory_by_method[method]),
                new_internal_cut_n=len(cut_by_method[method]), covered_bp=sum(lengths), coverage_fraction=sum(lengths)/size,
                mean_bp=float(np.mean(lengths)) if lengths else np.nan, median_bp=float(np.median(lengths)) if lengths else np.nan,
                added_coverage_bp=sum(r['end']-r['start'] for r in additions), new_territory_internal_pairs=support_new))
            for group, points in boundary_groups(old, additions, territory_by_method[method], size).items():
                block_groups.setdefault((method, group), {})[chrom] = points
        density = [dict(chrom=chrom, position=i*STEP+min(STEP, size-i*STEP)/2, raw=float(raw[i]), density25=float(smooth[0, i]),
                   density105=float(smooth[1, i]), density225=float(smooth[2, i]), background=float(bg[i]),
                   endpoint225=float(ep[i]), original_uncovered=int(available[i]),
                   hard_absolute=int(hs['hard_absolute'][i]), hard_relative=int(hs['hard_relative'][i])) for i in range(n)]
        write(root/f'f06_plot_data/{chrom}.density.tsv.gz', density)
        write(root/f'f06_plot_data/{chrom}.aglink.tsv.gz', [dict(chrom=chrom, left=e['left'], right=e['right'], acr_id=e['acr_id'], gene_id=e['gene_id']) for e in es], ['chrom', 'left', 'right', 'acr_id', 'gene_id'])
        for start in range(0, size, args.window_bp):
            windows.append(dict(chrom=chrom, start=start, end=min(size, start+args.window_bp)))
    # Pooled rows are computed from the existing per-chromosome products, not averaged fractions.
    for method in ['T0', 'T1', 'T2', 'T3']:
        ss = [r for r in stats if r['method'] == method]
        lengths = [r['end']-r['start'] for r in old_all if r['chrom'] in args.chroms]
        lengths += [r['end']-r['start'] for r in all_added if r['method'] == method]
        totals = {key: sum(r[key] for r in ss) for key in ['chromosome_bp', 'original_cru_n', 'cru_n', 'added_cru_n', 'new_territory_n', 'new_internal_cut_n', 'covered_bp', 'added_coverage_bp', 'new_territory_internal_pairs']}
        stats.append(dict(method=method, chrom='Selected_chromosomes', **totals,
                          coverage_fraction=totals['covered_bp']/totals['chromosome_bp'],
                          mean_bp=float(np.mean(lengths)), median_bp=float(np.median(lengths))))
    write(root/'f04_summary/01_chromosome_comparison.tsv', stats)
    write(root/'f04_summary/02_candidate_decisions.tsv', candidates)
    write(root/'f04_summary/03_original_excluded_regions.tsv', exclusion_rows)
    write(root/'f04_summary/04_criterion_bp_audit.tsv', bin_diagnostics)
    write(root/'f04_summary/05_added_CRU_detail.tsv', all_added)
    write(root/'f04_summary/06_excluded_positive_signal_islands.tsv', signal_islands)
    write(root/'f06_plot_data/01_windows.tsv', windows)
    for c in args.chroms:
        cc = [r for r in candidates if r['chrom'] == c]
        chosen = []
        for status in [1, 0]:
            pool = [r for r in cc if r['accepted'] == status]
            # Distinct locations, deterministic seed; avoid repeated T2/T3 candidate records.
            distinct = {(r['start'], r['end']): r for r in pool}
            pool = list(distinct.values())
            random.Random(args.seed).shuffle(pool)
            chosen.extend(pool[:args.examples_per_status])
        write(root/f'f06_plot_data/{c}.examples.tsv', chosen, list(candidates[0]) if candidates else ['chrom', 'start', 'end', 'accepted'])
    if not args.skip_block:
        validate_blocks(args, block_groups, sizes, root/'f05_block_validation')
    if any(md5(p) != h for p, h in hashes.items()):
        raise ValueError('An input changed during execution')
    write(root/'02_input_MD5.tsv', [dict(path=p, before_md5=h, after_md5=md5(p)) for p, h in hashes.items()])
    report = [
        '# Low-density recall sensitivity', '',
        'This is a diagnostic experiment, not a replacement V11 result.',
        'Frozen global input: {} CRUs; evaluated chromosomes: {}.'.format(len(old_all), ', '.join(args.chroms)),
        'All inherited selected-chromosome CRU records are retained verbatim in T0/T1/T2/T3.',
        'Background/absolute-core flags can overlap; their bp must not be added as exclusive causes.',
        'Legacy seed failure is NOT an active territory rejection criterion in this calling path.', '',
        '| Chromosome | Test | CRUs | New territories | New internal cuts | Coverage |',
        '| --- | --- | ---: | ---: | ---: | ---: |']
    for row in stats:
        report.append('| {chrom} | {method} | {cru_n} | {new_territory_n} | {new_internal_cut_n} | {pct:.2f}% |'.format(
            **row, pct=100*row['coverage_fraction']))
    report += ['', '## Candidate failures', '',
               'Failure counts below are nonexclusive and refer to proposed islands, not all chromosome bins.']
    for method in ['T1', 'T2', 'T3']:
        cc = [r for r in candidates if r['method'] == method]
        reasons = sorted({s for r in cc if not r['accepted'] for s in r['reject_reason'].split(',')})
        report.append('{}: {} candidate records, {} accepted territories.'.format(method, len(cc), sum(r['accepted'] for r in cc)))
        for reason in reasons:
            report.append('- {}: {}'.format(reason, sum(reason in r['reject_reason'].split(',') for r in cc if not r['accepted'])))
    report += ['', 'Block validation (if enabled) is independent and uses the original local null.',
               'Read f05_block_validation/01_block_enrichment.tsv for fold, P, numerators and denominators.',
               'No score or filtering of the fixed V11 CRUs has been applied.']
    (root/'04_REPORT.md').write_text('\n'.join(report)+'\n')
    (root/'03_PREPARATION.DONE').write_text('Frozen CRUs unchanged; only uncovered regions tested.\n')
    print('Outputs:', root, flush=True)
    for row in stats:
        if row['chrom'] == 'Selected_chromosomes':
            print(row['method'], 'CRUs', row['cru_n'], 'new territories', row['new_territory_n'], 'new internal cuts', row['new_internal_cut_n'], flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for k in ['source-dir', 'settings', 'aglink', 'gene-bed', 'acr-bed', 'genome-sizes', 'block-bed', 'block-helper', 'outdir']:
        p.add_argument('--'+k, type=Path, required=True)
    p.add_argument('--chroms', nargs='+', default=['4B', '4D', '5D', '7D', '3A', '5A', '1B', '2D'])
    p.add_argument('--expected-cru-n', type=int, default=3159)
    p.add_argument('--min-daughter-bp', type=int, default=300000)
    p.add_argument('--split-new', choices=['none', 'frozen'], default='frozen')
    p.add_argument('--window-bp', type=int, default=5000000)
    p.add_argument('--perm-n', type=int, default=1000)
    p.add_argument('--seed', type=int, default=20260805)
    p.add_argument('--local-window-bp', type=int, default=10000000)
    p.add_argument('--block-cutoffs', nargs='+', type=int, default=[25000, 50000, 100000, 200000])
    p.add_argument('--skip-block', action='store_true')
    p.add_argument('--resume', action='store_true', help='Retry matching interrupted preparation only before any test BED exists')
    p.add_argument('--examples-per-status', type=int, default=3)
    main(p.parse_args())
