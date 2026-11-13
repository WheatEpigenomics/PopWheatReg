# m14_support_recall: T4/T5/T6 support recall
# Frozen rule module, copied unchanged from the published calling chain
# (source: f19_low_density_support/01_test_support_recall.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
from ._load import load as _load
"""Append T4/T5/T6 recalls to saved T3; never reconstruct or edit inherited CRUs."""
import argparse
from collections import Counter, defaultdict
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import pickle
import random
import shutil
from types import SimpleNamespace

import numpy as np
from scipy.ndimage import uniform_filter1d


HERE = Path(__file__).resolve().parent
from . import m13_low_density_recall as P
STEP = P.STEP
METHODS = ['T3', 'T4', 'T5', 'T6']


def table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def load_edges(path, chroms):
    out = {c: [] for c in chroms}
    with path.open() as f:
        for r in csv.DictReader(f, delimiter='\t'):
            if r['chrom'] not in out:
                continue
            for k in ['p1', 'p2', 'left', 'right']:
                r[k] = int(r[k])
            if float(r['weight']) != 1:
                raise ValueError('Frozen input is not equal-weight agLINK')
            out[r['chrom']].append(r)
    return out


class EdgeIndex:
    """Endpoint ownership is half-open BED membership, not span intersection."""
    def __init__(self, edges):
        self.edges = edges
        self.p1 = np.array([r['p1'] for r in edges], dtype=np.int64)
        self.p2 = np.array([r['p2'] for r in edges], dtype=np.int64)
        seen = {}
        for e in edges:
            key = (e['acr_id'], e['gene_id'])
            if key in seen and (seen[key]['p1'], seen[key]['p2']) != (e['p1'], e['p2']):
                raise ValueError('Same ACR/gene pair has conflicting endpoint coordinates')
            seen[key] = e
        self.unique = list(seen.values())
        self.u1 = np.array([r['p1'] for r in self.unique], dtype=np.int64)
        self.u2 = np.array([r['p2'] for r in self.unique], dtype=np.int64)

    def describe(self, start, end):
        result = {}
        selected = []
        for prefix, edges, p1, p2 in [('raw', self.edges, self.p1, self.p2),
                                      ('unique', self.unique, self.u1, self.u2)]:
            a, b = (p1 >= start) & (p1 < end), (p2 >= start) & (p2 < end)
            inside, single = a & b, a ^ b
            bridge = ~a & ~b & (np.minimum(p1, p2) < start) & (np.maximum(p1, p2) >= end)
            outside = ~inside & ~single & ~bridge
            for label, mask in [('both_internal', inside), ('single_internal', single),
                                ('span_only', bridge), ('outside', outside)]:
                result[prefix+'_'+label+'_n'] = int(mask.sum())
            assert sum(result[prefix+'_'+label+'_n'] for label in ['both_internal', 'single_internal', 'span_only', 'outside']) == len(edges)
            if prefix == 'unique':
                selected = [edges[i] for i in np.flatnonzero(inside)]
        positions = sorted({e[k] for e in selected for k in ['p1', 'p2']})
        acrs = {e['acr_id'] for e in selected}
        genes = {e['gene_id'] for e in selected}
        spans = sorted((e['left'], e['right']+1) for e in selected)
        covered, lo, hi = 0, None, None
        for a, b in spans:
            if hi is None or a > hi:
                if hi is not None:
                    covered += hi-lo
                lo, hi = a, b
            else:
                hi = max(hi, b)
        if hi is not None:
            covered += hi-lo
        # Bipartite connected components: shared nodes are real ID identities.
        parent = {}
        def find(x):
            parent.setdefault(x, x)
            while parent[x] != x:
                parent[x] = parent[parent[x]]
                x = parent[x]
            return x
        for e in selected:
            a, b = ('ACR', e['acr_id']), ('gene', e['gene_id'])
            parent[find(a)] = find(b)
        nodes, counts = defaultdict(set), Counter()
        for e in selected:
            a, b = ('ACR', e['acr_id']), ('gene', e['gene_id'])
            key = find(a)
            nodes[key].update([a, b])
            counts[key] += 1
        component_edges = max(counts.values(), default=0)
        component_nodes = max((len(nodes[k]) for k in counts if counts[k] == component_edges), default=0)
        # Unique IDs in the densest 100-kb endpoint window; duplicated rows never inflate it.
        endpoint_nodes = sorted([(e['p1'], ('ACR', e['acr_id'])) for e in selected] +
                                [(e['p2'], ('gene', e['gene_id'])) for e in selected])
        window_max = 0
        active = Counter()
        j = 0
        for i, (pos, node) in enumerate(endpoint_nodes):
            while j < len(endpoint_nodes) and endpoint_nodes[j][0]-pos < 100000:
                active[endpoint_nodes[j][1]] += 1
                j += 1
            window_max = max(window_max, len(active))
            active[node] -= 1
            if not active[node]:
                del active[node]
        result.update(internal_aglink_n=result['raw_both_internal_n'],
            independent_pair_n=len(selected), acr_n=len(acrs), gene_n=len(genes),
            acr_position_n=len({e['p1'] for e in selected}), tss_position_n=len({e['p2'] for e in selected}),
            endpoint_extent_bp=positions[-1]-positions[0] if positions else 0,
            endpoint_25kb_bin_n=len({p//25000 for p in positions}),
            endpoint_100kb_max_unique_nodes=window_max,
            internal_span_covered_fraction=covered/(end-start),
            independent_pairs_per_100kb=len(selected)*100000/(end-start),
            largest_component_edges=component_edges, largest_component_nodes=component_nodes,
            local_500kb_internal_n=sum(e['right']-e['left'] <= 500000 for e in selected),
            local_1mb_internal_n=sum(e['right']-e['left'] <= 1000000 for e in selected),
            touching_internal_fraction=len(selected)/(result['unique_both_internal_n']+result['unique_single_internal_n'])
                if result['unique_both_internal_n']+result['unique_single_internal_n'] else 0.)
        return result


def assess(row, opt, rescue=False):
    requirements = [
        ('length_bp', 'min_territory_bp'), ('independent_pair_n', 'min_pairs'),
        ('acr_n', 'min_acrs'), ('gene_n', 'min_genes'),
        ('acr_position_n', 'min_acrs'), ('tss_position_n', 'min_genes'),
        ('independent_pairs_per_100kb', 'min_pairs_per_100kb'),
        ('internal_span_covered_fraction', 'min_internal_span_fraction'),
        ('endpoint_extent_bp', 'min_endpoint_extent_bp'),
        ('endpoint_25kb_bin_n', 'min_endpoint_bins'),
        ('endpoint_100kb_max_unique_nodes', 'min_local_endpoint_nodes'),
        ('largest_component_edges', 'min_component_edges'),
        ('largest_component_nodes', 'min_component_nodes'),
        ('local_1mb_internal_n', 'min_local_1mb_edges'),
        ('supported_bin_fraction', 'min_supported_fraction')]
    failed = [field for field, key in requirements if row[field] < opt[key]-1e-9]
    if row['max_zero_run_bp'] > opt['max_zero_run_bp']:
        failed.append('wide_zero_gap')
    if row['protected_zero_overlap_bp']:
        failed.append('protected_no_signal_hard_gap')
    if rescue:
        if row['background_near_zero_fraction'] < opt['min_near_zero_fraction']:
            failed.append('background_not_near_zero')
        if row['flank_background_median'] > opt['max_flank_background']:
            failed.append('flank_background_not_low')
        if row['fine_peak_excess'] < opt['min_peak_excess']:
            failed.append('weak_peak_excess')
        if row['peak_to_flank_ratio_regularized'] < opt['min_peak_to_flank_ratio']:
            failed.append('no_local_peak')
    return failed


def failed_margin(row, opt, rescue=False):
    """Descriptive ordering only; never affects candidate acceptance."""
    fail = assess(row, opt, rescue)
    fields = {'length_bp': 'min_territory_bp', 'independent_pair_n': 'min_pairs',
              'acr_n': 'min_acrs', 'gene_n': 'min_genes', 'acr_position_n': 'min_acrs',
              'tss_position_n': 'min_genes', 'independent_pairs_per_100kb': 'min_pairs_per_100kb',
              'internal_span_covered_fraction': 'min_internal_span_fraction',
              'endpoint_extent_bp': 'min_endpoint_extent_bp', 'endpoint_25kb_bin_n': 'min_endpoint_bins',
              'endpoint_100kb_max_unique_nodes': 'min_local_endpoint_nodes',
              'largest_component_edges': 'min_component_edges', 'largest_component_nodes': 'min_component_nodes',
              'local_1mb_internal_n': 'min_local_1mb_edges', 'supported_bin_fraction': 'min_supported_fraction',
              'background_not_near_zero': 'min_near_zero_fraction', 'weak_peak_excess': 'min_peak_excess',
              'no_local_peak': 'min_peak_to_flank_ratio'}
    actual = {'background_not_near_zero': 'background_near_zero_fraction', 'weak_peak_excess': 'fine_peak_excess',
              'no_local_peak': 'peak_to_flank_ratio_regularized'}
    return sum(max(0., 1-row[actual.get(field, field)]/opt[fields[field]]) if field in fields and opt[fields[field]] else 1.
               for field in fail)


def audit_previous(prior, out):
    rows = table(prior/'f04_summary/02_candidate_decisions.tsv')
    near, summary = [], []
    for method in ['T1', 'T2', 'T3']:
        selected = [r for r in rows if r['method'] == method and int(r['length_bp']) >= 300000]
        failures = Counter(s for r in selected if r['accepted'] == '0' for s in r['reject_reason'].split(','))
        for reason, count in failures.items():
            summary.append(dict(method=method, denominator=len(selected), criterion=reason, rejected_n=count))
        for r in selected:
            if r['accepted'] == '0' and len(r['reject_reason'].split(',')) == 1:
                near.append(r)
    P.write(out/'01_previous_failure_counts.tsv', summary)
    P.write(out/'02_previous_single_failure.tsv', near)


def all_summary(method, chrom, rows, t3, old, territories, cuts, size, new_pairs):
    lengths = [r['end']-r['start'] for r in rows]
    return dict(method=method, chrom=chrom, chromosome_bp=size, cru_n=len(rows),
        t3_cru_n=len(t3), original_v11_cru_n=len(old), increase_vs_t3=len(rows)-len(t3),
        new_territory_n=len(territories), new_internal_cut_n=cuts,
        covered_bp=sum(lengths), coverage_fraction=sum(lengths)/size,
        mean_bp=float(np.mean(lengths)) if lengths else np.nan,
        median_bp=float(np.median(lengths)) if lengths else np.nan,
        additional_coverage_vs_t3_bp=sum(lengths)-sum(r['end']-r['start'] for r in t3),
        new_territory_independent_pairs=new_pairs)


def main(args):
    prior, source, root = args.prior_dir.resolve(), args.source_dir.resolve(), args.outdir.resolve()
    if root in [prior, source] or prior in root.parents or source in root.parents:
        raise ValueError('Use an independent result directory')
    if not (prior/'03_PREPARATION.DONE').exists():
        raise ValueError('Previous T3 preparation is incomplete')
    if len(set(args.chroms)) != len(args.chroms):
        raise ValueError('Duplicate chromosomes')
    config = json.loads(args.settings.read_text())
    rules = json.loads((source/'f04_statistics/07_effective_rules.json').read_text())
    previous_cfg = json.loads((prior/'01_parameters.json').read_text())
    previous_stats = table(prior/'f04_summary/01_chromosome_comparison.tsv')
    signature = json.loads((source/'00_call_signature.json').read_text())
    for path, h in signature['inputs'].items():
        if not Path(path).is_file() or P.md5(path) != h:
            raise ValueError('Frozen input/function changed: '+path)
    for role in ['aglink', 'gene_bed', 'acr_bed', 'genome_sizes']:
        supplied = getattr(args, role)
        expected = [h for p, h in signature['inputs'].items() if Path(p).resolve() == supplied.resolve()]
        if expected != [P.md5(supplied)]:
            raise ValueError('Input differs from original V11: '+str(supplied))
    validation = json.loads((source.parent/'f04_block_enrichment/01_parameters.json').read_text())
    for path in [args.block_bed, args.block_helper, args.genome_sizes]:
        if validation['inputs'].get(str(path)) != P.md5(path):
            raise ValueError('Original block validation input changed: '+str(path))
    if (args.perm_n, args.seed, args.local_window_bp) != (validation['perm_n'], validation['seed'], validation['local_window_bp']):
        raise ValueError('Use original local-permutation settings')
    sizes = {r.split()[0]: int(r.split()[1]) for r in args.genome_sizes.read_text().splitlines() if r.strip()}
    if set(args.chroms)-set(sizes) or args.window_bp != 5000000:
        raise ValueError('Invalid chromosome / use 5-Mb atlas windows')
    for method in ['T4', 'T5', 'T6']:
        opt = config[method]
        if opt['min_territory_bp'] != 300000 or opt['min_pairs'] < 2 or opt['max_zero_run_bp'] > 10000:
            raise ValueError('Keep 300-kb length, >=2 distinct pairs and zero-gap protection')
    old_path = source/'f03_results/f06_V11_300kb/02_CRU.bed'
    original = P.read_bed(old_path)
    if len(original) != args.expected_cru_n:
        raise ValueError('Original population differs from expected count')
    data_files = [old_path, source/'00_call_signature.json', prior/'01_parameters.json',
                  prior/'f04_summary/01_chromosome_comparison.tsv', prior/'f04_summary/02_candidate_decisions.tsv',
                  source/'f04_statistics/07_effective_rules.json', args.settings, args.aglink, args.gene_bed,
                  args.acr_bed, args.genome_sizes, args.block_bed, args.block_helper,
                  source/'f01_shared_input/01_count_only_aglink.annotated.bedpe', Path(P.__file__)]
    for c in args.chroms:
        data_files += [prior/f'f03_tests/T3/{c}.{name}.bed' for name in ['CRU', 'added_territory']]
        data_files += [prior/f'f02_diagnostics/{c}/05_verified_refined_gaps.tsv', source/f'f02_inputs/{c}.profiles.npz',
                       prior/f'f06_plot_data/{c}.density.tsv.gz', prior/f'f06_plot_data/{c}.aglink.tsv.gz']
    data_files += [p for p in HERE.iterdir() if p.is_file()]
    hashes = {str(p): P.md5(p) for p in data_files}
    run_signature = dict(inputs=hashes, options={k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items() if k != 'resume'})
    marker = root/'01_parameters.json'
    if root.exists() and any(root.iterdir()):
        if not args.resume or not marker.exists() or json.loads(marker.read_text())['signature'] != run_signature:
            raise ValueError('Existing run signature differs; preserve results and use a new folder')
    for d in ['f01_scripts', 'f02_diagnostics', 'f03_tests', 'f04_summary', 'f05_block_validation', 'f06_plot_data', 'f07_figures', 'f08_checkpoints']:
        (root/d).mkdir(parents=True, exist_ok=True)
    for p in HERE.iterdir():
        if p.is_file():
            shutil.copy2(p, root/'f01_scripts'/p.name)
    shutil.copy2(Path(P.__file__), root/'f01_scripts/05_previous_recall_functions.py')
    (root/'f01_scripts/06_frozen_V11_rules.json').write_text(json.dumps(rules, indent=2)+'\n')
    (root/'f01_scripts/07_previous_T3_density_rules.json').write_text(json.dumps(previous_cfg['settings_content'], indent=2)+'\n')
    shutil.copy2(source/'00_call_signature.json', root/'f01_scripts/08_original_calling_signature.json')
    marker.write_text(json.dumps(dict(signature=run_signature, source_dir=str(source), block_bed=str(args.block_bed),
        settings=config, prior_dir=str(prior)), indent=2)+'\n')
    audit_previous(prior, root/'f02_diagnostics')
    engine = _load('m11_engine')
    M, H = engine.M, engine.H
    edges = load_edges(source/'f01_shared_input/01_count_only_aglink.annotated.bedpe', args.chroms)
    all_decisions, metrics, territory_details, invariants, block_groups, windows = [], [], [], [], {}, []
    all_beds = {}
    for c in args.chroms:
        print('Support sensitivity:', c, flush=True)
        checkpoint = root/f'f08_checkpoints/{c}.pickle'
        if args.resume and checkpoint.exists():
            with checkpoint.open('rb') as f:
                result = pickle.load(f)
            if result['signature'] != run_signature:
                raise ValueError('Checkpoint differs: '+c)
        else:
            size, es = sizes[c], edges[c]
            index = EdgeIndex(es)
            t3 = sorted(P.read_bed(prior/f'f03_tests/T3/{c}.CRU.bed'), key=lambda r: r['start'])
            old = [r for r in original if r['chrom'] == c]
            if not P.disjoint(t3) or any(r not in t3 for r in old):
                raise ValueError('T3 does not preserve original V11')
            expected_t3 = [int(r['cru_n']) for r in previous_stats if r['chrom'] == c and r['method'] == 'T3']
            if expected_t3 != [len(t3)]:
                raise ValueError('T3 BED count differs from saved previous summary: '+c)
            with np.load(source/f'f02_inputs/{c}.profiles.npz') as z:
                raw, smooth, bg = z['full_raw'].copy(), z['full_smooth'].copy(), z['full_background'].copy()
            endpoint = np.zeros(len(raw))
            np.add.at(endpoint, index.p1//STEP, 1)
            np.add.at(endpoint, index.p2//STEP, 1)
            ep = uniform_filter1d(endpoint, 45, mode='constant')*45
            sig = M.D.territory_signals(raw, smooth, ep, SimpleNamespace(**rules['primary']))
            if not np.allclose(sig['background'], bg, atol=1e-9, rtol=1e-9):
                raise ValueError('Saved background differs')
            protected_zero = np.zeros(len(raw), bool)
            zero_records = []
            for a, b in P.runs(raw <= 1e-9):
                if min(size, b*STEP)-a*STEP >= rules['hard_gap']['hard_gap_min_bp']:
                    protected_zero[a:b] = True
                    zero_records.append(dict(chrom=c, start=a*STEP, end=min(size,b*STEP),
                        length_bp=min(size,b*STEP)-a*STEP, reason='raw_zero_ge_original_hard_gap_duration',
                        **index.describe(a*STEP, min(size,b*STEP))))
            P.write(root/f'f02_diagnostics/{c}_protected_zero_hard_gaps.tsv', zero_records,
                    None if zero_records else ['chrom','start','end','length_bp','reason'])
            gaps = table(prior/f'f02_diagnostics/{c}/05_verified_refined_gaps.tsv')
            for r in gaps:
                r['start'], r['end'] = int(r['start']), int(r['end'])
            gap_mask = P.mask_intervals(gaps, size)
            acrs, genes = H.annotation_positions(c, args.gene_bed, args.acr_bed)
            result = dict(signature=run_signature, decisions=[], stats=[], territories=[], invariants=[], beds={}, groups={})
            inherited = list(t3)
            territories, cumulative_cuts = [], 0
            # Each method retains its predecessor exactly. T4/T5 use unchanged T3
            # density routes; only T6 introduces an additional lower-density route.
            for method in METHODS:
                if method != 'T3':
                    opt = config[method]
                    routes = [('T2_density', previous_cfg['settings_content']['T2'], False),
                              ('T3_low_background', previous_cfg['settings_content']['T3_rescue'], True)]
                    if method == 'T6':
                        routes.append(('T6_low_background_peak', config['T6_density'], True))
                    for route, density_opt, low_bg in routes:
                        available = ~P.mask_intervals(inherited, size) & ~protected_zero
                        proposed = P.proposal_mask(raw, smooth, bg, available, density_opt, low_bg)
                        # Runs within a route are disjoint; no greedy overlapping candidate selection.
                        for a, b in P.runs(proposed):
                            row = dict(method=method, stage=route, chrom=c, **P.describe(raw, smooth, bg, ep, a, b, size))
                            row.update(index.describe(row['start'], row['end']))
                            flank = 250000//STEP
                            flank_indices = np.r_[np.arange(max(0, a-flank), a), np.arange(b, min(len(raw), b+flank))]
                            bg_flank = bg[flank_indices]
                            fine_flank = smooth[0, flank_indices]
                            peak = float(smooth[0, a:b].max())
                            threshold = np.maximum(density_opt['absolute_floor'], bg[a:b]+density_opt['background_sigma']*np.sqrt(bg[a:b]+1))
                            sustained = (raw[a:b] >= density_opt['absolute_floor']-1e-9) & (np.sum(smooth[:, a:b] >= threshold-1e-9, axis=0) >= 2)
                            if low_bg:
                                sustained &= bg[a:b] <= density_opt['near_zero_background']
                            row.update(supported_bin_fraction=float(sustained.mean()),
                                background_near_zero_fraction=float(np.mean(bg[a:b] <= density_opt.get('near_zero_background', .1))),
                                background_near_zero_cutoff=density_opt.get('near_zero_background', .1),
                                original_hard_gap_fraction=float(gap_mask[a:b].mean()),
                                protected_zero_overlap_bp=int(protected_zero[a:b].sum()*STEP),
                                flank_background_median=float(np.median(bg_flank)) if len(bg_flank) else 0.,
                                fine_peak_excess=float((smooth[0, a:b]-bg[a:b]).max()),
                                peak_to_flank_ratio_regularized=(peak+.25)/(float(np.median(fine_flank))+.25) if len(fine_flank) else (peak+.25)/.25)
                            rescue = route == 'T6_low_background_peak'
                            failed = assess(row, opt, rescue)
                            row.update(accepted=int(not failed), reject_reason=','.join(failed) or 'qualified_internal_topology',
                                failed_condition_n=len(failed), failure_margin=failed_margin(row, opt, rescue),
                                original_hard_gap_status='supported_low_signal_overlap' if gap_mask[a:b].any() else 'outside_saved_hard_gap')
                            result['decisions'].append(row)
                            if failed:
                                continue
                            parent = dict(chrom=c, start=row['start'], end=row['end'], cru_id=f'{c}_{method}_territory_{len(territories)+1}')
                            if args.split_new == 'frozen':
                                local = [e for e in es if e['left'] < parent['end'] and e['right'] >= parent['start']]
                                rr, bb, cc, _ = M.solve(parent, sig, local, size, SimpleNamespace(**rules['primary']),
                                    M.options(rules['shape_settings'], 300000), acrs, genes, method, rules['max_rounds'])
                                P.write(root/f'f02_diagnostics/{c}_{method}_{len(territories)+1:04d}_split_candidates.tsv', cc)
                            else:
                                rr, bb = [parent.copy()], []
                            rr.sort(key=lambda r: r['start'])
                            assert rr[0]['start'] == parent['start'] and rr[-1]['end'] == parent['end']
                            assert all(x['end'] == y['start'] for x, y in zip(rr, rr[1:]))
                            assert all(r['end']-r['start'] >= 300000 for r in rr)
                            for j, r in enumerate(rr, 1):
                                r['cru_id'] = parent['cru_id']+f'_CRU_{j}'
                            inherited += rr
                            territories.append(parent)
                            cumulative_cuts += len(bb)
                            result['territories'].append(dict(**row, territory_id=parent['cru_id'], child_cru_n=len(rr), internal_cut_n=len(bb)))
                inherited.sort(key=lambda r: r['start'])
                assert P.disjoint(inherited) and all(r in inherited for r in t3)
                extra = [r for r in inherited if r not in t3]
                extra_mask = P.mask_intervals(extra, size)
                assert not (extra_mask & protected_zero).any()
                result['beds'][method] = list(inherited)
                pairs = sum(index.describe(r['start'], r['end'])['independent_pair_n'] for r in territories)
                result['stats'].append(all_summary(method, c, inherited, t3, old, territories, cumulative_cuts, size, pairs))
                old_points = {p for r in old for p in [r['start'], r['end']]}
                t3_points = {p for r in t3 for p in [r['start'], r['end']]}
                current_points = {p for r in inherited for p in [r['start'], r['end']]}
                result['invariants'].append(dict(chrom=c, method=method, t3_record_n=len(t3),
                    lost_or_changed_t3_records=sum(r not in inherited for r in t3), old_boundary_missing_n=len(old_points-current_points),
                    overlap_n=0, protected_zero_filled_bp=int((extra_mask & protected_zero).sum()*STEP)))
                groups = P.boundary_groups(old, [r for r in inherited if r not in old],
                    P.read_bed(prior/f'f03_tests/T3/{c}.added_territory.bed')+territories, size)
                # Keep the T0/T3 comparison units distinct from genuinely NEW this-round points.
                groups['added_since_T3'] = {p for p in current_points-t3_points if 0 < p < size}
                result['groups'][method] = groups
            temp = checkpoint.with_suffix('.tmp')
            with temp.open('wb') as f:
                pickle.dump(result, f, protocol=4)
            temp.replace(checkpoint)
        all_decisions += result['decisions']
        metrics += result['stats']
        territory_details += result['territories']
        invariants += result['invariants']
        for method in METHODS:
            rows = result['beds'][method]
            all_beds[(method, c)] = rows
            t3 = result['beds']['T3']
            P.save_bed(root/f'f03_tests/{method}/{c}.CRU.bed', rows)
            P.save_bed(root/f'f03_tests/{method}/{c}.added_since_T3.bed', [r for r in rows if r not in t3])
            for group, points in result['groups'][method].items():
                block_groups.setdefault((method, group), {})[c] = points
        for suffix in ['density.tsv.gz', 'aglink.tsv.gz']:
            shutil.copy2(prior/f'f06_plot_data/{c}.{suffix}', root/f'f06_plot_data/{c}.{suffix}')
        for start in range(0, sizes[c], args.window_bp):
            windows.append(dict(chrom=c, start=start, end=min(sizes[c], start+args.window_bp)))
    for method in METHODS:
        mm = [r for r in metrics if r['method'] == method]
        values = [r['end']-r['start'] for c in args.chroms for r in all_beds[(method, c)]]
        totals = {k: sum(r[k] for r in mm) for k in ['chromosome_bp', 'cru_n', 't3_cru_n', 'original_v11_cru_n', 'increase_vs_t3',
                  'new_territory_n', 'new_internal_cut_n', 'covered_bp', 'additional_coverage_vs_t3_bp', 'new_territory_independent_pairs']}
        metrics.append(dict(method=method, chrom='Selected_chromosomes', **totals,
            coverage_fraction=totals['covered_bp']/totals['chromosome_bp'], mean_bp=float(np.mean(values)), median_bp=float(np.median(values))))
    P.write(root/'f04_summary/01_chromosome_comparison.tsv', metrics)
    P.write(root/'f04_summary/02_candidate_decisions.tsv', all_decisions,
            None if all_decisions else ['method', 'chrom', 'start', 'end', 'accepted'])
    P.write(root/'f04_summary/03_new_territory_support.tsv', territory_details,
            None if territory_details else ['method', 'chrom', 'start', 'end', 'accepted'])
    P.write(root/'f04_summary/04_boundary_invariants.tsv', invariants)
    reasons, passing = [], []
    for method in ['T4', 'T5', 'T6']:
        for scope in args.chroms+['Selected_chromosomes']:
            rr = [r for r in all_decisions if r['method'] == method and (scope == 'Selected_chromosomes' or r['chrom'] == scope)]
            for subset, ss in [('all_proposals', rr), ('length_ge300kb', [r for r in rr if r['length_bp'] >= 300000])]:
                n = len(ss)
                passing.append(dict(method=method, scope=scope, subset=subset, candidate_n=n, accepted_n=sum(r['accepted'] for r in ss),
                    accepted_fraction=sum(r['accepted'] for r in ss)/n if n else np.nan))
                for cause, count in Counter(s for r in ss if not r['accepted'] for s in r['reject_reason'].split(',')).items():
                    reasons.append(dict(method=method, scope=scope, subset=subset, candidate_n=n, criterion=cause, reject_n=count))
    P.write(root/'f04_summary/05_rejection_counts.tsv', reasons)
    P.write(root/'f04_summary/06_acceptance_rates.tsv', passing)
    gap_records = [r for r in all_decisions if r['original_hard_gap_fraction'] > 0]
    P.write(root/'f04_summary/07_saved_hard_gap_candidates.tsv', gap_records,
            None if gap_records else ['method', 'chrom', 'start', 'end', 'accepted'])
    gap_summary = []
    for method in ['T4', 'T5', 'T6']:
        rr = [r for r in gap_records if r['method'] == method]
        gap_summary.append(dict(method=method, criterion='overlaps_saved_hard_gap', candidate_n=len(rr),
            accepted_n=sum(r['accepted'] for r in rr), rejected_n=sum(not r['accepted'] for r in rr)))
        for cause, count in Counter(s for r in rr if not r['accepted'] for s in r['reject_reason'].split(',')).items():
            gap_summary.append(dict(method=method, criterion=cause, candidate_n=len(rr), accepted_n=0, rejected_n=count))
    P.write(root/'f04_summary/09_hard_gap_rejection_summary.tsv', gap_summary)
    distributions = []
    for method in ['T4', 'T5', 'T6']:
        rr = [r for r in territory_details if METHODS.index(r['method']) <= METHODS.index(method)]
        for field in ['length_bp', 'independent_pair_n', 'acr_n', 'gene_n', 'unique_both_internal_n',
                      'unique_single_internal_n', 'unique_span_only_n', 'raw_both_internal_n', 'raw_single_internal_n', 'raw_span_only_n']:
            vv = [r[field] for r in rr]
            distributions.append(dict(method=method, field=field, territory_n=len(rr), median=float(np.median(vv)) if vv else np.nan,
                mean=float(np.mean(vv)) if vv else np.nan, p25=float(np.quantile(vv, .25)) if vv else np.nan,
                p75=float(np.quantile(vv, .75)) if vv else np.nan, max=max(vv) if vv else np.nan))
    P.write(root/'f04_summary/08_support_distributions.tsv', distributions)
    accepted = list(territory_details)
    random.Random(args.seed).shuffle(accepted)
    accepted = accepted[:args.example_n]
    # One record per physical candidate; rank first by failure count, then deficit.
    unique_rejected = {}
    for r in all_decisions:
        if r['accepted'] or r['length_bp'] < 300000:
            continue
        if any(t['chrom'] == r['chrom'] and t['start'] < r['end'] and r['start'] < t['end'] for t in territory_details):
            continue
        key = r['chrom'], r['start'], r['end']
        rank = r['failed_condition_n'], r['failure_margin'], r['method'], r['stage']
        if key not in unique_rejected or rank < unique_rejected[key][0]:
            unique_rejected[key] = rank, r
    rejected = [r for _, r in sorted(unique_rejected.values(), key=lambda x: x[0])][:args.example_n]
    single_failures = [r for _, r in unique_rejected.values() if r['failed_condition_n'] == 1]
    P.write(root/'f04_summary/10_still_rejected_single_failure.tsv', single_failures,
            None if single_failures else ['method', 'chrom', 'start', 'end', 'reject_reason'])
    examples = []
    for kind, rows in [('new_territory', accepted), ('near_pass_rejected', rejected)]:
        for i, r in enumerate(rows, 1):
            examples.append(dict(r, example_id=f'{kind}_{i:02d}', example_kind=kind))
    P.write(root/'f06_plot_data/02_examples.tsv', examples, None if examples else ['example_id', 'example_kind', 'chrom', 'start', 'end', 'accepted'])
    P.write(root/'f06_plot_data/01_windows.tsv', windows)
    if not args.skip_block:
        P.validate_blocks(args, block_groups, sizes, root/'f05_block_validation')
    if any(P.md5(p) != h for p, h in hashes.items()):
        raise ValueError('Source changed during analysis')
    P.write(root/'02_input_MD5.tsv', [dict(path=p, before_md5=h, after_md5=P.md5(p)) for p, h in hashes.items()])
    report = ['# T3/T4/T5/T6 support sensitivity', '', 'Selected chromosomes only: '+', '.join(args.chroms),
        'All inherited V11/T3 records retained; hard-gap identification unchanged.',
        'T4/T5 use unchanged T3 density routes. T6 adds a lower-density peak-rescue route.',
        'T4/T5/T6 are cumulative, so summaries and candidate-stage denominators are not independent.', '',
        '| Method | CRUs | Increase vs T3 | New territories | New internal cuts | Coverage |',
        '| --- | ---: | ---: | ---: | ---: | ---: |']
    for r in metrics:
        if r['chrom'] == 'Selected_chromosomes':
            report.append('| {method} | {cru_n} | {increase_vs_t3} | {new_territory_n} | {new_internal_cut_n} | {pct:.2f}% |'.format(**r, pct=100*r['coverage_fraction']))
    report += ['', f'Seeded new-territory examples: {len(accepted)}; near-pass rejected examples: {len(rejected)}.',
        'Near-pass means descriptive ranking by fewest failed conditions then normalized deficit, not a relaxed acceptance rule.',
        'Rejected examples exclude any overlap with territories accepted by any new method.',
        'See f04_summary/04_boundary_invariants.tsv for missing records/overlap/zero-gap checks.',
        'See f05_block_validation for original local-null observed/random data; no selection uses blocks.',
        'In block tables added_all is relative to original V11; added_since_T3 is new this round.']
    (root/'04_REPORT.md').write_text('\n'.join(report)+'\n')
    (root/'03_PREPARATION.DONE').write_text('Support sensitivity finished; no replacement of original results.\n')
    for r in metrics:
        if r['chrom'] == 'Selected_chromosomes':
            print(r['method'], r['cru_n'], 'increase vs T3', r['increase_vs_t3'], flush=True)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    for k in ['source-dir', 'prior-dir', 'settings', 'aglink', 'gene-bed', 'acr-bed', 'genome-sizes', 'block-bed', 'block-helper', 'outdir']:
        p.add_argument('--'+k, type=Path, required=True)
    p.add_argument('--chroms', nargs='+', default=['4B', '4D', '5D', '7D', '3A', '5A', '1B', '2D'])
    p.add_argument('--expected-cru-n', type=int, default=3159)
    p.add_argument('--split-new', choices=['none', 'frozen'], default='frozen')
    p.add_argument('--window-bp', type=int, default=5000000)
    p.add_argument('--perm-n', type=int, default=1000)
    p.add_argument('--seed', type=int, default=20260805)
    p.add_argument('--local-window-bp', type=int, default=10000000)
    p.add_argument('--block-cutoffs', nargs='+', type=int, default=[25000, 50000, 100000, 200000])
    p.add_argument('--example-n', type=int, default=20)
    p.add_argument('--skip-block', action='store_true')
    p.add_argument('--resume', action='store_true')
    main(p.parse_args())
