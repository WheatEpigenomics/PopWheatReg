# m01_density_signal: span/endpoint density signal
# Frozen rule module, copied unchanged from the published calling chain
# (source: f02_span_density_pilot/01_call_span_density_pilot.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Independent, count-only agLINK span-density pilot; blocks are evaluation only."""

import argparse
import bisect
import csv
import gzip
import hashlib
import json
import math
import random
import runpy
import shutil
import sys
from collections import Counter, defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np
from scipy.ndimage import uniform_filter1d
from scipy.signal import find_peaks, peak_widths


METHODS = {'V1': 'f02_V1_span_valley', 'V2': 'f03_V2_recursive_background',
           'V3': 'f04_V3_crossing_filter', 'M3': 'f05_M3_reference'}
COLORS = {'V1': '#96CEDD', 'V2': '#E89B6C', 'V3': '#5FB6BD', 'M3': '#B086CF'}


def log(text):
    print(f'[{datetime.now().isoformat(timespec="seconds")}] {text}', flush=True)


def write_tsv(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(dict.fromkeys(k for row in rows for k in row))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'wt', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def read_bed(path, chroms):
    result = defaultdict(list)
    with open(path) as handle:
        for lineno, line in enumerate(handle, 1):
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if fields[0] not in chroms:
                continue
            start, end = map(int, fields[1:3])
            if start < 0 or end <= start:
                raise ValueError(f'Invalid BED interval {path}:{lineno}')
            result[fields[0]].append({'chrom': fields[0], 'start': start, 'end': end,
                                      'cru_id': fields[3] if len(fields) > 3 else f'row{lineno}'})
    for chrom in chroms:
        result[chrom].sort(key=lambda r: (r['start'], r['end']))
    return result


def metadata(path, chroms):
    result = {}
    with open(path) as handle:
        for line in handle:
            if not line.strip() or line.startswith('#'):
                continue
            fields = line.split()
            if len(fields) >= 4 and fields[0] in chroms:
                result[fields[3]] = (fields[0], int(fields[1]), int(fields[2]), fields[5] if len(fields) > 5 else '.')
    return result


def load_count_only_edges(args, sizes, shared):
    genes, acrs = metadata(args.gene_bed, args.chroms), metadata(args.acr_bed, args.chroms)
    by_chrom, audit, pair_counts = defaultdict(list), Counter(), Counter()
    expected_bedpe = Counter()
    with open(args.aglink) as handle:
        for lineno, row in enumerate(csv.reader(handle, delimiter='\t'), 1):
            if not row or row[0].startswith('#'):
                continue
            audit['input_rows'] += 1
            if len(row) < 2:
                raise ValueError(f'agLINK row {lineno} has fewer than two columns')
            gene, acr = row[:2]
            if gene not in genes or acr not in acrs:
                audit['not_both_in_selected_metadata'] += 1
                continue
            gc, gs, ge, strand = genes[gene]
            ac, ast, ae, _ = acrs[acr]
            if ac != gc:
                audit['selected_interchromosomal_rows'] += 1
                continue
            tss, amid = (ge - 1 if strand == '-' else gs), (ast + ae) // 2
            p1, p2 = amid, tss
            left, right = sorted((p1, p2))
            if not 0 <= left <= right < sizes[gc]:
                raise ValueError(f'Out-of-range input anchors, line {lineno}')
            edge = {'chrom': gc, 'p1': p1, 'p2': p2, 'left': left, 'right': right,
                    'span_start': left, 'span_end': right + 1, 'distance_bp': right - left,
                    'node1': acr, 'node2': gene, 'edge_id': f'aglink_input_{lineno}',
                    'edge_type': 'aglink', 'source_line': lineno, 'weight': 1.0,
                    'acr_start': ast, 'acr_end': ae, 'gene_start': gs, 'gene_end': ge,
                    'gene_strand': strand}
            by_chrom[gc].append(edge)
            pair_counts[(gene, acr)] += 1
            anchors = sorted([(gc, gs, ge), (ac, ast, ae)])
            expected_bedpe[tuple(anchors[0] + anchors[1])] += 1
            audit[f'{gc}_retained'] += 1
    audit['duplicate_pair_excess_rows'] = sum(n - 1 for n in pair_counts.values())
    # The original six-column BEDPE lacks IDs. Verify its coordinate multiset,
    # then retain ID-bearing point anchors from the same significant-pair table.
    if args.bedpe_reference:
        actual = Counter()
        with open(args.bedpe_reference) as handle:
            for line in handle:
                fields = line.split()
                if len(fields) < 6 or fields[0] not in args.chroms or fields[3] not in args.chroms:
                    continue
                if fields[0] != fields[3]:
                    continue
                anchors = sorted([(fields[0], int(fields[1]), int(fields[2])),
                                  (fields[3], int(fields[4]), int(fields[5]))])
                actual[tuple(anchors[0] + anchors[1])] += 1
        audit['bedpe_reference_missing_rows'] = sum((expected_bedpe - actual).values())
        audit['bedpe_reference_extra_rows'] = sum((actual - expected_bedpe).values())
        if audit['bedpe_reference_missing_rows'] or audit['bedpe_reference_extra_rows']:
            raise ValueError(f'Original BEDPE does not match selected significant pairs: {audit}')
    rows = []
    for chrom in args.chroms:
        by_chrom[chrom].sort(key=lambda r: (r['left'], r['right'], r['source_line']))
        if not by_chrom[chrom]:
            raise ValueError(f'No agLINKs on {chrom}')
        for e in by_chrom[chrom]:
            rows.append({'chrom1': chrom, 'start1': e['p1'], 'end1': e['p1'] + 1,
                         'chrom2': chrom, 'start2': e['p2'], 'end2': e['p2'] + 1,
                         'acr_id': e['node1'], 'gene_id': e['node2'], **e})
    write_tsv(shared / '01_count_only_aglink.annotated.bedpe', rows)
    write_tsv(shared / '02_input_audit.tsv', [{'metric': k, 'value': v} for k, v in sorted(audit.items())])
    return by_chrom


def span_coverage(edges, size, step, predicate=lambda e: True):
    n = math.ceil(size / step)
    diff = np.zeros(n + 1)
    partial = np.zeros(n)
    for e in edges:
        if not predicate(e):
            continue
        start, end = e['span_start'], e['span_end']
        a, b = start // step, (end - 1) // step
        if a == b:
            partial[a] += end - start
        else:
            partial[a] += (a + 1) * step - start
            partial[b] += end - b * step
            diff[a + 1] += 1
            diff[b] -= 1
    lengths = np.minimum(step, size - np.arange(n) * step)
    return np.cumsum(diff[:-1]) + partial / lengths


def runs(mask):
    edges = np.diff(np.r_[False, mask, False].astype(int))
    return list(zip(np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)))


def signals(edges, size, args):
    raw = span_coverage(edges, size, args.bin_bp)
    endpoint = np.zeros(len(raw))
    for e in edges:
        endpoint[e['p1'] // args.bin_bp] += 1
        endpoint[e['p2'] // args.bin_bp] += 1
    widths = [max(1, int(round(bp / args.bin_bp)) | 1) for bp in args.smooth_bp]
    smooth = [uniform_filter1d(raw, width, mode='nearest') for width in widths]
    ep = uniform_filter1d(endpoint, widths[-1], mode='constant') * widths[-1]
    positive = smooth[-1][smooth[-1] > 0]
    threshold = max(args.min_span_count, float(np.quantile(positive, args.active_quantile)))
    # A single long edge has depth <=1; no gap filling can bridge its span.
    active = (raw >= args.min_span_count) & (smooth[-1] >= threshold) & (ep >= args.min_endpoint_count)
    return {'raw': raw, 'endpoint': endpoint, 'endpoint_broad': ep, 'smooth': smooth,
            'long_ge1mb': span_coverage(edges, size, args.bin_bp, lambda e: e['distance_bp'] >= 1_000_000),
            'widths': widths, 'active': active, 'threshold': threshold}


def arrays(edges):
    return np.array([e['left'] for e in edges]), np.array([e['right'] for e in edges])


def internal_n(lefts, rights, start, end):
    return int(np.count_nonzero((lefts >= start) & (rights < end)))


def candidate_valleys(y, start, end, args, scale):
    a, b = max(0, start // args.bin_bp), min(len(y), math.ceil(end / args.bin_bp))
    if b - a < 3:
        return []
    values = y[a:b]
    peaks, props = find_peaks(-values, prominence=args.min_prominence, distance=max(1, args.min_cru_bp // args.bin_bp))
    if len(peaks) == 0:
        return []
    widths, _, lip, rip = peak_widths(-values, peaks, rel_height=0.5)
    result = []
    flank = max(1, args.shoulder_bp // args.bin_bp)
    for j, p in enumerate(peaks):
        idx = a + int(p)
        pos = idx * args.bin_bp + args.bin_bp // 2
        lpeak = float(np.max(y[max(a, idx - flank):idx]))
        rpeak = float(np.max(y[idx + 1:min(b, idx + flank + 1)]))
        floor = float(y[idx])
        shoulder = min(lpeak, rpeak)
        result.append({'boundary_pos': pos, 'scale': scale, 'floor': floor,
                       'left_peak': lpeak, 'right_peak': rpeak,
                       'relative_depth': (shoulder - floor) / max(shoulder, 1e-12),
                       'peak_balance': min(lpeak, rpeak) / max(lpeak, rpeak, 1e-12),
                       'prominence': float(props['prominences'][j]),
                       'valley_start': int((a + lip[j]) * args.bin_bp),
                       'valley_end': int((a + rip[j]) * args.bin_bp),
                       'width_bp': float(widths[j] * args.bin_bp)})
    return result


def split_intervals(territory, cuts):
    points = [territory['start']] + sorted(cuts) + [territory['end']]
    return [(a, b) for a, b in zip(points, points[1:])]


def viable_cut(pos, cuts, territory, lefts, rights, args):
    parts = [territory['start']] + sorted(cuts) + [territory['end']]
    i = bisect.bisect_left(parts, pos)
    if i == 0 or i >= len(parts) or parts[i] == pos:
        return False
    start, end = parts[i - 1], parts[i]
    return (min(pos - start, end - pos) >= args.min_cru_bp
            and internal_n(lefts, rights, start, pos) >= args.min_internal
            and internal_n(lefts, rights, pos, end) >= args.min_internal)


def evaluate_crossing(pos, lefts, rights, args):
    flank = args.cross_flank_bp
    il = internal_n(lefts, rights, max(0, pos - flank), pos)
    ir = internal_n(lefts, rights, pos, pos + flank)
    global_c = int(np.count_nonzero((lefts < pos) & (rights > pos)))
    local_c = int(np.count_nonzero((lefts < pos) & (rights > pos)
                                 & (lefts >= pos - flank) & (rights < pos + flank)))
    weakest = min(il, ir)
    fraction = global_c / (global_c + weakest) if global_c + weakest else 1.0
    return {'left_internal': il, 'right_internal': ir, 'crossing_global': global_c,
            'crossing_local': local_c, 'weak_side_cross_fraction': fraction,
            'crossing_pass': int(weakest >= args.min_internal and fraction <= args.max_cross_fraction)}


def independent_background(candidates, y, args):
    chosen = []
    for c in sorted(candidates, key=lambda c: (c['floor'], c['boundary_pos'])):
        if all(abs(c['boundary_pos'] - old['boundary_pos']) >= args.background_sep_bp
               and (c['valley_end'] <= old['valley_start'] or old['valley_end'] <= c['valley_start']) for old in chosen):
            chosen.append(c)
        if len(chosen) == args.background_regions:
            break
    floors = []
    for c in chosen:
        lo = c['valley_start'] // args.bin_bp
        hi = min(len(y), max(lo + 1, math.ceil(c['valley_end'] / args.bin_bp)))
        floors.append(float(np.median(y[lo:hi])))
    return (float(np.median(floors)) if len(floors) >= args.min_background_regions else None), chosen


def build_methods(chrom, size, edges, sig, args):
    lefts, rights = arrays(edges)
    territories = []
    for a, b in runs(sig['active']):
        start, end = a * args.bin_bp, min(size, b * args.bin_bp)
        if end - start >= args.min_cru_bp and internal_n(lefts, rights, start, end) >= args.min_internal:
            territories.append({'chrom': chrom, 'start': start, 'end': end,
                                'territory_id': f'{chrom}_T{len(territories) + 1}'})
    results = {method: {'crus': [], 'boundaries': [], 'candidates': [], 'backgrounds': [], 'actions': []} for method in ['V1', 'V2', 'V3']}
    for territory in territories:
        tid = territory['territory_id']
        common = {'chrom': chrom, 'territory_id': tid, 'territory_start': territory['start'], 'territory_end': territory['end']}
        cuts1, cut_info = [], {}
        pool = candidate_valleys(sig['smooth'][0], territory['start'], territory['end'], args, 'fine')
        for c in sorted(pool, key=lambda c: (-c['relative_depth'], -c['prominence'], c['boundary_pos'])):
            basic = c['relative_depth'] >= args.v1_depth and c['width_bp'] >= args.valley_min_width_bp and c['peak_balance'] >= args.min_peak_balance
            passed = basic and viable_cut(c['boundary_pos'], cuts1, territory, lefts, rights, args)
            row = {**common, **c, 'stage': 'V1', 'accepted': int(passed),
                   'reject_reason': '' if passed else ('valley_shape' if not basic else 'daughter_support_or_spacing')}
            results['V1']['candidates'].append(row)
            if passed:
                cuts1.append(c['boundary_pos'])
                cut_info[c['boundary_pos']] = row
        cuts2 = list(cuts1)
        actions, bgrows, extra = [], [], []
        todo = [(start, end, 0) for start, end in split_intervals(territory, cuts1)]
        while todo:
            start, end, depth = todo.pop()
            if end - start < args.long_region_bp:
                continue
            pool2 = candidate_valleys(sig['smooth'][1], start, end, args, 'medium')
            eligible_bg = [c for c in pool2 if c['width_bp'] >= args.valley_min_width_bp]
            bg, chosen = independent_background(eligible_bg, sig['smooth'][1], args)
            for c in chosen:
                bgrows.append({**common, 'parent_start': start, 'parent_end': end, 'depth': depth,
                               'background_region_start': c['valley_start'], 'background_region_end': c['valley_end'],
                               'local_background': bg if bg is not None else '', 'background_region_n': len(chosen)})
            selected = []
            for c in pool2:
                shoulder = min(c['left_peak'], c['right_peak'])
                bgdepth = (shoulder - c['floor']) / max(shoulder - bg, 1.0) if bg is not None else -1
                valid_shape = (bg is not None and shoulder > bg + args.min_prominence and bgdepth >= args.v2_background_depth
                               and c['width_bp'] >= args.valley_min_width_bp and c['peak_balance'] >= args.min_peak_balance)
                valid = valid_shape and viable_cut(c['boundary_pos'], cuts2, territory, lefts, rights, args)
                row = {**common, **c, 'stage': 'V2', 'parent_start': start, 'parent_end': end,
                       'depth': depth, 'background_region_n': len(chosen), 'local_background': bg if bg is not None else '',
                       'background_depth': bgdepth, 'accepted': 0,
                       'reject_reason': 'missing_independent_background' if bg is None else 'shape_or_daughter_support'}
                extra.append(row)
                if valid:
                    selected.append(row)
            if selected and depth < args.max_recursion_depth:
                best = max(selected, key=lambda c: (c['background_depth'], c['prominence'], -c['boundary_pos']))
                best['accepted'], best['reject_reason'] = 1, ''
                pos = best['boundary_pos']
                cuts2.append(pos)
                cut_info[pos] = best
                todo.extend([(start, pos, depth + 1), (pos, end, depth + 1)])
                status = 'background_supported_split'
            else:
                status = 'kept_long_no_supported_valley' if depth < args.max_recursion_depth else 'kept_long_recursion_safety_limit'
            actions.append({**common, 'parent_start': start, 'parent_end': end, 'depth': depth,
                            'action': status, 'candidate_n': len(pool2), 'background_region_n': len(chosen)})
        cuts3 = []
        evidence = {}
        for pos in sorted(cuts2):
            ev = evaluate_crossing(pos, lefts, rights, args)
            evidence[pos] = ev
            if ev['crossing_pass']:
                cuts3.append(pos)
        for method, cuts in [('V1', cuts1), ('V2', cuts2), ('V3', cuts3)]:
            for start, end in split_intervals(territory, cuts):
                results[method]['crus'].append({'chrom': chrom, 'start': start, 'end': end, 'territory_id': tid})
            for pos in sorted(cuts):
                results[method]['boundaries'].append({**cut_info[pos], **evidence[pos], 'method': method})
            if method != 'V1':
                results[method]['candidates'].extend([{**r} for r in results['V1']['candidates'] if r['territory_id'] == tid] + [{**r} for r in extra])
                results[method]['backgrounds'].extend(bgrows)
                results[method]['actions'].extend(actions)
        for r in results['V3']['candidates']:
            pos = r['boundary_pos']
            if r['territory_id'] == tid and r['accepted']:
                r.update(evidence[pos])
                if not evidence[pos]['crossing_pass']:
                    r['accepted'], r['reject_reason'] = 0, 'global_crossing_exceeds_weak_internal_side'
    for method, result in results.items():
        for i, row in enumerate(sorted(result['crus'], key=lambda r: r['start']), 1):
            row['cru_id'] = f'{chrom}_{method}_CRU_{i}'
        result['territories'] = territories
    return results


def interval_summary(method, chrom, crus, edges, size):
    lengths = np.array([r['end'] - r['start'] for r in crus])
    positions = sorted((r['start'], r['end']) for r in crus)
    starts = [r[0] for r in positions]
    if any(b[0] < a[1] for a, b in zip(positions, positions[1:])):
        raise ValueError(f'Overlapping CRUs in {method}:{chrom}')
    edge_class = Counter()
    for e in edges:
        ids = []
        for pos in [e['p1'], e['p2']]:
            i = bisect.bisect_right(starts, pos) - 1
            ids.append(i if i >= 0 and positions[i][0] <= pos < positions[i][1] else -1)
        key = 'intra' if ids[0] == ids[1] and ids[0] >= 0 else 'inter' if min(ids) >= 0 else 'one_endpoint' if max(ids) >= 0 else 'outside'
        edge_class[key] += 1
    out = {'method': method, 'chrom': chrom, 'chromosome_bp': size, 'cru_n': len(crus),
           'covered_bp': int(lengths.sum()), 'coverage_fraction': float(lengths.sum()) / size,
           'aglink_n': len(edges), 'intra_aglink_n': edge_class['intra'],
           'intra_aglink_fraction': edge_class['intra'] / len(edges) if edges else 0,
           'inter_aglink_n': edge_class['inter'], 'one_endpoint_aglink_n': edge_class['one_endpoint'],
           'outside_aglink_n': edge_class['outside'], 'length_ge2mb_n': int(sum(lengths >= 2_000_000)),
           'length_ge5mb_n': int(sum(lengths >= 5_000_000))}
    for label, q in [('p25', .25), ('median', .5), ('p75', .75), ('p90', .9), ('p95', .95)]:
        out[label + '_bp'] = float(np.quantile(lengths, q)) if len(lengths) else ''
    out['mean_bp'] = float(np.mean(lengths)) if len(lengths) else ''
    out['max_bp'] = int(max(lengths)) if len(lengths) else ''
    return out


def boundary_sets(result, sizes):
    allpts = {(r['chrom'], pos) for r in result['crus'] for pos in [r['start'], r['end']] if 0 < pos < sizes[r['chrom']]}
    internal = {(r['chrom'], r['boundary_pos']) for r in result['boundaries']} & allpts
    return {'all_final_cru_boundaries': allpts, 'internal_boundaries': internal, 'territory_outer_edges': allpts - internal}


def block_evaluation(results, sizes, chroms, args, helper, comparison):
    blocks = helper['load_block_bed'](args.block_bed, chroms)
    index = helper['block_test_index'](blocks, chroms)
    metrics, observed, nearest, permutation_detail = [], [], [], []
    directions = {f'within_{k}kb_prop': 'higher_is_better' for k in [25, 50, 100, 200]}
    directions.update({'median_distance_bp': 'lower_is_better', 'mean_distance_bp': 'lower_is_better'})
    for method, result in results.items():
        sets = boundary_sets(result, sizes)
        for source, pointset in sets.items():
            for scope in chroms + ['Tested_chromosomes']:
                points = [{'chrom': c, 'pos': p} for c, p in sorted(pointset) if scope == 'Tested_chromosomes' or c == scope]
                if not points:
                    observed.append({'method': method, 'boundary_source': source, 'chrom': scope, 'boundary_n': 0})
                    continue
                obs = helper['boundary_block_metrics'](points, index)
                observed.append({'method': method, 'boundary_source': source, 'chrom': scope, **obs})
                if scope != 'Tested_chromosomes':
                    nearest.extend({'method': method, 'boundary_source': source, 'chrom': r['chrom'], 'position': r['pos'],
                                    'nearest_block_boundary_distance_bp': helper['nearest_block_distance'](r['pos'], index[r['chrom']])} for r in points)
                for null_name, func in [('local_boundary_points_10Mb', helper['local_random_points']),
                                        ('chromosome_100bin_boundary_points', helper['chr100bin_random_points'])]:
                    rng = random.Random(args.block_seed)
                    nullrows = []
                    for rep in range(args.block_perm_n):
                        randomized = func(points, sizes, args.block_local_window_bp, rng) if null_name.startswith('local_') else func(points, sizes, rng)
                        vals = helper['boundary_block_metrics'](randomized, index)
                        nullrows.append(vals)
                        permutation_detail.append({'method': method, 'boundary_source': source, 'chrom': scope,
                                                   'null_model': null_name, 'replicate': rep + 1, **vals})
                    rows = helper['summarize_permutation'](method, null_name, obs, nullrows, directions)
                    for row in rows:
                        row.update({'method': method, 'boundary_source': source, 'chrom': scope,
                                    'boundary_n': len(points), 'block_seed': args.block_seed})
                    metrics.extend(rows)
        log(f'{method}: block-boundary permutations complete')
    write_tsv(comparison / '02_boundary_block_observed.tsv', observed)
    write_tsv(comparison / '03_boundary_block_permutation_summary.tsv', metrics)
    write_tsv(comparison / '04_boundary_block_nearest_distance.tsv', nearest)
    write_tsv(comparison / '05_boundary_permutation_detail.tsv.gz', permutation_detail)
    for method in results:
        write_tsv(args.outdir / METHODS[method] / '06_boundary_block_permutation_summary.tsv', [r for r in metrics if r['method'] == method])
    return blocks, metrics


def plot_region(method, result, edges, sig, blocks, region, args, helper, outdir):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    chrom, start, end = region
    plt.rcParams.update({'font.family': args.font_family, 'font.size': 10, 'pdf.fonttype': 42, 'ps.fonttype': 42})
    fig, axes = plt.subplots(7, 1, figsize=(12, 9), sharex=True,
                             gridspec_kw={'height_ratios': [.7, 1.4, 1.25, .9, .7, .55, .45], 'hspace': .10})
    visible = [e for e in edges if e['right'] > start and e['left'] < end]
    helper['draw_oldstyle_arc_track'](axes[0], visible, start, end, '#2E8DC8', .28, .7, arc_points=36)
    axes[0].set_ylabel('agLINK', fontsize=10)
    axes[0].set_title(f'{method}   {chrom}:{start / 1e6:g}-{end / 1e6:g} Mb   ({len(visible):,} links)', fontsize=12, loc='left')
    idx = np.arange(len(sig['raw']))
    x = (idx + .5) * args.bin_bp
    keep = (x >= start) & (x <= end)
    axes[1].plot(x[keep], sig['raw'][keep], color='#C1C0BF', lw=.65, label='Raw span coverage')
    for y, color, width in zip(sig['smooth'], ['#96CEDD', '#5FB6BD', '#E89B6C'], sig['widths']):
        axes[1].plot(x[keep], y[keep], color=color, lw=1.1, label=f'{width * args.bin_bp / 1000:g}-kb mean')
    axes[1].axhline(sig['threshold'], color='#777777', ls='--', lw=.75, label='Activity threshold')
    axes[1].set_ylabel('Span\ncoverage')
    axes[1].legend(loc='upper left', ncol=3, fontsize=8, frameon=False)
    y = sig['smooth'][1]
    axes[2].plot(x[keep], y[keep], color='#5FB6BD', lw=1.1, label='Medium-scale density')
    for bg in result.get('backgrounds', []):
        if bg['chrom'] == chrom and bg['parent_end'] > start and bg['parent_start'] < end and bg['local_background'] != '':
            axes[2].hlines(bg['local_background'], max(start, bg['parent_start']), min(end, bg['parent_end']), color='#777777', ls='--', lw=.9)
            axes[2].axvspan(bg['background_region_start'], bg['background_region_end'], color='#EEC566', alpha=.15)
    candidates = [c for c in result.get('candidates', []) if c['chrom'] == chrom and start <= c['boundary_pos'] <= end]
    for pos in sorted({c['boundary_pos'] for c in candidates}):
        axes[2].plot(pos, np.interp(pos, x, y), marker='v', ms=4, color='#B5B5B5', linestyle='None')
    axes[2].set_ylabel('Local background\n/ valleys')
    axes[2].legend(loc='upper left', fontsize=8, frameon=False)
    axes[3].plot(x[keep], sig['endpoint_broad'][keep], color='#2E8DC8', lw=1.1, label='225-kb endpoint count')
    axes[3].set_ylabel('Endpoint\nactivity')
    axes[3].legend(loc='upper left', fontsize=8, frameon=False)
    endpoints = {(r['chrom'], p) for r in result['crus'] for p in [r['start'], r['end']]}
    internal = sorted({r['boundary_pos'] for r in result['boundaries'] if r['chrom'] == chrom
                       and (chrom, r['boundary_pos']) in endpoints and start <= r['boundary_pos'] <= end})
    # Build visible final endpoints without allowing another chromosome's sizes
    # to leak into this region's boundary set.
    internal_all = {(r['chrom'], r['boundary_pos']) for r in result['boundaries']}
    outer = {(r['chrom'], p) for r in result['crus'] for p in [r['start'], r['end']] if (r['chrom'], p) not in internal_all}
    for pos in internal:
        axes[4].vlines(pos, .25, .65, color=COLORS[method], lw=1.4)
        axes[1].axvline(pos, color=COLORS[method], ls=':', lw=.85, alpha=.7)
    for c, pos in outer:
        if c == chrom and start <= pos <= end:
            axes[4].vlines(pos, .25, .65, color='#AAAAAA', ls='--', lw=1)
    axes[4].set_ylim(0, 1)
    axes[4].set_yticks([])
    axes[4].set_ylabel('CRU\nboundary')
    from matplotlib.lines import Line2D
    axes[4].legend(handles=[Line2D([0], [0], color=COLORS[method], label='Internal cut'),
                            Line2D([0], [0], color='#AAAAAA', ls='--', label='Outer edge')],
                   loc='upper right', fontsize=8, ncol=2, frameon=False)
    for i, row in enumerate(r for r in result['crus'] if r['chrom'] == chrom and r['start'] < end and r['end'] > start):
        lo, hi = max(start, row['start']), min(end, row['end'])
        axes[5].add_patch(Rectangle((lo, .22), hi - lo, .56, facecolor=['#3F8FC1', '#92BCD8'][i % 2], edgecolor='#555555', lw=.85))
    axes[5].set_ylabel(method + '\nCRU')
    axes[5].set_ylim(0, 1)
    axes[5].set_yticks([])
    for i, row in enumerate(blocks.get(chrom, [])):
        if row['start'] < end and row['end'] > start:
            lo, hi = max(start, row['start']), min(end, row['end'])
            axes[6].add_patch(Rectangle((lo, 0), hi - lo, 1, facecolor=['#E6E6E6', '#CDCDCD'][i % 2], edgecolor='#777777', lw=.8))
    axes[6].set_ylim(0, 1)
    axes[6].set_yticks([])
    axes[6].set_ylabel('Block')
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    for i, ax in enumerate(axes):
        ax.set_xlim(start, end)
        ax.grid(False)
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_color('#222222')
            spine.set_linewidth(.9)
        ax.tick_params(width=.8, labelsize=9)
        if i in [1, 2, 3]:
            low, high = ax.get_ylim()
            ymax = max(high, 1)
            ax.set_ylim(-.07 * ymax, 1.35 * ymax)
        if i < 6:
            ax.tick_params(axis='x', bottom=False, labelbottom=False)
    axes[6].xaxis.set_major_locator(MaxNLocator(7))
    axes[6].xaxis.set_major_formatter(FuncFormatter(lambda v, _: f'{v / 1e6:g}'))
    axes[6].set_xlabel(f'{chrom} position (Mb)')
    fig.subplots_adjust(left=.11, right=.99, top=.96, bottom=.07)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = outdir / f'{chrom}_{start}_{end}'
    fig.savefig(str(stem) + '.pdf')
    fig.savefig(str(stem) + '.png', dpi=160)
    plt.close(fig)


def plot_comparison(results, edges, sig, blocks, region, args, helper, outdir):
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.ticker import FuncFormatter, MaxNLocator
    chrom, start, end = region
    plt.rcParams.update({'font.family': args.font_family, 'font.size': 10, 'pdf.fonttype': 42})
    fig, axes = plt.subplots(7, 1, figsize=(12, 7.8), sharex=True,
                             gridspec_kw={'height_ratios': [.8, 1.4, .6, .6, .6, .6, .4], 'hspace': .12})
    visible = [e for e in edges if e['right'] > start and e['left'] < end]
    helper['draw_oldstyle_arc_track'](axes[0], visible, start, end, '#2E8DC8', .28, .7, arc_points=36)
    axes[0].set_ylabel('agLINK')
    axes[0].set_title(f'{chrom}:{start / 1e6:g}-{end / 1e6:g} Mb   |   independent agLINK-only pilot vs saved M3', loc='left', fontsize=12)
    x = (np.arange(len(sig['raw'])) + .5) * args.bin_bp
    keep = (x >= start) & (x <= end)
    for y, color, width in zip(sig['smooth'], ['#96CEDD', '#5FB6BD', '#E89B6C'], sig['widths']):
        axes[1].plot(x[keep], y[keep], color=color, lw=1.15, label=f'{width * args.bin_bp / 1000:g}-kb mean')
    axes[1].set_ylabel('Span\ncoverage')
    axes[1].legend(loc='upper left', ncol=3, fontsize=8, frameon=False)
    high = max(axes[1].get_ylim()[1], 1)
    axes[1].set_ylim(-.06 * high, high * 1.3)
    for ax, (method, result) in zip(axes[2:6], results.items()):
        cuts = {r['boundary_pos'] for r in result['boundaries'] if r['chrom'] == chrom}
        for row in result['crus']:
            if row['chrom'] == chrom and row['start'] < end and row['end'] > start:
                lo, hi = max(start, row['start']), min(end, row['end'])
                ax.add_patch(Rectangle((lo, .18), hi - lo, .42, facecolor=COLORS[method], edgecolor='#555555', lw=.85))
                for pos in [row['start'], row['end']]:
                    if start <= pos <= end:
                        ax.vlines(pos, .65, .91, color=COLORS[method] if pos in cuts else '#999999',
                                  linestyle='-' if pos in cuts else '--', lw=1.2)
        ax.set_ylabel(method)
        ax.set_ylim(0, 1)
        ax.set_yticks([])
    for i, row in enumerate(blocks.get(chrom, [])):
        if row['start'] < end and row['end'] > start:
            lo, hi = max(start, row['start']), min(end, row['end'])
            axes[6].add_patch(Rectangle((lo, 0), hi - lo, 1, facecolor=['#E6E6E6', '#CDCDCD'][i % 2], edgecolor='#777777', lw=.8))
    axes[6].set_ylabel('Block')
    axes[6].set_ylim(0, 1)
    axes[6].set_yticks([])
    for i, ax in enumerate(axes):
        ax.set_xlim(start, end)
        for spine in ax.spines.values():
            spine.set_linewidth(.9)
            spine.set_visible(True)
            spine.set_color('#222222')
        ax.tick_params(width=.8, labelsize=9)
        if i < 6:
            ax.tick_params(axis='x', bottom=False, labelbottom=False)
    axes[6].xaxis.set_major_locator(MaxNLocator(7))
    axes[6].xaxis.set_major_formatter(FuncFormatter(lambda value, _: f'{value / 1e6:g}'))
    axes[6].set_xlabel(f'{chrom} position (Mb)')
    fig.subplots_adjust(left=.09, right=.99, top=.95, bottom=.08)
    outdir.mkdir(parents=True, exist_ok=True)
    stem = outdir / f'{chrom}_{start}_{end}.comparison'
    fig.savefig(str(stem) + '.pdf')
    fig.savefig(str(stem) + '.png', dpi=160)
    plt.close(fig)


def parser():
    ap = argparse.ArgumentParser(description=__doc__)
    for flag in ['aglink', 'gene-bed', 'acr-bed', 'genome-sizes', 'block-bed', 'm3-dir', 'm3-script', 'outdir']:
        ap.add_argument('--' + flag, type=Path, required=True)
    ap.add_argument('--bedpe-reference', type=Path)
    ap.add_argument('--chroms', default='1A,1B,2A')
    ap.add_argument('--regions', default='1A:34000000-36000000,2A:15000000-22000000,1A:538500000-542000000,2A:765000000-771000000,1B:164375000-166400000')
    ap.add_argument('--bin-bp', type=int, default=5000)
    ap.add_argument('--smooth-bp', default='25000,105000,225000')
    ap.add_argument('--active-quantile', type=float, default=.25)
    ap.add_argument('--min-span-count', type=float, default=2)
    ap.add_argument('--min-endpoint-count', type=float, default=2)
    ap.add_argument('--min-cru-bp', type=int, default=100000)
    ap.add_argument('--min-internal', type=int, default=2)
    ap.add_argument('--min-prominence', type=float, default=1)
    ap.add_argument('--valley-min-width-bp', type=int, default=10000)
    ap.add_argument('--shoulder-bp', type=int, default=500000)
    ap.add_argument('--min-peak-balance', type=float, default=.2)
    ap.add_argument('--v1-depth', type=float, default=.35)
    ap.add_argument('--long-region-bp', type=int, default=1500000)
    ap.add_argument('--background-regions', type=int, default=5)
    ap.add_argument('--min-background-regions', type=int, default=2)
    ap.add_argument('--background-sep-bp', type=int, default=100000)
    ap.add_argument('--v2-background-depth', type=float, default=.55)
    ap.add_argument('--max-recursion-depth', type=int, default=16)
    ap.add_argument('--cross-flank-bp', type=int, default=250000)
    ap.add_argument('--max-cross-fraction', type=float, default=.5)
    ap.add_argument('--block-perm-n', type=int, default=1000)
    ap.add_argument('--block-seed', type=int, default=20260805)
    ap.add_argument('--block-local-window-bp', type=int, default=10000000)
    ap.add_argument('--font-family', default='Liberation Sans')
    ap.add_argument('--skip-plots', action='store_true')
    return ap


def main():
    args = parser().parse_args()
    args.chroms = args.chroms.split(',')
    if len(set(args.chroms)) != len(args.chroms) or len(args.chroms) > 3 or 'all' in args.chroms:
        raise ValueError('This pilot is restricted to at most three explicitly named chromosomes.')
    args.smooth_bp = [int(x) for x in args.smooth_bp.split(',')]
    if len(args.smooth_bp) != 3 or args.smooth_bp != sorted(args.smooth_bp):
        raise ValueError('Supply exactly three increasing smoothing widths.')
    args.sizes = dict((c, int(n)) for c, n in (line.split()[:2] for line in open(args.genome_sizes) if line.strip()))
    if (args.outdir / '01_parameters.json').exists():
        raise ValueError('Output run exists. Choose an unused run directory; no historical files are overwritten.')
    args.outdir.mkdir(parents=True, exist_ok=True)
    payload = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()}
    payload['method_notes'] = {'pcc': 'not read, filtered or weighted; upstream significant agLINK input reused',
                               'coacr': 'not loaded', 'm3_boundaries': 'not used in calling',
                               'span': 'ACR midpoint to strand-aware gene TSS, inclusive integer positions represented as BED [left,right+1)',
                               'block_null': 'Reuse original M3 randomization and metric functions; reset seed per method/source/scope/null, not bit-identical historical random draws.'}
    (args.outdir / '01_parameters.json').write_text(json.dumps(payload, indent=2) + '\n')
    shared = args.outdir / 'f01_shared_input'
    shared.mkdir()
    provenance = []
    for key in ['aglink', 'gene_bed', 'acr_bed', 'genome_sizes', 'bedpe_reference', 'block_bed', 'm3_script']:
        path = getattr(args, key)
        if path:
            provenance.append({'role': key, 'path': str(path), 'sha256': sha256(path)})
    for filename in ['M3_CRU.bed', 'M3_primary_boundaries.bed', 'M3_secondary_boundaries.bed', 'M3_parameters.tsv']:
        path = args.m3_dir / filename
        provenance.append({'role': 'M3_reference', 'path': str(path), 'sha256': sha256(path)})
    write_tsv(shared / '03_input_provenance.tsv', provenance)
    edges_by_chrom = load_count_only_edges(args, args.sizes, shared)
    results = {method: {'crus': [], 'boundaries': [], 'candidates': [], 'backgrounds': [], 'actions': [], 'territories': []} for method in ['V1', 'V2', 'V3']}
    signals_by_chrom, influence = {}, []
    for chrom in args.chroms:
        edges = edges_by_chrom[chrom]
        log(f'{chrom}: {len(edges):,} equal-weight agLINKs; no PCC/coACR/M3 boundary inputs to construction')
        sig = signals(edges, args.sizes[chrom], args)
        signals_by_chrom[chrom] = sig
        called = build_methods(chrom, args.sizes[chrom], edges, sig, args)
        for method in results:
            for key in results[method]:
                results[method][key].extend(called[method][key])
        write_tsv(shared / f'{chrom}.profiles.tsv.gz', ({'chrom': chrom, 'start': i * args.bin_bp,
            'end': min(args.sizes[chrom], (i + 1) * args.bin_bp), 'span_coverage': sig['raw'][i],
            'smooth_fine': sig['smooth'][0][i], 'smooth_medium': sig['smooth'][1][i], 'smooth_broad': sig['smooth'][2][i],
            'endpoint_n': sig['endpoint'][i], 'endpoint_broad': sig['endpoint_broad'][i],
            'long_ge1mb_coverage': sig['long_ge1mb'][i], 'active': int(sig['active'][i]),
            'activity_threshold': sig['threshold']} for i in range(len(sig['raw']))))
        lens = np.array([e['distance_bp'] for e in edges])
        influence.append({'chrom': chrom, 'aglink_n': len(edges), 'ge1mb_aglink_n': int(sum(lens >= 1e6)),
                          'ge1mb_link_fraction': float(np.mean(lens >= 1e6)),
                          'ge1mb_span_mass_fraction': float(sum(lens[lens >= 1e6] + 1) / sum(lens + 1)),
                          'single_link_only_bin_n': int(sum((sig['raw'] > 0) & (sig['raw'] <= 1))),
                          'active_bin_n': int(sum(sig['active'])), 'activity_threshold': sig['threshold']})
        log(f'{chrom}: ' + '; '.join(f'{m}={len(called[m]["crus"])} CRUs/{len(called[m]["boundaries"])} internal cuts' for m in called))
    write_tsv(shared / '04_long_aglink_influence.tsv', influence)
    # Calling is now finished. M3 intervals and blocks are first exposed below.
    m3bed = read_bed(args.m3_dir / 'M3_CRU.bed', args.chroms)
    m3cuts = []
    for filename in ['M3_primary_boundaries.bed', 'M3_secondary_boundaries.bed']:
        for rows in read_bed(args.m3_dir / filename, args.chroms).values():
            m3cuts.extend({**r, 'boundary_pos': (r['start'] + r['end']) // 2} for r in rows)
    results['M3'] = {'crus': [r for c in args.chroms for r in m3bed[c]], 'boundaries': m3cuts,
                     'candidates': [], 'backgrounds': [], 'territories': [], 'actions': []}
    comparison = args.outdir / 'f06_comparison'
    comparison.mkdir()
    summary = []
    lengthrows = []
    for method, result in results.items():
        method_dir = args.outdir / METHODS[method]
        method_dir.mkdir()
        (method_dir / '01_parameters.json').write_text(json.dumps({**payload, 'method': method}, indent=2) + '\n')
        with (method_dir / '02_CRU.bed').open('w') as handle:
            for r in result['crus']:
                handle.write(f'{r["chrom"]}\t{r["start"]}\t{r["end"]}\t{r["cru_id"]}\n')
        write_tsv(method_dir / '03_internal_boundaries.tsv', result['boundaries'], None if result['boundaries'] else ['chrom', 'boundary_pos'])
        write_tsv(method_dir / '04_candidate_valleys.tsv', result['candidates'], None if result['candidates'] else ['chrom', 'boundary_pos', 'accepted'])
        write_tsv(method_dir / '07_local_backgrounds.tsv', result['backgrounds'], None if result['backgrounds'] else ['chrom', 'parent_start', 'parent_end'])
        write_tsv(method_dir / '08_recursive_actions.tsv', result['actions'], None if result['actions'] else ['chrom', 'action'])
        method_summary = []
        for chrom in args.chroms:
            method_summary.append(interval_summary(method, chrom, [r for r in result['crus'] if r['chrom'] == chrom], edges_by_chrom[chrom], args.sizes[chrom]))
        # Quantiles must be pooled over actual intervals, never averaged by chr.
        pooled = interval_summary(method, 'Tested_chromosomes', [], [], sum(args.sizes[c] for c in args.chroms))
        lengths = np.array([r['end'] - r['start'] for r in result['crus']])
        for key in ['cru_n', 'covered_bp', 'aglink_n', 'intra_aglink_n', 'inter_aglink_n', 'one_endpoint_aglink_n', 'outside_aglink_n', 'length_ge2mb_n', 'length_ge5mb_n']:
            pooled[key] = sum(r[key] for r in method_summary)
        pooled['coverage_fraction'] = pooled['covered_bp'] / pooled['chromosome_bp']
        pooled['intra_aglink_fraction'] = pooled['intra_aglink_n'] / pooled['aglink_n']
        for label, q in [('p25', .25), ('median', .5), ('p75', .75), ('p90', .9), ('p95', .95)]:
            pooled[label + '_bp'] = float(np.quantile(lengths, q)) if len(lengths) else ''
        pooled['mean_bp'], pooled['max_bp'] = (float(np.mean(lengths)), int(max(lengths))) if len(lengths) else ('', '')
        method_summary.append(pooled)
        write_tsv(method_dir / '05_summary.tsv', method_summary)
        summary.extend(method_summary)
        lengthrows.extend({'method': method, **r, 'length_bp': r['end'] - r['start']} for r in result['crus'])
    write_tsv(comparison / '01_CRU_summary.tsv', summary)
    write_tsv(comparison / '06_CRU_lengths.tsv', lengthrows)
    write_tsv(comparison / '07_long_CRUs_ge2mb.tsv', [r for r in lengthrows if r['length_bp'] >= 2_000_000])
    helper = runpy.run_path(str(args.m3_script), run_name='m3_diagnostics_and_arcs_only')
    blocks, metrics = block_evaluation(results, args.sizes, args.chroms, args, helper, comparison)
    regions = []
    for spec in args.regions.split(','):
        c, bounds = spec.split(':')
        start, end = map(int, bounds.split('-'))
        if c not in args.chroms or not 0 <= start < end <= args.sizes[c]:
            raise ValueError(f'Invalid plot region {spec}')
        regions.append((c, start, end))
    write_tsv(comparison / '08_plot_regions.tsv', [{'chrom': c, 'start': a, 'end': b} for c, a, b in regions])
    if not args.skip_plots:
        for region in regions:
            for method, result in results.items():
                plot_region(method, result, edges_by_chrom[region[0]], signals_by_chrom[region[0]], blocks, region, args, helper,
                            args.outdir / METHODS[method] / 'f01_publication_tracks')
            plot_comparison(results, edges_by_chrom[region[0]], signals_by_chrom[region[0]], blocks, region, args, helper,
                            comparison / 'f01_region_comparisons')
            log(f'Regional PDF/PNG rendered: {region}')
    frozen = args.outdir / 'f07_scripts'
    frozen.mkdir()
    for filename in ['01_call_span_density_pilot.py', '02_test_span_density_pilot.py', '03_plot_pilot_summary.R', '04_METHODS.md']:
        source = Path(__file__).with_name(filename)
        if source.exists():
            shutil.copy2(source, frozen / filename)
    write_tsv(comparison / '09_source_integrity_after_run.tsv', [{**r, 'sha256_after': sha256(r['path']),
               'unchanged': int(sha256(r['path']) == r['sha256'])} for r in provenance])
    log('TESTED CHROMOSOMES ONLY: ' + ','.join(args.chroms))
    for row in summary:
        if row['chrom'] == 'Tested_chromosomes':
            log(f'{row["method"]}: CRUs={row["cru_n"]}; coverage={row["coverage_fraction"]:.2%}; median={row["median_bp"]} bp; intra-agLINK={row["intra_aglink_fraction"]:.2%}')
    for row in metrics:
        if row['chrom'] == 'Tested_chromosomes' and row['boundary_source'] == 'all_final_cru_boundaries' and row['null_model'] == 'local_boundary_points_10Mb' and row['metric'] == 'within_25kb_prop':
            log(f'{row["method"]}: all-final <=25 kb observed={row["observed"]:.2%}, random={row["permutation_mean"]:.2%}, fold={row["observed_over_permutation"]:.3f}, P={row["empirical_p_ge"]:.4g}')
    (args.outdir / '02_ANALYSIS.DONE').write_text('Pilot calling, post-hoc evaluation and requested regional plots completed.\n')


if __name__ == '__main__':
    main()
