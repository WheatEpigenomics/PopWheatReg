# m07_changepoint: V6/V7 change point
# Frozen rule module, copied unchanged from the published calling chain
# (source: f11_B145_V6_V7_changepoint/01_call_V6_V7.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Frozen V5 partitions plus generic density/connection and local change-point tests."""
import argparse
from collections import Counter, defaultdict
import importlib.util
import json
from pathlib import Path
import shutil
from types import SimpleNamespace

import numpy as np
from scipy.signal import find_peaks


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


from . import m06_depth_sensitivity as A
F, C, P = A.F, A.C, A.P
STEP, CHROMS = A.STEP, A.CHROMS
METHODS = {'Baseline': 'f01_baseline', 'V5_reference': 'f02_V5_reference',
           'V6': 'f03_V6', 'V7': 'f04_V7'}
REGIONS = A.REGIONS
NEW_SOURCES = {'sharp', 'step', 'changepoint'}


def config(args):
    old = json.loads((args.reference_dir / '01_test_settings.json').read_text())
    opts = A.tier(SimpleNamespace(**old), 'V5')
    opts.sharp_depth = args.sharp_depth
    opts.medium_depth = args.medium_depth
    opts.sharp_max_width_bp = args.sharp_max_width_bp
    opts.step_relative = args.step_relative
    opts.max_cross_fraction = args.max_cross_fraction
    opts.min_daughter_bp = args.min_daughter_bp
    opts.enable_strong = False
    opts.cp_candidate = args.cp_candidate
    opts.cp_accept = args.cp_accept
    opts.cp_sse_gain = args.cp_sse_gain
    return opts


def finite(value, default=np.nan):
    try:
        return float(value)
    except (ValueError, TypeError):
        return default


def cp_proposals(sig, start, end, opts):
    """A local two-level contrast scan; no boundary or maximum-length prior."""
    y = sig['smooth'][1]
    a, b = F.bounds(start, end, len(y))
    pref = np.r_[0., np.cumsum(y)]
    grad = np.gradient(sig['smooth'][0])
    out = set()
    for flank in [200000, 400000, 800000]:
        w, skip = flank // STEP, 50000 // STEP
        ii = np.arange(a + w + skip, b - w - skip - 1)
        if len(ii) < 3:
            continue
        left = (pref[ii - skip] - pref[ii - skip - w]) / w
        right = (pref[ii + skip + 1 + w] - pref[ii + skip + 1]) / w
        rel = np.abs(right - left) / np.maximum(np.maximum(left, right), 1.)
        peaks, _ = find_peaks(rel, height=opts.cp_candidate, prominence=.03, distance=10)
        for j in peaks:
            i = int(ii[j]); direction = 1 if right[j] > left[j] else -1
            lo, hi = max(a + 1, i - 15), min(b - 1, i + 16)
            if hi > lo:
                refined = lo + int(np.argmax(direction * grad[lo:hi]))
                out.add(refined)
    return [('changepoint', i) for i in sorted(out)]


def cp_morphology(i, sig, start, end, opts):
    a, b = F.bounds(start, end, len(sig['raw']))
    tested = []
    for flank in [200000, 400000, 800000]:
        w, skip = flank // STEP, 50000 // STEP
        if i - skip - w < a or i + skip + 1 + w > b:
            continue
        row = {'source': 'changepoint', 'boundary_pos': i * STEP + STEP // 2,
               'parent_start': start, 'parent_end': end, 'cp_flank_bp': flank,
               'local_background': float(sig['background'][i]), 'raw_span': float(sig['raw'][i])}
        reason = []
        segments = []
        signs = []
        for k, label in enumerate(['25k', '105k', '225k']):
            y = sig['smooth'][k]
            ls, rs = y[i-skip-w:i-skip], y[i+skip+1:i+skip+1+w]
            lv, rv = float(np.median(ls)), float(np.median(rs))
            signs.append(1 if rv > lv else -1)
            row.update({label+'_left_level': lv, label+'_right_level': rv,
                        label+'_relative_change': abs(rv-lv)/max(lv, rv, 1.)})
            segments.append((ls, rs))
        direction = signs[1]; row['direction'] = 'rise' if direction > 0 else 'drop'
        lv, rv = row['105k_left_level'], row['105k_right_level']; high = max(lv, rv, 1.)
        ls, rs = segments[1]
        lb = np.array([np.median(z) for z in np.array_split(ls, 4)])
        rb = np.array([np.median(z) for z in np.array_split(rs, 4)])
        sustain = int(np.count_nonzero(direction*(rb-lb) >= .20*high))
        spread = max(float(np.ptp(lb)), float(np.ptp(rb))) / high
        values = np.r_[ls, rs]
        sse0 = float(np.sum((values-np.mean(values))**2))
        sse1 = float(np.sum((ls-np.mean(ls))**2)+np.sum((rs-np.mean(rs))**2))
        gain = 1-sse1/sse0 if sse0 > 0 else 0.
        delta = abs(row['25k_left_level']-row['25k_right_level'])
        yl = sig['smooth'][0]; low25 = min(row['25k_left_level'], row['25k_right_level'])
        high25 = max(row['25k_left_level'], row['25k_right_level'])
        pre = np.arange(max(a, i-40), i+1); post = np.arange(i, min(b, i+41))
        if direction > 0:
            pl = pre[yl[pre] <= low25+.2*delta]; pr = post[yl[post] >= high25-.2*delta]
        else:
            pl = pre[yl[pre] >= high25-.2*delta]; pr = post[yl[post] <= low25+.2*delta]
        width = max(STEP, int(pr[0]-pl[-1])*STEP) if len(pl) and len(pr) else 0
        if row['105k_relative_change'] < opts.cp_accept: reason.append('cp_105k_relative_change')
        if row['25k_relative_change'] < opts.cp_candidate or signs[0] != direction: reason.append('cp_fine_direction_or_magnitude')
        if sustain < 3: reason.append('cp_not_sustained')
        if spread > .40: reason.append('cp_unstable_platforms')
        if gain < opts.cp_sse_gain: reason.append('cp_small_two_level_SSE_gain')
        if delta < opts.min_absolute_change: reason.append('small_absolute_change')
        if not STEP <= width <= 200000: reason.append('cp_transition_width')
        row.update(transition_width_bp=width, absolute_change=delta, sustained_block_n=sustain,
                   persistence_blocks=4, cp_two_level_SSE_gain=gain, cp_platform_spread=spread,
                   scale_support_n=2+int(signs[2] == direction and row['225k_relative_change'] >= .10),
                   shape_pass=int(not reason), shape_reason=','.join(reason),
                   strength=row['105k_relative_change']*(2+int(gain >= .5)))
        tested.append(row)
    if not tested:
        return {'source': 'changepoint', 'boundary_pos': i*STEP+STEP//2,
                'parent_start': start, 'parent_end': end, 'shape_pass': 0,
                'shape_reason': 'cp_insufficient_flanks', 'strength': 0, 'scale_support_n': 0}
    return max(tested, key=lambda r: (r['shape_pass'], r['cp_two_level_SSE_gain']*r['105k_relative_change'], -r['cp_flank_bp']))


def evaluate(source, i, sig, parent, es, opts, acrs, genes, cache):
    row = cp_morphology(i, sig, parent['start'], parent['end'], opts) if source == 'changepoint' else A.morphology(source, i, sig, parent['start'], parent['end'], opts)
    row.update(F.graph(es, parent['start'], parent['end'], row['boundary_pos'], opts))
    row.update(chrom=parent['chrom'], baseline_id=parent['baseline_id'], v5_parent_id=parent['v5_parent_id'])
    key = (parent['start'], parent['end'], i)
    if key not in cache:
        cache[key] = [A.short_evidence(es, parent['start'], parent['end'], row['boundary_pos'], cap, acrs, genes) for cap in [500000, 1000000]]
    short = cache[key]
    for d in short:
        row.update({f'short_{d["cap_bp"]}_{k}': v for k, v in d.items()})
    hard = []
    if min(row['boundary_pos']-parent['start'], parent['end']-row['boundary_pos']) < opts.min_daughter_bp: hard.append('minimum_daughter_500kb')
    if min(row['whole_left_internal'], row['whole_right_internal']) < 3: hard.append('weak_daughter_internal')
    if min(row['whole_left_ACR_n'], row['whole_right_ACR_n']) < 2: hard.append('daughter_ACR_diversity')
    if row['global_cut_loss'] > .40: hard.append('global_cut_loss')
    contexts = []
    for flank in [250000, 500000]:
        prefix = f'context_{flank}_'
        ctx = {k[len(prefix):]: v for k, v in row.items() if k.startswith(prefix)}
        if min(ctx['left_internal'], ctx['right_internal']) >= 2 and min(ctx['left_ACR_n'], ctx['right_ACR_n']) >= 2:
            contexts.append(ctx)
    if not contexts: hard.append('necessary_local_module_support')
    direct = [q for q in contexts if q['local_cross_fraction'] <= opts.max_cross_fraction]
    rescues = []
    for d in short:
        enough = min(d['left_internal'], d['right_internal']) >= 3 and min(d['left_ACR_n'], d['right_ACR_n']) >= 2
        direction_row = {**row, 'source': 'step' if source == 'changepoint' else source}
        pc = A.profile_change(d, direction_row)
        D = d['relative_depletion']
        graph_vote = np.isfinite(D) and D >= .10
        profile_vote = np.isfinite(pc) and pc >= .20
        separation = d['cross_fraction'] <= .75
        # Graph/profile are correlated diagnostics of the same edges, not independent validation.
        strong_shape = row.get('25k_relative_change', 0) >= .50 and row.get('105k_relative_change', 0) >= .15
        supported = (graph_vote and profile_vote) or (graph_vote and separation) or (strong_shape and profile_vote and d['cross_fraction'] <= .65)
        if enough and supported:
            rescues.append((d, pc, int(graph_vote)+int(profile_vote)+int(separation)))
    route = 'rejected'; rescue = None
    if row['shape_pass'] and not hard:
        if direct:
            row.update(min(direct, key=lambda q: (q['support_context_bp'], q['local_cross_fraction'])))
            route = 'relaxed_touching_crossing'
        elif rescues:
            rescue, change, votes = max(rescues, key=lambda v: (v[2], finite(v[0]['relative_depletion'], -np.inf)))
            row.update(min(contexts, key=lambda q: q['local_cross_fraction']))
            row.update(rescue_cap_bp=rescue['cap_bp'], rescue_profile_change=change, connection_evidence_votes=votes)
            route = 'combined_short_connection_evidence'
    row.update(graph_pass=int(not hard and bool(direct or rescues)), eligible=int(route != 'rejected'), accepted=0,
               acceptance_route=route, structural_reason=','.join(hard),
               raw_crossing_limit=opts.max_cross_fraction,
               overridden_veto='touching_crossing_fraction' if rescue else '',
               reject_reason=','.join(q for q in [row['shape_reason'], ','.join(hard),
                   'touching_crossing_no_combined_evidence' if not direct and not rescues else ''] if q))
    return row


def call_parent(parent, sig, es, opts, method, acrs, genes, cache):
    todo = [(parent['start'], parent['end'], 0)]; cuts = []; leaves = []; rows = []
    while todo:
        start, end, depth = todo.pop(); child = {**parent, 'start': start, 'end': end}
        proposals = F.proposals(sig, start, end, opts)
        if method == 'V7': proposals += cp_proposals(sig, start, end, opts)
        candidates = []
        for source, i in sorted(set(proposals), key=lambda z: (z[1], z[0])):
            r = evaluate(source, i, sig, child, es, opts, acrs, genes, cache)
            r.update(method=method, recursive_depth=depth); candidates.append(r)
        valid = [r for r in candidates if r['eligible']]
        if valid and depth < 16:
            best = max(valid, key=lambda r: (r['strength']/(1+r['local_cross_fraction']), r['scale_support_n'], -r['boundary_pos'], r['source']))
            best.update(accepted=1, reject_reason='accepted_'+best['acceptance_route'])
            cuts.append(dict(best)); p = best['boundary_pos']
            for r in valid:
                if r is not best: r['reject_reason'] = 'eligible_rechecked_in_child'
            todo.extend([(p, end, depth+1), (start, p, depth+1)])
        else:
            for r in valid: r['reject_reason'] = 'recursion_safety_limit'
            leaves.append(child)
        rows.extend(candidates)
    return leaves, cuts, rows


def audit_v5(rows, edges, out):
    reason = Counter(); category = defaultdict(list); detail = []; longs = []
    groups = {'crossing': {'touching_crossing_and_no_short_rescue'},
              'valley_depth': {'sharp_bilateral_fine_depth', 'step_relative_change'},
              'valley_width': {'sharp_width_or_single_bin', 'transition_not_fast', 'insufficient_raw_valley_duration'},
              'multiscale_support': {'insufficient_105k_sharp_support', 'step_multiscale_direction_or_magnitude'},
              'minimum_length_150kb': {'near_fixed_or_new_parent_edge'},
              'module_support': {'necessary_local_module_support', 'weak_daughter_internal', 'daughter_ACR_diversity'},
              'persistence': {'not_sustained_step', 'unstable_post_transition_platform', 'insufficient_sustained_flanks'},
              'global_cut_loss': {'global_cut_loss'}}
    arrays = {c: (np.array([e['left'] for e in edges[c]]), np.array([e['right'] for e in edges[c]])) for c in CHROMS}
    for r in rows:
        failed = set(r['reject_reason'].split(',')) if not int(r['accepted']) else set()
        reason.update(failed)
        for name, keys in groups.items():
            if keys & failed: category[name].append(r)
        pos, start, end = (int(r[k]) for k in ['boundary_pos', 'parent_start', 'parent_end'])
        left, right = arrays[r['chrom']]; crossing = (left < pos) & (right >= pos)
        vals = {}
        for flank in [250000, 500000]:
            a, b = max(start, pos-flank), min(end, pos+flank)
            touch = crossing & (((left >= a) & (left < b)) | ((right >= a) & (right < b)))
            n = F.unique_support(edges[r['chrom']], touch)[0]
            assert n == int(r[f'context_{flank}_crossing_touching_context'])
            for cap in [500000, 1000000]:
                long = F.unique_support(edges[r['chrom']], touch & (right-left > cap))[0]
                vals[f'context_{flank}_touching_long_gt{cap}_n'] = long
                vals[f'context_{flank}_touching_long_gt{cap}_fraction'] = long/n if n else np.nan
        inherited_short = min(pos-start, end-pos) < 500000
        detail.append({**r, **vals, 'would_fail_new_500kb_daughter_rule': int(inherited_short)})
        longs.append({'chrom': r['chrom'], 'boundary_pos': pos, 'parent_start': start, 'parent_end': end,
                      'source': r['source'], 'accepted_in_V5': r['accepted'], **vals})
    denom = sum(not int(r['accepted']) for r in rows)
    C.write(out/'02_V5_rejection_reason_summary.tsv', [
        {'reason': k, 'rejected_evaluation_n': v, 'all_rejected_evaluation_n': denom,
         'fraction_of_rejected_evaluations': v/denom, 'nonexclusive': 1} for k, v in reason.most_common()])
    C.write(out/'03_V5_rejection_category_summary.tsv', [
        {'condition': k, 'rejected_evaluation_n': len(v),
         'rejected_unique_position_n': len({(q['chrom'], q['boundary_pos']) for q in v}),
         'all_rejected_evaluation_n': denom, 'nonexclusive': 1} for k, v in category.items()])
    C.write(out/'04_V5_long_link_audit.tsv', longs)
    C.write(out/'05_V5_candidate_audit.tsv.gz', detail)
    print('V5 audit: evaluations', len(rows), 'unique positions', len({(r['chrom'],r['boundary_pos']) for r in rows}), 'rejections', reason.most_common(8), flush=True)


def select_random(cuts, n, seed):
    unique = {}
    for r in cuts:
        c, p = r['chrom'], int(r['boundary_pos'])
        if any(c == cc and a <= p < b for cc, a, b, focus in REGIONS): continue
        unique[(c, p)] = {**r, 'outside_representative_regions': 1}
    keys = sorted(unique)
    rng = np.random.default_rng(seed)
    chosen = sorted(rng.choice(len(keys), min(n, len(keys)), replace=False)) if keys else []
    return [unique[keys[i]] for i in chosen]


def main(args):
    root = args.outdir
    settings = {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items() if k not in ['stage', 'preview_approved']}
    root.mkdir(parents=True, exist_ok=True)
    setting_path = root/'01_test_settings.json'
    if setting_path.exists(): assert json.loads(setting_path.read_text()) == settings
    else: setting_path.write_text(json.dumps(settings, indent=2)+'\n')
    phase = root/('f01_preview' if args.stage == 'preview' else 'f02_chromosomes')
    if phase.exists(): raise ValueError('Output stage exists; never overwrite')
    if args.stage == 'chromosomes':
        assert args.preview_approved and (root/'f01_preview/f04_figures/f01_region_tracks/03_R_PLOTTING.DONE').exists()
    phase.mkdir()
    for d in ['f01_results', 'f02_plot_data', 'f03_evaluation', 'f04_figures', 'f05_scripts', 'f06_checks']: (phase/d).mkdir()
    base = args.baseline_dir; ref = args.reference_dir/'f02_chromosomes'
    meta = json.loads((base/'f01_baseline/01_parameters.json').read_text()); hist = Path(meta['history_root'])
    opts = config(args)
    basebed = base/'f01_baseline/02_CRU.bed'; v5bed = ref/'f01_results/f04_V5/02_CRU.bed'
    sources = [basebed, base/'f01_baseline/03_internal_boundaries.tsv', v5bed,
               ref/'f01_results/f04_V5/03_internal_boundaries.tsv', ref/'f01_results/f04_V5/04_candidate_evaluations.tsv',
               Path(meta['aglink_bedpe']), Path(meta['A_gap_file']), args.gene_bed, args.acr_bed,
               Path(__file__), Path(A.__file__), Path(F.__file__),
               args.reference_dir/'01_test_settings.json', base/'f01_baseline/01_parameters.json']
    sources += [hist/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz' for c in CHROMS]
    hashes = [{'path': str(p), 'md5_before': C.md5(p)} for p in sources]
    if args.stage == 'chromosomes':
        oldhash = {r['path']: r['md5_before'] for r in C.rows(root/'f01_preview/f06_checks/01_input_MD5.tsv')}
        assert all(oldhash[r['path']] == r['md5_before'] for r in hashes)
    baseline = P.read_bed(basebed, CHROMS); v5 = P.read_bed(v5bed, CHROMS)
    assert sum(map(len, baseline.values())) == 365 and sum(map(len, v5.values())) == 414
    assert C.md5(basebed) == 'a4eb35cc75719f543fa8b21d6a65fc30'
    acrs, genes = A.annotation_positions(args.gene_bed, args.acr_bed)
    edges = defaultdict(list)
    for r in C.rows(meta['aglink_bedpe']):
        if r['chrom'] not in CHROMS: continue
        for k in ['p1', 'p2', 'left', 'right']: r[k] = int(r[k])
        assert float(r['weight']) == 1
        edges[r['chrom']].append(r)
    oldc = C.rows(ref/'f01_results/f04_V5/04_candidate_evaluations.tsv')
    if args.stage == 'preview':
        audit_v5(oldc, edges, phase/'f06_checks')
    else:
        for name in ['02_V5_rejection_reason_summary.tsv', '03_V5_rejection_category_summary.tsv',
                     '04_V5_long_link_audit.tsv', '05_V5_candidate_audit.tsv.gz']:
            shutil.copy2(root/'f01_preview/f06_checks'/name, phase/'f06_checks'/name)
    results = {m: {'crus': [], 'boundaries': [], 'candidates': []} for m in METHODS}
    original = C.rows(base/'f01_baseline/03_internal_boundaries.tsv')
    oldbounds = C.rows(ref/'f01_results/f04_V5/03_internal_boundaries.tsv')
    for m, bed, bb in [('Baseline', baseline, original), ('V5_reference', v5, oldbounds)]:
        results[m]['crus'] = [r for c in CHROMS for r in bed[c]]
        results[m]['boundaries'] = [{**r, 'boundary_pos': int(r['boundary_pos']), 'method': m} for r in bb]
    results['V5_reference']['candidates'] = [{**r, 'method': 'V5_reference'} for r in oldc]
    sigs = {}; selected = []
    for c in CHROMS:
        with np.load(hist/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz') as z: sig = {k: z[k] for k in z.files}
        sigs[c] = sig
        for parent in v5[c]:
            if args.stage == 'preview' and not any(c == cc and parent['start'] < b and parent['end'] > a for cc, a, b, f in REGIONS): continue
            baseparent = next(r for r in baseline[c] if r['start'] <= parent['start'] and parent['end'] <= r['end'])
            par = {**parent, 'baseline_id': baseparent['cru_id'], 'v5_parent_id': parent['cru_id']}; selected.append(par)
            es = [e for e in edges[c] if e['left'] < parent['end'] and e['right'] >= parent['start']]; cache = {}
            for method in ['V6', 'V7']:
                rr, bb, cs = call_parent(par, sig, es, opts, method, acrs[c], genes[c], cache)
                results[method]['crus'].extend(rr); results[method]['boundaries'].extend(bb); results[method]['candidates'].extend(cs)
        for method in ['V6', 'V7']:
            results[method]['boundaries'].extend({**r, 'boundary_pos': int(r['boundary_pos']),
                'original_source': r.get('source', ''), 'source': 'inherited_V5', 'method': method} for r in oldbounds
                if r['chrom'] == c and any(p['chrom'] == c and p['start'] <= int(r['boundary_pos']) <= p['end'] for p in selected))
        print(args.stage, c, [(m, sum(r['chrom'] == c for r in v['crus'])) for m, v in results.items()], flush=True)
    for method, res in results.items():
        d = phase/'f01_results'/METHODS[method]; d.mkdir()
        if method in ['V6', 'V7']:
            for c in CHROMS:
                for i, r in enumerate(sorted((r for r in res['crus'] if r['chrom'] == c), key=lambda r: r['start']), 1): r['cru_id'] = f'{c}_{method}_CRU_{i}'
            for par in selected:
                rows = sorted((r for r in res['crus'] if r['v5_parent_id'] == par['v5_parent_id']), key=lambda r: r['start'])
                assert rows[0]['start'] == par['start'] and rows[-1]['end'] == par['end']
                assert all(l['end'] == r['start'] for l, r in zip(rows, rows[1:]))
            if args.stage == 'chromosomes':
                preview = P.read_bed(root/'f01_preview/f01_results'/METHODS[method]/'02_CRU.bed', CHROMS)
                coords = {(r['chrom'], r['start'], r['end']) for r in res['crus']}
                assert {(r['chrom'], r['start'], r['end']) for c in CHROMS for r in preview[c]} <= coords
        C.savebed(d/'02_CRU.bed', res['crus']); C.write(d/'03_internal_boundaries.tsv', res['boundaries'])
        C.write(d/'04_candidate_evaluations.tsv', res['candidates'], None if res['candidates'] else ['chrom', 'boundary_pos'])
        (d/'01_parameters.json').write_text(json.dumps({'method': method, 'stage': args.stage, 'rules': vars(opts)}, indent=2)+'\n')
    allnew = [{'method': m, **r} for m in ['V6', 'V7'] for r in results[m]['boundaries'] if r['source'] in NEW_SOURCES]
    sample = select_random(allnew, args.random_sample_n, args.random_seed) if args.stage == 'chromosomes' else []
    sample_rows = []
    for r in sample:
        sample_rows.append({**r, 'manual_inspection_status': 'pending_R_visual_check'})
    C.write(phase/'f03_evaluation/11_random_new_cut_sample.tsv', sample_rows,
            None if sample_rows else ['chrom', 'boundary_pos', 'method', 'manual_inspection_status'])
    regions = [{'region': f'{c}_{a}_{b}', 'chrom': c, 'start': a, 'end': b, 'focus': f, 'kind': 'representative'} for c, a, b, f in REGIONS]
    for j, r in enumerate(sample, 1):
        c, p = r['chrom'], int(r['boundary_pos']); a, b = max(0, p-1000000), min(meta['sizes'][c], p+1000000)
        regions.append({'region': f'random{j:02d}_{c}_{a}_{b}', 'chrom': c, 'start': a, 'end': b, 'focus': p, 'kind': 'random_new_cut'})
    oldparams = json.loads((hist/'f08_span_only_control/01_parameters.json').read_text())
    blockpath = Path(meta['block_bed']); m3path = Path(oldparams['m3_dir'])/'M3_CRU.bed'
    blocks = P.read_bed(blockpath, CHROMS); m3 = P.read_bed(m3path, CHROMS)
    hashes += [{'path': str(p), 'md5_before': C.md5(p)} for p in [blockpath, m3path]]
    pd = phase/'f02_plot_data'; arcs = []; profiles = []; intervals = []; plotblocks = []; cand = []; pgaps = []; bounds = []
    gaps = C.rows(meta['A_gap_file'])
    C.write(pd/'01_regions.tsv', regions)
    for reg in regions:
        c, a, b, region = (reg[k] for k in ['chrom', 'start', 'end', 'region']); sig = sigs[c]
        arcs.extend({'region': region, **{k: e[k] for k in ['chrom', 'left', 'right', 'acr_id', 'gene_id']}} for e in edges[c] if e['left'] < b and e['right'] > a)
        x = (np.arange(len(sig['raw']))+.5)*STEP
        for i in np.flatnonzero((x >= a-STEP) & (x <= b+STEP)):
            profiles.append({'region': region, 'chrom': c, 'position': x[i], 'raw_span': sig['raw'][i],
                'span_25k': sig['smooth'][0, i], 'span_105k': sig['smooth'][1, i], 'span_225k': sig['smooth'][2, i], 'local_background': sig['background'][i]})
        for method, res in results.items():
            intervals.extend({'region': region, 'method': method, **r} for r in res['crus'] if r['chrom'] == c and r['start'] < b and r['end'] > a)
            cand.extend({'region': region, **r} for r in res['candidates'] if r['chrom'] == c and a <= int(r['boundary_pos']) < b)
            bounds.extend({'region': region, 'method': method, **r} for r in res['boundaries'] if r['chrom'] == c and a <= int(r['boundary_pos']) < b)
        intervals.extend({'region': region, 'method': 'M3_reference', **r} for r in m3[c] if r['start'] < b and r['end'] > a)
        plotblocks.extend({'region': region, **r} for r in blocks[c] if r['start'] < b and r['end'] > a)
        pgaps.extend({'region': region, **r} for r in gaps if r['chrom'] == c and int(r['start']) < b and int(r['end']) > a)
    for name, rows in [('02_aglink_arcs.tsv.gz', arcs), ('03_signals.tsv', profiles), ('04_intervals.tsv', intervals), ('05_blocks.tsv', plotblocks), ('06_candidates.tsv', cand), ('07_frozen_gaps.tsv', pgaps), ('08_boundaries.tsv', bounds)]: C.write(pd/name, rows)
    for r in hashes:
        r['md5_after'] = C.md5(r['path']); assert r['md5_before'] == r['md5_after']
    C.write(phase/'f06_checks/01_input_MD5.tsv', hashes)
    C.write(phase/'f06_checks/06_selected_V5_parents.tsv', selected)
    cfg = {**settings, 'stage': args.stage, 'history_root': str(hist), 'sizes': meta['sizes'],
           'block_bed': meta['block_bed'], 'm3_script': oldparams['m3_script'], 'aglink_bedpe': meta['aglink_bedpe'],
           'frozen_V5_count': 414, 'no_coACR': True, 'new_rules': vars(opts)}
    (phase/'01_parameters.json').write_text(json.dumps(cfg, indent=2)+'\n')
    for p in Path(__file__).parent.iterdir():
        if p.is_file(): shutil.copy2(p, phase/'f05_scripts'/p.name)
    (phase/'03_CALLING.DONE').write_text('Frozen V5 and B145 geometry. New daughters >=500 kb; no coACR or block-guided cuts.\n')


if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    for k in ['baseline-dir', 'reference-dir', 'outdir', 'gene-bed', 'acr-bed']: ap.add_argument('--'+k, type=Path, required=True)
    ap.add_argument('--stage', choices=['preview', 'chromosomes'], required=True)
    ap.add_argument('--preview-approved', action='store_true')
    for k, value in [('sharp-depth', .30), ('medium-depth', .06), ('step-relative', .30), ('max-cross-fraction', .80), ('cp-candidate', .25), ('cp-accept', .30), ('cp-sse-gain', .25)]: ap.add_argument('--'+k, type=float, default=value)
    ap.add_argument('--sharp-max-width-bp', type=int, default=300000)
    ap.add_argument('--min-daughter-bp', type=int, default=500000)
    ap.add_argument('--random-sample-n', type=int, default=12)
    ap.add_argument('--random-seed', type=int, default=20261107)
    args = ap.parse_args()
    if args.min_daughter_bp < 500000: ap.error('This pilot does not relax daughters below 500 kb')
    main(args)
