# m08_baseline_targets: baseline targets and core components
# Frozen rule module, copied unchanged from the published calling chain
# (source: f12_four_target_tracks/01_prepare_targets.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Replay frozen generic rules for one missing chromosome, then export target tracks."""
import argparse
import csv
import gzip
import importlib.util
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
import numpy as np


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
    return mod


from . import m07_changepoint as N
A, F, C, P = N.A, N.F, N.C, N.P
from . import m04_primary_refine as E
D = E.D


def savebed(path, rows):
    with Path(path).open('w') as fh:
        for r in sorted(rows, key=lambda r: (r['chrom'], r['start'], r['end'])):
            fh.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['cru_id']}\n")


def coords(rows):
    return {(r['chrom'], int(r['start']), int(r['end'])) for r in rows}


def target_window(row, size, width, padding):
    a, b = int(row['start_1based'])-1, int(row['end_1based'])
    if not 0 <= a < b <= size: raise ValueError('Invalid 1-based inclusive target')
    width = max(width, b-a+2*padding)
    if width > size: raise ValueError('Window longer than chromosome')
    start = max(0, min(size-width, (a+b)//2-width//2)); end = start+width
    if a-start < padding or end-b < padding: raise ValueError('Cannot provide requested target flanks')
    return {**row, 'target_start': a, 'target_end': b, 'start': start, 'end': end,
            'focus': (a+b)//2, 'kind': 'user_target', 'region': f"{row['target_id']}_{row['chrom']}_{start}_{end}",
            'left_flank_bp': a-start, 'right_flank_bp': end-b, 'window_bp': width}


def baseline(chrom, size, sig, es, hard, edgeopts, rule):
    hs, mask, gaps, cores = C.identify_gaps(chrom, size, sig, es, hard, D)
    parents, actions = C.base_components(chrom, size, mask, sig, es, hard, D)
    refined, events = E.refine_edges(chrom, sig, es, gaps, cores, edgeopts)
    original = [{'chrom': chrom, 'start': p['start'], 'end': p['end'], 'cru_id': p['territory_id']} for p in parents]
    expanded, moves = E.apply_edges(original, refined)
    opts = SimpleNamespace(**rule); opts.valley_depth = .45
    result = {k: [] for k in ['crus', 'boundaries', 'candidates', 'actions']}
    for par in expanded:
        rr = C.recursive(chrom, par['start'], par['end'], es, sig, opts, D, 'primary45')
        for k in result: result[k].extend(rr[k])
    result['crus'].sort(key=lambda r: r['start'])
    for i, r in enumerate(result['crus'], 1): r['cru_id'] = f'{chrom}_Aedge_B145_CRU_{i}'
    return result, {'hard_signals': hs, 'mask': mask, 'original_gaps': gaps, 'cores': cores,
                    'refined_gaps': refined, 'edge_events': events, 'parent_actions': actions, 'refined_parents': expanded}


def annotation_positions(chrom, gene_bed, acr_bed):
    g, a = P.metadata(gene_bed, [chrom]), P.metadata(acr_bed, [chrom])
    genes = np.sort([end-1 if strand == '-' else start for c, start, end, strand in g.values()])
    acrs = np.sort([(start+end)//2 for c, start, end, strand in a.values()])
    return acrs, genes


def additions(chrom, parents, original, es, sig, opts, method, acrs, genes, v5=False, baseparents=None):
    result = {'crus': [], 'boundaries': [], 'candidates': []}
    for par in parents:
        parent = {**par, 'baseline_id': par['cru_id']}
        if not v5:
            source = next(r for r in baseparents if r['start'] <= par['start'] and par['end'] <= r['end'])
            parent.update(baseline_id=source['cru_id'], v5_parent_id=par['cru_id'])
        local = [e for e in es if e['left'] < par['end'] and e['right'] >= par['start']]
        fn = A.call_parent if v5 else N.call_parent
        rr, bb, cs = fn(parent, sig, local, opts, method, acrs, genes, {})
        result['crus'].extend(rr); result['boundaries'].extend(bb); result['candidates'].extend(cs)
    result['boundaries'].extend({**r, 'boundary_pos': int(r['boundary_pos']), 'original_source': r.get('source', ''),
        'source': 'inherited_B145' if v5 else 'inherited_V5', 'method': method} for r in original)
    result['crus'].sort(key=lambda r: r['start'])
    for i, r in enumerate(result['crus'], 1): r['cru_id'] = f'{chrom}_{method}_CRU_{i}'
    oldpts = {(r['chrom'], p) for r in parents for p in [r['start'], r['end']]}
    pts = {(r['chrom'], p) for r in result['crus'] for p in [r['start'], r['end']]}
    assert oldpts <= pts and sum(r['end']-r['start'] for r in parents) == sum(r['end']-r['start'] for r in result['crus'])
    for r in result['boundaries']:
        if r['source'] in N.NEW_SOURCES:
            assert min(int(r['boundary_pos'])-int(r['parent_start']), int(r['parent_end'])-int(r['boundary_pos'])) >= opts.min_daughter_bp
    return result


def main(args):
    root, hist, newchrom = args.outdir, args.history_root, args.add_chrom
    if root.exists(): raise ValueError('Output exists; never overwrite')
    root.mkdir(parents=True)
    for f in ['f01_chromosome_extension', 'f02_plot_data', 'f03_checks', 'f04_figures', 'f05_scripts']: (root/f).mkdir()
    inputs = root/'f01_chromosome_extension/f01_inputs'; inputs.mkdir()
    prior = hist/'f20_V6_V7_local_changepoint/f02_chromosomes'
    cfg = json.loads((prior/'01_parameters.json').read_text())
    params = json.loads((hist/'f08_span_only_control/01_parameters.json').read_text())
    basemeta = json.loads((hist/'f15_Aedge_B145_coACR_support/f01_baseline/01_parameters.json').read_text())
    hard = SimpleNamespace(**json.loads((hist/'f13_B_route_hard_gap_first/01_parameters.json').read_text()))
    edgeopts = SimpleNamespace(**json.loads((hist/'f14_B_route_edge_and_long_split/01_parameters.json').read_text()))
    rule = basemeta['primary_rule']
    v5settings = json.loads((hist/'f19_B145_V4_V5_short_crossing/01_test_settings.json').read_text())
    v5opts = A.tier(SimpleNamespace(**v5settings), 'V5'); currentopts = SimpleNamespace(**cfg['new_rules'])
    selected = C.rows(args.targets); regions = [target_window(r, cfg['sizes'][r['chrom']], args.window_bp, args.padding_bp) for r in selected]
    assert newchrom not in ['1A', '1B', '2A'] and args.replay_chrom in ['1A', '1B', '2A']
    paths = [args.targets, Path(__file__), Path(N.__file__), Path(A.__file__), Path(F.__file__), Path(E.__file__), Path(D.__file__), Path(P.__file__), Path(C.__file__),
             Path(params['aglink']), Path(params['gene_bed']), Path(params['acr_bed']), Path(params['bedpe_reference']),
             prior/'01_parameters.json', hist/'f08_span_only_control/01_parameters.json',
             hist/'f13_B_route_hard_gap_first/01_parameters.json', hist/'f14_B_route_edge_and_long_split/01_parameters.json',
             hist/'f15_Aedge_B145_coACR_support/f01_baseline/01_parameters.json', hist/'f19_B145_V4_V5_short_crossing/01_test_settings.json',
             hist/'f15_Aedge_B145_coACR_support/f01_baseline/02_CRU.bed', Path(basemeta['A_gap_file']),
             hist/'f08_span_only_control/f01_shared_input/01_count_only_aglink.annotated.bedpe']
    for method, folder in N.METHODS.items():
        paths.extend(prior/'f01_results'/folder/name for name in ['02_CRU.bed', '03_internal_boundaries.tsv', '04_candidate_evaluations.tsv'])
    paths += [hist/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz' for c in set(r['chrom'] for r in regions) & set(['1A', '1B', '2A'])]
    before = [{'path': str(p), 'md5_before': C.md5(p)} for p in dict.fromkeys(paths)]
    cache = hist/'f08_span_only_control/f01_shared_input/01_count_only_aglink.annotated.bedpe'
    existing_edges = D.load_edges(cache, ['1A', '1B', '2A'])
    # Replay the generic upstream chain on an already verified full chromosome.
    replay = args.replay_chrom
    with np.load(hist/f'f13_B_route_hard_gap_first/f01_inputs/{replay}.B_signals.npz') as z: replay_sig = {k: z[k] for k in z.files}
    replay_result, replay_aux = baseline(replay, cfg['sizes'][replay], replay_sig, existing_edges[replay], hard, edgeopts, rule)
    oldbase = P.read_bed(hist/'f15_Aedge_B145_coACR_support/f01_baseline/02_CRU.bed', [replay])[replay]
    assert coords(replay_result['crus']) == coords(oldbase), 'Generic A-edge+B1-45 replay must exactly reproduce the frozen chromosome'
    C.write(root/'f03_checks/01_existing_chromosome_replay.tsv', [{'chrom': replay, 'saved_cru_n': len(oldbase), 'replayed_cru_n': len(replay_result['crus']), 'coordinate_set_equal': 'PASS'}])
    loadargs = SimpleNamespace(**{k: params[k] for k in ['gene_bed', 'acr_bed', 'aglink', 'bedpe_reference']}, chroms=[newchrom])
    loaded = P.load_count_only_edges(loadargs, cfg['sizes'], inputs)
    es = D.load_edges(inputs/'01_count_only_aglink.annotated.bedpe', [newchrom])[newchrom]
    sigargs = SimpleNamespace(**params); sigargs.bin_bp = 5000
    original_signals = P.signals(loaded[newchrom], cfg['sizes'][newchrom], sigargs)
    sig = D.territory_signals(original_signals['raw'], np.array(original_signals['smooth']), original_signals['endpoint_broad'], SimpleNamespace(**rule))
    sig['old_active'] = original_signals['active']
    result, aux = baseline(newchrom, cfg['sizes'][newchrom], sig, es, hard, edgeopts, rule)
    np.savez_compressed(inputs/f'{newchrom}.signals.npz', **sig, **aux['hard_signals'], hard_excluded=aux['mask'])
    for name, rows in [('03_original_hard_gaps.tsv', aux['original_gaps']), ('04_hard_gap_cores.tsv', aux['cores']),
                       ('05_A_refined_gaps.tsv', aux['refined_gaps']), ('06_A_edge_events.tsv', aux['edge_events']),
                       ('07_parent_actions.tsv', aux['parent_actions']), ('08_A_refined_parents.tsv', aux['refined_parents'])]: C.write(inputs/name, rows)
    called = {'Baseline': result}; acrs, genes = annotation_positions(newchrom, params['gene_bed'], params['acr_bed'])
    called['V5_reference'] = additions(newchrom, result['crus'], result['boundaries'], es, sig, v5opts, 'V5', acrs, genes, v5=True)
    for method in ['V6', 'V7']:
        called[method] = additions(newchrom, called['V5_reference']['crus'], called['V5_reference']['boundaries'], es, sig, currentopts, method, acrs, genes, baseparents=result['crus'])
    statistics = []
    for method, folder in [('Baseline', 'f02_B145'), ('V5_reference', 'f03_V5'), ('V6', 'f04_V6'), ('V7', 'f05_V7')]:
        d = root/'f01_chromosome_extension'/folder; d.mkdir(); r = called[method]
        savebed(d/'02_CRU.bed', r['crus']); C.write(d/'03_internal_boundaries.tsv', r['boundaries'])
        C.write(d/'04_candidate_evaluations.tsv', r['candidates'], None if r['candidates'] else ['chrom', 'boundary_pos'])
        (d/'01_parameters.json').write_text(json.dumps({'method': method, 'chrom': newchrom, 'whole_chromosome': True,
            'primary_rule': rule, 'hard_gap_rule': vars(hard), 'A_edge_rule': vars(edgeopts),
            'additional_rule': vars(v5opts) if method == 'V5_reference' else vars(currentopts) if method in ['V6','V7'] else None}, indent=2)+'\n')
        stat = P.interval_summary(method, newchrom, r['crus'], es, cfg['sizes'][newchrom]); stat['internal_cut_n'] = len(r['boundaries']); statistics.append(stat)
        print('Full chromosome', newchrom, method, 'CRUs', stat['cru_n'], 'coverage', stat['coverage_fraction'], flush=True)
    C.write(root/'f03_checks/02_added_chromosome_summary.tsv', statistics)
    available = ['1A', '1B', '2A', newchrom]; bed = {}; bounds = {}; cand = {}
    for method, folder in N.METHODS.items():
        d = prior/'f01_results'/folder
        bed[method] = P.read_bed(d/'02_CRU.bed', ['1A','1B','2A'])
        bed[method][newchrom] = called[method]['crus']
        bounds[method] = C.rows(d/'03_internal_boundaries.tsv') + called[method]['boundaries']
        cand[method] = C.rows(d/'04_candidate_evaluations.tsv') + called[method]['candidates']
    signals = {newchrom: sig}; edges = {**existing_edges, newchrom: es}
    for c in set(r['chrom'] for r in regions)-{newchrom}:
        with np.load(hist/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz') as z: signals[c] = {k:z[k] for k in z.files}
    blocks = P.read_bed(params['block_bed'], available); m3 = P.read_bed(Path(params['m3_dir'])/'M3_CRU.bed', available)
    before += [{'path': str(p), 'md5_before': C.md5(p)} for p in [Path(params['block_bed']), Path(params['m3_dir'])/'M3_CRU.bed', cache]]
    gaps = C.rows(basemeta['A_gap_file'])+aux['refined_gaps']
    arcs=[]; profiles=[]; intervals=[]; plotblocks=[]; candidates=[]; pgaps=[]; boundaries=[]; visible=[]
    for reg in regions:
        c, start, end, rid = (reg[k] for k in ['chrom','start','end','region']); s=signals[c]
        arcs.extend({'region':rid, **{k:e[k] for k in ['chrom','left','right','acr_id','gene_id']}} for e in edges[c] if e['left']<end and e['right']>start)
        x=(np.arange(len(s['raw']))+.5)*5000
        for i in np.flatnonzero((x>=start-5000)&(x<=end+5000)):
            profiles.append({'region':rid,'chrom':c,'position':x[i],'raw_span':s['raw'][i],
                'span_25k':s['smooth'][0,i],'span_105k':s['smooth'][1,i],'span_225k':s['smooth'][2,i],'local_background':s['background'][i]})
        for method in N.METHODS:
            rows=[r for r in bed[method][c] if r['start']<end and r['end']>start]
            intervals.extend({'region':rid,'method':method,**r} for r in rows)
            visible.extend({'region':rid,'method':method,'target_overlap_bp':max(0,min(r['end'],reg['target_end'])-max(r['start'],reg['target_start'])),**r} for r in rows)
            candidates.extend({'region':rid,**r,'method':method} for r in cand[method] if r['chrom']==c and start<=int(r['boundary_pos'])<end)
            boundaries.extend({'region':rid,**r,'method':method} for r in bounds[method] if r['chrom']==c and start<=int(r['boundary_pos'])<end)
        intervals.extend({'region':rid,'method':'M3_reference',**r} for r in m3[c] if r['start']<end and r['end']>start)
        visible.extend({'region':rid,'method':'M3_reference','target_overlap_bp':max(0,min(r['end'],reg['target_end'])-max(r['start'],reg['target_start'])),**r} for r in m3[c] if r['start']<end and r['end']>start)
        plotblocks.extend({'region':rid,**r} for r in blocks[c] if r['start']<end and r['end']>start)
        pgaps.extend({'region':rid,**r} for r in gaps if r['chrom']==c and int(r['start'])<end and int(r['end'])>start)
    pd=root/'f02_plot_data'; C.write(pd/'01_regions.tsv',regions)
    P.write_tsv(pd/'02_aglink_arcs.tsv.gz',arcs)
    for name, rows in [('03_signals.tsv',profiles),('04_intervals.tsv',intervals),('05_blocks.tsv',plotblocks),
                       ('06_candidates.tsv',candidates),('07_frozen_gaps.tsv',pgaps),('08_boundaries.tsv',boundaries)]: C.write(pd/name,rows)
    C.write(root/'f03_checks/03_visible_CRUs_and_target_overlap.tsv',visible)
    for r in before:
        r['md5_after']=C.md5(r['path']); assert r['md5_after']==r['md5_before']
    C.write(root/'f03_checks/04_input_MD5.tsv',before)
    settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}
    (root/'01_parameters.json').write_text(json.dumps({**settings,'target_coordinates':'1-based inclusive input -> BED [start-1,end)',
        'whole_chromosome_added':newchrom,'existing_1A_2A_reused':True,'calling_rules_unchanged':True,
        'block_enrichment_recomputed':False,'no_coACR':True},indent=2)+'\n')
    for p in Path(__file__).parent.iterdir():
        if p.is_file(): shutil.copy2(p,root/'f05_scripts'/p.name)
    (root/'02_DATA.DONE').write_text('Exact baseline replay passed; one missing chromosome called with frozen rules; target tracks exported without window-local calling.\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ['history-root','outdir','targets']: ap.add_argument('--'+k,type=Path,required=True)
    ap.add_argument('--add-chrom',default='7B'); ap.add_argument('--replay-chrom',default='1A')
    ap.add_argument('--window-bp',type=int,default=2000000); ap.add_argument('--padding-bp',type=int,default=500000)
    main(ap.parse_args())
