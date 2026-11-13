# m02_adaptive_background: adaptive background density
# Frozen rule module, copied unchanged from the published calling chain
# (source: f03_adaptive_span_pilot/01_adaptive_span_pilot.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Staged agLINK-only pilot: adaptive territories and parent-local deep valleys."""
import argparse
import bisect
import csv
import gzip
import hashlib
import importlib.util
import inspect
import json
import math
import runpy
import shutil
from collections import defaultdict
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.ndimage import maximum_filter1d, percentile_filter, uniform_filter1d
from scipy.signal import find_peaks


FOLDERS = {'Adaptive': 'f01_new_adaptive', 'Old_V1': 'f02_old_V1', 'Old_V2': 'f03_old_V2',
           'Old_V3': 'f04_old_V3', 'M3': 'f05_M3_reference'}
COLORS = {'Adaptive': '#5FB6BD', 'Old_V1': '#96CEDD', 'Old_V2': '#E89B6C', 'Old_V3': '#C1C0BF', 'M3': '#B086CF'}


def module(path):
    spec = importlib.util.spec_from_file_location('previous_span_utilities', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def read(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def write(path, rows, fields=None):
    rows = list(rows)
    fields = fields or list(dict.fromkeys(k for r in rows for k in r)) or ['chrom']
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


def sha(path):
    h = hashlib.sha256()
    with open(path, 'rb') as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


def runs(mask):
    changes = np.diff(np.r_[False, mask, False].astype(int))
    return list(zip(np.flatnonzero(changes == 1), np.flatnonzero(changes == -1)))


def persistent(mask, min_bins):
    out = np.zeros(len(mask), dtype=bool)
    for a, b in runs(mask):
        if b - a >= min_bins:
            out[a:b] = True
    return out


def local_background(fine, args):
    k = max(1, args.background_tile_bp // args.bin_bp)
    floors = np.array([np.quantile(fine[i:i+k], .1) for i in range(0, len(fine), k)])
    width = max(3, args.background_window_bp // args.background_tile_bp) | 1
    bg = percentile_filter(floors, percentile=20, size=width, mode='nearest')
    bg = np.repeat(bg, k)[:len(fine)]
    return uniform_filter1d(bg, max(1, args.smooth_bp[-1] // args.bin_bp) | 1, mode='nearest')


def side_maxima(values, flank):
    flank = max(1, flank)
    left_shift = np.r_[values[0], values[:-1]]
    right_shift = np.r_[values[1:], values[-1]]
    left = maximum_filter1d(left_shift, size=flank, origin=(flank-1)//2, mode='nearest')
    right = maximum_filter1d(right_shift, size=flank, origin=-(flank//2), mode='nearest')
    return left, right


def territory_signals(raw, smooth, endpoint, args):
    raw = np.maximum(raw, 0)
    smooth = np.maximum(smooth, 0)
    background = local_background(smooth[0], args)
    threshold = np.maximum(args.min_seed_depth, background + args.background_sigma * np.sqrt(background + 1))
    seed = (np.count_nonzero(smooth >= threshold, axis=0) >= 2) & (raw >= args.min_seed_depth) & (endpoint >= args.seed_endpoints)
    left, right = side_maxima(smooth[1], args.gap_flank_bp // args.bin_bp)
    weak = np.minimum(left, right)
    relative = np.count_nonzero(smooth <= weak * args.gap_fraction, axis=0) >= 2
    below_bg = np.count_nonzero(smooth <= threshold, axis=0) >= 2
    low_relative = relative & below_bg
    low_thin = (raw < args.bridge_min_links) & (smooth[0] < args.bridge_min_links) & (smooth[1] < args.bridge_min_links) & (endpoint < args.seed_endpoints)
    low_zero = raw <= 1e-9
    low = low_zero | low_relative | low_thin
    block = persistent(low, max(1, math.ceil(args.gap_min_bp / args.bin_bp)))
    return {'raw': raw, 'smooth': smooth, 'endpoint': endpoint, 'background': background,
            'seed_threshold': threshold, 'seed': seed, 'low_zero': low_zero,
            'low_relative': low_relative, 'low_thin': low_thin, 'blocked': block, 'domain': ~block,
            'left_reference': left, 'right_reference': right}


def internal_n(left, right, start, end):
    return int(np.count_nonzero((left >= start) & (right < end)))


def territories(chrom, size, edges, sig, args):
    left, right = np.array([e['left'] for e in edges]), np.array([e['right'] for e in edges])
    kept, actions = [], []
    for a, b in runs(sig['domain']):
        positives = np.flatnonzero(sig['raw'][a:b] > 1e-9)
        if not len(positives):
            continue
        lo, hi = a + int(positives[0]), a + int(positives[-1]) + 1
        start, end = lo * args.bin_bp, min(size, hi * args.bin_bp)
        seed_n = int(np.count_nonzero(sig['seed'][lo:hi]))
        support = internal_n(left, right, start, end)
        reasons = []
        if not seed_n: reasons.append('no_locally_enriched_seed')
        if end - start < args.min_cru_bp: reasons.append('too_short')
        if support < args.min_internal: reasons.append('insufficient_internal_links')
        row = {'chrom': chrom, 'start': start, 'end': end, 'domain_start': a*args.bin_bp,
               'domain_end': min(size,b*args.bin_bp), 'seed_bin_n': seed_n, 'aglink_internal_n': support,
               'passed': int(not reasons), 'reason': ','.join(reasons) or 'kept',
               'zero_only_edge_trim_bp': (lo-a+b-hi)*args.bin_bp}
        actions.append(row)
        if not reasons:
            kept.append({**row, 'territory_id': f'{chrom}_T{len(kept)+1}'})
    gaps = []
    for a,b in runs(sig['blocked']):
        causes = [name for name in ['low_zero','low_relative','low_thin'] if np.any(sig[name][a:b])]
        gaps.append({'chrom': chrom,'start': a*args.bin_bp,'end':min(size,b*args.bin_bp),
                     'cause':','.join(causes),'median_raw_depth':float(np.median(sig['raw'][a:b])),
                     'median_background':float(np.median(sig['background'][a:b]))})
    return kept, actions, gaps


def score_scale(y, idx, a, b, peaks, args):
    flank = args.shoulder_bp // args.bin_bp
    lp = peaks[(peaks < idx) & (peaks >= max(a, idx-flank))]
    rp = peaks[(peaks > idx) & (peaks < min(b, idx+flank+1))]
    if len(lp): li = int(lp[-1])
    else:
        start = max(a,idx-flank)
        li = start + int(np.argmax(y[start:idx]))
    if len(rp): ri = int(rp[0])
    else: ri = idx+1+int(np.argmax(y[idx+1:min(b,idx+flank+1)]))
    floor, lpeak, rpeak = float(y[idx]),float(y[li]),float(y[ri])
    dl,dr = (lpeak-floor)/max(lpeak,1e-12),(rpeak-floor)/max(rpeak,1e-12)
    level = floor + .5*max(0,min(lpeak,rpeak)-floor)
    lo,hi=idx,idx
    while lo>li and y[lo-1] <= level: lo-=1
    while hi<ri and y[hi+1] <= level: hi+=1
    return {'floor':floor,'left_peak':lpeak,'right_peak':rpeak,'left_depth':dl,'right_depth':dr,
            'depth':min(dl,dr),'width_bp':(hi-lo+1)*args.bin_bp,
            'left_peak_pos':li*args.bin_bp+args.bin_bp//2,'right_peak_pos':ri*args.bin_bp+args.bin_bp//2}


def parent_candidates(chrom, start, end, left, right, sig, args, depth=0):
    # Recursion cuts are bin centers, not bin edges. Search only centers in
    # [start,end), including a center at start but never one at end.
    a=max(0,math.ceil(start/args.bin_bp-.5))
    b=min(len(sig['raw']),math.ceil(end/args.bin_bp-.5))
    if b-a<3:return []
    indices, props = find_peaks(-sig['smooth'][0,a:b],prominence=args.min_prominence,
                                distance=max(1,args.min_cru_bp//args.bin_bp))
    peaks=[]
    for y in sig['smooth']:
        positions,properties=find_peaks(y[a:b],prominence=args.min_peak_prominence)
        absolute=positions+a
        supported=properties['prominences']>=args.min_relative_peak_prominence*np.maximum(y[absolute],1e-12)
        peaks.append(absolute[supported])
    result=[]
    for j,p in enumerate(indices):
        idx=a+int(p); pos=idx*args.bin_bp+args.bin_bp//2
        scores=[score_scale(y,idx,a,b,pk,args) for y,pk in zip(sig['smooth'],peaks)]
        scale_pass=[s['depth']>=args.valley_depth and s['width_bp']>=args.valley_min_width_bp for s in scores]
        cross_all=(left<pos)&(right>pos)
        global_n=int(np.count_nonzero(cross_all))
        contexts=[]
        for flank in args.support_flanks_bp:
            local_start,local_end=max(start,pos-flank),min(end,pos+flank)
            il,ir=internal_n(left,right,local_start,pos),internal_n(left,right,pos,local_end)
            local=int(np.count_nonzero(cross_all&(left>=local_start)&(right<local_end)))
            weak=min(il,ir);frac=local/(local+weak) if local+weak else 1
            contexts.append((flank,il,ir,local,frac))
        good=[v for v in contexts if min(v[1],v[2])>=args.min_internal and v[4]<=args.max_local_cross_fraction]
        flank,il,ir,local,frac=good[0] if good else contexts[-1]
        weak=min(il,ir)
        reasons=[]
        if pos-start<args.min_cru_bp or end-pos<args.min_cru_bp:reasons.append('daughter_length')
        if not scale_pass[0]:reasons.append('fine_depth_or_width')
        if sum(scale_pass)<args.min_scale_support:reasons.append('multiscale_depth_or_width')
        if weak<args.min_internal:reasons.append('local_internal_support')
        if frac>args.max_local_cross_fraction:reasons.append('local_crossing_fraction')
        if min(internal_n(left,right,start,pos),internal_n(left,right,pos,end))<args.min_internal:reasons.append('daughter_internal_support')
        row={'chrom':chrom,'parent_start':start,'parent_end':end,'recursive_depth':depth,'boundary_pos':pos,
             'prominence':float(props['prominences'][j]),'scale_support_n':sum(scale_pass),
             'left_internal':il,'right_internal':ir,'crossing_local':local,'crossing_global':global_n,
             'outside_local_crossing':global_n-local,'local_cross_fraction':frac,
             'support_context_bp':flank,
             'local_background':float(sig['background'][idx]),'eligible':int(not reasons),'accepted':0,
             'reject_reason':','.join(reasons),'score':sum(s['depth'] for s in scores if s['depth']>=args.valley_depth)}
        for label,s in zip(['25k','105k','225k'],scores):
            row.update({label+'_'+key:value for key,value in s.items()})
        for f,lval,rval,cval,ratio in contexts:
            row.update({f'context_{f}_left_internal':lval,f'context_{f}_right_internal':rval,
                        f'context_{f}_crossing':cval,f'context_{f}_cross_fraction':ratio})
        result.append(row)
    return result


def split_territory(territory, edges, sig, args):
    left,right=np.array([e['left'] for e in edges]),np.array([e['right'] for e in edges])
    todo=[(territory['start'],territory['end'],0)]
    candidates,boundaries,leaves,actions=[],[],[],[]
    while todo:
        start,end,depth=todo.pop()
        rows=parent_candidates(territory['chrom'],start,end,left,right,sig,args,depth)
        candidates.extend(rows)
        eligible=[r for r in rows if r['eligible']]
        if eligible and depth<args.max_recursion_depth:
            best=max(eligible,key=lambda r:(r['scale_support_n'],r['score'],r['prominence'],-r['boundary_pos']))
            best['accepted']=1
            for r in eligible:
                if r is not best:r['reject_reason']='eligible_not_selected_in_this_parent'
            boundaries.append(dict(best))
            pos=best['boundary_pos']
            todo.extend([(pos,end,depth+1),(start,pos,depth+1)])
            action='accepted_deep_multiscale_local_split'
        else:
            leaves.append({'chrom':territory['chrom'],'start':start,'end':end,'territory_id':territory['territory_id']})
            action='kept_no_qualified_valley' if depth<args.max_recursion_depth else 'recursion_safety_limit'
        actions.append({'chrom':territory['chrom'],'start':start,'end':end,'depth':depth,'candidate_n':len(rows),
                        'eligible_n':len(eligible),'action':action})
    return {'crus':leaves,'boundaries':boundaries,'candidates':candidates,'actions':actions,'territories':[territory]}


def load_edges(path, chroms):
    out=defaultdict(list)
    for row in read(path):
        if row['chrom'] not in chroms:continue
        for k in ['p1','p2','left','right','span_start','span_end','distance_bp','source_line']:
            row[k]=int(row[k])
        row['weight']=float(row['weight'])
        if row['weight']!=1:raise ValueError('Input is not count-only')
        out[row['chrom']].append(row)
    return out


def freeze(args, target):
    target.mkdir(parents=True,exist_ok=True)
    for path in Path(__file__).parent.iterdir():
        if path.is_file():shutil.copy2(path,target/path.name)
    shutil.copy2(args.pilot_script,target/'07_previous_pilot_utilities.py')


def config(args):
    keys=['bin_bp','smooth_bp','background_tile_bp','background_window_bp','background_sigma','min_seed_depth',
          'seed_endpoints','gap_min_bp','gap_fraction','gap_flank_bp','bridge_min_links','min_cru_bp','min_internal',
          'shoulder_bp','valley_depth','valley_min_width_bp','min_prominence','min_peak_prominence','min_relative_peak_prominence','min_scale_support',
          'support_flanks_bp','max_local_cross_fraction','max_recursion_depth']
    rules={k:getattr(args,k) for k in keys}
    algorithm='\n'.join(inspect.getsource(f) for f in [runs,persistent,local_background,side_maxima,territory_signals,internal_n,territories,score_scale,parent_candidates,split_territory])
    rules['algorithm_sha256']=hashlib.sha256(algorithm.encode()).hexdigest()
    rules['weighting']='all agLINKs weight=1; PCC/coACR not read; no distance decay'
    return rules


def prepare(args,pilot):
    shared=args.outdir/'f01_shared'
    if (shared/'01_rule_signature.json').exists():
        if json.loads((shared/'01_rule_signature.json').read_text())!=config(args):raise ValueError('Existing rule signature differs')
        return
    shared.mkdir(parents=True,exist_ok=True)
    source=json.loads((args.source_dir/'01_parameters.json').read_text())
    if not set(args.chroms)<=set(source['chroms']):raise ValueError('Cached chromosome missing')
    provenance=[]
    for row in read(args.source_dir/'f01_shared_input/03_input_provenance.tsv'):
        if sha(row['path'])!=row['sha256']:raise ValueError('Original source changed: '+row['path'])
        provenance.append(row)
    profile_fields={'raw':'span_coverage','endpoint':'endpoint_broad','fine':'smooth_fine','medium':'smooth_medium','broad':'smooth_broad','old_active':'active'}
    for chrom in args.chroms:
        path=args.source_dir/f'f01_shared_input/{chrom}.profiles.tsv.gz'
        rows=read(path)
        data={k:np.array([float(r[v]) for r in rows]) for k,v in profile_fields.items()}
        smooth=np.array([data.pop('fine'),data.pop('medium'),data.pop('broad')])
        sig=territory_signals(data['raw'],smooth,data['endpoint'],args)
        sig['old_active']=data['old_active']
        np.savez_compressed(shared/f'{chrom}.signals.npz',**sig)
        provenance.append({'role':'old_profile','path':str(path),'sha256':sha(path)})
    bedpe=args.source_dir/'f01_shared_input/01_count_only_aglink.annotated.bedpe'
    provenance.append({'role':'ID_preserving_aglink_cache','path':str(bedpe),'sha256':sha(bedpe)})
    write(shared/'02_provenance.tsv',provenance)
    (shared/'03_original_parameters.json').write_text(json.dumps(source,indent=2)+'\n')
    (shared/'01_rule_signature.json').write_text(json.dumps(config(args),indent=2)+'\n')


def references(args,pilot,regions_only=False):
    source=json.loads((args.source_dir/'01_parameters.json').read_text())
    out={}
    for method,folder in [('Old_V1','f02_V1_span_valley'),('Old_V2','f03_V2_recursive_background'),('Old_V3','f04_V3_crossing_filter')]:
        by_chr=pilot.read_bed(args.source_dir/folder/'02_CRU.bed',args.chroms)
        bounds=read(args.source_dir/folder/'03_internal_boundaries.tsv')
        bounds=[{**r,'boundary_pos':int(r['boundary_pos'])} for r in bounds if r['chrom'] in args.chroms]
        out[method]={'crus':[r for c in args.chroms for r in by_chr[c]],'boundaries':bounds,'candidates':[],'actions':[],'territories':[]}
    by_chr=pilot.read_bed(Path(source['m3_dir'])/'M3_CRU.bed',args.chroms)
    bounds=[]
    for name in ['M3_primary_boundaries.bed','M3_secondary_boundaries.bed']:
        for rows in pilot.read_bed(Path(source['m3_dir'])/name,args.chroms).values():
            bounds.extend({**r,'boundary_pos':(r['start']+r['end'])//2} for r in rows)
    out['M3']={'crus':[r for c in args.chroms for r in by_chr[c]],'boundaries':bounds,'candidates':[],'actions':[],'territories':[]}
    return out,source


def old_diagnosis(args,sigs,refs,terr,new,phase):
    output=[]
    oldroot=args.source_dir.parent
    guard={}
    for c in args.chroms:
        cached=read(oldroot/f'f01_shared_input/{c}.profiles.tsv.gz')
        guard[c]={'active':np.array([int(r['active']) for r in cached],dtype=bool),
                  'threshold':np.array([float(r['activity_threshold']) for r in cached])}
        del cached
    for c,start,end in args.regions:
        label=f'{c}:{start}-{end}'
        for setting,base in [('span_only',args.source_dir),('endpoint_guard_2',oldroot)]:
            params=json.loads((base/'01_parameters.json').read_text())
            sig=sigs[c]
            lo,hi=start//args.bin_bp,min(len(sig['raw']),math.ceil(end/args.bin_bp))
            oldactive=sig['old_active'].astype(bool) if setting=='span_only' else guard[c]['active']
            for a,b in runs(~oldactive[lo:hi]):
                a,b=a+lo,b+lo
                causes=[]
                if np.any(sig['raw'][a:b]<params['min_span_count']):causes.append('old_raw_span_floor')
                threshold=guard[c]['threshold'][a:b]
                if np.any(sig['smooth'][2,a:b]<threshold):causes.append('old_chromosome_broad_threshold')
                if np.any(sig['endpoint'][a:b]<params['min_endpoint_count']-1e-9):causes.append('endpoint_guard')
                output.append({'region':label,'setting':setting,'event':'excluded_activity_run',
                               'start':a*args.bin_bp,'end':b*args.bin_bp,'cause':','.join(causes),
                               'median_raw':float(np.median(sig['raw'][a:b])), 'median_endpoint':float(np.median(sig['endpoint'][a:b])),
                               'territory_threshold':float(np.median(threshold)), 'trim':'none in old span-only caller'})
            for method,folder in [('V1','f02_V1_span_valley'),('V2','f03_V2_recursive_background'),('V3','f04_V3_crossing_filter')]:
                for row in read(base/folder/'03_internal_boundaries.tsv'):
                    if row['chrom']==c and start<=int(row['boundary_pos'])<=end:
                        output.append({'region':label,'setting':setting,'event':method+'_accepted_cut',
                                       'start':int(row['boundary_pos']),'end':int(row['boundary_pos']),
                                       'cause':'inherited_V1_valley' if row.get('stage')=='V1' else 'recursive_background_split',
                                       'old_relative_depth':row.get('relative_depth',''),'old_background_depth':row.get('background_depth',''),
                                       'trim':'none'})
        for gap in terr[c]['gaps']:
            if gap['start']<end and gap['end']>start:output.append({'region':label,'setting':'new_adaptive','event':'sustained_gap',**gap})
        for r in new['candidates']:
            if r['chrom']==c and start<=r['boundary_pos']<=end:
                output.append({'region':label,'setting':'new_adaptive','event':'candidate_evaluation',
                               'start':r['boundary_pos'],'end':r['boundary_pos'],'cause':r['reject_reason'] or 'accepted',
                               'parent_start':r['parent_start'],'parent_end':r['parent_end'],
                               'fine_depth':r['25k_depth'],'medium_depth':r['105k_depth'],'broad_depth':r['225k_depth'],
                               'accepted':r['accepted']})
    write(phase/'f06_comparison/07_region_event_causes.tsv',output)


def render(region,sig,edges,results,terr,blocks,phase,args,helper):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle
    from matplotlib.lines import Line2D
    from matplotlib.ticker import FuncFormatter,MaxNLocator
    c,start,end=region
    plt.rcParams.update({'font.family':args.font_family,'font.size':9,'pdf.fonttype':42})
    fig,axes=plt.subplots(10,1,figsize=(13,11),sharex=True,
                         gridspec_kw={'height_ratios':[.75,1.5,1.3,.45,.45,.45,.4,.55,.45,.4],'hspace':.12})
    visible=[e for e in edges if e['left']<end and e['right']>start]
    alpha=min(.28,max(.025,100/max(len(visible),1)))
    helper['draw_oldstyle_arc_track'](axes[0],visible,start,end,'#2E8DC8',alpha,.65,arc_points=36)
    axes[0].set_ylabel('agLINK')
    axes[0].set_title(f'{c}:{start/1e6:g}-{end/1e6:g} Mb | adaptive territory + parent-local 45% multiscale valleys',loc='left',fontsize=12)
    x=(np.arange(len(sig['raw']))+.5)*args.bin_bp
    keep=(x>=start)&(x<=end)
    for y,color,name in zip(sig['smooth'],['#96CEDD','#5FB6BD','#E89B6C'],['25 kb','105 kb','225 kb']):
        axes[1].plot(x[keep],y[keep],color=color,lw=1.1,label=name)
    axes[1].plot(x[keep],sig['background'][keep],color='#777777',lw=.8,ls='--',label='Local background')
    axes[1].set_ylabel('Span\ncoverage')
    axes[1].legend(ncol=4,loc='upper left',frameon=False,fontsize=8)
    axes[2].plot(x[keep],sig['smooth'][1,keep],color='#5FB6BD',lw=1.1)
    for gap in terr['gaps']:
        if gap['start']<end and gap['end']>start:axes[2].axvspan(gap['start'],gap['end'],color='#C1C0BF',alpha=.25)
    cs=[r for r in results['Adaptive']['candidates'] if r['chrom']==c and start<=r['boundary_pos']<=end]
    unique={}
    for r in cs:
        pos=r['boundary_pos']
        if pos not in unique or r['accepted']>unique[pos]['accepted']:unique[pos]=r
    for pos,row in unique.items():
        axes[2].plot(pos,np.interp(pos,x,sig['smooth'][1]),marker='|' if row['accepted'] else 'v',ms=7 if row['accepted'] else 4,
                     color='#287F9C' if row['accepted'] else '#999999',linestyle='None')
    # A small, deterministic subset is labeled; every evaluation is in the TSV.
    chosen=sorted(unique.values(),key=lambda r:(-r['prominence'],r['boundary_pos']))[:4]
    for j,row in enumerate(sorted(chosen,key=lambda r:r['boundary_pos'])):
        reason='accepted' if row['accepted'] else row['reject_reason'].split(',')[0].replace('_',' ')
        axes[2].text(.02+j*.245,.95,f'{row["boundary_pos"]/1e6:.3f} Mb\n{reason}',transform=axes[2].transAxes,
                     va='top',fontsize=7.3,color='#222222')
    axes[2].set_ylabel('Valley / gap\ndiagnostics')
    axes[2].legend(handles=[Line2D([0],[0],marker='v',color='#999999',ls='',label='Rejected valley'),
                            Line2D([0],[0],marker='|',color='#287F9C',ls='',label='Accepted cut')],
                   loc='lower left',ncol=2,fontsize=7.5,frameon=False)
    rowspec=[('Old_V1',axes[3]),('Old_V2',axes[4]),('Old_V3',axes[5]),('Territory',axes[6]),('Adaptive',axes[7]),('M3',axes[8])]
    for method,ax in rowspec:
        arr=terr['territories'] if method=='Territory' else results[method]['crus']
        cuts=set() if method=='Territory' else {r['boundary_pos'] for r in results[method]['boundaries'] if r['chrom']==c}
        for row in arr:
            if row['chrom']!=c or row['start']>=end or row['end']<=start:continue
            lo,hi=max(start,row['start']),min(end,row['end'])
            color='#B8DCCD' if method=='Territory' else COLORS[method]
            ax.add_patch(Rectangle((lo,.18),hi-lo,.42,facecolor=color,edgecolor='#555555',lw=.85))
            for pos in [row['start'],row['end']]:
                if start<=pos<=end:ax.vlines(pos,.65,.91,color=color if pos in cuts else '#888888',lw=1.1,ls='-' if pos in cuts else '--')
        ax.set_ylim(0,1);ax.set_yticks([])
        ax.set_ylabel({'Old_V1':'Old V1','Old_V2':'Old V2','Old_V3':'Old V3','Adaptive':'New CRU'}.get(method,method))
    for i,row in enumerate(blocks.get(c,[])):
        if row['start']<end and row['end']>start:
            lo,hi=max(start,row['start']),min(end,row['end'])
            axes[9].add_patch(Rectangle((lo,0),hi-lo,1,facecolor=['#E6E6E6','#CDCDCD'][i%2],edgecolor='#777777',lw=.8))
    axes[9].set_ylim(0,1);axes[9].set_yticks([]);axes[9].set_ylabel('Block')
    for i,ax in enumerate(axes):
        ax.set_xlim(start,end);ax.grid(False);ax.tick_params(labelsize=8.5,width=.8)
        for spine in ax.spines.values():spine.set_linewidth(.9);spine.set_visible(True);spine.set_color('#222222')
        if i in [1,2]:
            high=max(ax.get_ylim()[1],1);ax.set_ylim(-.06*high,high*(1.6 if i==2 else 1.3))
        if i<9:ax.tick_params(axis='x',bottom=False,labelbottom=False)
    axes[9].xaxis.set_major_locator(MaxNLocator(7));axes[9].xaxis.set_major_formatter(FuncFormatter(lambda v,_:f'{v/1e6:g}'))
    axes[9].set_xlabel(f'{c} position (Mb)')
    fig.subplots_adjust(left=.105,right=.99,top=.96,bottom=.06)
    folder=phase/'f07_figures';folder.mkdir(exist_ok=True)
    path=folder/f'{c}_{start}_{end}.comparison'
    fig.savefig(str(path)+'.pdf');fig.savefig(str(path)+'.png',dpi=160);plt.close(fig)


def summaries(results,edges,sizes,args,pilot,phase):
    rows,lengthrows=[],[]
    for method,res in results.items():
        folder=phase/FOLDERS[method];folder.mkdir(exist_ok=True)
        (folder/'01_parameters.json').write_text(json.dumps(config(args),indent=2)+'\n')
        with (folder/'02_CRU.bed').open('w') as f:
            for r in sorted(res['crus'],key=lambda r:(r['chrom'],r['start'])):f.write(f'{r["chrom"]}\t{r["start"]}\t{r["end"]}\t{r["cru_id"]}\n')
        write(folder/'03_internal_boundaries.tsv',res['boundaries'],None if res['boundaries'] else ['chrom','boundary_pos'])
        write(folder/'04_candidate_valleys.tsv',res.get('candidates',[]),None if res.get('candidates') else ['chrom','boundary_pos'])
        write(folder/'08_recursive_actions.tsv',res.get('actions',[]),None if res.get('actions') else ['chrom','action'])
        stat=[]
        for c in args.chroms:stat.append(pilot.interval_summary(method,c,[r for r in res['crus'] if r['chrom']==c],edges[c],sizes[c]))
        lengths=np.array([r['end']-r['start'] for r in res['crus']])
        pooled={'method':method,'chrom':'Tested_chromosomes','chromosome_bp':sum(sizes[c] for c in args.chroms)}
        for k in ['cru_n','covered_bp','aglink_n','intra_aglink_n','inter_aglink_n','one_endpoint_aglink_n','outside_aglink_n','length_ge2mb_n','length_ge5mb_n']:pooled[k]=sum(r[k] for r in stat)
        pooled['coverage_fraction']=pooled['covered_bp']/pooled['chromosome_bp']
        pooled['intra_aglink_fraction']=pooled['intra_aglink_n']/pooled['aglink_n']
        for label,q in [('p25',.25),('median',.5),('p75',.75),('p90',.9),('p95',.95)]:pooled[label+'_bp']=float(np.quantile(lengths,q)) if len(lengths) else ''
        pooled.update({'mean_bp':float(np.mean(lengths)) if len(lengths) else '', 'max_bp':int(max(lengths)) if len(lengths) else ''})
        stat.append(pooled);rows.extend(stat);write(folder/'05_summary.tsv',stat)
        lengthrows.extend({'method':method,**r,'length_bp':r['end']-r['start']} for r in res['crus'])
    write(phase/'f06_comparison/01_CRU_summary.tsv',rows)
    write(phase/'f06_comparison/06_CRU_lengths.tsv',lengthrows)
    return rows


def parser():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--source-dir',type=Path,required=True)
    ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--pilot-script',type=Path,default=Path(str(Path(__file__).resolve().parent / 'm01_density_signal.py')))
    ap.add_argument('--stage',choices=['preview','chromosomes'],required=True)
    ap.add_argument('--preview-approved',action='store_true')
    ap.add_argument('--chroms',default='1A,1B,2A')
    ap.add_argument('--regions',default='1A:34000000-36000000,1A:538500000-542000000,1B:164375000-166400000,2A:765000000-771000000')
    ap.add_argument('--bin-bp',type=int,default=5000)
    ap.add_argument('--smooth-bp',default='25000,105000,225000')
    ap.add_argument('--background-tile-bp',type=int,default=100000)
    ap.add_argument('--background-window-bp',type=int,default=4100000)
    ap.add_argument('--background-sigma',type=float,default=1)
    ap.add_argument('--min-seed-depth',type=float,default=3)
    ap.add_argument('--seed-endpoints',type=float,default=2)
    ap.add_argument('--gap-min-bp',type=int,default=50000)
    ap.add_argument('--gap-fraction',type=float,default=.2)
    ap.add_argument('--gap-flank-bp',type=int,default=500000)
    ap.add_argument('--bridge-min-links',type=float,default=4)
    ap.add_argument('--min-cru-bp',type=int,default=100000)
    ap.add_argument('--min-internal',type=int,default=2)
    ap.add_argument('--shoulder-bp',type=int,default=500000)
    ap.add_argument('--valley-depth',type=float,default=.45)
    ap.add_argument('--valley-min-width-bp',type=int,default=25000)
    ap.add_argument('--min-prominence',type=float,default=1)
    ap.add_argument('--min-peak-prominence',type=float,default=.5)
    ap.add_argument('--min-relative-peak-prominence',type=float,default=.2)
    ap.add_argument('--min-scale-support',type=int,default=2)
    ap.add_argument('--cross-flank-bp',type=int,default=250000)
    ap.add_argument('--support-flanks-bp',default='250000,500000')
    ap.add_argument('--max-local-cross-fraction',type=float,default=.5)
    ap.add_argument('--max-recursion-depth',type=int,default=16)
    ap.add_argument('--block-perm-n',type=int,default=1000)
    ap.add_argument('--block-seed',type=int,default=20260805)
    ap.add_argument('--block-local-window-bp',type=int,default=10000000)
    ap.add_argument('--font-family',default='Liberation Sans')
    return ap


def main():
    args=parser().parse_args();args.chroms=args.chroms.split(',');args.smooth_bp=[int(v) for v in args.smooth_bp.split(',')]
    args.support_flanks_bp=[int(v) for v in args.support_flanks_bp.split(',')]
    if args.support_flanks_bp!=sorted(set(args.support_flanks_bp)):raise ValueError('Support flanks must be positive and increasing')
    if len(args.chroms)>3 or not set(args.chroms)<={'1A','1B','2A'}:raise ValueError('Pilot restricted to 1A,1B,2A')
    if args.bin_bp!=5000 or args.smooth_bp!=[25000,105000,225000]:raise ValueError('Cached profiles require original 5-kb grid and three scales')
    args.regions=[(r.split(':')[0],*map(int,r.split(':')[1].split('-'))) for r in args.regions.split(',')]
    pilot=module(args.pilot_script)
    args.outdir.mkdir(parents=True,exist_ok=True)
    prepare(args,pilot)
    phase=args.outdir/('f02_preview' if args.stage=='preview' else 'f03_chromosomes')
    if (phase/'01_run_started.json').exists():raise ValueError('Stage output already exists; use an unused output directory.')
    if args.stage=='chromosomes' and (not args.preview_approved or not (args.outdir/'f02_preview/02_PREVIEW.DONE').exists()):
        raise ValueError('Inspect the four preview PDFs first, then explicitly use --preview-approved.')
    phase.mkdir(exist_ok=True);(phase/'f06_comparison').mkdir(exist_ok=True)
    (phase/'01_run_started.json').write_text(json.dumps({'stage':args.stage,'rules':config(args),'source_code_sha256':sha(__file__),'chroms':args.chroms,'regions':args.regions},indent=2)+'\n')
    source=json.loads((args.source_dir/'01_parameters.json').read_text());sizes=source['sizes']
    edges=load_edges(args.source_dir/'f01_shared_input/01_count_only_aglink.annotated.bedpe',args.chroms)
    sigs,terr={},{}
    result={'crus':[],'boundaries':[],'candidates':[],'actions':[],'territories':[]}
    for c in args.chroms:
        with np.load(args.outdir/f'f01_shared/{c}.signals.npz') as data:sigs[c]={k:data[k] for k in data.files}
        ts,acts,gaps=territories(c,sizes[c],edges[c],sigs[c],args)
        terr[c]={'territories':ts,'actions':acts,'gaps':gaps}
        selected=ts if args.stage=='chromosomes' else [t for t in ts if any(cc==c and t['start']<e and t['end']>s for cc,s,e in args.regions)]
        for t in selected:
            called=split_territory(t,edges[c],sigs[c],args)
            for key in result:result[key].extend(called[key])
        pilot.log(f'{args.stage} {c}: territories={len(ts)}; selected complete territories={len(selected)}; current new CRUs={len(result["crus"])}')
    for c in args.chroms:
        for i,r in enumerate(sorted((r for r in result['crus'] if r['chrom']==c),key=lambda r:r['start']),1):r['cru_id']=f'{c}_ADAPTIVE45_CRU_{i}'
    refs,source=references(args,pilot)
    results={'Adaptive':result,**refs}
    old_diagnosis(args,sigs,refs,terr,result,phase)
    write(phase/'f06_comparison/08_new_territory_actions.tsv',[r for c in args.chroms for r in terr[c]['actions']])
    write(phase/'f06_comparison/09_sustained_gap_evidence.tsv',[r for c in args.chroms for r in terr[c]['gaps']])
    for method,res in results.items():
        folder=phase/FOLDERS[method];folder.mkdir(exist_ok=True)
        write(folder/'03_internal_boundaries.tsv',res['boundaries'],None if res['boundaries'] else ['chrom','boundary_pos'])
        write(folder/'04_candidate_valleys.tsv',res.get('candidates',[]),None if res.get('candidates') else ['chrom','boundary_pos'])
        with (folder/'02_CRU.bed').open('w') as f:
            for r in res['crus']:f.write(f'{r["chrom"]}\t{r["start"]}\t{r["end"]}\t{r["cru_id"]}\n')
    helper=runpy.run_path(source['m3_script'],run_name='m3_diagnostics_and_arcs_only')
    blocks=helper['load_block_bed'](source['block_bed'],args.chroms)
    if args.stage=='chromosomes':
        preview=pilot.read_bed(args.outdir/'f02_preview/f01_new_adaptive/02_CRU.bed',args.chroms)
        actual={(r['chrom'],r['start'],r['end']) for r in result['crus']}
        expected={(r['chrom'],r['start'],r['end']) for rows in preview.values() for r in rows}
        if not expected<=actual:raise AssertionError('Preview/full coordinate mismatch')
        summaries(results,edges,sizes,args,pilot,phase)
        pilot.METHODS=FOLDERS
        args.block_bed=source['block_bed']
        old_out=args.outdir;args.outdir=phase
        pilot.block_evaluation(results,sizes,args.chroms,args,helper,phase/'f06_comparison')
        args.outdir=old_out
    for region in args.regions:
        render(region,sigs[region[0]],edges[region[0]],results,terr[region[0]],blocks,phase,args,helper)
        pilot.log('Rendered '+str(region))
    freeze(args,phase/'f08_scripts')
    marker='02_PREVIEW.DONE' if args.stage=='preview' else '02_CHROMOSOMES.DONE'
    (phase/marker).write_text('Only the requested pilot scope was called.\n')
    pilot.log('Complete: '+str(phase))


if __name__=='__main__':main()
