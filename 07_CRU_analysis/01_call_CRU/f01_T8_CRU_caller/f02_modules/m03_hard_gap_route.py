# m03_hard_gap_route: hard-gap B route
# Frozen rule module, copied unchanged from the published calling chain
# (source: f04_hard_gap_pilot/01_run_B_route.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Independent hard-gap-first pilot. Blocks and fixed reference calls are evaluation only."""
import argparse
import copy
import csv
from collections import defaultdict
import gzip
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import runpy
import shutil
from types import SimpleNamespace

import numpy as np


CHROMS=['1A','1B','2A']
REGIONS=[('1B',0,5000000),('1A',5000000,10000000),('1A',10000000,15000000),
         ('1B',5000000,10000000),('1B',164375000,166400000),('1A',538500000,542000000),
         ('2A',765000000,771000000),('1A',34000000,36000000)]
FOLDERS={'B1_45':'f01_B1_45','B2_40':'f02_B2_40','B2_30':'f03_B2_30',
         'Adaptive_fixed':'f04_Adaptive_fixed','M3_fixed':'f05_M3_fixed'}


def module(path,name):
    spec=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m


def md5(path):
    h=hashlib.md5()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1<<20),b''):h.update(b)
    return h.hexdigest()


def rows(path):
    op=gzip.open if str(path).endswith('.gz') else open
    with op(path,'rt') as f:return list(csv.DictReader(f,delimiter='\t'))


def write(path,data,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True);data=list(data)
    fields=fields or list(dict.fromkeys(k for r in data for k in r)) or ['chrom']
    with path.open('w') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)


def savebed(path,data):
    with Path(path).open('w') as f:
        for r in sorted(data,key=lambda r:(CHROMS.index(r['chrom']),r['start'])):
            f.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['cru_id']}\n")


def hard_signals(sig,args,adaptive):
    raw=sig['raw'];fine=sig['smooth'][0]
    left,right=adaptive.side_maxima(sig['smooth'][1],args.hard_reference_bp//args.bin_bp)
    weak=np.minimum(left,right)
    bg=sig['background']+np.sqrt(sig['background']+1)
    limit=np.maximum(args.absolute_low_depth,np.minimum(bg,args.hard_core_relative*weak))
    absolute=(raw<=args.absolute_low_depth+1e-9)&(fine<=args.absolute_low_depth+1+1e-9)
    relative=(raw<=limit+1e-9)&(fine<=limit+1e-9)&(weak>=args.min_reference_depth)
    mask=absolute|relative
    closed=mask.copy()
    # Tolerate only short internal spikes, never fill a wide high-signal island.
    for a,b in adaptive.runs(~mask):
        if a>0 and b<len(mask) and (b-a)*args.bin_bp<=args.max_hole_bp:closed[a:b]=True
    return {'left_reference_B':left,'right_reference_B':right,'weak_reference_B':weak,
            'hard_core_limit':limit,'hard_absolute':absolute,'hard_relative':relative,
            'hard_core':mask,'hard_core_closed':closed}


def identify_gaps(chrom,size,sig,edge,args,adaptive):
    hs=hard_signals(sig,args,adaptive);candidates=[];grown=[]
    raw=sig['raw'];fine=sig['smooth'][0];medium=sig['smooth'][1]
    flank=args.hard_reference_bp//args.bin_bp
    for a,b in adaptive.runs(hs['hard_core_closed']):
        duration=(b-a)*args.bin_bp
        if duration<args.hard_gap_min_bp:
            candidates.append({'chrom':chrom,'core_start':a*args.bin_bp,'core_end':min(size,b*args.bin_bp),
                               'accepted':0,'reason':'core_duration_below_minimum','core_bp':duration});continue
        l=max(0,a-flank);r=min(len(raw),b+flank)
        lp=float(np.max(medium[l:a])) if a>l else 0.
        rp=float(np.max(medium[b:r])) if r>b else 0.
        weak=min(lp,rp)
        # Width detection is multiscale, but final edges grow on RAW bins.
        # No peak is searched inside a previously excluded parent.
        threshold=max(args.absolute_low_depth,args.hard_growth_relative*weak)
        lo,hi=a,b
        while lo>0 and raw[lo-1]<=threshold+1e-9:lo-=1
        while hi<len(raw) and raw[hi]<=threshold+1e-9:hi+=1
        row={'chrom':chrom,'core_start':a*args.bin_bp,'core_end':min(size,b*args.bin_bp),
             'start':lo*args.bin_bp,'end':min(size,hi*args.bin_bp),'core_bp':duration,
             'accepted':1,'reason':'sustained_hard_gap_no_internal_side_veto','left_reference':lp,'right_reference':rp,
             'raw_edge_threshold':threshold,'short_spike_bins':int(np.count_nonzero(~hs['hard_core'][a:b]))}
        candidates.append(row);grown.append((lo,hi))
    merged=[]
    for a,b in sorted(grown):
        if merged and a<=merged[-1][1]:merged[-1][1]=max(b,merged[-1][1])
        else:merged.append([a,b])
    excluded=np.zeros(len(raw),bool);gaps=[]
    left=np.array([e['left'] for e in edge]);right=np.array([e['right'] for e in edge])
    for i,(a,b) in enumerate(merged,1):
        excluded[a:b]=True;start,end=a*args.bin_bp,min(size,b*args.bin_bp)
        center=(start+end)//2
        gaps.append({'chrom':chrom,'start':start,'end':end,'gap_id':f'{chrom}_HARD_GAP_{i}',
             'length_bp':end-start,'median_raw':float(np.median(raw[a:b])),'p90_raw':float(np.quantile(raw[a:b],.9)),
             'max_raw':float(np.max(raw[a:b])),'mean_endpoint_225k':float(np.mean(sig['endpoint'][a:b])),
             'fully_internal_edges':int(np.count_nonzero((left>=start)&(right<end))),
             'edges_spanning_entire_gap':int(np.count_nonzero((left<start)&(right>=end))),
             'crossing_center':int(np.count_nonzero((left<center)&(right>center))),
             'low_le5_bin_fraction':float(np.mean(raw[a:b]<=5+1e-9))})
    return hs,excluded,gaps,candidates


def base_components(chrom,size,mask,sig,edge,args,adaptive):
    left=np.array([e['left'] for e in edge]);right=np.array([e['right'] for e in edge]);kept=[];actions=[]
    for a,b in adaptive.runs(~mask):
        positive=np.flatnonzero(sig['raw'][a:b]>1e-9)
        if not len(positive):continue
        lo,hi=a+int(positive[0]),a+int(positive[-1])+1
        start,end=lo*args.bin_bp,min(size,hi*args.bin_bp)
        support=adaptive.internal_n(left,right,start,end)
        reason='kept' if end-start>=args.min_cru_bp and support>=args.min_internal else 'short_or_insufficient_internal'
        r={'chrom':chrom,'start':start,'end':end,'territory_id':f'{chrom}_HARD_COMPLEMENT_{len(kept)+1}',
           'internal_n':support,'reason':reason,'passed':int(reason=='kept'),'source':'whole_chromosome_minus_hard_gaps_no_Adaptive_territory'}
        actions.append(r)
        if r['passed']:kept.append(r)
    return kept,actions


def recursive(chrom,start,end,edge,sig,options,adaptive,stage,long_min=0):
    left=np.array([e['left'] for e in edge]);right=np.array([e['right'] for e in edge]);todo=[(start,end,0)]
    result={'crus':[],'boundaries':[],'candidates':[],'actions':[]}
    while todo:
        a,b,depth=todo.pop()
        cs=adaptive.parent_candidates(chrom,a,b,left,right,sig,options,depth) if b-a>=long_min else []
        for r in cs:r.update({'stage':stage,'original_parent_start':start,'original_parent_end':end})
        result['candidates'].extend(cs);valid=[r for r in cs if r['eligible']]
        if valid and depth<options.max_recursion_depth:
            best=max(valid,key=lambda r:(r['scale_support_n'],r['score'],r['prominence'],-r['boundary_pos']))
            best['accepted']=1;result['boundaries'].append(dict(best));p=best['boundary_pos']
            for r in valid:
                if r is not best:r['reject_reason']='eligible_not_selected_in_this_parent'
            todo.extend([(p,b,depth+1),(a,p,depth+1)]);action='accepted_soft_valley'
        else:
            result['crus'].append({'chrom':chrom,'start':a,'end':b,'original_parent_start':start,'original_parent_end':end})
            action='kept_no_qualified_valley' if b-a>=long_min else 'secondary_parent_below_length_trigger'
        result['actions'].append({'chrom':chrom,'start':a,'end':b,'stage':stage,'depth':depth,'candidate_n':len(cs),
                                  'accepted':int(action=='accepted_soft_valley'),'action':action})
    return result


def call_chromosome(chrom,size,sig,edge,args,adaptive,rule,regions=None):
    hs,mask,gaps,candidates=identify_gaps(chrom,size,sig,edge,args,adaptive)
    components,component_actions=base_components(chrom,size,mask,sig,edge,args,adaptive)
    selected=components if regions is None else [r for r in components if any(c==chrom and r['start']<e and r['end']>s for c,s,e in regions)]
    primary=SimpleNamespace(**rule);primary.valley_depth=args.primary_depth
    b1={'crus':[],'boundaries':[],'candidates':[],'actions':[]}
    for r in selected:
        child=recursive(chrom,r['start'],r['end'],edge,sig,primary,adaptive,'primary45')
        for k in b1:b1[k].extend(child[k])
    output={'B1_45':b1}
    for name,depth in [('B2_40',.4),('B2_30',.3)]:
        opts=SimpleNamespace(**rule);opts.valley_depth=depth
        method={'crus':[],'boundaries':copy.deepcopy(b1['boundaries']),'candidates':copy.deepcopy(b1['candidates']),
                'actions':copy.deepcopy(b1['actions'])}
        for r in b1['crus']:
            child=recursive(chrom,r['start'],r['end'],edge,sig,opts,adaptive,'secondary'+str(int(depth*100)),args.secondary_parent_min_bp)
            for k in child:method[k].extend(child[k])
        output[name]=method
    for name,res in output.items():
        for i,r in enumerate(sorted(res['crus'],key=lambda r:r['start']),1):r['cru_id']=f'{chrom}_{name}_CRU_{i}'
        for r in res['crus']:
            assert not any(r['start']<g['end'] and r['end']>g['start'] for g in gaps)
        if name!='B1_45':assert {r['boundary_pos'] for r in b1['boundaries']}<={r['boundary_pos'] for r in res['boundaries']}
    return output,{'signals':hs,'mask':mask,'gaps':gaps,'candidates':candidates,'components':component_actions}


def area(mask,start,end,size,step=5000,prefix=None):
    if prefix is None:
        lengths=np.minimum(step,size-np.arange(len(mask))*step)
        prefix=np.r_[0,np.cumsum(mask*lengths)]
    def f(p):
        i=min(len(mask),p//step)
        return int(prefix[i]+((p-i*step)*mask[i] if i<len(mask) else 0))
    return f(end)-f(start)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--history-root',type=Path,required=True);ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--stage',choices=['diagnose','preview','chromosomes'],required=True)
    ap.add_argument('--preview-approved',action='store_true')
    ap.add_argument('--hard-gap-min-bp',type=int,default=100000);ap.add_argument('--hard-reference-bp',type=int,default=2000000)
    ap.add_argument('--absolute-low-depth',type=float,default=2);ap.add_argument('--hard-core-relative',type=float,default=.1)
    ap.add_argument('--hard-growth-relative',type=float,default=.2);ap.add_argument('--max-hole-bp',type=int,default=10000)
    ap.add_argument('--min-reference-depth',type=float,default=10);ap.add_argument('--primary-depth',type=float,default=.45)
    ap.add_argument('--secondary-parent-min-bp',type=int,default=1500000)
    ap.add_argument('--block-perm-n',type=int,default=1000);ap.add_argument('--block-seed',type=int,default=20260805)
    ap.add_argument('--block-local-window-bp',type=int,default=10000000)
    args=ap.parse_args();args.bin_bp=5000;args.min_cru_bp=100000;args.min_internal=2
    base=args.history_root;source=base/'f08_span_only_control';fixed=base/'f10_adaptive_territory_v2'
    out=args.outdir;out.mkdir(parents=True,exist_ok=True)
    settings={k:v for k,v in vars(args).items() if k not in ['stage','preview_approved','block_perm_n','block_seed','block_local_window_bp']}
    settings={k:str(v) if isinstance(v,Path) else v for k,v in settings.items()}
    settings['secondary_depths']=[.4,.3]
    signature=out/'01_parameters.json'
    if signature.exists():assert json.loads(signature.read_text())==settings,'Rules differ from this run'
    else:signature.write_text(json.dumps(settings,indent=2)+'\n')
    phase=out/{'diagnose':'f02_diagnosis','preview':'f03_preview','chromosomes':'f04_chromosomes'}[args.stage]
    if phase.exists():raise ValueError('Stage already exists; do not overwrite')
    if args.stage=='chromosomes':assert args.preview_approved and (out/'f03_preview/02_STAGE.DONE').exists()
    phase.mkdir();(out/'f01_inputs').mkdir(exist_ok=True);(out/'f05_scripts').mkdir(exist_ok=True)
    execution={'code_md5':md5(__file__),'stage':args.stage,'block_perm_n':args.block_perm_n,
               'block_seed':args.block_seed,'block_local_window_bp':args.block_local_window_bp}
    if args.stage=='chromosomes':
        assert json.loads((out/'f03_preview/08_execution_parameters.json').read_text())['code_md5']==execution['code_md5']
    (phase/'08_execution_parameters.json').write_text(json.dumps(execution,indent=2)+'\n')
    here=Path(__file__).parent
    for p in here.iterdir():
        if p.is_file():shutil.copy2(p,out/'f05_scripts'/p.name)
    adaptive_path=Path(str(Path(__file__).resolve().parent / 'm02_adaptive_background.py'))
    pilot_path=Path(str(Path(__file__).resolve().parent / 'm01_density_signal.py'))
    adaptive=module(adaptive_path,'adaptive_scoring_only');pilot=module(pilot_path,'historical_evaluation_only')
    params=json.loads((source/'01_parameters.json').read_text());sizes=params['sizes']
    rule=json.loads((fixed/'f01_shared/01_rule_signature.json').read_text())
    expected={r['path']:r['md5'] for r in rows(base/'f11_reference_provenance_check/f02_tables/01_input_inventory_MD5.tsv')}
    input_paths=[source/'f01_shared_input/01_count_only_aglink.annotated.bedpe',adaptive_path,pilot_path,
        fixed/'f03_chromosomes/f01_new_adaptive/02_CRU.bed',Path(params['m3_dir'])/'M3_CRU.bed',Path(params['block_bed'])]
    input_paths += [fixed/f'f01_shared/{c}.signals.npz' for c in CHROMS]
    before=[]
    for p in input_paths:
        value=md5(p)
        if str(p) in expected:assert value==expected[str(p)],str(p)
        before.append({'path':str(p),'md5_before':value})
    edges=adaptive.load_edges(source/'f01_shared_input/01_count_only_aglink.annotated.bedpe',CHROMS)
    sigs={};aux={};results={k:{'crus':[],'boundaries':[],'candidates':[],'actions':[]} for k in ['B1_45','B2_40','B2_30']}
    diagnostics=[]
    for c in CHROMS:
        with np.load(fixed/f'f01_shared/{c}.signals.npz') as z:sigs[c]={k:z[k] for k in z.files}
        s=sigs[c];recomputed=pilot.span_coverage(edges[c],sizes[c],5000)
        assert np.allclose(recomputed,s['raw'],atol=1e-8,rtol=1e-12),'Raw cache mismatch'
        endpoint=np.zeros(len(s['raw']),int)
        for e in edges[c]:endpoint[e['p1']//5000]+=1;endpoint[e['p2']//5000]+=1
        centers=(np.arange(len(s['raw']))+.5)*5000
        crossing=np.searchsorted(sorted(e['left'] for e in edges[c]),centers,side='left')-np.searchsorted(sorted(e['right'] for e in edges[c]),centers,side='right')
        hs=hard_signals(s,args,adaptive)
        for cc,start,end in REGIONS:
            if cc!=c:continue
            for i in range(start//5000,min(len(endpoint),math.ceil(end/5000))):
                diagnostics.append({'region':f'{c}:{start}-{end}','chrom':c,'start':i*5000,'end':min(sizes[c],(i+1)*5000),
                    'raw_span':s['raw'][i],'endpoint_bin_n':endpoint[i],'endpoint_225k':s['endpoint'][i],
                    'aglink_crossing_bin_center':int(crossing[i]),'local_background':s['background'][i],
                    'old_seed_threshold':s['seed_threshold'][i],'old_left_ref_500k':s['left_reference'][i],
                    'old_right_ref_500k':s['right_reference'][i],'old_low_zero':int(s['low_zero'][i]),
                    'old_low_relative':int(s['low_relative'][i]),'old_low_thin':int(s['low_thin'][i]),'old_blocked':int(s['blocked'][i]),
                    **{k:float(v[i]) for k,v in hs.items()},**{f'smooth_{w}k':s['smooth'][j,i] for j,w in enumerate([25,105,225])}})
        if args.stage=='diagnose':continue
        called,a=call_chromosome(c,sizes[c],s,edges[c],args,adaptive,rule,REGIONS if args.stage=='preview' else None);aux[c]=a
        for method,res in called.items():
            for k in res:results[method][k].extend(res[k])
        pilot.log(f'{args.stage} {c}: hard gaps={len(a["gaps"])}; B1={len(called["B1_45"]["crus"])}; B2_40={len(called["B2_40"]["crus"])}; B2_30={len(called["B2_30"]["crus"])}')
    write(phase/'01_regional_bin_diagnosis.tsv',diagnostics)
    reasons=[]
    for c,start,end in REGIONS:
        rs=[r for r in diagnostics if r['region']==f'{c}:{start}-{end}']
        reasons.append({'region':f'{c}:{start}-{end}','raw_min':min(r['raw_span'] for r in rs),
            'raw_median':float(np.median([r['raw_span'] for r in rs])),'raw_p90':float(np.quantile([r['raw_span'] for r in rs],.9)),
            'old_gap_bin_n':sum(r['old_blocked'] for r in rs),'new_core_bin_n':sum(r['hard_core'] for r in rs),
            'reason':'Old hard floor <4 plus 500kb shoulders can lie inside a wide nonzero valley; see per-bin flags, no cache mismatch'})
    write(phase/'03_regional_cause_summary.tsv',reasons)
    if args.stage=='diagnose':
        for r in before:r['md5_after']=md5(r['path']);r['unchanged']=int(r['md5_before']==r['md5_after']);assert r['unchanged']
        write(phase/'07_input_MD5_unchanged.tsv',before)
        (phase/'02_STAGE.DONE').write_text('Diagnosis only; raw profile independently verified.\n');return
    for name,folder,bedpath,cutpaths in [
        ('Adaptive_fixed','f04_Adaptive_fixed',fixed/'f03_chromosomes/f01_new_adaptive/02_CRU.bed',[fixed/'f03_chromosomes/f01_new_adaptive/03_internal_boundaries.tsv']),
        ('M3_fixed','f05_M3_fixed',Path(params['m3_dir'])/'M3_CRU.bed',[Path(params['m3_dir'])/'M3_primary_boundaries.bed',Path(params['m3_dir'])/'M3_secondary_boundaries.bed'])]:
        by=pilot.read_bed(bedpath,CHROMS);bounds=[]
        for p in cutpaths:
            if str(p).endswith('.tsv'):bounds += [{**r,'boundary_pos':int(r['boundary_pos'])} for r in rows(p)]
            else:
                for arr in pilot.read_bed(p,CHROMS).values():bounds += [{**r,'boundary_pos':(r['start']+r['end'])//2} for r in arr]
        results[name]={'crus':[r for c in CHROMS for r in by[c]],'boundaries':bounds,'candidates':[],'actions':[],
                       'reference_bed':str(bedpath),'reference_md5':md5(bedpath)}
    write(phase/'04_hard_gaps.tsv',[r for c in CHROMS for r in aux[c]['gaps']])
    write(phase/'05_hard_gap_candidates.tsv',[r for c in CHROMS for r in aux[c]['candidates']])
    write(phase/'06_hard_complement_actions.tsv',[r for c in CHROMS for r in aux[c]['components']])
    for c in CHROMS:
        short=pilot.span_coverage(edges[c],sizes[c],5000,lambda e:e['distance_bp']<=250000)
        np.savez_compressed(out/f'f01_inputs/{c}.B_signals.npz',**sigs[c],**aux[c]['signals'],hard_excluded=aux[c]['mask'],short_le250k=short)
    summaries=[];valleys=[]
    for method,res in results.items():
        folder=phase/FOLDERS[method];folder.mkdir()
        meta={'method':method,'construction_rules':settings} if method.startswith('B') else {'fixed_reference':True,'bed':res['reference_bed'],'md5':res['reference_md5']}
        (folder/'01_parameters.json').write_text(json.dumps(meta,indent=2)+'\n')
        savebed(folder/'02_CRU.bed',res['crus'])
        write(folder/'03_internal_boundaries.tsv',res['boundaries'],None if res['boundaries'] else ['chrom','boundary_pos'])
        write(folder/'04_candidate_valleys.tsv',res['candidates'],None if res['candidates'] else ['chrom','boundary_pos'])
        write(folder/'08_recursive_actions.tsv',res['actions'])
        if args.stage=='chromosomes':
            lengths=[];stats=[]
            for c in CHROMS:
                rr=[r for r in res['crus'] if r['chrom']==c];stat=pilot.interval_summary(method,c,rr,edges[c],sizes[c]);lengths.extend(r['end']-r['start'] for r in rr)
                for key,mask in [('raw_le2',sigs[c]['raw']<=2+1e-9),('raw_le5',sigs[c]['raw']<=5+1e-9),
                                 ('below_local_background',sigs[c]['raw']<=sigs[c]['background']+np.sqrt(sigs[c]['background']+1)),('in_B_hard_gap',aux[c]['mask'])]:
                    lens=np.minimum(5000,sizes[c]-np.arange(len(mask))*5000)
                    prefix=np.r_[0,np.cumsum(mask*lens)]
                    bp=sum(area(mask,r['start'],r['end'],sizes[c],prefix=prefix) for r in rr);stat[key+'_bp']=bp
                    stat[key+'_fraction']=bp/stat['covered_bp'] if stat['covered_bp'] else 0
                stat['length_ge10mb_n']=sum(r['end']-r['start']>=10000000 for r in rr);stats.append(stat)
            pooled={'method':method,'chrom':'Tested_chromosomes','chromosome_bp':sum(sizes[c] for c in CHROMS)}
            summed=['cru_n','covered_bp','aglink_n','intra_aglink_n','inter_aglink_n','one_endpoint_aglink_n','outside_aglink_n',
                    'length_ge2mb_n','length_ge5mb_n','length_ge10mb_n','raw_le2_bp','raw_le5_bp','below_local_background_bp','in_B_hard_gap_bp']
            for k in summed:pooled[k]=sum(r[k] for r in stats)
            pooled.update({'coverage_fraction':pooled['covered_bp']/pooled['chromosome_bp'],
                'intra_aglink_fraction':pooled['intra_aglink_n']/pooled['aglink_n'],'max_bp':max(lengths),'mean_bp':float(np.mean(lengths))})
            for label,q in [('p25',.25),('median',.5),('p75',.75),('p90',.9),('p95',.95)]:pooled[label+'_bp']=float(np.quantile(lengths,q))
            for k in ['raw_le2','raw_le5','below_local_background','in_B_hard_gap']:pooled[k+'_fraction']=pooled[k+'_bp']/pooled['covered_bp']
            stats.append(pooled);summaries.extend(stats);write(folder/'05_summary.tsv',stats)
        pos={(r['chrom'],r['boundary_pos']) for r in res['candidates']};accepted={(r['chrom'],r['boundary_pos']) for r in res['boundaries']}
        for c in CHROMS+['Tested_chromosomes']:
            choose=lambda r:c=='Tested_chromosomes' or r['chrom']==c
            valleys.append({'method':method,'chrom':c,'candidate_evaluation_n':sum(choose(r) for r in res['candidates']),
                'candidate_unique_n':sum(c=='Tested_chromosomes' or cc==c for cc,p in pos),
                'accepted_soft_cut_n':sum(choose(r) for r in res['boundaries']),
                'rejected_only_unique_n':sum(c=='Tested_chromosomes' or cc==c for cc,p in pos-accepted),
                'secondary_extra_cut_n':sum(choose(r) and r.get('stage','').startswith('secondary') for r in res['boundaries']),
                **{f'cuts_in_original_parent_ge{n}Mb':sum(choose(r) and int(r.get('original_parent_end',0))-int(r.get('original_parent_start',0))>=n*1000000 for r in res['boundaries']) for n in [2,5,10]}})
    comp=phase/'f06_comparison';comp.mkdir()
    write(comp/'07_soft_valley_counts.tsv',valleys)
    if args.stage=='chromosomes':
        write(comp/'01_CRU_summary.tsv',summaries)
        for method in ['B1_45','B2_40','B2_30']:
            preview=pilot.read_bed(out/'f03_preview'/FOLDERS[method]/'02_CRU.bed',CHROMS)
            actual={(r['chrom'],r['start'],r['end']) for r in results[method]['crus']}
            assert {(r['chrom'],r['start'],r['end']) for cc in CHROMS for r in preview[cc]}<=actual
        # Construction has finished. Only now read blocks and run the unchanged evaluation.
        helper=runpy.run_path(params['m3_script'],run_name='M3_evaluation_only')
        pilot.METHODS=FOLDERS;args.chroms=CHROMS;args.block_bed=params['block_bed']
        saved=args.outdir;args.outdir=phase;pilot.block_evaluation(results,sizes,CHROMS,args,helper,comp);args.outdir=saved
    for r in before:r['md5_after']=md5(r['path']);r['unchanged']=int(r['md5_before']==r['md5_after']);assert r['unchanged']
    write(phase/'07_input_MD5_unchanged.tsv',before)
    (phase/'02_STAGE.DONE').write_text('Fixed 1A/1B/2A scope; calling and reference provenance checks passed.\n')
    pilot.log('Complete '+str(phase))


if __name__=='__main__':main()
