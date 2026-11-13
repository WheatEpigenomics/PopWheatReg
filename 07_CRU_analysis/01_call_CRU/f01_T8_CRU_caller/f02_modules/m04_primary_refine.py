# m04_primary_refine: B1/B2 primary refinement
# Frozen rule module, copied unchanged from the published calling chain
# (source: f05_B_route_refinement/01_refine_B_route.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Factorial pilot on frozen B2-30: supported gap-edge recovery and conservative long-CRU cuts."""
import argparse
import copy
from collections import defaultdict
import importlib.util
import json
import math
from pathlib import Path
import runpy
import shutil
from types import SimpleNamespace

import numpy as np


def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

from . import m03_hard_gap_route as C
from . import m02_adaptive_background as D
from . import m01_density_signal as P
CHROMS=C.CHROMS
REGIONS=[('1B',0,5000000),('1B',15000000,20000000),('1B',20000000,25000000),
         ('1B',50000000,59000000),('1B',164375000,166400000),('2A',765000000,771000000),
         ('2A',773000000,787782082),('1A',5000000,10000000),('1A',34000000,36000000),
         ('1A',538500000,542000000)]
FOLDERS={'A_edge':'f01_A_edge','B_long':'f02_B_long','AB':'f03_AB',
         'B1_fixed':'f04_B1_fixed','B230_fixed':'f05_B230_fixed'}


def node_support(es):
    pairs={(e['acr_id'],e['gene_id']) for e in es}
    return {'unique_pair_n':len(pairs),'acr_n':len({e['acr_id'] for e in es}),
            'gene_n':len({e['gene_id'] for e in es}),'acr_position_n':len({e['p1'] for e in es}),
            'TSS_n':len({e['p2'] for e in es}),'edge_n':len(es)}


def internal(es,start,end):
    return [e for e in es if start<=e['left'] and e['right']<end]


def close_short(mask,width):
    out=mask.copy()
    for a,b in D.runs(~mask):
        if a>0 and b<len(mask) and b-a<=width:out[a:b]=True
    return out


def recover_band(sig,es,start,end,side,core_floor,args):
    a,b=start//5000,end//5000
    if b-a<args.recovery_bp//5000:
        return None,{'reason':'growth_band_shorter_than_recovery_duration','band_start':start,'band_end':end}
    threshold=np.maximum.reduce([sig['background'][a:b]+args.recovery_sigma*np.sqrt(sig['background'][a:b]+1),
                                 np.full(b-a,args.core_recovery_fold*core_floor),np.full(b-a,4.)])
    raw=sig['raw'][a:b];scales=sig['smooth'][:,a:b]
    supported=(raw>=threshold-1e-9)&(np.count_nonzero(scales>=threshold,axis=0)>=2)
    connected=close_short(supported,args.recovery_hole_bp//5000)
    runs=[(x,y) for x,y in D.runs(connected) if (x==0 if side=='left' else y==b-a)]
    if not runs:return None,{'reason':'no_density_recovery_connected_to_kept_CRU','band_start':start,'band_end':end}
    x,y=runs[0];duration=(y-x)*5000
    if duration<args.recovery_bp:
        return None,{'reason':'connected_recovery_too_short','recovery_bp':duration,'band_start':start,'band_end':end}
    window=args.recovery_bp//5000
    possible=range(y-window,x-1,-1) if side=='left' else range(x,y-window+1)
    last_support=None
    local=internal(es,start,end)
    for i in possible:
        lo,hi=(a+i)*5000,(a+i+window)*5000
        support=node_support(internal(local,lo,hi));last_support=support
        if support['unique_pair_n']<args.edge_min_pairs or min(support['acr_n'],support['gene_n'],support['acr_position_n'],support['TSS_n'])<2:continue
        position=hi if side=='left' else lo
        return position,{'reason':'sustained_density_plus_independent_internal_support','band_start':start,'band_end':end,
            'window_start':lo,'window_end':hi,'edge_position':position,'core_floor_q90':core_floor,
            'recovery_level_median':float(np.median(threshold[i:i+window])),
            'raw_median':float(np.median(raw[i:i+window])),**support}
    return None,{'reason':'density_recovered_but_no_independent_internal_window','band_start':start,'band_end':end,
                 **(last_support or {})}


def refine_edges(chrom,sig,es,gaps,cores,args):
    refined=[];events=[]
    for old in gaps:
        start,end=int(old['start']),int(old['end'])
        inside=[r for r in cores if int(r['accepted']) and int(r['core_start'])>=start and int(r['core_end'])<=end]
        assert inside,'Every original merged gap needs its frozen accepted cores'
        ca=min(int(r['core_start']) for r in inside);cb=max(int(r['core_end']) for r in inside)
        ix=np.concatenate([np.arange(int(r['core_start'])//5000,int(r['core_end'])//5000) for r in inside])
        floor=float(np.quantile(sig['raw'][np.unique(ix)],.9))
        left,le=recover_band(sig,es,start,ca,'left',floor,args)
        right,re=recover_band(sig,es,cb,end,'right',floor,args)
        a=start if left is None else left;b=end if right is None else right
        assert start<=a<=ca<cb<=b<=end
        refined.append({**old,'start':a,'end':b,'old_start':start,'old_end':end,'fixed_core_start':ca,'fixed_core_end':cb})
        for side,event in [('left',le),('right',re)]:
            events.append({'chrom':chrom,'gap_id':old['gap_id'],'side':side,'old_gap_start':start,'old_gap_end':end,
                'new_gap_start':a,'new_gap_end':b,'fixed_core_start':ca,'fixed_core_end':cb,
                'core_q90':floor,'changed':int((a!=start) if side=='left' else (b!=end)),**event})
    return refined,events


def apply_edges(base,gaps):
    starts={int(g['old_end']):int(g['end']) for g in gaps};ends={int(g['old_start']):int(g['start']) for g in gaps}
    out=[];events=[]
    for r in base:
        a=starts.get(r['start'],r['start']);b=ends.get(r['end'],r['end'])
        out.append({**r,'start':a,'end':b,'baseline_id':r['cru_id']})
        if (a,b)!=(r['start'],r['end']):events.append({'chrom':r['chrom'],'baseline_id':r['cru_id'],
            'old_start':r['start'],'old_end':r['end'],'new_start':a,'new_end':b,'added_bp':b-a-(r['end']-r['start'])})
    assert len(out)==len(base)
    return out,events


def graph_evidence(es,start,end,pos,args):
    local=internal(es,start,end)
    left=[e for e in local if e['p1']<pos and e['p2']<pos]
    right=[e for e in local if e['p1']>=pos and e['p2']>=pos]
    cross=[e for e in local if (e['p1']<pos)!=(e['p2']<pos)]
    L,R,X=node_support(left),node_support(right),node_support(cross)
    il,ir,c=L['unique_pair_n'],R['unique_pair_n'],X['unique_pair_n']
    mixing=c/(il+ir+c) if il+ir+c else 1.
    possible_left=L['acr_n']*L['gene_n'];possible_right=R['acr_n']*R['gene_n']
    possible_cross=L['acr_n']*R['gene_n']+R['acr_n']*L['gene_n']
    dl=il/possible_left if possible_left else 0.;dr=ir/possible_right if possible_right else 0.
    dc=c/possible_cross if possible_cross else float('inf')
    ratio=dc/min(dl,dr) if min(dl,dr)>0 else float('inf')
    reasons=[]
    if min(il,ir)<args.long_min_pairs:reasons.append('insufficient_independent_pairs')
    if min(L['acr_n'],L['gene_n'],R['acr_n'],R['gene_n'],L['acr_position_n'],L['TSS_n'],R['acr_position_n'],R['TSS_n'])<2:
        reasons.append('insufficient_independent_ACR_gene_modules')
    if mixing>args.max_total_cross_fraction:reasons.append('cross_fraction_exceeds_limit')
    if ratio>args.max_cross_density_ratio:reasons.append('cross_density_not_depleted')
    return {**{'graph_left_'+k:v for k,v in L.items()},**{'graph_right_'+k:v for k,v in R.items()},
        'graph_cross_n':c,'graph_cross_fraction':mixing,'graph_left_density':dl,'graph_right_density':dr,
        'graph_cross_density':dc,'graph_cross_density_ratio':ratio,'graph_reason':','.join(reasons),'graph_pass':int(not reasons)}


def long_scan(row,es,sig,rule,args,do_split):
    todo=[(row['start'],row['end'],0)];candidates=[];cuts=[];leaves=[]
    left=np.array([e['left'] for e in es]);right=np.array([e['right'] for e in es])
    options=SimpleNamespace(**rule);options.valley_depth=args.long_depth
    while todo:
        a,b,depth=todo.pop()
        if b-a<args.long_trigger_bp:
            leaves.append({**row,'start':a,'end':b,'baseline_id':row['cru_id']});continue
        cs=D.parent_candidates(row['chrom'],a,b,left,right,sig,options,depth)
        for r in cs:
            reasons=[]
            if not (r['25k_depth']>=args.long_depth and r['25k_width_bp']>=rule['valley_min_width_bp']):reasons.append('fine_depth_or_width')
            passing=[s for s in ['25k','105k','225k'] if r[s+'_depth']>=args.long_depth and r[s+'_width_bp']>=rule['valley_min_width_bp']]
            if len(passing)<2:reasons.append('multiscale_depth_or_width')
            if min(r['boundary_pos']-a,b-r['boundary_pos'])<rule['min_cru_bp']:reasons.append('daughter_length')
            # Peak-anchored context is fixed a priori, clipped inside THIS parent.
            scales=passing or ['25k','105k','225k']
            lo=max(a,min(r['boundary_pos']-250000,min(r[s+'_left_peak_pos'] for s in scales)-args.peak_padding_bp))
            hi=min(b,max(r['boundary_pos']+250000,max(r[s+'_right_peak_pos'] for s in scales)+args.peak_padding_bp))
            ev=graph_evidence(es,lo,hi,r['boundary_pos'],args)
            if not ev['graph_pass']:reasons += ev['graph_reason'].split(',')
            r.update(ev);r.update({'baseline_id':row['cru_id'],'baseline_start':row['start'],'baseline_end':row['end'],
                'graph_context_start':lo,'graph_context_end':hi,'old_reject_reason':r['reject_reason'],
                'eligible':int(not reasons),'accepted':0,'reject_reason':','.join(reasons),'stage':'new_conservative_long'})
        candidates.extend(cs);valid=[r for r in cs if r['eligible']]
        if valid and do_split and depth<rule['max_recursion_depth']:
            best=max(valid,key=lambda r:(r['scale_support_n'],r['score'],r['prominence'],-r['boundary_pos']))
            best['accepted']=1;cuts.append(copy.deepcopy(best));p=best['boundary_pos']
            todo.extend([(p,b,depth+1),(a,p,depth+1)])
            for r in valid:
                if r is not best:r['reject_reason']='eligible_not_selected_in_this_parent'
        else:
            leaves.append({**row,'start':a,'end':b,'baseline_id':row['cru_id']})
    return leaves,cuts,candidates


def split_at_saved(rows,cuts):
    positions=defaultdict(list)
    for r in cuts:positions[r['baseline_id']].append(r['boundary_pos'])
    out=[]
    for r in rows:
        points=[r['start']]+sorted(positions[r['baseline_id']])+[r['end']]
        out += [{**r,'start':a,'end':b} for a,b in zip(points,points[1:])]
    return out


def preview_ids(base,changes,chrom,size):
    views=[(a,b) for c,a,b in REGIONS if c==chrom]
    for r in changes:
        if r['chrom']!=chrom:continue
        for old,new in [('old_start','new_start'),('old_end','new_end')]:
            if r[old]!=r[new]:views.append((max(0,min(r[old],r[new])-500000),min(size,max(r[old],r[new])+500000)))
    return {r['cru_id'] for r in base if any(r['start']<b and r['end']>a for a,b in views)}


def mask_bp(sig,start,end,size):
    return C.area(sig,start,end,size)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--history-root',type=Path,required=True);ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--stage',choices=['diagnose','preview','chromosomes'],required=True)
    ap.add_argument('--preview-approved',action='store_true');ap.add_argument('--combine-approved',action='store_true')
    ap.add_argument('--recovery-bp',type=int,default=100000);ap.add_argument('--recovery-sigma',type=float,default=2.)
    ap.add_argument('--core-recovery-fold',type=float,default=2.);ap.add_argument('--recovery-hole-bp',type=int,default=10000)
    ap.add_argument('--edge-min-pairs',type=int,default=3);ap.add_argument('--long-trigger-bp',type=int,default=5000000)
    ap.add_argument('--long-depth',type=float,default=.45);ap.add_argument('--peak-padding-bp',type=int,default=100000)
    ap.add_argument('--long-min-pairs',type=int,default=5);ap.add_argument('--max-total-cross-fraction',type=float,default=.2)
    ap.add_argument('--max-cross-density-ratio',type=float,default=.5)
    ap.add_argument('--block-perm-n',type=int,default=1000);ap.add_argument('--block-seed',type=int,default=20260805)
    ap.add_argument('--block-local-window-bp',type=int,default=10000000)
    args=ap.parse_args();base=args.history_root;old=base/'f13_B_route_hard_gap_first/f04_chromosomes';out=args.outdir
    settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items() if k not in ['stage','preview_approved','combine_approved','block_perm_n','block_seed','block_local_window_bp']}
    out.mkdir(parents=True,exist_ok=True);signature=out/'01_parameters.json'
    if signature.exists():assert json.loads(signature.read_text())==settings
    else:signature.write_text(json.dumps(settings,indent=2)+'\n')
    phase=out/{'diagnose':'f02_diagnosis','preview':'f03_preview','chromosomes':'f04_chromosomes'}[args.stage]
    if phase.exists():raise ValueError('Stage exists; do not overwrite')
    if args.stage=='chromosomes':assert args.preview_approved and (out/'f03_preview/02_STAGE.DONE').exists()
    phase.mkdir();(out/'f01_inputs').mkdir(exist_ok=True);(out/'f05_scripts').mkdir(exist_ok=True)
    for p in Path(__file__).parent.iterdir():
        if p.is_file():shutil.copy2(p,out/'f05_scripts'/p.name)
    execution={'code_md5':C.md5(__file__),'combine_approved':args.combine_approved,'block_perm_n':args.block_perm_n,'block_seed':args.block_seed}
    if args.stage=='chromosomes':assert json.loads((out/'f03_preview/08_execution.json').read_text())['code_md5']==execution['code_md5']
    (phase/'08_execution.json').write_text(json.dumps(execution,indent=2)+'\n')
    params=json.loads((base/'f08_span_only_control/01_parameters.json').read_text());sizes=params['sizes']
    baseline=P.read_bed(old/'f03_B2_30/02_CRU.bed',CHROMS);b1=P.read_bed(old/'f01_B1_45/02_CRU.bed',CHROMS)
    oldcuts=[{**r,'boundary_pos':int(r['boundary_pos'])} for r in C.rows(old/'f03_B2_30/03_internal_boundaries.tsv')]
    b1cuts=[{**r,'boundary_pos':int(r['boundary_pos'])} for r in C.rows(old/'f01_B1_45/03_internal_boundaries.tsv')]
    edges=D.load_edges(base/'f08_span_only_control/f01_shared_input/01_count_only_aglink.annotated.bedpe',CHROMS)
    rule=json.loads((base/'f10_adaptive_territory_v2/f01_shared/01_rule_signature.json').read_text())
    oldgaps=C.rows(old/'04_hard_gaps.tsv');oldcores=C.rows(old/'05_hard_gap_candidates.tsv')
    protected=[old/'f03_B2_30/02_CRU.bed',old/'f01_B1_45/02_CRU.bed',old/'04_hard_gaps.tsv',old/'05_hard_gap_candidates.tsv',
               old/'f03_B2_30/03_internal_boundaries.tsv',old/'f01_B1_45/03_internal_boundaries.tsv',Path(params['block_bed']),
               base/'f08_span_only_control/f01_shared_input/01_count_only_aglink.annotated.bedpe']
    protected += [base/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz' for c in CHROMS]
    checks=[{'path':str(p),'md5_before':C.md5(p)} for p in protected]
    frozen=out/'f01_inputs/01_frozen_source_manifest.tsv'
    if frozen.exists():
        expected={r['path']:r['md5_before'] for r in C.rows(frozen)}
        assert all(expected[r['path']]==r['md5_before'] for r in checks),'Frozen source changed between stages'
    else:C.write(frozen,checks)
    result={m:{'crus':[],'boundaries':copy.deepcopy(oldcuts),'candidates':[]} for m in ['A_edge','B_long']}
    if args.stage=='chromosomes' and args.combine_approved:result['AB']={'crus':[],'boundaries':copy.deepcopy(oldcuts),'candidates':[]}
    sigs={};agaps=[];events=[];moves=[];audit=[];diagnostic=[]
    for c in CHROMS:
        with np.load(base/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz') as z:sig={k:z[k] for k in z.files}
        sigs[c]=sig;gg=[r for r in oldgaps if r['chrom']==c];cc=[r for r in oldcores if r['chrom']==c]
        newgaps,ee=refine_edges(c,sig,edges[c],gg,cc,args);agaps+=newgaps;events+=ee
        A,mm=apply_edges(baseline[c],newgaps);moves+=mm
        longrows=[];cuts=[];candidates=[]
        for r in baseline[c]:
            rr,bb,cs=long_scan(r,edges[c],sig,rule,args,args.stage!='diagnose')
            longrows+=rr;cuts+=bb;candidates+=cs
        diagnostic+=candidates
        for c0,start,end in REGIONS:
            if c0!=c:continue
            for i in range(start//5000,min(len(sig['raw']),math.ceil(end/5000))):
                audit.append({'region':f'{c}:{start}-{end}','chrom':c,'start':i*5000,'end':min(sizes[c],(i+1)*5000),
                    'raw':sig['raw'][i],'endpoint_225k':sig['endpoint'][i],'background':sig['background'][i],
                    **{f'smooth_{w}k':sig['smooth'][j,i] for j,w in enumerate([25,105,225])}})
        if args.stage=='diagnose':continue
        selected=preview_ids(baseline[c],mm,c,sizes[c])
        def chosen(r):return args.stage=='chromosomes' or r['baseline_id'] in selected
        result['A_edge']['crus'] += [r for r in A if chosen(r)]
        result['B_long']['crus'] += [r for r in longrows if chosen(r)]
        result['B_long']['boundaries'] += cuts;result['B_long']['candidates'] += candidates
        if 'AB' in result:
            result['AB']['crus'] += split_at_saved(A,cuts)
            result['AB']['boundaries'] += copy.deepcopy(cuts);result['AB']['candidates'] += copy.deepcopy(candidates)
    C.write(phase/'03_A_gap_edge_evidence.tsv',events);C.write(phase/'04_A_refined_gaps.tsv',agaps)
    C.write(phase/'05_A_CRU_edge_changes.tsv',moves);C.write(phase/'06_B_candidate_evidence.tsv',diagnostic)
    C.write(phase/'09_regional_bin_profiles.tsv',audit)
    if args.stage=='diagnose':
        for r in checks:r['md5_after']=C.md5(r['path']);assert r['md5_before']==r['md5_after']
        C.write(phase/'07_input_MD5.tsv',checks);(phase/'02_STAGE.DONE').write_text('Diagnosis only; no new BED emitted.\n');return
    for m,by,bc,path in [('B1_fixed',b1,b1cuts,old/'f01_B1_45/02_CRU.bed'),('B230_fixed',baseline,oldcuts,old/'f03_B2_30/02_CRU.bed')]:
        result[m]={'crus':[r for c in CHROMS for r in by[c]],'boundaries':bc,'candidates':[],'source':str(path),'md5':C.md5(path)}
    for method,res in result.items():
        if method in ['A_edge','B_long','AB']:
            for c in CHROMS:
                for i,r in enumerate(sorted((r for r in res['crus'] if r['chrom']==c),key=lambda r:r['start']),1):r['cru_id']=f'{c}_{method}_{i}'
        folder=phase/FOLDERS[method];folder.mkdir()
        C.savebed(folder/'02_CRU.bed',res['crus'])
        C.write(folder/'06_baseline_ID_mapping.tsv',[{'chrom':r['chrom'],'new_id':r['cru_id'],
            'baseline_id':r.get('baseline_id',r['cru_id']),'start':r['start'],'end':r['end']} for r in res['crus']])
        C.write(folder/'03_internal_boundaries.tsv',res['boundaries'])
        C.write(folder/'04_candidate_valleys.tsv',res['candidates'],None if res['candidates'] else ['chrom','boundary_pos'])
        meta={'fixed_reference':True,'source_BED':res['source'],'md5':res['md5']} if method.endswith('fixed') else {'method':method,'settings':settings,'baseline':'f13 B2_30','B_rules_applied':method in ['B_long','AB'],'A_edges_applied':method in ['A_edge','AB']}
        (folder/'01_parameters.json').write_text(json.dumps(meta,indent=2)+'\n')
    comp=phase/'f06_comparison';comp.mkdir();summary=[]
    if args.stage=='chromosomes':
        for method,res in result.items():
            stats=[];lengths=[]
            for c in CHROMS:
                rr=[r for r in res['crus'] if r['chrom']==c];s=P.interval_summary(method,c,rr,edges[c],sizes[c]);lengths += [r['end']-r['start'] for r in rr]
                for key,mask in [('raw_le2',sigs[c]['raw']<=2+1e-9),('raw_le5',sigs[c]['raw']<=5+1e-9)]:
                    lens=np.minimum(5000,sizes[c]-np.arange(len(mask))*5000);prefix=np.r_[0,np.cumsum(mask*lens)]
                    bp=sum(C.area(mask,r['start'],r['end'],sizes[c],prefix=prefix) for r in rr);s[key+'_bp']=bp;s[key+'_fraction']=bp/s['covered_bp']
                s['length_ge10mb_n']=sum(r['end']-r['start']>=10000000 for r in rr);stats.append(s)
            total={'method':method,'chrom':'Tested_chromosomes','chromosome_bp':sum(sizes[c] for c in CHROMS)}
            for k in ['cru_n','covered_bp','aglink_n','intra_aglink_n','inter_aglink_n','one_endpoint_aglink_n','outside_aglink_n','length_ge2mb_n','length_ge5mb_n','length_ge10mb_n','raw_le2_bp','raw_le5_bp']:total[k]=sum(r[k] for r in stats)
            total.update({'coverage_fraction':total['covered_bp']/total['chromosome_bp'],'intra_aglink_fraction':total['intra_aglink_n']/total['aglink_n'],
                          'mean_bp':float(np.mean(lengths)),'max_bp':max(lengths)})
            for label,q in [('p25',.25),('median',.5),('p75',.75),('p90',.9),('p95',.95)]:total[label+'_bp']=float(np.quantile(lengths,q))
            for k in ['raw_le2','raw_le5']:total[k+'_fraction']=total[k+'_bp']/total['covered_bp']
            stats.append(total);summary+=stats;C.write(phase/FOLDERS[method]/'05_summary.tsv',stats)
        C.write(comp/'01_CRU_summary.tsv',summary)
        for method in ['A_edge','B_long']:
            prev=P.read_bed(out/'f03_preview'/FOLDERS[method]/'02_CRU.bed',CHROMS)
            actual={(r['chrom'],r['start'],r['end']) for r in result[method]['crus']}
            assert {(r['chrom'],r['start'],r['end']) for c in CHROMS for r in prev[c]}<=actual
        helper=runpy.run_path(params['m3_script'],run_name='evaluation_only')
        P.METHODS={m:FOLDERS[m] for m in result};args.chroms=CHROMS;args.block_bed=params['block_bed']
        saved=args.outdir;args.outdir=phase;P.block_evaluation(result,sizes,CHROMS,args,helper,comp);args.outdir=saved
    for r in checks:r['md5_after']=C.md5(r['path']);assert r['md5_before']==r['md5_after']
    C.write(phase/'07_input_MD5.tsv',checks)
    (phase/'02_STAGE.DONE').write_text('Only 1A/1B/2A; no fixed baseline modified.\n')
    P.log('Complete '+str(phase))


if __name__=='__main__':main()
