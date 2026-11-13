# m05_density_transitions: agLINK density transitions
# Frozen rule module, copied unchanged from the published calling chain
# (source: f09_aglink_density_transitions/01_call_transitions.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Additional sharp/step candidates on frozen A-edge+B1-45; no coACR or blocks."""
import argparse
from collections import defaultdict
import csv
import importlib.util
import json
import math
from pathlib import Path
import shutil
from types import SimpleNamespace
import numpy as np
from scipy.signal import find_peaks

def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

from . import m03_hard_gap_route as C
from . import m01_density_signal as P
CHROMS=['1A','1B','2A'];STEP=5000
METHODS={'Baseline':'f01_baseline','V1_sharp':'f02_V1_sharp','V2_step':'f03_V2_step','V3_both':'f04_V3_both'}
REGIONS=[('1A',5000000,10000000,7500000),('1A',535000000,540000000,538400000),
 ('1B',20000000,25000000,23000000),('1B',15000000,20000000,18600000),
 ('1B',164000000,167000000,164400000),('1A',34000000,36000000,35000000),
 ('1B',50000000,59000000,55000000),('2A',765000000,771000000,766000000)]


def bounds(start,end,n):
    return max(0,math.ceil(start/STEP-.5)),min(n,math.ceil(end/STEP-.5))


def proposals(sig,start,end,args):
    fine=sig['smooth'][0];a,b=bounds(start,end,len(fine));out=[]
    if b-a<5:return out
    idx,_=find_peaks(-fine[a:b],prominence=args.min_absolute_change,distance=10)
    out.extend(('sharp',int(i+a)) for i in idx)
    w=args.persistence_bp//STEP;skip=args.transition_exclusion_bp//STEP
    ids=np.arange(a+w+skip,b-w-skip-1)
    if len(ids)<3:return out
    prefix=np.r_[0,np.cumsum(fine)]
    left=(prefix[ids-skip]-prefix[ids-skip-w])/w
    right=(prefix[ids+skip+1+w]-prefix[ids+skip+1])/w
    contrast=np.abs(left-right)/np.maximum(np.maximum(left,right),1e-9)
    peaks,_=find_peaks(contrast,height=.20,prominence=.03,distance=10)
    grad=np.gradient(fine);chosen={}
    for j in peaks:
        i=int(ids[j]);direction=1 if right[j]>left[j] else -1
        lo=max(a+1,i-15);hi=min(b-1,i+16)
        if lo>=hi:continue
        aligned=lo+int(np.argmax(direction*grad[lo:hi]))
        if aligned not in chosen or contrast[j]>chosen[aligned]:chosen[aligned]=float(contrast[j])
    out.extend(('step',i) for i in sorted(chosen))
    return sorted(set(out),key=lambda t:(t[1],t[0]))


def morph(source,i,sig,start,end,args):
    a,b=bounds(start,end,len(sig['raw']));pos=i*STEP+STEP//2
    row={'source':source,'boundary_pos':pos,'parent_start':start,'parent_end':end,
         'local_background':float(sig['background'][i]),'raw_span':float(sig['raw'][i])}
    reasons=[];support=0
    if not a<i<b-1:return {**row,'shape_pass':0,'shape_reason':'parent_edge','strength':0}
    if source=='sharp':
        flank=args.sharp_recovery_bp//STEP
        lo=max(a,i-flank);hi=min(b,i+flank+1)
        li=lo+int(np.argmax(sig['smooth'][0,lo:i]));ri=i+1+int(np.argmax(sig['smooth'][0,i+1:hi]))
        depths=[]
        for k,label,threshold in zip(range(3),['25k','105k','225k'],[args.sharp_depth,.20,.10]):
            y=sig['smooth'][k];floor=float(y[i]);lp=float(y[li]);rp=float(y[ri])
            depth=min((lp-floor)/max(lp,1e-9),(rp-floor)/max(rp,1e-9));depths.append(depth)
            support+=int(depth>=threshold)
            row.update({label+'_floor':floor,label+'_left_level':lp,label+'_right_level':rp,label+'_relative_change':depth})
        fine=sig['smooth'][0];floor=float(fine[i]);weak=float(min(fine[li],fine[ri]));level=floor+.5*max(0,weak-floor)
        l,r=i,i
        while l>li and fine[l-1]<=level:l-=1
        while r<ri and fine[r+1]<=level:r+=1
        width=(r-l+1)*STEP
        delta=weak-floor
        row.update(direction='valley',left_peak_pos=li*STEP+STEP//2,right_peak_pos=ri*STEP+STEP//2,
            transition_width_bp=width,absolute_change=delta,floor_above_background=int(floor>sig['background'][i]),
            recovery_left_bp=(i-li)*STEP,recovery_right_bp=(ri-i)*STEP)
        if depths[0]<args.sharp_depth:reasons.append('sharp_bilateral_fine_depth')
        if not STEP<=width<=args.sharp_max_width_bp:reasons.append('not_a_sharp_width')
        if delta<args.min_absolute_change:reasons.append('small_absolute_change')
        if support<2:reasons.append('sharp_multiscale_support')
        strength=max(0,depths[0])*support
    else:
        w=args.persistence_bp//STEP;skip=args.transition_exclusion_bp//STEP
        if i-skip-w<a or i+skip+1+w>b:
            return {**row,'shape_pass':0,'shape_reason':'insufficient_sustained_flanks','strength':0}
        depths=[];direction=None
        for k,label,threshold in zip(range(3),['25k','105k','225k'],[args.step_relative,.20,.10]):
            y=sig['smooth'][k];ls=y[i-skip-w:i-skip];rs=y[i+skip+1:i+skip+1+w]
            lv=float(np.median(ls));rv=float(np.median(rs));sign=1 if rv>lv else -1
            if direction is None:direction=sign
            rel=abs(lv-rv)/max(lv,rv,1e-9);depths.append(rel)
            support+=int(sign==direction and rel>=threshold)
            row.update({label+'_left_level':lv,label+'_right_level':rv,label+'_relative_change':rel})
        lv,rv=row['25k_left_level'],row['25k_right_level'];delta=abs(lv-rv)
        y=sig['smooth'][0];ls=y[i-skip-w:i-skip];rs=y[i+skip+1:i+skip+1+w]
        blocks=4;left_blocks=np.array([np.median(z) for z in np.array_split(ls,blocks)])
        right_blocks=np.array([np.median(z) for z in np.array_split(rs,blocks)])
        persist=int(np.count_nonzero(direction*(right_blocks-left_blocks)>=.20*max(lv,rv)))
        high,low=max(lv,rv),min(lv,rv);width=None
        search=args.max_transition_bp//STEP
        pre=np.arange(max(a,i-search),i+1);post=np.arange(i,min(b,i+search+1))
        if direction<0:
            l=pre[y[pre]>=high-.20*delta];r=post[y[post]<=low+.20*delta]
        else:
            l=pre[y[pre]<=low+.20*delta];r=post[y[post]>=high-.20*delta]
        if len(l) and len(r):width=max(STEP,int(r[0]-l[-1])*STEP)
        row.update(direction='rise' if direction>0 else 'drop',absolute_change=delta,
            transition_width_bp=width if width is not None else '',persistence_blocks=blocks,
            sustained_block_n=persist,low_high_ratio=(low+1)/(high+1),
            change_per_kb=delta/(max(STEP,width or STEP)/1000))
        if depths[0]<args.step_relative:reasons.append('step_relative_change')
        if delta<args.min_absolute_change:reasons.append('small_absolute_change')
        if support<2:reasons.append('step_multiscale_direction_or_magnitude')
        if persist<3:reasons.append('not_sustained_step')
        if width is None or width>args.max_transition_bp:reasons.append('transition_not_fast')
        strength=depths[0]*support
    row.update(scale_support_n=support,shape_pass=int(not reasons),shape_reason=','.join(reasons),strength=strength)
    return row


def unique_support(es,mask):
    selected=[e for e,m in zip(es,mask) if m]
    return len({(e['acr_id'],e['gene_id']) for e in selected}),len({e['acr_id'] for e in selected}),len({e['gene_id'] for e in selected})


def graph(es,start,end,pos,args):
    left=np.array([e['left'] for e in es]);right=np.array([e['right'] for e in es])
    wholeL=unique_support(es,(left>=start)&(right<pos));wholeR=unique_support(es,(left>=pos)&(right<end))
    crossing=(left<pos)&(right>=pos);global_cross=unique_support(es,crossing)[0]
    contexts=[]
    for flank in [250000,500000]:
        a,b=max(start,pos-flank),min(end,pos+flank)
        il=unique_support(es,(left>=a)&(right<pos));ir=unique_support(es,(left>=pos)&(right<b))
        local=unique_support(es,crossing&(left>=a)&(right<b))[0]
        touch=unique_support(es,crossing&(((left>=a)&(left<b))|((right>=a)&(right<b))))[0]
        frac=touch/(touch+min(il[0],ir[0])) if touch+min(il[0],ir[0]) else 1.
        contexts.append({'support_context_bp':flank,'left_internal':il[0],'right_internal':ir[0],
            'left_ACR_n':il[1],'right_ACR_n':ir[1],'left_gene_n':il[2],'right_gene_n':ir[2],
            'crossing_fully_local':local,'crossing_touching_context':touch,'local_cross_fraction':frac})
    good=[c for c in contexts if min(c['left_internal'],c['right_internal'])>=args.min_local_internal and
          min(c['left_ACR_n'],c['right_ACR_n'])>=2 and c['local_cross_fraction']<=args.max_cross_fraction]
    chosen=min(good,key=lambda c:(c['support_context_bp'],c['local_cross_fraction'])) if good else contexts[-1]
    reasons=[]
    if pos-start<args.min_daughter_bp or end-pos<args.min_daughter_bp:reasons.append('near_fixed_or_new_parent_edge')
    if min(wholeL[0],wholeR[0])<args.min_daughter_internal:reasons.append('weak_or_edge_like_daughter')
    if min(wholeL[1],wholeR[1])<2:reasons.append('daughter_ACR_diversity')
    if min(chosen['left_internal'],chosen['right_internal'])<args.min_local_internal:reasons.append('local_internal_support')
    if min(chosen['left_ACR_n'],chosen['right_ACR_n'])<2:reasons.append('local_ACR_diversity')
    if chosen['local_cross_fraction']>args.max_cross_fraction:reasons.append('touching_crossing_fraction')
    loss=global_cross/max(1,global_cross+wholeL[0]+wholeR[0])
    if loss>args.max_global_cut_loss:reasons.append('global_cut_loss')
    result={**chosen,'whole_left_internal':wholeL[0],'whole_right_internal':wholeR[0],
        'whole_left_ACR_n':wholeL[1],'whole_right_ACR_n':wholeR[1],
        'crossing_global':global_cross,'global_cut_loss':loss,'graph_pass':int(not reasons),
        'graph_reason':','.join(reasons)}
    for c in contexts:
        result.update({f'context_{c["support_context_bp"]}_{k}':v for k,v in c.items()})
    return result


def evaluate(source,i,sig,parent,es,args):
    row=morph(source,i,sig,parent['start'],parent['end'],args)
    row.update(graph(es,parent['start'],parent['end'],row['boundary_pos'],args))
    reason=','.join(z for z in [row['shape_reason'],row['graph_reason']] if z)
    row.update(chrom=parent['chrom'],baseline_id=parent['baseline_id'],eligible=int(not reason),
        accepted=0,reject_reason=reason,classification='internal_candidate' if not row['graph_reason'] else 'unsupported_or_edge_like')
    return row


def call_parent(parent,sig,es,sources,args,method):
    todo=[(parent['start'],parent['end'],0)];rows=[];cuts=[];leaves=[]
    while todo:
        start,end,depth=todo.pop();child={**parent,'start':start,'end':end}
        candidates=[]
        for source,i in proposals(sig,start,end,args):
            if source not in sources:continue
            row=evaluate(source,i,sig,child,es,args);row.update(method=method,recursive_depth=depth)
            candidates.append(row)
        valid=[r for r in candidates if r['eligible']]
        if valid and depth<16:
            best=max(valid,key=lambda r:(r['strength']/(1+r['local_cross_fraction']),r['scale_support_n'],-r['boundary_pos'],r['source']))
            best['accepted']=1;best['reject_reason']='accepted_density_and_connection_evidence';cuts.append(dict(best))
            for row in valid:
                if row is not best:row['reject_reason']='eligible_rechecked_in_child'
            p=best['boundary_pos'];todo.extend([(p,end,depth+1),(start,p,depth+1)])
        else:
            if valid:
                for row in valid:row['reject_reason']='recursion_safety_limit'
            leaves.append({**child,'recursive_depth':depth})
        rows.extend(candidates)
    return leaves,cuts,rows


def main(args):
    root=args.outdir;baseline=args.baseline_dir
    if root.exists():raise ValueError('Output exists; never overwrite history')
    meta=json.loads((baseline/'f01_baseline/01_parameters.json').read_text());hist=Path(meta['history_root'])
    bed=baseline/'f01_baseline/02_CRU.bed';assert C.md5(bed)==(baseline/'f01_baseline/12_BASELINE.DONE').read_text().strip()
    old=P.read_bed(bed,CHROMS);assert sum(map(len,old.values()))==365
    edges=defaultdict(list)
    for r in C.rows(meta['aglink_bedpe']):
        if r['chrom'] not in CHROMS:continue
        assert float(r['weight'])==1
        for key in ['p1','p2','left','right']:r[key]=int(r[key])
        edges[r['chrom']].append(r)
    root.mkdir()
    for d in ['f01_results','f02_plot_data','f03_evaluation','f04_figures','f05_scripts','f06_checks']:(root/d).mkdir()
    original_cuts=C.rows(baseline/'f01_baseline/03_internal_boundaries.tsv')
    old_candidates=C.rows(baseline/'f01_baseline/04_candidate_valleys.tsv')
    sources=[bed,Path(meta['aglink_bedpe']),Path(meta['A_gap_file']),baseline/'f01_baseline/03_internal_boundaries.tsv',
             baseline/'f01_baseline/04_candidate_valleys.tsv',baseline/'f01_baseline/01_parameters.json',Path(__file__)]
    sources += [hist/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz' for c in CHROMS]
    hashes=[{'path':str(p),'md5_before':C.md5(p)} for p in sources]
    results={m:{'crus':[],'boundaries':[],'candidates':[]} for m in METHODS}
    inventory=[];signals={};focus=[];region_signals=[];plot_arcs=[]
    for c in CHROMS:
        with np.load(hist/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz') as z:sig={k:z[k] for k in z.files}
        signals[c]=sig
        for r in old[c]:r['baseline_id']=r['cru_id']
        results['Baseline']['crus'].extend(old[c])
        results['Baseline']['boundaries'].extend({**r,'boundary_pos':int(r['boundary_pos'])} for r in original_cuts if r['chrom']==c)
        x=(np.arange(len(sig['raw']))+.5)*STEP
        gl=np.array([e['left'] for e in edges[c]]);gr=np.array([e['right'] for e in edges[c]])
        cross=np.searchsorted(np.sort(gl),x,side='left')-np.searchsorted(np.sort(gr),x,side='left')
        # Inventory is chromosome-wide; acceptance always re-evaluates within the frozen parent.
        for source,i in proposals(sig,0,meta['sizes'][c],args):
            pos=i*STEP+STEP//2
            parent=next((r for r in old[c] if r['start']<=pos<r['end']),None)
            if parent:
                es=[e for e in edges[c] if e['left']<parent['end'] and e['right']>=parent['start']]
                row=evaluate(source,i,sig,parent,es,args)
            else:
                row={**morph(source,i,sig,0,meta['sizes'][c],args),'chrom':c,'baseline_id':'',
                     'eligible':0,'accepted':0,'classification':'outside_frozen_CRU_gap',
                     'reject_reason':'outside_fixed_CRU_no_boundary_change_allowed'}
            inventory.append(row)
        for method,wanted in [('V1_sharp',{'sharp'}),('V2_step',{'step'}),('V3_both',{'sharp','step'})]:
            for parent in old[c]:
                es=[e for e in edges[c] if e['left']<parent['end'] and e['right']>=parent['start']]
                rr,bb,cc=call_parent(parent,sig,es,wanted,args,method)
                results[method]['crus'].extend(rr);results[method]['boundaries'].extend(bb);results[method]['candidates'].extend(cc)
            for r in original_cuts:
                if r['chrom']==c:results[method]['boundaries'].append({**r,'boundary_pos':int(r['boundary_pos']),'source':'inherited_B145','method':method})
        print(c,'complete',[(m,sum(r['chrom']==c for r in res['crus'])) for m,res in results.items()],flush=True)
        for chrom,a,b,target in REGIONS:
            if chrom!=c:continue
            region=f'{c}_{a}_{b}'
            for e in edges[c]:
                if e['left']<b and e['right']>a:plot_arcs.append({'region':region,**{k:e[k] for k in ['chrom','left','right','acr_id','gene_id']}})
            for i in np.flatnonzero((x>=a-STEP)&(x<=b+STEP)):
                region_signals.append({'region':region,'chrom':c,'position':x[i],'raw_span':sig['raw'][i],
                    'span_25k':sig['smooth'][0,i],'span_105k':sig['smooth'][1,i],'span_225k':sig['smooth'][2,i],
                    'local_background':sig['background'][i],'ag_cross_global':cross[i]})
            near=sorted((r for r in inventory if r['chrom']==c and abs(r['boundary_pos']-target)<=150000),
                        key=lambda r:abs(r['boundary_pos']-target))[:6]
            previous=sorted((r for r in old_candidates if r['chrom']==c),key=lambda r:abs(int(r['boundary_pos'])-target))[:2]
            for r in near:focus.append({'region':region,'focus_position':target,'record_type':'new_inventory',**r})
            for r in previous:focus.append({'region':region,'focus_position':target,'record_type':'saved_B145_evaluation',**r})
    for method,res in results.items():
        folder=root/'f01_results'/METHODS[method];folder.mkdir()
        if method!='Baseline':
            for c in CHROMS:
                rr=sorted((r for r in res['crus'] if r['chrom']==c),key=lambda r:r['start'])
                for i,r in enumerate(rr,1):r['cru_id']=f'{c}_{method}_CRU_{i}'
        C.savebed(folder/'02_CRU.bed',res['crus'])
        C.write(folder/'03_internal_boundaries.tsv',res['boundaries'])
        C.write(folder/'04_candidate_evaluations.tsv',res['candidates'],None if res['candidates'] else ['chrom','boundary_pos'])
        (folder/'01_parameters.json').write_text(json.dumps({'method':method,'baseline_BED':str(bed),
            'baseline_md5':C.md5(bed),'baseline_internal_cuts_locked':True,'outer_edges_locked':True,
            'no_coACR':True,'rules':vars(args)},default=str,indent=2)+'\n')
        for c in CHROMS:
            rr=sorted((r for r in res['crus'] if r['chrom']==c),key=lambda r:r['start'])
            assert all(a['end']<=b['start'] for a,b in zip(rr,rr[1:]))
            for parent in old[c]:
                child=sorted((r for r in rr if r['baseline_id']==parent['cru_id']),key=lambda r:r['start'])
                assert child[0]['start']==parent['start'] and child[-1]['end']==parent['end']
                assert all(a['end']==b['start'] for a,b in zip(child,child[1:]))
    C.write(root/'f06_checks/02_chromosome_candidate_inventory.tsv',inventory)
    C.write(root/'f06_checks/03_focus_and_original_B145_reasons.tsv',focus)
    # Plot data are exported only after the caller has finished; block positions cannot guide cuts.
    oldparams=json.loads((hist/'f08_span_only_control/01_parameters.json').read_text())
    block_path=Path(meta['block_bed']);m3_path=Path(oldparams['m3_dir'])/'M3_CRU.bed'
    blocks=P.read_bed(block_path,CHROMS);m3=P.read_bed(m3_path,CHROMS)
    hashes += [{'path':str(p),'md5_before':C.md5(p)} for p in [block_path,m3_path]]
    pd=root/'f02_plot_data'
    C.write(pd/'01_regions.tsv',[{'region':f'{c}_{a}_{b}','chrom':c,'start':a,'end':b,'focus':t} for c,a,b,t in REGIONS])
    C.write(pd/'02_aglink_arcs.tsv.gz',plot_arcs)
    C.write(pd/'03_signals.tsv',region_signals)
    ints=[];bbs=[];cand=[];gaps=[];boundaries=[]
    frozen_gaps=C.rows(meta['A_gap_file'])
    for c,a,b,t in REGIONS:
        region=f'{c}_{a}_{b}'
        for method,res in results.items():
            for r in res['crus']:
                if r['chrom']==c and r['start']<b and r['end']>a:ints.append({'region':region,'method':method,**r})
            for r in res['candidates']:
                if r['chrom']==c and a<=r['boundary_pos']<b:cand.append({'region':region,**r})
            for r in res['boundaries']:
                if r['chrom']==c and a<=int(r['boundary_pos'])<b:boundaries.append({'region':region,'method':method,**r})
        for r in m3[c]:
            if r['start']<b and r['end']>a:ints.append({'region':region,'method':'M3_reference',**r})
        for r in blocks[c]:
            if r['start']<b and r['end']>a:bbs.append({'region':region,**r})
        for r in frozen_gaps:
            if r['chrom']==c and int(r['start'])<b and int(r['end'])>a:gaps.append({'region':region,**r})
    C.write(pd/'04_intervals.tsv',ints);C.write(pd/'05_blocks.tsv',bbs)
    C.write(pd/'06_candidates.tsv',cand);C.write(pd/'07_frozen_gaps.tsv',gaps);C.write(pd/'08_boundaries.tsv',boundaries)
    for h in hashes:
        h['md5_after']=C.md5(h['path']);assert h['md5_before']==h['md5_after']
    C.write(root/'f06_checks/01_input_MD5.tsv',hashes)
    config={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}
    config.update(sizes=meta['sizes'],block_bed=meta['block_bed'],m3_script=oldparams['m3_script'],
        aglink_bedpe=meta['aglink_bedpe'],baseline_MD5=C.md5(bed),no_coACR=True,block_positions_used_for_calling=False)
    (root/'01_parameters.json').write_text(json.dumps(config,indent=2)+'\n')
    for p in Path(__file__).parent.iterdir():
        if p.is_file():shutil.copy2(p,root/'f05_scripts'/p.name)
    (root/'03_CALLING.DONE').write_text('Only 1A/1B/2A; inherited boundaries and coverage locked; no coACR.\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--baseline-dir',type=Path,required=True);ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--sharp-depth',type=float,default=.40);ap.add_argument('--sharp-recovery-bp',type=int,default=250000)
    ap.add_argument('--sharp-max-width-bp',type=int,default=150000);ap.add_argument('--step-relative',type=float,default=.35)
    ap.add_argument('--persistence-bp',type=int,default=200000);ap.add_argument('--transition-exclusion-bp',type=int,default=25000)
    ap.add_argument('--max-transition-bp',type=int,default=100000);ap.add_argument('--min-absolute-change',type=float,default=2)
    ap.add_argument('--min-daughter-bp',type=int,default=150000);ap.add_argument('--min-local-internal',type=int,default=2)
    ap.add_argument('--min-daughter-internal',type=int,default=3);ap.add_argument('--max-cross-fraction',type=float,default=.50)
    ap.add_argument('--max-global-cut-loss',type=float,default=.40)
    main(ap.parse_args())
