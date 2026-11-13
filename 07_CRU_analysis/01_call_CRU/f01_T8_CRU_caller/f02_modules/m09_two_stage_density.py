# m09_two_stage_density: V8/V9 two-stage density
# Frozen rule module, copied unchanged from the published calling chain
# (source: f13_two_stage_density/01_call_two_stage.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Two-stage density sensitivity from frozen B145; no target counts or block input."""
import argparse
from collections import defaultdict
import importlib.util
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
import numpy as np
from scipy.ndimage import percentile_filter


def load(path, name):
    spec=importlib.util.spec_from_file_location(name,path)
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);return mod


from . import m07_changepoint as N
A,F,C,P=N.A,N.F,N.C,N.P
write=P.write_tsv
CHROMS=['1A','1B','2A']; STEP=5000
METHODS={'Baseline':'f01_baseline','V7_reference':'f02_V7_reference','V8':'f03_V8','V9':'f04_V9'}
NEW_SOURCES={'sharp','step','changepoint'}
REGIONS=[('1A',5000000,10000000,7500000),('1B',15000000,20000000,17700000),
 ('1B',20000000,25000000,23000000),('1A',535000000,540000000,538400000),
 ('1A',34000000,36000000,35200000),('1B',50000000,59000000,55180000),
 ('1B',164000000,167000000,165000000),('2A',778000000,784000000,781000000)]


def tier(args, method):
    p=method.lower()
    return SimpleNamespace(high_depth=getattr(args,p+'_high_depth'), low_depth=getattr(args,p+'_low_depth'),
        high_snr=args.high_snr, low_snr=getattr(args,p+'_low_snr'), linear_gain=getattr(args,p+'_linear_gain'),
        max_touch=getattr(args,p+'_max_touch'), max_short_fraction=.80 if method=='V8' else .85,
        profile_min=.10 if method=='V8' else .05, dominant_snr=4. if method=='V8' else 3.,
        dominant_global_loss=.65 if method=='V8' else .75, max_platform_spread=.35 if method=='V8' else .45,
        max_transition_bp=150000 if method=='V8' else 200000, sharp_recovery_bp=500000,
        sharp_max_width_bp=400000, min_absolute_change=2., min_raw_valley_bp=10000,
        min_daughter_bp=args.min_daughter_bp, min_local_internal=2, min_daughter_internal=3,
        max_cross_fraction=getattr(args,p+'_max_touch'), max_global_cut_loss=1.,
        sharp_depth=.15, step_relative=.15, persistence_bp=200000, transition_exclusion_bp=25000,
        cp_candidate=.15, cp_accept=.15, cp_sse_gain=.20)


def density_states(sig,start,end):
    a,b=F.bounds(start,end,len(sig['raw']))
    if end>(len(sig['raw'])-1)*STEP:b=len(sig['raw'])
    med=sig['smooth'][1,a:b]; broad=sig['smooth'][2,a:b]
    bg=sig['background'][a:b]; width=min(401,max(3,len(med)//2*2+1))
    tm=np.maximum.reduce([percentile_filter(med,60,size=width,mode='nearest'),bg+1.5*np.sqrt(bg+1),np.full(len(bg),3.)])
    tb=np.maximum.reduce([percentile_filter(broad,50,size=width,mode='nearest'),bg+np.sqrt(bg+1),np.full(len(bg),2.)])
    high=(med>=tm-1e-9)&(broad>=tb-1e-9)&(sig['raw'][a:b]>=2)
    state=np.where(high,2,np.where(sig['raw'][a:b]>0,1,0))
    return {'offset':a,'high':high,'state':state,'threshold105':tm,'threshold225':tb}


def regime_at(states,i,flank=100000):
    j=i-states['offset']; w=flank//STEP
    l=states['high'][max(0,j-w):j];r=states['high'][j+1:min(len(states['high']),j+w+1)]
    lf=float(np.mean(l)) if len(l) else 0.;rf=float(np.mean(r)) if len(r) else 0.
    return {'high_density_context':int(max(lf,rf)>=.30),'left_high_fraction':lf,'right_high_fraction':rf,
            'local_high_threshold105':float(states['threshold105'][j]),'local_high_threshold225':float(states['threshold225'][j])}


def basic_row(source,i,sig,start,end,states,opts):
    r={'source':source,'boundary_pos':i*STEP+STEP//2,'parent_start':start,'parent_end':end,
       'raw_span':float(sig['raw'][i]),'local_background':float(sig['background'][i]),**regime_at(states,i)}
    high=r['high_density_context'];r['effective_depth_threshold']=opts.high_depth if high else opts.low_depth
    r['effective_snr_threshold']=opts.high_snr if high else opts.low_snr
    return r


def sharp_shape(i,sig,start,end,states,opts):
    row=basic_row('sharp',i,sig,start,end,states,opts); a,b=F.bounds(start,end,len(sig['raw']))
    flank=opts.sharp_recovery_bp//STEP;lo=max(a,i-flank);hi=min(b,i+flank+1)
    if i<=lo or i+1>=hi:return {**row,'shape_pass':0,'shape_reason':'insufficient_peak_flanks','strength':0,'scale_support_n':0}
    for k,label in enumerate(['25k','105k','225k']):
        y=sig['smooth'][k];li=lo+int(np.argmax(y[lo:i]));ri=i+1+int(np.argmax(y[i+1:hi]))
        lp,rp,floor=map(float,[y[li],y[ri],y[i]])
        depth=min((lp-floor)/max(lp,1.),(rp-floor)/max(rp,1.))
        row.update({label+'_left_level':lp,label+'_right_level':rp,label+'_floor':floor,label+'_relative_change':depth,
                    label+'_left_peak_pos':li*STEP+STEP//2,label+'_right_peak_pos':ri*STEP+STEP//2})
    weak=min(row['25k_left_level'],row['25k_right_level']);floor=row['25k_floor'];delta=weak-floor
    level=floor+.5*max(0,delta);fine=sig['smooth'][0]
    li=(row['25k_left_peak_pos']-STEP//2)//STEP;ri=(row['25k_right_peak_pos']-STEP//2)//STEP
    l=r=i
    while l>li and fine[l-1]<=level:l-=1
    while r<ri and fine[r+1]<=level:r+=1
    width=(r-l+1)*STEP
    raw=sig['raw'];l=r=i
    if raw[i]<=level:
        while l>lo and raw[l-1]<=level:l-=1
        while r+1<hi and raw[r+1]<=level:r+=1
        duration=(r-l+1)*STEP
    else:duration=0
    snr=delta/np.sqrt(max(1.,weak+floor+2*row['local_background']+1))
    reason=[];depth=row['effective_depth_threshold']; medium=max(.05,.25*depth)
    if row['25k_relative_change']<depth:reason.append('valley_depth')
    if row['105k_relative_change']<medium:reason.append('valley_105k_support')
    if not opts.min_raw_valley_bp<=width<=opts.sharp_max_width_bp:reason.append('valley_width')
    if duration<opts.min_raw_valley_bp:reason.append('single_bin_or_short_raw_valley')
    if delta<max(opts.min_absolute_change,.8*np.sqrt(row['local_background']+1)):reason.append('absolute_change_vs_background')
    if snr<row['effective_snr_threshold']:reason.append('low_absolute_density_SNR')
    support=int(row['25k_relative_change']>=depth)+int(row['105k_relative_change']>=medium)+int(row['225k_relative_change']>=.05)
    row.update(direction='valley',absolute_change=delta,density_snr=snr,transition_width_bp=width,
               raw_valley_duration_bp=duration,scale_support_n=support,shape_pass=int(not reason),
               shape_reason=','.join(reason),strength=snr)
    return row


def level_shape(source,i,sig,start,end,states,opts):
    base=basic_row(source,i,sig,start,end,states,opts);a,b=F.bounds(start,end,len(sig['raw']));tested=[]
    for flank in [100000,200000,400000,800000]:
        w,skip=flank//STEP,25000//STEP
        if i-skip-w<a or i+skip+1+w>b:continue
        row=dict(base);segments=[];signs=[]
        for k,label in enumerate(['25k','105k','225k']):
            y=sig['smooth'][k];left=y[i-skip-w:i-skip];right=y[i+skip+1:i+skip+1+w]
            lv,rv=float(np.median(left)),float(np.median(right));signs.append(1 if rv>lv else -1)
            row.update({label+'_left_level':lv,label+'_right_level':rv,label+'_relative_change':abs(rv-lv)/max(lv,rv,1.)})
            segments.append((left,right))
        lv,rv=row['105k_left_level'],row['105k_right_level'];direction=signs[1];high=max(lv,rv,1.);delta=abs(lv-rv)
        left,right=segments[1];lb=np.array([np.median(x) for x in np.array_split(left,4)]);rb=np.array([np.median(x) for x in np.array_split(right,4)])
        depth=row['effective_depth_threshold'];sustain=int(np.count_nonzero(direction*(rb-lb)>=max(.10,.6*depth)*high))
        spread=max(float(np.ptp(lb)),float(np.ptp(rb)))/high
        yy=np.r_[left,right];xx=np.r_[np.arange(i-skip-w,i-skip),np.arange(i+skip+1,i+skip+1+w)].astype(float)
        xx=(xx-i)*STEP/1000.;X=np.column_stack([np.ones(len(xx)),xx]);coef=np.linalg.lstsq(X,yy,rcond=None)[0]
        sse0=float(np.sum((yy-np.mean(yy))**2));sse1=float(np.sum((left-np.mean(left))**2)+np.sum((right-np.mean(right))**2))
        linear=float(np.sum((yy-X@coef)**2));gain=1-sse1/sse0 if sse0>0 else 0.
        advantage=1-sse1/max(linear,1e-12)
        delta25=abs(row['25k_left_level']-row['25k_right_level']);low25=min(row['25k_left_level'],row['25k_right_level']);high25=max(row['25k_left_level'],row['25k_right_level'])
        y=sig['smooth'][0];search=opts.max_transition_bp//STEP
        pre=np.arange(max(a,i-search),i+1);post=np.arange(i,min(b,i+search+1))
        if direction>0:pl=pre[y[pre]<=low25+.2*delta25];pr=post[y[post]>=high25-.2*delta25]
        else:pl=pre[y[pre]>=high25-.2*delta25];pr=post[y[post]<=low25+.2*delta25]
        width=max(STEP,int(pr[0]-pl[-1])*STEP) if len(pl) and len(pr) else 0
        snr=delta/np.sqrt(max(1.,lv+rv+2*row['local_background']+1))
        reason=[]
        if row['105k_relative_change']<depth:reason.append('medium_relative_change')
        if row['25k_relative_change']<.7*depth or signs[0]!=direction:reason.append('fine_direction_or_magnitude')
        if sustain<3:reason.append('unstable_step_persistence')
        if spread>opts.max_platform_spread:reason.append('variable_platform_levels')
        if gain<.20:reason.append('small_two_level_SSE_gain')
        if advantage<opts.linear_gain:reason.append('continuous_slope_better_than_two_levels')
        if not STEP<=width<=opts.max_transition_bp:reason.append('transition_width')
        if delta<max(opts.min_absolute_change,.8*np.sqrt(row['local_background']+1)):reason.append('absolute_change_vs_background')
        if snr<row['effective_snr_threshold']:reason.append('low_absolute_density_SNR')
        support=2+int(signs[2]==direction and row['225k_relative_change']>=.05)
        row.update(direction='rise' if direction>0 else 'drop',absolute_change=delta,density_snr=snr,
            transition_width_bp=width,cp_flank_bp=flank,cp_two_level_SSE_gain=gain,step_vs_linear_gain=advantage,
            platform_spread=spread,sustained_block_n=sustain,persistence_blocks=4,scale_support_n=support,
            shape_pass=int(not reason),shape_reason=','.join(reason),strength=snr)
        tested.append(row)
    if not tested:return {**base,'shape_pass':0,'shape_reason':'insufficient_sustained_flanks','strength':0,'scale_support_n':0}
    return max(tested,key=lambda r:(r['shape_pass'],r['density_snr'],r['step_vs_linear_gain'],-r['cp_flank_bp']))


def evaluate(source,i,sig,parent,es,opts,acrs,genes,cache,states):
    row=sharp_shape(i,sig,parent['start'],parent['end'],states,opts) if source=='sharp' else level_shape(source,i,sig,parent['start'],parent['end'],states,opts)
    row.update(F.graph(es,parent['start'],parent['end'],row['boundary_pos'],opts))
    key=(parent['start'],parent['end'],i)
    if key not in cache:
        cache[key]=[A.short_evidence(es,parent['start'],parent['end'],row['boundary_pos'],cap,acrs,genes) for cap in [250000,500000]]
    short=cache[key]
    for d in short:row.update({f'short_{d["cap_bp"]}_{k}':v for k,v in d.items()})
    hard=[];p=row['boundary_pos']
    if min(p-parent['start'],parent['end']-p)<opts.min_daughter_bp:hard.append('minimum_daughter_500kb')
    if min(row['whole_left_internal'],row['whole_right_internal'])<3:hard.append('weak_whole_daughter_internal')
    if min(row['whole_left_ACR_n'],row['whole_right_ACR_n'])<2:hard.append('daughter_ACR_diversity')
    supported=[];density_only=[];all_votes=[]
    for d in short:
        enough=min(d['left_internal'],d['right_internal'])>=2 and min(d['left_ACR_n'],d['right_ACR_n'])>=2
        change=A.profile_change(d,{**row,'source':'step' if source=='changepoint' else source})
        raw_vote=row['local_cross_fraction']<=opts.max_touch
        local_vote=d['cross_fraction']<=opts.max_short_fraction
        depletion_vote=np.isfinite(d['relative_depletion']) and d['relative_depletion']>=0
        profile_vote=np.isfinite(change) and change>=opts.profile_min
        votes={'raw_touching':int(raw_vote),'short_fraction':int(local_vote),'short_depletion':int(depletion_vote),'short_profile':int(profile_vote)}
        all_votes.append({'cap':d['cap_bp'],**votes,'local_support':int(enough)})
        if enough and any(votes.values()):supported.append((d,change,votes))
        dominant=(enough and min(d['left_internal'],d['right_internal'])>=3 and row.get('high_density_context',0) and
                  row.get('density_snr',0)>=opts.dominant_snr and row.get('scale_support_n',0)>=2 and row['global_cut_loss']<=opts.dominant_global_loss)
        if dominant:density_only.append((d,change,votes))
    if not any(q['local_support'] for q in all_votes):hard.append('necessary_short_local_module_support')
    route='rejected';chosen=None
    if row['shape_pass'] and not hard:
        if supported:chosen=max(supported,key=lambda v:(sum(v[2].values()),-v[0]['cap_bp']));route='density_plus_local_connection'
        elif density_only:chosen=density_only[0];route='strong_density_with_two_sided_modules'
    if chosen:
        d,pc,votes=chosen;row.update(evidence_cap_bp=d['cap_bp'],evidence_profile_change=pc,
            connection_evidence_categories=','.join(k for k,v in votes.items() if v),
            evidence_left_internal=d['left_internal'],evidence_right_internal=d['right_internal'])
    reason=','.join(x for x in [row['shape_reason'],','.join(hard),'no_combined_connection_evidence' if not supported and not density_only else ''] if x)
    row.update(chrom=parent['chrom'],baseline_id=parent['baseline_id'],eligible=int(route!='rejected'),accepted=0,
               graph_pass=int(not hard and bool(supported or density_only)),structural_reason=','.join(hard),
               acceptance_route=route,reject_reason=reason,raw_touching_is_a_veto=False,
               short250_votes=json.dumps(all_votes[0],sort_keys=True),short500_votes=json.dumps(all_votes[1],sort_keys=True))
    return row


def call_parent(parent,sig,es,opts,method,acrs,genes):
    todo=[(parent['start'],parent['end'],0)];rows=[];cuts=[];leaves=[];cache={}
    while todo:
        start,end,depth=todo.pop();child={**parent,'start':start,'end':end};states=density_states(sig,start,end)
        proposals=sorted(set(F.proposals(sig,start,end,opts)+N.cp_proposals(sig,start,end,opts)),key=lambda r:(r[1],r[0]))
        cc=[]
        for source,i in proposals:
            r=evaluate(source,i,sig,child,es,opts,acrs,genes,cache,states);r.update(method=method,recursive_depth=depth);cc.append(r)
        valid=[r for r in cc if r['eligible']]
        if valid and depth<16:
            best=max(valid,key=lambda r:(r['high_density_context'],r['density_snr'],r.get('105k_relative_change',0),-r['local_cross_fraction'],-r['boundary_pos'],r['source']))
            best.update(accepted=1,reject_reason='accepted_'+best['acceptance_route']);cuts.append(dict(best));p=best['boundary_pos']
            todo.extend([(p,end,depth+1),(start,p,depth+1)])
            for r in valid:
                if r is not best:r['reject_reason']='eligible_rechecked_in_child'
        else:
            for r in valid:r['reject_reason']='recursion_safety_limit'
            leaves.append(child)
        rows.extend(cc)
    return leaves,cuts,rows


def main(args):
    root=args.outdir;hist=args.history_root
    settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items() if k not in ['stage','preview_approved']}
    root.mkdir(parents=True,exist_ok=True);signature=root/'01_test_settings.json'
    if signature.exists():assert json.loads(signature.read_text())==settings
    else:signature.write_text(json.dumps(settings,indent=2)+'\n')
    phase=root/('f01_preview' if args.stage=='preview' else 'f02_chromosomes')
    if phase.exists():raise ValueError('Stage exists; never overwrite')
    if args.stage=='chromosomes':assert args.preview_approved and (root/'f01_preview/f04_figures/f01_region_tracks/03_R_PLOTTING.DONE').exists()
    phase.mkdir()
    for d in ['f01_results','f02_plot_data','f03_evaluation','f04_figures','f05_scripts','f06_checks','f07_density_states']:(phase/d).mkdir()
    baseroot=hist/'f15_Aedge_B145_coACR_support/f01_baseline';prior=hist/'f20_V6_V7_local_changepoint/f02_chromosomes'
    cfg=json.loads((prior/'01_parameters.json').read_text());meta=json.loads((baseroot/'01_parameters.json').read_text())
    base=P.read_bed(baseroot/'02_CRU.bed',CHROMS);reference=P.read_bed(prior/'f01_results/f04_V7/02_CRU.bed',CHROMS)
    assert sum(map(len,base.values()))==365 and sum(map(len,reference.values()))==462
    original=C.rows(baseroot/'03_internal_boundaries.tsv');oldbounds=C.rows(prior/'f01_results/f04_V7/03_internal_boundaries.tsv')
    oldc=C.rows(prior/'f01_results/f04_V7/04_candidate_evaluations.tsv')
    paths=[Path(__file__),Path(N.__file__),Path(A.__file__),Path(F.__file__),Path(P.__file__),Path(C.__file__),
           baseroot/'02_CRU.bed',baseroot/'03_internal_boundaries.tsv',baseroot/'01_parameters.json',
           prior/'01_parameters.json',prior/'f01_results/f04_V7/02_CRU.bed',prior/'f01_results/f04_V7/03_internal_boundaries.tsv',
           prior/'f01_results/f04_V7/04_candidate_evaluations.tsv',Path(meta['aglink_bedpe']),Path(meta['A_gap_file']),Path(cfg['gene_bed']),Path(cfg['acr_bed'])]
    paths += [hist/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz' for c in CHROMS]
    hashes=[{'path':str(p),'md5_before':C.md5(p)} for p in paths]
    if args.stage=='chromosomes':
        oldhash={r['path']:r['md5_before'] for r in C.rows(root/'f01_preview/f06_checks/01_input_MD5.tsv')}
        assert all(oldhash[r['path']]==r['md5_before'] for r in hashes)
    acrs,genes=A.annotation_positions(cfg['gene_bed'],cfg['acr_bed']);edges=defaultdict(list)
    for r in C.rows(meta['aglink_bedpe']):
        if r['chrom'] not in CHROMS:continue
        for k in ['p1','p2','left','right']:r[k]=int(r[k])
        assert float(r['weight'])==1;edges[r['chrom']].append(r)
    results={m:{'crus':[],'boundaries':[],'candidates':[]} for m in METHODS};sigs={};globalstates={};selected=[]
    for m,bed,bb in [('Baseline',base,original),('V7_reference',reference,oldbounds)]:
        results[m]['crus']=[r for c in CHROMS for r in bed[c]];results[m]['boundaries']=[{**r,'method':m,'boundary_pos':int(r['boundary_pos'])} for r in bb]
    results['V7_reference']['candidates']=[{**r,'method':'V7_reference'} for r in oldc]
    for c in CHROMS:
        with np.load(hist/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz') as z:sig={k:z[k] for k in z.files}
        sigs[c]=sig;states=density_states(sig,0,cfg['sizes'][c]);globalstates[c]=states
        np.savez_compressed(phase/f'f07_density_states/{c}.density_states.npz',**states)
        parents=[{**r,'baseline_id':r['cru_id']} for r in base[c] if args.stage=='chromosomes' or any(c==cc and r['start']<b and r['end']>a for cc,a,b,f in REGIONS)]
        selected.extend(parents)
        for parent in parents:
            es=[e for e in edges[c] if e['left']<parent['end'] and e['right']>=parent['start']]
            for method in ['V8','V9']:
                rr,bb,cs=call_parent(parent,sig,es,tier(args,method),method,acrs[c],genes[c])
                results[method]['crus'].extend(rr);results[method]['boundaries'].extend(bb);results[method]['candidates'].extend(cs)
        for method in ['V8','V9']:
            results[method]['boundaries'].extend({**r,'source':'inherited_B145','method':method,'boundary_pos':int(r['boundary_pos'])} for r in original if r['chrom']==c and any(p['start']<=int(r['boundary_pos'])<=p['end'] for p in parents))
        print(args.stage,c,[(m,sum(r['chrom']==c for r in res['crus'])) for m,res in results.items()],flush=True)
    for method,folder in METHODS.items():
        d=phase/'f01_results'/folder;d.mkdir();res=results[method]
        if method in ['V8','V9']:
            for c in CHROMS:
                for i,r in enumerate(sorted((r for r in res['crus'] if r['chrom']==c),key=lambda r:r['start']),1):r['cru_id']=f'{c}_{method}_CRU_{i}'
            for parent in selected:
                rr=sorted((r for r in res['crus'] if r['baseline_id']==parent['baseline_id']),key=lambda r:r['start'])
                assert rr[0]['start']==parent['start'] and rr[-1]['end']==parent['end'] and all(l['end']==r['start'] for l,r in zip(rr,rr[1:]))
            if args.stage=='chromosomes':
                prev=P.read_bed(root/'f01_preview/f01_results'/folder/'02_CRU.bed',CHROMS)
                assert {(r['chrom'],r['start'],r['end']) for c in CHROMS for r in prev[c]}<={(r['chrom'],r['start'],r['end']) for r in res['crus']}
        C.savebed(d/'02_CRU.bed',res['crus']);write(d/'03_internal_boundaries.tsv',res['boundaries'])
        write(d/'04_candidate_evaluations.tsv',res['candidates'],None if res['candidates'] else ['chrom','boundary_pos'])
        (d/'01_parameters.json').write_text(json.dumps({'method':method,'stage':args.stage,'rules':vars(tier(args,method)) if method in ['V8','V9'] else {'frozen_reference':True}},indent=2)+'\n')
    cuts=[{'method':m,**r} for m in ['V8','V9'] for r in results[m]['boundaries'] if r['source'] in NEW_SOURCES]
    unique={(r['chrom'],int(r['boundary_pos'])):r for r in cuts if not any(r['chrom']==c and a<=int(r['boundary_pos'])<b for c,a,b,f in REGIONS)}
    keys=sorted(unique);rng=np.random.default_rng(args.random_seed)
    take=sorted(rng.choice(len(keys),min(args.random_sample_n,len(keys)),replace=False)) if args.stage=='chromosomes' and keys else []
    sample=[{**unique[keys[i]],'manual_inspection_status':'pending_visual_review'} for i in take]
    write(phase/'f03_evaluation/11_random_new_cut_sample.tsv',sample,None if sample else ['chrom','boundary_pos','method'])
    regions=[{'region':f'{c}_{a}_{b}','chrom':c,'start':a,'end':b,'focus':f,'kind':'representative'} for c,a,b,f in REGIONS]
    for j,r in enumerate(sample,1):
        c,p=r['chrom'],int(r['boundary_pos']);a=max(0,p-1000000);b=min(cfg['sizes'][c],p+1000000)
        regions.append({'region':f'random{j:02d}_{c}_{a}_{b}','chrom':c,'start':a,'end':b,'focus':p,'kind':'random_new_cut'})
    blocks=P.read_bed(cfg['block_bed'],CHROMS);m3path=Path(hist.parent/'f01_M3_genomewide/M3_CRU.bed');m3=P.read_bed(m3path,CHROMS)
    hashes += [{'path':str(p),'md5_before':C.md5(p)} for p in [Path(cfg['block_bed']),m3path]]
    pd=phase/'f02_plot_data';arcs=[];profiles=[];intervals=[];plotblocks=[];candidate=[];pgaps=[];boundaries=[];gaps=C.rows(meta['A_gap_file'])
    for reg in regions:
        c,a,b,rid=(reg[k] for k in ['chrom','start','end','region']);sig=sigs[c];states=globalstates[c]
        arcs.extend({'region':rid,**{k:e[k] for k in ['chrom','left','right','acr_id','gene_id']}} for e in edges[c] if e['left']<b and e['right']>a)
        x=(np.arange(len(sig['raw']))+.5)*STEP
        for i in np.flatnonzero((x>=a-STEP)&(x<=b+STEP)):
            profiles.append({'region':rid,'chrom':c,'position':x[i],'raw_span':sig['raw'][i],'span_25k':sig['smooth'][0,i],
                'span_105k':sig['smooth'][1,i],'span_225k':sig['smooth'][2,i],'local_background':sig['background'][i],
                'density_state':int(states['state'][i]),'high_threshold105':states['threshold105'][i],'high_threshold225':states['threshold225'][i]})
        for method,res in results.items():
            intervals.extend({'region':rid,**r,'method':method} for r in res['crus'] if r['chrom']==c and r['start']<b and r['end']>a)
            candidate.extend({'region':rid,**r,'method':method} for r in res['candidates'] if r['chrom']==c and a<=int(r['boundary_pos'])<b)
            boundaries.extend({'region':rid,**r,'method':method} for r in res['boundaries'] if r['chrom']==c and a<=int(r['boundary_pos'])<b)
        intervals.extend({'region':rid,**r,'method':'M3_reference'} for r in m3[c] if r['start']<b and r['end']>a)
        plotblocks.extend({'region':rid,**r} for r in blocks[c] if r['start']<b and r['end']>a)
        pgaps.extend({'region':rid,**r} for r in gaps if r['chrom']==c and int(r['start'])<b and int(r['end'])>a)
    write(pd/'01_regions.tsv',regions)
    for name,rr in [('02_aglink_arcs.tsv.gz',arcs),('03_signals.tsv',profiles),('04_intervals.tsv',intervals),('05_blocks.tsv',plotblocks),('06_candidates.tsv',candidate),('07_frozen_gaps.tsv',pgaps),('08_boundaries.tsv',boundaries)]:write(pd/name,rr)
    for r in hashes:r['md5_after']=C.md5(r['path']);assert r['md5_after']==r['md5_before']
    write(phase/'f06_checks/01_input_MD5.tsv',hashes);write(phase/'f06_checks/02_selected_baseline_parents.tsv',selected)
    (phase/'01_parameters.json').write_text(json.dumps({**settings,'stage':args.stage,'baseline_dir':str(baseroot.parent),'history_root':str(hist),
        'sizes':cfg['sizes'],'block_bed':cfg['block_bed'],'m3_script':cfg['m3_script'],'aglink_bedpe':meta['aglink_bedpe'],
        'gene_bed':cfg['gene_bed'],'acr_bed':cfg['acr_bed'],'no_coACR':True,'count_target_used':False},indent=2)+'\n')
    for p in Path(__file__).parent.iterdir():
        if p.is_file():shutil.copy2(p,phase/'f05_scripts'/p.name)
    (phase/'03_CALLING.DONE').write_text('Fixed B145 gaps, covered bases and original cuts. Optional V5-V7 cuts reevaluated. No count target or block input.\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--history-root',type=Path,required=True);ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--stage',choices=['preview','chromosomes'],required=True);ap.add_argument('--preview-approved',action='store_true')
    for prefix,values in [('v8',[.25,.30,2.5,.10,.85]),('v9',[.20,.25,2.0,.00,.90])]:
        for name,value in zip(['high-depth','low-depth','low-snr','linear-gain','max-touch'],values):ap.add_argument('--'+prefix+'-'+name,type=float,default=value)
    ap.add_argument('--high-snr',type=float,default=2.0);ap.add_argument('--min-daughter-bp',type=int,default=500000)
    ap.add_argument('--random-sample-n',type=int,default=20);ap.add_argument('--random-seed',type=int,default=20261108)
    args=ap.parse_args()
    if args.min_daughter_bp<500000:ap.error('Do not relax new daughters below 500 kb')
    main(args)
