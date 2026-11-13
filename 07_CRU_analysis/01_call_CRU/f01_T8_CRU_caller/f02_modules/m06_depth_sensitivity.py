# m06_depth_sensitivity: B1-45 depth sensitivity
# Frozen rule module, copied unchanged from the published calling chain
# (source: f10_B145_sensitivity/01_test_sensitivity.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Frozen B145 plus two sensitivity tiers; metadata-normalized short-link rescue."""
import argparse
from collections import defaultdict, Counter
import importlib.util
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
import numpy as np

def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

from . import m05_density_transitions as F
C,P=F.C,F.P
CHROMS=F.CHROMS;STEP=F.STEP
METHODS={'Baseline':'f01_baseline','V3_reference':'f02_V3_reference','V4':'f03_V4',
         'V5':'f04_V5','V5_no_strong':'f05_V5_no_strong'}
REGIONS=[('1A',5000000,10000000,7500000),('1B',15000000,20000000,18600000),
 ('1B',20000000,25000000,23000000),('1A',535000000,540000000,538400000),
 ('1A',34000000,36000000,35090000),('1B',50000000,59000000,55000000),
 ('1B',164000000,167000000,165000000),('2A',778000000,784000000,781000000)]


def tier(args,name):
    aggressive=name!='V4';prefix='v5' if aggressive else 'v4'
    return SimpleNamespace(min_absolute_change=2.,sharp_recovery_bp=250000,
        persistence_bp=200000,transition_exclusion_bp=25000,max_transition_bp=100000,
        sharp_depth=getattr(args,prefix+'_sharp_depth'),sharp_max_width_bp=getattr(args,prefix+'_max_width'),
        medium_depth=getattr(args,prefix+'_medium_depth'),step_relative=getattr(args,prefix+'_step_relative'),
        post_range=getattr(args,prefix+'_post_range'),depletion_min=getattr(args,prefix+'_depletion'),
        max_cross_fraction=.5*getattr(args,prefix+'_cross_multiplier'),
        min_daughter_bp=150000,min_local_internal=2,min_daughter_internal=3,max_global_cut_loss=.40,
        min_raw_valley_bp=args.min_raw_valley_bp,short_min_internal=args.short_min_internal,
        enable_strong=name=='V5',strong_change=args.strong_change,strong_depletion=args.strong_depletion,
        strong_profile=args.strong_profile,strong_min_internal=args.strong_min_internal,
        strong_max_cross=args.strong_max_cross)


def morphology(source,i,sig,start,end,opts):
    old=F.morph(source,i,sig,start,end,opts);r=dict(old)
    reasons=old['shape_reason'].split(',') if old['shape_reason'] else []
    if source=='sharp' and '25k_relative_change' in old:
        reasons=[z for z in reasons if z not in ['sharp_multiscale_support','not_a_sharp_width']]
        medium=old['105k_relative_change']>=opts.medium_depth
        if not medium:reasons.append('insufficient_105k_sharp_support')
        if not opts.min_raw_valley_bp<=old['transition_width_bp']<=opts.sharp_max_width_bp:reasons.append('sharp_width_or_single_bin')
        a,b=F.bounds(start,end,len(sig['raw']));lo=max(a,i-opts.sharp_max_width_bp//STEP);hi=min(b,i+opts.sharp_max_width_bp//STEP+1)
        level=old['25k_floor']+.5*(min(old['25k_left_level'],old['25k_right_level'])-old['25k_floor'])
        low=sig['raw'][lo:hi]<=level
        # Only a run containing the candidate can protect against one-bin artifacts.
        local=i-lo;l=h=local
        if low[local]:
            while l>0 and low[l-1]:l-=1
            while h+1<len(low) and low[h+1]:h+=1
            duration=(h-l+1)*STEP
        else:duration=0
        if duration<opts.min_raw_valley_bp:reasons.append('insufficient_raw_valley_duration')
        r.update(raw_valley_duration_bp=duration,medium_support_required=opts.medium_depth,
                 broad_support_auxiliary=int(old['225k_relative_change']>=.05))
        r['scale_support_n']=int(old['25k_relative_change']>=opts.sharp_depth)+int(medium)+int(old['225k_relative_change']>=.05)
    elif source=='step' and '25k_relative_change' in old:
        # A transient low bin cannot be rescued: all original sustained-step gates remain.
        w=opts.persistence_bp//STEP;skip=opts.transition_exclusion_bp//STEP
        post=sig['smooth'][0,i+skip+1:i+skip+1+w]
        medians=np.array([np.median(z) for z in np.array_split(post,4)])
        spread=float(np.ptp(medians)/max(old['25k_left_level'],old['25k_right_level'],1))
        if spread>opts.post_range:reasons.append('unstable_post_transition_platform')
        r['post_platform_range_over_high']=spread
    r.update(shape_pass=int(not reasons),shape_reason=','.join(reasons))
    r['strength']=max(0,r.get('25k_relative_change',0))*r.get('scale_support_n',0)
    return r


def pair_opportunities(acrs,genes,uppers):
    total=[]
    for d in uppers:
        total.append(int(np.sum(np.searchsorted(genes,acrs+d,side='right')-
                                np.searchsorted(genes,acrs-d,side='left'))))
    return np.diff([0]+total).astype(np.int64)


def edge_hist(es,mask,uppers):
    pairs={(e['acr_id'],e['gene_id']):e['right']-e['left'] for e,m in zip(es,mask) if m}
    d=np.array(list(pairs.values()),dtype=np.int64)
    bins=np.searchsorted(uppers,d,side='left');bins=bins[bins<len(uppers)]
    return np.bincount(bins,minlength=len(uppers))


def short_evidence(es,start,end,p,cap,acr_pos,gene_pos):
    flank=cap;a,b=max(start,p-flank),min(end,p+flank)
    left=np.array([e['left'] for e in es]);right=np.array([e['right'] for e in es])
    distance=right-left;valid=(left>=a)&(right<b)&(distance<=cap)
    il=valid&(right<p);ir=valid&(left>=p);cc=valid&(left<p)&(right>=p)
    L=F.unique_support(es,il);R=F.unique_support(es,ir);cross=F.unique_support(es,cc)[0]
    AL=acr_pos[(acr_pos>=a)&(acr_pos<p)];AR=acr_pos[(acr_pos>=p)&(acr_pos<b)]
    GL=gene_pos[(gene_pos>=a)&(gene_pos<p)];GR=gene_pos[(gene_pos>=p)&(gene_pos<b)]
    upper=np.array([25000,50000,100000,250000,500000,1000000]);upper=upper[upper<=cap]
    PL=pair_opportunities(AL,GL,upper);PR=pair_opportunities(AR,GR,upper)
    PC=pair_opportunities(AL,GR,upper)+pair_opportunities(AR,GL,upper)
    OL=edge_hist(es,il,upper);OR=edge_hist(es,ir,upper)
    ratesL=np.divide(OL,PL,out=np.zeros(len(PL),dtype=float),where=PL>0)
    ratesR=np.divide(OR,PR,out=np.zeros(len(PR),dtype=float),where=PR>0)
    expected=float(np.sum(PC*np.minimum(ratesL,ratesR)))
    depletion=1-cross/expected if expected>0 else np.nan
    short=(distance<=cap)&(left>=start)&(right<end)
    sl=np.sort(left[short]);sr=np.sort(right[short])
    def level(lo,hi):
        xx=np.arange(max(start,lo)+STEP//2,min(end,hi),STEP)
        if len(xx)<5:return np.nan
        return float(np.median(np.searchsorted(sl,xx,side='left')-np.searchsorted(sr,xx,side='left')))
    cl=level(p-150000,p-50000);cr=level(p+50000,p+150000);c0=level(p-12500,p+12500)
    valley=1-c0/min(cl,cr) if min(cl,cr)>0 else np.nan
    step=(cr-cl)/max(cl,cr) if max(cl,cr)>0 else np.nan
    return {'cap_bp':cap,'context_start':a,'context_end':b,'left_internal':L[0],'right_internal':R[0],
        'left_ACR_n':L[1],'right_ACR_n':R[1],'left_gene_n':L[2],'right_gene_n':R[2],
        'crossing':cross,'cross_fraction':cross/(cross+min(L[0],R[0])) if cross+min(L[0],R[0]) else 1.,
        'expected_cross_from_internal_opportunities':expected,'relative_depletion':depletion,
        'profile_left':cl,'profile_right':cr,'profile_center':c0,'profile_valley_depletion':valley,
        'profile_signed_step':step,'opportunities_left':int(PL.sum()),'opportunities_right':int(PR.sum()),
        'opportunities_cross':int(PC.sum()),'status':'available' if expected>0 and min(L[0],R[0])>=2 else 'insufficient_short_evidence'}


def profile_change(short,row):
    if row['source']=='sharp':return short['profile_valley_depletion']
    signed=short['profile_signed_step']
    return signed if row.get('direction')=='rise' else -signed


def evaluate(source,i,sig,parent,es,opts,acr_pos,gene_pos,cache):
    row=morphology(source,i,sig,parent['start'],parent['end'],opts)
    g=F.graph(es,parent['start'],parent['end'],row['boundary_pos'],opts)
    row.update(g)
    key=(parent['start'],parent['end'],i)
    if key not in cache:
        cache[key]=[short_evidence(es,parent['start'],parent['end'],row['boundary_pos'],cap,acr_pos,gene_pos)
                    for cap in [500000,1000000]]
    short=cache[key]
    for evidence in short:row.update({f'short_{evidence["cap_bp"]}_{k}':v for k,v in evidence.items()})
    # Length, module support and global-cut-loss vetoes are never overridden.
    structural=[]
    p=row['boundary_pos'];start,end=parent['start'],parent['end']
    if min(p-start,end-p)<opts.min_daughter_bp:structural.append('near_fixed_or_new_parent_edge')
    if min(g['whole_left_internal'],g['whole_right_internal'])<opts.min_daughter_internal:structural.append('weak_daughter_internal')
    if min(g['whole_left_ACR_n'],g['whole_right_ACR_n'])<2:structural.append('daughter_ACR_diversity')
    if g['global_cut_loss']>opts.max_global_cut_loss:structural.append('global_cut_loss')
    contexts=[]
    for f in [250000,500000]:
        ctx={k[len(f'context_{f}_'):]:v for k,v in g.items() if k.startswith(f'context_{f}_')}
        if min(ctx['left_internal'],ctx['right_internal'])>=opts.min_local_internal and min(ctx['left_ACR_n'],ctx['right_ACR_n'])>=2:contexts.append(ctx)
    if not contexts:structural.append('necessary_local_module_support')
    direct=[r for r in contexts if r['local_cross_fraction']<=opts.max_cross_fraction+1e-12]
    rescue=[];strong=[]
    if row['shape_pass'] and not structural:
        for d in short:
            enough=min(d['left_internal'],d['right_internal'])>=opts.short_min_internal and min(d['left_ACR_n'],d['right_ACR_n'])>=2
            change=profile_change(d,row)
            if enough and np.isfinite(change) and change>=opts.depletion_min and d['relative_depletion']>=opts.depletion_min and d['cross_fraction']<=opts.max_cross_fraction:
                rescue.append(d)
            isstrong=(source=='step' and row.get('25k_relative_change',0)>=opts.strong_change and
                      row.get('105k_relative_change',0)>=.40)
            enoughstrong=min(d['left_internal'],d['right_internal'])>=opts.strong_min_internal and min(d['left_ACR_n'],d['right_ACR_n'])>=2
            if opts.enable_strong and isstrong and enoughstrong and np.isfinite(change) and change>=opts.strong_profile and d['relative_depletion']>=opts.strong_depletion and d['cross_fraction']<=.75 and min(c['local_cross_fraction'] for c in contexts)<=opts.strong_max_cross:
                strong.append(d)
    route='rejected';chosen=None
    if row['shape_pass'] and not structural:
        if direct:route='direct_relaxed_crossing';selected=direct[0];row.update(selected)
        elif rescue:route='short_depletion_rescue';chosen=rescue[0];row.update(min(contexts,key=lambda c:c['local_cross_fraction']))
        elif strong:route='strong_drop_rescue';chosen=strong[0];row.update(min(contexts,key=lambda c:c['local_cross_fraction']))
    if chosen:row['rescue_cap_bp']=chosen['cap_bp']
    reasons=[r for r in [row['shape_reason'],','.join(structural)] if r]
    if not direct and not rescue and not strong:reasons.append('touching_crossing_and_no_short_rescue')
    row.update(chrom=parent['chrom'],baseline_id=parent['baseline_id'],eligible=int(route!='rejected'),accepted=0,
        raw_graph_pass=g['graph_pass'],graph_pass=int(not structural and bool(direct or rescue or strong)),
        acceptance_route=route,overridden_veto='touching_crossing_fraction' if route.endswith('rescue') else '',
        raw_crossing_limit=opts.max_cross_fraction,reject_reason=','.join(reasons),
        unrelaxed_crossing_would_fail=int(row['local_cross_fraction']>.50),
        strong_rescue_available=int(bool(strong)),short_rescue_available=int(bool(rescue)))
    return row


def call_parent(parent,sig,es,opts,method,acr_pos,gene_pos,cache):
    todo=[(parent['start'],parent['end'],0)];rows=[];cuts=[];leaves=[]
    while todo:
        start,end,depth=todo.pop();child={**parent,'start':start,'end':end}
        candidates=[]
        for source,i in F.proposals(sig,start,end,opts):
            r=evaluate(source,i,sig,child,es,opts,acr_pos,gene_pos,cache)
            r.update(method=method,recursive_depth=depth);candidates.append(r)
        valid=[r for r in candidates if r['eligible']]
        if valid and depth<16:
            best=max(valid,key=lambda r:(r['strength']/(1+r['local_cross_fraction']),r['scale_support_n'],-r['boundary_pos'],r['source']))
            best['accepted']=1;best['reject_reason']='accepted_'+best['acceptance_route'];cuts.append(dict(best))
            for r in valid:
                if r is not best:r['reject_reason']='eligible_rechecked_in_child'
            p=best['boundary_pos'];todo.extend([(p,end,depth+1),(start,p,depth+1)])
        else:
            if valid:
                for r in valid:r['reject_reason']='recursion_safety_limit'
            leaves.append(child)
        rows.extend(candidates)
    return leaves,cuts,rows


def annotation_positions(gene_bed,acr_bed):
    genes=defaultdict(list);acrs=defaultdict(list);known={}
    for path,isgene in [(gene_bed,True),(acr_bed,False)]:
        with Path(path).open() as f:
            for line in f:
                s=line.split()
                if not s or s[0] not in CHROMS:continue
                c,a,b,k=s[:4];a,b=int(a),int(b)
                p=(b-1 if s[5]=='-' else a) if isgene else (a+b)//2
                ident=(isgene,k)
                if ident in known:assert known[ident]==(c,p);continue
                known[ident]=(c,p)
                (genes if isgene else acrs)[c].append(p)
    return {c:np.sort(v) for c,v in acrs.items()},{c:np.sort(v) for c,v in genes.items()}


def main(args):
    root=args.outdir;baseline=args.baseline_dir;reference=args.reference_dir
    settings={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items() if k not in ['stage','preview_approved']}
    root.mkdir(parents=True,exist_ok=True)
    param=root/'01_test_settings.json'
    if param.exists():assert json.loads(param.read_text())==settings
    else:param.write_text(json.dumps(settings,indent=2)+'\n')
    phase=root/('f01_preview' if args.stage=='preview' else 'f02_chromosomes')
    if phase.exists():raise ValueError('Stage exists; do not overwrite')
    if args.stage=='chromosomes':assert args.preview_approved and (root/'f01_preview/f04_figures/f01_region_tracks/03_R_PLOTTING.DONE').exists()
    phase.mkdir()
    for d in ['f01_results','f02_plot_data','f03_evaluation','f04_figures','f05_scripts','f06_checks']:(phase/d).mkdir()
    meta=json.loads((baseline/'f01_baseline/01_parameters.json').read_text());hist=Path(meta['history_root'])
    bed=baseline/'f01_baseline/02_CRU.bed';assert C.md5(bed)==(baseline/'f01_baseline/12_BASELINE.DONE').read_text().strip()
    base=P.read_bed(bed,CHROMS);assert sum(map(len,base.values()))==365
    sources=[bed,Path(meta['aglink_bedpe']),Path(meta['A_gap_file']),baseline/'f01_baseline/03_internal_boundaries.tsv',
             reference/'f01_results/f04_V3_both/02_CRU.bed',reference/'f01_results/f04_V3_both/03_internal_boundaries.tsv',
             reference/'f01_results/f04_V3_both/04_candidate_evaluations.tsv',Path(args.gene_bed),Path(args.acr_bed),Path(__file__),Path(F.__file__)]
    sources += [hist/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz' for c in CHROMS]
    hashes=[{'path':str(p),'md5_before':C.md5(p)} for p in sources]
    previous=root/'f01_preview/f06_checks/01_input_MD5.tsv'
    if args.stage=='chromosomes':
        expected={r['path']:r['md5_before'] for r in C.rows(previous)}
        assert all(expected[h['path']]==h['md5_before'] for h in hashes)
    acrs,genes=annotation_positions(args.gene_bed,args.acr_bed)
    edges=defaultdict(list)
    for r in C.rows(meta['aglink_bedpe']):
        if r['chrom'] not in CHROMS:continue
        for k in ['p1','p2','left','right']:r[k]=int(r[k])
        assert float(r['weight'])==1;edges[r['chrom']].append(r)
    original=C.rows(baseline/'f01_baseline/03_internal_boundaries.tsv')
    oldc=C.rows(reference/'f01_results/f04_V3_both/04_candidate_evaluations.tsv')
    audit=[]
    for r in oldc:
        C0=int(r['crossing_touching_context']);W=min(int(r['left_internal']),int(r['right_internal']))
        value=C0/(C0+W) if C0+W else 1
        assert np.isclose(value,float(r['local_cross_fraction']))
        audit.append({**r,'recomputed_touching_fraction':value,'passes_cross_0625':int(value<=.625),
                      'passes_cross_075':int(value<=.75),'cross_over_weak':C0/W if W else 'NA'})
    C.write(phase/'f06_checks/02_original_V3_crossing_audit.tsv',audit)
    results={m:{'crus':[],'boundaries':[],'candidates':[]} for m in METHODS}
    v3=P.read_bed(reference/'f01_results/f04_V3_both/02_CRU.bed',CHROMS)
    v3bounds=C.rows(reference/'f01_results/f04_V3_both/03_internal_boundaries.tsv')
    for m,by,bb in [('Baseline',base,original),('V3_reference',v3,v3bounds)]:
        results[m]['crus']=[r for c in CHROMS for r in by[c]]
        results[m]['boundaries']=[{**r,'boundary_pos':int(r['boundary_pos']),'method':m} for r in bb]
    results['V3_reference']['candidates']=[{**r,'method':'V3_reference'} for r in oldc]
    sigs={};selected=[]
    for c in CHROMS:
        with np.load(hist/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz') as z:sig={k:z[k] for k in z.files}
        sigs[c]=sig
        parents=[{**r,'baseline_id':r['cru_id']} for r in base[c] if args.stage=='chromosomes' or
                 any(c==cc and r['start']<b and r['end']>a for cc,a,b,t in REGIONS)]
        selected.extend(parents)
        for parent in parents:
            es=[e for e in edges[c] if e['left']<parent['end'] and e['right']>=parent['start']];cache={}
            for method in ['V4','V5','V5_no_strong']:
                rr,bb,cs=call_parent(parent,sig,es,tier(args,method),method,acrs[c],genes[c],cache)
                results[method]['crus'].extend(rr);results[method]['boundaries'].extend(bb);results[method]['candidates'].extend(cs)
        for method in ['V4','V5','V5_no_strong']:
            results[method]['boundaries'].extend({**r,'boundary_pos':int(r['boundary_pos']),'source':'inherited_B145','method':method}
                for r in original if r['chrom']==c and any(p['start']<=int(r['boundary_pos'])<=p['end'] for p in parents))
        print(args.stage,c,[(m,sum(r['chrom']==c for r in res['crus'])) for m,res in results.items()],flush=True)
    for method,res in results.items():
        folder=phase/'f01_results'/METHODS[method];folder.mkdir()
        if method in ['V4','V5','V5_no_strong']:
            for c in CHROMS:
                for i,r in enumerate(sorted((r for r in res['crus'] if r['chrom']==c),key=lambda r:r['start']),1):r['cru_id']=f'{c}_{method}_CRU_{i}'
            for parent in selected:
                rr=sorted((r for r in res['crus'] if r['baseline_id']==parent['baseline_id']),key=lambda r:r['start'])
                assert rr[0]['start']==parent['start'] and rr[-1]['end']==parent['end']
                assert all(a['end']==b['start'] for a,b in zip(rr,rr[1:]))
        C.savebed(folder/'02_CRU.bed',res['crus']);C.write(folder/'03_internal_boundaries.tsv',res['boundaries'])
        C.write(folder/'04_candidate_evaluations.tsv',res['candidates'],None if res['candidates'] else ['chrom','boundary_pos'])
        (folder/'01_parameters.json').write_text(json.dumps({'method':method,'stage':args.stage,
            'rules':vars(tier(args,method)) if method in ['V4','V5','V5_no_strong'] else {'frozen_reference':True}},indent=2)+'\n')
        if args.stage=='chromosomes' and method in ['V4','V5','V5_no_strong']:
            prev=P.read_bed(root/'f01_preview/f01_results'/METHODS[method]/'02_CRU.bed',CHROMS)
            coords={(r['chrom'],r['start'],r['end']) for r in res['crus']}
            assert {(r['chrom'],r['start'],r['end']) for c in CHROMS for r in prev[c]}<=coords
    C.write(phase/'f06_checks/03_selected_full_parents.tsv',selected)
    # Only after all algorithmic decisions, load fixed block/M3 tracks for R export.
    oldparams=json.loads((hist/'f08_span_only_control/01_parameters.json').read_text())
    blockpath=Path(meta['block_bed']);m3path=Path(oldparams['m3_dir'])/'M3_CRU.bed'
    blocks=P.read_bed(blockpath,CHROMS);m3=P.read_bed(m3path,CHROMS)
    hashes += [{'path':str(p),'md5_before':C.md5(p)} for p in [blockpath,m3path]]
    pd=phase/'f02_plot_data';arcs=[];profiles=[];intervals=[];plotblocks=[];candidates=[];pgaps=[];bounds=[]
    gaps=C.rows(meta['A_gap_file'])
    C.write(pd/'01_regions.tsv',[{'region':f'{c}_{a}_{b}','chrom':c,'start':a,'end':b,'focus':t} for c,a,b,t in REGIONS])
    for c,a,b,t in REGIONS:
        region=f'{c}_{a}_{b}';sig=sigs[c];x=(np.arange(len(sig['raw']))+.5)*STEP
        arcs.extend({'region':region,**{k:e[k] for k in ['chrom','left','right','acr_id','gene_id']}} for e in edges[c] if e['left']<b and e['right']>a)
        for i in np.flatnonzero((x>=a-STEP)&(x<=b+STEP)):
            profiles.append({'region':region,'chrom':c,'position':x[i],'raw_span':sig['raw'][i],
                'span_25k':sig['smooth'][0,i],'span_105k':sig['smooth'][1,i],'span_225k':sig['smooth'][2,i],'local_background':sig['background'][i]})
        for method,res in results.items():
            intervals.extend({'region':region,'method':method,**r} for r in res['crus'] if r['chrom']==c and r['start']<b and r['end']>a)
            candidates.extend({'region':region,**r} for r in res['candidates'] if r['chrom']==c and a<=int(r['boundary_pos'])<b)
            bounds.extend({'region':region,'method':method,**r} for r in res['boundaries'] if r['chrom']==c and a<=int(r['boundary_pos'])<b)
        intervals.extend({'region':region,'method':'M3_reference',**r} for r in m3[c] if r['start']<b and r['end']>a)
        plotblocks.extend({'region':region,**r} for r in blocks[c] if r['start']<b and r['end']>a)
        pgaps.extend({'region':region,**r} for r in gaps if r['chrom']==c and int(r['start'])<b and int(r['end'])>a)
    C.write(pd/'02_aglink_arcs.tsv.gz',arcs);C.write(pd/'03_signals.tsv',profiles);C.write(pd/'04_intervals.tsv',intervals)
    C.write(pd/'05_blocks.tsv',plotblocks);C.write(pd/'06_candidates.tsv',candidates);C.write(pd/'07_frozen_gaps.tsv',pgaps);C.write(pd/'08_boundaries.tsv',bounds)
    for h in hashes:
        h['md5_after']=C.md5(h['path']);assert h['md5_before']==h['md5_after']
    C.write(phase/'f06_checks/01_input_MD5.tsv',hashes)
    cfg={**settings,'stage':args.stage,'sizes':meta['sizes'],'block_bed':meta['block_bed'],
         'm3_script':oldparams['m3_script'],'aglink_bedpe':meta['aglink_bedpe'],'baseline_MD5':C.md5(bed),
         'no_coACR':True,'block_positions_used_for_calling':False}
    (phase/'01_parameters.json').write_text(json.dumps(cfg,indent=2)+'\n')
    for p in Path(__file__).parent.iterdir():
        if p.is_file():shutil.copy2(p,phase/'f05_scripts'/p.name)
    (phase/'03_CALLING.DONE').write_text('Frozen B145 geometry; no coACR; no block-guided calling.\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    for k in ['baseline-dir','reference-dir','outdir','gene-bed','acr-bed']:ap.add_argument('--'+k,type=Path,required=True)
    ap.add_argument('--stage',choices=['preview','chromosomes'],required=True);ap.add_argument('--preview-approved',action='store_true')
    for prefix,values in [('v4',[1.25,.40,.10,200000,.40,.30,.25]),('v5',[1.50,.35,.075,250000,.325,.40,.15])]:
        for name,value in zip(['cross-multiplier','sharp-depth','medium-depth','max-width','step-relative','post-range','depletion'],values):
            ap.add_argument('--'+prefix+'-'+name,type=int if name=='max-width' else float,default=value)
    ap.add_argument('--min-raw-valley-bp',type=int,default=10000);ap.add_argument('--short-min-internal',type=int,default=3)
    ap.add_argument('--strong-change',type=float,default=.60);ap.add_argument('--strong-depletion',type=float,default=.10)
    ap.add_argument('--strong-profile',type=float,default=.40);ap.add_argument('--strong-min-internal',type=int,default=5)
    ap.add_argument('--strong-max-cross',type=float,default=.90)
    main(ap.parse_args())
