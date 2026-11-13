# m10_local_subdivision: V11 local subdivision
# Frozen rule module, copied unchanged from the published calling chain
# (source: f14_local_density_subdivision/01_test_local_density.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Local agLINK profiles and simultaneous compatible cuts; bounded exploratory run."""
import argparse
from collections import Counter,defaultdict
import importlib.util
import json
from pathlib import Path
import shutil
from types import SimpleNamespace
import numpy as np
from scipy.ndimage import uniform_filter1d


def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m

from . import m09_two_stage_density as T
A,F,C,P,N=T.A,T.F,T.C,T.P,T.N
from . import m02_adaptive_background as D
write=P.write_tsv
CHROMS=['1A','1B','2A'];STEP=5000
METHODS={'Baseline':'f01_baseline','V7':'f02_V7_reference','V9':'f03_V9_reference',
         'V10_500kb':'f04_V10_500kb','V11_500kb':'f05_V11_500kb','V11_300kb':'f06_V11_300kb'}
REGIONS=[('1A',34000000,36000000),('1B',15000000,20000000),('1B',50000000,59000000)]


def weighted_selection(rows,start,end,minimum):
    """Maximum summed log(1+shape effect) with interval compatibility, not count."""
    best={}
    for r in rows:
        p=int(r['boundary_pos'])
        if min(p-start,end-p)<minimum:continue
        if p not in best or float(r['density_snr'])>float(best[p]['density_snr']):best[p]=r
    positions=sorted(best)
    if not positions:return []
    weights=np.array([np.log1p(max(0.,float(best[p]['density_snr']))) for p in positions])
    previous=np.searchsorted(positions,np.array(positions)-minimum,side='right')-1
    dp=np.zeros(len(positions)+1);take=np.zeros(len(positions),bool)
    for i in range(len(positions)):
        candidate=weights[i]+dp[previous[i]+1]
        if candidate>dp[i]+1e-12:dp[i+1]=candidate;take[i]=True
        else:dp[i+1]=dp[i]
    selected=[];i=len(positions)-1
    while i>=0:
        if take[i]:selected.append(best[positions[i]]);i=previous[i]
        else:i-=1
    return list(reversed(selected))


def span_profile(es,size,origin,weights):
    n=(size-origin+STEP-1)//STEP;diff=np.zeros(n+1);partial=np.zeros(n)
    if es:
        left=np.array([e['left']-origin for e in es]);end=np.array([e['right']+1-origin for e in es])
        aa=left//STEP;bb=(end-1)//STEP;w=np.asarray(weights);same=aa==bb
        np.add.at(partial,aa[same],w[same]*(end[same]-left[same]))
        k=~same
        np.add.at(partial,aa[k],w[k]*((aa[k]+1)*STEP-left[k]));np.add.at(partial,bb[k],w[k]*(end[k]-bb[k]*STEP))
        np.add.at(diff,aa[k]+1,w[k]);np.add.at(diff,bb[k],-w[k])
    lengths=np.minimum(STEP,size-origin-np.arange(n)*STEP)
    return np.maximum(0.,np.cumsum(diff[:-1])+partial/lengths)


def shifted_edges(es,offset):
    return [{**e,**{k:e[k]-offset for k in ['p1','p2','left','right']}} for e in es]


def profile_context(sig,es,chromsize,start,end,rule):
    lo=start//STEP;hi=min(len(sig['raw']),(end+STEP-1)//STEP);origin=lo*STEP;size=min(chromsize,hi*STEP)
    inside=[e for e in es if start<=e['left'] and e['right']<end]
    distances=np.array([e['right']-e['left'] for e in inside])
    profiles={'full':{'raw':sig['raw'][lo:hi], 'smooth':sig['smooth'][:,lo:hi], 'background':sig['background'][lo:hi]}}
    for name,cap in [('short500',500000),('short1M',1000000),('weighted',None)]:
        use=inside if cap is None else [e for e in inside if e['right']-e['left']<=cap]
        w=1/np.sqrt(np.maximum(1.,distances/500000.)) if cap is None else np.ones(len(use))
        raw=span_profile(use,size,origin,w);smooth=np.array([uniform_filter1d(raw,k,mode='nearest') for k in [5,21,45]])
        profiles[name]={'raw':raw,'smooth':smooth,'background':D.local_background(smooth[0],rule)}
    return {'origin':origin,'start':start-origin,'end':end-origin,'profiles':profiles,
            'edges':shifted_edges([e for e in es if e['left']<end and e['right']>=start],origin)}


def options(settings,minimum):
    opt=T.tier(SimpleNamespace(**settings),'V9');opt.min_daughter_bp=minimum;return opt


def proposals(context,opt):
    out=set()
    for name,sig in context['profiles'].items():
        local=set(F.proposals(sig,context['start'],context['end'],opt)+N.cp_proposals(sig,context['start'],context['end'],opt))
        for source,i in local:out.add((i*STEP+STEP//2+context['origin'],source))
    return out


def evaluate(position,source,context,opt,acrs,genes,method):
    off=context['origin'];p=position-off;i=(p-STEP//2)//STEP;start,end=context['start'],context['end']
    row={'boundary_pos':position,'source':source,'parent_start':start+off,'parent_end':end+off,'method':method,
         'accepted':0,'eligible':0,'density_snr':0.,'local_profile':'none','reject_reason':''}
    if min(p-start,end-p)<opt.min_daughter_bp:
        return {**row,'reject_reason':'minimum_daughter_length'}
    shapes=[]
    for name,sig in context['profiles'].items():
        if 'states' not in sig:sig['states']=T.density_states(sig,start,end)
        states=sig['states']
        shape=T.sharp_shape(i,sig,start,end,states,opt) if source=='sharp' else T.level_shape(source,i,sig,start,end,states,opt)
        shapes.append({**shape,'local_profile':name})
    eligible=[r for r in shapes if r['shape_pass']]
    if not eligible:
        best=max(shapes,key=lambda r:r.get('density_snr',0));return {**row,**best,'boundary_pos':position,
            'parent_start':start+off,'parent_end':end+off,'reject_reason':'shape:'+best['shape_reason']}
    local=[r for r in eligible if r['local_profile']!='full']
    chosen=max(local or eligible,key=lambda r:(r['density_snr'],r.get('105k_relative_change',0)))
    graph=F.graph(context['edges'],start,end,p,opt)
    hard=[]
    if min(graph['whole_left_internal'],graph['whole_right_internal'])<3:hard.append('whole_internal_support')
    if min(graph['whole_left_ACR_n'],graph['whole_right_ACR_n'])<2:hard.append('whole_ACR_diversity')
    short=[A.short_evidence(context['edges'],start,end,p,cap,acrs-off,genes-off) for cap in [250000,500000]]
    connected=[];modules=[]
    for d in short:
        support=min(d['left_internal'],d['right_internal'])>=2 and min(d['left_ACR_n'],d['right_ACR_n'])>=2
        change=A.profile_change(d,{'source':'step' if source=='changepoint' else source,'direction':chosen.get('direction')})
        if support:modules.append(d)
        if support and (d['cross_fraction']<=.85 or d['relative_depletion']>=0 or (np.isfinite(change) and change>=.05)):connected.append(d)
    if not modules:hard.append('short_local_module_support')
    # Local shape can pass even when full-span shape fails: that is the structural test.
    route='local_density_with_short_separation' if local and connected else ''
    if not route and not local and connected:route='full_reference_with_short_separation'
    if not route and method.startswith('V11') and local and any(min(d['left_internal'],d['right_internal'])>=3 for d in modules):
        route='local_density_with_two_sided_modules'
    if not route:hard.append('connection_evidence')
    out={**row,**chosen,**graph,'boundary_pos':position,'parent_start':start+off,'parent_end':end+off,
         'shape_pass':1,'eligible':int(not hard),'accepted':0,'acceptance_route':route,
         'reject_reason':','.join(hard),'local_profiles_passing':','.join(r['local_profile'] for r in eligible),
         'full_shape_pass':int(any(r['local_profile']=='full' for r in eligible))}
    for d in short:
        out.update({f'short_{d["cap_bp"]}_{k}':v for k,v in d.items()})
    for key,value in list(out.items()):
        if (key.endswith('_peak_pos') or key.endswith('_context_start') or key.endswith('_context_end')) and np.isfinite(value):out[key]=value+off
    out['profile_origin']=off
    return out


def solve(parent,sig,edges,chromsize,rule,opt,acrs,genes,method,max_rounds):
    start,end=parent['start'],parent['end'];current=set();pool=set();records=[];seen=set();iterations=[];final=[]
    if end-start<2*opt.min_daughter_bp:return [dict(parent)],[],[],[]
    for iteration in range(max_rounds):
        points=[start]+sorted(current)+[end];contexts={};evaluations={}
        def context(a,b):
            if (a,b) not in contexts:contexts[(a,b)]=profile_context(sig,edges,chromsize,a,b,rule)
            return contexts[(a,b)]
        for a,b in zip(points,points[1:]):
            if b-a>=2*opt.min_daughter_bp:pool.update(proposals(context(a,b),opt))
        for p,source in sorted(pool):
            if p in current:
                j=points.index(p);a,b=points[j-1],points[j+1]
            else:
                j=np.searchsorted(points,p)-1;a,b=points[j],points[j+1]
            r=evaluate(p,source,context(a,b),opt,acrs,genes,method)
            r.update(chrom=parent['chrom'],baseline_id=parent['cru_id'],iteration=iteration)
            evaluations[(p,source)]=r;records.append(r)
        valid=[r for r in evaluations.values() if r['eligible']]
        chosen=weighted_selection(valid,start,end,opt.min_daughter_bp)
        # Refresh all selected neighborhoods simultaneously; retry alternatives, not greedy cuts.
        invalid=set();stable=[]
        for retry in range(12):
            candidates=[r for r in valid if (int(r['boundary_pos']),r['source']) not in invalid]
            chosen=weighted_selection(candidates,start,end,opt.min_daughter_bp)
            pp=[start]+[int(r['boundary_pos']) for r in chosen]+[end];stable=[];failed=[]
            for j,r in enumerate(chosen,1):
                q=evaluate(int(r['boundary_pos']),r['source'],context(pp[j-1],pp[j+1]),opt,acrs,genes,method)
                q.update(chrom=parent['chrom'],baseline_id=parent['cru_id'],iteration=iteration)
                if q['eligible']:stable.append(q)
                else:failed.append((int(r['boundary_pos']),r['source']))
            if not failed:break
            invalid.update(failed)
        # Rare retry exhaustion keeps only individually revalidated cuts, then rechecks again.
        while True:
            pp=[start]+[int(r['boundary_pos']) for r in stable]+[end];updated=[]
            for j,r in enumerate(stable,1):
                q=evaluate(int(r['boundary_pos']),r['source'],context(pp[j-1],pp[j+1]),opt,acrs,genes,method)
                q.update(chrom=parent['chrom'],baseline_id=parent['cru_id'],iteration=iteration)
                if q['eligible']:updated.append(q)
            if len(updated)==len(stable):stable=updated;break
            stable=updated
        new={int(r['boundary_pos']) for r in stable}
        iterations.append({'chrom':parent['chrom'],'baseline_id':parent['cru_id'],'method':method,'iteration':iteration,
                           'candidate_pairs':len(pool),'eligible_evaluations':len(valid),'selected_n':len(new),'joint_rejected_n':len(invalid)})
        final=stable
        if new==current:break
        key=tuple(sorted(new))
        if key in seen:break
        seen.add(key);current=new
    for r in final:r.update(accepted=1,reject_reason='accepted_'+r['acceptance_route'])
    pp=[start]+[int(r['boundary_pos']) for r in final]+[end]
    leaves=[{**parent,'start':a,'end':b} for a,b in zip(pp,pp[1:])]
    return leaves,final,records,iterations


def diagnose(hist,out):
    source=hist/'f22_two_stage_density_V8_V9/f02_chromosomes/f01_results/f04_V9'
    rows=C.rows(source/'04_candidate_evaluations.tsv');stages=[];reasons=[]
    for kind in ['ALL','sharp','step','changepoint']:
        rr=rows if kind=='ALL' else [r for r in rows if r['source']==kind]
        stages.append({'source':kind,'evaluations':len(rr),'unique_positions':len({(r['chrom'],r['boundary_pos']) for r in rr}),
            'shape_pass':sum(int(r['shape_pass']) for r in rr),'graph_pass':sum(int(r['graph_pass']) for r in rr),
            'eligible':sum(int(r['eligible']) for r in rr),'accepted':sum(int(r['accepted']) for r in rr)})
        counter=Counter(z for r in rr if not int(r['accepted']) for z in r['reject_reason'].split(',') if z)
        reasons.extend({'source':kind,'reason':k,'rejected_evaluations':v,'nonexclusive':1} for k,v in counter.most_common())
    write(out/'01_V9_stage_counts.tsv',stages);write(out/'02_V9_rejection_counts.tsv',reasons)
    capacities=[]
    for name,path in [('B145',hist/'f15_Aedge_B145_coACR_support/f01_baseline/02_CRU.bed'),('V9',source/'02_CRU.bed')]:
        bed=P.read_bed(path,CHROMS)
        for minimum in [500000,300000]:
            for c in CHROMS+['Tested_chromosomes']:
                rr=[r for cc in CHROMS for r in bed[cc] if c=='Tested_chromosomes' or cc==c]
                upper=sum(max(1,(r['end']-r['start'])//minimum) for r in rr)
                capacities.append({'reference':name,'chrom':c,'minimum_bp':minimum,'current_n':len(rr),'length_only_upper_n':upper,'extra_cuts_upper_n':upper-len(rr)})
    write(out/'03_length_capacity.tsv',capacities)
    order=[]
    for bid in sorted({r['baseline_id'] for r in rows}):
        rr=[r for r in rows if r['baseline_id']==bid and int(r['recursive_depth'])==0 and int(r['eligible'])]
        if not rr:continue
        a,b=int(rr[0]['parent_start']),int(rr[0]['parent_end']);dp=weighted_selection(rr,a,b,500000)
        greedy=[]
        for r in sorted(rr,key=lambda r:float(r['density_snr']),reverse=True):
            p=int(r['boundary_pos'])
            if all(abs(p-int(q['boundary_pos']))>=500000 for q in greedy):greedy.append(r)
        order.append({'baseline_id':bid,'root_eligible_unique_n':len({r['boundary_pos'] for r in rr}),
            'fixed_pool_greedy_n':len(greedy),'fixed_pool_DP_n':len(dp),
            'counterfactual_only':'Same root candidates; no child re-evaluation, not actual CRU counts'})
    write(out/'04_V9_fixed_pool_order_diagnostic.tsv',order)
    (out/'05_structural_diagnosis.txt').write_text('V9 density states only alter morphology thresholds and priority. There is no independent high-density local-span generator; all three morphology families use the same full-span input. Repeated eligible_rechecked_in_child records are not unique suppressed boundaries. Fixed-pool order contrasts are descriptive, not proof that greedy alone caused every rejection.\n')


def main(args):
    root=args.outdir;hist=args.history_root
    if root.exists():raise ValueError('Output exists; never overwrite')
    root.mkdir(parents=True)
    for d in ['f01_diagnosis','f02_inputs','f03_results','f04_statistics','f05_figures','f06_scripts']:(root/d).mkdir()
    diagnose(hist,root/'f01_diagnosis')
    cfg=json.loads((hist/'f22_two_stage_density_V8_V9/f02_chromosomes/01_parameters.json').read_text())
    settings=json.loads((hist/'f22_two_stage_density_V8_V9/01_test_settings.json').read_text())
    rule=SimpleNamespace(**json.loads((hist/'f15_Aedge_B145_coACR_support/f01_baseline/01_parameters.json').read_text())['primary_rule'])
    base=P.read_bed(hist/'f15_Aedge_B145_coACR_support/f01_baseline/02_CRU.bed',CHROMS)
    original=C.rows(hist/'f15_Aedge_B145_coACR_support/f01_baseline/03_internal_boundaries.tsv')
    paths=[Path(__file__),Path(T.__file__),Path(A.__file__),Path(N.__file__),Path(F.__file__),Path(D.__file__),Path(P.__file__),
        hist/'f15_Aedge_B145_coACR_support/f01_baseline/02_CRU.bed',hist/'f15_Aedge_B145_coACR_support/f01_baseline/03_internal_boundaries.tsv',
        Path(cfg['aglink_bedpe']),Path(cfg['gene_bed']),Path(cfg['acr_bed'])]
    hashes=[{'path':str(p),'md5_before':C.md5(p)} for p in paths]
    edges=D.load_edges(cfg['aglink_bedpe'],CHROMS);acrs,genes=A.annotation_positions(cfg['gene_bed'],cfg['acr_bed'])
    results={m:{'crus':[],'boundaries':[],'candidates':[],'iterations':[]} for m in METHODS}
    results['Baseline']['crus']=[r for c in CHROMS for r in base[c]];results['Baseline']['boundaries']=original
    for name,path in [('V7',hist/'f20_V6_V7_local_changepoint/f02_chromosomes/f01_results/f04_V7'),('V9',hist/'f22_two_stage_density_V8_V9/f02_chromosomes/f01_results/f04_V9')]:
        bed=P.read_bed(path/'02_CRU.bed',CHROMS);results[name]['crus']=[r for c in CHROMS for r in bed[c]];results[name]['boundaries']=C.rows(path/'03_internal_boundaries.tsv')
        hashes.append({'path':str(path/'02_CRU.bed'),'md5_before':C.md5(path/'02_CRU.bed')})
    signals={};influence=[]
    for c in CHROMS:
        signalpath=hist/f'f13_B_route_hard_gap_first/f01_inputs/{c}.B_signals.npz'
        with np.load(signalpath) as z:sig={k:z[k] for k in z.files}
        hashes.append({'path':str(signalpath),'md5_before':C.md5(signalpath)})
        context=profile_context(sig,edges[c],cfg['sizes'][c],0,cfg['sizes'][c],rule);signals[c]=context['profiles']
        np.savez_compressed(root/f'f02_inputs/{c}.profiles.npz',**{f'{k}_{field}':v[field] for k,v in context['profiles'].items() for field in ['raw','smooth','background']})
        for par in base[c]:
            a,b=par['start']//STEP,(par['end']+STEP-1)//STEP
            fullarea=float(np.sum(context['profiles']['full']['raw'][a:b]));shortarea=float(np.sum(context['profiles']['short500']['raw'][a:b]))
            influence.append({'chrom':c,'baseline_id':par['cru_id'],'start':par['start'],'end':par['end'],
                'full_density_bin_sum':fullarea,'short500_density_bin_sum':shortarea,'long_gt500_span_fraction':1-shortarea/fullarea if fullarea else np.nan})
            local=[e for e in edges[c] if e['left']<par['end'] and e['right']>=par['start']]
            for method,minimum in [('V10_500kb',500000),('V11_500kb',500000),('V11_300kb',300000)]:
                rr,bb,cs,it=solve(par,sig,local,cfg['sizes'][c],rule,options(settings,minimum),acrs[c],genes[c],method,args.max_rounds)
                results[method]['crus'].extend(rr);results[method]['boundaries'].extend(bb);results[method]['candidates'].extend(cs);results[method]['iterations'].extend(it)
        for method in ['V10_500kb','V11_500kb','V11_300kb']:
            results[method]['boundaries'].extend({**r,'source':'inherited_B145','method':method} for r in original if r['chrom']==c)
        print('COUNT',c,{m:sum(r['chrom']==c for r in q['crus']) for m,q in results.items()},flush=True)
    write(root/'f01_diagnosis/06_long_span_contribution.tsv',influence)
    for method,folder in METHODS.items():
        out=root/'f03_results'/folder;out.mkdir();r=results[method]
        if method.startswith('V1'):
            for c in CHROMS:
                for j,q in enumerate(sorted((q for q in r['crus'] if q['chrom']==c),key=lambda q:q['start']),1):q['cru_id']=f'{c}_{method}_CRU_{j}'
            for c in CHROMS:
                for par in base[c]:
                    rr=sorted((q for q in r['crus'] if q['chrom']==c and par['start']<=q['start'] and q['end']<=par['end']),key=lambda q:q['start'])
                    assert rr[0]['start']==par['start'] and rr[-1]['end']==par['end'] and all(l['end']==u['start'] for l,u in zip(rr,rr[1:]))
        C.savebed(out/'02_CRU.bed',r['crus']);write(out/'03_boundaries.tsv',r['boundaries'],None if r['boundaries'] else ['chrom','boundary_pos'])
        write(out/'04_candidates.tsv.gz',r['candidates'],None if r['candidates'] else ['chrom','boundary_pos'])
        write(out/'05_selection_iterations.tsv',r['iterations'],None if r['iterations'] else ['chrom','baseline_id'])
    summary=[]
    for method,q in results.items():
        for c in CHROMS:
            rr=[r for r in q['crus'] if r['chrom']==c];row=P.interval_summary(method,c,rr,edges[c],cfg['sizes'][c])
            row['length_ge10mb_n']=sum(r['end']-r['start']>=10000000 for r in rr);summary.append(row)
        ll=[r['end']-r['start'] for r in q['crus']];rr=[r for r in summary if r['method']==method]
        total={'method':method,'chrom':'Tested_chromosomes','cru_n':len(ll),'covered_bp':sum(ll),
            'coverage_fraction':sum(ll)/sum(cfg['sizes'][c] for c in CHROMS),'median_bp':float(np.median(ll)),
            'mean_bp':float(np.mean(ll)),'max_bp':max(ll),'length_ge5mb_n':sum(v>=5000000 for v in ll),'length_ge10mb_n':sum(v>=10000000 for v in ll),
            'intra_aglink_n':sum(r['intra_aglink_n'] for r in rr),'aglink_n':sum(r['aglink_n'] for r in rr)}
        total['intra_aglink_fraction']=total['intra_aglink_n']/total['aglink_n'];summary.append(total)
        print('TOTAL',method,len(ll),flush=True)
    write(root/'f04_statistics/01_CRU_summary.tsv',summary)
    for r in hashes:r['md5_after']=C.md5(r['path']);assert r['md5_before']==r['md5_after']
    write(root/'f02_inputs/01_source_MD5.tsv',hashes)
    (root/'01_parameters.json').write_text(json.dumps({'history_root':str(hist),'max_rounds':args.max_rounds,'seed':args.seed,
        'sizes':cfg['sizes'],'block_bed':cfg['block_bed'],'m3_script':cfg['m3_script'],'aglink_bedpe':cfg['aglink_bedpe'],
        'methods':METHODS,'no_coACR':True,'count_target_used':False,
        'density_weights':'1/sqrt(max(1,distance/500kb))','morphology':'unchanged V9 shape gates applied to each profile',
        'local_density':'both endpoints inside current child; <=500kb, <=1Mb, mild weighted',
        'selection':'simultaneous maximum sum log(1+shape effect), compatible spacing; child-local joint revalidation'},indent=2)+'\n')
    for p in Path(__file__).parent.iterdir():
        if p.is_file():shutil.copy2(p,root/'f06_scripts'/p.name)
    (root/'02_CALLING.DONE').write_text('BEDs and counts saved. Only 1A/1B/2A; main 500 kb and 300 kb sensitivity are separate.\n')


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--history-root',type=Path,required=True);ap.add_argument('--outdir',type=Path,required=True)
    ap.add_argument('--max-rounds',type=int,default=3);ap.add_argument('--seed',type=int,default=20261108);main(ap.parse_args())
