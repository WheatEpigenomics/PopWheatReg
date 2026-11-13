#!/usr/bin/env python3
"""The complete CRU calling chain: one in-memory pass over the frozen rule set.

Every boundary is localized before its candidate is accepted or a parent is
split; no coordinate is post-processed afterwards.
"""
from types import SimpleNamespace as NS
import importlib.util
from pathlib import Path
import sys

import numpy as np
from scipy.signal import find_peaks


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0, str(HERE))

from f02_modules import m11_engine as E
from f02_modules import m13_low_density_recall as P
from f02_modules import m14_support_recall as S
from f02_modules import m15_island as T7

M, H = E.M, E.H


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, str(path))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


O = _load(HERE / '03_fine_boundary.py', 'fine_boundary')
BIRTH = _load(HERE / '04_boundary_birth.py', 'boundary_birth')


def evaluate_local(position,source,context,opt,acrs,genes,method):
    # Same frozen M.evaluate gates, with the bin containing the true coordinate.
    off=context['origin'];p=position-off;i=p//5000;start,end=context['start'],context['end']
    row=dict(boundary_pos=position,source=source,parent_start=start+off,parent_end=end+off,
        method=method,accepted=0,eligible=0,density_snr=0.,local_profile='none',reject_reason='')
    if min(p-start,end-p)<opt.min_daughter_bp:return dict(row,reject_reason='minimum_daughter_length')
    shapes=[]
    for name,sig in context['profiles'].items():
        if 'states' not in sig:sig['states']=M.T.density_states(sig,start,end)
        shape=M.T.sharp_shape(i,sig,start,end,sig['states'],opt) if source=='sharp' else M.T.level_shape(source,i,sig,start,end,sig['states'],opt)
        shapes.append(dict(shape,local_profile=name))
    eligible=[r for r in shapes if r['shape_pass']]
    if not eligible:
        best=max(shapes,key=lambda r:r.get('density_snr',0))
        return {**row,**best,'reject_reason':'shape:'+best['shape_reason']}
    local=[r for r in eligible if r['local_profile']!='full']
    chosen=max(local or eligible,key=lambda r:(r['density_snr'],r.get('105k_relative_change',0)))
    graph=M.F.graph(context['edges'],start,end,p,opt);hard=[]
    if min(graph['whole_left_internal'],graph['whole_right_internal'])<3:hard.append('whole_internal_support')
    if min(graph['whole_left_ACR_n'],graph['whole_right_ACR_n'])<2:hard.append('whole_ACR_diversity')
    short=[M.A.short_evidence(context['edges'],start,end,p,cap,acrs-off,genes-off) for cap in [250000,500000]]
    connected=[];modules=[]
    for d in short:
        support=min(d['left_internal'],d['right_internal'])>=2 and min(d['left_ACR_n'],d['right_ACR_n'])>=2
        change=M.A.profile_change(d,dict(source='step' if source=='changepoint' else source,direction=chosen.get('direction')))
        if support:modules.append(d)
        if support and (d['cross_fraction']<=.85 or d['relative_depletion']>=0 or (np.isfinite(change) and change>=.05)):connected.append(d)
    if not modules:hard.append('short_local_module_support')
    route='local_density_with_short_separation' if local and connected else ''
    if not route and not local and connected:route='full_reference_with_short_separation'
    if not route and method.startswith('V11') and local and any(min(d['left_internal'],d['right_internal'])>=3 for d in modules):
        route='local_density_with_two_sided_modules'
    if not route:hard.append('connection_evidence')
    out={**row,**chosen,**graph,'boundary_pos':position,'parent_start':start+off,'parent_end':end+off,
        'shape_pass':1,'eligible':int(not hard),'accepted':0,'acceptance_route':route,'reject_reason':','.join(hard),
        'local_profiles_passing':','.join(r['local_profile'] for r in eligible),
        'full_shape_pass':int(any(r['local_profile']=='full' for r in eligible)),'shape_bin_center':(position//5000)*5000+2500}
    for d in short:out.update({f'short_{d["cap_bp"]}_{k}':v for k,v in d.items()})
    return out


def evaluate_primary(position,parent_start,parent_end,es,sig,opt,chrom,depth,template):
    D=H.D;a=max(0,int(np.ceil(parent_start/5000-.5)));b=min(len(sig['raw']),int(np.ceil(parent_end/5000-.5)))
    idx=position//5000
    peaks=[]
    for y in sig['smooth']:
        pp,props=find_peaks(y[a:b],prominence=opt.min_peak_prominence)
        absolute=pp+a
        peaks.append(absolute[props['prominences']>=opt.min_relative_peak_prominence*np.maximum(y[absolute],1e-12)])
    scores=[D.score_scale(y,idx,a,b,pk,opt) for y,pk in zip(sig['smooth'],peaks)]
    scale_pass=[s['depth']>=opt.valley_depth and s['width_bp']>=opt.valley_min_width_bp for s in scores]
    left=np.array([e['left'] for e in es]);right=np.array([e['right'] for e in es]);cross=(left<position)&(right>position)
    contexts=[]
    for flank in opt.support_flanks_bp:
        lo,hi=max(parent_start,position-flank),min(parent_end,position+flank)
        il,ir=D.internal_n(left,right,lo,position),D.internal_n(left,right,position,hi)
        local=int(np.count_nonzero(cross&(left>=lo)&(right<hi)));weak=min(il,ir)
        contexts.append((flank,il,ir,local,local/(local+weak) if local+weak else 1))
    good=[v for v in contexts if min(v[1],v[2])>=opt.min_internal and v[4]<=opt.max_local_cross_fraction]
    flank,il,ir,local,fraction=good[0] if good else contexts[-1]
    fail=[]
    if min(position-parent_start,parent_end-position)<opt.min_cru_bp:fail.append('daughter_length')
    if not scale_pass[0]:fail.append('fine_depth_or_width')
    if sum(scale_pass)<opt.min_scale_support:fail.append('multiscale_depth_or_width')
    if min(il,ir)<opt.min_internal:fail.append('local_internal_support')
    if fraction>opt.max_local_cross_fraction:fail.append('local_crossing_fraction')
    if min(D.internal_n(left,right,parent_start,position),D.internal_n(left,right,position,parent_end))<opt.min_internal:fail.append('daughter_internal_support')
    return dict(template,boundary_pos=position,eligible=int(not fail),reject_reason=','.join(fail),
        scale_support_n=sum(scale_pass),score=sum(s['depth'] for s in scores if s['depth']>=opt.valley_depth),
        left_internal=il,right_internal=ir,crossing_local=local,local_cross_fraction=fraction)


def primary(chrom,start,end,es,sig,opt,refiner):
    if not refiner.enabled:return H.C.recursive(chrom,start,end,es,sig,opt,H.D,'primary45')
    todo=[(start,end,0)];result={k:[] for k in ['crus','boundaries','candidates','actions']}
    left=np.array([e['left'] for e in es]);right=np.array([e['right'] for e in es])
    while todo:
        a,b,depth=todo.pop();original=H.D.parent_candidates(chrom,a,b,left,right,sig,opt,depth);rows=[]
        for old in original:
            center=int(old['boundary_pos'])
            row=refiner.choose(center,a,b,opt.min_cru_bp,es,
                lambda p: evaluate_primary(p,a,b,es,sig,opt,chrom,depth,old),
                'B145_primary',old.get('source','valley'),f'depth={depth}')
            rows.append(row)
        result['candidates'].extend(rows);valid=[r for r in rows if r['eligible']]
        if valid and depth<opt.max_recursion_depth:
            best=max(valid,key=lambda r:(r['scale_support_n'],r['score'],r['prominence'],-r['boundary_pos']))
            best=dict(best,accepted=1);result['boundaries'].append(best);p=int(best['boundary_pos'])
            todo.extend([(p,b,depth+1),(a,p,depth+1)])
            result['actions'].append(dict(chrom=chrom,start=a,end=b,action='fine_point_used_before_recursion',boundary_pos=p))
        else:result['crus'].append(dict(chrom=chrom,start=a,end=b,original_parent_start=start,original_parent_end=end))
    return result


def solve(parent,sig,edges,size,rule,opt,acrs,genes,method,max_rounds,refiner):
    if not refiner.enabled:return M.solve(parent,sig,edges,size,rule,opt,acrs,genes,method,max_rounds)
    start,end=parent['start'],parent['end'];pool=set();current={};seen=set();records=[];iterations=[];final=[]
    if end-start<2*opt.min_daughter_bp:return [dict(parent)],[],[],[]
    for iteration in range(max_rounds):
        points=[start]+sorted(int(r['boundary_pos']) for r in current.values())+[end];contexts={}
        def ctx(a,b):
            if (a,b) not in contexts:contexts[(a,b)]=M.profile_context(sig,edges,size,a,b,rule)
            return contexts[(a,b)]
        def refine(key,a,b):
            center,source=key;context=ctx(a,b)
            r=refiner.choose(center,a,b,opt.min_daughter_bp,edges,
                lambda p:evaluate_local(p,source,context,opt,acrs,genes,method),method,source,f'iteration={iteration}')
            return dict(r,chrom=parent['chrom'],baseline_id=parent['cru_id'],iteration=iteration,
                        candidate_center=center,proposal_source=source)
        for a,b in zip(points,points[1:]):
            if b-a>=2*opt.min_daughter_bp:pool.update(M.proposals(ctx(a,b),opt))
        evaluations={}
        for key in sorted(pool):
            anchor=int(current[key]['boundary_pos']) if key in current else key[0]
            if anchor in points[1:-1]:
                j=points.index(anchor);a,b=points[j-1],points[j+1]
            else:
                j=int(np.searchsorted(points,anchor))-1;j=max(0,min(j,len(points)-2));a,b=points[j],points[j+1]
            row=refine(key,a,b);evaluations[key]=row;records.append(row)
        valid={k:r for k,r in evaluations.items() if r['eligible']};blocked=set();stable=[]
        for retry in range(12):
            chosen=M.weighted_selection([r for k,r in valid.items() if k not in blocked],start,end,opt.min_daughter_bp)
            pp=[start]+[int(r['boundary_pos']) for r in chosen]+[end]
            updated=[];changed=False
            for j,row in enumerate(chosen,1):
                key=(int(row['candidate_center']),row['proposal_source'])
                q=refine(key,pp[j-1],pp[j+1]);records.append(q)
                if not q['eligible']:blocked.add(key);changed=True;continue
                valid[key]=q;updated.append(q)
                if q['boundary_pos']!=row['boundary_pos']:changed=True
            if not changed:
                stable=updated;break
        else:
            # Reject unresolved simultaneous fine coordinates; never silently
            # accept a point with stale neighboring-region constraints.
            stable=[]
            iterations.append(dict(chrom=parent['chrom'],baseline_id=parent['cru_id'],iteration=iteration,
                                   action='joint_fine_coordinates_did_not_converge'))
        final=stable;new={(int(r['candidate_center']),r['proposal_source']):r for r in final}
        state=tuple(sorted((k,int(r['boundary_pos'])) for k,r in new.items()))
        oldstate=tuple(sorted((k,int(r['boundary_pos'])) for k,r in current.items()))
        iterations.append(dict(chrom=parent['chrom'],baseline_id=parent['cru_id'],iteration=iteration,
            candidate_n=len(pool),eligible_n=len(valid),selected_n=len(final),action='joint_fine_revalidation'))
        if state==oldstate:break
        if state in seen:
            final=[];iterations.append(dict(chrom=parent['chrom'],baseline_id=parent['cru_id'],iteration=iteration,action='fine_coordinate_cycle_rejected'));break
        seen.add(state);current=new
    final=sorted([dict(r,accepted=1,reject_reason='accepted_'+r['acceptance_route']) for r in final],key=lambda r:r['boundary_pos'])
    points=[start]+[int(r['boundary_pos']) for r in final]+[end]
    assert all(b-a>=opt.min_daughter_bp for a,b in zip(points,points[1:]))
    return [dict(parent,start=a,end=b) for a,b in zip(points,points[1:])],final,records,iterations


def run_chain(chrom,size,es,rules,t3,t6,islands,acrs,genes,acr_iv,gene_iv,enabled):
    raw=M.P.signals(es,size,NS(**rules['signal']))
    sig=M.D.territory_signals(raw['raw'],np.asarray(raw['smooth']),raw['endpoint_broad'],NS(**rules['primary']))
    sig['old_active']=raw['active'];protected=np.zeros(len(sig['raw']),bool);hard=[]
    for a,b in P.runs(sig['raw']<=1e-9):
        if min(size,b*5000)-a*5000>=100000:protected[a:b]=True;hard.append((a*5000,min(size,b*5000)))
    refiner=O.Refiner(acr_iv,gene_iv,hard,acrs,genes,enabled=enabled)
    hs,mask,gaps,cores=H.C.identify_gaps(chrom,size,sig,es,NS(**rules['hard_gap']),H.D)
    parents,actions=H.C.base_components(chrom,size,mask,sig,es,NS(**rules['hard_gap']),H.D)
    refined,events=H.E.refine_edges(chrom,sig,es,gaps,cores,NS(**rules['A_edge']))
    original=[dict(chrom=chrom,start=r['start'],end=r['end'],cru_id=r['territory_id']) for r in parents]
    expanded,_=H.E.apply_edges(original,refined);coarse_parents=[dict(r) for r in expanded]
    primary_rows=[];boundaries=[];origins=[];candidate_records=[];cut_checks={}
    index=S.EdgeIndex(es)
    def metric(a,b,opt,low_bg=False,island=False):
        lo=max(0,a//5000);hi=min(len(sig['raw']),(b+4999)//5000)
        row=(T7.shape(lo,hi,sig['raw'],sig['smooth'],sig['background'],size,islands) if island else
             P.describe(sig['raw'],sig['smooth'],sig['background'],sig['endpoint'],lo,hi,size))
        row.update(start=a,end=b,length_bp=b-a);row.update(index.describe(a,b))
        if island:
            inside=[e for e in index.unique if a<=e['left'] and e['right']<b]
            row.update(single_aglink_span_bp=inside[0]['right']-inside[0]['left'] if len(inside)==1 else 0,
                protected_zero_overlap_bp=int(protected[lo:hi].sum()*5000),coacr_internal_n=0,coacr_acr_n=0,
                coacr_linked_component_edges=0,coacr_linked_component_nodes=0,coacr_short500_n=0,coacr_endpoint_extent_bp=0)
            return row
        threshold=np.maximum(opt['absolute_floor'],sig['background'][lo:hi]+opt['background_sigma']*np.sqrt(sig['background'][lo:hi]+1))
        supported=(sig['raw'][lo:hi]>=opt['absolute_floor']-1e-9)&(np.sum(sig['smooth'][:,lo:hi]>=threshold-1e-9,axis=0)>=2)
        if low_bg:supported &= sig['background'][lo:hi]<=opt['near_zero_background']
        ii=np.r_[np.arange(max(0,lo-50),lo),np.arange(hi,min(len(sig['raw']),hi+50))]
        row.update(supported_bin_fraction=float(supported.mean()),
            background_near_zero_fraction=float(np.mean(sig['background'][lo:hi]<=opt.get('near_zero_background',.1))),
            protected_zero_overlap_bp=int(protected[lo:hi].sum()*5000),
            flank_background_median=float(np.median(sig['background'][ii])) if len(ii) else 0.,
            fine_peak_excess=float(np.max(sig['smooth'][0,lo:hi]-sig['background'][lo:hi])),
            peak_to_flank_ratio_regularized=(float(sig['smooth'][0,lo:hi].max())+.25)/(float(np.median(sig['smooth'][0,ii]))+.25) if len(ii) else (float(sig['smooth'][0,lo:hi].max())+.25)/.25)
        return row
    def gate_for(method,opt,low_bg=False,rescue=False):
        def check(a,b):
            if a>=b:return False,'invalid_interval'
            row=metric(a,b,opt,low_bg,method=='T7_A')
            if method=='T3':
                # Frozen T3 supported fraction does not add a separate low-bg mask.
                lo=a//5000;hi=(b+4999)//5000
                th=np.maximum(opt['absolute_floor'],sig['background'][lo:hi]+opt['background_sigma']*np.sqrt(sig['background'][lo:hi]+1))
                fraction=float(np.mean((sig['raw'][lo:hi]>=opt['absolute_floor']-1e-9)&(np.sum(sig['smooth'][:,lo:hi]>=th-1e-9,axis=0)>=2)))
                fail=P.reasons_for(row,opt,fraction)
            elif method=='T7_A':fail,_=T7.decide(row,islands)
            else:fail=S.assess(row,t6[method],rescue)
            return not fail,','.join(fail)
        return check
    def base_gate(a,b):
        n=sum(a<=e['left'] and e['right']<b for e in es)
        fail=[]
        if b-a<rules['hard_gap']['min_cru_bp']:fail.append('component_min_length')
        if n<rules['hard_gap']['min_internal']:fail.append('component_internal_links')
        return not fail,','.join(fail)
    def check_internal(rows,affected):
        for left,right in zip(rows,rows[1:]):
            q=left['end']
            if q!=right['start'] or q not in cut_checks or not ({left['cru_id'],right['cru_id']} & affected):continue
            info=cut_checks[q];a,b=left['start'],right['end']
            if info['kind']=='primary':
                r=evaluate_primary(q,a,b,es,sig,NS(**rules['primary']),chrom,0,info['template'])
            else:
                context=M.profile_context(sig,es,size,a,b,NS(**rules['primary']))
                r=evaluate_local(q,info['source'],context,M.options(rules['shape_settings'],300000),acrs,genes,info['method'])
            if not r['eligible']:return False,'neighbor_internal_cut:'+r.get('reject_reason','unknown')
        return True,''
    birth=BIRTH.BirthRefiner(refiner,O,es,sig['raw'],size,enabled,check_internal)
    def lineage(rr,par,bb,kind,method):
        centers={int(r['boundary_pos']):int(r.get('candidate_center',r['boundary_pos'])) for r in bb}
        for r in rr:
            r['_gate']=par.get('_gate',base_gate)
            r['_territory_id']=par.get('_territory_id',par['cru_id'])
            r['_minimum_bp']=300000 if kind!='primary' and len(rr)>1 else min(par.get('_minimum_bp',100000),r['end']-r['start'])
            r['_start_center']=par.get('_start_center',par['start']) if r['start']==par['start'] else centers.get(r['start'],r['start'])
            r['_end_center']=par.get('_end_center',par['end']) if r['end']==par['end'] else centers.get(r['end'],r['end'])
        for r in bb:cut_checks[int(r['boundary_pos'])]=dict(kind=kind,source=r.get('source','valley'),method=method,template=r)
        return rr
    root_parents=[]
    for original_par in expanded:
        par,root_parents=birth.accept(original_par,root_parents,base_gate,'Aedge_territory',100000,
            hard=O.IntervalIndex(hard+[(r['start'],r['end']) for r in refined]))
        root_parents.append(par)
    for par in root_parents:
        origins.extend(dict(chrom=chrom,position=p,source='Aedge_territory_outer',refinement_at_birth=1) for p in [par['start'],par['end']])
        out=primary(chrom,par['start'],par['end'],es,sig,NS(**rules['primary']),refiner)
        primary_rows.extend(lineage(out['crus'],par,out['boundaries'],'primary','B145_primary'))
        boundaries.extend(out['boundaries']);candidate_records.extend(out['candidates'])
    primary_rows.sort(key=lambda r:r['start'])
    print('CRU' if enabled else 'CONTROL','primary components',len(primary_rows),flush=True)
    for i,r in enumerate(primary_rows,1):r['cru_id']=f'{chrom}_Aedge_B145_CRU_{i}'
    rows=[];iteration_records=[]
    def subdivide(par,method):
        local=[e for e in es if e['left']<par['end'] and e['right']>=par['start']]
        rr,bb,cc,it=solve(par,sig,local,size,NS(**rules['primary']),M.options(rules['shape_settings'],300000),acrs,genes,
            method,rules['max_rounds'],refiner)
        boundaries.extend(bb);candidate_records.extend(cc);iteration_records.extend(it)
        return lineage(rr,par,bb,'local',method)
    for j,par in enumerate(primary_rows,1):
        rows.extend(subdivide(par,'V11_300kb'))
        if enabled and j%20==0:print('CRU local parents completed',j,'/',len(primary_rows),flush=True)
    rows.sort(key=lambda r:r['start'])
    for i,r in enumerate(rows,1):r['cru_id']=f'{chrom}_V11_300kb_CRU_{i}'
    baseline_coverage=sum(r['end']-r['start'] for r in rows)
    print('CRU' if enabled else 'CONTROL','local subdivision CRUs',len(rows),flush=True)
    recall_records=[]
    for method in ['T3','T4','T5','T6']:
        if method=='T3':routes=[('T2',t3['T2'],False),('T3_rescue',t3['T3_rescue'],True)]
        else:
            routes=[('T2_density',t3['T2'],False),('T3_low_background',t3['T3_rescue'],True)]
            if method=='T6':routes.append(('T6_low_background_peak',t6['T6_density'],True))
        new_n=0
        for route,opt,low_bg in routes:
            available=~P.mask_intervals(rows,size)
            if method!='T3':available &= ~protected
            proposed=P.proposal_mask(sig['raw'],sig['smooth'],sig['background'],available,opt,low_bg)
            for a,b in P.runs(proposed):
                row=dict(method=method,stage=route,chrom=chrom,**P.describe(sig['raw'],sig['smooth'],sig['background'],sig['endpoint'],a,b,size))
                threshold=np.maximum(opt['absolute_floor'],sig['background'][a:b]+opt['background_sigma']*np.sqrt(sig['background'][a:b]+1))
                supported=(sig['raw'][a:b]>=opt['absolute_floor']-1e-9)&(np.sum(sig['smooth'][:,a:b]>=threshold-1e-9,axis=0)>=2)
                if method=='T3':
                    row.update(P.support(es,row['start'],row['end']));row['supported_bin_fraction']=float(supported.mean())
                    fail=P.reasons_for(row,opt,row['supported_bin_fraction'])
                else:
                    if low_bg:supported &= sig['background'][a:b]<=opt['near_zero_background']
                    row.update(index.describe(row['start'],row['end']))
                    flank=250000//5000;ii=np.r_[np.arange(max(0,a-flank),a),np.arange(b,min(len(sig['raw']),b+flank))]
                    row.update(supported_bin_fraction=float(supported.mean()),
                        background_near_zero_fraction=float(np.mean(sig['background'][a:b]<=opt.get('near_zero_background',.1))),
                        protected_zero_overlap_bp=int(protected[a:b].sum()*5000),
                        flank_background_median=float(np.median(sig['background'][ii])) if len(ii) else 0.,
                        fine_peak_excess=float(np.max(sig['smooth'][0,a:b]-sig['background'][a:b])),
                        peak_to_flank_ratio_regularized=(float(sig['smooth'][0,a:b].max())+.25)/(float(np.median(sig['smooth'][0,ii]))+.25) if len(ii) else (float(sig['smooth'][0,a:b].max())+.25)/.25)
                    fail=S.assess(row,t6[method],route=='T6_low_background_peak')
                row.update(accepted=int(not fail),reject_reason=','.join(fail));recall_records.append(row)
                if fail:continue
                new_n+=1;par=dict(chrom=chrom,start=row['start'],end=row['end'],cru_id=f'{chrom}_{method}_recalled_{route}_{new_n}')
                gate=gate_for(method,opt,low_bg,route=='T6_low_background_peak')
                par,rows=birth.accept(par,rows,gate,method+':'+route,300000)
                origins.extend(dict(chrom=chrom,position=p,source=method+'_recall_outer',refinement_at_birth=1) for p in [par['start'],par['end']])
                rows.extend(subdivide(par,method))
        rows.sort(key=lambda r:r['start']);assert P.disjoint(rows)
        print('CRU' if enabled else 'CONTROL',method,'CRUs',len(rows),flush=True)
    # T7-A adds islands only; no coACR/B rescue is read or executed.
    available=P.mask_intervals(rows,size);added=[];island_n=0
    for a,b in T7.generate(sig['raw'],sig['smooth'],sig['background'],size,available,protected,islands):
        row=dict(chrom=chrom,**T7.shape(a,b,sig['raw'],sig['smooth'],sig['background'],size,islands))
        row.update(index.describe(row['start'],row['end']));inside=[e for e in index.unique if row['start']<=e['left'] and e['right']<row['end']]
        row.update(single_aglink_span_bp=inside[0]['right']-inside[0]['left'] if len(inside)==1 else 0,
            protected_zero_overlap_bp=int(protected[a:b].sum()*5000),
            coacr_internal_n=0,coacr_acr_n=0,coacr_linked_component_edges=0,coacr_linked_component_nodes=0,coacr_short500_n=0,coacr_endpoint_extent_bp=0)
        fail,_=T7.decide(row,islands);row.update(method='T7_A',accepted=int(not fail),reject_reason=','.join(fail));recall_records.append(row)
        if fail:continue
        island_n+=1;par=dict(chrom=chrom,start=row['start'],end=row['end'],cru_id=f'{chrom}_T7_A_island_{island_n}')
        par,rows=birth.accept(par,rows+added,gate_for('T7_A',islands), 'T7_A_island',300000)
        # Existing island additions are now in the transactional working list.
        added=[];rows.append(par)
        origins.extend(dict(chrom=chrom,position=p,source='T7_A_island_outer',refinement_at_birth=1) for p in [par['start'],par['end']])
    rows=sorted(rows+added,key=lambda r:r['start']);assert P.disjoint(rows)
    return dict(crus=rows,boundaries=boundaries,origins=origins,audit=refiner.audit,outer_audit=birth.audit,candidates=candidate_records,
        iterations=iteration_records,recall=recall_records,sig=sig,refined_gaps=refined,parents=coarse_parents,
        baseline_coverage=baseline_coverage,protected=protected)

