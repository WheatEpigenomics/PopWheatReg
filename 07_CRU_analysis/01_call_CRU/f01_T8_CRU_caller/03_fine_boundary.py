#!/usr/bin/env python3
"""Fine localization of every boundary before it is accepted; never edits a BED."""
import bisect
from collections import OrderedDict
import numpy as np


class IntervalIndex:
    def __init__(self, intervals):
        self.intervals=sorted((int(a),int(b)) for a,b in intervals if b>a)
        self.starts=[a for a,b in self.intervals]
        self.prefix_ends=np.maximum.accumulate([b for a,b in self.intervals]).tolist() if self.intervals else []

    def overlapping(self,lo,hi):
        j=bisect.bisect_left(self.prefix_ends,lo)
        while j<len(self.intervals) and self.intervals[j][0]<=hi:
            a,b=self.intervals[j]
            if b>=lo:yield a,b
            j+=1

    def split_body_n(self,p):
        return sum(a<p<b for a,b in self.overlapping(p,p))


def local_edges(edges,start,end):
    unique={}
    for e in edges:
        if start<=e['p1']<end and start<=e['p2']<end:
            unique[(e['acr_id'],e['gene_id'])]=e
    return list(unique.values())


class EdgeQueryIndex:
    """Exact half-open endpoint queries, retaining original order and pair deduplication."""
    def __init__(self,edges,cache_size=1024):
        self.edges=edges;self.cache_size=cache_size;self.cache=OrderedDict()
        self.order=np.argsort([min(e['p1'],e['p2']) for e in edges],kind='stable')
        self.left=np.array([min(edges[i]['p1'],edges[i]['p2']) for i in self.order])

    def query(self,start,end):
        key=(start,end)
        if key in self.cache:
            self.cache.move_to_end(key)
            return list(self.cache[key])
        lo,hi=np.searchsorted(self.left,[start,end],side='left')
        indexes=sorted(int(i) for i in self.order[lo:hi]
                       if max(self.edges[i]['p1'],self.edges[i]['p2'])<end)
        unique={}
        for i in indexes:
            e=self.edges[i];unique[(e['acr_id'],e['gene_id'])]=e
        result=tuple(unique.values());self.cache[key]=result
        if len(self.cache)>self.cache_size:self.cache.popitem(last=False)
        return list(result)


def crossing(edges,p):
    return sum(min(e['p1'],e['p2'])<p<=max(e['p1'],e['p2']) for e in edges)


def legal_ranges(center,start,end,minimum,radius,acr,gene,hard):
    lo=max(start+minimum,center-radius);hi=min(end-minimum,center+radius)
    if lo>hi:return []
    forbidden=[]
    for index in [acr,gene]:
        forbidden.extend((a+1,b-1) for a,b in index.overlapping(lo,hi) if b-a>1)
    for a,b in hard.overlapping(min(center,lo),max(center,hi)):
        if a<center<b:return []
        if b<=center:lo=max(lo,b)
        elif a>=center:hi=min(hi,a)
    for a,b in sorted(forbidden):
        if b<lo:continue
        if a>hi:break
        if a>lo:yield_lo=lo;yield_hi=min(hi,a-1);yield yield_lo,yield_hi
        lo=max(lo,b+1)
        if lo>hi:return
    if lo<=hi:yield lo,hi


def equivalence_points(center,ranges,edges,start,end,acr_pos,gene_pos):
    """Partition all 1-bp choices at endpoint/profile/constraint change events.

    Only one nearest-to-center representative per constant-evidence cell is
    needed. Events include the frozen 250/500-kb contexts and crossing-profile
    sampling offsets; q and q+1 retain half-open endpoint transitions.
    """
    if not ranges:return []
    lo=min(a for a,b in ranges);hi=max(b for a,b in ranges)
    nodes=np.array([e[k] for e in edges for k in ['p1','p2']],dtype=np.int64)
    shifts={0,250000,-250000,500000,-500000,150000,-150000,50000,-50000,12500,-12500}
    samples=list(range(-147500,-49999,5000))+list(range(52500,150001,5000))+list(range(-10000,10001,5000))
    shifts.update(-x for x in samples)
    events={lo,hi+1}
    for shift in shifts:
        for delta in [0,1]:
            v=nodes+shift+delta
            events.update(int(p) for p in v[(v>=lo)&(v<=hi+1)])
    annotated=np.r_[acr_pos,gene_pos]
    for shift in [0,250000,-250000,500000,-500000]:
        for delta in [0,1]:
            v=annotated+shift+delta
            events.update(int(p) for p in v[(v>=lo)&(v<=hi+1)])
    for q in [start,end]:
        for shift in shifts:
            for delta in [0,1]:
                if lo<=q+shift+delta<=hi+1:events.add(q+shift+delta)
    events.update(range((lo//5000)*5000,hi+5001,5000))
    points=set()
    for a,b in ranges:
        cut=sorted({a,b+1}|{p for p in events if a<p<=b})
        for x,y in zip(cut,cut[1:]):points.add(max(x,min(center,y-1)))
    return sorted(points)


class Refiner:
    def __init__(self,acr,gene,hard,acr_pos,gene_pos,radius=10000,enabled=True):
        self.acr=IntervalIndex(acr);self.gene=IntervalIndex(gene);self.hard=IntervalIndex(hard)
        self.acr_pos=np.asarray(acr_pos);self.gene_pos=np.asarray(gene_pos)
        self.radius=radius;self.enabled=enabled;self.audit=[]
        self._edge_indexes=OrderedDict()

    def query_edges(self,edges,start,end):
        key=id(edges)
        if key not in self._edge_indexes:
            self._edge_indexes[key]=EdgeQueryIndex(edges)
            if len(self._edge_indexes)>8:self._edge_indexes.popitem(last=False)
        self._edge_indexes.move_to_end(key)
        return self._edge_indexes[key].query(start,end)

    def choose(self,center,start,end,minimum,edges,evaluate,stage,source,context_tag=''):
        original=evaluate(center)
        if not self.enabled:return original
        local=self.query_edges(edges,start,end)
        info=dict(stage=stage,source=source,candidate_center=center,parent_start=start,parent_end=end,
            minimum_daughter_bp=minimum,context_tag=context_tag,cross_before=crossing(local,center),
            local_unique_aglink_n=len(local),original_ACR_body_overlap_n=self.acr.split_body_n(center),
            original_gene_body_overlap_n=self.gene.split_body_n(center),original_eligible=int(original.get('eligible',0)))
        shape_failure=original.get('reject_reason','').startswith('shape:') or any(
            x in original.get('reject_reason','').split(',') for x in ['fine_depth_or_width','multiscale_depth_or_width'])
        if shape_failure:
            r=dict(original);r['eligible']=0
            self.audit.append(dict(info,refined_pos='',cross_after='',status='original_candidate_shape_rejected',trials=0,
                                   original_reject_reason=original.get('reject_reason',''),refined_reject_reason=r.get('reject_reason','')))
            return r
        ranges=list(legal_ranges(center,start,end,minimum,self.radius,self.acr,self.gene,self.hard))
        positions=equivalence_points(center,ranges,edges,start,end,self.acr_pos,self.gene_pos)
        l=np.sort([min(e['p1'],e['p2']) for e in local]);rr=np.sort([max(e['p1'],e['p2']) for e in local])
        ps=np.array(positions,dtype=np.int64)
        counts=np.searchsorted(l,ps,side='left')-np.searchsorted(rr,ps,side='left')
        order=sorted(zip(positions,counts),key=lambda v:(int(v[1]),abs(v[0]-center),v[0]))
        trials=0;first_failure='';selected=None
        for p,count in order:
            trials+=1;r=evaluate(p)
            if not r.get('eligible'):
                if not first_failure:first_failure=r.get('reject_reason','unknown')
                continue
            assert self.acr.split_body_n(p)==0 and self.gene.split_body_n(p)==0
            selected=dict(r,candidate_center=center,boundary_pos=p,fine_cross_before=info['cross_before'],
                fine_cross_after=int(count),fine_shift_bp=p-center,fine_annotation_legal=1)
            self.audit.append(dict(info,refined_pos=p,cross_after=int(count),status='accepted_refined' if p!=center else 'accepted_original_legal_optimum',
                trials=trials,geometric_cell_n=len(order),shift_bp=p-center,refined_ACR_body_overlap_n=0,
                refined_gene_body_overlap_n=0,original_reject_reason=original.get('reject_reason',''),refined_reject_reason=''))
            break
        if selected is not None:return selected
        if original.get('eligible') and min(center-start,end-center)>=minimum and not any(a<=center<b for a,b in self.hard.overlapping(center,center)):
            kept=dict(original,candidate_center=center,boundary_pos=center,fine_cross_before=info['cross_before'],
                fine_cross_after=info['cross_before'],fine_shift_bp=0,fine_annotation_legal=int(
                    not self.acr.split_body_n(center) and not self.gene.split_body_n(center)))
            self.audit.append(dict(info,refined_pos=center,cross_after=info['cross_before'],shift_bp=0,
                status='kept_original_no_suitable_fine_point',trials=trials,geometric_cell_n=len(order),
                refined_ACR_body_overlap_n=self.acr.split_body_n(center),refined_gene_body_overlap_n=self.gene.split_body_n(center),
                original_reject_reason=original.get('reject_reason',''),refined_reject_reason=first_failure or 'annotation_window_blocked'))
            return kept
        reason='no_annotation_legal_position' if not order else 'no_legal_position_passing_frozen_split_rules'
        rejected=dict(original,eligible=0,accepted=0,reject_reason=reason,candidate_center=center)
        self.audit.append(dict(info,refined_pos='',cross_after='',status=reason,trials=trials,
            geometric_cell_n=len(order),original_reject_reason=original.get('reject_reason',''),refined_reject_reason=first_failure))
        return rejected
