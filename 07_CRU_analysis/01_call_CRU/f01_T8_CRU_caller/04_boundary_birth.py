#!/usr/bin/env python3
"""Refine outer/shared endpoints transactionally while a new territory is accepted.

No final BED is read or rewritten. Existing neighboring regions can be adjusted
transactionally only at the birth of a new shared junction, with revalidation.
"""

class BirthRefiner:
    def __init__(self,optimizer,functions,edges,raw,size,enabled,check_internal):
        self.opt=optimizer;self.edges=edges;self.raw=raw;self.size=size
        self.enabled=enabled;self.check_internal=check_internal;self.audit=[];self.functions=functions
        self._edge_index=functions.EdgeQueryIndex(edges)

    def pairs(self,a,b):
        return {(e['acr_id'],e['gene_id']) for e in self._edge_index.query(a,b)}

    def valid_rows(self,rows,affected):
        rows=sorted(rows,key=lambda r:r['start'])
        if any(a['end']>b['start'] for a,b in zip(rows,rows[1:])):return False,'overlap'
        checked_groups=set()
        for r in rows:
            if r['end']-r['start']<r.get('_minimum_bp',1):return False,'minimum_length'
            if '_gate' in r and r['cru_id'] in affected:
                group=r.get('_territory_id',r['cru_id'])
                if group in checked_groups:continue
                checked_groups.add(group)
                parts=[s for s in rows if s.get('_territory_id',s['cru_id'])==group]
                ok,why=r['_gate'](min(s['start'] for s in parts),max(s['end'] for s in parts))
                if not ok:return False,'frozen_territory_gate:'+why
        return self.check_internal(rows,affected)

    def accept(self,parent,existing,gate,stage,minimum,hard=None):
        """Return (new parent, updated neighbors) before subdivision/commit."""
        parent=dict(parent,_gate=gate,_minimum_bp=minimum)
        parent.setdefault('_territory_id',parent['cru_id'])
        parent.setdefault('_start_center',parent['start']);parent.setdefault('_end_center',parent['end'])
        if not self.enabled:return parent,list(existing)
        working=[dict(r) for r in existing]
        original=dict(parent);new_keep=self.pairs(original['start'],original['end'])
        for side in ['start','end']:
            p0=parent[side];near=[]
            for j,r in enumerate(working):
                d=p0-r['end'] if side=='start' else r['start']-p0
                if 0<=d<5000:near.append((d,j))
            neighbor_index=min(near)[1] if near else None
            neighbor=working[neighbor_index] if neighbor_index is not None else None
            neighbor_side='end' if side=='start' else 'start'
            center=int(neighbor.get('_'+neighbor_side+'_center',neighbor[neighbor_side])) if neighbor is not None else int(parent['_'+side+'_center'])
            fallback=int(neighbor[neighbor_side]) if neighbor is not None else int(parent[side])
            neighbor_keep=self.pairs(neighbor['start'],neighbor['end']) if neighbor is not None else set()
            neighbor_group=neighbor.get('_territory_id',neighbor['cru_id']) if neighbor is not None else None
            group_before=[r for r in working if r.get('_territory_id',r['cru_id'])==neighbor_group] if neighbor is not None else []
            group_keep=self.pairs(min(r['start'] for r in group_before),max(r['end'] for r in group_before)) if group_before else set()
            low=0 if side=='start' else parent['start']+minimum
            high=parent['end']-minimum if side=='start' else self.size
            if neighbor is not None:
                if side=='start':low=max(low,neighbor['start']+neighbor.get('_minimum_bp',1))
                else:high=min(high,neighbor['end']-neighbor.get('_minimum_bp',1))
            else:
                left=[r['end'] for r in working if r['end']<=parent['start']]
                right=[r['start'] for r in working if r['start']>=parent['end']]
                low=max(low,max(left) if left else 0);high=min(high,min(right) if right else self.size)
            hard_index=hard or self.opt.hard
            functions=vars(self.functions)
            ranges=list(functions['legal_ranges'](center,low,high,0,self.opt.radius,self.opt.acr,self.opt.gene,hard_index))
            context_start=neighbor['start'] if neighbor is not None and side=='start' else max(0,parent['start']-500000)
            context_end=neighbor['end'] if neighbor is not None and side=='end' else min(self.size,parent['end']+500000)
            local=self._edge_index.query(context_start,context_end)
            def count(point,candidate):
                if neighbor is not None:return functions['crossing'](local,point)
                return sum(min(e['p1'],e['p2'])<point<=max(e['p1'],e['p2']) and
                    ((candidate['start']<=e['p1']<candidate['end']) or (candidate['start']<=e['p2']<candidate['end'])) for e in local)
            def tentative(point):
                par=dict(parent);par[side]=int(point);others=[dict(r) for r in working]
                if neighbor_index is not None:others[neighbor_index][neighbor_side]=int(point)
                if not new_keep<=self.pairs(par['start'],par['end']):return None,None,'would_lose_original_internal_pairs'
                if neighbor_index is not None and not neighbor_keep<=self.pairs(others[neighbor_index]['start'],others[neighbor_index]['end']):
                    return None,None,'would_lose_neighbor_internal_pairs'
                if group_before:
                    group_after=[r for r in others if r.get('_territory_id',r['cru_id'])==neighbor_group]
                    if not group_keep<=self.pairs(min(r['start'] for r in group_after),max(r['end'] for r in group_after)):
                        return None,None,'would_lose_neighbor_territory_internal_pairs'
                affected={par['cru_id']}
                if neighbor_index is not None:affected.add(others[neighbor_index]['cru_id'])
                ok,why=self.valid_rows(others+[par],affected)
                if not ok:return None,None,why
                return par,others,''
            positions=functions['equivalence_points'](center,ranges,self.edges,context_start,context_end,self.opt.acr_pos,self.opt.gene_pos)
            ranking=[]
            for p in positions:
                candidate=dict(parent);candidate[side]=p
                ranking.append((count(p,candidate),abs(p-center),p))
            selected=None;trials=0;failure=''
            for cross,_,p in sorted(ranking):
                trials+=1;par,others,why=tentative(p)
                if par is not None:selected=par,others,cross;break
                if not failure:failure=why
            status='accepted_fine_outer_or_shared'
            if selected is None:
                par,others,why=tentative(fallback)
                if par is not None:selected=par,others,count(fallback,par);status='kept_previous_coordinate_no_suitable_fine_point'
                else:
                    # A tiny availability-bin join may fail a frozen gate.
                    # Keep the original coarse candidate and do not edit its neighbor.
                    par=dict(parent);others=list(working);selected=par,others,count(par[side],par)
                    status='kept_coarse_coordinate_join_not_legal';failure=failure or why
                    center=int(parent['_'+side+'_center'])
            par,others,cross=selected
            q=int(par[side])
            self.audit.append(dict(stage=stage,source='shared_at_territory_birth' if neighbor is not None else 'outer_at_territory_birth',
                side=side,candidate_center=center,previous_coordinate=p0,refined_pos=q,shift_bp=q-center,
                cross_before=count(p0,parent),cross_after=cross,local_context_start=context_start,local_context_end=context_end,
                internal_pairs_coarse=len(new_keep),internal_pairs_after=len(self.pairs(par['start'],par['end'])),
                neighbor_internal_pairs_before=len(neighbor_keep),neighbor_internal_pairs_after=len(self.pairs(others[neighbor_index]['start'],others[neighbor_index]['end'])) if neighbor_index is not None else 0,
                original_ACR_body_overlap_n=self.opt.acr.split_body_n(p0),original_gene_body_overlap_n=self.opt.gene.split_body_n(p0),
                refined_ACR_body_overlap_n=self.opt.acr.split_body_n(q),refined_gene_body_overlap_n=self.opt.gene.split_body_n(q),
                status=status,trials=trials,first_rejection_reason=failure,minimum_new_region_bp=minimum))
            par['_'+side+'_center']=center
            if neighbor_index is not None and status!='kept_coarse_coordinate_join_not_legal':others[neighbor_index]['_'+neighbor_side+'_center']=center
            parent,working=par,others
        ok,why=gate(parent['start'],parent['end'])
        assert ok,'Birth fallback must retain a valid coarse region: '+why
        return parent,working
