# m11_engine: V11 engine (M/H namespaces)
# Frozen rule module, copied unchanged from the published calling chain
# (source: f17_V11_pipeline/05_call_V11.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Whole-chromosome wrapper around the frozen A-edge+B145 and V11-300kb functions."""
import argparse
import contextlib
import csv
import hashlib
import importlib.util
import json
import multiprocessing
from pathlib import Path
import pickle
import shutil
from types import SimpleNamespace
import numpy as np


def load(path,name):
    s=importlib.util.spec_from_file_location(name,path);m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


HERE=Path(__file__).resolve().parent
from . import m10_local_subdivision as M
from . import m08_baseline_targets as H
CTX={}


def md5(path):
    h=hashlib.md5()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def savebed(path,rows):
    with path.open('w') as f:
        for r in rows:f.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['cru_id']}\n")


def worker(chrom):
    root=CTX['root'];d=root/'chromosomes'/chrom;d.mkdir(parents=True,exist_ok=True)
    done=d/'DONE.json';statefile=d/'state.pickle'
    if done.exists() and CTX['resume']:
        old=json.loads(done.read_text())
        if old['signature']==CTX['signature'] and all(Path(p).is_file() and md5(p)==v for p,v in old['artifacts'].items()):
            return chrom,'reused'
    with (d/'run.log').open('a',buffering=1) as log,contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
        rule=CTX['rules'];size=CTX['sizes'][chrom];es=CTX['edges'][chrom]
        signal_args=SimpleNamespace(**rule['signal'])
        raw=M.P.signals(es,size,signal_args)
        sig=M.D.territory_signals(raw['raw'],np.asarray(raw['smooth']),raw['endpoint_broad'],SimpleNamespace(**rule['primary']))
        sig['old_active']=raw['active']
        baseline,aux=H.baseline(chrom,size,sig,es,SimpleNamespace(**rule['hard_gap']),SimpleNamespace(**rule['A_edge']),rule['primary'])
        acrs,genes=H.annotation_positions(chrom,CTX['gene_bed'],CTX['acr_bed'])
        rows=[];bounds=[];candidates=[];iterations=[]
        options=M.options(rule['shape_settings'],rule['minimum_daughter_bp'])
        for parent in baseline['crus']:
            local=[e for e in es if e['left']<parent['end'] and e['right']>=parent['start']]
            rr,bb,cc,it=M.solve(parent,sig,local,size,SimpleNamespace(**rule['primary']),options,acrs,genes,'V11_300kb',rule['max_rounds'])
            rr.sort(key=lambda r:r['start'])
            assert rr[0]['start']==parent['start'] and rr[-1]['end']==parent['end']
            assert all(a['end']==b['start'] for a,b in zip(rr,rr[1:]))
            rows.extend(rr);bounds.extend(bb);candidates.extend(cc);iterations.extend(it)
        rows.sort(key=lambda r:r['start'])
        for i,r in enumerate(rows,1):r['cru_id']=f'{chrom}_V11_300kb_CRU_{i}'
        bounds.extend({**r,'source':'inherited_B145','method':'V11_300kb'} for r in baseline['boundaries'])
        assert sum(r['end']-r['start'] for r in rows)==sum(r['end']-r['start'] for r in baseline['crus'])
        profiles=M.profile_context(sig,es,size,0,size,SimpleNamespace(**rule['primary']))['profiles']
        np.savez_compressed(root/f'f02_inputs/{chrom}.profiles.npz',**{f'{k}_{field}':v[field] for k,v in profiles.items() for field in ['raw','smooth','background']})
        savebed(d/'CRU.bed',rows);savebed(d/'Aedge_B145.bed',baseline['crus'])
        M.write(d/'boundaries.tsv',bounds,['chrom','boundary_pos','source'] if not bounds else None)
        M.write(d/'candidates.tsv.gz',candidates,['chrom','boundary_pos'] if not candidates else None)
        M.write(d/'selection_iterations.tsv',iterations,['chrom','baseline_id'] if not iterations else None)
        M.write(d/'refined_hard_gaps.tsv',aux['refined_gaps'],['chrom','start','end'] if not aux['refined_gaps'] else None)
        stat=M.P.interval_summary('V11_300kb',chrom,rows,es,size)
        payload={'crus':rows,'boundaries':bounds,'stat':stat}
        temp=Path(str(statefile)+'.tmp')
        with temp.open('wb') as f:pickle.dump(payload,f,protocol=4)
        temp.replace(statefile)
        paths=[statefile,d/'CRU.bed',d/'Aedge_B145.bed',d/'boundaries.tsv',d/'candidates.tsv.gz',d/'selection_iterations.tsv',d/'refined_hard_gaps.tsv',root/f'f02_inputs/{chrom}.profiles.npz']
        temp=Path(str(done)+'.tmp');temp.write_text(json.dumps({'signature':CTX['signature'],'artifacts':{str(p):md5(p) for p in paths}},indent=2)+'\n');temp.replace(done)
        print('CALL DONE',chrom,len(rows),'CRUs',flush=True)
    return chrom,'called'


def main(args):
    root=args.outdir.resolve();rules=json.loads(args.rules.read_text())
    sizes={}
    with args.genome_sizes.open() as f:
        for line in f:
            if line.strip():c,n=line.split()[:2];sizes[c]=int(n)
    chroms=args.chroms
    if set(chroms)-set(sizes) or len(set(chroms))!=len(chroms):raise ValueError('Invalid chromosome list')
    if args.workers<1:raise ValueError('--workers must be positive')
    for path,isgene in [(args.gene_bed,True),(args.acr_bed,False)]:
        identities={}
        with path.open() as f:
            for line in f:
                r=line.split()
                if not r or r[0] not in chroms:continue
                if len(r)<4:raise ValueError('BED requires chromosome/start/end/ID')
                c,a,b,k=r[:4];a,b=int(a),int(b)
                if not 0<=a<b<=sizes[c]:raise ValueError('BED out of genome bounds: '+str(path))
                if isgene and (len(r)<6 or r[5] not in ['+','-']):raise ValueError('gene BED must have reliable +/- strand in column 6')
                value=(c,a,b,r[5] if isgene else '.')
                if k in identities and identities[k]!=value:raise ValueError('Conflicting annotation ID: '+k)
                identities[k]=value
    sourcefiles=[Path(__file__).resolve(),args.rules,args.aglink,args.gene_bed,args.acr_bed,args.genome_sizes,
                 Path(M.__file__),Path(H.__file__),Path(H.E.__file__),Path(M.T.__file__),Path(M.A.__file__),Path(M.N.__file__),Path(M.F.__file__),Path(M.C.__file__),Path(M.D.__file__),Path(M.P.__file__)]
    hashes={str(p):md5(p) for p in sourcefiles}
    signature={'inputs':hashes,'chroms':chroms}
    if root.exists() and any(root.iterdir()) and not (root/'00_call_signature.json').exists():raise ValueError('Unknown nonempty calling output directory')
    root.mkdir(parents=True,exist_ok=True);marker=root/'00_call_signature.json'
    if marker.exists():
        if not args.resume or json.loads(marker.read_text())!=signature:raise ValueError('Calling inputs/rules differ; use a new directory')
    else:marker.write_text(json.dumps(signature,indent=2)+'\n')
    for name in ['f01_shared_input','f02_inputs','f03_results','f04_statistics','f06_scripts','chromosomes']:(root/name).mkdir(exist_ok=True)
    loader=SimpleNamespace(aglink=args.aglink,gene_bed=args.gene_bed,acr_bed=args.acr_bed,chroms=chroms,bedpe_reference=None)
    # ID-bearing ACR midpoint / strand-aware gene TSS input, identical to historical calling.
    M.P.load_count_only_edges(loader,sizes,root/'f01_shared_input')
    edges=M.D.load_edges(root/'f01_shared_input/01_count_only_aglink.annotated.bedpe',chroms)
    CTX.update(root=root,rules=rules,sizes=sizes,edges=edges,resume=args.resume,signature=signature,gene_bed=args.gene_bed,acr_bed=args.acr_bed)
    with multiprocessing.get_context('fork').Pool(min(args.workers,len(chroms))) as pool:
        for chrom,state in pool.imap_unordered(worker,chroms):print('CALL',state,chrom,flush=True)
    merged=[];boundaries=[];statistics=[]
    for c in chroms:
        with (root/'chromosomes'/c/'state.pickle').open('rb') as f:state=pickle.load(f)
        merged.extend(state['crus']);boundaries.extend(state['boundaries']);statistics.append(state['stat'])
    final=root/'f03_results/f06_V11_300kb';final.mkdir(exist_ok=True)
    savebed(final/'02_CRU.bed',merged);M.write(final/'03_boundaries.tsv',boundaries,['chrom','boundary_pos','source'] if not boundaries else None)
    lengths=[r['end']-r['start'] for r in merged]
    total={'method':'V11_300kb','chrom':'Tested_chromosomes','cru_n':len(merged),'covered_bp':sum(lengths),
           'coverage_fraction':sum(lengths)/sum(sizes[c] for c in chroms),'median_bp':float(np.median(lengths)) if lengths else np.nan,
           'intra_aglink_n':sum(r['intra_aglink_n'] for r in statistics),'aglink_n':sum(r['aglink_n'] for r in statistics)}
    if lengths:
        total.update(mean_bp=float(np.mean(lengths)),max_bp=max(lengths))
        for label,q in [('p25',.25),('p75',.75),('p90',.9),('p95',.95)]:total[label+'_bp']=float(np.quantile(lengths,q))
    M.write(root/'f04_statistics/08_complete_CRU_summary.tsv',statistics+[total])
    (root/'f04_statistics/07_effective_rules.json').write_text(json.dumps(rules,indent=2)+'\n')
    parameters={'chroms':chroms,'sizes':sizes,'max_rounds':rules['max_rounds'],'aglink_bedpe':str(root/'f01_shared_input/01_count_only_aglink.annotated.bedpe'),
                'gene_bed':str(args.gene_bed),'acr_bed':str(args.acr_bed),'block_bed':str(args.block_bed),'rules':str(args.rules),
                'methods':{'V11_300kb':'f06_V11_300kb'},'expected_intra_aglink_n':total['intra_aglink_n'],'no_coACR_for_calling':True,
                'block_for_plotting_only':True}
    (root/'01_parameters.json').write_text(json.dumps(parameters,indent=2)+'\n')
    assert all(md5(p)==h for p,h in hashes.items())
    M.write(root/'f02_inputs/01_source_MD5.tsv',[{'path':p,'md5_before':h,'md5_after':md5(p)} for p,h in hashes.items()])
    for p in HERE.iterdir():
        if p.is_file():shutil.copy2(p,root/'f06_scripts'/p.name)
    (root/'02_CALLING.DONE').write_text(str(len(merged))+' CRUs; fixed generic rules; only selected chromosomes.\n')
    print('MERGED V11-300kb:',len(merged),'CRUs',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['aglink','gene-bed','acr-bed','genome-sizes','block-bed','outdir','rules']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--chroms',nargs='+',required=True);p.add_argument('--workers',type=int,default=7);p.add_argument('--resume',action='store_true')
    main(p.parse_args())
