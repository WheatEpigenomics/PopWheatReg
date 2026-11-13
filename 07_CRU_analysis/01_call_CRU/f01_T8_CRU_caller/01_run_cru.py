#!/usr/bin/env python3
"""CRU calling: one chromosome per process, optional 5-Mb atlas and block test.

Plotting never blocks calling: a failed R atlas is recorded in
07_plot_failures.tsv while the calls are still merged and reported.
Redraw later with --mode plot --plot all --resume.
"""
import argparse
from collections import Counter
from concurrent.futures import ThreadPoolExecutor, as_completed
import contextlib
import csv
from datetime import datetime
import fcntl
import gc
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
from types import SimpleNamespace as NS


HERE=Path(__file__).resolve().parent
STANDARD=[f'{i}{s}' for i in range(1,8) for s in 'ABD']


def md5(path):
    h=hashlib.md5()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(1<<20),b''):h.update(block)
    return h.hexdigest()


def atomic_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.tmp')
    temp.write_text(json.dumps(data,indent=2,sort_keys=True)+'\n');temp.replace(path)


def write_table(path,rows,fields=None):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]),delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)


def load_core():
    s=importlib.util.spec_from_file_location('cru_caller',HERE/'02_call_chain.py')
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def font_module():
    s=importlib.util.spec_from_file_location('t8_fonts',HERE/'08_font_environment.py')
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def inventory(base,paths):
    return {str(Path(p).relative_to(base)):md5(p) for p in paths}


def complete(marker,signature,base):
    if not marker.is_file():return False
    try:
        old=json.loads(marker.read_text())
        return old['signature']==signature and bool(old['artifacts']) and all(
            (base/p).is_file() and md5(base/p)==h for p,h in old['artifacts'].items())
    except (ValueError,KeyError,OSError):return False


def read_sizes(path):
    sizes={}
    for n,line in enumerate(Path(path).read_text().splitlines(),1):
        if not line.strip() or line.startswith('#'):continue
        r=line.split()
        if len(r)<2 or r[0] in sizes or int(r[1])<=0:raise ValueError(f'Invalid genome size at line {n}')
        sizes[r[0]]=int(r[1])
    return sizes


def partition_inputs(args,chroms,sizes,root,signature):
    """Stream the raw gene/ACR-pair table once; workers never open other chromosomes."""
    base=root/'f01_inputs';done=base/'01_INPUTS.DONE.json'
    if complete(done,signature,root):return
    base.mkdir(exist_ok=True)
    metadata=[];artifacts=[]
    for c in chroms:(base/c).mkdir(exist_ok=True)
    with contextlib.ExitStack() as stack:
        for source,filename in [(args.gene_bed,'01_genes.bed'),(args.acr_bed,'02_ACRs.bed')]:
            handles={c:stack.enter_context((base/c/filename).open('w')) for c in chroms};meta={}
            with source.open() as f:
                for n,line in enumerate(f,1):
                    if not line.strip() or line.startswith('#'):continue
                    row=line.split()
                    if len(row)<4 or row[0] not in handles:continue
                    a,b=int(row[1]),int(row[2]);c=row[0]
                    if not 0<=a<b<=sizes[c]:raise ValueError(f'Invalid annotation coordinate {source}:{n}')
                    if row[3] in meta and meta[row[3]]!=c:raise ValueError(f'Ambiguous cross-chromosome annotation ID {row[3]}')
                    meta[row[3]]=c;handles[c].write(line if line.endswith('\n') else line+'\n')
            metadata.append(meta)
            artifacts.extend(base/c/filename for c in chroms)
        handles={c:stack.enter_context((base/c/'03_aglink.tsv').open('w',newline='')) for c in chroms}
        writers={c:csv.writer(f,delimiter='\t',lineterminator='\n') for c,f in handles.items()}
        maps={c:stack.enter_context((base/c/'04_source_lines.tsv').open('w')) for c in chroms}
        audit=Counter();genes,acrs=metadata
        with args.aglink.open(newline='') as f:
            for n,row in enumerate(csv.reader(f,delimiter='\t'),1):
                if not row or row[0].startswith('#'):continue
                if len(row)<2:raise ValueError(f'agLINK row {n} has fewer than two columns')
                audit['input_rows']+=1;g,a=row[:2]
                if g not in genes or a not in acrs:audit['not_both_in_selected_metadata']+=1;continue
                if genes[g]!=acrs[a]:audit['selected_interchromosomal_rows']+=1;continue
                c=genes[g];writers[c].writerow(row);maps[c].write(str(n)+'\n');audit[c+'_retained']+=1
        artifacts.extend(base/c/name for c in chroms for name in ['03_aglink.tsv','04_source_lines.tsv'])
        for c in chroms:
            if not audit[c+'_retained']:raise ValueError(f'No agLINKs on {c}; cannot run frozen density quantiles')
        for source,filename in [(args.block_bed,'05_blocks.bed'),(args.reference_bed,'07_reference.bed')]:
            if source is None:continue
            outs={c:stack.enter_context((base/c/filename).open('w')) for c in chroms}
            with source.open() as f:
                for line in f:
                    if not line.strip() or line.startswith('#'):continue
                    row=line.split()
                    if row[0] not in outs:continue
                    if len(row)<3 or not 0<=int(row[1])<int(row[2])<=sizes[row[0]]:
                        raise ValueError(f'Invalid BED row in {source}')
                    outs[row[0]].write(line if line.endswith('\n') else line+'\n')
            artifacts.extend(base/c/filename for c in chroms)
    for c in chroms:
        path=base/c/'06_chromosome.sizes';path.write_text(f'{c}\t{sizes[c]}\n');artifacts.append(path)
    path=base/'02_partition_audit.tsv'
    write_table(path,[dict(metric=k,value=v) for k,v in sorted(audit.items())]);artifacts.append(path)
    atomic_json(done,dict(signature=signature,artifacts=inventory(root,artifacts)))


def partition_coacr(args,chroms,root):
    """Plot-only ACR midpoint pairs; separate signature cannot alter calling."""
    if args.coacr is None:return
    signature=dict(coacr_md5=md5(args.coacr),acr_md5=md5(args.acr_bed),chroms=chroms,
                   min_abs_pcc=args.coacr_min_abs_pcc,max_distance_bp=args.coacr_max_distance_bp)
    base=root/'f01_inputs';done=base/'03_COACR.DONE.json'
    if complete(done,signature,root):return
    meta={}
    with args.acr_bed.open() as f:
        for line in f:
            r=line.split()
            if not r or r[0] not in chroms or len(r)<4:continue
            value=(r[0],(int(r[1])+int(r[2]))//2)
            if r[3] in meta and meta[r[3]]!=value:raise ValueError('Conflicting ACR ID in coACR display metadata: '+r[3])
            meta[r[3]]=value
    pairs={c:{} for c in chroms};audit=Counter()
    with args.coacr.open() as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not {'id1','id2','correlation'}<=set(reader.fieldnames or []):raise ValueError('coACR requires id1,id2,correlation columns')
        for r in reader:
            audit['input_rows']+=1;a,b=r['id1'],r['id2']
            if a==b or a not in meta or b not in meta:continue
            c,p=meta[a];d,q=meta[b];value=float(r['correlation'])
            if c!=d or not math.isfinite(value) or abs(value)<args.coacr_min_abs_pcc or abs(p-q)>args.coacr_max_distance_bp:continue
            key=tuple(sorted((a,b)))
            if key not in pairs[c] or abs(value)>abs(pairs[c][key]):pairs[c][key]=value
    artifacts=[]
    for c in chroms:
        path=base/c/'08_coacr.tsv';rows=[]
        for (a,b),value in pairs[c].items():
            _,p=meta[a];_,q=meta[b]
            rows.append(dict(chrom=c,left=min(p,q),right=max(p,q),id1=a,id2=b,correlation=value))
        write_table(path,rows,fields=['chrom','left','right','id1','id2','correlation'])
        artifacts.append(path);audit[c+'_unique_display_pairs']=len(rows)
    path=base/'04_coacr_display_audit.tsv'
    write_table(path,[dict(metric=k,value=v) for k,v in sorted(audit.items())]);artifacts.append(path)
    atomic_json(done,dict(signature=signature,artifacts=inventory(root,artifacts)))


def archive_partial(work):
    paths=[p for p in work.iterdir() if p.name!='f06_failed_attempts'] if work.exists() else []
    if not paths:return
    archive=work/'f06_failed_attempts';archive.mkdir(exist_ok=True)
    n=1
    while (archive/f'attempt{n:03d}').exists():n+=1
    dest=archive/f'attempt{n:03d}';dest.mkdir()
    for p in paths:shutil.move(str(p),dest/p.name)


def bed_rows(path,chrom,size):
    rows=[];ids=set()
    for line in path.read_text().splitlines():
        r=line.split();a,b=int(r[1]),int(r[2])
        if r[0]!=chrom or not 0<=a<b<=size or r[3] in ids:raise ValueError(f'Invalid or duplicated CRU in {path}')
        rows.append(dict(chrom=chrom,start=a,end=b,cru_id=r[3]));ids.add(r[3])
    if rows!=sorted(rows,key=lambda r:(r['start'],r['end'])):raise ValueError(f'Unsorted BED: {path}')
    if any(a['end']>b['start'] for a,b in zip(rows,rows[1:])):raise ValueError(f'Overlapping CRUs: {path}')
    return rows


def perform_call(plan,chrom,work):
    import numpy as np
    C=load_core();source=Path(plan['outdir'])/'f01_inputs'/chrom;size=plan['sizes'][chrom]
    raw=work/'f01_raw_input';out=work/'f02_calling';qc=work/'f03_diagnostics';plot=work/'f04_plot_data'
    for d in [raw,out,qc,plot]:d.mkdir(parents=True,exist_ok=True)
    inputs=NS(aglink=source/'03_aglink.tsv',gene_bed=source/'01_genes.bed',acr_bed=source/'02_ACRs.bed',
              chroms=[chrom],bedpe_reference=None)
    C.M.P.load_count_only_edges(inputs,{chrom:size},raw)
    # Restore original table line IDs after partitioning; equal-position ordering is unchanged.
    original=[int(x) for x in (source/'04_source_lines.tsv').read_text().splitlines()]
    path=raw/'01_count_only_aglink.annotated.bedpe'
    with path.open() as f:annotated=list(csv.DictReader(f,delimiter='\t'))
    for row in annotated:
        line=original[int(row['source_line'])-1];row['source_line']=line;row['edge_id']=f'aglink_input_{line}'
    C.P.write(path,annotated)
    edges=C.M.D.load_edges(path,[chrom])[chrom]
    acrs,genes=C.H.annotation_positions(chrom,inputs.gene_bed,inputs.acr_bed)
    def bodies(p):return [(int(r[1]),int(r[2])) for r in (x.split() for x in p.read_text().splitlines()) if r and r[0]==chrom]
    acr_iv,gene_iv=bodies(inputs.acr_bed),bodies(inputs.gene_bed)
    config=json.loads(Path(plan['config']).read_text())
    parameters=[config[k] for k in ['baseline','low_density_recall','support_recall','island_rescue']]
    common=(chrom,size,edges,*parameters,acrs,genes,acr_iv,gene_iv)
    started=time.monotonic();control=None
    if plan['control']=='replay':
        print('Replaying unmodified control from raw chromosome inputs',flush=True)
        control=C.run_chain(*common,False)
        if plan['reference_bed']:
            historical=C.P.read_bed(source/'07_reference.bed')
            coords=lambda rr:[(int(r['start']),int(r['end'])) for r in sorted(rr,key=lambda r:r['start'])]
            if coords(control['crus'])!=coords(historical):raise ValueError('Raw control differs from supplied historical reference')
        C.P.save_bed(out/f'{chrom}.T7A_control.CRU.bed',control['crus'])
    print('CRU calling:',chrom,flush=True)
    result=C.run_chain(*common,True)
    if control is not None:
        for name in ['parents','refined_gaps']:
            assert [(r['start'],r['end']) for r in result[name]]==[(r['start'],r['end']) for r in control[name]]
    assert C.P.disjoint(result['crus'])
    for r in result['crus']:
        assert not result['protected'][r['start']//5000:(r['end']+4999)//5000].any()
    assert all(abs(r.get('shift_bp',0))<=10000 for r in result['audit']+result['outer_audit'] if r.get('shift_bp') is not None)
    for i,r in enumerate(result['crus'],1):r['cru_id']=f'{chrom}_CRU_{i}'
    C.P.save_bed(out/f'{chrom}.CRU.bed',result['crus'])
    tables=[('01_fine_evaluations.tsv.gz','audit'),('02_final_internal_cuts.tsv','boundaries'),
            ('03_outer_boundary_sources.tsv','origins'),('04_candidate_evaluations.tsv.gz','candidates'),
            ('05_iterations.tsv','iterations'),('06_recall_decisions.tsv.gz','recall'),
            ('07_preserved_hard_gaps.tsv','refined_gaps'),('08_outer_birth_evaluations.tsv','outer_audit')]
    for filename,key in tables:C.P.write(qc/filename,result[key])
    def stats(method,data):
        lengths=[r['end']-r['start'] for r in data['crus']]
        points={p for r in data['crus'] for p in [r['start'],r['end']]}
        return dict(chrom=chrom,method=method,cru_n=len(lengths),chromosome_bp=size,covered_bp=sum(lengths),
            coverage_fraction=sum(lengths)/size,median_bp=float(np.median(lengths)) if lengths else None,
            mean_bp=float(np.mean(lengths)) if lengths else None,p90_bp=float(np.percentile(lengths,90)) if lengths else None,
            max_bp=max(lengths,default=0),ge5Mb_n=sum(x>=5000000 for x in lengths),ge10Mb_n=sum(x>=10000000 for x in lengths),
            internal_cut_n=len(data['boundaries']),unique_boundary_n=len(points),
            annotation_split_boundary_n=0)
    summary=[stats('CRU',result)]
    if control is not None:summary.insert(0,stats('T7_A',control))
    ai,gi=C.O.IntervalIndex(acr_iv),C.O.IntervalIndex(gene_iv)
    for row,data in [(summary[-1],result)]+([(summary[0],control)] if control is not None else []):
        points={p for r in data['crus'] for p in [r['start'],r['end']]}
        row['annotation_split_boundary_n']=sum(bool(ai.split_body_n(p) or gi.split_body_n(p)) for p in points)
    C.P.write(qc/'09_summary.tsv',summary)
    positions=np.arange(len(result['sig']['raw']))*5000+2500;sig=result['sig']
    C.P.write(plot/'01_density.tsv.gz',(
        dict(position=int(p),raw=float(sig['raw'][i]),density25=float(sig['smooth'][0,i]),
             density105=float(sig['smooth'][1,i]),density225=float(sig['smooth'][2,i]),background=float(sig['background'][i]))
        for i,p in enumerate(positions)),fields=['position','raw','density25','density105','density225','background'])
    C.P.write(plot/'02_aglink.tsv.gz',(
        dict(left=e['left'],right=e['right'],acr_id=e['acr_id'],gene_id=e['gene_id']) for e in edges),
        fields=['left','right','acr_id','gene_id'])
    C.P.write(plot/'03_ACR_bodies.tsv',[dict(start=a,end=b) for a,b in acr_iv],fields=['start','end'])
    C.P.write(plot/'04_gene_bodies.tsv',[dict(start=a,end=b) for a,b in gene_iv],fields=['start','end'])
    # Checkpoint validation is based on files, not the worker exit code alone.
    bed_rows(out/f'{chrom}.CRU.bed',chrom,size)
    C.P.write(qc/'10_stage_counts.tsv',[dict(stage=k,accepted_candidate_n=v) for k,v in sorted(Counter(
        r['method'] for r in result['recall'] if r.get('accepted')).items())],fields=['stage','accepted_candidate_n'])
    print('CALL DONE',chrom,summary,'elapsed_seconds',round(time.monotonic()-started,2),flush=True)


def run_worker(plan,chrom):
    root=Path(plan['outdir']);work=root/'f02_chromosomes'/chrom;work.mkdir(parents=True,exist_ok=True)
    done=work/'01_CALL.DONE.json';signature=dict(run=plan['signature'],chrom=chrom)
    valid=complete(done,signature,work)
    if plan['mode'] in ['plot','block'] and not valid:raise ValueError('This mode requires a validated calling checkpoint: '+chrom)
    if not valid:
        archive_partial(work)
        with (root/'f03_logs'/f'{chrom}.call.log').open('a',buffering=1) as log:
            with contextlib.redirect_stdout(log),contextlib.redirect_stderr(log):
                print('\nSTART',datetime.now().isoformat(),chrom,flush=True);perform_call(plan,chrom,work)
        artifacts=[p for d in ['f01_raw_input','f02_calling','f03_diagnostics','f04_plot_data'] for p in (work/d).rglob('*') if p.is_file()]
        atomic_json(done,dict(signature=signature,artifacts=inventory(work,artifacts)))
        gc.collect()
    if plan['plot']=='all' and plan['mode'] in ['all','plot']:
        opts={k:plan[k] for k in ['window_bp','font_family','dpi','max_aglink_arcs','max_coacr_arcs','seed']}
        opts['coacr_md5']=md5(root/'f01_inputs'/chrom/'08_coacr.tsv')
        opts['font_bundle']=font_module().signature()
        psig=dict(call_checkpoint=md5(done),options=opts,script_md5=md5(HERE/'06_plot_windows.R'),
                  pdf_checker_md5=md5(HERE/'09_pdf_qc.py'))
        pdone=work/'02_PLOT.DONE.json'
        if not complete(pdone,psig,work):
            cmd=[plan['rscript'],'--vanilla',str(HERE/'06_plot_windows.R'),'--chrom-dir',str(work),
                 '--input-dir',str(root/'f01_inputs'/chrom),'--chrom',chrom,'--size',str(plan['sizes'][chrom]),
                 '--window-bp',str(plan['window_bp']),'--font-family',plan['font_family'],'--dpi',str(plan['dpi']),
                 '--max-aglink-arcs',str(plan['max_aglink_arcs']),'--max-coacr-arcs',str(plan['max_coacr_arcs']),
                 '--seed',str(plan['seed'])]
            env=dict(os.environ,CRU_PYTHON=sys.executable);env.pop('R_LIBS_USER',None)
            failed=work/'03_PLOT.FAILED.json'
            try:
                with (root/'f03_logs'/f'{chrom}.plot.log').open('a') as log:
                    with font_module().plotting_environment(env) as render_env:
                        subprocess.run(cmd,env=render_env,stdout=log,stderr=subprocess.STDOUT,check=True)
                artifacts=[p for p in (work/'f05_figures').rglob('*') if p.is_file()]
                atomic_json(pdone,dict(signature=psig,artifacts=inventory(work,artifacts)))
                if failed.exists():failed.unlink()
            except Exception as exc:
                # Figures are optional: the validated calls are already on disk.
                atomic_json(failed,dict(signature=psig,error=f'{type(exc).__name__}: {exc}',
                    note='calling is intact; redraw with --mode plot --plot all --resume'))
                print('PLOT FAILED',chrom,'- calling is intact; redraw later with --mode plot --plot all --resume:',
                      exc,flush=True)
    print('CHROMOSOME DONE',chrom,flush=True)


def merge_results(plan):
    import numpy as np
    root=Path(plan['outdir']);out=root/'f04_merged';out.mkdir(exist_ok=True)
    allrows=[];summaries=[];indexes=[]
    for c in plan['chroms']:
        work=root/'f02_chromosomes'/c
        if not complete(work/'01_CALL.DONE.json',dict(run=plan['signature'],chrom=c),work):raise ValueError('Missing or invalid call: '+c)
        rows=bed_rows(work/'f02_calling'/f'{c}.CRU.bed',c,plan['sizes'][c]);allrows.extend(rows)
        with (work/'f03_diagnostics/09_summary.tsv').open() as f:
            s=next(r for r in csv.DictReader(f,delimiter='\t') if r['method']=='CRU')
        assert int(s['cru_n'])==len(rows) and int(s['covered_bp'])==sum(r['end']-r['start'] for r in rows)
        summaries.append(s)
        path=work/'f05_figures/01_window_index.tsv'
        if plan['plot']=='all' and plan['mode'] in ['all','plot'] and path.is_file():
            with path.open() as f:indices=list(csv.DictReader(f,delimiter='\t'))
            assert len(indices)==(plan['sizes'][c]+plan['window_bp']-1)//plan['window_bp']
            indexes.extend(indices)
    path=out/'CRU.bed'
    with path.open('w') as f:
        for r in allrows:f.write(f"{r['chrom']}\t{r['start']}\t{r['end']}\t{r['cru_id']}\n")
    lengths=[r['end']-r['start'] for r in allrows];den=sum(plan['sizes'][c] for c in plan['chroms'])
    total=dict(summaries[0],chrom='Genomewide' if len(plan['chroms'])==21 else 'Selected_chromosomes',
               cru_n=len(allrows),chromosome_bp=den,covered_bp=sum(lengths),coverage_fraction=sum(lengths)/den,
               median_bp=float(np.median(lengths)) if lengths else None,mean_bp=float(np.mean(lengths)) if lengths else None,
               p90_bp=float(np.percentile(lengths,90)) if lengths else None,max_bp=max(lengths,default=0))
    for k in ['ge5Mb_n','ge10Mb_n','internal_cut_n','unique_boundary_n','annotation_split_boundary_n']:
        total[k]=sum(int(s[k]) for s in summaries)
    write_table(out/'01_chromosome_summary.tsv',summaries+[total])
    if indexes:write_table(out/'02_all_window_index.tsv',indexes)
    print('SUMMARY',json.dumps(total),flush=True)
    print('BED:',path,'\nFigures:',root/'f02_chromosomes/<chrom>/f05_figures',flush=True)


def main(args):
    if args.worker_plan:
        run_worker(json.loads(args.worker_plan.read_text()),args.worker_chrom);return
    if not args.outdir:raise ValueError('--outdir is required')
    for k in ['aglink','gene_bed','acr_bed','genome_sizes','config','block_bed','reference_bed','coacr']:
        p=getattr(args,k)
        if p is not None:
            p=p.resolve();setattr(args,k,p)
            if not p.is_file():raise ValueError('Missing input --'+k.replace('_','-')+': '+str(p))
        elif k not in ['block_bed','reference_bed','coacr']:raise ValueError('--'+k.replace('_','-')+' is required')
    if args.reference_bed and args.control!='replay':raise ValueError('--reference-bed requires --control replay')
    cfg=json.loads(args.config.read_text())
    if not all(isinstance(cfg.get(k),dict) for k in ['baseline','low_density_recall','support_recall','island_rescue']):
        raise ValueError('config requires baseline, low_density_recall, support_recall, island_rescue objects')
    if args.workers<1 or args.window_bp<1 or args.dpi<1 or args.max_aglink_arcs<0:raise ValueError('Invalid worker/render setting')
    if args.max_coacr_arcs<0 or args.coacr_min_abs_pcc<0 or args.coacr_max_distance_bp<0:raise ValueError('Invalid coACR display setting')
    if args.mode=='plot' and args.plot!='all':raise ValueError('--mode plot requires --plot all')
    if args.mode=='block' and args.block=='none':raise ValueError('--mode block requires --block all')
    if args.block=='all' and (args.block_perm_n<1 or args.local_window_bp<5000):raise ValueError('Invalid block-permutation settings')
    sizes=read_sizes(args.genome_sizes)
    chroms=[c for c in STANDARD if c in sizes] if args.chroms==['all'] else args.chroms
    if not chroms or len(chroms)!=len(set(chroms)) or any(c not in sizes for c in chroms):raise ValueError('Invalid chromosome list')
    chroms=sorted(chroms,key=lambda c:(STANDARD.index(c) if c in STANDARD else len(STANDARD),c))
    if args.block=='all' and args.block_bed is None:raise ValueError('--block-bed is required for block validation')
    if args.plot=='all' and args.mode!='call':
        if args.block_bed is None:raise ValueError('--block-bed is required for the block reference track')
        if args.mode!='block' and args.coacr is None:raise ValueError('--coacr is required for the coACR display track')
        if shutil.which(args.rscript) is None:raise ValueError('Rscript is not available')
        print('Plot font:',font_module().validate_family(args.font_family),flush=True)
    root=args.outdir.resolve()
    if root==HERE or HERE in root.parents:raise ValueError('Keep output outside the source-code folder')
    existing=root.exists() and any(root.iterdir())
    if existing and not args.resume:raise ValueError('Nonempty output; use --resume or a new output directory')
    if existing and not (root/'01_run_signature.json').is_file():raise ValueError('Unmanaged output directory; refusing to change it')
    root.mkdir(parents=True,exist_ok=True)
    with (root/'05_run.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        sources=[HERE/n for n in ['01_run_cru.py','02_call_chain.py','03_fine_boundary.py',
                                  '04_boundary_birth.py','05_block_enrichment.py']]
        sources+=sorted((HERE/'f02_modules').rglob('*.py'))
        input_paths=[getattr(args,k) for k in ['aglink','gene_bed','acr_bed','genome_sizes','config','block_bed','reference_bed'] if getattr(args,k) is not None]
        # Block data never enter calling, so the block BED cannot invalidate a call.
        calling_inputs=[getattr(args,k) for k in ['aglink','gene_bed','acr_bed','genome_sizes','config','reference_bed'] if getattr(args,k) is not None]
        signature=dict(chroms=chroms,control=args.control,inputs={str(p):md5(p) for p in calling_inputs},
                       code={str(p.relative_to(HERE)):md5(p) for p in sources})
        marker=root/'01_run_signature.json'
        if marker.exists() and json.loads(marker.read_text())!=signature:raise ValueError('Calling inputs/code/settings changed; use a new output directory')
        atomic_json(marker,signature)
        for d in ['f01_inputs','f02_chromosomes','f03_logs','f04_merged']:(root/d).mkdir(exist_ok=True)
        if not (root/'f05_scripts').exists():
            shutil.copytree(HERE,root/'f05_scripts',ignore=shutil.ignore_patterns('__pycache__','*.pyc'))
        write_table(root/'02_input_MD5.tsv',[dict(path=str(p),md5=md5(p)) for p in input_paths])
        partition_inputs(args,chroms,sizes,root,signature)
        if args.plot=='all' and args.mode in ['all','plot']:partition_coacr(args,chroms,root)
        plan={k:str(v) if isinstance(v,Path) else v for k,v in vars(args).items()}
        plan.update(outdir=str(root),signature=signature,chroms=chroms,sizes={c:sizes[c] for c in chroms})
        atomic_json(root/'06_execution_plan.json',plan)
        if (root/'04_RUN.DONE.json').exists():(root/'04_RUN.DONE.json').unlink()
        status={c:dict(chrom=c,status='pending',error='') for c in chroms}
        def dispatch(c):
            env=dict(os.environ,PYTHONDONTWRITEBYTECODE='1',OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1',NUMEXPR_NUM_THREADS='1')
            with (root/'f03_logs'/f'{c}.worker.log').open('a') as f:
                subprocess.run([sys.executable,'-B',str(HERE/'01_run_cru.py'),'--worker-plan',str(root/'06_execution_plan.json'),
                                '--worker-chrom',c],env=env,stdout=f,stderr=subprocess.STDOUT,check=True)
        print('Chromosomes:',','.join(chroms),'workers:',min(args.workers,len(chroms)),'plot:',args.plot,'control:',args.control,flush=True)
        with ThreadPoolExecutor(max_workers=min(args.workers,len(chroms))) as pool:
            pending={pool.submit(dispatch,c):c for c in chroms}
            for future in as_completed(pending):
                c=pending[future]
                try:future.result();status[c]['status']='done';print('DONE',c,flush=True)
                except Exception as exc:status[c].update(status='failed',error=str(exc));print('FAILED',c,'see',root/'f03_logs'/f'{c}.worker.log',flush=True)
                write_table(root/'03_status.tsv',[status[c] for c in chroms])
        if any(s['status']!='done' for s in status.values()):raise RuntimeError('Some chromosomes failed; valid checkpoints are preserved. Repeat with --resume.')
        merge_results(plan)
        plot_failures=[]
        if args.block=='all':
            block_out=args.block_outdir.resolve() if args.block_outdir else root/'f06_block_enrichment'
            cmd=[sys.executable,'-B',str(HERE/'05_block_enrichment.py'),
                 '--cru-bed',str(root/'f04_merged/CRU.bed'),'--block-bed',str(args.block_bed),
                 '--genome-sizes',str(args.genome_sizes),'--outdir',str(block_out),
                 '--chroms',*chroms,'--perm-n',str(args.block_perm_n),'--seed',str(args.seed),
                 '--local-window-bp',str(args.local_window_bp),'--rscript',args.rscript,
                 '--font-family',args.font_family,'--plot','none']
            env=dict(os.environ,CRU_PYTHON=sys.executable);env.pop('R_LIBS_USER',None)
            with (root/'f03_logs/block_enrichment.log').open('a') as log:
                subprocess.run(cmd,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
                if args.plot=='all' and args.mode!='call':
                    try:
                        subprocess.run(cmd[:-2]+['--plot','all','--plot-only'],env=env,
                                       stdout=log,stderr=subprocess.STDOUT,check=True)
                    except subprocess.CalledProcessError as exc:
                        plot_failures.append(('block_figures',str(exc)))
                        print('PLOT FAILED block figures - block statistics are complete;',flush=True)
            print('Block summary:',block_out/'f01_statistics/02_enrichment_summary.tsv',flush=True)
        for c in chroms:
            record=root/'f02_chromosomes'/c/'03_PLOT.FAILED.json'
            if record.is_file():
                try: plot_failures.append((c,json.loads(record.read_text()).get('error','')))
                except ValueError: plot_failures.append((c,'unreadable plot failure record'))
        if plot_failures:
            write_table(root/'07_plot_failures.tsv',[dict(target=x,error=e) for x,e in plot_failures])
            print('PLOTS: %d failed, %d ok. Calling, merging and block statistics are complete.'
                  %(len(plot_failures),len(chroms)-sum(1 for x,_ in plot_failures if x in chroms)),flush=True)
            print('Redraw chromosome figures later with: --mode plot --plot all --resume',flush=True)
        atomic_json(root/'04_RUN.DONE.json',dict(signature=signature,chroms=chroms,plot=args.plot,block=args.block,
                    artifacts=inventory(root,[p for p in (root/'f04_merged').iterdir() if p.is_file()])))


def parser():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['aglink','gene-bed','acr-bed','genome-sizes','block-bed','reference-bed','coacr','outdir']:
        p.add_argument('--'+k,type=Path)
    p.add_argument('--config',type=Path,default=HERE/'f01_rules/CRU_rules.json')
    p.add_argument('--chroms',nargs='+',default=['all'])
    p.add_argument('--workers',type=int,default=7)
    p.add_argument('--control',choices=['none','replay'],default='none',help='optional raw-input historical replay; not required')
    p.add_argument('--plot',choices=['none','all'],default='none')
    p.add_argument('--mode',choices=['all','call','plot','block'],default='all')
    p.add_argument('--block',choices=['none','all'],default='none',help='block permutation test (needs --block-bed)')
    p.add_argument('--block-perm-n',type=int,default=1000)
    p.add_argument('--local-window-bp',type=int,default=10000000)
    p.add_argument('--block-outdir',type=Path)
    p.add_argument('--resume',action='store_true')
    p.add_argument('--rscript',default='Rscript')
    p.add_argument('--font-family',default='Liberation Sans')
    p.add_argument('--window-bp',type=int,default=5000000)
    p.add_argument('--dpi',type=int,default=120)
    p.add_argument('--max-aglink-arcs',type=int,default=1200)
    p.add_argument('--max-coacr-arcs',type=int,default=1200)
    p.add_argument('--coacr-min-abs-pcc',type=float,default=.6)
    p.add_argument('--coacr-max-distance-bp',type=int,default=1000000)
    p.add_argument('--seed',type=int,default=20260805)
    p.add_argument('--worker-plan',type=Path,help=argparse.SUPPRESS)
    p.add_argument('--worker-chrom',help=argparse.SUPPRESS)
    return p


if __name__=='__main__':main(parser().parse_args())
