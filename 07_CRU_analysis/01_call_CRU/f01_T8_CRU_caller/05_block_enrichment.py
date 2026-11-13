#!/usr/bin/env python3
"""Haplotype-block permutation test with the unchanged V11/M3 local null."""
import argparse
import csv
import hashlib
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
from types import SimpleNamespace as NS


HERE=Path(__file__).resolve().parent
if str(HERE) not in sys.path: sys.path.insert(0, str(HERE))
FROZEN=HERE/'f02_modules'
from f02_modules import m12_block_test as B
from f02_modules import m16_junctions as J


def load(path,name):
    s=importlib.util.spec_from_file_location(name,path)
    m=importlib.util.module_from_spec(s);s.loader.exec_module(m);return m


def md5(path):
    h=hashlib.md5()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda:f.read(1<<20),b''):h.update(chunk)
    return h.hexdigest()


def main(args):
    if args.perm_n<1 or args.local_window_bp<5000:raise ValueError('Require perm-n >=1 and local-window-bp >=5000')
    container=args.outdir.resolve();container.mkdir(parents=True,exist_ok=True)
    root=container/'f01_statistics'
    junction=container/'01_shared_junctions.tsv'
    J.main(NS(cru_bed=args.cru_bed,genome_sizes=args.genome_sizes,chroms=args.chroms,label='CRU',out=junction))
    if args.plot_only and not (root/'07_ENRICHMENT.DONE.json').is_file():
        raise ValueError('Plot-only requires completed block statistics; use --mode block first')
    B.main(NS(cru_bed=args.cru_bed,boundaries=junction,block_bed=args.block_bed,genome_sizes=args.genome_sizes,
              caller=FROZEN/'m17_block_reference_null.py',outdir=root,chroms=args.chroms,
              perm_n=args.perm_n,seed=args.seed,local_window_bp=args.local_window_bp,resume=True))
    if args.plot=='all':
        script=HERE/'07_plot_block.R'
        figures=root/'f01_figures'
        if figures.is_dir() and not (figures/'03_PLOT.DONE.rds').exists():
            backup=root/'f02_failed_plots';backup.mkdir(exist_ok=True);n=1
            while (backup/f'attempt{n:03d}').exists():n+=1
            shutil.move(str(figures),backup/f'attempt{n:03d}')
        env=dict(os.environ);env.pop('R_LIBS_USER',None)
        fonts=load(HERE/'08_font_environment.py','t8_fonts')
        with fonts.plotting_environment(env) as render_env:
            subprocess.run([args.rscript,'--vanilla',str(script),'--analysis-dir',str(root),
                            '--font-family',args.font_family,'--line-scale','1.3'],env=render_env,check=True)
        checker=load(HERE/'09_pdf_qc.py','pdf_qc')
        checks=[checker.check_pdf(p,expected=1,rscript=args.rscript) for p in sorted(figures.glob('*.pdf'))]
        with (root/'08_pdf_page_QC.tsv').open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=['pdf','expected_pages','actual_pages','page_count_status','page_count_backend'],delimiter='\t')
            w.writeheader();w.writerows(checks)
    with (root/'02_enrichment_summary.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    for r in rows:
        if r['scope']=='All_selected' and r['boundary_group']=='all_final':
            print('ALL-FINAL <=%dkb: real %.3f%%; random %.3f%%; fold %.3f; P %.6g; n=%s' %
                  (int(r['cutoff_bp'])//1000,100*float(r['observed_fraction']),100*float(r['random_fraction_mean']),
                   float(r['observed_over_random']),float(r['empirical_p_ge']),r['boundary_n']),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['cru-bed','block-bed','genome-sizes','outdir']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--chroms',nargs='+',required=True)
    p.add_argument('--perm-n',type=int,default=1000);p.add_argument('--seed',type=int,default=20260805)
    p.add_argument('--local-window-bp',type=int,default=10000000)
    p.add_argument('--plot',choices=['none','all'],default='none');p.add_argument('--plot-only',action='store_true')
    p.add_argument('--rscript',default='/home/jzzhang/anaconda3/bin/Rscript')
    p.add_argument('--font-family',default='Liberation Sans')
    main(p.parse_args())
