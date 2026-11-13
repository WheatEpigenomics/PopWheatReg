# m12_block_test: haplotype-block permutation test
# Frozen rule module, copied unchanged from the published calling chain
# (source: f17_V11_pipeline/07_block_enrichment.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Independent outer/internal/all-final block validation of a completed CRU BED."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import random
import runpy
import numpy as np


def md5(path):
    h=hashlib.md5()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()


def write(path,rows):
    with path.open('w') as f:
        w=csv.DictWriter(f,list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)


def distances(points,index,chroms):
    p=np.array([r['pos'] for r in points],dtype=np.int64);distance=np.empty(len(p),dtype=np.int64);targets=np.empty(len(p),dtype=np.int64)
    for c in chroms:
        ii=np.array([i for i,r in enumerate(points) if r['chrom']==c],dtype=int)
        if not len(ii):continue
        bb=np.asarray(index[c]['boundaries'],dtype=np.int64)
        if not len(bb):raise ValueError('Missing block boundaries for '+c)
        j=np.searchsorted(bb,p[ii]);left=bb[np.maximum(j-1,0)];right=bb[np.minimum(j,len(bb)-1)]
        target=np.where(abs(left-p[ii])<=abs(right-p[ii]),left,right)
        targets[ii]=target;distance[ii]=abs(target-p[ii])
    return targets,distance


def main(args):
    inputs=[args.cru_bed,args.boundaries,args.block_bed,args.genome_sizes,args.caller,Path(__file__).resolve()]
    hashes={str(p):md5(p) for p in inputs}
    signature={'inputs':hashes,'chroms':args.chroms,'perm_n':args.perm_n,'seed':args.seed,'local_window_bp':args.local_window_bp}
    root=args.outdir.resolve();root.mkdir(parents=True,exist_ok=True);meta=root/'01_parameters.json'
    if meta.exists():
        if not args.resume or json.loads(meta.read_text())!=signature:raise ValueError('Existing enrichment inputs/options changed')
        if (root/'07_ENRICHMENT.DONE.json').exists():
            saved=json.loads((root/'07_ENRICHMENT.DONE.json').read_text())
            if all(Path(p).exists() and md5(p)==h for p,h in saved.items()):print('REUSE block enrichment',flush=True);return
    elif any(root.iterdir()):raise ValueError('Unknown enrichment output directory')
    else:meta.write_text(json.dumps(signature,indent=2)+'\n')
    helper=runpy.run_path(str(args.caller));sizes=helper['load_genome_sizes'](str(args.genome_sizes))
    starts=set();ends=set()
    with args.cru_bed.open() as f:
        for line in f:
            if not line.strip():continue
            c,a,b=line.split()[:3];a,b=int(a),int(b)
            if c not in args.chroms:continue
            assert 0<=a<b<=sizes[c];starts.add((c,a));ends.add((c,b))
    total={(c,p) for c,p in starts|ends if 0<p<sizes[c]};internal=starts&ends&total
    with args.boundaries.open() as f:saved={(r['chrom'],int(r['boundary_pos'])) for r in csv.DictReader(f,delimiter='\t')}
    assert saved&total==internal,'Saved cut coordinates and shared-CRU endpoints disagree'
    groups={'outer':total-internal,'internal':internal,'all_final':total}
    blocks=helper['load_block_bed'](str(args.block_bed),args.chroms);index=helper['block_test_index'](blocks,args.chroms)
    summary=[];reps=[];details=[];cutoffs=[25000,50000,100000,200000]
    for group,points in groups.items():
        for scope in args.chroms+['All_selected']:
            pp=[{'chrom':c,'pos':p} for c,p in sorted(points) if scope=='All_selected' or c==scope]
            n=len(pp)
            if not n:continue
            nearest,dd=distances(pp,index,args.chroms)
            observed=[int(np.count_nonzero(dd<=k)) for k in cutoffs]
            if scope=='All_selected':details.extend({'group':group,'chrom':r['chrom'],'position':r['pos'],'nearest_block_boundary':int(t),'distance_bp':int(d)} for r,t,d in zip(pp,nearest,dd))
            null=np.zeros((args.perm_n,len(cutoffs)),dtype=int);rng=random.Random(args.seed)
            for rep in range(args.perm_n):
                shuffled=helper['local_random_points'](pp,sizes,args.local_window_bp,rng);assert len(shuffled)==n
                _,dist=distances(shuffled,index,args.chroms)
                row={'boundary_group':group,'scope':scope,'replicate':rep+1,'boundary_n':n}
                for j,k in enumerate(cutoffs):
                    null[rep,j]=np.count_nonzero(dist<=k);row[f'match_{k//1000}kb_n']=int(null[rep,j]);row[f'match_{k//1000}kb_fraction']=null[rep,j]/n
                reps.append(row)
            for j,k in enumerate(cutoffs):
                values=null[:,j];mean=float(values.mean());lo,hi=np.quantile(values,[.025,.975]);obs=observed[j]
                summary.append({'boundary_group':group,'scope':scope,'cutoff_bp':k,'boundary_n':n,'observed_match_n':obs,'observed_fraction':obs/n,
                    'random_match_n_mean':mean,'random_match_n_p025':float(lo),'random_match_n_p975':float(hi),
                    'random_fraction_mean':mean/n,'random_fraction_p025':float(lo/n),'random_fraction_p975':float(hi/n),
                    'observed_over_random':obs/mean if mean else np.nan,'empirical_p_ge':float((1+np.count_nonzero(values>=obs))/(args.perm_n+1)),
                    'permutation_n':args.perm_n})
            if scope=='All_selected':print('BLOCK:',group,'n=',n,'<=25kb fold=',round(summary[-4]['observed_over_random'],3),'P=',summary[-4]['empirical_p_ge'],flush=True)
    write(root/'02_enrichment_summary.tsv',summary);write(root/'03_random_replicates.tsv',reps);write(root/'04_boundary_distances.tsv',details)
    assert all(md5(p)==h for p,h in hashes.items())
    write(root/'05_input_MD5.tsv',[{'path':p,'md5_before':h,'md5_after':md5(p)} for p,h in hashes.items()])
    (root/'06_METHODS.md').write_text('Independent validation only: no calling/score changes. Outer=nonshared endpoints; internal=shared CRU end/start verified against saved cuts; all-final=deduplicated union. Physical chromosome 0/size endpoints excluded from CRU queries. Block start/end coordinates use the original M3 helper, deduplicated per chromosome. Original +/-10-Mb local point null excludes distances <5 kb from each source and clips to 1..size-1. Seed resets per group/scope; random coincidences are not deduplicated so query N is preserved. 95% intervals are empirical null percentiles, not a confidence interval on the mean. P=(1+#null>=observed)/(B+1). All_selected includes exactly the requested chromosomes.\n')
    outputs=[root/name for name in ['02_enrichment_summary.tsv','03_random_replicates.tsv','04_boundary_distances.tsv','05_input_MD5.tsv','06_METHODS.md']]
    (root/'07_ENRICHMENT.DONE.json').write_text(json.dumps({str(p):md5(p) for p in outputs},indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for k in ['cru-bed','boundaries','block-bed','genome-sizes','caller','outdir']:p.add_argument('--'+k,type=Path,required=True)
    p.add_argument('--chroms',nargs='+',required=True);p.add_argument('--perm-n',type=int,default=1000)
    p.add_argument('--seed',type=int,default=20260805);p.add_argument('--local-window-bp',type=int,default=10000000)
    p.add_argument('--resume',action='store_true');main(p.parse_args())
