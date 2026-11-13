#!/usr/bin/env python3
"""Neutral Wright--Fisher simulation using the original prepared biological inputs.

Scientific functions are consolidated from the source project. Input staging,
CLI orchestration and chromosome-count-dependent labels are adapted for this package.
"""
from __future__ import annotations
import argparse, bisect, configparser, contextlib, csv, datetime, fcntl, gzip
import hashlib, json, math, os, random, re, shlex, shutil, socket, subprocess, sys
import tempfile, time
from pathlib import Path
from typing import Dict, List
from collections import defaultdict, Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
try:
    import numpy as np
except ImportError:
    raise SystemExit('NumPy is required. Install numpy or set PYTHON_BIN to a Python with NumPy.')

CHROMOSOMES = [f'{number}{genome}' for number in range(1, 8) for genome in 'ABD']
PLOTTED_GENERATIONS = [20, 50, 100, 200]
DISTANCE_BINS = [0, 10000, 20000, 50000, 100000, 250000, 500000, 750000, 1000000]
FEATURE_NAMES = ('gene_density', 'acr_density', 'te_fraction', 'recomb_rank')


# Original functions: common.py

def load_config(path):
    cp = configparser.ConfigParser(interpolation=configparser.BasicInterpolation())
    with open(path) as fh:
        cp.read_file(fh)
    return cp


def csv_list(s: str) -> List[str]:
    return [x.strip() for x in s.split(',') if x.strip()]


def int_list(s: str) -> List[int]:
    return [int(x.strip()) for x in s.split(',') if x.strip()]


def open_text(path):
    path = str(path)
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path, 'rt')


def norm_chr(x: str) -> str:
    x = x.strip()
    return x[3:] if x.lower().startswith('chr') else x


def parse_gff_attrs(s: str) -> Dict[str,str]:
    out = {}
    for field in s.split(';'):
        if '=' in field:
            k,v = field.split('=',1)
            out[k] = v
    return out


def interval_overlap(a0:int,a1:int,b0:int,b1:int)->bool:
    return a0 < b1 and b0 < a1


def stable_seed(base:int, *parts) -> int:
    h = hashlib.sha256((str(base)+'|'+'|'.join(map(str,parts))).encode()).digest()
    return int.from_bytes(h[:8],'little') & 0x7fffffffffffffff


def ensure_dir(p):
    Path(p).mkdir(parents=True, exist_ok=True)


def read_blocks_bed(path, chrom=None, max_block_bp=None, coordinate_mode='bed0'):
    if coordinate_mode not in ('bed0','one_based_inclusive'):
        raise ValueError(f'Unsupported block coordinate mode: {coordinate_mode}')
    rows=[]; large=[]
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'): continue
            f=line.rstrip('\n').split('\t')
            if len(f)<3: continue
            c=norm_chr(f[0]); s=int(f[1]); e=int(f[2])
            if coordinate_mode == 'one_based_inclusive':
                s -= 1
            bid=f[3] if len(f)>3 else f'{c}:{s}-{e}'
            if chrom and c!=chrom: continue
            if s < 0 or e <= s:
                raise ValueError(f'Invalid block interval after coordinate conversion: {line.rstrip()}')
            rec=(c,s,e,bid)
            if max_block_bp is not None and e-s>max_block_bp:
                large.append(rec)
            else:
                rows.append(rec)
    rows.sort(key=lambda x:(x[0],x[1],x[2]))
    large.sort(key=lambda x:(x[0],x[1],x[2]))
    return rows, large


def distance_bin_label(d:int, edges:List[int]):
    for a,b in zip(edges[:-1],edges[1:]):
        if a <= d < b:
            def fmt(x):
                if x==0: return '0'
                if x%1000000==0: return f'{x//1000000}Mb'
                return f'{x//1000}kb'
            return f'{fmt(a)}-{fmt(b)}'
    return None


def write_json(path,obj):
    with open(path,'w') as fh:
        json.dump(obj,fh,indent=2,sort_keys=True)



# Original functions: prepare_inputs.py

def read_genes(gff, chroms):
    genes={}; chrom_len=defaultdict(int)
    with open_text(gff) as fh:
        for line in fh:
            if not line.strip() or line.startswith('#'): continue
            f=line.rstrip('\n').split('\t')
            if len(f)!=9: continue
            c=norm_chr(f[0]); chrom_len[c]=max(chrom_len[c], int(f[4]))
            if c not in chroms or f[2] != 'gene': continue
            start1,end1=int(f[3]),int(f[4]); strand=f[6]
            attrs=parse_gff_attrs(f[8]); gid=attrs.get('ID') or attrs.get('Name')
            if not gid: continue
            if gid.startswith('gene:'): gid=gid[5:]
            s0=start1-1; e0=end1
            tss0=s0 if strand=='+' else e0-1
            genes[gid]=(c,s0,e0,strand,tss0)
    return genes,chrom_len


def promoter_interval(gene, promoter_bp):
    c,s,e,strand,tss=gene
    if strand=='+': return max(0,tss-promoter_bp), tss
    return tss+1, tss+1+promoter_bp


def standardized_difference(values, reference, scales):
    difference = np.asarray(values, dtype=float) - np.asarray(reference, dtype=float)
    result = np.zeros_like(difference)
    np.divide(difference, scales, out=result, where=scales > 0)
    constant = scales == 0
    if np.any(constant):
        result[..., constant] = np.where(difference[..., constant] == 0, 0.0,
                                         np.copysign(np.inf, difference[..., constant]))
    return result


def background_candidate_mask(candidate_values, real_values, candidate_valid, scales, caliper):
    return (np.asarray(candidate_valid, dtype=bool)
            & np.all(np.isfinite(candidate_values), axis=1)
            & np.all(np.abs(standardized_difference(candidate_values, real_values, scales)) <= caliper + 1e-12, axis=1))


def distance_candidate_indices(d, tss, gid, acr_ids, acr_start, acr_end, acr_mid,
                               genes, known_partners, cfg):
    """Original distance/promoter/known-partner rules, before background filtering."""
    cal = cfg.getfloat('random_aglink_null', 'distance_caliper_fraction')
    threshold = cfg.getint('random_aglink_null', 'short_distance_threshold_bp', fallback=0)
    bins = int_list(cfg['analysis']['distance_bins_bp'])
    short = threshold > 0 and d < threshold
    if short:
        bounds = next(((a, b) for a, b in zip(bins[:-1], bins[1:]) if a <= d < b), None)
        if bounds is None:
            raise ValueError(f'No distance stratum contains {d} bp')
        lo, hi = bounds
    else:
        lo = max(1, int(math.floor(d * (1.0 - cal))))
        hi = max(lo, int(math.ceil(d * (1.0 + cal))))
    ps, pe = promoter_interval(genes[gid], cfg.getint('analysis', 'promoter_bp'))
    exclude_known = cfg.getboolean('random_aglink_null', 'exclude_known_partners')
    result = []
    seen = set()
    for left, right in ((tss-hi, tss-lo), (tss+lo, tss+hi)):
        low = int(np.searchsorted(acr_mid, left, side='left'))
        high = int(np.searchsorted(acr_mid, right, side='right'))
        for j in range(low, high):
            if j in seen:
                continue
            if exclude_known and acr_ids[j] in known_partners.get(gid, set()):
                continue
            if interval_overlap(int(acr_start[j]), int(acr_end[j]), ps, pe):
                continue
            rd = abs(int(acr_mid[j]) - tss)
            if short:
                if not lo <= rd < hi:
                    continue
            elif d <= 0 or not 1.0-cal-1e-12 <= rd/d <= 1.0+cal+1e-12:
                continue
            seen.add(j)
            result.append(j)
    if short:
        result.sort(key=lambda j: (abs(abs(int(acr_mid[j])-tss)-d), j))
    else:
        result.sort()
    return result, short, (lo, hi)


def empty_balance():
    return {key: np.zeros(4, dtype=np.int64 if key.startswith('n_') else float)
            for key in ('n_before', 'real_before', 'random_before', 'abs_before', 'max_before',
                        'n_after', 'real_after', 'random_after', 'abs_after', 'max_after')}


def make_null_design(chrom, eval_rows, acr_chr, genes, known_partners, cfg, out_npz, out_audit):
    mincand=cfg.getint('random_aglink_null','min_candidates_per_edge')
    short_mincand=cfg.getint('random_aglink_null','short_min_candidates_per_edge',fallback=mincand)
    short_pool=cfg.getint('random_aglink_null','short_candidate_pool_size',fallback=max(short_mincand,mincand))
    distance_bins=int_list(cfg['analysis']['distance_bins_bp'])
    bg_caliper=cfg.getfloat('background','match_caliper_sd',fallback=0.2)
    min_coverage=cfg.getfloat('background','min_recomb_coverage',fallback=0.5)
    if not math.isfinite(bg_caliper) or bg_caliper < 0:
        raise ValueError('Background matching caliper must be finite and nonnegative')
    # local ACR catalogue sorted by midpoint
    acr_ids=np.array([x[0] for x in acr_chr],dtype='U128')
    acr_start=np.array([x[2] for x in acr_chr],dtype=np.int64)
    acr_end=np.array([x[3] for x in acr_chr],dtype=np.int64)
    acr_mid=np.array([x[4] for x in acr_chr],dtype=np.int64)
    acr_index={aid:i for i,aid in enumerate(acr_ids.tolist())}
    bgpath=Path(cfg['paths']['work_dir'])/'02_background'/chrom/'acr_background.npz'
    with np.load(bgpath,allow_pickle=False) as background:
        if not np.array_equal(background['acr_ids'],acr_ids) or not np.array_equal(background['acr_mid0'],acr_mid):
            raise ValueError(f'Background ACR identity/order mismatch: {chrom}')
        features=np.column_stack([background[name] for name in FEATURE_NAMES])
        feature_valid=np.all(np.isfinite(features),axis=1)&(background['recomb_coverage']>=min_coverage)
    scales=features[feature_valid].std(axis=0,ddof=1) if feature_valid.sum()>1 else np.zeros(4)
    if feature_valid.any(): scales[np.ptp(features[feature_valid],axis=0)==0]=0
    if not np.all(np.isfinite(scales)):
        raise ValueError(f'Invalid background feature scales on {chrom}')
    edge_indices=[]; gene_pos=[]; real_acr_idx=[]; dist=[]; ga=[]; counts=[]; supports=[]; choice_rows=[]
    before_counts=[]; reasons=[]
    audits=[]; B=cfg.getint('random_aglink_null','n_random_sets'); seed0=cfg.getint('random_aglink_null','base_seed')
    strata=['ALL']+[distance_bin_label((a+b)//2,distance_bins) for a,b in zip(distance_bins[:-1],distance_bins[1:])]
    support_audit={s:Counter() for s in strata}
    balance={s:empty_balance() for s in strata}
    for x in eval_rows:
        if x['chrom']!=chrom: continue
        d=int(x['distance_bp']); tss=int(x['gene_tss0']); gid=x['gene_id']; aid=x['acr_id']
        ri=acr_index[aid]
        distance_candidates,is_short,bounds=distance_candidate_indices(d,tss,gid,acr_ids,acr_start,acr_end,acr_mid,genes,known_partners,cfg)
        minimum=short_mincand if is_short else mincand
        before_n=len(distance_candidates)
        before_pool=distance_candidates[:short_pool] if is_short else distance_candidates
        ci=np.asarray(distance_candidates,dtype=np.int32)
        if feature_valid[ri] and before_n:
            passing=background_candidate_mask(features[ci],features[ri],feature_valid[ci],scales,bg_caliper)
            candidate_idxs=ci[passing].tolist()
        else:
            candidate_idxs=[]
        n=len(candidate_idxs)
        sampling_idxs=candidate_idxs[:short_pool] if is_short else candidate_idxs
        if not feature_valid[ri]: reason='true_background_missing'
        elif before_n<minimum: reason='insufficient_distance_candidates'
        elif n<minimum: reason='insufficient_background_candidates'
        else: reason='supported'
        support=int(reason=='supported')
        edge_indices.append(int(x['edge_index'])); gene_pos.append(tss); real_acr_idx.append(acr_index[aid]); dist.append(d)
        ga.append(int(x['group_all_significant']))
        counts.append(n); supports.append(support); before_counts.append(before_n); reasons.append(reason)
        edge_strata=['ALL']
        label=distance_bin_label(d,distance_bins)
        if label is not None: edge_strata.append(label)
        for s in edge_strata:
            support_audit[s]['n_real_edges']+=1
            support_audit[s]['n_distance_eligible']+=int(before_n>=minimum)
            support_audit[s]['n_'+reason]+=1
        if support:
            rng=np.random.default_rng(stable_seed(seed0,chrom,int(x['edge_index'])))
            carr=np.asarray(sampling_idxs,dtype=np.int32)
            chosen=carr[rng.integers(0,len(carr),size=B)]
            selected_d=np.abs(acr_mid[chosen]-tss)
            if is_short:
                distance_ok=np.all((selected_d>=bounds[0])&(selected_d<bounds[1]))
            else:
                cal=cfg.getfloat('random_aglink_null','distance_caliper_fraction')
                distance_ok=np.all((selected_d/d>=1-cal-1e-12)&(selected_d/d<=1+cal+1e-12))
            z_after=standardized_difference(features[chosen],features[ri],scales)
            if not distance_ok or not feature_valid[chosen].all() or not np.all(np.abs(z_after)<=bg_caliper+1e-12):
                raise RuntimeError(f'Random pair failed distance/background validation: {chrom} {x["edge_id"]}')
            choice_rows.append(chosen)
            pre=features[np.asarray(before_pool,dtype=np.int32)]
            pre_finite=np.isfinite(pre)
            pre_count=pre_finite.sum(axis=0)
            pre_mean=np.full(4,np.nan)
            np.divide(np.where(pre_finite,pre,0).sum(axis=0),pre_count,out=pre_mean,where=pre_count>0)
            usable=np.isfinite(pre_mean)
            pre_z=standardized_difference(pre,features[ri],scales)
            pre_max=np.max(np.where(pre_finite,np.abs(pre_z),0),axis=0)
            after_mean=features[chosen].mean(axis=0)
            before_abs=np.abs(standardized_difference(pre_mean,features[ri],scales))
            after_abs=np.abs(standardized_difference(after_mean,features[ri],scales))
            for s in edge_strata:
                b=balance[s]
                b['n_before']+=usable.astype(np.int64)
                b['real_before']+=np.where(usable,features[ri],0)
                b['random_before']+=np.where(usable,pre_mean,0)
                b['abs_before']+=np.where(usable,before_abs,0)
                b['max_before']=np.maximum(b['max_before'],pre_max)
                b['n_after']+=1; b['real_after']+=features[ri]; b['random_after']+=after_mean
                b['abs_after']+=after_abs; b['max_after']=np.maximum(b['max_after'],np.max(np.abs(z_after),axis=0))
        else:
            choice_rows.append(np.full(B,-1,dtype=np.int32))
        audits.append((x['edge_id'],chrom,gid,aid,d,n,support,before_n,n,reason,int(feature_valid[ri])))
    np.savez_compressed(out_npz,
        acr_ids=acr_ids, acr_start0=acr_start, acr_end0=acr_end, acr_mid0=acr_mid,
        edge_indices=np.asarray(edge_indices,dtype=np.int64), gene_pos0=np.asarray(gene_pos,dtype=np.int64),
        real_acr_idx=np.asarray(real_acr_idx,dtype=np.int32), distance_bp=np.asarray(dist,dtype=np.int64),
        group_all_significant=np.asarray(ga,dtype=np.int8),
        candidate_count=np.asarray(counts,dtype=np.int32), support=np.asarray(supports,dtype=np.int8),
        distance_candidate_count=np.asarray(before_counts,dtype=np.int32),support_reason=np.asarray(reasons,dtype='U40'),
        background_feature_names=np.asarray(FEATURE_NAMES,dtype='U32'),background_scales=scales,
        background_match_caliper_sd=np.asarray(bg_caliper))
    choice_mat=np.vstack(choice_rows).astype(np.int32,copy=False) if choice_rows else np.empty((0,B),dtype=np.int32)
    np.save(Path(out_npz).with_name('random_choice_idx.npy'),choice_mat,allow_pickle=False)
    with open(out_audit,'w') as out:
        out.write('edge_id\tchrom\tgene_id\treal_acr_id\tdistance_bp\tn_random_acr_candidates\tsupported\tn_distance_candidates\tn_background_candidates\tsupport_reason\ttrue_background_valid\n')
        for r in audits: out.write('\t'.join(map(str,r))+'\n')
    with open(Path(out_npz).with_name('match_support_by_distance.tsv'),'w',newline='') as out:
        names=['chrom','distance_stratum','n_real_edges','n_distance_eligible','n_true_background_missing','n_insufficient_distance_candidates','n_insufficient_background_candidates','n_supported','support_fraction']
        writer=csv.DictWriter(out,fieldnames=names,delimiter='\t',lineterminator='\n'); writer.writeheader()
        for s in strata:
            a=support_audit[s]
            row={k:int(a[k]) for k in names[2:-1]}
            row.update(chrom=chrom,distance_stratum=s,support_fraction=a['n_supported']/a['n_real_edges'] if a['n_real_edges'] else 'NA')
            writer.writerow(row)
    with open(Path(out_npz).with_name('background_balance.tsv'),'w',newline='') as out:
        writer=csv.writer(out,delimiter='\t',lineterminator='\n')
        writer.writerow(['chrom','distance_stratum','stage','feature','n_supported_edges_with_values','real_mean','random_mean','paired_standardized_mean_difference_real_minus_random','mean_absolute_expected_pair_standardized_difference','max_absolute_individual_pair_standardized_difference','chromosome_feature_sample_sd','caliper_sd'])
        for s in strata:
            for stage in ('before','after'):
                b=balance[s]
                for k,name in enumerate(FEATURE_NAMES):
                    n=int(b['n_'+stage][k])
                    real=b['real_'+stage][k]/n if n else float('nan')
                    rand=b['random_'+stage][k]/n if n else float('nan')
                    delta=(real-rand)/scales[k] if scales[k]>0 else (0.0 if real==rand else float('nan'))
                    writer.writerow([chrom,s,stage,name,n,real,rand,delta,b['abs_'+stage][k]/n if n else 'NA',b['max_'+stage][k] if n else 'NA',scales[k],bg_caliper])
    return dict(n_edges=len(edge_indices), n_supported=sum(supports), min_candidates=min(counts) if counts else 0,
                median_candidates=float(np.median(counts)) if counts else 0.0,
                n_distance_eligible=support_audit['ALL']['n_distance_eligible'],support_reasons=dict(Counter(reasons)),
                background_feature_scales=dict(zip(FEATURE_NAMES,scales.tolist())),background_match_caliper_sd=bg_caliper,
                background_feature_transform='raw values; sample SD across chromosome ACRs with all features finite and sufficient observed map coverage',
                background_n_valid_acrs=int(feature_valid.sum()),random_choices_validated=True)



# Original functions: blocks_utils.py

def plink_blocks_to_bed(det,chrom,outbed,max_block):
    kept=large=0; lengths=[]
    with open(outbed,'w') as out:
        if not Path(det).exists() or Path(det).stat().st_size==0:
            return kept,large,lengths
        with open(det) as fh:
            header=fh.readline().strip().split()
            if not header:
                return kept,large,lengths
            idx={k:i for i,k in enumerate(header)}
            if 'BP1' not in idx or 'BP2' not in idx:
                raise RuntimeError(f'Cannot find BP1/BP2 in {det}: {header}')
            for n,line in enumerate(fh,1):
                if not line.strip(): continue
                f=line.split(); bp1=int(f[idx['BP1']]); bp2=int(f[idx['BP2']]); s=bp1-1; e=bp2
                if e<=s: continue
                L=e-s
                if L>max_block: large+=1; continue
                out.write(f'{chrom}\t{s}\t{e}\tsim_block{n:07d}:{chrom}_{s}_{e}\n'); kept+=1; lengths.append(L)
    return kept,large,lengths


def call_plink(vcf,outprefix,chrom,cfg):
    staged=Path(cfg['paths']['work_dir'])/'01_bin'/'plink1.9'
    plink=str(staged if staged.exists() else Path(cfg['paths']['plink']))
    allow=['--allow-extra-chr'] if cfg.getboolean('plink','allow_extra_chr') else []
    threads=cfg.getint('plink','threads',fallback=cfg.getint('simulation','threads',fallback=1))
    if os.environ.get('SLURM_CPUS_PER_TASK'):
        threads=min(threads,int(os.environ['SLURM_CPUS_PER_TASK']))
    if threads<1: raise ValueError('PLINK threads must be positive')
    memory=cfg.getint('plink','memory_mb',fallback=3072)
    if memory<64: raise ValueError('PLINK memory_mb must be at least 64')
    resources=['--threads',str(threads),'--memory',str(memory)]
    def run_logged(cmd,stage):
        # PLINK overwrites its native .log on its second invocation.  Preserve
        # both command streams and native logs, including all restart attempts.
        stamp=datetime.datetime.now(datetime.timezone.utc).isoformat()
        with open(str(outprefix)+'.commands.jsonl','a') as log:
            log.write(json.dumps(dict(started_at=stamp,stage=stage,command=cmd))+'\n')
        with open(str(outprefix)+'.'+stage+'.stdout.log','a') as log:
            log.write('\n'+stamp+' '+shlex.join(cmd)+'\n'); log.flush()
            try:
                subprocess.check_call(cmd,stdout=log,stderr=subprocess.STDOUT)
            finally:
                native=Path(str(outprefix)+'.log')
                if native.exists():
                    with open(str(outprefix)+'.'+stage+'.plink.log','a') as saved:
                        saved.write('\n'+stamp+'\n'+native.read_text())
    cmd=[plink,'--vcf',str(vcf),'--double-id']+allow+resources+shlex.split(cfg['plink']['makebed_args'])+['--out',str(outprefix)]
    run_logged(cmd,'makebed')
    block_args=shlex.split(cfg['plink']['block_args'])
    if '--blocks-max-kb' in block_args:
        raise ValueError('Put chromosome-specific maxima in block_max_kb_<CHR>, not in block_args')
    max_kb=cfg.getint('plink',f'block_max_kb_{chrom}')
    cmd=[plink,'--bfile',str(outprefix)]+allow+resources+block_args+['--blocks-max-kb',str(max_kb),'--out',str(outprefix)]
    run_logged(cmd,'blocks')
    return Path(str(outprefix)+'.blocks.det')



# Original functions: score_blocks.py

def shared_block_span(pos0, blocks, starts, max_block_bp):
    if not blocks: return (-1,-1)
    i=bisect.bisect_right(starts,pos0)
    left=None; right=None
    j=i-1
    floor=pos0-max_block_bp
    while j>=0 and blocks[j][1]>=floor:
        _,s,e,_=blocks[j]
        if s<=pos0<e:
            left=s if left is None else min(left,s); right=e if right is None else max(right,e)
        j-=1
    return (-1,-1) if left is None else (left,right)


def score_blockset(cfg, chrom, blocks_path, out_prefix, source='sim'):
    maxb=cfg.getint('analysis','max_block_bp'); bins=int_list(cfg['analysis']['distance_bins_bp'])
    groups=csv_list(cfg['analysis']['evaluation_groups']); strata=['ALL']+[distance_bin_label((a+b)//2,bins) for a,b in zip(bins[:-1],bins[1:])]
    # deduplicate labels defensively
    strata=list(dict.fromkeys([x for x in strata if x]))
    work=Path(cfg['paths']['work_dir']); design=np.load(work/'02_prepared'/chrom/'random_null_design.npz',allow_pickle=False)
    mode=cfg.get('analysis','real_blocks_coordinate_mode',fallback='one_based_inclusive') if source.startswith('observed') else 'bed0'
    blocks,large=read_blocks_bed(blocks_path,chrom,maxb,mode); starts=[b[1] for b in blocks]
    edge_idx=design['edge_indices']; gene=design['gene_pos0']; real_idx=design['real_acr_idx']; dist=design['distance_bp']; support=design['support'].astype(bool)
    acr_mid=design['acr_mid0']; random_choice=np.load(work/'02_prepared'/chrom/'random_choice_idx.npy',mmap_mode='r')
    membership={g:design['group_'+g.replace('.','_')].astype(bool) for g in groups}
    B=cfg.getint('random_aglink_null','n_random_sets')
    G=len(groups); S=len(strata)
    real_all_same=np.zeros((G,S),dtype=np.int64); n_all_edges=np.zeros((G,S),dtype=np.int64)
    real_same=np.zeros((G,S),dtype=np.int64); n_edges=np.zeros((G,S),dtype=np.int64); random_same=np.zeros((G,S,B),dtype=np.int64)
    # Cache block spans by gene TSS, then score all edges in vectorized chunks.
    span_cache={}
    left=np.empty(len(edge_idx),dtype=np.int64); right=np.empty(len(edge_idx),dtype=np.int64)
    for ei,pos in enumerate(gene):
        key=int(pos)
        if key not in span_cache: span_cache[key]=shared_block_span(key,blocks,starts,maxb)
        left[ei],right[ei]=span_cache[key]
    real_pos=acr_mid[real_idx]
    real_hit=(left>=0)&(real_pos>=left)&(real_pos<right)
    stratum_masks=[np.ones(len(edge_idx),dtype=bool)]
    for low,high in zip(bins[:-1],bins[1:]): stratum_masks.append((dist>=low)&(dist<high))
    chunk_size=2048
    for gi,g in enumerate(groups):
        gm=membership[g]
        for si,sm in enumerate(stratum_masks):
            all_mask=gm&sm
            n_all_edges[gi,si]=int(all_mask.sum()); real_all_same[gi,si]=int(real_hit[all_mask].sum())
            matched=all_mask&support
            idx=np.flatnonzero(matched)
            n_edges[gi,si]=len(idx); real_same[gi,si]=int(real_hit[idx].sum())
            for start in range(0,len(idx),chunk_size):
                ii=idx[start:start+chunk_size]
                chosen=np.asarray(random_choice[ii,:],dtype=np.int32)
                if np.any(chosen<0): raise RuntimeError('Supported edge has an invalid random choice')
                positions=acr_mid[chosen]
                hit=(left[ii,None]>=0)&(positions>=left[ii,None])&(positions<right[ii,None])
                random_same[gi,si,:]+=hit.sum(axis=0,dtype=np.int64)
    out_prefix=Path(out_prefix); ensure_dir(out_prefix.parent)
    np.savez_compressed(str(out_prefix)+'.npz',group_names=np.asarray(groups,dtype='U32'),strata_names=np.asarray(strata,dtype='U32'),real_all_same=real_all_same,n_all_edges=n_all_edges,real_same=real_same,n_edges=n_edges,random_same=random_same)
    with open(str(out_prefix)+'.tsv','w') as out:
        out.write('source\tchrom\tgroup\tdistance_stratum\tn_all_edges\treal_all_same\treal_all_fraction\tn_supported_edges\treal_supported_same\treal_supported_fraction\trandom_mean\trandom_q025\trandom_q975\tdelta_supported_real_minus_random\tempirical_p_random_ge_real\n')
        for gi,g in enumerate(groups):
            for si,s in enumerate(strata):
                na=int(n_all_edges[gi,si]); n=int(n_edges[gi,si])
                if na==0: continue
                rfa=real_all_same[gi,si]/na
                if n>0:
                    rf=real_same[gi,si]/n; rv=random_same[gi,si,:]/n; p=(1+int(np.sum(rv>=rf)))/(B+1)
                    out.write(f'{source}\t{chrom}\t{g}\t{s}\t{na}\t{real_all_same[gi,si]}\t{rfa:.10g}\t{n}\t{real_same[gi,si]}\t{rf:.10g}\t{rv.mean():.10g}\t{np.quantile(rv,.025):.10g}\t{np.quantile(rv,.975):.10g}\t{(rf-rv.mean()):.10g}\t{p:.10g}\n')
                else:
                    out.write(f'{source}\t{chrom}\t{g}\t{s}\t{na}\t{real_all_same[gi,si]}\t{rfa:.10g}\t0\t0\tnan\tnan\tnan\tnan\tnan\tnan\n')
    return dict(n_valid_blocks=len(blocks),n_large_blocks_ignored=len(large),npz=str(out_prefix)+'.npz',tsv=str(out_prefix)+'.tsv')


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def atomic_text(path, text):
    path=Path(path)
    tmp=path.with_name(path.name+f'.tmp.{os.getpid()}')
    with tmp.open('w') as handle:
        handle.write(text)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp,path)


def atomic_json(path, obj):
    atomic_text(path,json.dumps(obj,sort_keys=True,indent=2)+'\n')


def read_row(path, trajectory_id):
    with Path(path).open() as handle:
        for row in csv.DictReader(handle,delimiter='\t'):
            if int(row['trajectory_id'])==trajectory_id:
                return row
    raise KeyError(f'No trajectory {trajectory_id} in {path}')


def sha256(path):
    digest=hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda:handle.read(1024*1024),b''):
            digest.update(chunk)
    return digest.hexdigest()


def scientific_signature(cfg, row, inputs):
    # Resource allocations and output retention may change without invalidating
    # completed scientific results.  Scientific parameters and input bytes may not.
    sections={name:dict(cfg[name]) for name in
              ('analysis','random_aglink_null','simulation','plink','ld_qc','background') if name in cfg}
    for name,keys in {'simulation':('threads','task_table_filename'),
                      'plink':('threads','memory_mb')}.items():
        for key in keys: sections.get(name,{}).pop(key,None)
    payload=dict(parameters=sections,trajectory=row,
                 input_sha256={name:sha256(path) for name,path in inputs.items()})
    digest=hashlib.sha256(json.dumps(payload,sort_keys=True).encode()).hexdigest()
    return digest,payload


def complete(path, signature):
    path=Path(path)
    marker=path/'DONE'
    if not marker.exists(): return False
    try: saved=json.loads(marker.read_text())
    except (ValueError,OSError) as exc:
        raise RuntimeError(f'Invalid completion marker: {marker}') from exc
    if saved.get('signature')!=signature:
        raise RuntimeError(f'Completed output has different inputs/parameters: {path}')
    gen=int(saved['generation'])
    required=[f'score_G{gen}.npz',f'score_G{gen}.tsv','block_architecture.tsv',
              'ld_qc.tsv','maf_qc.tsv',f'G{gen}.plink.blocks.le10Mb.bed']
    if not all((path/name).is_file() for name in required):
        raise RuntimeError(f'Completed output is missing required files: {path}')
    return True


@contextlib.contextmanager
def timed_step(log_path, step, **metadata):
    started=utc_now(); tick=time.perf_counter(); status='ok'; error=None
    try:
        yield
    except BaseException as exc:
        status='failed'; error=repr(exc)
        raise
    finally:
        event=dict(step=step,started_at=started,finished_at=utc_now(),
                   elapsed_seconds=time.perf_counter()-tick,status=status,**metadata)
        if error is not None: event['error']=error
        with Path(log_path).open('a') as log:
            log.write(json.dumps(event,sort_keys=True)+'\n')
        print(json.dumps(event,sort_keys=True),flush=True)


def validate_vcf(path, expected_variants, expected_samples):
    """A successful exit alone does not detect an incomplete disk write."""
    n=0; samples=None
    with Path(path).open() as handle:
        for line in handle:
            if not line.endswith('\n'):
                raise RuntimeError(f'Truncated VCF line in {path}')
            if line.startswith('#CHROM\t'):
                samples=len(line.rstrip('\n').split('\t'))-9
            elif not line.startswith('#'):
                n+=1
    if n!=expected_variants or samples!=expected_samples:
        raise RuntimeError(f'Invalid VCF {path}: {n}/{expected_variants} variants, '
                           f'{samples}/{expected_samples} samples')


def split_qc(source, destination, generation, allow_empty=False):
    with Path(source).open() as handle:
        reader=csv.DictReader(handle,delimiter='\t')
        fields=reader.fieldnames
        if not fields or 'generation' not in fields:
            raise RuntimeError(f'Invalid QC header: {source}')
        rows=[row for row in reader if int(row['generation'])==generation]
    if not rows and not allow_empty:
        raise RuntimeError(f'No QC for generation {generation}: {source}')
    text='\t'.join(fields)+'\n'
    text+=''.join('\t'.join(row[key] for key in fields)+'\n' for row in rows)
    atomic_text(destination,text)


def validate_score(path, n_random_sets):
    with np.load(path,allow_pickle=False) as score:
        real=score['real_same']; n=score['n_edges']; random=score['random_same']
        if real.shape!=n.shape or random.shape!=n.shape+(n_random_sets,):
            raise RuntimeError(f'Unexpected score dimensions: {path}')
        for values,denominator in ((real,n),(random,n[...,None])):
            if np.any(values<0) or np.any(values>denominator):
                raise RuntimeError(f'Invalid counts in {path}')


def cleanup_checkpoint(cfg, directory, generation):
    prefix=directory/f'G{generation}.plink'
    if not cfg.getboolean('output','keep_checkpoint_vcf',fallback=False):
        (directory/f'G{generation}.vcf').unlink(missing_ok=True)
    if not cfg.getboolean('output','keep_plink_binary',fallback=False):
        for suffix in ('.bed','.bim','.fam','.nosex'):
            Path(str(prefix)+suffix).unlink(missing_ok=True)
    # Preserve all command/native/stdout logs even when compact transient
    # genotype and block-caller files are removed.
    if not cfg.getboolean('output','keep_plink_log',fallback=False):
        for suffix in ('.hh','.blocks','.blocks.det'):
            Path(str(prefix)+suffix).unlink(missing_ok=True)


def allocate_vcf_directory(trajectory, shared_directory, signature):
    """Keep large disposable genotypes off the shared filesystem where possible."""
    candidates=[]
    if os.environ.get('SLURM_TMPDIR'): candidates.append(Path(os.environ['SLURM_TMPDIR']))
    candidates.append(Path('/tmp'))
    failures=[]
    for root in dict.fromkeys(candidates):
        try:
            if not root.is_dir() or shutil.disk_usage(root).free<2*1024**3:
                failures.append(dict(root=str(root),reason='missing directory or less than 2 GiB available'))
                continue
            path=Path(tempfile.mkdtemp(prefix='neutral_genotypes_',dir=str(root)))
            owner=dict(signature=signature,trajectory=str(trajectory),host=socket.gethostname(),
                       created_at=utc_now(),job_id=os.environ.get('SLURM_JOB_ID'))
            atomic_json(path/'OWNER.json',owner)
            return path,True,failures
        except OSError as exc:
            failures.append(dict(root=str(root),reason=repr(exc)))
    return shared_directory,False,failures


def remove_private_scratch(path, trajectory, signature, recorded_host):
    """Only remove the private directory created and marked by this trajectory."""
    if recorded_host!=socket.gethostname() or not path.exists(): return
    owner=path/'OWNER.json'
    if not owner.is_file(): raise RuntimeError(f'Refusing to clean unmarked scratch: {path}')
    saved=json.loads(owner.read_text())
    if saved.get('signature')!=signature or saved.get('trajectory')!=str(trajectory) or saved.get('host')!=recorded_host:
        raise RuntimeError(f'Refusing to clean scratch owned by another trajectory: {path}')
    shutil.rmtree(path)


def run_trajectory(cfg, trajectory_id):
    work=Path(cfg['paths']['work_dir']).resolve()
    prep=work/'02_prepared'
    row=read_row(prep/'trajectory_task_table.tsv',trajectory_id)
    chrom=row['chrom']; rep=int(row['replicate']); seed=int(row['seed'])
    generations=sorted(set(int_list(cfg['simulation']['checkpoints'])))
    if not generations: raise ValueError('No checkpoints configured')
    trajectory=work/'03_trajectories'/f'chr{chrom}_rep{rep:03d}'
    ensure_dir(trajectory); ensure_dir(work/'03_neutral')
    lock=(trajectory/'LOCK').open('a')
    try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError as exc:
        raise RuntimeError(f'Trajectory is already running: {trajectory}') from exc
    # Keep the lock alive for the whole simulation and all checkpoint scoring.
    with lock:
        threads=cfg.getint('simulation','threads')
        if threads<1: raise ValueError('Thread count must be positive')
        cfg['simulation']['threads']=str(threads)
        cfg['plink']['threads']=str(threads)
        os.environ['OMP_NUM_THREADS']=str(threads)
        os.environ['OPENBLAS_NUM_THREADS']='1'
        os.environ['MKL_NUM_THREADS']='1'
        exe=work/'01_bin'/'forward_neutral'
        map_path=work/'02_background'/chrom/'recombination_map.tsv'
        snps=prep/chrom/'snps.tsv'
        for path in (exe,map_path,snps):
            if not path.is_file() or path.stat().st_size==0:
                raise FileNotFoundError(f'Missing required input: {path}')
        inputs=dict(executable=exe,recombination_map=map_path,snps=snps,
                    chromosome_length=prep/chrom/'chromosome_length.txt',
                    null_design=prep/chrom/'random_null_design.npz',
                    random_choices=prep/chrom/'random_choice_idx.npy',
                    analysis_code=Path(__file__),
                    plink_executable=work/'01_bin'/'plink1.9')
        if (prep/chrom/'ld_pairs.tsv').exists(): inputs['ld_pairs']=prep/chrom/'ld_pairs.tsv'
        signature,provenance=scientific_signature(cfg,row,inputs)
        marker=trajectory/'SIMULATOR_DONE.json'
        saved_marker=json.loads(marker.read_text()) if marker.exists() else {}
        if marker.exists() and saved_marker.get('signature')!=signature:
            raise RuntimeError(f'Simulator completion marker has changed inputs: {marker}')
        atomic_json(trajectory/'provenance.json',dict(signature=signature,**provenance))
        audit=trajectory/'timing.jsonl'
        finals={g:work/'03_neutral'/f'chr{chrom}_rep{rep:03d}_G{g:03d}' for g in generations}
        pending=[g for g in generations if not complete(finals[g],signature)]
        if not pending:
            for g in generations: cleanup_checkpoint(cfg,finals[g],g)
            if saved_marker.get('local_scratch'):
                remove_private_scratch(Path(saved_marker['vcf_directory']),trajectory,signature,saved_marker['scratch_host'])
            atomic_json(trajectory/'DONE',dict(signature=signature,finished_at=utc_now(),checkpoints=generations))
            print(f'Already complete: {trajectory}',flush=True)
            return
        # Small checkpoint QC remains shared even when the large genotype VCFs
        # live on the current compute node.  Local paths are never assumed to
        # identify the same files on a different host after requeue/resubmission.
        simdir=trajectory/'simulator_output'; ensure_dir(simdir)
        vcfdir=Path(saved_marker.get('vcf_directory',str(simdir)))
        local_scratch=bool(saved_marker.get('local_scratch',False))
        scratch_host=saved_marker.get('scratch_host',socket.gethostname())
        def available_vcf(g):
            local_accessible=(not local_scratch or scratch_host==socket.gethostname())
            candidates=([vcfdir/f'G{g}.vcf'] if local_accessible else [])
            candidates.extend([finals[g]/f'G{g}.vcf',
                               finals[g].with_name('.'+finals[g].name+'.staging')/f'G{g}.vcf'])
            return next((p for p in candidates if p.exists() and (local_accessible or not p.is_symlink())),None)
        # Partial VCFs are never reused without a successful, validated whole
        # simulator completion marker.  Missing unfinished checkpoints force a
        # deterministic regeneration, leaving completed checkpoint results intact.
        simulate=(not marker.exists() or any(available_vcf(g) is None for g in pending))
        if simulate:
            marker.unlink(missing_ok=True)
            vcfdir,local_scratch,scratch_failures=allocate_vcf_directory(trajectory,simdir,signature)
            scratch_host=socket.gethostname()
            scratch_record=dict(vcf_directory=str(vcfdir),local_scratch=local_scratch,
                                scratch_host=scratch_host,created_at=utc_now(),fallback_reasons=scratch_failures)
            # Preserve every attempt's scratch location, including failed runs.
            with (trajectory/'scratch_history.jsonl').open('a') as handle:
                handle.write(json.dumps(scratch_record,sort_keys=True)+'\n')
            chrlen=int((prep/chrom/'chromosome_length.txt').read_text().strip())
            command=[str(exe),'--snps',str(snps),'--outdir',str(vcfdir),'--chrom',chrom,
                     '--N',cfg['simulation']['population_size'],'--outN',cfg['simulation']['output_sample_size'],
                     '--recomb',cfg['simulation']['recomb_rate_per_bp'],'--recomb-map',str(map_path),
                     '--selfing',cfg['simulation']['selfing_rate'],'--chrlen',str(chrlen),
                     '--seed',str(seed),'--threads',str(threads),
                     '--checkpoints',','.join(map(str,generations)),'--emit-vcf','1',
                     '--ld-out',str(simdir/'ld_qc.tsv'),'--maf-out',str(simdir/'maf_qc.tsv')]
            if 'ld_pairs' in inputs: command+=['--ldpairs',str(inputs['ld_pairs'])]
            with timed_step(audit,'simulate',command=command):
                with (trajectory/'simulator.stdout.log').open('a') as log:
                    log.write('\n'+utc_now()+' '+json.dumps(command)+'\n'); log.flush()
                    subprocess.check_call(command,stdout=log,stderr=subprocess.STDOUT)
            with timed_step(audit,'validate_simulator_outputs'):
                with snps.open() as handle: n_variants=sum(1 for _ in handle)-1
                n_samples=min(cfg.getint('simulation','population_size'),cfg.getint('simulation','output_sample_size'))
                for g in generations:
                    validate_vcf(vcfdir/f'G{g}.vcf',n_variants,n_samples)
                if not (simdir/'ld_qc.tsv').is_file() or not (simdir/'maf_qc.tsv').is_file():
                    raise RuntimeError('Simulator did not emit the required QC files')
            atomic_json(marker,dict(signature=signature,completed_at=utc_now(),checkpoints=generations,
                                   vcf_directory=str(vcfdir),local_scratch=local_scratch,scratch_host=scratch_host))
        for g in pending:
            final=finals[g]
            stage=final.with_name('.'+final.name+'.staging'); ensure_dir(stage)
            vcf=available_vcf(g)
            if vcf is None: raise RuntimeError(f'Missing validated VCF for generation {g}')
            staged_vcf=stage/f'G{g}.vcf'
            if vcf!=staged_vcf:
                if local_scratch and vcf.resolve().parent==vcfdir.resolve():
                    staged_vcf.unlink(missing_ok=True)
                    staged_vcf.symlink_to(vcf.resolve())
                else:
                    os.replace(vcf,staged_vcf)
            split_qc(simdir/'ld_qc.tsv',stage/'ld_qc.tsv',g,allow_empty=True)
            split_qc(simdir/'maf_qc.tsv',stage/'maf_qc.tsv',g)
            prefix=stage/f'G{g}.plink'
            with timed_step(audit,'plink',generation=g):
                det=call_plink(staged_vcf,prefix,chrom,cfg)
            with timed_step(audit,'score',generation=g):
                bed=Path(str(prefix)+'.blocks.le10Mb.bed')
                kept,large,lengths=plink_blocks_to_bed(det,chrom,bed,cfg.getint('analysis','max_block_bp'))
                score_blockset(cfg,chrom,bed,stage/f'score_G{g}',source=f'neutral_rep{rep}_G{g}')
                validate_score(stage/f'score_G{g}.npz',cfg.getint('random_aglink_null','n_random_sets'))
                arr=np.asarray(lengths,dtype=float)
                header='chrom\treplicate\tgeneration\tn_blocks_le10Mb\tn_blocks_gt10Mb_excluded\tmedian_block_bp\tmean_block_bp\tq95_block_bp\ttotal_block_bp\n'
                data=(f'{chrom}\t{rep}\t{g}\t{kept}\t{large}\t'
                      f'{(np.median(arr) if len(arr) else 0):.10g}\t'
                      f'{(np.mean(arr) if len(arr) else 0):.10g}\t'
                      f'{(np.quantile(arr,.95) if len(arr) else 0):.10g}\t'
                      f'{(np.sum(arr) if len(arr) else 0):.10g}\n')
                atomic_text(stage/'block_architecture.tsv',header+data)
            if cfg.getboolean('output','keep_checkpoint_vcf',fallback=False) and staged_vcf.is_symlink():
                # A retained output must be a real shared file, never a dangling
                # reference to a node-local directory removed after completion.
                materialized=stage/f'G{g}.vcf.materialized'
                shutil.copy2(staged_vcf,materialized)
                os.replace(materialized,staged_vcf)
            ensure_dir(final)
            for item in stage.iterdir():
                if item.is_file(): os.replace(item,final/item.name)
            stage.rmdir()
            # Publish DONE last; summaries never see partially published results.
            atomic_json(final/'DONE',dict(signature=signature,generation=g,completed_at=utc_now()))
            cleanup_checkpoint(cfg,final,g)
            if local_scratch: (vcfdir/f'G{g}.vcf').unlink(missing_ok=True)
            print(f'Completed chromosome={chrom} replicate={rep} generation={g}: {final}',flush=True)
        for g in generations:
            if not complete(finals[g],signature): raise RuntimeError(f'Incomplete checkpoint {g}')
            # A deterministic regeneration after a partial failure may have
            # re-emitted VCFs for checkpoints that were already scored.
            if not cfg.getboolean('output','keep_checkpoint_vcf',fallback=False):
                (vcfdir/f'G{g}.vcf').unlink(missing_ok=True)
        if local_scratch: remove_private_scratch(vcfdir,trajectory,signature,scratch_host)
        atomic_json(trajectory/'DONE',dict(signature=signature,finished_at=utc_now(),checkpoints=generations))
        print(f'Completed trajectory: {trajectory}',flush=True)



# Original functions: summarize_results.py

def read_task_table(path):
    with open(path) as fh: return list(csv.DictReader(fh,delimiter='\t'))


def q(v,p): return float(np.quantile(np.asarray(v,dtype=float),p)) if len(v) else float('nan')


def sample_sd(v): return float(np.std(np.asarray(v,dtype=float),ddof=1)) if len(v)>1 else (0.0 if len(v)==1 else float('nan'))


def load_score(path):
    z=np.load(path,allow_pickle=False)
    return dict(groups=z['group_names'].tolist(),strata=z['strata_names'].tolist(),real_all=z['real_all_same'],n_all=z['n_all_edges'],real=z['real_same'],n=z['n_edges'],rand=z['random_same'])



# Original pooled counts and statistics.
def summarize(cfg):
    work=Path(cfg['paths']['work_dir']); prep=work/'02_prepared'; out=work/'04_summary'; ensure_dir(out)
    chroms=csv_list(cfg.get('analysis','run_chromosomes',fallback=cfg['analysis']['chromosomes'])); generations=int_list(cfg['simulation']['checkpoints']); B=cfg.getint('random_aglink_null','n_random_sets')
    task_table=cfg.get('simulation','task_table_filename',fallback='neutral_task_table.tsv')
    tasks=read_task_table(prep/task_table); expected=len(tasks); done=0; scores={}; block_rows=[]; ld_rows=[]; maf_rows=[]
    for row in tasks:
        c=row['chrom']; rep=int(row['replicate']); gen=int(row['generation']); d=work/'03_neutral'/f'chr{c}_rep{rep:03d}_G{gen:03d}'
        if not (d/'DONE').exists(): continue
        done+=1; scores[(c,rep,gen)]=load_score(d/f'score_G{gen}.npz')
        with open(d/'block_architecture.tsv') as fh: block_rows += list(csv.DictReader(fh,delimiter='\t'))
        p=d/'ld_qc.tsv'
        if p.exists():
            with open(p) as fh:
                for x in csv.DictReader(fh,delimiter='\t'): x.update(chrom=c,replicate=str(rep)); ld_rows.append(x)
        p=d/'maf_qc.tsv'
        if p.exists():
            with open(p) as fh:
                for x in csv.DictReader(fh,delimiter='\t'): x.update(chrom=c,replicate=str(rep)); maf_rows.append(x)
    if done!=expected: raise SystemExit(f'Neutral tasks incomplete: {done}/{expected} DONE')
    # Observed real-block reference scored on exactly the same supported real/random agLink sets.
    obs={}
    for c in chroms:
        prefix=out/f'observed_{c}'
        score_blockset(cfg,c,cfg['paths']['real_blocks_bed'],prefix,source='observed_real_blocks')
        obs[c]=load_score(str(prefix)+'.npz')
    groups=obs[chroms[0]]['groups']; strata=obs[chroms[0]]['strata']; G=len(groups); S=len(strata)
    # Observed reference, chromosome-specific and pooled.
    with open(out/'observed_reference.tsv','w') as fh:
        fh.write('chrom\tgroup\tdistance_stratum\tn_all_edges\treal_all_fraction\tn_supported_edges\treal_supported_fraction\trandom_aglink_mean\trandom_q025\trandom_q975\tdelta_supported_real_minus_random\tempirical_p\n')
        for ckey in chroms+['POOLED']:
            if ckey=='POOLED':
                rra=sum(obs[c]['real_all'] for c in chroms); nna=sum(obs[c]['n_all'] for c in chroms); rr=sum(obs[c]['real'] for c in chroms); nn=sum(obs[c]['n'] for c in chroms); rand=sum(obs[c]['rand'] for c in chroms)
            else: rra,nna,rr,nn,rand=obs[ckey]['real_all'],obs[ckey]['n_all'],obs[ckey]['real'],obs[ckey]['n'],obs[ckey]['rand']
            for gi,g in enumerate(groups):
                for si,s in enumerate(strata):
                    na=int(nna[gi,si]); n=int(nn[gi,si])
                    if na==0: continue
                    real_all=float(rra[gi,si]/na)
                    if n>0:
                        realf=float(rr[gi,si]/n); rv=rand[gi,si,:]/n; p=(1+int(np.sum(rv>=realf)))/(B+1)
                        fh.write(f'{ckey}\t{g}\t{s}\t{na}\t{real_all:.10g}\t{n}\t{realf:.10g}\t{rv.mean():.10g}\t{np.quantile(rv,.025):.10g}\t{np.quantile(rv,.975):.10g}\t{(realf-rv.mean()):.10g}\t{p:.10g}\n')
                    else:
                        fh.write(f'{ckey}\t{g}\t{s}\t{na}\t{real_all:.10g}\t0\tnan\tnan\tnan\tnan\tnan\tnan\n')
    # Neutral forward summary. Pooled combines the same replicate and random-set indices across chromosomes.
    with open(out/'neutral_forward_summary.tsv','w') as fh:
        fh.write('generation\tchrom\tgroup\tdistance_stratum\tn_all_edges\tn_supported_edges\tn_sim_replicates\treal_all_mean\treal_all_q025\treal_all_q975\treal_supported_mean\treal_supported_q025\treal_supported_q975\trandom_aglink_mean\trandom_repmean_q025\trandom_repmean_q975\trandom_set_q025\trandom_set_q975\tdelta_supported_real_minus_random_mean\tdelta_q025\tdelta_q975\tempirical_p_random_ge_real\treal_supported_sd\trandom_repmean_sd\n')
        reps=range(1,cfg.getint('simulation','replicates')+1)
        for gen in generations:
            for ckey in chroms+['POOLED']:
                for gi,g in enumerate(groups):
                    for si,s in enumerate(strata):
                        realallv=[]; realv=[]; randmean=[]; deltav=[]; allrand=[]; ge=0; total=0; n_ref=None; na_ref=None
                        for rep in reps:
                            if ckey=='POOLED':
                                ss=[scores[(c,rep,gen)] for c in chroms]
                                rra=sum(x['real_all'][gi,si] for x in ss); nna=sum(x['n_all'][gi,si] for x in ss); rr=sum(x['real'][gi,si] for x in ss); nn=sum(x['n'][gi,si] for x in ss); rand=sum(x['rand'][gi,si,:] for x in ss)
                            else:
                                x=scores[(ckey,rep,gen)]; rra=x['real_all'][gi,si]; nna=x['n_all'][gi,si]; rr=x['real'][gi,si]; nn=x['n'][gi,si]; rand=x['rand'][gi,si,:]
                            na=int(nna); n=int(nn)
                            if na>0:
                                na_ref=na if na_ref is None else na_ref; realallv.append(float(rra/na))
                            if n==0: continue
                            n_ref=n if n_ref is None else n_ref
                            rf=float(rr/n); rv=rand/n; rm=float(rv.mean())
                            realv.append(rf); randmean.append(rm); deltav.append(rf-rm); allrand.extend(rv.tolist()); ge+=int(np.sum(rv>=rf)); total+=len(rv)
                        if not realallv: continue
                        p=(1+ge)/(1+total) if total else float('nan')
                        fh.write(f'{gen}\t{ckey}\t{g}\t{s}\t{na_ref or 0}\t{n_ref or 0}\t{len(realallv)}\t{np.mean(realallv):.10g}\t{q(realallv,.025):.10g}\t{q(realallv,.975):.10g}\t{(np.mean(realv) if realv else float("nan")):.10g}\t{q(realv,.025):.10g}\t{q(realv,.975):.10g}\t{(np.mean(randmean) if randmean else float("nan")):.10g}\t{q(randmean,.025):.10g}\t{q(randmean,.975):.10g}\t{q(allrand,.025):.10g}\t{q(allrand,.975):.10g}\t{(np.mean(deltav) if deltav else float("nan")):.10g}\t{q(deltav,.025):.10g}\t{q(deltav,.975):.10g}\t{p:.10g}\t{sample_sd(realv):.10g}\t{sample_sd(randmean):.10g}\n')
    # Block architecture summary.
    with open(out/'block_architecture_summary.tsv','w') as fh:
        fh.write('generation\tchrom\tn_reps\tmean_n_blocks_le10Mb\tmean_n_blocks_gt10Mb_excluded\tmedian_of_median_block_bp\tmean_total_block_bp\n')
        for gen in generations:
            for ckey in chroms+['POOLED']:
                vals=[x for x in block_rows if int(x['generation'])==gen and (ckey=='POOLED' or x['chrom']==ckey)]
                if not vals: continue
                if ckey=='POOLED':
                    # Aggregate per replicate first.
                    byrep=defaultdict(list)
                    for x in vals: byrep[int(x['replicate'])].append(x)
                    nb=[]; nl=[]; med=[]; tot=[]
                    for rep,zs in byrep.items():
                        nb.append(sum(int(z['n_blocks_le10Mb']) for z in zs)); nl.append(sum(int(z['n_blocks_gt10Mb_excluded']) for z in zs)); med.append(float(np.median([float(z['median_block_bp']) for z in zs]))); tot.append(sum(float(z['total_block_bp']) for z in zs))
                else:
                    nb=[int(z['n_blocks_le10Mb']) for z in vals]; nl=[int(z['n_blocks_gt10Mb_excluded']) for z in vals]; med=[float(z['median_block_bp']) for z in vals]; tot=[float(z['total_block_bp']) for z in vals]
                fh.write(f'{gen}\t{ckey}\t{len(nb)}\t{np.mean(nb):.10g}\t{np.mean(nl):.10g}\t{np.median(med):.10g}\t{np.mean(tot):.10g}\n')
    # Null-design support audit.
    with open(out/'random_null_support.tsv','w') as fh:
        fh.write('chrom\tgroup\tn_edges\tn_supported\tsupport_fraction\tmedian_candidates_supported\tmin_candidates_supported\n')
        for c in chroms:
            z=np.load(prep/c/'random_null_design.npz',allow_pickle=False); support=z['support'].astype(bool); cnt=z['candidate_count']; mem={g:z['group_'+g.replace('.','_')].astype(bool) for g in groups}
            for g in groups:
                m=mem[g]; n=int(m.sum()); ms=m&support; ns=int(ms.sum()); cc=cnt[ms]
                fh.write(f'{c}\t{g}\t{n}\t{ns}\t{(ns/n if n else float("nan")):.10g}\t{(np.median(cc) if len(cc) else float("nan")):.10g}\t{(np.min(cc) if len(cc) else float("nan")):.10g}\n')
    # LD QC summary: compare simulated r2 trajectory with observed r2 for the same sampled pairs.
    if cfg.getboolean('ld_qc','enabled',fallback=True):
        real_ld={}
        for c in chroms:
            p=prep/c/'real_ld.tsv'
            if p.exists():
                with open(p) as fh:
                    for x in csv.DictReader(fh,delimiter='\t'): real_ld[(c,int(x['bin_index']))]=float(x['mean_r2'])
        with open(out/'ld_qc_summary.tsv','w') as fh:
            fh.write('generation\tchrom\tbin_index\treal_mean_r2\tsim_mean_r2\tsim_q025\tsim_q975\tn_reps\n')
            for gen in generations:
                for c in chroms:
                    bins=sorted({int(x['bin_index']) for x in ld_rows if x['chrom']==c and int(x['generation'])==gen})
                    for b in bins:
                        v=[float(x['mean_r2']) for x in ld_rows if x['chrom']==c and int(x['generation'])==gen and int(x['bin_index'])==b and x['mean_r2'].lower()!='nan']
                        if not v: continue
                        fh.write(f'{gen}\t{c}\t{b}\t{real_ld.get((c,b),float("nan")):.10g}\t{np.mean(v):.10g}\t{q(v,.025):.10g}\t{q(v,.975):.10g}\t{len(v)}\n')
    with open(out/'maf_qc_summary.tsv','w') as fh:
        fh.write('generation\tchrom\taf_rmse_mean\tfixed_fraction_mean\tmean_het_target\tmean_het_simulated\tn_reps\n')
        for gen in generations:
            for c in chroms:
                v=[x for x in maf_rows if x['chrom']==c and int(x['generation'])==gen]
                if not v: continue
                fh.write(f'{gen}\t{c}\t{np.mean([float(x["af_rmse"]) for x in v]):.10g}\t{np.mean([float(x["fixed_fraction"]) for x in v]):.10g}\t{np.mean([float(x["mean_het_G0_target"]) for x in v]):.10g}\t{np.mean([float(x["mean_het_simulated"]) for x in v]):.10g}\t{len(v)}\n')
    print(out/'neutral_forward_summary.tsv')


# Original functions: prepare_plot_data.py

def read_tsv(path):
    with open(path) as handle:
        return list(csv.DictReader(handle,delimiter='\t'))


def extract(rows,chrom,generation,strata,expected_replicates=None):
    relevant=[row for row in rows if row['chrom']==chrom
              and row['group']=='all_significant' and int(row['generation'])==generation]
    indexed={row['distance_stratum']:row for row in relevant}
    if len(indexed)!=len(relevant):
        raise RuntimeError(f'Duplicate summary strata for {chrom} G{generation}')
    values=[]
    for label,low,high,mid_kb in strata:
        row=indexed.get(label)
        n=int(row['n_supported_edges']) if row is not None else 0
        n_all=int(row['n_all_edges']) if row is not None else 0
        if not 0<=n<=n_all: raise RuntimeError(f'Invalid denominator for {chrom} G{generation} {label}')
        value=dict(chrom=chrom,generation=generation,distance_stratum=label,
                   distance_low_bp=low,distance_high_bp=high,distance_midpoint_kb=mid_kb,
                   n_supported=n,n_all=n_all,
                   support_status=('supported' if n else 'no_background_matched_pairs' if n_all else 'no_eligible_real_pairs'))
        columns={'real':'real_supported_mean','real_sd':'real_supported_sd',
                 'real_low':'real_supported_q025','real_high':'real_supported_q975',
                 'random':'random_aglink_mean','random_sd':'random_repmean_sd',
                 'random_low':'random_repmean_q025','random_high':'random_repmean_q975'}
        for key,column in columns.items():
            value[key]=100*float(row[column]) if n else float('nan')
        if n:
            if not np.all(np.isfinite([value[k] for k in columns])):
                raise RuntimeError(f'Non-finite mean, SD or quantile for {chrom} G{generation} {label}')
            if expected_replicates is not None and int(row['n_sim_replicates'])!=expected_replicates:
                raise RuntimeError(f'Incomplete replicate count for {chrom} G{generation} {label}')
            if any(value[k]<0 for k in ('real_sd','random_sd')):
                raise RuntimeError(f'Negative SD for {chrom} G{generation} {label}')
            if any(not 0<=value[k]<=100 for k in ('real','random','real_low','real_high','random_low','random_high')):
                raise RuntimeError(f'Percentage outside 0--100 for {chrom} G{generation} {label}')
        values.append(value)
    return values


def write_source(path,selected):
    fields=['chrom','generation','distance_stratum','distance_low_bp','distance_high_bp',
            'distance_midpoint_kb','n_supported','real','real_sd','real_low','real_high',
            'random','random_sd','random_low','random_high','n_all','support_status']
    with open(path,'w',newline='') as out:
        writer=csv.DictWriter(out,fieldnames=fields,delimiter='\t')
        writer.writerow(dict(zip(fields,['chrom','generation','distance_stratum','distance_low_bp',
            'distance_high_bp','distance_midpoint_kb','n_supported_real_aglinks',
            'real_mean_percent','real_sd_percent','real_q025_percent','real_q975_percent',
            'random_mean_percent','random_sd_percent','random_q025_percent','random_q975_percent',
            'n_all_real_aglinks','support_status'])))
        for chrom in selected:
            for generation in (20,50,100,200):
                for row in selected[chrom][generation]:
                    writer.writerow({key:(('NA' if not np.isfinite(row[key]) else f'{row[key]:.10g}')
                                          if isinstance(row[key],float) else row[key]) for key in fields})



def prepare_plot_data(cfg):
    work=Path(cfg['paths']['work_dir'])
    rows=read_tsv(work/'04_summary'/'neutral_forward_summary.tsv')
    chromosomes=csv_list(cfg.get('analysis','run_chromosomes',fallback=cfg['analysis']['chromosomes']))
    bins=int_list(cfg['analysis']['distance_bins_bp'])
    if bins!=[0,10000,20000,50000,100000,250000,500000,750000,1000000]:
        raise RuntimeError('Expected the eight predefined distance strata')
    n_reps=cfg.getint('simulation','replicates')
    strata=[(distance_bin_label((low+high)//2,bins),low,high,(low+high)/2000.0)
            for low,high in zip(bins[:-1],bins[1:])]
    outdir=work/'05_figures'/'source_data'; ensure_dir(outdir)
    selected={}
    for chrom in chromosomes:
        selected[chrom]={generation:extract(rows,chrom,generation,strata,n_reps)
                         for generation in (20,50,100,200)}
    write_source(outdir/'all_chromosome_distance_source_data.tsv',selected)
    pooled={f'ALL_{len(chromosomes)}':{generation:extract(rows,'POOLED',generation,strata,n_reps)
                      for generation in (20,50,100,200)}}
    write_source(outdir/'pooled_distance_source_data.tsv',pooled)
    print(f'Created final plotting source tables for {len(chromosomes)} chromosomes and the pooled result in {outdir}')



# Original functions: finalize_outputs.py

def required_file(path, allow_empty=False):
    path=Path(path)
    if not path.is_file() or (not allow_empty and path.stat().st_size==0):
        raise RuntimeError(f'Missing or empty required output: {path}')
    return path


def assert_close(actual,expected,description,atol=2e-8):
    actual=float(actual); expected=float(expected)
    if np.isnan(expected):
        if not np.isnan(actual): raise RuntimeError(f'Expected NA: {description}, got {actual}')
    elif not np.isfinite(actual) or not np.isclose(actual,expected,rtol=1e-8,atol=atol):
        raise RuntimeError(f'Mismatched {description}: {actual} versus {expected}')


def fractions(numerators,denominators):
    numerators=np.asarray(numerators,dtype=float)
    denominators=np.asarray(denominators,dtype=float)
    result=np.full(np.broadcast_shapes(numerators.shape,denominators.shape),np.nan)
    np.divide(numerators,denominators,out=result,where=denominators>0)
    return result


def verify_statistics(cfg,work):
    if cfg.getint('simulation','replicates')!=30:
        raise RuntimeError('This run requires exactly 30 simulation replicates per chromosome')
    if cfg.getint('random_aglink_null','n_random_sets')!=1000:
        raise RuntimeError('This run requires exactly 1000 fixed random control sets')
    chroms=csv_list(cfg.get('analysis','run_chromosomes',fallback=cfg['analysis']['chromosomes']))
    if not chroms or len(set(chroms))!=len(chroms) or not set(chroms)<=set(CHROMOSOMES):
        raise RuntimeError('Invalid chromosome list')
    generations=int_list(cfg['simulation']['checkpoints'])
    if generations!=[10,20,50,100,200]: raise RuntimeError('Original five checkpoints must be retained')
    bins=int_list(cfg['analysis']['distance_bins_bp'])
    if bins!=DISTANCE_BINS: raise RuntimeError('Expected the original eight distance strata')
    strata=['ALL']+[distance_bin_label((a+b)//2,bins) for a,b in zip(bins[:-1],bins[1:])]
    prep=work/'02_prepared'
    tasks=read_tsv(prep/cfg.get('simulation','task_table_filename',fallback='neutral_task_table.tsv'))
    expected={(c,r,g) for c in chroms for r in range(1,31) for g in generations}
    task_keys=[(row['chrom'],int(row['replicate']),int(row['generation'])) for row in tasks]
    if len(tasks)!=len(expected) or len(set(task_keys))!=len(expected) or set(task_keys)!=expected:
        raise RuntimeError('Checkpoint task table does not contain all chromosome × 30 × 5 unique tasks')
    trajectories=read_tsv(prep/'trajectory_task_table.tsv')
    trajectory_keys={(row['chrom'],int(row['replicate'])):row for row in trajectories}
    if len(trajectories)!=len(chroms)*30 or set(trajectory_keys)!={(c,r) for c in chroms for r in range(1,31)}:
        raise RuntimeError('Trajectory table does not contain all chromosome × 30 unique trajectories')
    for row in tasks:
        if row['seed']!=trajectory_keys[(row['chrom'],int(row['replicate']))]['seed']:
            raise RuntimeError('Checkpoint and trajectory random seeds differ')
    design_counts={}
    for chrom in chroms:
        with np.load(prep/chrom/'random_null_design.npz',allow_pickle=False) as design:
            group=design['group_all_significant'].astype(bool)
            support=design['support'].astype(bool)
            distances=design['distance_bp']
            memberships=[np.ones(len(distances),dtype=bool)]+[
                (distances>=low)&(distances<high) for low,high in zip(bins[:-1],bins[1:])]
            n_all=np.asarray([int((group&mask).sum()) for mask in memberships],dtype=np.int64)
            n=np.asarray([int((group&support&mask).sum()) for mask in memberships],dtype=np.int64)
            random=np.load(prep/chrom/'random_choice_idx.npy',mmap_mode='r',allow_pickle=False)
            if random.shape!=(len(distances),1000):
                raise RuntimeError(f'Expected 1000 random sets per real pair: {chrom}')
            if np.any(random[support]<0): raise RuntimeError(f'Invalid supported random choice: {chrom}')
            design_counts[chrom]=(n_all,n)
    design_counts['POOLED']=(sum(design_counts[c][0] for c in chroms),sum(design_counts[c][1] for c in chroms))
    # Counts are first accumulated within a replicate across chromosomes.
    counts={}
    for chrom in chroms+['POOLED']:
        for generation in generations:
            counts[(chrom,generation)]=dict(real=np.zeros((30,9),dtype=np.int64),
                                           real_all=np.zeros((30,9),dtype=np.int64),
                                           random=np.zeros((30,9,1000),dtype=np.int64))
    for chrom,rep,generation in task_keys:
        directory=work/'03_neutral'/f'chr{chrom}_rep{rep:03d}_G{generation:03d}'
        required_file(directory/'DONE')
        required_file(directory/'block_architecture.tsv')
        required_file(directory/'ld_qc.tsv'); required_file(directory/'maf_qc.tsv')
        required_file(directory/f'G{generation}.plink.blocks.le10Mb.bed',allow_empty=True)
        with np.load(directory/f'score_G{generation}.npz',allow_pickle=False) as score:
            if score['group_names'].tolist()!=['all_significant'] or score['strata_names'].tolist()!=strata:
                raise RuntimeError(f'Unexpected score strata/groups: {directory}')
            n_all,n=design_counts[chrom]
            if not np.array_equal(score['n_all_edges'],n_all[None,:]) or not np.array_equal(score['n_edges'],n[None,:]):
                raise RuntimeError(f'Unequal or inconsistent real/random denominators: {directory}')
            if score['real_same'].shape!=(1,9) or score['real_all_same'].shape!=(1,9) or score['random_same'].shape!=(1,9,1000):
                raise RuntimeError(f'Unexpected score dimensions: {directory}')
            arrays=dict(real=score['real_same'][0],real_all=score['real_all_same'][0],random=score['random_same'][0])
            for name,array in arrays.items():
                denominator=n_all if name=='real_all' else n[:,None] if name=='random' else n
                if np.any(array<0) or np.any(array>denominator): raise RuntimeError(f'Invalid {name} counts: {directory}')
                counts[(chrom,generation)][name][rep-1]=array
                counts[('POOLED',generation)][name][rep-1]+=array
    summary_rows=read_tsv(work/'04_summary/neutral_forward_summary.tsv')
    summary={(row['chrom'],int(row['generation']),row['distance_stratum']):row for row in summary_rows}
    if len(summary)!=len(summary_rows): raise RuntimeError('Duplicate neutral summary rows')
    calculated={}; expected_summary=set()
    for (chrom,generation),arrays in counts.items():
        n_all,n=design_counts[chrom]
        real=fractions(arrays['real'],n[None,:])
        real_all=fractions(arrays['real_all'],n_all[None,:])
        random=fractions(arrays['random'].mean(axis=2),n[None,:])
        for j,stratum in enumerate(strata):
            key=(chrom,generation,stratum)
            metrics=dict(n_all=int(n_all[j]),n=int(n[j]),
                         real=float(real[:,j].mean()) if n[j] else np.nan,
                         real_sd=float(real[:,j].std(ddof=1)) if n[j] else np.nan,
                         random=float(random[:,j].mean()) if n[j] else np.nan,
                         random_sd=float(random[:,j].std(ddof=1)) if n[j] else np.nan)
            calculated[key]=metrics
            if n_all[j]==0: continue
            expected_summary.add(key)
            row=summary.get(key)
            if row is None: raise RuntimeError(f'Missing neutral summary: {key}')
            if int(row['n_all_edges'])!=n_all[j] or int(row['n_supported_edges'])!=n[j] or int(row['n_sim_replicates'])!=30:
                raise RuntimeError(f'Incorrect summary denominator/replicate count: {key}')
            for output,value in [('real_all_mean',real_all[:,j].mean()),
                                 ('real_supported_mean',metrics['real']),('real_supported_sd',metrics['real_sd']),
                                 ('random_aglink_mean',metrics['random']),('random_repmean_sd',metrics['random_sd']),
                                 ('delta_supported_real_minus_random_mean',metrics['real']-metrics['random'])]:
                assert_close(row[output],value,f'{key} {output}')
    if set(summary)!=expected_summary: raise RuntimeError('Unexpected extra neutral summary rows')
    return calculated,design_counts,dict(checkpoints_verified=len(tasks),trajectories=len(trajectories),
                                       replicates_per_chromosome=30,random_sets=1000,
                                       chromosomes=chroms,generations=generations,summary_rows=len(summary_rows))


def verify_plots(cfg,work,calculated):
    chromosomes=csv_list(cfg['analysis']['chromosomes'])
    pooled_label=f'ALL_{len(chromosomes)}'
    source=work/'05_figures/source_data'
    filenames=['all_chromosome_distance_source_data.tsv','pooled_distance_source_data.tsv']
    sizes=[len(chromosomes)*32,32]; all_rows=[]
    for filename,size in zip(filenames,sizes):
        rows=read_tsv(source/filename)
        if len(rows)!=size: raise RuntimeError(f'Expected {size} plot rows: {filename}, found {len(rows)}')
        all_rows.extend(rows)
    keys=set()
    for row in all_rows:
        chrom='POOLED' if row['chrom']==pooled_label else row['chrom']
        generation=int(row['generation']); stratum=row['distance_stratum']; key=(chrom,generation,stratum)
        if key in keys: raise RuntimeError(f'Duplicate plotting row: {key}')
        keys.add(key)
        if generation not in PLOTTED_GENERATIONS or stratum=='ALL' or key not in calculated:
            raise RuntimeError(f'Unexpected plotting row: {key}')
        expected=calculated[key]
        if int(row['n_supported_real_aglinks'])!=expected['n'] or int(row['n_all_real_aglinks'])!=expected['n_all']:
            raise RuntimeError(f'Incorrect plotting denominators: {key}')
        for column,metric in [('real_mean_percent','real'),('real_sd_percent','real_sd'),
                              ('random_mean_percent','random'),('random_sd_percent','random_sd')]:
            value=np.nan if row[column] in ('NA','NaN','nan','') else float(row[column])
            assert_close(value,100*expected[metric],f'{key} {column}',atol=2e-7)
        expected_status='supported' if expected['n'] else 'no_background_matched_pairs' if expected['n_all'] else 'no_eligible_real_pairs'
        if row['support_status']!=expected_status: raise RuntimeError(f'Incorrect support status: {key}')
    expected_keys={key for key in calculated if key[1] in PLOTTED_GENERATIONS and key[2]!='ALL'}
    if keys!=expected_keys: raise RuntimeError('Incomplete plot grid')
    figures=work/'05_figures/final'
    expected_pdfs=[f'Chromosome_{c}_distance_generations.pdf' for c in chromosomes]+[f'All_{len(chromosomes)}_chromosomes_pooled_distance_generations.pdf']
    for name in expected_pdfs:
        with required_file(figures/name).open('rb') as handle:
            if handle.read(5)!=b'%PDF-': raise RuntimeError(f'Invalid PDF: {figures/name}')
    limits=read_tsv(figures/'plot_limits_audit.tsv')
    if len(limits)!=(len(chromosomes)+1)*4 or len({(x['chrom'],int(x['generation'])) for x in limits})!=(len(chromosomes)+1)*4:
        raise RuntimeError('Expected an axis audit for all chromosome and pooled panels')
    expected_panels={(c,g) for c in chromosomes+[pooled_label] for g in PLOTTED_GENERATIONS}
    if {(x['chrom'],int(x['generation'])) for x in limits}!=expected_panels: raise RuntimeError('Unexpected panel axis audit')
    for row in limits:
        if row['error_bars_clipped'].upper()!='FALSE': raise RuntimeError('Error bars were clipped')
        if int(row['n_supported_strata'])+int(row['n_unsupported_strata'])!=8: raise RuntimeError('Missing categorical bins')
        if int(row['n_supported_strata']):
            if float(row['min_mean_minus_sd'])<float(row['y_min'])-1e-7 or float(row['max_mean_plus_sd'])>float(row['y_max'])+1e-7:
                raise RuntimeError('Panel limits do not include all mean ± SD values')
    return dict(chromosome_plot_rows=len(chromosomes)*32,pooled_plot_rows=32,pdfs=len(chromosomes)+1,panel_limits=(len(chromosomes)+1)*4)


# Package input staging and execution. Scientific routines above retain source logic.
def table_rows(path):
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def copy_input(source, destination):
    source=Path(source); destination=Path(destination)
    ensure_dir(destination.parent)
    if source.suffix=='.gz' and destination.suffix!='.gz':
        with gzip.open(source,'rb') as inp, destination.open('wb') as out:
            shutil.copyfileobj(inp,out)
    else:
        shutil.copy2(source,destination)
        if sha256(source)!=sha256(destination):
            raise RuntimeError(f'Input copy differs: {source}')


def input_spec(source, cfg):
    chromosomes=csv_list(cfg['analysis']['run_chromosomes'])
    spec=dict(chromosomes=chromosomes,edges=source/'02_prepared/evaluation_edges.tsv',
              genes=source/Path(cfg['input']['gff3']).expanduser(),
              blocks=source/Path(cfg['input']['real_blocks_bed']).expanduser(),chrom={})
    for chrom in chromosomes:
        d=source/'02_prepared'/chrom; b=source/'02_background'/chrom
        spec['chrom'][chrom]=dict(snps=d/'snps.tsv',length=d/'chromosome_length.txt',
            ld_pairs=d/'ld_pairs.tsv',real_ld=d/'real_ld.tsv',
            catalog=d/'random_null_design.npz',background=b/'acr_background.npz',map=b/'recombination_map.tsv')
    for key in ('chromosomes','run_chromosomes'):
        cfg.set('analysis',key,','.join(chromosomes))
    paths=[spec[key] for key in ('edges','genes','blocks')]
    paths.extend(path for files in spec['chrom'].values() for path in files.values())
    missing=[str(path) for path in paths if not path.is_file()]
    if missing:
        raise FileNotFoundError('Missing prepared input(s): ' + ', '.join(missing[:5]) +
            '. Set INPUT_ROOT to your prepared data directory and check [input] in settings.ini; or use demo mode.')
    return spec


def validate_parameters(cfg):
    if (cfg.getint('simulation','replicates'),cfg.getint('random_aglink_null','n_random_sets'),
            cfg.getint('simulation','population_size'),cfg.getint('simulation','output_sample_size'))!=(30,1000,203,203):
        raise ValueError('This package uses 30 replicates, 1000 background sets and 203 diploid individuals/samples.')
    if int_list(cfg['simulation']['checkpoints'])!=[10,20,50,100,200]:
        raise ValueError('Retain the original five generation checkpoints.')
    if cfg.getint('simulation','threads')<1:
        raise ValueError('threads must be positive')


def prepare_inputs(cfg, spec):
    work=Path(cfg['paths']['work_dir']); prep=work/'02_prepared'; ensure_dir(prep)
    chromosomes=spec['chromosomes']
    edges=table_rows(spec['edges'])
    edges=[row for row in edges if row['chrom'] in chromosomes]
    numeric=('edge_index','gene_tss0','acr_mid0','distance_bp','input_order','group_all_significant')
    known=defaultdict(set); per_chrom=defaultdict(list)
    for row in edges:
        for name in numeric: row[name]=int(row[name])
        known[row['gene_id']].add(row['acr_id']); per_chrom[row['chrom']].append(row)
    if len({row['edge_index'] for row in edges})!=len(edges):
        raise ValueError('Duplicate input edge indices')
    genes,_=read_genes(spec['genes'],set(chromosomes))
    for row in edges:
        if genes[row['gene_id']][4]!=row['gene_tss0']:
            raise ValueError('Gene TSS differs from frozen evaluation edges')
    copy_input(spec['blocks'],prep/'observed_blocks.bed')
    cfg.set('paths','real_blocks_bed',str(prep/'observed_blocks.bed'))
    audits={}
    for chrom in chromosomes:
        inputs=spec['chrom'][chrom]; d=prep/chrom; bg=work/'02_background'/chrom
        ensure_dir(d); ensure_dir(bg)
        for key,name in [('snps','snps.tsv'),('length','chromosome_length.txt'),
                         ('ld_pairs','ld_pairs.tsv'),('real_ld','real_ld.tsv')]:
            copy_input(inputs[key],d/name)
        copy_input(inputs['background'],bg/'acr_background.npz')
        copy_input(inputs['map'],bg/'recombination_map.tsv')
        with np.load(inputs['catalog'],allow_pickle=False) as z:
            acrs=[(str(aid),chrom,int(start),int(end),int(mid)) for aid,start,end,mid in
                  zip(z['acr_ids'],z['acr_start0'],z['acr_end0'],z['acr_mid0'])]
        if any(a[4]>b[4] for a,b in zip(acrs,acrs[1:])):
            raise ValueError('ACR catalogue is not sorted')
        if len({row[0] for row in acrs})!=len(acrs):
            raise ValueError('Duplicate catalogue ACR IDs')
        snps=table_rows(d/'snps.tsv'); length=int((d/'chromosome_length.txt').read_text())
        positions=[int(row['pos1']) for row in snps]
        if not positions or any(a>b for a,b in zip(positions,positions[1:])):
            raise ValueError(f'Empty or unsorted SNP input: {chrom}')
        if any(int(row['snp_index'])!=i or int(row['pos0'])!=positions[i]-1 or not 0<float(row['alt_freq'])<1 for i,row in enumerate(snps)):
            raise ValueError(f'Invalid SNP index, coordinate or frequency: {chrom}')
        if min(positions)<1 or max(positions)>length:
            raise ValueError('SNPs outside chromosome bounds')
        for row in table_rows(d/'ld_pairs.tsv'):
            i,j=int(row['i']),int(row['j'])
            if not 0<=i<j<len(snps) or abs(positions[j]-positions[i])!=int(row['distance_bp']):
                raise ValueError('Invalid LD pair input')
        audit=make_null_design(chrom,per_chrom[chrom],acrs,genes,known,cfg,
                               d/'random_null_design.npz',d/'random_null_audit.tsv')
        audits[chrom]=dict(audit,n_snps=len(snps),n_catalogue_acrs=len(acrs),chromosome_length=length)
        print(f'Prepared {chrom}: {len(snps)} SNPs, {audit["n_supported"]}/{audit["n_edges"]} supported pairs, 1000 control sets',flush=True)
    with (prep/'evaluation_edges.tsv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fieldnames=list(edges[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader(); writer.writerows(edges)
    tasks=[]; trajectories=[]
    for rep in range(1,31):
        for chrom in chromosomes:
            seed=stable_seed(cfg.getint('simulation','base_seed'),chrom,rep)
            trajectories.append((len(trajectories),chrom,rep,seed))
            for generation in int_list(cfg['simulation']['checkpoints']):
                tasks.append((len(tasks),chrom,rep,generation,seed))
    for filename,header,rows in [
        ('trajectory_task_table.tsv','trajectory_id\tchrom\treplicate\tseed',trajectories),
        ('neutral_task_table.tsv','task_id\tchrom\treplicate\tgeneration\tseed',tasks)]:
        (prep/filename).write_text(header+'\n'+''.join('\t'.join(map(str,row))+'\n' for row in rows))
    write_json(prep/'prepare_audit.json',dict(chromosomes=audits,replicates=30,n_random_sets=1000,
               population_size=203,checkpoints=[10,20,50,100,200],n_trajectories=len(trajectories),n_checkpoints=len(tasks)))
    return len(trajectories)


def find_plink(work, source):
    requested=os.environ.get('PLINK_BIN')
    candidates=([shutil.which(requested) or requested] if requested else
                [work/'01_bin/plink1.9', source/'01_bin/plink1.9',
                 shutil.which('plink1.9'), shutil.which('plink')])
    for candidate in candidates:
        if not candidate: continue
        path=Path(candidate).resolve()
        if not path.is_file(): continue
        try:
            version=subprocess.run([str(path),'--version'],capture_output=True,text=True,check=True).stdout.strip()
        except (OSError,subprocess.CalledProcessError):
            continue
        if re.search(r'PLINK v1\.9',version): return path,version
    raise RuntimeError('PLINK 1.9 was not found. Set PLINK_BIN to a working executable.')


def compile_tools(work, package, plink):
    bindir=work/'01_bin'; ensure_dir(bindir)
    source=package/'02_forward_neutral.cpp'; executable=bindir/'forward_neutral'
    marker=bindir/'simulator_build.json'; built={}
    if marker.is_file():
        try: built=json.loads(marker.read_text())
        except ValueError: pass
    source_hash=sha256(source)
    valid=(executable.is_file() and os.access(executable,os.X_OK)
           and built.get('source_sha256')==source_hash
           and built.get('executable_sha256')==sha256(executable))
    if not valid:
        fd,name=tempfile.mkstemp(prefix='.neutral_build_',dir=bindir); os.close(fd)
        temporary=Path(name)
        cmd=[os.environ.get('CXX','g++'),'-O3','-std=c++11','-fopenmp',str(source),'-o',str(temporary)]
        try:
            with (bindir/'compile.log').open('w') as log:
                try: subprocess.run(cmd[:4]+['-static']+cmd[4:],check=True,stdout=log,stderr=subprocess.STDOUT)
                except subprocess.CalledProcessError:
                    subprocess.run(cmd,check=True,stdout=log,stderr=subprocess.STDOUT)
            with temporary.open('rb') as handle:
                if handle.read(4)!=b'\x7fELF': raise RuntimeError('Compiler did not produce a valid Linux executable')
            temporary.chmod(0o755)
            os.replace(temporary,executable)
            atomic_json(marker,dict(source_sha256=source_hash,executable_sha256=sha256(executable)))
        finally:
            temporary.unlink(missing_ok=True)
    dest=bindir/'plink1.9'
    # Keep a real copy: /opt may be absent when a run resumes on another node.
    expected=sha256(plink)
    if dest.is_symlink() or not dest.is_file() or sha256(dest)!=expected:
        fd,name=tempfile.mkstemp(prefix='.plink_copy_',dir=bindir); os.close(fd)
        temporary=Path(name)
        try:
            shutil.copy2(plink,temporary)
            if sha256(temporary)!=expected: raise RuntimeError('PLINK copy differs from selected executable')
            os.replace(temporary,dest)
        finally:
            temporary.unlink(missing_ok=True)


def trajectory_worker(config_path, trajectory_id):
    cfg=load_config(config_path); work=Path(cfg['paths']['work_dir'])
    with (work/'00_logs'/f'trajectory_{trajectory_id:04d}.log').open('a') as log:
        with contextlib.redirect_stdout(log): run_trajectory(cfg,trajectory_id)
    return trajectory_id


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--jobs',type=int,default=1)
    parser.add_argument('--demo',action='store_true',help='Use only the bundled chromosome 1A input subset')
    args=parser.parse_args(); package=Path(__file__).resolve().parent
    source=args.input.resolve(); work=args.output.resolve()
    if args.jobs<1: parser.error('--jobs must be positive')
    inside_package=package in work.parents
    if work in (source,package) or work in source.parents or (source in work.parents and not inside_package):
        parser.error('Use an output directory inside this package or outside the original source run.')
    cfg=load_config(package/'settings.ini'); validate_parameters(cfg)
    if args.demo:
        cfg.set('analysis','run_chromosomes','1A')
        cfg.set('input','gff3','genes.gff3')
        cfg.set('input','real_blocks_bed','observed_blocks.bed')
    try:
        spec=input_spec(source,cfg)
    except (OSError,KeyError,ValueError) as error:
        parser.error(str(error))
    try:
        plink,version=find_plink(work,source)
    except RuntimeError as error:
        parser.error(str(error))
    rscript=os.environ.get('RSCRIPT_BIN','Rscript')
    if not shutil.which(os.environ.get('CXX','g++')):
        parser.error('C++ compiler not found. Install g++ with OpenMP or set CXX.')
    try:
        subprocess.run([rscript,'-e','stopifnot(capabilities("cairo"))'],check=True)
    except (OSError,subprocess.CalledProcessError):
        parser.error('R with Cairo is required. Install R or set RSCRIPT_BIN to a working Rscript.')
    inputs={key:spec[key] for key in ['edges','genes','blocks']}
    for chrom,files in spec['chrom'].items():
        inputs.update({f'{chrom}/{key}':path for key,path in files.items()})
    metadata=dict(mode='demo' if args.demo else 'full',chromosomes=spec['chromosomes'],
        parameters={section:dict(cfg[section]) for section in cfg.sections()},
        input_sha256={name:sha256(path) for name,path in inputs.items()},
        code_sha256={name:sha256(package/name) for name in ['01_run_neutral.py','02_forward_neutral.cpp','03_plot_results.R']},
        plink_version=version,plink_sha256=sha256(plink))
    signature=hashlib.sha256(json.dumps(metadata,sort_keys=True).encode()).hexdigest()
    marker=work/'RUN.json'
    if marker.exists():
        if json.loads(marker.read_text()).get('signature')!=signature:
            parser.error('Existing output has different inputs, code or parameters; choose a new output directory.')
    elif work.exists() and any(work.iterdir()):
        parser.error('Output directory is not empty and has no matching RUN.json; choose a new directory.')
    ensure_dir(work); ensure_dir(work/'00_logs'); ensure_dir(work/'.cache')
    lock=(work/'RUN.lock').open('a')
    try: fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError: parser.error('This output directory is already running.')
    with lock:
        atomic_json(marker,dict(signature=signature,**metadata))
        cfg.add_section('paths')
        cfg.set('paths','work_dir',str(work)); cfg.set('paths','plink',str(plink))
        cfg.set('paths','real_blocks_bed',str(work/'02_prepared/observed_blocks.bed'))
        os.environ['XDG_CACHE_HOME']=str(work/'.cache')
        compile_tools(work,package,plink)
        prepared=work/'02_prepared/PREPARED.json'
        if not prepared.exists():
            n=prepare_inputs(cfg,spec)
            atomic_json(prepared,dict(signature=signature,trajectories=n))
        else:
            state=json.loads(prepared.read_text())
            if state['signature']!=signature: raise RuntimeError('Prepared inputs differ from current run')
            n=int(state['trajectories'])
        config_path=work/'settings.ini'
        with config_path.open('w') as handle: cfg.write(handle)
        print(f'Design: {len(spec["chromosomes"])} chromosomes x 30 replicates x 5 checkpoints; 1000 controls; N=203',flush=True)
        if args.jobs==1:
            for task in range(n):
                trajectory_worker(str(config_path),task)
                print(f'Completed trajectory {task+1}/{n}',flush=True)
        else:
            with ProcessPoolExecutor(max_workers=min(args.jobs,n)) as pool:
                futures=[pool.submit(trajectory_worker,str(config_path),task) for task in range(n)]
                for count,future in enumerate(as_completed(futures),1):
                    future.result(); print(f'Completed trajectories {count}/{n}',flush=True)
        summarize(cfg); prepare_plot_data(cfg)
        plots=work/'05_figures'
        subprocess.run([rscript,str(package/'03_plot_results.R'),
            str(plots/'source_data/all_chromosome_distance_source_data.tsv'),
            str(plots/'source_data/pooled_distance_source_data.tsv'),str(plots/'final')],check=True)
        calculated,_,validation=verify_statistics(cfg,work)
        validation.update(verify_plots(cfg,work,calculated))
        validation.update(status='PASS',population_size=203,output_sample_size=203)
        atomic_json(work/'VALIDATION.json',validation)
        print('PASS:',json.dumps(validation,sort_keys=True)); print('Output:',work)


if __name__=='__main__':
    main()
