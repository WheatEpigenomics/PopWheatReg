"""Independent nominal cis scan; old global/subgroup implementations stay untouched."""
import csv
import gzip
import json
import os
from pathlib import Path
import re
import time

import numpy as np
import pandas as pd

FIELDS = ['variant_id','phenotype_id','pval','beta','beta_se','af','n_samples',
          'ref','alt','effect_allele','other_allele','allele_source',
          'chrom','variant_pos','phenotype_center','distance','df','n_genotyped','fit_status']
QC_FIELDS = ['phenotype_id','chrom','center','candidate_variants','output_pairs','status']
LEGACY_FIELDS = ['phenotype_id','variant_id','start_distance','end_distance','af',
                 'ma_samples','ma_count','pval_nominal','slope','slope_se']


def identity(path):
    p=Path(path).resolve();s=p.stat()
    return dict(path=str(p),size=s.st_size,mtime_ns=s.st_mtime_ns)


def canonical(values):
    raw=pd.Series(list(map(str,values)))
    normalized=raw.str.replace(r'^(?i:chr)','',regex=True)
    if normalized.eq('').any() or not normalized.str.fullmatch(r'[A-Za-z0-9_.-]+').all():
        raise ValueError('Invalid chromosome name')
    mapping=pd.DataFrame(dict(raw=raw,chrom=normalized)).drop_duplicates()
    if mapping['chrom'].duplicated().any():
        raise ValueError('Ambiguous chromosome aliases within one input (e.g. 1A and chr1A)')
    return normalized.to_numpy()


def positive_integer(values,name,minimum=1):
    a=pd.to_numeric(values,errors='raise').to_numpy(dtype=np.float64)
    if not np.isfinite(a).all() or (a<minimum).any() or (a!=np.floor(a)).any():
        raise ValueError(name+' must contain finite integer coordinates')
    return a.astype(np.int64)


def positions(args,phenotype,core):
    """Use explicit BED coordinates; never guess gene strand from an ID."""
    ids=phenotype.index
    annotation=args.phenotype_positions
    if annotation:
        table=pd.read_csv(annotation,sep='\t',comment='#',header=None,dtype=str)
        if table.shape[1]<4:raise ValueError('--phenotype_positions needs BED4 or BED6')
        table.index=table[3];core._check_unique(table.index,'position phenotype IDs')
        missing=ids.difference(table.index)
        if len(missing):raise ValueError('Phenotypes absent from position BED: '+core._preview(missing))
        table=table.loc[ids]
        chrom=table[0];start=positive_integer(table[1],'BED start',0);end=positive_integer(table[2],'BED end')
        strand=table[5].to_numpy() if table.shape[1]>=6 else None
        source=str(Path(annotation).resolve())
    elif str(args.phenotype).lower().endswith(('.bed','.bed.gz')):
        table=pd.read_csv(args.phenotype,sep='\t',usecols=[0,1,2,3],dtype=str)
        table.index=table.iloc[:,3];core._check_unique(table.index,'BED phenotype IDs');table=table.loc[ids]
        chrom=table.iloc[:,0];start=positive_integer(table.iloc[:,1],'BED start',0);end=positive_integer(table.iloc[:,2],'BED end')
        strand=None;source='phenotype BED coordinates'
    elif args.phenotype_kind=='caqtl':
        matches=[re.fullmatch(r'[^:\s]+:([^:\s]+)_(\d+)_(\d+)',str(g)) for g in ids]
        if not all(matches):raise ValueError('TSV/parquet caQTL needs --phenotype_positions BED4, or ACR IDs encoding chr_start_end')
        chrom=[m[1] for m in matches];start=np.array([int(m[2]) for m in matches]);end=np.array([int(m[3]) for m in matches])
        strand=None;source='ACR IDs: BED0 half-open chr_start_end'
    else:
        raise ValueError('eQTL needs --phenotype_positions BED6 for strand-aware TSS (or a single-base phenotype BED)')
    if (end<=start).any():raise ValueError('BED coordinates require 0 <= start < end')
    if args.phenotype_kind=='eqtl':
        if strand is None:
            if not (end==start+1).all():raise ValueError('Gene intervals require BED6 strand information; refusing to guess TSS')
            center=end
        else:
            if not np.isin(strand,['+','-']).all():raise ValueError('Gene strand must be + or -')
            center=np.where(strand=='+',start+1,end)
    else:
        center=(start+end)//2+1
    return pd.DataFrame(dict(chrom=canonical(chrom),center=center,
                             start=start+1,end=end),index=ids),source


def variant_metadata(args,genotype,variant_pos,core):
    prefix=args.genotype
    if all(Path(prefix+'.'+s).exists() for s in ('pgen','pvar','psam')):
        from tensorqtl.pgen import read_pvar
        table=read_pvar(prefix+'.pvar').rename(columns={'id':'variant_id'}).set_index('variant_id')
        table=table[['chrom','pos','ref','alt']].copy()
        if table[['ref','alt']].isna().any().any():raise ValueError('PVAR has missing REF/ALT')
        # Pgenlib hardcall read() defaults to allele_idx=1: counts first ALT.
        table['effect_allele']=table['alt'].str.split(',').str[0]
        table['other_allele']=table['ref'].where(~table['alt'].str.contains(','),'NA')
        table['allele_source']='PVAR_first_ALT_count'
    elif all(Path(prefix+'.'+s).exists() for s in ('bed','bim','fam')):
        table=pd.read_csv(prefix+'.bim',sep=r'\s+',header=None,dtype=str)
        table.columns=['chrom','variant_id','cm','pos','effect_allele','other_allele'];table=table.set_index('variant_id')
        # Verified against VCF GT, PGEN and BED using installed TensorQTL 1.0.10.
        # BIM A1/A2 alone do not establish genomic REF/ALT.
        table['ref']='NA';table['alt']='NA';table['allele_source']='BIM_A1_count_REF_ALT_unknown'
    else:
        if not args.variant_metadata:
            raise ValueError('Matrix genotypes require --variant_metadata TSV: variant_id,chrom,pos,effect_allele,other_allele (REF/ALT optional)')
        table=pd.read_csv(args.variant_metadata,sep='\t',dtype=str).set_index('variant_id')
        required={'chrom','pos','effect_allele','other_allele'}
        if not required.issubset(table):raise ValueError('Variant metadata missing required columns: '+str(sorted(required-set(table))))
        for col in ('ref','alt'):
            if col not in table:table[col]='NA'
        table['allele_source']='explicit_matrix_counted_allele_metadata'
    core._check_unique(genotype.index,'genotype variant IDs');core._check_unique(table.index,'variant metadata IDs')
    missing=genotype.index.difference(table.index)
    if len(missing):raise ValueError('Variants lack position/allele metadata: '+core._preview(missing))
    table=table.loc[genotype.index].copy()
    table['chrom']=canonical(table['chrom']);table['pos']=positive_integer(table['pos'],'SNP position')
    for col in ('effect_allele','other_allele'):
        bad=table[col].isna()|table[col].isin(['','.','NA','0'])
        # A multiallelic PGEN has a known ALT1 count; non-effect alleles are a mixture.
        if col=='other_allele':bad &= ~table['alt'].fillna('').str.contains(',')
        if bad.any():raise ValueError('Missing/unreliable '+col+' in variant metadata')
    if variant_pos is not None:
        vp=variant_pos.loc[genotype.index]
        if not np.array_equal(canonical(vp.iloc[:,0]),table['chrom'].values) or not np.array_equal(
                positive_integer(vp.iloc[:,1],'loaded SNP position'),table['pos'].values):
            raise ValueError('Variant metadata differs from loaded genotype positions')
    return table


def run(args,core):
    from tensorqtl.core import get_allele_stats
    started=time.monotonic();out=Path(args.out).resolve()
    if not args.phenotype_kind:raise ValueError('cis_all requires --phenotype_kind eqtl|caqtl')
    if args.window<0:raise ValueError('--window must be >=0')
    if args.resume:raise ValueError('cis_all does not support --resume; use a new empty output directory')
    if not out.is_dir() or any(out.iterdir()):raise ValueError('cis_all --out must be an existing empty directory; existing results are preserved')
    core.validate_parameters(1.,args.maf_threshold,args.batch_size,args.pheno_chunk)
    paths=[args.phenotype]+([args.covariates] if args.covariates else [])
    paths += [p for p in (args.phenotype_positions,args.variant_metadata) if p]
    paths += ([args.genotype+'.'+s for s in ('pgen','pvar','psam')]
              if Path(args.genotype+'.pgen').exists() else
              [args.genotype+'.'+s for s in ('bed','bim','fam')] if Path(args.genotype+'.bim').exists() else [args.genotype])
    inputs=[identity(p) for p in paths]
    y=core.load_phenotypes(args.phenotype);core._check_unique(y.index,'phenotype IDs');core._check_unique(y.columns,'phenotype sample IDs')
    pos,position_source=positions(args,y,core)
    cov=core.load_covariates(args.covariates) if args.covariates else None
    genotype_samples=core.plink_sample_ids(args.genotype)
    if genotype_samples is not None:
        core._check_unique(genotype_samples,'genotype sample IDs')
        selected=y.columns.intersection(genotype_samples,sort=False).tolist()
        if not selected:raise ValueError('Genotype and phenotype have 0 common sample IDs')
    else:selected=y.columns.tolist()
    if args.genotype.endswith('.gz'):core.check_raw_sample_header(args.genotype,4 if args.genotype.endswith('.bed.gz') else 1,'genotype')
    g,vp=core.load_genotypes(args.genotype,select_samples=selected)
    core.validate_matrices(g,y,cov)
    common=y.columns.intersection(g.columns,sort=False);y=y.loc[:,common];g=g.loc[:,common]
    if cov is not None:
        missing=common.difference(cov.index)
        if len(missing):raise ValueError('Analysis samples missing from covariates: '+core._preview(missing))
        cov=cov.loc[common]
    core.validate_alignment(g,y,cov)
    n=len(common);q=0 if cov is None else cov.shape[1];df=n-2-q
    if df<=0:raise ValueError('cis_all requires n_samples - 2 - n_covariates > 0')
    design=np.ones((n,q+1),dtype=np.float64)
    if cov is not None:design[:,1:]=cov.values
    if np.linalg.matrix_rank(design)!=q+1:raise ValueError('Covariates plus intercept must have full rank')
    Q=np.linalg.qr(design,mode='reduced')[0]
    meta=variant_metadata(args,g,vp,core)
    torch=core.torch;device=core._get_device()
    residualizer=core.Residualizer(torch.tensor(cov.values,dtype=torch.float32,device=device)) if cov is not None else None
    logger=core.SimpleLogger(logfile=args.logfile,verbose=True)
    logger.write('cis_all: '+args.phenotype_kind+'; device='+str(device))
    logger.write(f'  * {n} samples, {q} covariates, df={df}; mean-imputed genotype counts')
    logger.write('  * all computable cis pairs; --pval_threshold ignored; no nominal P filter')
    logger.write('  * host input loading matches existing modes; only cis SNP batches reach regression')
    chrom_indexes={}
    for chrom,rows in meta.groupby('chrom',sort=False):
        ix=g.index.get_indexer(rows.index);order=np.argsort(rows['pos'].values,kind='stable')
        chrom_indexes[chrom]=(rows['pos'].values[order],ix[order])
    chromosomes=sorted(pos['chrom'].unique());qc_path=out/f'{len(chromosomes)+2:02d}_phenotype_qc.tsv'
    legacy_path=out/f'{len(chromosomes)+3:02d}_cis_nominal_legacy.tsv.gz'
    legacy_partial=Path(str(legacy_path)+'.partial');legacy_pairs=0
    counters=dict(candidate_pairs=0,output_pairs=0,invalid_genotype_pairs=0,maf_filtered_pairs=0,
                  nonestimable_pairs=0,invalid_statistic_pairs=0,zero_correlation_pairs=0,perfect_fit_pairs=0)
    outputs=[];qc_status={}
    with Path(str(qc_path)+'.partial').open('w') as qc_file, gzip.open(legacy_partial,'wt') as legacy_file:
        qc_writer=csv.writer(qc_file,delimiter='\t',lineterminator='\n');qc_writer.writerow(QC_FIELDS)
        for file_no,chrom in enumerate(chromosomes,1):
            target=out/f'{file_no:02d}_{chrom}.cis_all.tsv.gz';partial=Path(str(target)+'.partial');chrom_pairs=0
            subset=pos[pos['chrom']==chrom].sort_values('center',kind='stable')
            with gzip.open(partial,'wt') as f:
                f.write('\t'.join(FIELDS)+'\n')
                for pid,row in subset.iterrows():
                    center=int(row['center']);candidates=np.array([],dtype=np.int64)
                    if chrom in chrom_indexes:
                        pp,ii=chrom_indexes[chrom]
                        lo=np.searchsorted(pp,max(1,center-args.window),'left');hi=np.searchsorted(pp,center+args.window,'right')
                        candidates=ii[lo:hi]
                    counters['candidate_pairs']+=len(candidates);written=0
                    values=y.loc[pid].to_numpy(dtype=np.float64)
                    yr=values-Q@(Q.T@values)
                    y_estimable=np.linalg.norm(yr)>np.finfo(float).eps*max(design.shape)*max(1.,np.linalg.norm(values))*10
                    if y_estimable:
                        yt=torch.tensor(values[None,:],dtype=torch.float32,device=device)
                        for b in range(0,len(candidates),args.batch_size):
                            ix=candidates[b:b+args.batch_size]
                            raw=g.iloc[ix].to_numpy(dtype=np.float32)
                            observed=(raw!=-9)&np.isfinite(raw);observed_n=observed.sum(1)
                            valid=np.isfinite(raw).all(1)&(((raw>=0)&(raw<=2))|(raw==-9)).all(1)&(observed_n>0)
                            counters['invalid_genotype_pairs']+=int((~valid).sum());ix=ix[valid];raw=raw[valid];observed_n=observed_n[valid]
                            if not len(ix):continue
                            gt=torch.tensor(raw,dtype=torch.float32,device=device);core.impute_mean(gt)
                            ids=g.index[ix].to_numpy();gt,kept_ids,af=core.filter_maf(gt,ids,args.maf_threshold)
                            local=pd.Index(ids).get_indexer(kept_ids);observed_n=observed_n[local];ix=ix[local]
                            counters['maf_filtered_pairs']+=len(ids)-len(kept_ids)
                            if not len(ix):continue
                            # Float64 full-design rank check is QC only; statistics use the unchanged TensorQTL float32 core.
                            gd=gt.cpu().numpy().astype(np.float64);gr=gd-(gd@Q)@Q.T
                            estimable=np.linalg.norm(gr,axis=1)>np.finfo(float).eps*max(design.shape)*np.maximum(1.,np.linalg.norm(gd,axis=1))*10
                            counters['nonestimable_pairs']+=int((~estimable).sum())
                            r,vg,vy=core.calculate_corr(gt,yt,residualizer=residualizer,return_var=True)
                            r=r[:,0].double();ratio=torch.sqrt(vy[0]/vg).double()
                            eligible=torch.as_tensor(estimable,device=device)&torch.isfinite(r)&(r.abs()<=1+1e-6)&(vg>0)&(vy[0]>0)
                            r=r.clamp(-1,1);r2=r.square();t=r*torch.sqrt(df/(1-r2))
                            beta=(r*ratio).float();se=(ratio*torch.sqrt((1-r2)/df)).float()
                            pval=core.get_t_pval(t.cpu().numpy(),df,log=False)
                            beta_np=beta.cpu().numpy();se_np=se.cpu().numpy();af_np=af.cpu().numpy();r_np=r.cpu().numpy()
                            valid_stats=eligible.cpu().numpy()&np.isfinite(pval)&np.isfinite(beta_np)&np.isfinite(se_np)&(se_np>=0)
                            counters['invalid_statistic_pairs']+=int((estimable&~valid_stats).sum())
                            take=np.flatnonzero(valid_stats)
                            if len(take):
                                m=meta.iloc[ix[take]]
                                data=dict(variant_id=m.index.to_numpy(),phenotype_id=pid,pval=pval[take],beta=beta_np[take],
                                    beta_se=se_np[take],af=af_np[take],n_samples=n,
                                    ref=m['ref'].values,alt=m['alt'].values,effect_allele=m['effect_allele'].values,
                                    other_allele=m['other_allele'].values,allele_source=m['allele_source'].values,
                                    chrom=chrom,variant_pos=m['pos'].values,phenotype_center=center,distance=m['pos'].values-center,
                                    df=df,n_genotyped=observed_n[take],fit_status=np.where(np.abs(r_np[take])==1,'perfect_fit','regular'))
                                pd.DataFrame(data,columns=FIELDS).to_csv(f,sep='\t',header=False,index=False,float_format='%.17g',na_rep='NA')
                                # Native TensorQTL minor-allele counts, calculated from
                                # the same mean-imputed genotypes used in regression.
                                # Do not infer these integer fields from rounded AF.
                                _,ma_samples,ma_count=get_allele_stats(gt)
                                legacy_data=dict(phenotype_id=pid,variant_id=m.index.to_numpy(),
                                    start_distance=m['pos'].values-int(row.start),
                                    end_distance=m['pos'].values-int(row.end),af=af_np[take],
                                    ma_samples=ma_samples.cpu().numpy()[take],ma_count=ma_count.cpu().numpy()[take],
                                    pval_nominal=pval[take],slope=beta_np[take],slope_se=se_np[take])
                                pd.DataFrame(legacy_data,columns=LEGACY_FIELDS).to_csv(legacy_file,sep='\t',
                                    header=False,index=False,float_format='%.17g',na_rep='NA')
                                legacy_pairs+=len(take)
                                written+=len(take);counters['zero_correlation_pairs']+=int((r_np[take]==0).sum())
                                counters['perfect_fit_pairs']+=int((np.abs(r_np[take])==1).sum())
                    else:counters['nonestimable_pairs']+=len(candidates)
                    status='nonestimable_phenotype' if not y_estimable else 'no_cis_variants' if not len(candidates) else 'done' if written else 'no_computable_pairs'
                    qc_writer.writerow([pid,chrom,center,len(candidates),written,status]);qc_status[status]=qc_status.get(status,0)+1
                    counters['output_pairs']+=written;chrom_pairs+=written
            os.replace(partial,target);outputs.append(dict(**identity(target),pairs=chrom_pairs))
            logger.write(f'  * {chrom}: {len(subset)} phenotypes, {chrom_pairs:,} cis pairs -> {target.name}')
    os.replace(str(qc_path)+'.partial',qc_path)
    if legacy_pairs!=counters['output_pairs']:
        raise RuntimeError('Legacy output pair count mismatch; no completion manifest written')
    os.replace(legacy_partial,legacy_path)
    if inputs!=[identity(p) for p in paths]:raise ValueError('Inputs changed during cis scan; no completion manifest written')
    accounted=sum(counters[k] for k in ('output_pairs','invalid_genotype_pairs','maf_filtered_pairs',
                                       'nonestimable_pairs','invalid_statistic_pairs'))
    if counters['candidate_pairs']!=accounted:
        raise RuntimeError('Internal cis QC accounting mismatch; no completion manifest written')
    manifest=dict(completed=True,mode='cis_all',phenotype_kind=args.phenotype_kind,window=args.window,
        center_source=position_source,coordinate_policy='BED0 half-open; SNP and centers 1-based; +TSS=start+1, -TSS=end; ACR center=floor((start+end)/2)+1; inclusive window',
        maf_threshold=args.maf_threshold,pval_threshold=None,n_samples=n,df=df,sample_ids=common.tolist(),
        n_phenotypes=len(y),n_variants_loaded=len(g),phenotype_qc_status=qc_status,counters=counters,
        genotype_missing='-9 mean imputation within actual analysis samples; no additional missing-rate cutoff',
        se='same residualized variances as global_pc; direct SE formula fixes beta=t=0 only in cis_all',
        host_genotypes='same full host genotype loader as existing modes; cis candidate restriction before regression',
        covariance_rank='full rank with intercept required; nonestimable phenotype/genotype pairs excluded',
        inputs=inputs,outputs=outputs,phenotype_qc=identity(qc_path),
        legacy_output=dict(**identity(legacy_path),pairs=legacy_pairs,columns=LEGACY_FIELDS,header=False,
            distance_policy='SNP_POS - BED_start_1based; SNP_POS - BED_end_1based; no strand flip',
            counts='native tensorqtl.core.get_allele_stats on mean-imputed analysis genotypes'),
        elapsed_seconds=time.monotonic()-started)
    path=out/f'{len(chromosomes)+1:02d}_cis_all_manifest.json';tmp=Path(str(path)+'.partial')
    tmp.write_text(json.dumps(manifest,indent=2)+'\n');os.replace(tmp,path)
    logger.write(f'cis_all complete: {counters["output_pairs"]:,} pairs; {path}')
    logger.write(f'  * legacy nominal: {legacy_pairs:,} pairs, 10 columns without header -> {legacy_path.name}')
    return manifest
