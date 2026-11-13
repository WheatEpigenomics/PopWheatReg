"""Native selected-phenotype all-pair workflow; legacy modes are untouched."""
import gzip
import hashlib
import importlib.util
import json
import os
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

ROUTES = ('global', 'pw', 'sg1', 'sg2', 'sg3', 'sg5')
DIRECTORIES = ('f01_input', 'f02_raw', 'f03_merged', 'f04_manha', 'f05_logs')


def log(message):
    print(time.strftime('%Y-%m-%d %H:%M:%S'), message, flush=True)


def identity(path):
    path = Path(path).resolve()
    st = path.stat()
    return dict(path=str(path), size=st.st_size, mtime_ns=st.st_mtime_ns)


def requests(args):
    acr = bool(args.acr_id or args.acr_list)
    gene = bool(args.gene_id or args.gene_list)
    if acr == gene:
        raise ValueError('Specify ACR selectors OR gene selectors: --acr-id/--acr-list or --gene-id/--gene-list.')
    kind = 'caqtl' if acr else 'eqtl'
    if args.phenotype_kind and args.phenotype_kind != kind:
        raise ValueError('--phenotype_kind disagrees with the ID selectors.')
    values = list(args.acr_id if acr else args.gene_id)
    path = args.acr_list if acr else args.gene_list
    if path:
        with open(path) as handle:
            values += [line.strip() for line in handle if line.strip() and not line.lstrip().startswith('#')]
    if not values or any(not v or re.search(r'\s', v) for v in values):
        raise ValueError('Selectors must contain nonempty IDs, one ID per list line (no header).')
    return kind, list(dict.fromkeys(values))


def select_phenotypes(path, wanted, kind, core):
    """Scan once, retain only requested rows; audit duplicate/ambiguous matches."""
    short = {s for s in wanted if kind == 'caqtl' and re.fullmatch(r'acr\d+', s)}
    exact = set(wanted) - short
    rows, loc, matches = {}, {}, {s: [] for s in wanted}

    def retain(name, values, position=None):
        tokens = ([name] if name in exact else []) + ([name.split(':', 1)[0]] if name.split(':', 1)[0] in short else [])
        if not tokens:
            return
        if name in rows:
            raise ValueError('Duplicate requested phenotype row: '+name)
        try:
            array = np.asarray(values, dtype=float)
        except (ValueError, TypeError) as exc:
            raise ValueError('Selected phenotype contains nonnumeric/missing/nonfinite values: '+name) from exc
        if len(array) != len(samples) or not np.isfinite(array).all():
            raise ValueError('Selected phenotype has wrong width or nonfinite values: '+name)
        rows[name] = array
        loc[name] = position
        for token in tokens:
            matches[token].append(name)

    path = str(path)
    if path.endswith('.parquet'):
        pf = pq.ParquetFile(path)
        metadata = json.loads((pf.schema_arrow.metadata or {}).get(b'pandas', b'{}'))
        index_cols = metadata.get('index_columns', [])
        if len(index_cols) != 1 or not isinstance(index_cols[0], str):
            raise ValueError('Phenotype Parquet must have one stored string index (phenotype ID).')
        idcol = index_cols[0]
        samples = [c for c in pf.schema_arrow.names if c != idcol]
        core._check_unique(samples, 'phenotype sample IDs')
        for batch in pf.iter_batches(batch_size=2048):
            frame = batch.to_pandas()
            # iloc preserves duplicate-index auditing without .loc expansion.
            for i, name in enumerate(frame.index.astype(str)):
                if name in exact or name.split(':', 1)[0] in short:
                    retain(name, frame.iloc[i].to_numpy())
    else:
        bed = path.lower().endswith(('.bed', '.bed.gz'))
        offset = 4 if bed else 1
        opener = gzip.open if path.endswith('.gz') else open
        with opener(path, 'rt') as handle:
            fields = handle.readline().rstrip('\r\n').split('\t')
            samples = fields[offset:]
            core._check_unique(samples, 'phenotype sample IDs')
            if not samples:
                raise ValueError('No phenotype sample columns.')
            for number, line in enumerate(handle, 2):
                head = line.rstrip('\r\n').split('\t', offset)
                if len(head) != offset+1:
                    raise ValueError(f'Malformed phenotype row {number}.')
                name = head[3] if bed else head[0]
                if name not in exact and name.split(':', 1)[0] not in short:
                    continue
                position = (head[0], head[1], head[2]) if bed else None
                retain(name, head[offset].split('\t'), position)
    bad = {s: v for s, v in matches.items() if len(v) != 1}
    if bad:
        raise ValueError('Missing/ambiguous phenotype IDs: '+str(bad))
    order = list(dict.fromkeys(matches[s][0] for s in wanted))
    pheno = pd.DataFrame([rows[s] for s in order], index=order, columns=samples)
    pheno.index.name = 'phenotype_id'
    manifest = []
    for i, name in enumerate(order, 1):
        position = loc[name]
        if position is None and kind == 'caqtl':
            m = re.fullmatch(r'[^:]+:([^:]+)_(\d+)_(\d+)', name)
            if m:
                position = m.groups()
        if position is not None:
            ch, start, end = position
            if not str(start).isdigit() or not str(end).isdigit() or not (0 <= int(start) < int(end)):
                raise ValueError('Invalid BED0 interval for '+name)
        else:
            ch, start, end = '', '', ''
        stem = re.sub(r'[^A-Za-z0-9_.-]', '_', name.split(':', 1)[0])
        manifest.append(dict(order=i, phenotype_id=name, chromosome=ch, start=start, end=end,
                             manha_file=f'{stem}_result.txt'))
    core._check_unique([r['manha_file'] for r in manifest], 'Manhattan filenames')
    return pheno, pd.DataFrame(manifest)


def global_table(core, genotype, phenotype, covariates, maf_threshold):
    """Same installed TensorQTL kernels/formulas; include estimable r=0/P=1.

    Legacy global_pc excludes exact zero correlations through its t-threshold
    shortcut and computes 0/0 SE there. Only this new all-output path uses the
    continuous SE limit at r=0; existing modes/functions are not changed.
    """
    torch = core.torch
    device = core._get_device()
    g = torch.tensor(genotype.to_numpy(dtype=np.float32), device=device)
    y = torch.tensor(phenotype.to_numpy(dtype=np.float32), device=device)
    core.impute_mean(g)
    g, ids, af = core.filter_maf(g, genotype.index.to_numpy(), maf_threshold)
    columns = ['variant_id','phenotype_id','pval','beta','beta_se','af']
    if not len(ids):
        return pd.DataFrame(columns=columns)
    residualizer = None
    if covariates is not None:
        residualizer = core.Residualizer(torch.tensor(covariates.values, dtype=torch.float32, device=device))
    r, gv, pv = core.calculate_corr(g, y, residualizer=residualizer, return_var=True)
    valid = (torch.isfinite(r) & (r.abs() <= 1) & (gv[:, None] > 0) & (pv[None, :] > 0)
             & torch.isfinite(gv[:, None]) & torch.isfinite(pv[None, :]))
    ii, jj = valid.nonzero(as_tuple=True)
    if not len(ii):
        return pd.DataFrame(columns=columns)
    rr = r[ii, jj].to(torch.float64)
    dof = phenotype.shape[1]-2-(0 if covariates is None else covariates.shape[1])
    t = rr*torch.sqrt(dof/(1-rr.square()))
    ratio = torch.sqrt(pv[jj]/gv[ii])
    beta = rr*ratio
    se = beta/t
    zero = rr == 0
    se[zero] = torch.sqrt(pv[jj[zero]]/(dof*gv[ii[zero]])).to(se.dtype)
    return pd.DataFrame(dict(variant_id=np.asarray(ids)[ii.cpu().numpy()],
        phenotype_id=phenotype.index.to_numpy()[jj.cpu().numpy()],
        pval=core.get_t_pval(t.cpu().numpy(), dof, log=False),
        beta=beta.to(torch.float32).cpu().numpy(), beta_se=se.to(torch.float32).cpu().numpy(),
        af=af[ii].cpu().numpy()))


class Quiet:
    def write(self, *args):
        pass


def draw(root, regions):
    """Use only exported four-column data; never open original inputs."""
    complete = root/'f05_logs/03_complete.json'
    if not complete.is_file() or not json.loads(complete.read_text()).get('successful'):
        raise ValueError('--plot-only requires a completed selected run (f05_logs/03_complete.json).')
    parsed = []
    for text in regions:
        m = re.fullmatch(r'([^:]+):(\d+)-(\d+)', text)
        if not m or not (1 <= int(m[2]) < int(m[3])):
            raise ValueError('Invalid --plot-region; use CHR:START-END, 1-based inclusive bp.')
        parsed.append((m[1], int(m[2]), int(m[3])))
    os.environ.setdefault('MPLCONFIGDIR', '/tmp/tensorqtl_selected_mpl')
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'pdf.fonttype':42, 'font.size':10, 'axes.spines.top':False, 'axes.spines.right':False})
    manifest = pd.read_csv(root/'f01_input/02_matched_phenotypes.tsv', sep='\t')
    configuration = json.loads((root/'f01_input/04_run.json').read_text())
    chromosomes = pd.read_csv(root/'f01_input/05_chromosome_lengths.tsv', sep='\t', dtype={'CHR':str})
    offsets, cursor = {}, 0
    for r in chromosomes.itertuples():
        offsets[r.CHR] = cursor
        cursor += r.length
    ticks = [(offsets[r.CHR]+r.length/2)/1e6 for r in chromosomes.itertuples()]
    n = 6
    while (root/f'f{n:02d}_figures').exists():
        n += 1
    destination = root/f'f{n:02d}_figures'
    destination.mkdir()
    number = 1
    records = []
    for row in manifest.itertuples():
        points = pd.read_csv(root/'f04_manha'/row.manha_file, sep='\t', dtype={'CHR':str})
        for region in [None]+parsed:
            local = points if region is None else points[(points.CHR == region[0]) & points.BP.between(region[1],region[2])]
            fig, ax = plt.subplots(figsize=(14,4))
            yy = -np.log10(np.maximum(local.P.to_numpy(),1e-300))
            if region is None:
                for i, ch in enumerate(chromosomes.CHR):
                    use = local.CHR == ch
                    ax.scatter((offsets[ch]+local.loc[use,'BP'])/1e6, yy[use], s=4, linewidths=0,
                               color=['#436d9f','#aaaaaa'][i%2], rasterized=True)
                ax.set_xticks(ticks)
                ax.set_xticklabels(chromosomes.CHR, fontsize=8)
                ax.set_xlim(0, cursor/1e6)
                ax.set_xlabel('Source chromosome')
                suffix = 'genome'
            else:
                ax.scatter(local.BP/1e6, yy, s=8, linewidths=0, color='#436d9f', rasterized=True)
                ax.set_xlim(region[1]/1e6,region[2]/1e6)
                ax.set_xlabel(region[0]+' position (Mb)')
                suffix = f'{region[0]}_{region[1]}_{region[2]}'
            if not len(local):
                ax.text(.5,.5,'No estimable SNPs',transform=ax.transAxes,ha='center')
            ax.axhline(5, color='#bb5555', ls='--', lw=.8)
            ax.set_ylim(0, max(6, float(np.max(yy))*1.05 if len(yy) else 6))
            ax.set_ylabel('-log10(P)')
            ax.set_title(f"{row.phenotype_id} | {configuration['method']} | raw associations",fontsize=11)
            fig.tight_layout()
            for extension in ('pdf','png'):
                name = f'{number:02d}_{Path(row.manha_file).stem}_{suffix}.{extension}'
                fig.savefig(destination/name,dpi=160)
                records.append(dict(file=name,phenotype_id=row.phenotype_id,points=len(local),region=suffix))
                number += 1
            plt.close(fig)
            log(f'Plot {row.phenotype_id} {suffix}: {len(local):,} SNPs')
    (destination/f'{number:02d}_plot_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
    log('Plots complete: '+str(destination))


def run(args, core):
    started = time.monotonic()
    root = Path(args.out).resolve()
    if args.plot_only:
        draw(root, args.plot_region)
        return
    kind, wanted = requests(args)
    if args.method not in ROUTES:
        raise ValueError('selected requires --method global|pw|sg1|sg2|sg3|sg5.')
    if args.pval_threshold != 1:
        raise ValueError('selected exports all associations: explicitly set --pval_threshold 1 (P<1e-5 is exported separately).')
    if args.resume or args.omnibus_p_threshold is not None or args.pmin_threshold is not None:
        raise ValueError('selected does not use --resume or PW/Pmin output gates; use a new batch directory.')
    if args.method != 'global' and args.covariates:
        raise ValueError('PW/single-sg use the formal within-group no-PC model; omit --covariates.')
    if args.method != 'global' and not args.subgroups:
        raise ValueError('PW/single-sg require --subgroups.')
    if root.suffix == '.parquet':
        raise ValueError('selected --out is a batch DIRECTORY, not a .parquet file.')
    if root.exists() and (not root.is_dir() or any(root.iterdir())):
        raise ValueError('selected output directory must be new or empty; existing results are preserved.')
    if not root.parent.is_dir():
        raise ValueError('Output parent directory does not exist: '+str(root.parent))
    core.validate_parameters(1,args.maf_threshold,args.batch_size,args.pheno_chunk)
    log('Select requested phenotypes (one scan; no INT transformation)')
    phenotype, manifest = select_phenotypes(args.phenotype, wanted, kind, core)
    original_samples = phenotype.columns.tolist()
    core._check_unique(phenotype.index,'selected phenotype IDs')
    cov = core.load_covariates(args.covariates) if args.covariates else None
    subgroup = core.load_subgroups(args.subgroups) if args.subgroups else None
    genotype_ids = core.plink_sample_ids(args.genotype)
    if genotype_ids is not None:
        core._check_unique(genotype_ids,'genotype sample IDs')
        common = phenotype.columns.intersection(genotype_ids,sort=False)
        phenotype = phenotype.loc[:,common]
    if not phenotype.shape[1]:
        raise ValueError('No common genotype/phenotype samples.')
    if args.method != 'global':
        counts = core.validate_subgroups(subgroup,phenotype.columns)
        if args.method == 'pw' and set(counts) != {'sg1','sg2','sg3','sg5'}:
            raise ValueError('PW requires all four analysis groups sg1,sg2,sg3,sg5; got '+str(counts))
        if args.method != 'pw':
            if args.method not in counts:
                raise ValueError('Requested subgroup absent: '+args.method)
            phenotype = phenotype.loc[:,subgroup.reindex(phenotype.columns).eq(args.method)]
        subgroup = subgroup.reindex(phenotype.columns)
    if cov is not None:
        core._check_unique(cov.index,'covariate sample IDs')
        missing = phenotype.columns.difference(cov.index)
        if len(missing):
            raise ValueError('Analysis samples missing covariates: '+core._preview(missing))
        cov = cov.loc[phenotype.columns]
        design = np.column_stack([np.ones(len(cov)),cov.values])
        if not np.isfinite(design).all() or np.linalg.matrix_rank(design) != design.shape[1]:
            raise ValueError('Covariates are nonfinite or rank deficient with the intercept.')
    dof = phenotype.shape[1]-2-(0 if cov is None else cov.shape[1])
    if dof <= 0:
        raise ValueError('Not enough samples for the selected model.')
    if args.genotype.endswith('.gz'):
        core.check_raw_sample_header(args.genotype,4 if args.genotype.endswith('.bed.gz') else 1,'genotype')
    log(f'Load genotype ONCE: {phenotype.shape[0]} phenotypes; {phenotype.shape[1]} samples; route={args.method}')
    genotype, variant_pos = core.load_genotypes(args.genotype,select_samples=phenotype.columns.tolist())
    core._check_unique(genotype.index,'genotype variant IDs')
    common = phenotype.columns.intersection(genotype.columns,sort=False)
    phenotype = phenotype.loc[:,common]
    if cov is not None:
        cov = cov.loc[common]
    if subgroup is not None:
        subgroup = subgroup.loc[common]
    core.validate_matrices(genotype,phenotype,cov)
    core.validate_alignment(genotype,phenotype,cov)
    if phenotype.shape[1]-2-(0 if cov is None else cov.shape[1]) <= 0:
        raise ValueError('Not enough common samples for the selected model.')
    if cov is not None:
        design = np.column_stack([np.ones(len(cov)),cov.values])
        if np.linalg.matrix_rank(design) != design.shape[1]:
            raise ValueError('Covariates are rank deficient after selecting common samples.')
    if args.method != 'global':
        core.validate_subgroups(subgroup,phenotype.columns)
    # Existing coordinate/allele reader is shared with cis_all; no position
    # filtering or cis definition is applied in this all-pair workflow.
    spec = importlib.util.spec_from_file_location('selected_metadata',Path(__file__).with_name('01_cis_all_mode_run.py'))
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    variants = helper.variant_metadata(args,genotype,variant_pos,core)
    root.mkdir(exist_ok=True)
    for name in DIRECTORIES:
        (root/name).mkdir()
    progress_path = Path(args.logfile).resolve() if getattr(args,'logfile',None) else root/'f05_logs/04_progress.log'
    if not progress_path.parent.is_dir():
        raise ValueError('Log parent directory does not exist: '+str(progress_path.parent))
    if progress_path.exists():
        raise ValueError('Log already exists; existing logs are preserved: '+str(progress_path))
    progress_path.write_text(f'{time.strftime("%Y-%m-%d %H:%M:%S")} route={args.method}; '
                             f'phenotypes={len(phenotype)}; samples={phenotype.shape[1]}; variants={len(genotype)}\n')
    (root/'f01_input/01_requested_ids.txt').write_text('\n'.join(wanted)+'\n')
    manifest.to_csv(root/'f01_input/02_matched_phenotypes.tsv',sep='\t',index=False)
    phenotype.to_csv(root/'f01_input/03_selected_INT.tsv',sep='\t',float_format='%.17g')
    if cov is not None:
        cov.T.to_csv(root/'f01_input/06_covariates.tsv',sep='\t',float_format='%.17g')
    samples = pd.DataFrame(dict(sample_id=original_samples))
    samples['included'] = samples.sample_id.isin(phenotype.columns).astype(int)
    if args.subgroups:
        full_groups = core.load_subgroups(args.subgroups)
        samples['subgroup'] = samples.sample_id.map(full_groups)
    samples.to_csv(root/'f01_input/07_samples.tsv',sep='\t',index=False)
    paths = [args.phenotype]+[p for p in (args.covariates,args.subgroups,args.acr_list,args.gene_list,args.variant_metadata) if p]
    paths += [args.genotype+'.'+s for s in ('pgen','pvar','psam') if Path(args.genotype+'.'+s).is_file()]
    if not any(p.startswith(args.genotype+'.') for p in paths):
        paths += [args.genotype+'.'+s for s in ('bed','bim','fam') if Path(args.genotype+'.'+s).is_file()]
    if Path(args.genotype).is_file():
        paths.append(args.genotype)
    configuration = dict(mode='selected',method=args.method,kind=kind,maf_threshold=args.maf_threshold,
        pval_threshold=1,batch_size=args.batch_size,pheno_chunk=args.pheno_chunk,
        n_samples=phenotype.shape[1],n_covariates=0 if cov is None else cov.shape[1],
        samples=phenotype.columns.tolist(),input_sources=[identity(p) for p in paths],
        code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in
                     (Path(core.__file__),Path(__file__),Path(__file__).with_name('01_cis_all_mode_run.py'))},
        genotype_variants=len(genotype),INT_recomputed=False,progress_log=str(progress_path),
        subgroup_counts={} if args.method=='global' else core.validate_subgroups(subgroup,phenotype.columns))
    (root/'f01_input/04_run.json').write_text(json.dumps(configuration,indent=2)+'\n')
    lengths = variants.groupby('chrom').pos.max().to_dict()
    chromosomes = sorted(lengths,key=lambda s:[int(v) if v.isdigit() else v for v in re.split(r'(\d+)',s)])
    pd.DataFrame([dict(CHR=c,length=lengths[c]) for c in chromosomes]).to_csv(root/'f01_input/05_chromosome_lengths.tsv',sep='\t',index=False)
    summary = {s:dict(phenotype_id=s,n_pairs=0,n_significant=0,min_p=None,lead_snp=None) for s in phenotype.index}
    manha = dict(zip(manifest.phenotype_id,manifest.manha_file))
    for name in manha.values():
        (root/'f04_manha'/name).write_text('SNP\tCHR\tBP\tP\n')
    raw_path = root/'f02_raw/01_associations.parquet'
    text_path = root/'f03_merged/01_all_associations.txt'
    significant_path = root/'f03_merged/02_p_lt_1e5.txt'
    writer, text_header, raw_rows = None, True, 0
    # Fixed Arrow types make empty and nonempty batches schema-compatible.
    if args.method == 'global':
        schema = pa.schema([('variant_id',pa.string()),('phenotype_id',pa.string()),('pval',pa.float64()),
                            ('beta',pa.float32()),('beta_se',pa.float32()),('af',pa.float32())])
    else:
        schema = core._subgroup_parquet_schema(sorted(configuration['subgroup_counts']))
    writer = pq.ParquetWriter(raw_path,schema,compression='zstd')
    try:
        for begin in range(0,len(genotype),args.batch_size):
            g = genotype.iloc[begin:begin+args.batch_size]
            for a in range(0,len(phenotype),args.pheno_chunk or len(phenotype)):
                y = phenotype.iloc[a:a+(args.pheno_chunk or len(phenotype))]
                if args.method == 'global':
                    table = global_table(core,g,y,cov,args.maf_threshold)
                else:
                    table = core.run_subgroup_omnibus(g,y,subgroup,pval_threshold=1,maf_threshold=args.maf_threshold,
                        batch_size=len(g),pheno_chunk=len(y),prefilter=False,logger=Quiet(),verbose=False)
                writer.write_table(pa.Table.from_pandas(table,schema=schema,preserve_index=False))
                raw_rows += len(table)
                pcol = 'pval' if args.method=='global' else 'P_W' if args.method=='pw' else args.method+'_p'
                organized = table.copy()
                organized['pval'] = organized[pcol]
                if args.method.startswith('sg'):
                    organized['beta'] = organized[args.method+'_beta']
                    organized['beta_se'] = organized[args.method+'_se']
                organized['chr'] = organized.variant_id.map(variants.chrom)
                organized['pos'] = organized.variant_id.map(variants.pos)
                organized['method'] = args.method
                if organized[['chr','pos']].isna().any().any():
                    raise ValueError('Association lacks SNP coordinates.')
                if not organized.pval.between(0,1).all():
                    raise ValueError('Nonfinite or invalid P value in an estimable association.')
                organized.to_csv(text_path,sep='\t',index=False,mode='a',header=text_header,na_rep='NA',float_format='%.17g')
                organized[organized.pval<1e-5].to_csv(significant_path,sep='\t',index=False,mode='a',header=text_header,na_rep='NA',float_format='%.17g')
                text_header = False
                for name, part in organized.groupby('phenotype_id',sort=False):
                    s = summary[name]
                    s['n_pairs'] += len(part)
                    s['n_significant'] += int((part.pval<1e-5).sum())
                    if len(part):
                        lead = part.sort_values(['pval','chr','pos','variant_id'],kind='stable').iloc[0]
                        if s['min_p'] is None or (float(lead.pval),str(lead.variant_id)) < (s['min_p'],s['lead_snp']):
                            s['min_p'],s['lead_snp'] = float(lead.pval),str(lead.variant_id)
                    part[['variant_id','chr','pos','pval']].to_csv(root/'f04_manha'/manha[name],sep='\t',index=False,
                        header=False,mode='a',float_format='%.17g')
            message = f'Association progress: SNPs {min(begin+args.batch_size,len(genotype)):,}/{len(genotype):,}; pairs={raw_rows:,}; elapsed={time.monotonic()-started:.1f}s'
            log(message)
            with progress_path.open('a') as handle:
                handle.write(time.strftime('%Y-%m-%d %H:%M:%S')+' '+message+'\n')
    finally:
        writer.close()
    rows = []
    for s in summary.values():
        s['status'] = 'OK' if s['n_pairs'] else 'no_estimable_pairs'
        s['n_input_variants'] = len(genotype)
        s['n_not_estimable_or_MAF_filtered'] = len(genotype)-s['n_pairs']
        rows.append(s)
    pd.DataFrame(rows).to_csv(root/'f05_logs/01_summary.tsv',sep='\t',index=False,na_rep='NA',float_format='%.17g')
    assert sum(s['n_pairs'] for s in rows) == pq.ParquetFile(raw_path).metadata.num_rows == raw_rows
    (root/'f05_logs/02_run_notes.txt').write_text(
        'No re-INT, peak filtering, region QC, Intro-QC, rescue or union.\n'
        'P<1e-5 counts refer to raw associations, not final filtered QTL.\n'
        'Global uses installed TensorQTL kernels; exact r=0 emits P=1 with finite continuous-limit SE.\n'
        'Single-sg reuses the subgroup exact path with one group; raw preserves its native columns.\n'
        'Chromosome lengths in plots are maximum loaded SNP positions, fixed across phenotypes.\n')
    (root/'f05_logs/03_complete.json').write_text(json.dumps(dict(successful=True,phenotypes=len(phenotype),
        pairs=raw_rows,method=args.method,n_samples=phenotype.shape[1],elapsed_seconds=time.monotonic()-started),indent=2)+'\n')
    with progress_path.open('a') as handle:
        handle.write('SELECTED COMPLETE\n')
    log('SELECTED COMPLETE: '+str(root))
    if args.plot:
        draw(root,args.plot_region)
