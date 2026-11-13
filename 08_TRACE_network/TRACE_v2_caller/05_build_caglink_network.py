#!/usr/bin/env python3
"""Fresh C-agLINK V2.2 build from standard edges; no old network/cache inputs.

SNP batches are disjoint. The functional dependency (SNP,TF)->tf_snp_meta
allows wide Local metadata to be reattached only when publishing evidence.
No evidence paths are removed by this normalization. Global QC operates on
weighted narrow records, retaining exact unique-pair set calculations.
"""
import argparse
import hashlib
import importlib.util
import json
import os
import resource
import threading
import time
from datetime import datetime
from pathlib import Path

import duckdb

HERE=Path(__file__).parent
spec=importlib.util.spec_from_file_location('trace_network_reference',HERE/'04_network_rules.py')
N=importlib.util.module_from_spec(spec)
spec.loader.exec_module(N)
VERSION="C-agLINK-v2.2"
NETWORK_ID_FORMAT='TFN%09d'
FLAGS=['pass_coordinates','pass_no_self','pass_tf_gene','pass_tf_acr','pass_snp_gene']
INPUTS=['eqtl','caqtl','aglink','coacr','binding','tf_snp','gene_bed']


def log(message):
    print(f'[Caglink_network {time.strftime("%H:%M:%S")}] {message}',flush=True)


def identity(path):
    p=Path(path).resolve(); st=p.stat()
    return dict(path=str(p),bytes=st.st_size,mtime_ns=st.st_mtime_ns)


def saved(con,key):
    row=con.execute('SELECT value FROM progress WHERE key=?',[key]).fetchone()
    return json.loads(row[0]) if row else None


def save(con,key,value):
    con.execute('INSERT OR REPLACE INTO progress VALUES (?,?)',[key,json.dumps(value)])


def stage(con,key,callback):
    previous=saved(con,key)
    if previous is not None:
        for info in previous.get('files',[]):
            if identity(info['path'])!=info:
                raise ValueError('Completed output changed: '+info['path'])
        log(f'Resume: {key} completed; skip')
        return previous
    log('Start: '+key)
    started=time.monotonic(); stopped=threading.Event()
    def beat():
        while not stopped.wait(30):
            peak=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2
            log(f'{key}: elapsed={time.monotonic()-started:.0f}s; peak_RSS={peak:.1f} GiB')
    worker=threading.Thread(target=beat,daemon=True); worker.start()
    try:
        con.execute('BEGIN')
        try:
            value=callback() or {}
            save(con,key,value)
            con.execute('COMMIT')
        except BaseException:
            con.execute('ROLLBACK')
            raise
        con.execute('CHECKPOINT')
    finally:
        stopped.set(); worker.join()
    log(f'Completed: {key}; elapsed={time.monotonic()-started:.1f}s')
    return value


def sql_stage(con,key,sql,validate=None):
    def build():
        for query in sql: con.execute(query)
        if validate: validate()
    return stage(con,key,build)



def binding_route(con,binding='raw_binding',dest='route_binding'):
    # Gate ONLY Route3. The full tf_snp remains available to Route1 and Route2.
    con.execute(f"""CREATE VIEW {dest} AS
        SELECT t.src_id snp,b.src_id tf,h.gene,t.tf_snp_type,'binding' route,
          h.acr1,NULL::VARCHAR acr2,'none' coacr_label,'aglink' link_label,
          concat_ws(';',replace(replace(b.meta,'PCC=','tf_acr_pcc='),'effect=','legacy_PCC_sign_label='),
                    'PCC_interpretation=correlation_only',concat('aglink_',h.gene_meta)) support,
          t.meta tf_snp_meta
        FROM {binding} b JOIN tf_snp t ON trim(b.src_id)=t.dst_id
        JOIN route3_gate g ON t.src_id=g.snp AND t.dst_id=g.tf AND g.keep_C_agLINK
        JOIN acr_gene_direct h ON b.dst_id=h.acr1 WHERE h.link_type='agLINK'""")


N.binding_route=binding_route


def prepare_route3_gate(con,args):
    sql_stage(con,'gate_locations',[r"""CREATE TABLE gate_locations AS WITH located AS (
      SELECT t.src_id snp,t.dst_id tf,t.tf_snp_type sources,
        contains(t.tf_snp_type,'local_unit') has_local_unit,
        contains(t.tf_snp_type,'tf_prox_snpeff') has_snpeff_prox,
        g.chr tf_chr,g.start0,g.finish,g.tss,
        nullif(regexp_extract(t.src_id,'^([^:_\-]+)[:_\-]([0-9]+)$',1),'') snp_chr,
        try_cast(regexp_extract(t.src_id,'^([^:_\-]+)[:_\-]([0-9]+)$',2) AS BIGINT) snp_pos1
      FROM tf_snp_all t LEFT JOIN gene_loc g ON t.dst_id=g.gene_id)
      SELECT *,CASE WHEN snp_chr=tf_chr THEN
          greatest(start0+1-snp_pos1,snp_pos1-finish,0) END snp_tf_distance_bp,
        coalesce(snp_chr=tf_chr AND snp_pos1 IS NOT NULL AND
          greatest(start0+1-snp_pos1,snp_pos1-finish,0)<=5000,false) within_5kb
      FROM located"""])
    sql_stage(con,'gate_cis_eqtl',["""CREATE TABLE gate_eq AS
      SELECT DISTINCT p.snp,p.tf FROM gate_locations p JOIN raw_eqtl_all e
        ON p.snp=trim(e.src_id) AND p.tf=trim(e.dst_id)
      WHERE lower(trim(e.edge_type))='cis-eqtl'"""])
    # Actual direct agLINK row required. Keep historical Local coordinates:
    # SNP 1-based; floor(ACR BED midpoint); '+' TSS=start0, '-' TSS=end1.
    sql_stage(con,'gate_cis_caqtl',[r"""CREATE TABLE gate_ca_raw AS WITH joined AS (
      SELECT DISTINCT p.snp,p.tf,c.dst_id acr,p.snp_chr,p.snp_pos1,p.tf_chr,p.tss,
        nullif(regexp_extract(split_part(c.dst_id,':',2),'^([^_]+)_([0-9]+)_([0-9]+)$',1),'') acr_chr,
        try_cast(regexp_extract(split_part(c.dst_id,':',2),'^([^_]+)_([0-9]+)_([0-9]+)$',2) AS BIGINT) acr_start0,
        try_cast(regexp_extract(split_part(c.dst_id,':',2),'^([^_]+)_([0-9]+)_([0-9]+)$',3) AS BIGINT) acr_end1
      FROM raw_caqtl_all c JOIN gate_locations p ON trim(c.src_id)=p.snp
      JOIN raw_aglink a ON c.dst_id=a.src_id AND p.tf=trim(a.dst_id)
      WHERE lower(trim(c.edge_type))='cis-caqtl' AND a.edge_type='agLINK'),
      centered AS (SELECT *,cast(floor((acr_start0+acr_end1)/2.0) AS BIGINT) acr_midpoint FROM joined)
      SELECT *,greatest(snp_pos1,acr_midpoint,tss)-least(snp_pos1,acr_midpoint,tss) local_span_bp
      FROM centered""",
      f"""CREATE TABLE gate_ca AS SELECT snp,tf,count(DISTINCT acr) n_aglink_cis_acrs
      FROM gate_ca_raw WHERE snp_chr=tf_chr AND snp_chr=acr_chr
        AND snp_pos1 IS NOT NULL AND acr_midpoint IS NOT NULL AND tss IS NOT NULL
        AND local_span_bp<={args.local_span} GROUP BY snp,tf"""])
    def validate_gate():
        bad=con.execute("""SELECT count(*) FROM route3_gate
            WHERE has_local_unit AND NOT keep_C_agLINK""").fetchone()[0]
        duplicates=con.execute('SELECT count(*)-count(DISTINCT(snp,tf)) FROM route3_gate').fetchone()[0]
        if bad or duplicates:
            raise ValueError(f'Invalid C-agLINK gate: lost_LU={bad}, duplicates={duplicates}')
    sql_stage(con,'route3_gate',["""CREATE TABLE route3_gate AS WITH annotated AS (
      SELECT p.*,coalesce(e.snp IS NOT NULL,false) has_cis_eqtl,
        coalesce(c.n_aglink_cis_acrs,0) n_aglink_cis_acrs
      FROM gate_locations p LEFT JOIN gate_eq e USING(snp,tf) LEFT JOIN gate_ca c USING(snp,tf))
      SELECT *,coalesce(has_local_unit OR (has_snpeff_prox AND within_5kb AND
        (has_cis_eqtl OR n_aglink_cis_acrs>0)),false) keep_C_agLINK FROM annotated"""],validate_gate)
    def export_gate():
        folder=Path(args.out_prefix).parent
        files=[copy_atomic(con,'SELECT * FROM route3_gate ORDER BY snp,tf',
            folder/'01_Caglink_SNP_TF_entrances.tsv',"HEADER,DELIMITER '\t',NULL 'NA'")]
        qc="""SELECT 'all_input_pairs' item,count(*) AS "value" FROM route3_gate
          UNION ALL SELECT 'C_agLINK_pairs',count(*) FROM route3_gate WHERE keep_C_agLINK
          UNION ALL SELECT 'C_agLINK_snps',count(DISTINCT snp) FROM route3_gate WHERE keep_C_agLINK
          UNION ALL SELECT 'C_agLINK_tfs',count(DISTINCT tf) FROM route3_gate WHERE keep_C_agLINK
          UNION ALL SELECT 'local_unit_kept',count(*) FROM route3_gate WHERE has_local_unit AND keep_C_agLINK
          UNION ALL SELECT 'additional_prox_eqtl',count(*) FROM route3_gate
            WHERE keep_C_agLINK AND NOT has_local_unit AND has_cis_eqtl
          UNION ALL SELECT 'additional_prox_caqtl_only',count(*) FROM route3_gate
            WHERE keep_C_agLINK AND NOT has_local_unit AND NOT has_cis_eqtl AND n_aglink_cis_acrs>0
          UNION ALL SELECT 'prox_outside_5kb',count(*) FROM route3_gate WHERE has_snpeff_prox AND NOT within_5kb"""
        files.append(copy_atomic(con,qc,folder/'02_Caglink_entrance_QC.tsv',"HEADER,DELIMITER '\t',NULL 'NA'"))
        return dict(files=files)
    stage(con,'gate_outputs',export_gate)


def prepare(con,args):
    def validate_tf():
        if con.execute('SELECT count(*)-count(DISTINCT(src_id,dst_id)) FROM tf_snp_all').fetchone()[0]:
            raise ValueError('Duplicate SNP-TF pairs; rebuild with 06')
        bad=con.execute("""SELECT * FROM tf_snp_all WHERE tf_snp_type IS NULL OR
            NOT regexp_full_match(tf_snp_type,'(local_unit|tf_prox_snpeff)(\\|(local_unit|tf_prox_snpeff))?') LIMIT 1""").fetchone()
        if bad: raise ValueError(f'Unexpected/obsolete SNP-TF source: {bad}')
    sql_stage(con,'tf_snp',[f"""CREATE TABLE tf_snp_all AS SELECT trim(src_id) src_id,
        trim(dst_id) dst_id,score tf_snp_type,meta FROM {N.read(args.tf_snp)}"""],validate_tf)
    sql_stage(con,'snps',[f"""CREATE TABLE snps AS SELECT src_id,
        cast(floor((row_number() OVER(ORDER BY src_id)-1)/{args.batch_snps}.0) AS BIGINT) batch
        FROM (SELECT DISTINCT src_id FROM tf_snp_all)"""])
    def validate_genes():
        if con.execute('SELECT count(*)-count(DISTINCT gene_id) FROM gene_loc').fetchone()[0]:
            raise ValueError('Duplicate gene IDs in gene BED')
    sql_stage(con,'genes',[f"""CREATE TABLE gene_loc AS SELECT column0 chr,
        cast(column1 AS BIGINT) start0,cast(column2 AS BIGINT) finish,column3 gene_id,
        (cast(column1 AS BIGINT)+cast(column2 AS BIGINT))/2 center,
        CASE WHEN column5='+' THEN cast(column1 AS BIGINT)
             ELSE cast(column2 AS BIGINT) END tss
        FROM {N.read(args.gene_bed,False)}"""],validate_genes)
    for name in ('eqtl','caqtl'):
        # Scan each QTL input once: retain cis for the new entrance gate.
        # Route1/2 queries still select only their original trans records.
        sql_stage(con,'input_'+name,[f"""CREATE TABLE raw_{name}_all AS SELECT q.*,s.batch
            FROM {N.read(getattr(args,name))} q JOIN snps s ON trim(q.src_id)=s.src_id
            WHERE q.edge_type LIKE 'trans%' OR lower(trim(q.edge_type)) IN ('cis-eqtl','cis-caqtl')
            ORDER BY s.batch"""])
    for name in ('aglink','coacr','binding'):
        def validate_binding():
            bad=con.execute("""SELECT src_id,dst_id,meta FROM raw_binding
                WHERE try_cast(regexp_extract(meta,'(^|;)PCC=([^;]+)',2) AS DOUBLE) IS NULL
                  OR NOT isfinite(try_cast(regexp_extract(meta,'(^|;)PCC=([^;]+)',2) AS DOUBLE))
                  OR abs(try_cast(regexp_extract(meta,'(^|;)PCC=([^;]+)',2) AS DOUBLE))<0.35 LIMIT 1""").fetchone()
            if bad: raise ValueError(f'Formal Binding must remain signed |PCC|>=0.35: {bad}')
        sql_stage(con,'input_'+name,[f'CREATE TABLE raw_{name} AS SELECT * FROM {N.read(getattr(args,name))}'],
                  validate_binding if name=='binding' else None)
    prepare_route3_gate(con,args)
    stage(con,'acr_hub',lambda: N.acr_hub(con))
    log('Disk inputs prepared: '+str({t:con.execute(f'SELECT count(*) FROM {t}').fetchone()[0]
        for t in ('tf_snp_all','raw_eqtl_all','raw_caqtl_all','acr_gene_hub')}))


def copy_atomic(con,query,path,options):
    path=Path(path); partial=Path(str(path)+'.partial')
    # Only this checkpoint's explicit output/part paths are used here; never
    # remove production inputs or unrelated files. The fresh prefix is reserved.
    if partial.exists(): partial.unlink()
    con.execute(f'COPY ({query}) TO {N.lit(partial)} ({options})')
    os.replace(partial,path)
    return identity(path)


def batch(con,args,batch_id,part_dir):
    path=part_dir/f'01_batch_{batch_id+1:06d}.parquet'
    def build():
        # The true meta remains in tf_snp_all, unique by (SNP,TF). NULL here
        # avoids duplicating long Local evidence in joins, DISTINCT and QC.
        con.execute(f"""CREATE OR REPLACE VIEW tf_snp AS
            SELECT t.src_id,t.dst_id,t.tf_snp_type,NULL::VARCHAR meta
            FROM tf_snp_all t JOIN snps s USING(src_id) WHERE s.batch={batch_id}""")
        for name in ('eqtl','caqtl'):
            con.execute(f'CREATE OR REPLACE VIEW raw_{name} AS SELECT * EXCLUDE(batch) '
                        f'FROM raw_{name}_all WHERE batch={batch_id}')
        for view in ('route_eqtl','route_caqtl','route_binding'):
            con.execute(f'DROP VIEW IF EXISTS {view}')
        log(f'batch_{batch_id}: build distinct Route evidence')
        N.routes(con)
        log(f'batch_{batch_id}: annotate final-filter flags')
        N.flags(con,'all_evidence','flagged',args.same_chr_trans_min_distance)
        log(f'batch_{batch_id}: aggregate network (distinct lists, no per-group ordered string sort)')
        N.save_network(con,args.out_prefix,write_outputs=False,low_memory=True)
        log(f'batch_{batch_id}: append compact network/QC; export evidence batch')
        con.execute('CREATE TABLE IF NOT EXISTS network_rows AS SELECT *,0::BIGINT batch '
                    'FROM network_main LIMIT 0')
        con.execute(f'INSERT INTO network_rows SELECT *,{batch_id}::BIGINT FROM network_main')
        keys='snp,tf,gene,route,'+','.join(FLAGS)
        con.execute(f'CREATE TABLE IF NOT EXISTS qc_rows AS SELECT {keys},count(*)::BIGINT n_records '
                    f'FROM flagged WHERE false GROUP BY {keys}')
        con.execute(f'INSERT INTO qc_rows SELECT {keys},count(*)::BIGINT FROM flagged GROUP BY {keys}')
        # This compact artifact retains all final evidence and distance fields;
        # the final publisher reattaches tf_snp_meta without changing the schema.
        info=copy_atomic(con,'SELECT * EXCLUDE(tf_snp_meta) FROM filtered_evidence',path,
                         "FORMAT PARQUET,COMPRESSION ZSTD,ROW_GROUP_SIZE 8192")
        metrics=dict(raw_evidence=con.execute('SELECT count(*) FROM all_evidence').fetchone()[0],
            retained_evidence=con.execute('SELECT count(*) FROM filtered_evidence').fetchone()[0],
            network_rows=con.execute('SELECT count(*) FROM network_main').fetchone()[0])
        for table in ('all_evidence','flagged','filtered_evidence','network_main'):
            con.execute(f'DROP TABLE {table}')
        return dict(files=[info],batch=batch_id,**metrics)
    return stage(con,f'batch_{batch_id}',build)


def finish_index(con):
    def build():
        # network_rows is already grouped within disjoint SNP batches; no second
        # global evidence GROUP BY or wide DISTINCT is needed.
        con.execute("""CREATE TABLE network_index AS
            SELECT row_number() OVER(ORDER BY snp,tf,gene) net_id,* EXCLUDE(net_id)
            FROM network_rows""")
    stage(con,'network_index',build)
    con.execute('CREATE OR REPLACE VIEW network_main AS SELECT * EXCLUDE(batch) FROM network_index')
    con.execute('CREATE OR REPLACE VIEW flagged AS SELECT * FROM qc_rows')


def publish(con,args,part_dir,batches):
    prefix=args.out_prefix
    def qc():
        path=Path(prefix+'.filter_qc.tsv'); partial=Path(str(path)+'.partial')
        if partial.exists(): partial.unlink()
        N.filter_qc(con,str(partial),weight_column='n_records')
        os.replace(partial,path)
        return dict(files=[identity(path)])
    stage(con,'filter_qc',qc)
    stage(con,'network_output',lambda: dict(files=[copy_atomic(con,
        "SELECT printf('TFN%09d',net_id) net_id,* EXCLUDE(net_id) FROM network_main ORDER BY network_main.net_id",prefix+'.network.tsv',
        "HEADER,DELIMITER '\t',NULL 'NA'")]))
    n=con.execute('SELECT count(*) FROM network_main').fetchone()[0]
    for number,start in enumerate(range(1,max(n,1)+1,args.evidence_part_net_ids),1):
        end=start+args.evidence_part_net_ids-1
        def evidence(number=number,start=start,end=end):
            selected=[x[0] for x in con.execute(
                f'SELECT DISTINCT batch FROM network_index WHERE net_id BETWEEN {start} AND {end}').fetchall()]
            # Empty networks still publish one correctly typed empty evidence part.
            paths=[str(part_dir/f'01_batch_{i+1:06d}.parquet') for i in (selected or batches[:1])]
            files='['+','.join(N.lit(p) for p in paths)+']'
            query=f"""SELECT printf('TFN%09d',n.net_id) net_id,e.snp,e.tf,e.gene,e.tf_snp_type,e.route,
                coalesce(e.acr1,'NA') acr1,coalesce(e.acr2,'NA') acr2,e.coacr_label,e.link_label,e.support,
                t.meta tf_snp_meta,e.qtl_methods,e.qtl_p_source,e.qtl_locus_id,e.qtl_p,e.qtl_beta,e.qtl_beta_se,
                e.tf_gene_center_distance_bp,e.tf_acr_center_distance_bp,e.snp_gene_interval_distance_bp
                FROM read_parquet({files}) e JOIN network_main n USING(snp,tf,gene)
                JOIN tf_snp_all t ON e.snp= t.src_id AND e.tf=t.dst_id
                WHERE n.net_id BETWEEN {start} AND {end}"""
            # e.snp is a raw QTL ID in R1/2, while the original join used trim.
            # Reattach using the identical trim rule, without altering output IDs.
            query=query.replace('e.snp= t.src_id','trim(e.snp)=t.src_id')
            info=copy_atomic(con,query,prefix+f'.evidence.part{number}.parquet',
                             'FORMAT PARQUET,COMPRESSION ZSTD,ROW_GROUP_SIZE 8192')
            return dict(files=[info])
        stage(con,f'evidence_output_{number}',evidence)
        log(f'Evidence part {number}/{max(1,(n+args.evidence_part_net_ids-1)//args.evidence_part_net_ids)} published')


def run(args):
    prefix=Path(args.out_prefix).resolve(); args.out_prefix=str(prefix)
    db=Path(args.checkpoint_db or str(prefix.parent/'f01_network_checkpoint/01_network.duckdb')).resolve()
    spill=Path(args.tmp_dir or str(db.parent/'f02_spill')).resolve()
    if min(args.threads,args.batch_snps,args.evidence_part_net_ids)<1 or min(args.same_chr_trans_min_distance,args.local_span)<0:
        raise ValueError('Invalid resource/batch/distance parameter')
    infos={key:identity(getattr(args,key)) for key in INPUTS}
    fingerprint=dict(version=VERSION,inputs=infos,out_prefix=str(prefix),
        distance=args.same_chr_trans_min_distance,local_span=args.local_span,batch_snps=args.batch_snps,
        evidence_part_net_ids=args.evidence_part_net_ids,
        network_id_format=NETWORK_ID_FORMAT,
        script_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                       for p in (Path(__file__),HERE/'04_network_rules.py')})
    if db.exists() and not args.resume:
        raise FileExistsError(f'Checkpoint exists; use --resume: {db}')
    for folder in (prefix.parent,db.parent,spill): folder.mkdir(parents=True,exist_ok=True)
    con=duckdb.connect(str(db))
    try:
        con.execute(f'SET memory_limit={N.lit(args.ram_limit)}')
        con.execute(f'SET threads={args.threads}')
        con.execute(f'SET temp_directory={N.lit(spill)}')
        con.execute('CREATE TABLE IF NOT EXISTS progress(key VARCHAR PRIMARY KEY,value VARCHAR)')
        previous=saved(con,'identity')
        if previous is not None and previous!=fingerprint:
            raise ValueError('Input/code/semantic parameters changed; use a new checkpoint and output prefix')
        if previous is None:
            existing=list(prefix.parent.glob(prefix.name+'.*'))
            if existing: raise FileExistsError(f'Refusing to overwrite existing network files: {existing[:3]}')
            save(con,'identity',fingerprint); con.execute('CHECKPOINT')
        # Validate committed artifacts before resuming, including already exported
        # Parquet parts. Manual edits are never silently overwritten.
        for key,value in con.execute('SELECT key,value FROM progress').fetchall():
            for info in json.loads(value).get('files',[]):
                if identity(info['path'])!=info:
                    raise ValueError('Completed output changed: '+info['path'])
        if saved(con,'manifest_output') is not None:
            log('All formal outputs completed; nothing to rerun')
            return
        prepare(con,args)
        batches=[r[0] for r in con.execute('SELECT DISTINCT batch FROM snps ORDER BY batch').fetchall()] or [0]
        part_dir=db.parent/'f01_evidence_batches'; part_dir.mkdir(exist_ok=True)
        total_rows=0; total_evidence=0
        for index,batch_id in enumerate(batches,1):
            done=batch(con,args,batch_id,part_dir)
            total_rows+=done['network_rows']; total_evidence+=done['retained_evidence']
            log(f'Batch {index}/{len(batches)} completed; network rows={total_rows:,}; final evidence={total_evidence:,}')
        finish_index(con)
        publish(con,args,part_dir,batches)
        result=N.counts(con,'flagged',N.PASS,weight_column='n_records')
        def manifest():
            cfg=vars(args).copy()
            # 08 reads this original schema and imports 07's reference functions.
            data=dict(created=datetime.now().isoformat(),parameters=cfg,
                script=str(Path(__file__).resolve()),
                script_sha256=fingerprint['script_sha256'][Path(__file__).name],
                execution_script=str(Path(__file__).resolve()),execution_fingerprint=fingerprint,
                network_id_format=NETWORK_ID_FORMAT,
                basis='fresh Route1/2 unchanged; Route3 C-agLINK entrance; downstream direct agLINK',
                route3_entrance='local_unit OR (allANN_prox AND body_distance<=5000 AND (same_SNP_cis_eQTL_TF OR same_SNP_cis_caQTL_ACR_direct_agLINK_TF))',
                local_span=args.local_span,
                binding_policy='abs(PCC)>=0.35; signed PCC retained; not activation/repression inference',
                snp_gene_distance='SNP 1-based position to inclusive gene bases [BED.start+1,BED.end]',
                center_distance='unchanged v4 (BED.start+BED.end)/2; strictly > threshold',inputs=infos,
                execution=dict(checkpoint_db=str(db),completed_batches=len(batches),final_counts=list(result),
                    peak_RSS_GiB=round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,3),
                    note='Disk-backed SNP batches; Local meta reattached only at final evidence export; memory budget is not hard RSS cap'))
            path=Path(str(prefix)+'.manifest.json'); partial=Path(str(path)+'.partial')
            partial.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n'); os.replace(partial,path)
            return dict(files=[identity(path)])
        stage(con,'manifest_output',manifest)
        log(f'DONE: evidence/SNP-TF-gene/TF-gene/TF = {result}')
    finally:
        con.close()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('eqtl','caqtl','aglink','coacr','binding','tf-snp','gene-bed','out-prefix'):
        p.add_argument('--'+key,required=True)
    p.add_argument('--same-chr-trans-min-distance',type=int,default=500000)
    p.add_argument('--local-span',type=int,default=500000)
    p.add_argument('--ram-limit',default='64GB')
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--tmp-dir')
    p.add_argument('--checkpoint-db')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--batch-snps',type=int,default=500)
    p.add_argument('--evidence-part-net-ids',type=int,default=200000)
    run(p.parse_args())


if __name__=='__main__':
    main()

