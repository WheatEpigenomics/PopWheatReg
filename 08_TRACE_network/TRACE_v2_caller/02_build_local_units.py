#!/usr/bin/env python3
"""Disk-backed legacy L1/L2 construction; no changes to biological rules.

One input scan per committed file; trans edges are never materialized. All
relationships for one SNP stay in the same batch, including direct-path tests.
The old pandas quicksort did not define tie ordering: equal-priority duplicate
evidence now uses input order deterministically. Unit IDs/order are regenerated.
"""
import argparse
import csv
import gzip
import hashlib
import json
import os
import resource
import threading
import time
from pathlib import Path

import duckdb

VERSION = 1
HEADER = ['unit_id', 'template', 'snp', 'acr1', 'acr2', 'gene', 'span_bp',
          'acr_gene_type', 'support']
EDGE_HEADER = ['src_type', 'src_id', 'dst_type', 'dst_id', 'edge_type', 'score', 'meta']


def log(message):
    print(f'[local_disk {time.strftime("%H:%M:%S")}] {message}', flush=True)


def lit(value):
    return "'" + str(value).replace("'", "''") + "'"


def identity(path):
    path = Path(path).resolve()
    st = path.stat()
    return dict(path=str(path), size=st.st_size, mtime_ns=st.st_mtime_ns)


def saved(con, key):
    row = con.execute('SELECT value FROM progress WHERE key=?', [key]).fetchone()
    return json.loads(row[0]) if row else None


def save(con, key, value):
    con.execute('INSERT OR REPLACE INTO progress VALUES (?, ?)', [key, json.dumps(value)])


def stage(con, key, sql):
    if saved(con, key) is not None:
        log(f'Resume: {key} completed; skip')
        return
    log(f'Start: {key}')
    started = time.monotonic()
    stop = threading.Event()

    def beat():
        while not stop.wait(30):
            rss = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024**2
            log(f'{key}: elapsed={time.monotonic()-started:.0f}s; peak_RSS={rss:.1f} GiB')

    worker = threading.Thread(target=beat, daemon=True)
    worker.start()
    try:
        con.execute('BEGIN')
        try:
            for query in sql:
                con.execute(query)
            save(con, key, True)
            con.execute('COMMIT')
        except BaseException:
            con.execute('ROLLBACK')
            raise
        con.execute('CHECKPOINT')
    finally:
        stop.set()
        worker.join()
    log(f'Completed: {key}; elapsed={time.monotonic()-started:.1f}s')


def read_edges(path):
    return (f"read_csv({lit(path)}, delim='\\t', header=true, all_varchar=true, "
            "auto_detect=true, nullstr='NA')")


def snp_chr(col):
    return f"nullif(regexp_extract({col}, '^([^:_\\-]+)[:_\\-]([0-9]+)$', 1), '')"


def snp_pos(col):
    return f"try_cast(regexp_extract({col}, '^([^:_\\-]+)[:_\\-]([0-9]+)$', 2) AS BIGINT)"


def acr_coord(col, part):
    # Legacy strips everything through the first colon, then uses floor midpoint.
    text = f"CASE WHEN contains({col}, ':') THEN substr({col}, strpos({col}, ':')+1) ELSE {col} END"
    pattern = '^([^_]+)_([0-9]+)_([0-9]+)$'
    if part == 'chr':
        return f"nullif(regexp_extract({text}, {lit(pattern)}, 1), '')"
    start = f"try_cast(regexp_extract({text}, {lit(pattern)}, 2) AS BIGINT)"
    end = f"try_cast(regexp_extract({text}, {lit(pattern)}, 3) AS BIGINT)"
    return f'cast(floor(({start}+{end})/2.0) AS BIGINT)'


def prepare(con, args):
    for i, path in enumerate(args.edges):
        opener = gzip.open if str(path).endswith('.gz') else open
        with opener(path, 'rt') as handle:
            names = next(csv.reader(handle, delimiter='\t'))
        if not set(EDGE_HEADER).issubset(names):
            raise ValueError(f'Not a seven-column standard edge table: {path}')
        # Single reader preserves deterministic source order. Filtering/projection
        # happens in DuckDB, not in an all-input pandas DataFrame.
        con.execute('SET threads=1')
        stage(con, f'input_{i+1}', [f"""
            CREATE TABLE input_{i} AS
            SELECT {i}::INTEGER file_id, row_number() OVER () ord,
                   src_id, dst_id, trim(edge_type) edge_type, coalesce(meta,'nan') meta
            FROM {read_edges(path)}
            WHERE lower(trim(edge_type)) IN ('cis-eqtl','cis-caqtl','co-acr')
               OR regexp_matches(lower(trim(edge_type)), 'aglink|distance')
        """])
        con.execute(f'SET threads={args.threads}')
        if args.report_edge_types:
            log(f'{path}: retained edge types '+str(con.execute(
                f'SELECT edge_type,count(*) FROM input_{i} GROUP BY edge_type ORDER BY edge_type').fetchall()))
    union = ' UNION ALL '.join(f'SELECT * FROM input_{i}' for i in range(len(args.edges)))
    con.execute(f'CREATE OR REPLACE VIEW raw AS {union}')
    bed = (f"read_csv({lit(args.gene_bed)}, delim='\\t', header=false, "
           "all_varchar=true, auto_detect=true)")
    stage(con, 'genes', [f"""CREATE TABLE genes AS SELECT column3 gene, column0 gene_chr,
        CASE WHEN column5='+' THEN cast(column1 AS BIGINT) ELSE cast(column2 AS BIGINT) END tss
        FROM {bed}"""])
    # These early conditions are necessary conditions of the final span check;
    # they remove no L1/L2 unit that the original script could retain.
    stage(con, 'eq', [f"""CREATE TABLE eq AS
        SELECT r.*,g.gene_chr,g.tss FROM
        (SELECT *, {snp_chr('src_id')} chr, {snp_pos('src_id')} pos
         FROM raw WHERE lower(edge_type)='cis-eqtl') r
        JOIN genes g ON r.dst_id=g.gene
        WHERE r.chr=g.gene_chr AND abs(r.pos-g.tss)<={args.span}"""])
    stage(con, 'snps', [f"""CREATE TABLE snps AS SELECT src_id snp,
        cast(floor((row_number() OVER (ORDER BY chr,pos,src_id)-1)/{args.batch_snps}.0) AS BIGINT) batch
        FROM (SELECT DISTINCT src_id,chr,pos FROM eq)"""])
    stage(con, 'ca', [f"""CREATE TABLE ca AS
        SELECT r.*,s.batch FROM
        (SELECT *, {snp_chr('src_id')} chr, {snp_pos('src_id')} pos,
         {acr_coord('dst_id','chr')} acr_chr, {acr_coord('dst_id','pos')} acr_pos
         FROM raw WHERE lower(edge_type)='cis-caqtl') r
        JOIN snps s ON r.src_id=s.snp
        WHERE r.chr=r.acr_chr AND abs(r.pos-r.acr_pos)<={args.span}
        ORDER BY s.batch"""])
    stage(con, 'ag', ["""CREATE TABLE ag AS
        SELECT * FROM raw WHERE regexp_matches(lower(edge_type),'aglink|distance')
        AND dst_id IN (SELECT DISTINCT dst_id FROM eq)"""])
    stage(con, 'co', [f"""CREATE TABLE co AS SELECT * FROM
        (SELECT *, {acr_coord('src_id','chr')} chr1, {acr_coord('src_id','pos')} pos1,
         {acr_coord('dst_id','chr')} chr2, {acr_coord('dst_id','pos')} pos2
         FROM raw WHERE lower(edge_type)='co-acr') r
        WHERE chr1=chr2 AND abs(pos1-pos2)<={args.span}
        AND src_id IN (SELECT DISTINCT dst_id FROM ca)
        AND dst_id IN (SELECT DISTINCT src_id FROM ag)"""])
    # ca and eq keep all evidence. No DISTINCT on SNP/ACR pairs here: duplicate
    # evidence is resolved only at the original unit-level deduplication step.
    log('Prepared disk tables: '+str({t: con.execute(f'SELECT count(*) FROM {t}').fetchone()[0]
                                    for t in ('eq','snps','ca','ag','co')}))


def build_batch(con, args, batch):
    span1 = 'greatest(c.pos,c.acr_pos,e.tss)-least(c.pos,c.acr_pos,e.tss)'
    span2 = 'greatest(c.pos,c.acr_pos,o.pos2,e.tss)-least(c.pos,c.acr_pos,o.pos2,e.tss)'
    pri = "CASE a.edge_type WHEN 'agLINK' THEN 2 WHEN 'distance' THEN 1 ELSE 0 END"
    support1 = "'caqtl_'||c.meta||';coACR_NA;eqtl_'||e.meta||';acr_gene_'||a.meta"
    support2 = "'caqtl_'||c.meta||';coACR_'||o.meta||';eqtl_'||e.meta||';acr_gene_'||a.meta"
    stage(con, f'batch_{batch}', [
        'DROP TABLE IF EXISTS batch_l1', 'DROP TABLE IF EXISTS batch_l2',
        f"""CREATE TABLE batch_l1 AS SELECT * EXCLUDE(rn) FROM (
          SELECT 'L1' AS unit_template,c.src_id snp,c.dst_id acr1,NULL::VARCHAR acr2,e.dst_id gene,
          {span1} span_bp,a.edge_type acr_gene_type,{support1} support,
          row_number() OVER (PARTITION BY c.src_id,e.dst_id,c.dst_id
              ORDER BY {pri} DESC,c.file_id,c.ord,a.file_id,a.ord,e.file_id,e.ord) rn
          FROM ca c JOIN eq e ON c.src_id=e.src_id
          JOIN ag a ON a.src_id=c.dst_id AND a.dst_id=e.dst_id
          WHERE c.batch={batch} AND {span1}<={args.span}) WHERE rn=1""",
        f"""CREATE TABLE batch_l2 AS SELECT * EXCLUDE(rn) FROM (
          SELECT 'L2' AS unit_template,c.src_id snp,c.dst_id acr1,o.dst_id acr2,e.dst_id gene,
          {span2} span_bp,a.edge_type acr_gene_type,{support2} support,
          row_number() OVER (PARTITION BY c.src_id,e.dst_id,c.dst_id,o.dst_id
              ORDER BY c.file_id,c.ord,o.file_id,o.ord,a.file_id,a.ord,e.file_id,e.ord) rn
          FROM ca c JOIN eq e ON c.src_id=e.src_id
          JOIN ag a ON a.dst_id=e.dst_id
          JOIN co o ON o.src_id=c.dst_id AND o.dst_id=a.src_id
          WHERE c.batch={batch} AND c.chr=o.chr2 AND {span2}<={args.span}
            AND NOT EXISTS (SELECT 1 FROM batch_l1 d WHERE d.snp=c.src_id
                AND d.gene=e.dst_id AND (d.acr1=c.dst_id OR d.acr1=o.dst_id))) WHERE rn=1"""
    ])


def export_batch(con, args, batch, folder):
    path = folder / f'01_part_{batch+1:06d}.tsv'
    key = f'part_{batch}'
    done = saved(con, key)
    if done:
        if identity(path) != done['file']:
            raise ValueError(f'Checkpoint part changed/missing: {path}')
        return done
    # A committed batch can only be reused while its tables still contain that
    # batch. Export follows build immediately, before the next batch overwrites.
    build_batch(con, args, batch)
    counts = dict(con.execute('SELECT unit_template,count(*) FROM (SELECT * FROM batch_l1 UNION ALL '
                              'SELECT * FROM batch_l2) GROUP BY unit_template').fetchall())
    partial = Path(str(path)+'.partial')
    with partial.open('w', newline='') as handle:
        writer = csv.writer(handle, delimiter='\t', lineterminator='\n')
        cur = con.execute('SELECT * FROM (SELECT * FROM batch_l1 UNION ALL SELECT * FROM batch_l2) '
                          'ORDER BY unit_template,snp,gene,acr1,acr2')
        while True:
            rows = cur.fetchmany(args.fetch_rows)
            if not rows:
                break
            writer.writerows([['NA' if v is None else v for v in row] for row in rows])
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(partial, path)
    done = dict(file=identity(path), counts=counts)
    save(con, key, done)
    con.execute('CHECKPOINT')
    return done


def write_json(path, value):
    temp = Path(str(path)+'.partial')
    temp.write_text(json.dumps(value, indent=2, ensure_ascii=False)+'\n')
    os.replace(temp, path)


def run(args):
    out = Path(args.out).resolve()
    db = Path(args.checkpoint_db or str(out.parent/'f01_local_checkpoint/01_local.duckdb')).resolve()
    qc = Path(args.qc or str(out.parent/'02_local_QC.tsv')).resolve()
    manifest = Path(args.manifest or str(out.parent/'03_local_manifest.json')).resolve()
    spill = Path(args.tmp_dir or str(db.parent/'f02_spill')).resolve()
    if args.span < 0 or min(args.threads,args.batch_snps,args.fetch_rows) < 1:
        raise ValueError('Invalid span/threads/batch size')
    if len({str(out),str(db),str(qc),str(manifest)}) != 4:
        raise ValueError('Output/checkpoint/QC/manifest paths must be distinct')
    sig = dict(version=VERSION, script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
               edges=[identity(p) for p in args.edges], gene_bed=identity(args.gene_bed),
               span=args.span, batch_snps=args.batch_snps, out=str(out),qc=str(qc),manifest=str(manifest))
    if db.exists() and not args.resume:
        raise FileExistsError(f'Checkpoint exists; use --resume: {db}')
    for path in (out.parent,db.parent,spill,qc.parent,manifest.parent):
        path.mkdir(parents=True,exist_ok=True)
    con = duckdb.connect(str(db))
    try:
        con.execute(f'SET memory_limit={lit(args.memory_limit)}')
        con.execute(f'SET threads={args.threads}')
        con.execute(f'SET temp_directory={lit(spill)}')
        con.execute('CREATE TABLE IF NOT EXISTS progress(key VARCHAR PRIMARY KEY,value VARCHAR)')
        previous = saved(con,'identity')
        if previous and previous != sig:
            raise ValueError('Inputs/script/semantic parameters changed; use a new checkpoint/output')
        if previous is None:
            for path in (out,qc,manifest):
                if path.exists():
                    raise FileExistsError(f'Refusing to overwrite existing output: {path}')
            save(con,'identity',sig)
            con.execute('CHECKPOINT')
        for key,path in [('output',out),('qc_output',qc),('manifest_output',manifest)]:
            done = saved(con,key)
            if done and identity(path) != done:
                raise ValueError(f'Completed output changed: {path}')
        if saved(con,'manifest_output'):
            log('All outputs completed; nothing to rerun')
            return
        prepare(con,args)
        batches = [r[0] for r in con.execute('SELECT DISTINCT batch FROM ca ORDER BY batch').fetchall()]
        parts = db.parent/'f01_export_parts'
        parts.mkdir(exist_ok=True)
        totals = {'L1':0,'L2':0}
        for index,batch in enumerate(batches,1):
            done = export_batch(con,args,batch,parts)
            for k,v in done['counts'].items():
                totals[k] += v
            log(f'Batch {index}/{len(batches)} completed; cumulative L1={totals["L1"]:,}; L2={totals["L2"]:,}')
        if saved(con,'output') is None:
            if out.exists():
                raise FileExistsError(f'Refusing to overwrite untracked output: {out}')
            partial = Path(str(out)+'.partial')
            with partial.open('w',newline='') as handle:
                writer = csv.writer(handle,delimiter='\t',lineterminator='\n')
                writer.writerow(HEADER)
                n = 0
                for batch in batches:
                    with (parts/f'01_part_{batch+1:06d}.tsv').open() as src:
                        for row in csv.reader(src,delimiter='\t'):
                            n += 1
                            writer.writerow([f'LU{n:09d}']+row)
                    log(f'Assembling output: batch {batch+1}; units={n:,}')
                handle.flush()
                os.fsync(handle.fileno())
            os.replace(partial,out)
            save(con,'output',identity(out))
            con.execute('CHECKPOINT')
        if saved(con,'qc_output') is None:
            partial = Path(str(qc)+'.partial')
            with partial.open('w',newline='') as handle:
                writer=csv.writer(handle,delimiter='\t',lineterminator='\n')
                writer.writerow(['metric','value'])
                for i,path in enumerate(args.edges):
                    for kind,n in con.execute(f'SELECT edge_type,count(*) FROM input_{i} GROUP BY edge_type').fetchall():
                        writer.writerow([f'input_{i+1}_retained_{kind}',n])
                for t in ('eq','ca','ag','co','snps'):
                    writer.writerow([f'prepared_{t}_rows',con.execute(f'SELECT count(*) FROM {t}').fetchone()[0]])
                writer.writerows([['L1_units',totals['L1']],['L2_units',totals['L2']],
                    ['total_units',sum(totals.values())],['completed_batches',len(batches)],
                    ['span_bp',args.span],['coordinate_rule','legacy SNP / ACR floor midpoint / gene BED strand TSS'],
                    ['tie_rule','agLINK priority for L1; equal-priority evidence uses first input order'],
                    ['peak_RSS_GiB',round(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024**2,3)]])
            os.replace(partial,qc)
            save(con,'qc_output',identity(qc))
        write_json(manifest,dict(inputs=sig,outputs=dict(network=identity(out),qc=identity(qc)),
            checkpoint_db=str(db),memory_limit=args.memory_limit,threads=args.threads,
            unit_order='batch,template,SNP,gene,ACR1,ACR2; unit IDs regenerated',
            biological_rules='unchanged legacy L1/L2; coACR direction preserved; both ACRs checked for direct L1',
            limitations='DuckDB memory limit is a buffer budget, not a hard process RSS limit; disk space required'))
        save(con,'manifest_output',identity(manifest))
        con.execute('CHECKPOINT')
        log(f'DONE: {sum(totals.values()):,} local units -> {out}')
    finally:
        con.close()


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--edges',nargs='+',required=True)
    p.add_argument('--gene-bed',required=True)
    p.add_argument('--out',required=True)
    p.add_argument('--span',type=int,default=500000)
    p.add_argument('--memory-limit',default='64GB')
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--batch-snps',type=int,default=10000)
    p.add_argument('--fetch-rows',type=int,default=10000)
    p.add_argument('--checkpoint-db')
    p.add_argument('--tmp-dir')
    p.add_argument('--qc')
    p.add_argument('--manifest')
    p.add_argument('--resume',action='store_true')
    p.add_argument('--verbose',action='store_true')
    p.add_argument('--report-edge-types',action='store_true')
    run(p.parse_args())


if __name__=='__main__':
    main()
