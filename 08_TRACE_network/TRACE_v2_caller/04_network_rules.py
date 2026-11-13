#!/usr/bin/env python3
"""Shared TRACE V2 rule library: routes, distance flags, ACR hub and network assembly.

Route1/2 and the ACR hub keep their original semantics. The C-agLINK builder imports these
functions unchanged; this module is a library, not an entry point.
"""
import csv
from pathlib import Path

import duckdb


def lit(x):
    return "'" + str(x).replace("'", "''") + "'"


def read(path, header=True):
    return (f"read_csv({lit(path)}, delim='\\t', header={str(header).lower()}, all_varchar=true, "
            "auto_detect=true, quote='', nullstr='NA')")


def connection(args):
    con = duckdb.connect()
    con.execute(f'SET memory_limit={lit(args.ram_limit)}')
    con.execute(f'SET threads={args.threads}')
    if args.tmp_dir:
        Path(args.tmp_dir).mkdir(parents=True, exist_ok=True)
        con.execute(f'SET temp_directory={lit(args.tmp_dir)}')
    con.execute('PRAGMA enable_progress_bar')
    return con


def inputs(con, args):
    con.execute(f"""CREATE TABLE gene_loc AS SELECT column0 chr, cast(column1 AS BIGINT) start0,
        cast(column2 AS BIGINT) finish, column3 gene_id,
        (cast(column1 AS BIGINT)+cast(column2 AS BIGINT))/2 center
        FROM {read(args.gene_bed, False)}""")
    if con.execute('SELECT count(*)-count(DISTINCT gene_id) FROM gene_loc').fetchone()[0]:
        raise ValueError('Duplicate gene IDs in gene BED')
    con.execute(f"""CREATE TABLE tf_snp AS SELECT trim(src_id) src_id, trim(dst_id) dst_id,
        score tf_snp_type, meta FROM {read(args.tf_snp)}""")
    if con.execute('SELECT count(*)-count(DISTINCT (src_id,dst_id)) FROM tf_snp').fetchone()[0]:
        raise ValueError('Duplicate SNP-TF pairs; rebuild with 06_build_TF_SNP_edges.py')
    bad = con.execute("""SELECT * FROM tf_snp WHERE tf_snp_type IS NULL OR
        NOT regexp_full_match(tf_snp_type, '(local_unit|tf_prox_snpeff)(\\|(local_unit|tf_prox_snpeff))?') LIMIT 1""").fetchone()
    if bad:
        raise ValueError(f'Unexpected/obsolete SNP-TF source: {bad}')
    for name in ('eqtl', 'caqtl', 'aglink', 'coacr'):
        if getattr(args, name, None):
            if name in ('eqtl','caqtl'):
                # One scan of each new QTL EDGE input. Retain only SNPs that can join the network;
                # no candidate/P filtering is added. Local used the COMPLETE cis+trans edge tables.
                con.execute(f'''CREATE TABLE raw_{name} AS SELECT q.* FROM {read(getattr(args,name))} q
                    WHERE EXISTS (SELECT 1 FROM tf_snp t WHERE t.src_id=trim(q.src_id))''')
            else:
                con.execute(f'CREATE VIEW raw_{name} AS SELECT * FROM {read(getattr(args, name))}')
    if getattr(args, 'binding', None):
        # Formal input must remain signed legacy |PCC|>=0.35; don't trust unsigned score.
        con.execute(f'CREATE TABLE raw_binding AS SELECT * FROM {read(args.binding)}')
        bad = con.execute("""SELECT src_id,dst_id,meta FROM raw_binding
          WHERE try_cast(regexp_extract(meta, '(^|;)PCC=([^;]+)', 2) AS DOUBLE) IS NULL
             OR NOT isfinite(try_cast(regexp_extract(meta, '(^|;)PCC=([^;]+)', 2) AS DOUBLE))
             OR abs(try_cast(regexp_extract(meta, '(^|;)PCC=([^;]+)', 2) AS DOUBLE)) < 0.35 LIMIT 1""").fetchone()
        if bad:
            raise ValueError(f'Formal binding input is not the signed |PCC|>=0.35 table: {bad}')


def acr_hub(con, direct_only=False):
    # Original v4 hub: prefer direct agLINK over distance; coACR only without an existing direct pair.
    con.execute("""CREATE TABLE acr_gene_direct AS
        SELECT acr1,gene,link_type,gene_meta FROM (
          SELECT src_id acr1,dst_id gene,edge_type link_type,meta gene_meta,
            row_number() OVER(PARTITION BY src_id,dst_id ORDER BY
               CASE WHEN edge_type='agLINK' THEN 0 ELSE 1 END) rn
          FROM raw_aglink WHERE edge_type IN ('agLINK','distance')) WHERE rn=1""")
    if direct_only:
        return
    con.execute("""CREATE TABLE acr_gene_hub AS
        SELECT acr1,NULL::VARCHAR acr2,gene,link_type,'direct' path_class,
               NULL::VARCHAR coacr_meta,gene_meta FROM acr_gene_direct
        UNION ALL
        SELECT c.src_id,c.dst_id,d.gene,d.link_type,'coacr',c.meta,d.gene_meta
        FROM raw_coacr c JOIN acr_gene_direct d ON c.dst_id=d.acr1
        WHERE NOT EXISTS (SELECT 1 FROM acr_gene_direct x WHERE x.acr1=c.src_id AND x.gene=d.gene)""")


def binding_route(con, binding='raw_binding', dest='route_binding'):
    con.execute(f"""CREATE VIEW {dest} AS
        SELECT t.src_id snp,b.src_id tf,h.gene,t.tf_snp_type,'binding' route,
          h.acr1,NULL::VARCHAR acr2,'none' coacr_label,'aglink' link_label,
          concat_ws(';',replace(replace(b.meta,'PCC=','tf_acr_pcc='),'effect=','legacy_PCC_sign_label='),
                    'PCC_interpretation=correlation_only',concat('aglink_',h.gene_meta)) support,
          t.meta tf_snp_meta
        FROM {binding} b JOIN tf_snp t ON trim(b.src_id)=t.dst_id
        JOIN acr_gene_direct h ON b.dst_id=h.acr1 WHERE h.link_type='agLINK'""")


def routes(con):
    # No change to the Route1/2 topology, direct/coACR preference, or trans labels.
    con.execute("""CREATE VIEW route_eqtl AS
        SELECT e.src_id snp,t.dst_id tf,e.dst_id gene,t.tf_snp_type,'eqtl_trans' route,
          NULL::VARCHAR acr1,NULL::VARCHAR acr2,'none' coacr_label,'NA' link_label,
          concat('eqtl_',e.meta) support,t.meta tf_snp_meta
        FROM raw_eqtl e JOIN tf_snp t ON trim(e.src_id)=t.src_id WHERE e.edge_type LIKE 'trans%'""")
    con.execute("""CREATE VIEW route_caqtl AS
        SELECT c.src_id snp,t.dst_id tf,h.gene,t.tf_snp_type,'caqtl_trans' route,h.acr1,h.acr2,
          CASE WHEN h.path_class='coacr' THEN 'all' ELSE 'none' END coacr_label,
          CASE WHEN h.link_type='agLINK' THEN 'aglink' ELSE 'dist' END link_label,
          concat_ws(';',concat('caqtl_',c.meta),
            CASE WHEN h.path_class='coacr' THEN concat('coacr_',h.coacr_meta) END,
            concat(CASE WHEN h.link_type='agLINK' THEN 'aglink_' ELSE 'acr_gene_' END,h.gene_meta)) support,
          t.meta tf_snp_meta
        FROM raw_caqtl c JOIN tf_snp t ON trim(c.src_id)=t.src_id
        JOIN acr_gene_hub h ON c.dst_id=h.acr1 WHERE c.edge_type LIKE 'trans%'""")
    binding_route(con)
    con.execute("""CREATE TABLE all_evidence AS SELECT DISTINCT * FROM (
        SELECT * FROM route_eqtl UNION ALL SELECT * FROM route_caqtl UNION ALL SELECT * FROM route_binding)""")


def flags(con, source, dest, distance=500000):
    # TF/ACR center distances are intentionally identical to v4. For SNP-gene use
    # 1-based SNP p and the inclusive gene bases [BED.start+1, BED.end].
    con.execute(f"""CREATE TABLE {dest} AS WITH located AS (
        SELECT e.*,tf.chr tf_chr,g.chr gene_chr,tf.center tf_center,g.center gene_center,
          g.start0 gene_start0,g.finish gene_end1,
          regexp_extract(e.snp,'^([^-:_]+)[-:_]([0-9]+)$',1) snp_chr,
          try_cast(regexp_extract(e.snp,'^([^-:_]+)[-:_]([0-9]+)$',2) AS BIGINT) snp_pos1,
          split_part(split_part(e.acr1,':',2),'_',1) acr_chr,
          (try_cast(split_part(split_part(e.acr1,':',2),'_',2) AS BIGINT)+
           try_cast(split_part(split_part(e.acr1,':',2),'_',3) AS BIGINT))/2 acr_center
        FROM {source} e LEFT JOIN gene_loc tf ON e.tf=tf.gene_id
                           LEFT JOIN gene_loc g ON e.gene=g.gene_id), distances AS (
        SELECT *,CASE WHEN tf_chr=gene_chr THEN abs(tf_center-gene_center) END tf_gene_center_distance_bp,
          CASE WHEN tf_chr=acr_chr THEN abs(tf_center-acr_center) END tf_acr_center_distance_bp,
          CASE WHEN snp_chr=gene_chr THEN greatest(gene_start0+1-snp_pos1,snp_pos1-gene_end1,0) END snp_gene_interval_distance_bp
        FROM located)
        SELECT *, (tf_chr IS NOT NULL AND gene_chr IS NOT NULL AND snp_pos1 IS NOT NULL
            AND (acr1 IS NULL OR (acr_chr!='' AND acr_center IS NOT NULL))) pass_coordinates,
          tf!=gene pass_no_self,
          coalesce(tf_chr!=gene_chr OR tf_gene_center_distance_bp>{distance},false) pass_tf_gene,
          (acr1 IS NULL OR coalesce(tf_chr!=acr_chr OR tf_acr_center_distance_bp>{distance},false)) pass_tf_acr,
          coalesce(snp_chr!=gene_chr OR snp_gene_interval_distance_bp>{distance},false) pass_snp_gene
        FROM distances""")


PASS = 'pass_coordinates AND pass_no_self AND pass_tf_gene AND pass_tf_acr AND pass_snp_gene'


def counts(con, table, where='true', weight_column=None):
    rows = 'count(*)' if weight_column is None else f'coalesce(sum({weight_column}),0)::BIGINT'
    return con.execute(f"""SELECT {rows},count(DISTINCT(snp,tf,gene)),
      count(DISTINCT(tf,gene)),count(DISTINCT tf) FROM {table} WHERE {where}""").fetchone()


def filter_qc(con, path, weight_column=None):
    fields = ['route','mode','rule','evidence_rows','snp_tf_gene_relations','tf_gene_pairs','TFs',
              'fraction_of_route_evidence','fraction_of_route_snp_relations','fraction_of_route_tf_gene_pairs',
              'fully_lost_snp_tf_gene_relations','fully_lost_tf_gene_pairs',
              'fully_lost_snp_relation_fraction','fully_lost_tf_gene_fraction']
    with open(path, 'x') as h:
        w = csv.writer(h, delimiter='\t', lineterminator='\n'); w.writerow(fields)
        for route in ['ALL','eqtl_trans','caqtl_trans','binding']:
            scope = 'true' if route == 'ALL' else f'route={lit(route)}'
            base = counts(con, 'flagged', scope, weight_column)
            def emit(mode, rule, pred, before=None, after=None):
                values = counts(con, 'flagged', f'({scope}) AND ({pred})', weight_column)
                fractions = [v/b if b else 0 for v,b in zip(values[:3],base[:3])]
                lost = ['NA','NA']
                lost_fractions = ['NA','NA']
                if before is not None:
                    lost = [con.execute(f"""SELECT count(*) FROM (
                        SELECT DISTINCT {key} FROM flagged WHERE ({scope}) AND ({before}) EXCEPT
                        SELECT DISTINCT {key} FROM flagged WHERE ({scope}) AND ({after}))""").fetchone()[0]
                            for key in ('snp,tf,gene','tf,gene')]
                    lost_fractions = [v/b if b else 0 for v,b in zip(lost,base[1:3])]
                w.writerow([route, mode, rule, *values, *fractions, *lost, *lost_fractions])
            emit('baseline','before_final_filters','true')
            preceding = 'true'
            for flag in ('pass_coordinates','pass_no_self','pass_tf_gene','pass_tf_acr','pass_snp_gene'):
                emit('independent_rejected',flag,f'NOT {flag}','true',flag)
                emit('sequential_rejected',flag,f'({preceding}) AND NOT {flag}',preceding,preceding+' AND '+flag)
                preceding += f' AND {flag}'
                emit('remaining_after',flag,preceding)
            # New SNP-gene rule's extra loss against the old center-only pruning, excluding self separately.
            emit('additional_loss','snp_gene_vs_center_only',
                 'pass_coordinates AND pass_tf_gene AND pass_tf_acr AND NOT pass_snp_gene',
                 'pass_coordinates AND pass_tf_gene AND pass_tf_acr',
                 'pass_coordinates AND pass_tf_gene AND pass_tf_acr AND pass_snp_gene')
            emit('final','retained',PASS)


def save_network(con, prefix, write_outputs=True, low_memory=False):
    # Append explicit provenance columns to old evidence schema; support retains ALL original QTL metadata.
    def value(key):
        raw = f"nullif(regexp_extract(support,'(^|;)(?:eqtl_|caqtl_)?{key}=([^;]+)',2),'')"
        return f"replace(replace(replace({raw},'%3B',';'),'%3D','='),'%25','%')"
    extra = ','.join(value(k)+' AS '+name for k,name in
        [('methods','qtl_methods'),('p_source','qtl_p_source'),('locus_id','qtl_locus_id'),
         ('p','qtl_p'),('beta','qtl_beta'),('beta_se','qtl_beta_se')])
    con.execute(f'CREATE TABLE filtered_evidence AS SELECT *,{extra} FROM flagged WHERE {PASS}')
    def aggregate(field):
        # Avoid one ordered string-aggregation sort state per network group in
        # DuckDB 0.9.2. DISTINCT list uses hash aggregation; sort only each final
        # short unique list. array_to_string ignores NULL, as string_agg does.
        if low_memory:
            return f"array_to_string(list_sort(list(DISTINCT {field})),'|')"
        return f"string_agg(DISTINCT {field},'|' ORDER BY {field})"
    con.execute(f"""CREATE TABLE network_main AS
      SELECT row_number() OVER(ORDER BY snp,tf,gene) net_id,snp,tf,gene,
        {aggregate('tf_snp_type')} tf_snp_type,
        {aggregate('route')} routes,count(*) n_support,
        CASE WHEN bool_and(coacr_label='all') THEN 'all' WHEN bool_and(coacr_label='none') THEN 'none' ELSE 'mix' END coacr_type,
        CASE WHEN bool_and(link_label='NA') THEN 'NA'
          WHEN bool_and(link_label IN ('aglink','NA')) AND bool_or(link_label='aglink') THEN 'aglink'
          WHEN bool_and(link_label IN ('dist','NA')) AND bool_or(link_label='dist') THEN 'dist' ELSE 'mix' END link_type,
        {aggregate('qtl_methods')} qtl_methods,
        {aggregate('qtl_p_source')} qtl_p_source,
        {aggregate('qtl_locus_id')} qtl_locus_id
      FROM filtered_evidence GROUP BY snp,tf,gene""")
    if not write_outputs:
        return
    con.execute(f"COPY network_main TO {lit(prefix+'.network.tsv')} (HEADER,DELIMITER '\t',NULL 'NA')")
    # Original first 11 evidence columns remain unchanged.
    con.execute("""CREATE VIEW evidence_final AS SELECT n.net_id,e.snp,e.tf,e.gene,e.tf_snp_type,e.route,
      coalesce(e.acr1,'NA') acr1,coalesce(e.acr2,'NA') acr2,e.coacr_label,e.link_label,e.support,
      e.tf_snp_meta,e.qtl_methods,e.qtl_p_source,e.qtl_locus_id,e.qtl_p,e.qtl_beta,e.qtl_beta_se,
      e.tf_gene_center_distance_bp,e.tf_acr_center_distance_bp,e.snp_gene_interval_distance_bp
      FROM filtered_evidence e JOIN network_main n USING(snp,tf,gene)""")
    n = con.execute('SELECT count(*) FROM network_main').fetchone()[0]
    for part,start in enumerate(range(1,max(n,1)+1,200000),1):
        con.execute(f"COPY (SELECT * FROM evidence_final WHERE net_id BETWEEN {start} AND {start+199999}) "
                    f"TO {lit(prefix+f'.evidence.part{part}.parquet')} (FORMAT PARQUET)")
