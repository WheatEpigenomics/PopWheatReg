#!/usr/bin/env python3
"""V2 simple Route-count score, with compact provenance in the human-facing main table.

Keep the legacy score DESC / n_support DESC order and all original TFN IDs.
Only qtl_methods is tokenized/deduplicated; qtl_p_source is untouched and
qtl_locus_id is omitted. Original network and full evidence remain read-only.
"""
import argparse
from functools import lru_cache
import os
from pathlib import Path
import re
import tempfile
import time

import duckdb
import pandas as pd


METHOD_ORDER=('Global','PW','sg1','sg2','sg3','sg5')


def log(message):
    print('[score_compact '+time.strftime('%H:%M:%S')+'] '+message,flush=True)


def lit(value):
    return "'"+str(value).replace("'","''")+"'"


@lru_cache(maxsize=8192)
def compact_methods(value):
    tokens={s.strip() for s in re.split(r'[;|]',value)}-{'','NA'}
    known=[m for m in METHOD_ORDER if m in tokens]
    # Preserve unexpected method labels instead of silently discarding them.
    return ';'.join(known+sorted(tokens-set(METHOD_ORDER))) or 'NA'


def calculate_chunk_score(frame):
    # Same legacy definition: number of Route labels, not number of evidence paths.
    routes=frame['routes'].fillna('')
    return routes.str.count(r'\|').add(1).where(routes!='',0).astype('int64')


def run(args):
    source=Path(args.input).resolve(); target=Path(args.output).resolve()
    if source==target:
        raise ValueError('Input network must not be overwritten')
    if args.chunk_rows<1 or args.threads<1:
        raise ValueError('Chunk rows and threads must be positive')
    columns=pd.read_csv(source,sep='\t',nrows=0).columns.tolist()
    if not {'routes','n_support','net_id'}.issubset(columns):
        raise ValueError('Required network columns: net_id, routes, n_support')
    target.parent.mkdir(parents=True,exist_ok=True)
    # A failed rerun must not truncate the previous scored table.
    with tempfile.TemporaryDirectory(prefix='score_compact_',dir=target.parent) as tmp:
        chunks=Path(tmp)/'01_scored_chunks.tsv'
        result=Path(tmp)/'02_sorted.tsv'
        rows=0; header=True
        for frame in pd.read_csv(source,sep='\t',dtype=str,na_filter=False,
                usecols=lambda col:col!='qtl_locus_id',chunksize=args.chunk_rows):
            frame['score']=calculate_chunk_score(frame)
            if 'qtl_methods' in frame:
                frame['qtl_methods']=frame['qtl_methods'].map(compact_methods)
            # Do not split, normalize, or otherwise rewrite qtl_p_source.
            frame.to_csv(chunks,sep='\t',index=False,mode='w' if header else 'a',header=header)
            header=False; rows+=len(frame)
            log(f'Processed {rows:,} rows')
        if header:
            fields=[c for c in columns if c not in ('qtl_locus_id','score')]+['score']
            pd.DataFrame(columns=fields).to_csv(chunks,sep='\t',index=False)
        con=duckdb.connect()
        try:
            con.execute(f'SET memory_limit={lit(args.ram_limit)}')
            con.execute(f'SET threads={args.threads}')
            con.execute(f'SET temp_directory={lit(Path(tmp)/"f01_spill")}')
            log('Sort by score DESC, n_support DESC (legacy order)')
            con.execute(f"""COPY (
                SELECT * FROM read_csv({lit(chunks)},delim='\t',header=true,
                    auto_detect=true,all_varchar=true,nullstr='')
                ORDER BY TRY_CAST(score AS INTEGER) DESC,TRY_CAST(n_support AS BIGINT) DESC
            ) TO {lit(result)} (HEADER,DELIMITER '\t')""")
        finally:
            con.close()
        os.replace(result,target)
    log(f'DONE: {rows:,} rows; {target}')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input',required=True)
    p.add_argument('--output',required=True)
    p.add_argument('--ram-limit',default='64GB')
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--chunk-rows',type=int,default=500000)
    run(p.parse_args())


if __name__=='__main__':
    main()
