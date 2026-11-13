#!/usr/bin/env python3
"""Stream native-TFN evidence parts to TXT, merge/gzip, verify, then optionally delete parts.

Only the explicit new C-agLINK run is touched. The main network, manifest, QC,
checkpoints and all other directories are preserved. No production task is started
by importing this module. A small state file permits safe conversion/cleanup resume.
"""
import argparse
import csv
import fcntl
import gzip
import hashlib
import json
import os
import re
from pathlib import Path

import pyarrow.parquet as pq


def log(message):
    print('[evidence_text] '+message, flush=True)


def identity(path):
    path=Path(path).resolve(); st=path.stat()
    return dict(path=str(path),bytes=st.st_size,mtime_ns=st.st_mtime_ns)


def atomic_json(path,value):
    temp=Path(str(path)+'.partial')
    temp.write_text(json.dumps(value,indent=2,ensure_ascii=False)+'\n')
    os.replace(temp,path)


def verified_file(info):
    if identity(info['path'])!=info:
        raise ValueError('Tracked file changed: '+info['path'])


def archive_hash(path):
    digest=hashlib.sha256()
    # Reading to EOF also validates every gzip member's CRC and truncated output.
    with gzip.open(path,'rb') as h:
        for chunk in iter(lambda:h.read(8*1024*1024),b''):
            digest.update(chunk)
    return digest.hexdigest()


def run(args):
    prefix=Path(args.prefix).resolve()
    out=Path(args.out).resolve(); qc=Path(args.qc).resolve()
    manifest=Path(str(prefix)+'.manifest.json')
    main=Path(str(prefix)+'.network.tsv')
    if out.parent!=prefix.parent or qc.parent!=prefix.parent:
        raise ValueError('Archive and QC must be in the explicit network run directory')
    if args.batch_rows<1 or not 0<=args.compresslevel<=9:
        raise ValueError('Invalid batch size or gzip level')
    cfg=json.loads(manifest.read_text())
    if cfg.get('execution_fingerprint',{}).get('version')!='C-agLINK-v2.2':
        raise ValueError('Expected a completed fresh C-agLINK V2.2 manifest; old runs are not modified')
    state_path=Path(str(out)+'.checkpoint.json')
    with open(str(out)+'.lock','a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if state_path.exists():
            state=json.loads(state_path.read_text())
            if state['out']!=str(out) or state['qc']!=str(qc) or state['prefix']!=str(prefix):
                raise ValueError('Conversion state belongs to a different output')
            verified_file(state['manifest']); verified_file(state['main'])
        else:
            if out.exists() or qc.exists():
                raise FileExistsError('Untracked archive/QC exists; choose a new output name')
            files=list(prefix.parent.glob(prefix.name+'.evidence.part*.parquet'))
            pattern=re.compile(re.escape(prefix.name)+r'\.evidence\.part([1-9][0-9]*)\.parquet')
            if not files or any(not pattern.fullmatch(p.name) or p.is_symlink() for p in files):
                raise ValueError('No parts, nonstandard part name, or symlink part')
            files.sort(key=lambda p:int(pattern.fullmatch(p.name).group(1)))
            numbers=[int(pattern.fullmatch(p.name).group(1)) for p in files]
            n_main=cfg['execution']['final_counts'][1]
            step=cfg['parameters']['evidence_part_net_ids']
            if numbers!=list(range(1,max(1,(n_main+step-1)//step)+1)):
                raise ValueError('Missing/extra evidence parts')
            parts=[]; columns=None
            for path in files:
                pf=pq.ParquetFile(path)
                if columns is None: columns=pf.schema_arrow.names
                if columns!=pf.schema_arrow.names or 'net_id' not in columns:
                    raise ValueError('Evidence schemas differ or native TFN column is absent')
                text=path.with_suffix('.txt')
                if text.exists():
                    raise FileExistsError('Untracked TXT exists: '+str(text))
                parts.append(dict(source=identity(path),text=str(text),rows=pf.metadata.num_rows))
            expected=cfg['execution']['final_counts'][0]
            if sum(p['rows'] for p in parts)!=expected:
                raise ValueError('Parquet footer row total differs from completed network manifest')
            state=dict(prefix=str(prefix),out=str(out),qc=str(qc),manifest=identity(manifest),
                       main=identity(main),parts=parts,columns=columns,expected_rows=expected)
            atomic_json(state_path,state)
        # Resume must not adopt newly added parts or delete paths outside this prefix.
        allowed={p['source']['path'] for p in state['parts']}
        actual={str(p.resolve()) for p in prefix.parent.glob(prefix.name+'.evidence.part*.parquet')}
        if actual-allowed:
            raise ValueError('Additional untracked evidence parts appeared')
        for part in state['parts']:
            source=Path(part['source']['path']); text=Path(part['text'])
            if source.parent!=prefix.parent or text!=source.with_suffix('.txt'):
                raise ValueError('Unsafe part path in conversion state')
            if source.exists(): verified_file(part['source'])
            elif 'archive' not in state:
                raise FileNotFoundError('Source disappeared before archive validation: '+str(source))
            if 'converted' in part and text.exists(): verified_file(part['converted'])
        if 'archive' not in state:
            for number,part in enumerate(state['parts'],1):
                if 'converted' in part:
                    verified_file(part['converted'])
                    log(f'TXT {number}/{len(state["parts"])} already verified; skip')
                    continue
                source=Path(part['source']['path']); text=Path(part['text'])
                if text.exists(): raise FileExistsError('Uncommitted TXT exists: '+str(text))
                partial=Path(str(text)+'.partial')
                rows=0
                log(f'Convert TXT {number}/{len(state["parts"])}: {source.name}')
                with open(partial,'w',newline='') as h:
                    w=csv.writer(h,delimiter='\t',lineterminator='\n')
                    w.writerow(state['columns'])
                    pf=pq.ParquetFile(source)
                    for batch in pf.iter_batches(batch_size=args.batch_rows):
                        values=batch.to_pydict()
                        for record in zip(*(values[c] for c in state['columns'])):
                            w.writerow(['NA' if v is None else v for v in record])
                            rows+=1
                verified_file(part['source'])
                if rows!=part['rows']: raise ValueError('TXT row count differs from Parquet footer')
                os.replace(partial,text)
                part['converted']=identity(text); atomic_json(state_path,state)
                log(f'TXT {number} complete; records={rows:,}')
            if out.exists():
                raise FileExistsError('Uncommitted archive exists; preserve it and choose a new output')
            digest=hashlib.sha256(); partial=Path(str(out)+'.partial')
            header=None
            with gzip.open(partial,'wb',compresslevel=args.compresslevel) as target:
                for number,part in enumerate(state['parts'],1):
                    verified_file(part['converted'])
                    with open(part['text'],'rb') as source:
                        first=source.readline()
                        if header is None:
                            header=first; target.write(first); digest.update(first)
                        elif first!=header:
                            raise ValueError('TXT headers differ')
                        for chunk in iter(lambda:source.read(8*1024*1024),b''):
                            target.write(chunk); digest.update(chunk)
                    log(f'Merged {number}/{len(state["parts"])} TXT parts')
            log('Validate gzip CRC and SHA256 of the complete uncompressed merged content')
            if archive_hash(partial)!=digest.hexdigest():
                raise ValueError('Merged gzip content validation failed; no source part is deleted')
            os.replace(partial,out)
            state['archive']=identity(out); state['uncompressed_sha256']=digest.hexdigest()
            atomic_json(state_path,state)
        verified_file(state['archive'])
        if args.delete_parts:
            # Validation succeeded. Delete only this manifest's exact part files,
            # never glob-delete, touch the old run, or remove checkpoint batches.
            for part in state['parts']:
                for info in (part['source'],part['converted']):
                    path=Path(info['path'])
                    if path.exists():
                        verified_file(info); path.unlink()
                part['deleted']=True
                atomic_json(state_path,state)
            log('Validated archive retained; this run\'s published Parquet/TXT parts deleted')
        with open(str(qc)+'.partial','w',newline='') as h:
            w=csv.writer(h,delimiter='\t',lineterminator='\n')
            w.writerow(['item','value'])
            w.writerows([['evidence_records',state['expected_rows']],['parts',len(state['parts'])],
                         ['gzip_crc_and_content_verified','TRUE'],
                         ['uncompressed_sha256',state['uncompressed_sha256']],
                         ['parts_deleted',all(p.get('deleted',False) for p in state['parts'])],
                         ['archive',out],['archive_bytes',state['archive']['bytes']]])
        os.replace(str(qc)+'.partial',qc)
        log(f'DONE: {out}; records={state["expected_rows"]:,}')


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('prefix','out','qc'): p.add_argument('--'+key,required=True)
    p.add_argument('--batch-rows',type=int,default=256)
    p.add_argument('--compresslevel',type=int,default=1)
    p.add_argument('--delete-parts',action='store_true')
    run(p.parse_args())


if __name__=='__main__':
    main()
