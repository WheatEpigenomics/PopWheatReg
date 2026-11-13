#!/usr/bin/env python3
"""Optional independent PDF page counting; reports "unavailable" instead of guessing."""
import argparse
import csv
import importlib.util
import os
from pathlib import Path
import shutil
import subprocess
import sys
import warnings


def find_pdfinfo(rscript=None):
    if os.environ.get('CRU_PDFINFO')=='none':return None
    candidates=[os.environ.get('CRU_PDFINFO'),shutil.which('pdfinfo')]
    if rscript:candidates.append(str(Path(rscript).resolve().parent/'pdfinfo'))
    candidates.extend(['/usr/bin/pdfinfo','/bin/pdfinfo'])
    for p in candidates:
        if p and Path(p).is_file() and os.access(p,os.X_OK):return str(Path(p).resolve())
    return None


def check_pdf(path,expected,rscript=None):
    path=Path(path)
    # These structural checks are not a substitute for independent page counting.
    with path.open('rb') as f:
        if not f.read(8).startswith(b'%PDF-'):raise ValueError(f'Not a PDF: {path}')
        f.seek(max(0,path.stat().st_size-1024))
        if b'%%EOF' not in f.read():raise ValueError(f'Incomplete PDF (no EOF): {path}')
    pages=None;backend='unavailable';tool=find_pdfinfo(rscript)
    if tool:
        result=subprocess.run([tool,str(path)],capture_output=True,text=True,check=True,
                              env=dict(os.environ,LC_ALL='C'))
        rows=[s for s in result.stdout.splitlines() if s.startswith('Pages:')]
        if len(rows)!=1:raise ValueError('Cannot parse pdfinfo page count: '+str(path))
        pages=int(rows[0].split(':',1)[1].strip());backend=tool
    else:
        for name in ['pypdf','PyPDF2','fitz']:
            if importlib.util.find_spec(name) is None:continue
            if name=='fitz':
                import fitz
                with fitz.open(str(path)) as doc:pages=len(doc)
            else:
                module=__import__(name)
                reader=getattr(module,'PdfReader',None) or getattr(module,'PdfFileReader')
                with path.open('rb') as f:pages=len(reader(f).pages)
            backend=name;break
    if pages is not None and pages!=expected:
        raise ValueError(f'PDF page count mismatch: {path}; expected={expected}, actual={pages}')
    status='verified' if pages is not None else 'unavailable'
    if pages is None:
        warnings.warn(f'No pdfinfo/pypdf/PyPDF2/PyMuPDF available: page count NOT independently verified for {path}. '
                      f'Expected {expected} pages; see render QC. Calling is unaffected.',RuntimeWarning)
    return dict(pdf=str(path),expected_pages=expected,actual_pages=pages,page_count_status=status,page_count_backend=backend)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--pdf',type=Path,required=True);p.add_argument('--expected',type=int,required=True)
    p.add_argument('--rscript');a=p.parse_args()
    row=check_pdf(a.pdf,a.expected,a.rscript)
    w=csv.DictWriter(sys.stdout,fieldnames=list(row),delimiter='\t',lineterminator='\n')
    w.writeheader();w.writerow({k:'NA' if v is None else v for k,v in row.items()})
