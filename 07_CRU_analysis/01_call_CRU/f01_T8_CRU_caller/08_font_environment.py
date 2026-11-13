#!/usr/bin/env python3
"""Bundled Liberation Sans plus a private fontconfig for each R process (no substitution)."""
import argparse
from contextlib import contextmanager
import csv
import hashlib
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import xml.etree.ElementTree as ET


HERE=Path(__file__).resolve().parent
FONTS=HERE/'f04_fonts'
NAMES=['LiberationSans-Regular.ttf','LiberationSans-Bold.ttf',
       'LiberationSans-Italic.ttf','LiberationSans-BoldItalic.ttf']


def md5(path):
    return hashlib.md5(Path(path).read_bytes()).hexdigest()


def bundle(source,license_path):
    FONTS.mkdir(exist_ok=True);rows=[]
    for name in NAMES:
        original=source/name;target=FONTS/name
        if target.exists() and md5(target)!=md5(original):raise ValueError('Different font already bundled: '+str(target))
        if not target.exists():shutil.copy2(original,target)
        rows.append(dict(source=str(original),file=name,md5=md5(target)))
    target=FONTS/'02_copyright.txt'
    if not target.exists():shutil.copy2(license_path,target)
    with (FONTS/'01_source_manifest.tsv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=['source','file','md5'],delimiter='\t');w.writeheader();w.writerows(rows)
    print('Bundled Liberation Sans regular/bold/italic/bold-italic and license:',FONTS)


def signature():
    return {p.name:md5(p) for p in [Path(__file__).resolve()]+[FONTS/n for n in NAMES]}


def validate_family(family):
    with plotting_environment() as env:
        row=subprocess.check_output(['fc-match','-f','%{family}\t%{file}',family],env=env,text=True)
    fields=row.split('\t',1)
    if len(fields)!=2 or family.lower() not in fields[0].lower() or not Path(fields[1]).is_file():
        raise ValueError('Requested plot font is unavailable: '+family+'; fontconfig returned '+row)
    return row


@contextmanager
def plotting_environment(base=None,include_system=True):
    for name in NAMES:
        if not (FONTS/name).is_file():raise ValueError('Missing bundled font: '+str(FONTS/name))
    with tempfile.TemporaryDirectory(prefix='t8_fontconfig_') as tmp:
        root=ET.Element('fontconfig')
        if include_system:
            ET.SubElement(root,'include',{'ignore_missing':'yes'}).text='/etc/fonts/fonts.conf'
        ET.SubElement(root,'dir').text=str(FONTS)
        ET.SubElement(root,'cachedir').text=str(Path(tmp)/'cache')
        config=Path(tmp)/'fonts.conf'
        ET.ElementTree(root).write(config,encoding='utf-8',xml_declaration=True)
        env=dict(os.environ if base is None else base)
        env['FONTCONFIG_FILE']=str(config)
        if shutil.which('fc-match',path=env.get('PATH','')) is None:
            for directory in [Path(sys.executable).resolve().parent,Path('/usr/bin'),Path('/bin')]:
                tool=directory/'fc-match'
                if tool.is_file() and os.access(tool,os.X_OK):
                    env['PATH']=str(directory)+os.pathsep+env.get('PATH','');break
        env.pop('R_LIBS_USER',None)
        yield env


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--bundle-from',type=Path,required=True)
    p.add_argument('--license',type=Path,required=True)
    a=p.parse_args();bundle(a.bundle_from,a.license)
