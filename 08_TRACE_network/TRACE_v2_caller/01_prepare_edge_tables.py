#!/usr/bin/env python3
"""V2 input adapters: retain all rows and essential QTL provenance; no QTL selection."""
import argparse
import csv
import gzip
import math
import re
from collections import Counter
from pathlib import Path

EDGE_HEADER = ['src_type', 'src_id', 'dst_type', 'dst_id', 'edge_type', 'score', 'meta']
CAQTL_PROVENANCE_FIELDS = ('methods', 'p_source', 'locus_id')
GENE = re.compile(r'TraesCS[0-9A-Z]+G[0-9]+(?![0-9])')


def text_open(path, mode='rt'):
    return gzip.open(path, mode, newline='') if str(path).endswith('.gz') else open(path, mode, newline='')


def output(path):
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    return text_open(p, 'xt')


def canonical_gene(*values):
    # Full match extraction, never split on '-' (intergenic annotations contain two IDs).
    ids = {m.group(0) for value in values for m in GENE.finditer(value or '')}
    return next(iter(ids)) if len(ids) == 1 else 'NA'


def ann_fields(names):
    prefix = 'ANN.' if 'ANN.GENE' in names else 'ANN[0].'
    required = ['ID', prefix + 'GENE', prefix + 'EFFECT', prefix + 'IMPACT']
    if not set(required) <= set(names):
        raise ValueError('Missing ANN fields: ' + str(required))
    return prefix


def annotation_gene(row, prefix):
    effect = row[prefix + 'EFFECT']
    if effect == 'intergenic_region':
        return 'NA'
    return canonical_gene(row[prefix + 'GENE'], row.get(prefix + 'GENEID', ''))


def write_qc(path, counts):
    with output(path) as h:
        w = csv.writer(h, delimiter='\t', lineterminator='\n')
        w.writerow(['metric', 'value'])
        w.writerows(sorted(counts.items()))


def normalize(args):
    counts = Counter()
    with text_open(args.input) as src, output(args.out) as dst:
        r = csv.DictReader(src, delimiter='\t')
        prefix = ann_fields(r.fieldnames)
        names = list(r.fieldnames)
        if 'canonical_gene_id' not in names:
            names.append('canonical_gene_id')
        w = csv.DictWriter(dst, names, delimiter='\t', lineterminator='\n')
        w.writeheader()
        for row in r:
            row['canonical_gene_id'] = annotation_gene(row, prefix)
            w.writerow(row)
            counts['annotation_rows'] += 1
            counts['canonical_' + ('missing_or_intergenic' if row['canonical_gene_id'] == 'NA' else 'resolved')] += 1
            if counts['annotation_rows'] % 1000000 == 0:
                print(f"[allANN] {counts['annotation_rows']:,} rows", flush=True)
    write_qc(args.qc, counts)


def edges(args):
    counts = Counter()
    expected_rows = None
    if args.preparation_summary:
        # Optional upstream QC table; when given, the row count must match exactly.
        with open(args.preparation_summary) as h:
            prep = list(csv.DictReader(h, delimiter='\t'))
        if args.kind == 'caqtl':
            items = {r['item']: r['value'] for r in prep if r['scope'] == 'all'}
            if items.get('rows_preserved') != 'True':
                raise ValueError('caQTL beta preparation has not completed with rows_preserved=True')
            expected_rows = int(items['output_rows'])
        else:
            items = {r['item']: r['value'] for r in prep}
            expected_rows = int(items['query_rows'])
    genes = {}
    with open(args.gene_bed) as h:
        for fields in csv.reader(h, delimiter='\t'):
            # The checked fil5 producer uses BED start (NOT strand-aware TSS).
            # Preserve its input cis/trans contract; Local's strand-aware TSS stays unchanged separately.
            genes[fields[3]] = (fields[0], int(fields[1]))
    with text_open(args.input) as src, output(args.out) as dst:
        w = csv.writer(dst, delimiter='\t', lineterminator='\n')
        w.writerow(EDGE_HEADER)
        if args.kind == 'eqtl':
            reader = csv.reader(src, delimiter='\t')
        else:
            reader = csv.DictReader(src, delimiter='\t')
            required = {'variant_id', 'phenotype_id', 'pvalue', 'cis_trans', 'beta', 'beta_se',
                        'methods', 'p_source', 'locus_id'}
            if not required <= set(reader.fieldnames):
                raise ValueError('caQTL must use the SOURCE-ANNOTATED beta output, not the 9-column BED')
        for n, row in enumerate(reader, 1):
            if args.kind == 'eqtl':
                if len(row) != 9:
                    raise ValueError(f'eQTL line {n}: expected 9 columns, got {len(row)}')
                chrom, start, end, snp, target, p, label, beta, se = row
                extra = {}
            else:
                snp, target, p, label, beta, se = [row[k] for k in
                    ('variant_id', 'phenotype_id', 'pvalue', 'cis_trans', 'beta', 'beta_se')]
                # beta/se already match p_source. Keep multi-method support and
                # source locus, not raw lookup diagnostics or other SG effects.
                # The full annotated input is retained unchanged for backtracking.
                extra = {k: row[k] for k in CAQTL_PROVENANCE_FIELDS}
            if label not in ('cis', 'tra'):
                raise ValueError(f'Invalid cis/tra label at line {n}: {label}')
            m = re.fullmatch(r'([^-:_]+)[-:_]([0-9]+)', snp)
            if not m:
                raise ValueError('Invalid SNP coordinate ID: ' + snp)
            snp_chr, snp_pos = m.group(1), int(m.group(2))
            if args.kind == 'eqtl':
                if target not in genes:
                    raise ValueError('QTL target missing from gene BED: ' + target)
                target_chr, target_pos = genes[target]
                input_target_distance = abs(snp_pos-target_pos)
            else:
                ac = re.fullmatch(r'[^:]+:([^_]+)_([0-9]+)_([0-9]+)', target)
                if not ac:
                    raise ValueError('Invalid ACR coordinate ID: ' + target)
                target_chr = ac.group(1)
                # Match union_outputs_legacy_20261016.classify exactly (do not reclassify using midpoint).
                input_target_distance = max(int(ac.group(2))-snp_pos,snp_pos-int(ac.group(3)),0)
            if label == 'tra' and snp_chr == target_chr and input_target_distance <= 1000000:
                raise ValueError(f'Trans label violates existing >1 Mb input rule at row {n}: {snp}, {target}; inspect upstream annotation')
            pv = float(p)
            if not math.isfinite(pv) or not 0 <= pv <= 1:
                raise ValueError(f'Invalid P at line {n}: {p}')
            meta = {'p': p, 'beta': beta or 'NA', 'beta_se': se or 'NA', **extra}
            if any(';' in str(v) or '=' in str(v) for v in meta.values()):
                # methods commonly uses ';'; encode separators without changing methods/p_source semantics.
                meta = {k: str(v).replace('%', '%25').replace(';', '%3B').replace('=', '%3D') for k, v in meta.items()}
            etype = ('cis-' if label == 'cis' else 'trans-') + ('eQTL' if args.kind == 'eqtl' else 'caQTL')
            w.writerow(['SNP', snp, 'Gene' if args.kind == 'eqtl' else 'ACR', target,
                        etype, 'inf' if pv == 0 else f'{-math.log10(pv):.12g}',
                        ';'.join(f'{k}={v}' for k, v in meta.items())])
            counts[etype] += 1
            counts['input_rows'] += 1
            counts['output_rows'] += 1
            if n % 1000000 == 0:
                print(f'[{args.kind} edges] {n:,} rows', flush=True)
    # A trans-only candidate table must not silently become Local input.
    if expected_rows is not None and counts['output_rows'] != expected_rows:
        raise ValueError(f'Incomplete/wrong QTL input: expected {expected_rows}, read {counts["output_rows"]}')
    if counts['cis-' + ('eQTL' if args.kind == 'eqtl' else 'caQTL')] == 0:
        raise ValueError('No cis QTL! Do not build Local from a trans-only input. Inspect this output before use.')
    write_qc(args.qc, counts)


def aglink(args):
    """agLINK and distance ACR-gene pairs -> ACR->Gene edge table."""
    counts = Counter()
    with text_open(args.input) as src, output(args.out) as dst:
        w = csv.writer(dst, delimiter='\t', lineterminator='\n')
        w.writerow(EDGE_HEADER)
        for n, line in enumerate(src, 1):
            row = line.rstrip('\n').split('\t')
            if len(row) < 5:
                raise ValueError(f'agLINK line {n}: expected gene, ACR, value, p, type; got {len(row)} columns')
            gene, acr, value, _, kind = row[:5]
            if kind == 'dist':
                w.writerow(['ACR', acr, 'Gene', gene, 'distance', 'NA', 'dist=' + value])
                counts['distance'] += 1
                counts['output_rows'] += 1
            elif kind in ('aglink', 'calink'):
                w.writerow(['ACR', acr, 'Gene', gene, 'agLINK', value, 'PCC=' + value])
                counts['agLINK'] += 1
                counts['output_rows'] += 1
            else:
                raise ValueError(f'agLINK line {n}: unknown pair type {kind!r}; expected aglink, calink or dist')
            counts['input_rows'] += 1
    write_qc(args.qc, counts)


def coacr(args):
    """ACR-ACR co-accessibility -> ACR->ACR edge table."""
    counts = Counter()
    with text_open(args.input) as src, output(args.out) as dst:
        w = csv.writer(dst, delimiter='\t', lineterminator='\n')
        w.writerow(EDGE_HEADER)
        for n, line in enumerate(src, 1):
            row = line.rstrip('\n').split('\t')
            if n == 1 and row and row[0] == 'id1':
                continue
            if len(row) < 3:
                raise ValueError(f'coACR line {n}: expected id1, id2, correlation; got {len(row)} columns')
            a, b, corr = row[:3]
            w.writerow(['ACR', a, 'ACR', b, 'co-ACR', corr, 'corr=' + corr])
            counts['input_rows'] += 1
            counts['output_rows'] += 1
    write_qc(args.qc, counts)


def binding(args):
    """TF-motif-ACR correlation -> TF->ACR binding table, signed |PCC| >= cutoff."""
    counts = Counter()
    with text_open(args.input) as src, output(args.out) as dst:
        w = csv.writer(dst, delimiter='\t', lineterminator='\n')
        w.writerow(EDGE_HEADER)
        for n, line in enumerate(src, 1):
            row = line.rstrip('\n').split('\t')
            if n == 1 and row and row[0].lower() in ('tf', 'tf_id'):
                continue
            if len(row) < 4:
                raise ValueError(f'binding line {n}: expected TF, motif, ACR, PCC[, effect]; got {len(row)} columns')
            counts['input_rows'] += 1
            tf, motif, acr, pcc = row[:4]
            effect = row[4] if len(row) > 4 and row[4] else 'NA'
            value = float(pcc)
            if not math.isfinite(value):
                counts['skipped_non_finite'] += 1
                continue
            if abs(value) < args.min_abs_pcc:
                counts['skipped_below_abs_pcc'] += 1
                continue
            if args.exclude_un and ':Un_' in acr:
                counts['skipped_un'] += 1
                continue
            score = pcc[1:] if pcc.startswith('-') else pcc
            w.writerow(['TF', tf, 'ACR', acr, 'bind_coacc', score,
                        f'PCC={pcc};effect={effect};motif={motif}'])
            counts['output_rows'] += 1
    write_qc(args.qc, counts)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    sub = p.add_subparsers(dest='mode', required=True)
    for mode in ('normalize-ann', 'edges', 'aglink', 'coacr', 'binding'):
        q = sub.add_parser(mode)
        q.add_argument('--input', required=True)
        q.add_argument('--out', required=True)
        q.add_argument('--qc', required=True)
        if mode == 'binding':
            q.add_argument('--min-abs-pcc', type=float, default=0.35,
                           help='keep signed correlations with abs(PCC) >= this cutoff (default 0.35)')
            q.add_argument('--keep-un', dest='exclude_un', action='store_false',
                           help='keep ACRs on the Un fraction (dropped by default, matching the published table)')
        if mode == 'edges':
            q.add_argument('--kind', choices=['eqtl', 'caqtl'], required=True)
            q.add_argument('--gene-bed', required=True, help='Check existing trans input rule: gene BED start (fil5 eQTL) / interval (six-route caQTL union)')
            q.add_argument('--preparation-summary', help='Optional upstream QC table; when given, rows must match exactly')
    args = p.parse_args()
    {'normalize-ann': normalize, 'edges': edges, 'aglink': aglink,
     'coacr': coacr, 'binding': binding}[args.mode](args)


if __name__ == '__main__':
    main()
