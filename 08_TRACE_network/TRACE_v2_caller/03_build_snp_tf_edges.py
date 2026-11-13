#!/usr/bin/env python3
"""V2 SNP-TF: original Local units + all non-intergenic ANN; no ACR-overlap source."""
import argparse
import csv
import importlib.util
from collections import Counter, defaultdict
from pathlib import Path

spec = importlib.util.spec_from_file_location('trace_inputs', Path(__file__).with_name('01_prepare_edge_tables.py'))
common = importlib.util.module_from_spec(spec)
spec.loader.exec_module(common)


def proximal(path, tfs, old=False):
    pairs = defaultdict(lambda: {'effects': set(), 'impacts': set()})
    counts = Counter()
    with common.text_open(path) as h:
        reader = csv.DictReader(h, delimiter='\t')
        prefix = common.ann_fields(reader.fieldnames)
        previous = None
        # allANN export is VCF/ANN ordered; count SNPs without storing all annotations.
        for row in reader:
            counts['annotation_rows'] += 1
            if counts['annotation_rows'] % 1000000 == 0:
                print(f"[proximal] read {counts['annotation_rows']:,} ANN rows", flush=True)
            if row['ID'] != previous:
                counts['ordered_snp_records'] += 1
                previous = row['ID']
            effect = row[prefix + 'EFFECT']
            if effect == 'intergenic_region':
                continue
            if old:
                match = common.GENE.search(row[prefix + 'GENE'])
                tf = match.group(0) if match else 'NA'  # reproduce historical ANN[0] matching for audit
            else:
                tf = common.annotation_gene(row, prefix)
                if 'canonical_gene_id' in row and row['canonical_gene_id'] != tf:
                    raise ValueError('Stale/incorrect canonical_gene_id: ' + row['ID'])
                if tf == 'NA':
                    counts['unresolved_nonintergenic_annotations'] += 1
            if tf in tfs:
                entry = pairs[(row['ID'], tf)]
                entry['effects'].add(effect)
                entry['impacts'].add(row[prefix + 'IMPACT'] or 'NA')
    counts.update(unique_proximal_pairs=len(pairs), unique_proximal_snps=len({s for s, t in pairs}),
                  unique_proximal_tfs=len({t for s, t in pairs}))
    return pairs, counts


def build(args):
    with common.text_open(args.tf_list) as h:
        tfs = {r['TF'] for r in csv.DictReader(h, delimiter='\t')}
    prox, counts = proximal(args.snpeff, tfs)
    if args.old_snpeff:
        old, _ = proximal(args.old_snpeff, tfs, old=True)
        lost = old.keys() - prox.keys()
        added = prox.keys() - old.keys()
        counts.update(old_proximal_pairs=len(old), added_proximal_pairs=len(added), lost_old_pairs=len(lost))
        counts['added_over_old_fraction'] = len(added) / len(old) if old else 'NA'
        for effect, n in Counter(e for k in added for e in prox[k]['effects']).items():
            counts['added_pair_effect:' + effect] = n
        if lost:
            raise ValueError(f'Old valid Proximal pairs lost: {len(lost)}, examples={list(lost)[:5]}')
    sources = defaultdict(set)
    meta = defaultdict(set)
    for key, entry in prox.items():
        sources[key].add('tf_prox_snpeff')
        meta[key].add('EFFECT=' + '|'.join(sorted(entry['effects'])) + ';IMPACT=' +
                      '|'.join(sorted(entry['impacts'])) + ';local_path=NA')
    if args.local_units:
        with common.text_open(args.local_units) as h:
            for row in csv.DictReader(h, delimiter='\t'):
                if row['gene'] not in tfs:
                    continue
                key = (row['snp'], row['gene'])
                sources[key].add('local_unit')
                path = 'direct' if row['acr2'] in ('', 'NA', 'NaN', 'nan') else 'coACR'
                meta[key].add(row['support'] + ';local_path=' + path)
                counts['local_unit_support_rows'] += 1
    with common.output(args.out) as h:
        w = csv.writer(h, delimiter='\t', lineterminator='\n')
        w.writerow(common.EDGE_HEADER)
        for (snp, tf), support in sorted(sources.items()):
            w.writerow(['SNP', snp, 'TF', tf, 'TF-SNP', '|'.join(sorted(support)),
                        ';'.join(sorted(meta[(snp, tf)]))])
            counts['source:' + '|'.join(sorted(support))] += 1
    counts.update(final_snp_tf_pairs=len(sources), final_snps=len({s for s, t in sources}),
                  final_tfs=len({t for s, t in sources}))
    common.write_qc(args.qc, counts)
    print(dict(counts), flush=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--local-units', help='Required for formal output; omit only for isolated Proximal audit')
    p.add_argument('--proximal-only', action='store_true')
    p.add_argument('--snpeff', required=True)
    p.add_argument('--old-snpeff', help='Optional ANN[0] regression reference, never used to add edges')
    p.add_argument('--tf-list', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--qc', required=True)
    args = p.parse_args()
    if not args.local_units and not args.proximal_only:
        p.error('Formal SNP-TF output requires --local-units; audit only: --proximal-only')
    build(args)


if __name__ == '__main__':
    main()
