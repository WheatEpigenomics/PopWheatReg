# m16_junctions: shared-junction export
# Frozen rule module, copied unchanged from the published calling chain
# (source: f21_T7A_pipeline/09_make_junctions.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
#!/usr/bin/env python3
"""Derive the shared-CRU junction table used by f17_V11_pipeline/07_block_enrichment.py.

The f17 block test needs a boundaries file whose positions restricted to the CRU
endpoints equal the "internal" set (starts that are also ends). T7-A has no
calling-style boundary table of its own, so we rebuild exactly that set from the
final BED: a junction exists where CRU_i ends exactly where CRU_{i+1} starts.
Using this file cannot change any CRU coordinate; it only tells the block test
which endpoints are shared cuts.
"""
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main(args):
    sizes = {}
    with open(args.genome_sizes) as f:
        for line in f:
            if line.strip():
                c, n = line.split()[:2]
                sizes[c] = int(n)
    by_chrom = defaultdict(list)
    with open(args.cru_bed) as f:
        for line in f:
            if not line.strip() or line.startswith('#'):
                continue
            p = line.split()
            by_chrom[p[0]].append((int(p[1]), int(p[2])))
    rows = []
    for c in args.chroms:
        rr = sorted(by_chrom.get(c, []))
        if not rr:
            continue
        starts = {a for a, b in rr}
        ends = {b for a, b in rr}
        for pos in sorted(starts & ends):
            if 0 < pos < sizes[c]:
                rows.append({'chrom': c, 'boundary_pos': pos, 'source': 'shared_CRU_junction',
                             'parent_start': '', 'parent_end': '', 'method': args.label, 'accepted': 1})
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, 'w', newline='') as f:
        w = csv.DictWriter(f, list(rows[0]) if rows else ['chrom', 'boundary_pos'], delimiter='\t',
                           lineterminator='\n')
        w.writeheader()
        w.writerows(rows)
    print('junctions', len(rows), '->', args.out)


if __name__ == '__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cru-bed', type=Path, required=True)
    p.add_argument('--genome-sizes', type=Path, required=True)
    p.add_argument('--chroms', nargs='+', required=True)
    p.add_argument('--label', default='T7')
    p.add_argument('--out', type=Path, required=True)
    main(p.parse_args())
