# m17_block_reference_null: M3 local null generator
# Frozen rule module, copied unchanged from the published calling chain
# (source: f23_block_reference/01_M3_local_null.py). Only the module-reference literals were rewritten to
# relative imports so the whole rule set ships as one package.
"""Verbatim M3 block readers/index and local point null; validation only."""
from collections import defaultdict


def load_genome_sizes(path):
    sizes = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            sizes[parts[0]] = int(float(parts[1]))
    return sizes


def load_block_bed(path, chroms):
    if not path:
        return {}
    wanted = set(chroms)
    by_chrom = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] not in wanted:
                continue
            start = int(float(parts[1]))
            end = int(float(parts[2]))
            if end <= start:
                continue
            by_chrom[parts[0]].append({
                "chrom": parts[0],
                "start": start,
                "end": end,
                "block_id": (
                    parts[3] if len(parts) >= 4
                    else f"{parts[0]}_BLOCK_{len(by_chrom[parts[0]]) + 1}"
                ),
            })
    for chrom in by_chrom:
        by_chrom[chrom].sort(key=lambda x: (x["start"], x["end"]))
    return by_chrom


def block_test_index(blocks_by_chrom, chroms):
    index = {}
    for chrom in chroms:
        blocks = sorted(
            blocks_by_chrom.get(chrom, []),
            key=lambda x: (x["start"], x["end"]),
        )
        boundaries = sorted({
            position
            for block in blocks
            for position in (block["start"], block["end"])
        })
        index[chrom] = {
            "blocks": blocks,
            "starts": [block["start"] for block in blocks],
            "boundaries": boundaries,
        }
    return index


def local_random_points(points, genome_sizes, window_bp, rng):
    randomized = []
    for row in points:
        chrom = row["chrom"]
        original = row["pos"]
        lo = max(1, original - window_bp)
        hi = min(genome_sizes[chrom] - 1, original + window_bp)
        exclusion = 5_000
        ranges = []
        if lo <= original - exclusion:
            ranges.append((lo, original - exclusion))
        if original + exclusion <= hi:
            ranges.append((original + exclusion, hi))
        widths = [right - left + 1 for left, right in ranges]
        if not widths or sum(widths) <= 0:
            continue
        ticket = rng.randrange(sum(widths))
        chosen = None
        for (left, right), width in zip(ranges, widths):
            if ticket < width:
                chosen = left + ticket
                break
            ticket -= width
        randomized.append({"chrom": chrom, "pos": chosen})
    return randomized
