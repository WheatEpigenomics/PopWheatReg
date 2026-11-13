#!/usr/bin/env python3
import gzip
import numpy as np
import matplotlib.pyplot as plt

def compute_ld_decay(ld_file, bin_size=10000, cutoff=0.2, out_prefix="ld_decay"):
    sum_r2 = {}
    count = {}

    with gzip.open(ld_file, "rt") as f:
        header = next(f)  # skip header
        for line in f:
            chr_a, bp_a, snp_a, chr_b, bp_b, snp_b, r2 = line.strip().split()
            if chr_a != chr_b:  # only intra-chr
                continue
            bp_a, bp_b, r2 = int(bp_a), int(bp_b), float(r2)
            d = abs(bp_b - bp_a)
            bin = d // bin_size
            sum_r2[bin] = sum_r2.get(bin, 0.0) + r2
            count[bin] = count.get(bin, 0) + 1

    bins = sorted(sum_r2.keys())
    x = np.array([b * bin_size for b in bins])
    y = np.array([sum_r2[b] / count[b] for b in bins])

    # 输出结果
    with open(f"{out_prefix}.tsv", "w") as out:
        out.write("distance_bp\tmean_r2\tcount\n")
        for xi, yi, ci in zip(x, y, [count[b] for b in bins]):
            out.write(f"{xi}\t{yi}\t{ci}\n")

    # 画图
    plt.figure(figsize=(6,4))
    plt.plot(x, y, marker="o", markersize=3, lw=1)
    plt.axhline(y=cutoff, color="r", ls="--")
    plt.xlabel("Distance (bp)")
    plt.ylabel("Mean r²")
    plt.title("LD decay")
    plt.tight_layout()
    plt.savefig(f"{out_prefix}.png", dpi=300)

    # 估算r²=0.2时的距离
    idx = np.where(y <= cutoff)[0]
    if len(idx) > 0:
        d_cut = x[idx[0]]
        print(f"r² drops below {cutoff} at ~{d_cut/1000:.1f} kb")
    else:
        print(f"r² never drops below {cutoff}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) < 2:
        print("Usage: python ld_decay.py input.ld.gz [bin_size] [cutoff]")
        sys.exit(1)
    ld_file = sys.argv[1]
    bin_size = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
    cutoff = float(sys.argv[3]) if len(sys.argv) > 3 else 0.2
    compute_ld_decay(ld_file, bin_size, cutoff)

