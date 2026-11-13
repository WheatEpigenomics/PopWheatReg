#!/usr/bin/env python3
"""Generate a tiny synthetic dataset for the demo.

No real data is used. The generator writes:

    demo/data/genotype.parquet   variants x samples (int8, -9 = missing)
    demo/data/phenotype.tsv      phenotypes x samples (inverse-normal transformed)
    demo/data/covariates.tsv     covariates x samples (3 PCs, 0-based layout)
    demo/data/subgroups.tsv      headerless: sample_id <TAB> subgroup

A few effects are planted so that both `global_pc` and `subgroup_omnibus`
produce non-empty output. Run `run_demo.sh` for the end-to-end example.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import ndtri

SEED = 20240101
N_VARIANTS = 200
N_PHENO = 60
GROUPS = ["sg1", "sg2", "sg3", "sg5"]
N_PER_GROUP = 40
N_PCS = 3


def main():
    here = Path(__file__).resolve().parent
    data = here / "data"
    data.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)
    n_samples = len(GROUPS) * N_PER_GROUP
    samples = [f"{g}_S{i:02d}" for g in GROUPS for i in range(N_PER_GROUP)]
    group_of = np.repeat(np.arange(len(GROUPS)), N_PER_GROUP)

    # ---------------- genotypes ----------------
    G = rng.binomial(2, 0.30, size=(N_VARIANTS, n_samples)).astype(np.int8)
    G[rng.random(G.shape) < 0.02] = -9           # 2% missing
    G[3, group_of == 2] = 0                      # monomorphic in sg3 -> K_est < Kmax

    # ---------------- phenotypes (INT) ----------------
    Y = rng.normal(size=(N_PHENO, n_samples)).astype(np.float64)
    planted = [
        # (variant index, phenotype index, {group: beta})
        (0, 0, {"sg1": 0.90}),                                   # one strong subgroup
        (1, 1, {"sg1": 0.45, "sg2": 0.45, "sg3": 0.45, "sg5": 0.45}),  # all subgroups
        (2, 2, {"sg1": 0.55, "sg2": 0.55, "sg3": 0.55}),         # multiple moderate
        (3, 3, {"sg1": 0.80, "sg2": 0.80, "sg5": 0.80}),         # sg3 monomorphic
    ]
    for vi, pi, effects in planted:
        g = G[vi].astype(np.float64).copy()
        g[g == -9] = 1.0                          # demo-only fill for simulation
        for gi, lab in enumerate(GROUPS):
            beta = effects.get(lab, 0.0)
            if beta:
                idx = group_of == gi
                Y[pi, idx] += beta * g[idx]

    # rank-based inverse-normal transform, per phenotype
    for pi in range(N_PHENO):
        rank = np.argsort(np.argsort(Y[pi]))
        Y[pi] = ndtri((rank + 0.5) / n_samples)

    # ---------------- covariates (3 PCs) ----------------
    C = rng.normal(size=(N_PCS, n_samples)).astype(np.float32)

    # ---------------- write tables ----------------
    variant_ids = [f"chr1:{1000 + i * 1000}" for i in range(N_VARIANTS)]
    pheno_ids = [f"acr{i:04d}" for i in range(N_PHENO)]

    pd.DataFrame(G, index=variant_ids, columns=samples).to_parquet(data / "genotype.parquet")
    pd.DataFrame(Y.astype(np.float32), index=pheno_ids, columns=samples).to_csv(
        data / "phenotype.tsv", sep="\t", float_format="%.6g")
    pd.DataFrame(C, index=[f"PC{i + 1}" for i in range(N_PCS)], columns=samples).to_csv(
        data / "covariates.tsv", sep="\t", float_format="%.6g")

    with (data / "subgroups.tsv").open("w") as fh:
        for sample, gi in zip(samples, group_of):
            fh.write(f"{sample}\t{GROUPS[gi]}\n")

    print(f"demo data written to {data}")
    print(f"  genotypes : {N_VARIANTS} variants x {n_samples} samples "
          f"(2% missing, variant 4 monomorphic in sg3)")
    print(f"  phenotypes: {N_PHENO} x {n_samples} (INT)")
    print(f"  covariates: {N_PCS} PCs")
    print(f"  subgroups : {', '.join(f'{g}={N_PER_GROUP}' for g in GROUPS)}")
    print("  planted pairs (variant, phenotype, subgroup betas):")
    for vi, pi, effects in planted:
        print(f"    {variant_ids[vi]}  x  {pheno_ids[pi]}  {effects}")


if __name__ == "__main__":
    main()
