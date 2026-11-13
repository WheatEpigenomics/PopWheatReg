# ACR-gene link

ACR-gene association (agLink), ACR co-accessibility, and the SMR/HEIDI test that
connects caQTL to eQTL. Four scripts; `00_aglink_main.sh` holds the runnable
blocks and a demo.

## Scripts

| File | Role |
|---|---|
| `01_calcu_acr_gene_pcc.py` | Pearson correlation between one gene expression profile and each candidate ACR profile |
| `02_fill_column_for_ztest.sh` | pads a control row to a fixed number of columns |
| `03_aglink_ztest.py` | one-step z-test of the real PCC against the permutation null |
| `04_co_access_analysis.py` | ACR-ACR co-accessibility within a window, GPU version |

## Pipeline

1. Candidate ACRs within +/-1 Mb of each gene TSS (`bedtools intersect`).
2. Real PCC for every gene-ACR pair (`01`).
3. Null distribution: 10,000 ACRs drawn from other chromosomes per gene, same PCC
   (`01`), padded to a fixed width (`02`).
4. z-test: `z = (|r_real| - |mean_null|) / sd_null`, one-sided `p = 1 - Phi(z)` (`03`).
5. Significant links: `|r| > 0.4` and `p < 0.05`.
6. Co-accessibility between ACRs within 1 Mb (`04`, needs CUDA).
7. SMR/HEIDI on caQTL as exposure and eQTL as outcome, keeping `p_SMR < 0.05`,
   `p_HEIDI > 0.05` and `nsnp_HEIDI >= 3`.

## Inputs

| File | Format |
|---|---|
| gene expression | genes x samples, space-delimited, no header |
| ACR signal | ACRs x samples, first column ACR id |
| gene BED | BED6 |
| ACR BED | BED4 |
| co-accessibility matrix | header `ACR_ID` then sample columns |
| `.smr` result | columns 23/24/25 = `p_SMR`, `p_HEIDI`, `nsnp_HEIDI` |

## Demo

`demo/` holds a four-ACR toy matrix and a five-ACR co-accessibility set; the demo
block in `00_aglink_main.sh` runs all four scripts in a few seconds.

## Requirements and runtime

Python 3 with numpy, scipy and pandas; bedtools and shuf; SMR 1.4.1 with a PLINK
bed for the final block; CUDA only for `04`.

The demo runs all four scripts in under a minute; a full run scales with the
number of genes and candidate ACRs.
