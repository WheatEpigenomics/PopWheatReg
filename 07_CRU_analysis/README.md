# CRU analysis

Two analyses of cis-regulatory units (CRUs):

| Module | What it does | Entry point |
|---|---|---|
| `01_call_CRU` | build CRUs from the agLink and co-accessibility networks, then score them | `call_cru_main.sh` |
| `02_CRU_perturbation` | sentinel perturbation analysis inside fixed CRUs | `00_RUN.sh` |

## `01_call_CRU`

The code lives in `01_call_CRU/f01_T8_CRU_caller/` with its own README; the
entry point is `01_call_CRU/call_cru_main.sh`, organised as copy-paste blocks:

1. calling (default: no figures, no block test)
2. haplotype-block permutation test
3. 5-Mb region figures (optional, needs R)
4. scoring with the frozen T7-A anchors
5. score figure (optional)
6. demo

```sh
cd 01_call_CRU
bash ./f01_T8_CRU_caller/demo/run_demo.sh     # toy genome, a few seconds
```

Output: `out_cru/f04_merged/CRU.bed` and `out_score/04_CRU_scores_simple.tsv`.

A failed figure is reported as `PLOT FAILED` and never blocks calling, merging or
the block test; redraw later with `--mode plot --plot all --resume`.

## `02_CRU_perturbation`

Independent of the caller: it reads prepared CRU, ACR signal and gene expression
matrices. `bash 00_RUN.sh demo` runs the bundled 8-CRU dataset and passes 147
internal checks.

## Requirements and runtime

Python 3 with numpy and scipy; R with ggplot2 and Cairo for figures only; CUDA is
not needed.

The bundled demos finish in under a minute each. A full run is dominated by the
per-chromosome calling step, which parallelises over chromosomes; scoring and the
block test are small additions on top of it.
