# Population simulations

Three independent simulation and regression analyses. Each subdirectory has its
own `README.md`, its own `00_RUN.sh` and a `demo/` dataset.

| Module | What it does | Entry point |
|---|---|---|
| `01_Multivariate_regression` | nested and SNP-included models of how ACR, gene expression and SNP density explain the phenotype | `00_RUN.sh` |
| `02_Neutral_Wright_Fisher_simulations` | neutral forward simulations matched to the observed LD and block structure | `00_RUN.sh` |
| `03_Forward-simulation` | forward simulations of the per-block PCC top 1% selection model | `00_RUN.sh` |

Every runner takes the same shape:

```sh
bash 00_RUN.sh                    # bundled demo -> output_demo/
PROJECT_ROOT=/path/to/project bash 00_RUN.sh project
INPUT_ROOT=/path/to/data     bash 00_RUN.sh full
```

`demo` needs no external data; `project` and `full` need the real inputs and stop
with an explicit message if a required directory or environment variable is
missing.

## Requirements and runtime

Python 3 with numpy (and matplotlib for the forward simulation), R with ggplot2
and Cairo. Demo runs finish in seconds to a couple of minutes; full runs are
documented in each module README.
