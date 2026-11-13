# Neutral Wright–Fisher simulations

Run the bundled real-data demo from any working directory:

```bash
bash 00_RUN.sh                         # demo -> output_demo/
bash 00_RUN.sh demo /path/to/output    # optional output directory
INPUT_ROOT=/path/to/prepared_data bash 00_RUN.sh full
```

The demo uses chromosome 1A: 293 real SNPs, 5 real agLinks, and 16 ACRs.
It runs the same **30 replicates, 1,000 matched background sets, N=203**, and
G10/G20/G50/G100/G200 checkpoints as the full analysis, followed by statistics
and R plotting. One of the five demo pairs meets the original matching criteria;
unsupported pairs/strata are recorded explicitly. It demonstrates execution on
a small subset and is not a replacement for the full chromosome analysis.

Full/project mode requires `INPUT_ROOT` (`SOURCE_RUN` is accepted as an alias).
The input directory must contain `02_prepared/`, `02_background/`,
`genes.gff3`, and `observed_blocks.bed`, following the bundled `demo/` layout.
The `[input]` section in `settings.ini` sets the gene/block files, relative to
`INPUT_ROOT` or as absolute paths. No external `config_all21.ini` is required.
The full chromosome list and scientific settings remain in `settings.ini`.
Missing inputs are reported before compilation or simulation.

- `00_RUN.sh`: demo/full entry point; defaults to demo.
- `01_run_neutral.py`: input preparation, simulation orchestration and statistics.
- `02_forward_neutral.cpp`: original simulation model.
- `03_plot_results.R`: independent plotting script.
- `settings.ini`: input paths and analysis parameters.
- `demo/`: self-contained real input subset; NPZ files contain named arrays.

Requirements: Python >=3.8 with NumPy >=1.20, g++/OpenMP, PLINK 1.9,
and R with Cairo. Install these dependencies before running; set `PYTHON_BIN`,
`PLINK_BIN`, `RSCRIPT_BIN`, or `CXX` to select installed executables.
`JOBS` controls concurrent trajectories (default 1).
Results are in `04_summary/`, figures in `05_figures/`, and checks in
`VALIDATION.json` under the output directory. Matching runs resume; changed
inputs/code/settings require a new output directory.

## Runtime

The bundled demo runs one chromosome and finishes in a few minutes. A full run
covers all 21 chromosomes with `JOBS` parallel workers; the wall-clock scales
roughly with `21 / JOBS` times the per-chromosome time.
