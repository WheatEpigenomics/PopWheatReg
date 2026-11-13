# CRU perturbation analysis

Data and algorithms come from `csv2_T8_CRU_20261109`. Original grouping rules,
distance matching, statistical methods and random seeds are preserved.

- `00_RUN.sh`: run the analysis.
- `01_prepare_cru.py`: select sentinels, match controls and calculate low/high effects.
- `02_analyze_cru.R`: calculate statistics, draw figures and validate results.
- `demo/`: prepared real data for 8 CRUs and all 200 samples, plus reference results.

Requirements: Bash, standard Unix tools and R >= 3.6. Raw-data processing also
requires Python >= 3.8 and NumPy >= 1.24 (`python3 -m pip install 'numpy>=1.24'`).

Run from this directory:

```bash
# Default: bundled real-data demo -> output_demo/
bash 00_RUN.sh

# Reuse a prepared project (data/ and docs/)
PROJECT_ROOT=/path/to/prepared_project bash 00_RUN.sh project

# Recompute from complete raw matrices
INPUT_ROOT=/path/to/01_data bash 00_RUN.sh full
```

Each mode accepts an optional output directory, e.g. `bash 00_RUN.sh demo /path/to/output`.
`PROJECT_ROOT` is required for project mode; `INPUT_ROOT` is required for full mode.
Project inputs include `data/cru_input.bed`; if absent, provide `INPUT_ROOT` to
read `14_CRU/csv2_T8_CRU.bed`. Full mode expects the original `01_ref/`,
`03_exp/`, `04_acr/`, and `14_CRU/` layout. Missing inputs are reported before analysis.
Use `PYTHON_BIN` and `RSCRIPT_BIN` to select interpreters.

Results are written to `results/` and three PDFs to `figures/` under the output
directory. Repeated runs replace generated files. The demo passes 147 checks;
all 10 full-project statistics/plot tables have been verified against the original.

The demo starts from prepared effects, so it does not run the Python preparation
step. Files in `demo/expected/` are reference tables for comparison; the compact
runner performs internal validation but does not automatically compare them.

## Runtime

The bundled demo passes 147 checks in about a minute. Full mode re-reads the
complete raw matrices and is dominated by I/O.
