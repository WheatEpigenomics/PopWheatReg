# Multivariate regression

Data and methods come from `publication_submission`. The package retains the
original 2-Mb/4-Mb analyses without SNP density, the 2-Mb analysis with SNP density,
the sensitivity models, and all six figures.

- `00_RUN.sh`: run the demo or the original project.
- `01_fit_models.R`: fit models, chromosome-cluster Webb bootstrap, and drop-one contributions.
- `02_plot_validate.R`: draw figures and independently validate the results.
- `demo/`: real input rows, parameter records, block IDs, source hashes and reference results.

Requirements: Bash, standard Unix tools, and R >= 3.6 with Cairo PDF support.
Only standard R packages are used, including `splines`. Arial is used for figures.
Set `RSCRIPT_BIN=/path/to/Rscript` to select an R installation.

Run from this directory:

```bash
# Default: refit the bundled real-data demo -> output_demo/
bash 00_RUN.sh

# Full analysis requires your prepared data directory
PROJECT_ROOT=/path/to/prepared_project bash 00_RUN.sh full
```

All modes accept an optional output directory, e.g. `bash 00_RUN.sh demo /path/to/output`.
`PROJECT_ROOT` is required for project/full mode and must contain the listed `data/` and `docs/` inputs. Missing inputs stop with a clear message before outputs are created. `full` is an alias for `project`.
Results go to `results/`, plotting tables to `plot_data/`, six PDFs to `figures/`,
and settings/checks to `docs/`. Repeated runs replace generated outputs.

The demo contains 252 original blocks: 12 per chromosome across all 21 chromosomes,
selected with Python `random.Random(20261113)` in sorted chromosome order. Both
contexts retain the same blocks and original row order, giving 504 rows in
`analysis_blocks.tsv.gz`; `SNP_density_2Mb.tsv.gz` contains the matching 252 SNP rows.
Input measurements are unchanged. `block_ids.txt` and `SOURCE_DATA.tsv` record the
subset; the SNP provenance table describes the original full-data archive recovery.
`analysis_settings.tsv` describes the SNP-excluded primary analysis. Settings for
the SNP-adjusted analysis are generated in `<OUTPUT_DIR>/docs/SNP_included_settings.tsv`.

Both modes start from these prepared block-level features and refit the models.
Raw genome-feature construction is outside this package. Each bootstrap test uses
the original 9,999 replicates and seeds. Demo standardization and results are
recomputed on the subset, so its values differ from the full 3,641-block analysis.
Reference tables in `demo/expected/` come from the original scripts on the same
subset; they are not used as fitting inputs or automatically compared by the runner.

The response is `log1p(block length)`. Original predictor transformations,
chromosome fixed effects, sequential models and drop-one definitions are retained.
Only fixed block-count checks and summary labels were adapted to the input size.
R-squared describes in-sample fit; drop-one contributions are not additive.

Demo and full-data runs each pass 381 checks with R 4.4.3. Demo tables match the
original scripts exactly; full-data tables agree with archived results within 1e-10.

## Runtime

The bundled demo runs all 381 internal checks in about a minute. A full project
run scales with the number of blocks and bootstrap replicates.
