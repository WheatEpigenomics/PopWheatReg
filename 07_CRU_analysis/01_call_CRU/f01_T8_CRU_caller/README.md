# CRU analysis

Construction and scoring of cis-regulatory units (CRUs) from agLINK and
co-accessibility networks. A CRU is a contiguous interval of the genome that is
densely covered by ACR-gene links and behaves as one regulatory block.

## Principle

The unit is defined on agLINK density alone. Every agLINK contributes one count
regardless of PCC, distance or co-accessibility, so block boundaries reflect
where the network stops being dense, not how strong individual links are.
co-accessibility is used for the score and for display only; it never moves a
boundary.

Calling is one in-memory pass over the frozen rule set in `f02_modules/`:

1. agLINK endpoints and spans are binned at 5 kb and smoothed at 25/105/225 kb.
2. Genome stretches with no link at all over 100 kb or more are protected and
   are never filled.
3. Hard gaps and A-edge recovery define the core territories.
4. Each territory is split recursively at low-density valleys (depth 0.45).
5. Each part is subdivided again with the 300-kb daughter rule.
6. Four recall rounds re-admit low-density regions that have enough support
   (T3 density floor, T4/T5/T6 support and shape).
7. Isolated dense islands that passed the agLINK evidence gates are added back.

Boundaries are localized at the moment a candidate is accepted or a parent is
split, within a fixed +/-10 kb window around the original candidate centre.
There is no final BED sweep and no stage-end optimiser, so the emitted interval
is exactly the one that was accepted.

## Layout

```
01_call_CRU/
├── call_cru_main.sh          entry point, copy the blocks you need
└── f01_T8_CRU_caller/        all code and assets
    ├── 01_run_cru.py         calling entry, one chromosome per process
    ├── 02_call_chain.py      the in-memory rule chain
    ├── 03_fine_boundary.py   +/-10 kb boundary localization
    ├── 04_boundary_birth.py  transactional new-territory acceptance
    ├── 05_block_enrichment.py haplotype-block permutation test
    ├── 06_plot_windows.R     5-Mb atlas
    ├── 07_plot_block.R       block enrichment figure
    ├── 08_font_environment.py bundled-font fontconfig
    ├── 09_pdf_qc.py          PDF page QC
    ├── 10_score_cru.py       scoring
    ├── 11_plot_score.R       score distribution
    ├── f01_rules/            CRU_rules.json
    ├── f02_modules/          frozen rule modules
    ├── f03_scoring/          frozen f25 scorer and anchors
    ├── f04_fonts/            Liberation Sans
    └── demo/                 toy genome and runner
```

## Requirements

- Python 3.8 or newer with `numpy` and `scipy`
- `R` with `ggplot2` and a cairo-capable device, only for figures and block plots

The code lives in this directory; `call_cru_main.sh` one level up is the entry
point a user copies from. Run the commands from that parent directory. Nothing
is installed and nothing outside these two directories is touched.

## Inputs

| File | Format |
|---|---|
| agLINK | `gene_id  ACR_id  PCC  p` per line, no header |
| gene BED | BED6 `chr start end gene_id score strand` |
| ACR BED | BED4 `chr start end acr_id` |
| genome sizes | `chr <TAB> length` |
| haplotype blocks | BED3, optional, for the permutation test and the block track |
| co-accessibility | `id1  id2  correlation` with header, optional |

## Running

`../call_cru_main.sh` holds the commands as copy-paste blocks: calling, block
test, region figures, scoring, score figure, demo. The short version, run from
the parent directory:

```sh
python3 -B ./f01_T8_CRU_caller/01_run_cru.py \
  --aglink aglink.tsv --gene-bed genes.bed --acr-bed acrs.bed \
  --genome-sizes genome.sizes --config ./f01_T8_CRU_caller/f01_rules/CRU_rules.json \
  --chroms all --workers 21 --mode all --plot none --block none --outdir out_cru

python3 -B ./f01_T8_CRU_caller/10_score_cru.py \
  --cru-bed out_cru/f04_merged/CRU.bed --aglink aglink.tsv --acr-bed acrs.bed \
  --gene-bed genes.bed --genome-sizes genome.sizes --coacr coacr.tsv \
  --outdir out_score
```

One chromosome per worker process, 21 in parallel by default. `--resume` skips
any step whose inputs and outputs still match their recorded MD5, so an
interrupted run is restarted with the same command. Failed or damaged calling
files are moved to `f06_failed_attempts/` rather than deleted.

Modes: `--mode all` (call), `--mode block` (permutation test only),
`--mode plot` (figures only), plus `--block none` and `--plot none`.

### Figures never block calling

Plotting runs after a chromosome's calls are already validated and stored. If R
is missing, `ggplot2`/cairo is unavailable or any window fails, the run records
`f02_chromosomes/<chr>/03_PLOT.FAILED.json`, prints `PLOT FAILED <chr>`, and
carries on: merging, the block test and the summary all complete. The failures
are listed in `07_plot_failures.tsv` and redrawn later with

```sh
--mode plot --plot all --resume
```

without recomputing any call.

## Outputs

| Path | Content |
|---|---|
| `f04_merged/CRU.bed` | final CRUs, BED4, one row per interval |
| `f04_merged/01_chromosome_summary.tsv` | count, coverage, length quantiles, cut counts |
| `f02_chromosomes/<chr>/f02_calling/<chr>.CRU.bed` | per-chromosome calls |
| `f02_chromosomes/<chr>/f03_diagnostics/` | boundary evaluations, cuts, recall decisions, preserved hard gaps |
| `f02_chromosomes/<chr>/f05_figures/` | 5-Mb atlas, five tracks: agLINK, co-accessibility, density, CRU, blocks |
| `f06_block_enrichment/f01_statistics/` | observed vs permuted block-boundary matching |
| `f05_scripts/` | snapshot of the exact code that produced the run |
| `07_plot_failures.tsv` | only when a figure failed |

## CRU score

`10_score_cru.py` runs the frozen f25 scorer (`f03_scoring/02_f25_scorer.py`)
unmodified. The only override is the normalization callback: every density
metric is rescaled with the T7-A log1p q05/q95 anchors in
`f03_scoring/01_norm_anchors.tsv`, so scores stay comparable across cohorts and
the quantiles are not re-estimated on the data being scored.

For each network X in {agLINK, co-accessibility}:

```
X_per_mb          = internal_edges * 1e6 / length_bp
X_per_node_degree = 2 * internal_edges / nodes
                    nodes = acr_n + gene_n (agLINK), acr_n (co-accessibility)
z                 = clip( (log1p(x) - q05) / (q95 - q05), 0, 1 )
X_density         = mean( z(per_mb), z(per_node_degree) )
X_independence    = (internal + 1) / (touching + 2)      # NaN when touching = 0
score             = 100 * sum(w * component) / sum(w of available components)
```

Weights: agLINK density 0.40, co-accessibility density 0.25, agLINK
independence 0.20, co-accessibility independence 0.15. A point belongs to a CRU
when `start <= point < end`; merely spanning the interval does not count.
Independence is not quantile-normalized. CRUs whose evidence is thin are still
scored, with the available weights renormalized, and are flagged in
`evidence_status` rather than dropped.

Outputs: `01_CRU_scores.tsv` (all metrics), `02_CRU_scored.bed` (BED6, column 5
= `round(score*10)`, range 0-1000) and `04_CRU_scores_simple.tsv`
(`chr start end cru_id score`).

## Frozen rules

`f02_modules/` holds the rule modules exactly as they were used for the reported
calls. They are imported as a package and their algorithm code is untouched; the
only edit is turning the old hard-coded file paths into relative imports. See
`f02_modules/README.md` for what each module contributes.

## Demo

```sh
bash demo/run_demo.sh
```

A 20-Mb toy chromosome with 508 ACRs and 1,467 agLinks: calling and scoring take
a few seconds and write to `$TMPDIR/cru_demo`.
