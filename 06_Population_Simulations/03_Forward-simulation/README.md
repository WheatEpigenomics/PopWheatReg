# Per-block PCC top 1% forward simulations

The default command runs a self-contained demo, including preparation,
simulation, checks, and all six PDF figures:

```bash
bash 00_RUN.sh                         # demo -> output_demo/
bash 00_RUN.sh demo /path/to/output    # optional new output directory
RAW_DIR=/path/to/raw_data bash 00_RUN.sh full
```

The same original algorithm is used in both modes: after promoter filtering,
representative-SNP mapping and SNP-pair deduplication, keep
`max(1, ceil(n * 0.01))` links ranked by raw PCC **within each ACR block**.
Compare background only, distance-matched random selection, and PCC selection.

| Setting | Demo | Full/project |
| --- | --- | --- |
| Chromosomes | 1A input excerpt | 1A–7D |
| Population | 203 | 10,000 |
| Replicates | 3 | 30 |
| Generations | G0, G1, G2 | G0–G300 |
| LD sample size | 203 | 203 |

Demo inputs are 6 ACRs, 3 real links, 6 genes, 6 blocks, and 12 SNPs with
203 original samples. Demo settings illustrate the workflow, not full-study
inference. Selfing, recombination, selection strength, matching rules and seed
are unchanged. Resolved settings are saved to `<OUTPUT>/config.json`.

Full mode requires `RAW_DIR` with `05_acr_with_name_fil110.bed`,
`06_aglink_sig_only.tsv`, `haplotype_block.bed`, `Ta_CS_v2.1.gff3`, and
`vcf_file/` containing VCFs for every configured chromosome. No parent-project
path is assumed. Missing inputs produce a message before simulation starts.

- `00_RUN.sh`: demo/full entry point; `project` aliases `full`.
- `01_forward_simulation.py`: prepare, simulate, summarize, and audit.
- `02_plot_results.py`: independent plotting and combined tables.
- `config.json`: full-analysis settings; demo overrides only the table above and paths.
- `demo/`: real BED/TSV/GFF3/VCF inputs. BED is 0-based half-open; GFF3/VCF are 1-based.

Requirements: Bash, Python >=3.8, NumPy >=1.23, Matplotlib >=3.6.
Use `python3 -m pip install 'numpy>=1.23' 'matplotlib>=3.6'` in your environment,
or select an installed environment with `PYTHON_BIN` for the runner.
Outputs contain prepared pairs, simulation arrays, TSV summaries, audit reports,
and PDF figures. Use a new output directory for another run.

Redraw completed results without repeating simulations:

```bash
python3 02_plot_results.py --config output_demo/config.json \
  --input output_demo --output replot
```

Use the saved config from the corresponding full output for full-result plots.
Add `--chrom 1A` to draw only that chromosome's four PDFs. Plotting checks saved
metadata and arrays before drawing; `combined/` contains equal-weight chromosome means.

## Runtime

The bundled demo runs three replicates on chromosome 1A in a couple of minutes.
A full run covers 21 chromosomes x 30 replicates and is two to three orders of
magnitude longer.
