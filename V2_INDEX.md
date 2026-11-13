# PopWheatReg analysis code (v2)

Custom scripts and analysis pipelines for the population-scale wheat regulatory
genome study. The code follows the order of the paper: variant calling, chromatin
accessibility, regulatory links, molecular QTLs, phenotype association, and the
downstream CRU and TRACE network analyses.

Only code is shipped. Input data, intermediate files and results are not.

## Modules

| Module | What it does | Entry point | Demo |
|---|---|---|---|
| `01_Population_Genetics` | GATK variant calling, hard filtering, population structure, phylogenetic tree, LD decay and haplotype blocks | `01_call_SNP_and_filter.sh`, `02_pop_structure_and_tree.sh`, `03_LD_decay_and_genome_blocks.sh` | — |
| `02_Chromatin_Accessibility` | MACS2 peak calling, peak merging, 200-bp window counting, peak x sample matrix | `01_call_peak_and_norm.sh` | — |
| `03_ACR_Gene_Link` | agLink (ACR-gene Pearson + permutation z-test), ACR co-accessibility, SMR/HEIDI | `00_aglink_main.sh` | `demo/` |
| `04_Molecular_QTLs` | eQTL and caQTL association with TensorQTL, including the subgroup omnibus mode | `01_caQTL_eQTL_main.sh` | `TensorQTL_subgroup_omnibus/demo/` |
| `05_Phenotype_Association` | GWAS with EMMAX, caWAS and TWAS with PrediXcan | `01_GWAS_caWAS_TWAS_main.sh` | — |
| `06_Population_Simulations` | Multivariate regression, neutral Wright-Fisher simulations, forward simulations | `*/00_RUN.sh` | `*/demo/` |
| `07_CRU_analysis` | CRU construction and scoring, plus the CRU perturbation analysis | `01_call_CRU/call_cru_main.sh`, `02_CRU_perturbation/00_RUN.sh` | both |
| `08_TRACE_network` | TRACE V2-2 (C-agLINK) regulatory network construction and scoring | `call_trace_main.sh` | `demo/` |

Each module has its own `README.md` with the input formats, the commands and the
runtime; start there.

## How to run

Every entry point is a shell script organised as copy-paste blocks. Read the
block, set the paths at the top, then run the blocks in order. Nothing is
installed and no state is shared between modules; each one reads the outputs of
the previous module from disk.

Modules 01, 02 and 05 wrap external tools (GATK, bcftools, PLINK, MACS2, bedtools,
EMMAX, PrediXcan) and ship no demo, because a meaningful run needs those tools and
full-size inputs. Modules 03, 04, 06, 07 and 08 ship a toy dataset that exercises
every step in seconds.

## Software

| Module | Requirements |
|---|---|
| 01 | GATK 4.1, bcftools, PLINK 1.9, ADMIXTURE, VCF2Dis, FastME, Python 3 with numpy/matplotlib |
| 02 | MACS2, bedtools |
| 03 | Python 3 with numpy/scipy/pandas, bedtools, SMR 1.4.1, PLINK bed (GPU for co-accessibility) |
| 04 | TensorQTL 1.0.10, PyTorch, pandas, numpy; scipy for the subgroup omnibus |
| 05 | EMMAX, PrediXcan (Python 2) |
| 06 | Python 3 with numpy (matplotlib for the forward simulation), R with ggplot2 and Cairo |
| 07 | Python 3 with numpy/scipy; R for figures only |
| 08 | Python 3 with duckdb and pandas |

## Runtime

Demo runtimes for the custom code in this package. Public tools (GATK, MACS2,
PLINK, EMMAX, PrediXcan) are excluded, because their cost depends on the local
installation. Full runs scale with the size of the input.

| Analysis | Demo runtime |
|---|---|
| agLink | under a minute |
| ACR co-accessibility | under a minute |
| eQTL / caQTL association | under a minute on CPU |
| Population simulations | about 5 min |
| CRU construction and scoring | under a minute |
| TRACE network | under a minute |

## Contact

zhang_jiazi@163.com, luzefu@caas.cn
