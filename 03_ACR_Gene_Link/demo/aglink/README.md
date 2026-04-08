# 03_ACR_Gene_Link

This directory contains key custom scripts used for ACR–gene association analysis in the study.

## Overview

The scripts in this directory were used to quantify associations between gene expression and chromatin accessibility across the wheat population and to identify significant ACR–gene links (agLinks). The core logic includes:

1. calculating Pearson correlation coefficients (PCCs) between gene expression and candidate ACR accessibility profiles;
2. organizing correlation results for downstream statistical testing;
3. performing significance assessment for agLink identification;
4. supporting downstream co-accessibility-related analyses.

These scripts represent core custom analysis modules used in the study. In the full project workflow, some intermediate files were generated dynamically during runtime and removed after use; therefore, they are not included in this repository.

## Repository contents

- `01_calcu_acr_gene_pcc.py`  
  Calculates PCCs between one gene expression profile and multiple candidate ACR accessibility profiles.

- `02_fill_column_for_ztest.sh`  
  Prepares intermediate tabular output for downstream statistical testing.

- `03_aglink_ztest.py`  
  Performs downstream statistical assessment for agLink identification.

- `demo_gene_real_matrix.txt`  
  A minimal toy input file for illustrating the core correlation analysis.

- `demo_acr_good.list`  
  A minimal example ACR identifier list corresponding to the ACR rows in the toy matrix.

## Input format for the toy example

The file `demo_gene_real_matrix.txt` is a minimal example input for the core correlation analysis.

- The **first row** represents the expression profile of one gene across samples.
- Each **subsequent row** represents the accessibility profile of one candidate ACR across the same samples.
- Columns correspond to the same set and order of samples across all rows.

Example structure:

- row 1: gene expression
- row 2+: candidate ACR accessibility

The file `demo_acr_good.list` contains the ACR identifiers corresponding to rows 2 onward in the matrix.

## System requirements

The scripts were developed and tested in a standard Linux environment.

Recommended environment:
- Linux
- Python 3
- Bash shell

No non-standard hardware is required for the toy example.

## Dependencies

The core scripts use standard Python libraries and Unix command-line utilities.

Typical dependencies include:
- Python 3
- `numpy`
- `scipy`
- standard Unix tools such as `awk`, `paste`, `grep`, and `shuf`

Please ensure these are available in your environment before running the scripts.

## Installation

No formal installation is required beyond ensuring that the required Python environment and standard command-line tools are available.

Typical setup time on a normal desktop computer is less than 10 minutes, assuming Python and the required packages are already installed.

## Demo

A minimal toy input is provided to illustrate the core correlation analysis.

### Example command

```bash
python 01_calcu_acr_gene_pcc.py demo_gene_real_matrix.txt > demo_pcc_output.txt