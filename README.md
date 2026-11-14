# PopWheatReg
Population-scale integrative analysis of chromatin accessibility, gene expression and regulatory variation underlying spike development and agronomic traits in polyploid wheat.

Noncoding regulatory variation is an important contributor to phenotypic diversity in crops, yet its regulatory mechanisms remain difficult to resolve in large and polyploid genomes such as bread wheat (*Triticum aestivum*). In particular, how genetic variants influence chromatin accessibility, gene expression and agronomic traits through local and long-range regulatory interactions remains poorly understood at the population scale.

In this project, we generated a population-scale chromatin accessibility atlas of wheat spikes by profiling over 200 globally diverse wheat accessions using ATAC-seq, together with RNA-seq and additional whole-genome resequencing data. By integrating chromatin accessibility, gene expression, genetic variation and population structure, we systematically characterized the regulatory landscape underlying spike development and agronomic trait variation in polyploid wheat.
<p align="center">
  <img width="750" alt="PopWheatReg overview" src="https://github.com/user-attachments/assets/81542193-3e3d-4783-a934-f5c39dde999a" />
</p>





This repository contains custom scripts and analysis pipelines used to:

- Identify accessible chromatin regions (ACRs) and quantify their population-level variation.
- Map chromatin accessibility QTLs (caQTLs) and expression QTLs (eQTLs), accounting for population structure and subgroup-specific associations.
- Identify ACR–gene regulatory links (agLinks) and chromatin co-accessibility relationships (coACRs).
- Construct *cis*-regulatory units (CRUs) and characterize their genomic organization, boundary properties and robustness.
- Apply TRACE to reconstruct putative regulatory pathways linking genetic variants to target genes through direct SNP–ACR–gene associations or transcription factor-mediated downstream ACR–gene regulation.
- Integrate regulatory mapping results with GWAS, TWAS and caWAS analyses to prioritize candidate causal regulatory variants and trait-associated genes.


Together, these analyses provide a framework for investigating how noncoding regulatory variation shapes gene expression and complex traits in polyploid crop genomes. The scripts in this repository support reproducible population-scale regulatory analyses and may be adapted for similar studies in other plant species.

**For detailed code organization, usage instructions and software requirements, see [V2 Code Index](V2_INDEX.md).**

If you have any questions, requests, or comments, we kindly invite you to contact us at zhang_jiazi@163.com, luzefu@caas.cn.
