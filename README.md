# integrative-genomics-gwas-eqtl-rnaseq
Reproducible R workflow integrating GWAS, liver eQTL, and RNA-seq for candidate gene prioritization and disease-context interpretation.

# Integrative genomics case study: GWAS, liver eQTL, and RNA-seq

## Overview

This repository presents a reproducible R-based integrative genomics workflow combining:

- Genome-wide association study (GWAS) summary statistics
- Liver eQTL data (GTEx)
- Bulk RNA-seq data (TCGA-style workflow)

The objective is to determine whether a disease-associated locus can be linked to a candidate gene via regulatory evidence and whether that gene shows transcriptional changes in the relevant tissue.

This project is presented as a portfolio case study demonstrating applied bioinformatics and data integration skills.

---

## Key question

Can a GWAS locus be connected to a biologically plausible target gene using eQTL evidence, and does that gene show differential expression in disease-relevant tissue?

---

## Workflow

### 1. GWAS quality control and visualization
- Data cleaning and filtering
- Manhattan plot
- QQ plot

### 2. Signal prioritization
- Genome-wide significant variants
- Distance-based selection of approximately independent loci

### 3. Regional GWAS–eQTL integration
- Focus on a lead locus
- Overlay of GWAS and liver eQTL signals
- Candidate gene: **PNPLA3**

### 4. Colocalisation analysis
- Bayesian colocalisation using `coloc`
- Estimation of posterior probabilities for shared causal variants

### 5. RNA-seq context (optional module)
- Differential expression analysis (DESeq2-style workflow)
- Candidate gene expression visualization

---

## Repository structure

```text
.
├── README.md
├── .gitignore
├── LICENSE
├── data/
│   ├── raw/          # input files (not tracked)
│   └── processed/    # intermediate files
├── results/
│   ├── figures/      # plots and visual outputs
│   └── tables/       # processed tables
├── scripts/          # analysis scripts
└── docs/             # project documentation
