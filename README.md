# integrative-genomics-gwas-eqtl-rnaseq
Integrative genomics case study in R linking GWAS, liver eQTL, and bulk RNA-seq to explore locus-to-gene relationships and disease-associated transcriptional change.

# Integrative genomics case study: GWAS, liver eQTL, and bulk RNA-seq

## Overview

This repository contains a reproducible R-based case study integrating:

- GWAS summary statistics
- liver eQTL data
- bulk RNA-seq data

The goal is to explore whether a disease-associated locus can be linked to a candidate gene through regulatory evidence and tissue-level expression changes.

## Main workflow

1. GWAS quality control and visualization
2. Identification of significant loci
3. Regional comparison of GWAS and liver eQTL signals
4. Colocalisation analysis
5. Bulk RNA-seq differential expression
6. Candidate gene expression visualization

## Repository structure

```text
.
├── README.md
├── .gitignore
├── LICENSE
├── data/
│   ├── raw/
│   └── processed/
├── results/
│   ├── figures/
│   └── tables/
├── scripts/
└── docs/
