# Project Summary

## Overview

This project demonstrates an integrative genomics workflow combining GWAS, liver eQTL, and RNA-seq data to explore how genetic variation may influence gene expression in disease.

The analysis is structured as a reproducible R-based pipeline and presented as a portfolio case study.

---

## Objective

To investigate whether a disease-associated GWAS locus can be linked to a candidate gene using regulatory evidence (eQTL), and whether that gene shows differential expression in relevant tissue.

---

## Approach

The workflow consists of five main steps:

1. **GWAS processing and visualization**
   - Data cleaning and filtering
   - Manhattan and QQ plots

2. **Signal prioritization**
   - Identification of genome-wide significant variants
   - Distance-based selection of approximately independent loci

3. **Regional GWAS–eQTL integration**
   - Extraction of variants around a lead locus
   - Comparison of GWAS and liver eQTL signals
   - Focus on candidate gene: **PNPLA3**

4. **Colocalisation analysis**
   - Bayesian colocalisation using the `coloc` package
   - Estimation of posterior probabilities for shared genetic signals

5. **RNA-seq analysis (optional module)**
   - Differential expression analysis
   - Visualization of candidate gene expression

---

## Key Concept

This project follows a **locus-to-gene prioritization framework**:

- GWAS identifies disease-associated variants  
- eQTL links variants to gene expression  
- RNA-seq provides biological context  

---

## Outputs

The analysis generates:

- Manhattan and QQ plots for GWAS data
- Independent locus summary table
- Regional GWAS–eQTL comparison plot
- Colocalisation posterior probability results
- RNA-seq differential expression outputs (optional)

---

## Tools and Methods

- R programming
- data.table, dplyr (data processing)
- qqman (GWAS visualization)
- ggplot2 (visualization)
- coloc (colocalisation analysis)
- DESeq2 (RNA-seq analysis)

---

## Limitations

- Distance-based locus selection does not replace LD-based clumping
- Colocalisation results depend on data harmonization and assumptions
- eQTL and GWAS signals may not reflect identical causal variants
- RNA-seq results reflect bulk tissue and do not establish causality

---

## Conclusion

This case study demonstrates how integrating multiple genomic data types can help move from genetic association signals toward biological interpretation and candidate gene prioritization.

---

## Purpose

This project is intended as a portfolio piece to demonstrate practical skills in bioinformatics, statistical genomics, and reproducible analysis workflows.