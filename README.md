# Differential Expression Analysis and Biological Insights

This repository contains the code and resources for performing differential expression analysis on SARS-CoV-2 and lung cancer cell line data. The analysis focuses on identifying differentially expressed genes (DEGs) in blood and leukocyte samples, as well as differentially expressed microRNAs (DEMs) in lung cancer cell lines.

## Project Overview

### Datasets
- **GSE157103**: Expression profiles from blood samples of SARS-CoV-2 positive and negative patients.
- **GSE148729**: Small RNA expression profiles from SARS-CoV-2 infected and non-infected Calu-3 human epithelial cell lines.

### Objectives
1. Identify DEGs in plasma and leukocyte samples of COVID-19 patients.
2. Identify DEMs in lung cancer cell lines.
3. Perform pathway enrichment analysis to reveal biological insights related to viral infection.

## Results Summary

### Principal Component Analysis (PCA)
- PCA was conducted on GSE157103 to ensure factors like gender, age, and Charlson severity score did not contribute to variation in gene expression.

### Differential Expression Analysis
- **GSE157103**: Identified 464 DEGs (389 upregulated, 75 downregulated) with an adjusted p-value < 0.05.
- **GSE148729**: Identified 37 DEMs (15 upregulated, 22 downregulated).

### Pathway Enrichment
- Reactome pathway analysis revealed 97 enriched pathways with a False Discovery Rate (FDR) cutoff of 0.01. Top pathways include:
  - Cell cycle regulation
  - DNA replication
  - RHO GTPase signaling

## Methods

### Data Retrieval
- The datasets were retrieved from the Gene Expression Omnibus (GEO) repository.
  - **GSE157103**: 126 blood samples (100 COVID-19 positive, 26 negative).
  - **GSE148729**: RNA profiles from Calu-3 cell lines (mock and infected conditions).

### Differential Expression Analysis
- **Tool**: DESeq2 Bioconductor package was used for analysis of DEGs and DEMs.
- **Thresholds**: Adjusted p-value < 0.05 and log2FoldChange > 0.5 were used to identify significant genes.

### Pathway Enrichment Analysis
- **Tool**: Reactome pathway analysis tool was used to identify enriched pathways.
- **Cutoff**: FDR < 0.01.

## Repository Structure
- `data/`: Contains the retrieved datasets from GEO.
- `notebooks/`: Jupyter notebooks with data analysis and visualization.
- `scripts/`: Python scripts used for preprocessing, analysis, and plotting.
- `results/`: Results from differential expression and pathway enrichment analyses.

## Dependencies
- R (version ≥ 4.0.0)
- DESeq2
- ReactomePA
- Python (version ≥ 3.8)
- Pandas
- Matplotlib


## References
- [GSE157103 Dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE157103)
- [GSE148729 Dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE148729)
- [Reactome Pathway Database](https://reactome.org/)
- [SARS-CoV-2 potential drugs, drug targets, and biomarkers: a viral-host interaction network-based analysis](https://www.nature.com/articles/s41598-022-15898-w#Sec9)

