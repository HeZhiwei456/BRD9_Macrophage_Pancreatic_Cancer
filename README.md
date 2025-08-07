# BRD9-Mediated Epigenetic Regulation of M2 Macrophage Polarization

This repository contains the complete bioinformatics pipeline and supporting analysis files used in the study:

**He Zhiwei, et al. (2025)**  
*"BRD9-Mediated Epigenetic Regulation of M2 Macrophage Polarization Drives Immune Evasion and Cancer Cachexia in Pancreatic Cancer"*

---

## 🔬 Project Overview

This project integrates ATAC-seq, Cut&Tag, RNA-seq, and TCGA-based immunogenomic data to investigate the epigenetic role of BRD9 in regulating M2 macrophage polarization and immune evasion in pancreatic cancer. It includes raw data processing scripts, differential expression/peak analysis, motif enrichment, survival analysis, and figure-generating outputs.

---

## 📁 Repository Structure

scripts/ # R and shell scripts for ATAC-seq, RNA-seq, Cut&Tag analysis
results/ # Processed output files (e.g., GSEA, KM, ROC, TCGA correlations)
requirements_R.txt # Required R packages
LICENSE # MIT License for open-source use
README.md # Project documentation
.gitignore # Files excluded from version control


---

## 🔧 Software Requirements

- **R version ≥ 4.1**
- **Linux/Unix shell environment**

### Required tools:
- STAR
- Bowtie2
- Samtools
- MACS2
- Picard
- HOMER

### Required R packages:
- DESeq2
- edgeR
- ggplot2
- pheatmap
- clusterProfiler
- org.Mm.eg.db
- TxDb.Mmusculus.UCSC.mm10.knownGene
- ChIPseeker
- DiffBind
- biomaRt
- EnhancedVolcano
- ggrepel

> All packages are listed in `requirements_R.txt`.

📊 Results
The results/ directory contains processed output files used to generate key figures in the manuscript, including:

TCGA gene expression correlations

Cox regression and Kaplan-Meier survival analysis

GSEA enrichment results

ROC curve data

Immune signature and infiltration correlation

Single-cell summary tables

📄 Citation
If you use this repository or any part of the workflow, please cite the original publication:

He Zhiwei, et al. (2025)
"BRD9-Mediated Epigenetic Regulation of M2 Macrophage Polarization Drives Immune Evasion and Cancer Cachexia in Pancreatic Cancer"
