# Flor-Yeast-Differentially-Expressed-Genes

## Introduction

This project aims to identify differentially expressed genes (DEGs) between  developmental stages of Saccharomyces cerevisiae (yeast) velum development. This analysis identifies genes expression from RNA-sequencing (RNA-seq) by quantifying gene expression changes across different stages. 

## Objectives

- Identify DEGs using edgeR 
- Implement quality control and preprocessing of raw sequence data
- Perform statistical analysis to determinesignificant gene expression changes
- Visualize and interpret results 

## Methods

Tools Used:

- STAR: Read alignment to reference genome 
- SAMtools: Processing and manipulation of alignment files
- Subread (featureCounts): Gene expression quantification 
- edgeR (RStudio): Differential expression analysis and statistical modeling 

## Workflow Overview:

- Data Preparation: Retrieve RNA-Seq data and preprocess FASTQ files
- Sequence Alignment: Align reads to the reference genome with STAR, followed by BAM file processing using SAMtools
- Gene Quantification: Use featureCounts to generate a count matrix for differential expression analysis

### Differential Expression Analysis:
- MDS for intial clustering
- Filter low-expression genes
- Normalize count data
- Estimate dispersion and fit statistical models in edgeR

### Visualization:
Volcano Plots: Show significant DEGs based on log fold change and p-values.
Venn Diagrams: Display overlap of DEGs between developmental stages.
Smear Plots: Illustrate distribution of gene expression changes.

Results:

A comprehensive list of genes that were found to be differentially expressed across all pairwise comparisons
Top DEGs were identified with log fold change (logFC) and statistical significance
Key findings highlight core genes involved in biofilm development and potential regulators of yeast velum formation

Discussion
The results indicate dynamic gene expression changes during velum development, with core regulatory genes identified as potential candidates for further experimental validation. The study highlights:

The reliability of RNA-Seq for transcriptomic analysis
The importance of stringent statistical analysis for DEG identification
Future directions, including functional validation and pathway enrichment analysis

Repository Structure
├── featureGenecount.txt/     # Bash Output Feature Count data
├── DGE_analysis.sh           # Bash Script file
├── Visualizations.R          # R Script file
├── Results/                  # Output files including DEG lists
├── Figures/                  # Plots (volcano, smear, venn diagrams)
├── README.md                 # Project documentation
├── Assignment2.pdf           # Full assignment write-up

