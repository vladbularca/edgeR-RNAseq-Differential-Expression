# Flor-Yeast-Differentially-Expressed-Genes

## Introduction

This project aims to identify differentially expressed genes (DEGs) between  developmental stages of Saccharomyces cerevisiae (yeast) velum development. This analysis identifies genes expression from RNA-sequencing (RNA-seq) by quantifying gene expression changes across different stages. 

## Objectives

- Identify DEGs using edgeR 
- Implement quality control and preprocessing of raw sequence data
- Perform statistical analysis to determine significant gene expression changes
- Visualize and interpret results 

## Methods

Tools Used:

- STAR: Read alignment to reference genome 
- SAMtools: Processing and manipulation of alignment files
- Subread (featureCounts): Gene expression quantification 
- edgeR (RStudio): Differential expression analysis and statistical modeling 

### Workflow Overview:

- Data Preparation: Retrieve RNA-Seq data and preprocess FASTQ files
- Sequence Alignment: Align reads to the reference genome with STAR, followed by BAM file processing using SAMtools
- Gene Quantification: Use featureCounts to generate a count matrix for differential expression analysis

### Differential Expression Analysis:
- MDS for intial clustering
- Filter low-expression genes
- Normalize count data
- Estimate dispersion and fit statistical models in edgeR

## Results

- A comprehensive list of genes that were found to be differentially expressed across all pairwise comparisons was identified. Top 10 signficant DEGs for each pairwise comparison were also determined, highlighting stage specific expression patterns. 
- Top DEGs were identified with log fold change (logFC) and statistical significance (FDR)
- Key findings highlight core genes involved in velum development at different stages

## Repository Structure 

```text ├── gene_counts.txt #  Output of featureCounts 
        ├── STAR_alignment.sh # Bash Script file
        ├── *genes.txt/ # Output files with DEGs
        ├── README.md # Project documentation
        ├── Assignment2.pdf # Full assignment write-up
        ├── edgeR_analysis.R # R Script for edgeR analysis
 ```
