# Mito-Nuclear Association Analysis

## Overview

This repository contains an R workflow for identifying nuclear loci associated with mitochondrial haplotypes (mito-nuclear interactions) using genome-wide SNP data.

The pipeline:

1. Converts VCF genotype data to GDS format.
2. Identifies SNPs located within a predefined set of mitochondrial-function genes.
3. Performs LD pruning.
4. Tests each SNP for association with mitochondrial haplotype.
5. Generates an empirical null distribution through repeated sampling of background SNPs.
6. Evaluates enrichment of significant mito-gene SNPs relative to the genomic background.
7. Identifies candidate genes associated with significant SNPs.

## Biological Motivation

Mitochondrial and nuclear genomes jointly encode components required for mitochondrial function. Because these genomes coevolve, mitochondrial haplotypes may show non-random associations with nuclear alleles involved in oxidative phosphorylation, metabolism, and other mitochondrial processes.

This workflow tests whether SNPs in mitochondrial-function genes exhibit stronger associations with mitochondrial haplotypes than expected by chance.

## Required Input Files

- `mallard.CDS_biallelic_snps.vcf.gz` – Biallelic nuclear SNP dataset
- `mito_haplotypes.csv` – Sample IDs and mitochondrial haplotype assignments
- `admixture_Qmatrix.csv` – Admixture ancestry coefficients
- `sample_metadata.csv` – Additional sample metadata
- `mallard_omics_numt_CDS.gff` – GFF file containing mitochondrial-function gene annotations

## Dependencies

### CRAN

- tidyverse
- lme4

### Bioconductor

- SNPRelate
- SeqArray
- GenomicRanges
- rtracklayer

## Workflow

### 1. Convert VCF to GDS

The VCF file is converted into GDS format using SNPRelate for efficient genotype access.

### 2. Define Candidate Mitochondrial Genes

The GFF annotation file is parsed to extract CDS coordinates for mitochondrial-function genes.

### 3. Identify Candidate and Background SNPs

- **Mito SNPs**: SNPs located within mitochondrial-function genes
- **Background SNPs**: all other SNPs

### 4. LD Pruning

LD pruning is performed separately on mito-gene SNPs and background SNPs using a threshold of 0.2.

### 5. Association Testing

Genotypes are tested for association with mitochondrial haplotype after merging haplotype, ancestry, and metadata tables.

### 6. Multiple Testing Correction

P-values are corrected using the Benjamini–Hochberg FDR procedure.

### 7. Empirical Null Distribution

The script generates 500 random background SNP subsets to estimate the expected number of significant associations under the null hypothesis.

### 8. Fisher Enrichment Test

A Fisher's Exact Test evaluates whether mito-gene SNPs are enriched for significant associations compared with background SNPs.

### 9. Candidate Gene Identification

Significant SNPs are mapped back to annotated genes.

Outputs:

- `candidate_genes.csv`
- `top_candidate_genes.csv`

## Running the Analysis

```r
source("mito_nuclear_analysis_CLEAN.R")
```

## Outputs

- `nuclear.gds`
- `candidate_genes.csv`
- `top_candidate_genes.csv`
- Histogram of empirical null distribution
- Console summary statistics and enrichment results
