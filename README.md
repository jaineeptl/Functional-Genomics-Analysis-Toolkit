# Functional-Genomics-Analysis-Toolkit
Reproducible R workflow for differential expression and Gene Ontology (GO) enrichment from Seurat objects. Automates gene ranking, BP enrichment (clusterProfiler), optional GO subtree filtering, and generates publication-ready barplots and chord diagrams!
This repository is designed to turn statistical differential expression into a clean and scalable interepretation. Instead of running ad hoc enrichment scripts inside notebooks, this workflow standardizes ranking, enrichment, subtree filtering, and figure generation into a configurable pipeline.

# What This Workflow Does
This pipeline takes a Seurat object and:
- Runs differential expression between two groups
- Ranks genes by combined effect size × statistical significance
- Performs GO Biological Process enrichment using clusterProfiler
- Optionally filters results to a GO subtree 
- Generates publication-ready:
    - GO barplots with term-level mean logFC coloring
    - Chord diagrams linking GO terms to contributing genes
    - Gene-level logFC visualizations for specific GO terms
- Exports structured result tables (DE + GO) for downstream use
- Everything is controlled through a single config.yml file

# Why This Exists
In many transcriptomic projects, GO enrichment becomes an inconsistent, copy-pasted step at the end of analysis. This workflow formalizes that step:
- Removes hard-coded dataset logic
- Avoids manual threshold tweaking in notebooks
- Standardizes ranking fallback when few genes pass p-value thresholds
- Produces consistent figures across projects
- Separates analysis logic from dataset-specific biology
- It is built to be reusable across datasets and conditions.

# Repository Structure
.
├── config.yml
├── 00_setup.R
├── 01_load_seurat.R
├── 02_subset_filters.R
├── 03_differential_expression.R
├── 04_rank_genes.R
├── 05_go_enrichment.R
├── 06_visualize_go_barplots.R
├── 07_visualize_go_chordplots.R
├── 08_query_specific_go_terms.R
├── 99_run_all.R
└── R/
    ├── utils_io.R
    ├── utils_go.R
    └── utils_plot.R

# Requirements
R ≥ 4.2 recommended
Core packages:
Seurat
clusterProfiler
org.Hs.eg.db (or org.Mm.eg.db if adapting for mouse)
GO.db
dplyr
tidyr
ggplot2
circlize
RColorBrewer
yaml

You can install missing packages via:
install.packages(c("dplyr","tidyr","ggplot2","yaml","RColorBrewer"))
BiocManager::install(c("Seurat","clusterProfiler","org.Hs.eg.db","GO.db"))

# Input Requirements 
A Seurat object (.rds) containing:
- A metadata column to compare (e.g., Cell_Line, orig.ident)
- Two groups defined in that metadata column
- A normalized assay (RNA or SCT)
No raw FASTQ or count matrices required.

# How to Run
1. Edit config.yml:
    Set the Seurat RDS path
    Set ident_col, group_a, and group_b
    Adjust optional filters or GO subtree root if you'd like
2. Run
   source("99_run_all.R")
   
Outputs will be written to:
results/ → DE + GO tables
figures/ → barplots, chord diagrams, gene-level plots

# Example Use Case
Compare two experimental conditions:
Run DE
Extract top-ranked genes
Perform GO enrichment on upregulated genes
Restrict to a specific biological process subtree
Visualize enriched terms with effect size context

# Extensibility
This workflow can be extended to:
Mouse datasets (swap OrgDb)
KEGG enrichment
GSEA
Multi-condition comparisons
Integration into Snakemake or targets pipelines
Containerized execution (Docker)


