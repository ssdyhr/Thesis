# Thesis analysis code and figures

Code for the thesis "Impact of age, sex and obesity on cancer biology" using TCGA data across 33 cancer types.

## Scripts

- `cibersortx_tcellextrect_plots.rmd` : immune cell deconvolution using CIBERSORTx and T cell fraction estimation using TcellExTRECT
- `driver_mutations_landscape.rmd` : oncoprint visualisation of driver mutation landscape across cancer types and clinical groups
- `driver_mutations_testing.rmd` : driver weight correlation, total mutation burden, and logistic regression analyses
- `pathway_analysis_final.Rmd` : visualisation of pathway enrichment results for RNA-seq and methylation data
- `spearman_correlation_meth.r` : spearman correlation analysis between methylation and clinical variables
- `spearman_correlation_rna.Rmd` : spearman correlation analysis between gene expression and clinical variables
- `missmeth.r` : overrepresentation analysis for methylation data using missMethyl, run per cancer type for BMI, age, and sex
- `GO_seq_enrich.r` : overrepresentation analysis for RNA-seq DEGs using GOseq, accounting for gene length bias

## Supplementary figures

- `CibersortX_supplementary` : all CIBERSORTx results across cell types and clinical variables
- `Oncoplots_supplementary` : oncoprints for all cancer types and clinical groups
