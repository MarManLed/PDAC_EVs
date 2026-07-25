# PDAC secretome–NK cell integrative analysis

This repository documents the bioinformatic workflow used to integrate proteomic secretome data from patient-derived pancreatic ductal adenocarcinoma cultures with single-cell RNA-sequencing profiles of NK cells from tumor and adjacent tissues.

## Script

`PDAC_secretome_NK_illustrative_analysis.R` is an intentionally non-executable, condensed representation of the exploratory analyses originally implemented in the project R Markdown notebooks.

The script illustrates:

1. secretome preprocessing and missing-value handling;
2. protein-level univariate analysis, effect-size calculation, PCA and GSEA;
3. exploratory BDA, random forest, SVM, XGBoost and SHAP analyses;
4. Seurat processing of GSE212966 PDAC and adjacent-tissue samples;
5. SingleR annotation and NK-cell selection;
6. filtering of expressed NK-cell plasma-membrane genes;
7. ExoCarta and GO-based selection of theoretical EV-associated secretome proteins;
8. STRING interaction queries and intergroup edge filtering;
9. Cytoscape visualization;
10. CytoHubba-based local-hub, bridge and global-influencer prioritization.

## Important note

The code is provided for methodological transparency and does not constitute a turnkey reproducible pipeline. Local paths, source data organization, intermediate `.rda` objects, Cytoscape exports and custom functions must be adapted before execution.

The original notebooks remain the authoritative exploratory code record.
