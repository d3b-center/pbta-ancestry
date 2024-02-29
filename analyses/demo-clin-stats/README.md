# Generate demographic and clinical summary statistics in PBTA ancestry cohort

__Module Authors: Ryan Corbett__ ([@rjcorb](https://github.com/rjcorb))

This analysis module calculates demographic and clinical summary statistics in PBTA ancestry cohort by genetic ancestry superpopulation and cancer group

## Usage
`bash run_module.sh`

## Folder content

- `01-demo-clin-stats.Rmd`; Generates summary statistics and p-value tables by genetic ancestry superpopulation and cancer group

- `results/` files: 
  - `demo-clin-stats-all.tsv`; table of summary stats by genetic ancestry superpopulation
  - `demo-clin-pvals-all.tsv`; corresponding summary stat p-values 
  - `demo-clin-pvalues-by-histology.xlsx`; summary stat tables by genetic ancestry superpopulation within cancer groups
  - `demo-clin-stats-by-histology.xlsx`; corresponding summary stat p-values by genetic ancestry superpopulation within cancer groups

- `util/` files: 
  - `summary_functions.R`; function script to generate frequency tables 
