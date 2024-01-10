# Merge patient ancestry and histology data

__Module Authors: Ryan Corbett__ ([@rjcorb](https://github.com/rjcorb))

This analysis merges patient ancestry prediction data with OpenPedCan histologies data for downstream analyses


## Usage
`Rscript -e "rmarkdown::render('01-add_histologies.Rmd')"`
`Rscript 02-summary_stats.R`

## Folder content

- `01-add_histologies.Rmd`; merges ancestry data from Somalier with OpenPedCan histologies

- `02-summary_stats.R`; plots summary statistics including PCA plots, concordance between ancestry and reported race and ethnicity, plot group by ancestry, etc. 

- `input/` files: 
  - `all_pnoc_normal_wgs.tsv`; biospecimen and participant IDs of PNOC normal samples
  - `Broad_to_BS.txt`; key for X01 patient biospecimen IDs
  - `DEI_CBTN-PNOC_rerun.somalier-ancestry.tsv`; Somalier ancestry output for all WGS normal samples
  - `histologies.tsv`; OpenPedCan v12 histologies file
  - `Meyer_Tableau_Pull_11.8.2022.xlsm`; contains data on EU patients that are to be filtered out for this project
  - `plot-mapping.tsv` table of `plot_group` assignments by `broad histology` and `cancer_group`

- `results/` files: 

` `plots/` files: 
