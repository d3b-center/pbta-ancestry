# Merge patient ancestry and histology data

__Module Authors: Ryan Corbett__ ([@rjcorb](https://github.com/rjcorb))

This analysis merges patient ancestry prediction data with OpenPedCan histologies data for downstream analyses


## Usage
`bash run_module.sh`

## Folder content

- `01-add_histologies.Rmd`; merges ancestry data from Somalier with OpenPedCan histologies

- `02-summary_stats.R`; plots summary statistics including PCA plots, concordance between ancestry and reported race and ethnicity, plot group by ancestry, etc. 

- `03-aluvial_plot.R`; generates alluvial plot showing relationship between genetic ancestry superpopulation, reported race, and reported ethnicity across cohort

- `input/` files: 
  - `plot-mapping.tsv` table of `plot_group` assignments by `broad histology` and `cancer_group`

- `results/` files: 
  - `merged_ancestry_histology_data.tsv`; merged ancestry and matched tumor histology data

- `plots/` files: 
  - `ancestry-pcs.pdf`; plot of somalier PCs 1-4
  - `ancestry-race-ethnicity-alluvial.pdf`; alluvial plot of cohort predicted ancestry, reported race, and reported ethnicity
  - `ethnicity_ancestry_ct_enr_heatmap.pdf`; enrichment heatmap of genetic ancestry subpopulations among reported ethnicity groups
  - `lgg_subtype_by_predicted_ancestry.pdf`
  - `lgg_tumor_location_by_predicted_ancestry.pdf`
  - `lgg_tumor_resection_by_predicted_ancestry.pdf`
  - `low_major_ancestry_heatmap.pdf`; heatmap of predicted ancestries probabilities for patients with major (predicted) ancestry probability < 0.9
  - `major_predicted_ancestry_hist.pdf`; histogram of major (predicted) ancestry probability
  - `plot_group_ancestry_ct_enr_heatmap.pdf`; 
  - `plot_group_ancestry_enrichment_heatmap.pdf`
  - `plot_group_by_ancestry.pdf`; plot group count by predicted ancestry 
  - `predicted_ancestry_counts_by_ethnicity.pdf`
  - `predicted_ancestry_counts_by_race.pdf`
  - `predicted_ancestry_percent_by_ethnicity.pdf`
  - `predicted_ancestry_percent_by_race.pdf`
  - `race_ancestry_ct_enr_heatmap.pdf`; enrichment heatmap of genetic ancestry subpopulations among reported race groups

