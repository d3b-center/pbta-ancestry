# Run survival analyses in PBTA ancestry cohort

This module formats data for survival analyses and assesses the significance of predicted ancestry, reported race, and reported ethnicity on overall (OS) and event-free survival (EFS) by histology group and molecular subtype. 

## Usage

`bash run_module.sh`

## Folder contents

1. `01-prepare-survival.Rmd	` pulls latest survival data from base histologies file, formats OS and EFS, and groups molecular subtypes into broader molecular subgroups for survival analyses.  

2. `02-run-survival.Rmd` Runs Kaplan Meier and cox proportional hazards models using predicted ancestry, reported race, and reported ethnicity as predictors, and including covariates molecular subgroup, extent of tumor resection (when applicable), and age at diagnosis. 

3. `03-plot-survival.R` plots survival models 

4. `04-survival-summary.R`generates summary table displaying mean OS and EFS by predicted ancestry and histology group/molecular subtype. 

5. `results/` directory contains the following files: 
- All survival models are saved to the following histology-specific folders: 
  - `atrt/` 
  - `cranio/`
  - `CPT/`
  - `DMG/`
  - `epn/`
  - `GERM/`
  - `gng/`
  - `hgg/`
  - `lgg/`
  - `mb/`
  - `mes/`
  - `mng/`
  - `nfp/`
  - `oligo/`
  - `swn/`
- `merged_ancestry_histology_survival.tsv` contains most up-to-date survival data for each patient
- `median-survival-by-ancestry-cancer-group.tsv` contains OS and EFS summary stats by predicted ancestry and histology group

6. `plots/` directory contains the following folders and files: 
- All survival plots are saved to the following histology-specific folders: 
  - `atrt/` 
  - `cranio/`
  - `CPT/`
  - `dmg/`
  - `epn/`
  - `GERM/`
  - `GNG/`
  - `HGG/`
  - `LGG/`
  - `MB/`
  - `MES/`
  - `MNG/`
  - `NFP/`
  - `OLIGO/`
  - `SWN/`
  - `*_subtype_by_predicted_ancestry.pdf`; heatmaps of no. samples by predicted ancestry and molecular subtype within histology groups. 
  
7. `util/` directory contains `survival_models.R` script that contains functions for generating and plotting survival models. 
