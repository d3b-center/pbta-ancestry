# BRAF fusion breakpoint distribution analysis

__Module Authors: Ryan Corbett__ ([@rjcorb](https://github.com/rjcorb))

This module assesses distribution of BRAF fusion breakpoint types across predicted ancestries, and models event-free survival by breakpoint type and predicted ancestry in the PBTA ancestry cohort. 

## Usage
`bash run_module.sh`

## Folder content

- `01-get-breakpoint-exons.Rmd`; annotates putative oncogenic BRAF fusions with exon id and rank in canonincal transcripts

- `02-fusion-breakpoints-by-ancestry.Rmd`; assesses BRAF fusion breakpoint distribution by predicted ancestries, methylation subtypes, and identifies associations between breakpoint type and CNS region and degree of tumor resection

- `03-survival-by-braf-breakpoints.Rmd`; generate event-free survival models to determine affect of breakpoint type and predicted ancestry

- `results/` files: 
  - `coxph_braf_efs*`; Cox proportional hazards model RDS file assessing EFS in all braf fusion or specific breakpoint type patients, and including covariates specified in file name. 
  - `braf-fusions-exon-annotation.tsv`; putative oncogenic BRAF fusion file with exon annotation 
  - `lgg-braf-fusion-breakpoint-annotation.tsv`; LGG BRAF fusion-subset of PBTA ancestry cohort with rare and common breakpoint types reported
  - `lgg-braf-fusion-common-breakpoint-freq.tsv`; frequency of common breakpoint types by predicted ancestry

- `plots/` files:
  - `pbta_ancestry_km_efs*`; Kaplan-meier survival curves in all BRAF fusion patients by breakpoint type or group, or in specific breakpoint types by predicted ancestrty.
  - `pbta_ancestry_forest_efs*`; Forest plots of EFS in all BRAF fusion patients or in specific breakpoint types, and including covariates specified in file name. 


```
.
├── 01-get-breakpoint-exons.Rmd
├── 01-get-breakpoint-exons.html
├── 02-fusion-breakpoints-by-ancestry.Rmd
├── 02-fusion-breakpoints-by-ancestry.html
├── 03-survival-by-braf-breakpoints.Rmd
├── 03-survival-by-braf-breakpoints.html
├── README.md
├── plots
│   ├── pbta_ancestry_forest_efs_breakpoint15_9_resection_ancestry.pdf
│   ├── pbta_ancestry_forest_efs_breakpoint16_9_resection_ancestry.pdf
│   ├── pbta_ancestry_forest_efs_resection_ancestry_breakpoint_group.pdf
│   ├── pbta_ancestry_forest_efs_resection_ancestry_breakpoint_type.pdf
│   ├── pbta_ancestry_forest_efs_resection_breakpoint_group.pdf
│   ├── pbta_ancestry_forest_efs_resection_breakpoint_type.pdf
│   ├── pbta_ancestry_km_efs_15_9_predicted_ancestry.pdf
│   ├── pbta_ancestry_km_efs_16_9_predicted_ancestry.pdf
│   ├── pbta_ancestry_km_efs_breakpoint_group.pdf
│   └── pbta_ancestry_km_efs_breakpoint_type.pdf
├── results
│   ├── braf-fusions-exon-annotation.tsv
│   ├── coxph_braf_efs_15_9_resection_ancestry_breakpoint.RDS
│   ├── coxph_braf_efs_16_9_resection_ancestry_breakpoint.RDS
│   ├── coxph_braf_efs_resection_ancestry_breakpoint_group.RDS
│   ├── coxph_braf_efs_resection_ancestry_breakpoint_type.RDS
│   ├── coxph_braf_efs_resection_breakpoint_group.RDS
│   ├── coxph_braf_efs_resection_breakpoint_type.RDS
│   ├── lgg-braf-fusion-breakpoint-annotation.tsv
│   └── lgg-braf-fusion-common-breakpoint-freq.tsv
└── run_module.sh
```
