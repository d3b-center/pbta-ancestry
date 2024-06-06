# BRAF fusion breakpoint distribution analysis

__Module Authors: Ryan Corbett__ ([@rjcorb](https://github.com/rjcorb))

This module assesses distribution of BRAF fusion breakpoint types across predicted ancestries, and models event-free survival by breakpoint type and predicted ancestry in the PBTA ancestry cohort. 

## Usage
`bash run_module.sh`

## Folder content

- `01-get-breakpoint-exons.Rmd`; annotates putative oncogenic BRAF fusions with exon id and rank in canonincal transcripts

- `02-fusion-breakpoints-by-ancestry.Rmd`; assesses BRAF fusion breakpoint distribution by predicted ancestries, methylation subtypes, and identifies associations between breakpoint type and CNS region and degree of tumor resection

- `03-survival-by-braf-breakpoints.Rmd`; generate event-free survival models to determine affect of breakpoint type and predicted ancestry

## Directory structure

```
.
├── 01-get-breakpoint-exons.Rmd
├── 01-get-breakpoint-exons.html
├── 02-fusion-breakpoints-by-ancestry.Rmd
├── 02-fusion-breakpoints-by-ancestry.html
├── 03-survival-by-braf-breakpoints.Rmd
├── 03-survival-by-braf-breakpoints.html
├── README.md
├── input
│   └── lgg_braf_fusion_final_diagnosis.txt
├── plots
│   ├── breakpoint_group_ancestry_ct_enr_heatmap.pdf
│   ├── breakpoint_group_methyl_subtype_ct_enr_heatmap.pdf
│   ├── breakpoint_group_region_ct_enr_heatmap.pdf
│   ├── breakpoint_group_resection_ct_enr_heatmap.pdf
│   ├── lgg_braf_fusion_breakpoint_group_diagnosis_enr_heatmap.pdf
│   ├── pbta_ancestry_forest_efs_add_breakpoint15_9_resection_ancestry.pdf
│   ├── pbta_ancestry_forest_efs_add_breakpoint16_9_resection_ancestry.pdf
│   ├── pbta_ancestry_forest_efs_add_resection_ancestry_breakpoint_group.pdf
│   ├── pbta_ancestry_forest_efs_add_resection_ancestry_breakpoint_type.pdf
│   ├── pbta_ancestry_forest_efs_add_resection_breakpoint_group.pdf
│   ├── pbta_ancestry_forest_efs_add_resection_breakpoint_type.pdf
│   ├── pbta_ancestry_forest_efs_int_resection_breakpoint_group.pdf
│   ├── pbta_ancestry_forest_efs_int_resection_breakpoint_type.pdf
│   ├── pbta_ancestry_km_efs_15_9_predicted_ancestry.pdf
│   ├── pbta_ancestry_km_efs_16_9_predicted_ancestry.pdf
│   ├── pbta_ancestry_km_efs_breakpoint_group.pdf
│   ├── pbta_ancestry_km_efs_breakpoint_type.pdf
│   └── rare_breakpoints_by_ancestry.pdf
├── results
│   ├── braf-fusion-breakpoints-by-patient.tsv
│   ├── braf-fusions-exon-annotation.tsv
│   ├── coxph_braf_efs_add_15_9_resection_ancestry_breakpoint.RDS
│   ├── coxph_braf_efs_add_16_9_resection_ancestry_breakpoint.RDS
│   ├── coxph_braf_efs_add_resection_ancestry_breakpoint_group.RDS
│   ├── coxph_braf_efs_add_resection_ancestry_breakpoint_type.RDS
│   ├── coxph_braf_efs_add_resection_breakpoint_group.RDS
│   ├── coxph_braf_efs_add_resection_breakpoint_type.RDS
│   ├── coxph_braf_efs_int_resection_breakpoint_group.RDS
│   ├── coxph_braf_efs_int_resection_breakpoint_type.RDS
│   ├── lgg-braf-fusion-breakpoint-annotation.tsv
│   └── lgg-braf-fusion-common-breakpoint-freq.tsv
├── run_module.sh
└── util
    └── heatmap_function.R
```
