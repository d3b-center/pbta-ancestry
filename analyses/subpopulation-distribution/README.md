# Assess genetic ancestry subpopulation probability distributions in PBTA ancestry cohort

This module utilizes Somalier-derived [genetic ancestry subpopulation](https://www.nature.com/articles/nature15393) probabilities to assess their distribution across the PBTA ancestry cohort

## Usage

`bash run_module.sh`

## Folder content 

1. `01-subpopulation_distribution.R` plots genetic ancestry subpopulation probabilities for each superpopulation group (AFR, AMR, EAS, EUR, SAS)

##Analysis module directory structure

```
.
├── 01-subpopulation_distribution.R
├── README.md
├── input
│   └── PBTA_population.somalier-ancestry.tsv
├── plots
│   ├── AFR_subpopulation_probablilties.pdf
│   ├── AMR_subpopulation_probablilties.pdf
│   ├── EAS_subpopulation_probablilties.pdf
│   ├── EUR_subpopulation_probablilties.pdf
│   └── SAS_subpopulation_probablilties.pdf
├── results
│   └── pbta-somalier-subpopulations.tsv
└── run_module.sh
```