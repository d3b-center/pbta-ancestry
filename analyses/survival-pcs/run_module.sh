#!/bin/bash

set -e
set -o pipefail

# prepare survival data
Rscript -e "rmarkdown::render('01-run_survival_pcs.Rmd')"

# Run survival models
Rscript --vanilla 02-survival-summary.R

# Plot survival models
Rscript --vanilla 03-plot-survival.R

# Generate survival summary table
Rscript --vanilla 04-pc-ancestry-correlations.R	
