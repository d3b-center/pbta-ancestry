#!/bin/bash

set -e
set -o pipefail

# prepare survival data
Rscript -e "rmarkdown::render('01-prepare-survival.Rmd')"

# Run survival models
Rscript -e "rmarkdown::render('02-run-survival.Rmd')"

# Plot survival models
Rscript --vanilla 03-plot-survival.R

# Generate survival summary table
Rscript --vanilla 04-survival-summary.R	
