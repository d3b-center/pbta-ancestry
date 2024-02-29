#!/bin/bash

set -e
set -o pipefail

# prepare cohort-specific histologies file
Rscript -e "rmarkdown::render('01-add_histologies.Rmd')"

# Get summary stats
Rscript --vanilla 02-summary_stats.R

# Generate alluvial plot
Rscript 03-alluvial_plot.R