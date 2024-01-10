#!/bin/bash

set -e
set -o pipefail

# prepare cohort-specific histologies file
Rscript -e "rmarkdown::render('01-add_histologies.Rmd')"

# Get summary stats
Rscript 02-summary_stats.R

