#!/bin/bash

set -e
set -o pipefail

# perform exon annotation of BRAF fusion breakpoints
Rscript -e "rmarkdown::render('01-get-breakpoint-exons.Rmd')"

# Breakpoint distribution analysis by predicted ancestry
Rscript -e "rmarkdown::render('02-fusion-breakpoints-by-ancestry.Rmd')"

# Survival by breakpoint type
Rscript -e "rmarkdown::render('03-survival-by-braf-breakpoints.Rmd')"