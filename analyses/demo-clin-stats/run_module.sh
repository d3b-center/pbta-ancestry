#!/bin/bash

set -e
set -o pipefail

# generate demo clin summary stat tables
Rscript -e "rmarkdown::render('01-demo-clin-stats.Rmd')"


