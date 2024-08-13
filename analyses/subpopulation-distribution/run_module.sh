#!/bin/bash

set -e
set -o pipefail

# Generate subpopulation barplots
Rscript 01-subpopulation_distribution.R
