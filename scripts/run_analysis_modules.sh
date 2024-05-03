#!/bin/sh

set -e
set -o pipefail

# Set the working directory to the directory of this file
cd "$(dirname "${BASH_SOURCE[0]}")"

# Get base directory of project
cd ..
BASEDIR="$(pwd)"
cd -

analyses_dir="$BASEDIR/analyses"

# Run add-histologies analysis module
echo "Run add-histologies"
cd ${analyses_dir}/add-histologies
bash run_module.sh

# Run demo-clin-stats analysis module
echo "Run demo-clin-stats"
cd ${analyses_dir}/demo-clin-stats
bash run_module.sh

# Run survival analysis module
echo "Run survival"
cd ${analyses_dir}/survival
bash run_module.sh

# Run BRAF breakpoint analysis module
echo "Run BRAF fusion breakpoint analysis"
cd ${analyses_dir}/braf-fusions
bash run_module.sh

printf "\nDone running analysis modules!\n"