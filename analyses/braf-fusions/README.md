# BRAF fusion breakpoint distribution analysis

__Module Authors: Ryan Corbett__ ([@rjcorb](https://github.com/rjcorb))

This module assesses distribution of BRAF fusion breakpoint types across predicted ancestries, and models event-free survival by breakpoint type and predicted ancestry in the PBTA ancestry cohort. 

## Usage
`bash run_module.sh`

## Folder content

- `01-get-breakpoint-exons.Rmd`; annotates putative oncogenic BRAF fusions with exon id and rank in canonincal transcripts

- `02-fusion-breakpoints-by-ancestry.Rmd`; assesses BRAF fusion breakpoint distribution by predicted ancestries, and identifies associations between breakpoint type and CNS region and degree of tumor resection

- `results/` files: 
  - `braf-fusions-exon-annotation.tsv`; putative oncogenic BRAF fusion file with exon annotation 
  - `lgg-braf-fusion-breakpoint-annotation.tsv`; LGG BRAF fusion-subset of PBTA ancestry cohort with rare and common breakpoint types reported
  - `lgg-braf-fusion-common-breakpoint-freq.tsv`; frequency of common breakpoint types by predicted ancestry


