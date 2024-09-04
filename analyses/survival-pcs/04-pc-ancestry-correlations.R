# R Corbett 2024
#
# Generate heatmap of somalier-derived PCs and genetic ancestry superpopulation probabilities 

# Load packages
library(ComplexHeatmap)
library(tidyverse)
library(circlize)

# Define directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "survival-pcs")
input_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

ancestry_file <- file.path(root_dir, "analyses",
                           "add-histologies", "results",
                           "merged_ancestry_histology_data.tsv")

# Wrangle data
ancestry <- read_tsv(ancestry_file)

# create empty correlation matrix
cor_mat <- matrix(0, 5, 5,
                  dimnames = list(c("AFR_prob", "AMR_prob", "EAS_prob", "EUR_prob", "SAS_prob"),
                                  c("PC1", "PC2", "PC3", "PC4", "PC5")))

# loop through PCs and superpopulations to calculate pearson correlation coefficients
for (anc in rownames(cor_mat)){
  
  for (pc in colnames(cor_mat)){
    
    cor_mat[anc, pc] = round(cor.test(unlist(ancestry[,anc]),
                            unlist(ancestry[,pc]))$estimate, 2)
    
  }
  
}

# plot correlation heatmap
pdf(file.path(plots_dir, "pc-ancestry-prob-correlation-mat.pdf"),
    width = 6, height = 4)

Heatmap(cor_mat,
        name = "Pearson correlation",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "black", lwd = 2),
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", cor_mat[i, j]), x, y, gp = gpar(fontsize = 12))
        })

dev.off()
