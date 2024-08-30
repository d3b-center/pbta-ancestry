# R Corbett 2023
#
# Generate prinicpal component survival model summary stats

# Load packages
library(tidyverse)
library(survival)
library(patchwork)
library(colorblindr)
library(circlize)
library(ComplexHeatmap)

# Establish directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "survival-pcs")
input_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

# Load survival functions
source(file.path(root_dir, "analyses", "survival", "util", "survival_models_additive.R"))

ancestry_file <- file.path(root_dir, "analyses", "survival", "results", "merged_ancestry_histology_survival.tsv")

# Define df of histology names and folder names
hist_df <- data.frame(group = c("Atypical Teratoid Rhabdoid Tumor",
                                "Craniopharyngioma", "High-grade glioma", 
                                "Ependymoma", "Germ cell tumor", "Mixed neuronal-glial tumor", 
                                "Low-grade glioma", "Medulloblastoma", "Meningioma", 
                                "Mesenchymal tumor",
                                "Schwannoma"),
                      hist = c("ATRT", "CRANIO", "HGG", "EPN", "GERM", 
                               "GNG", "LGG", "MB", "MNG", "MES", "SWN"))

# Load ancestry file
ancestry <- read_tsv(ancestry_file)

# Create subtype df of molecular subgroup names and folder names
subtype_df <- ancestry %>%
  count(mol_sub_group, plot_group) %>%
  filter(n >=20 & !grepl("classified", mol_sub_group) & !is.na(mol_sub_group)) %>%
  dplyr::mutate(hist = unlist(lapply(strsplit(mol_sub_group, ", "), function(x) x[1])),
                subtype = case_when(
                  grepl(",", mol_sub_group) ~ unlist(lapply(strsplit(mol_sub_group, ", "), function(x) x[2])),
                  TRUE ~ NA_character_)) %>%
  dplyr::mutate(subtype = str_replace_all(subtype, " ", "-")) %>%
  dplyr::mutate(hist = case_when(
    hist == "GNG/GNT" ~ "GNG",
    TRUE ~ hist
  ))

# merge histology and subtype names 
group_df <- subtype_df %>%
  rename(group = mol_sub_group) %>%
  dplyr::mutate(subtype = glue::glue("{hist}_{subtype}")) %>%
  dplyr::select(group, hist, subtype) %>%
  bind_rows(hist_df) %>%
  dplyr::mutate(subtype = case_when(
    is.na(subtype) ~ hist,
    TRUE ~ subtype
  )) %>%
  arrange(group)

# Create empty list to store results
pc_res_list <- list()

# define PC variables
pcs <- c("PC1", "PC2", "PC3", "PC4", "PC5")

for (pc in pcs){

  # Create empty df to store survival model summary stats
  pc_res_list[[pc]] <- data.frame(group = c(group_df$group),
                                   OS_HR = rep(0, length(group_df$group)),
                                   OS_p = rep(0, length(group_df$group)),
                                   EFS_HR = rep(0, length(group_df$group)),
                                   EFS_p = rep(0, length(group_df$group)))
  
  # Loop through histologies and molecular subtypes and calculate extract OS and EFS hazard ratios and associated p-values for PC covariates
  for (i in 1:nrow(group_df)){
    
    input_dir <- file.path(analysis_dir, "results", group_df$hist[i])
    
    if (grepl("Low-grade glioma|LGG", group_df$group[i])){
      survival_os_anc <- read_rds(
        file.path(input_dir,
                  glue::glue("cox_{group_df$subtype[i]}_OS_additive_terms_subtype_resection_{pc}.RDS")
        ))
    }else{
      survival_os_anc <- read_rds(
        file.path(input_dir,
                  glue::glue("cox_{group_df$subtype[i]}_OS_additive_terms_subtype_{pc}.RDS")
        ))
    }

    os_df <- broom::tidy(survival_os_anc)

    # extract OS HR and p value
    pc_res_list[[pc]][pc_res_list[[pc]]$group == group_df$group[i],]$OS_HR <- round(exp(os_df$estimate[os_df$term == pc]), 2)
    pc_res_list[[pc]][pc_res_list[[pc]]$group == group_df$group[i],]$OS_p <- os_df$p.value[os_df$term == pc]
    
    if (grepl("Low-grade glioma|LGG", group_df$group[i])){
      survival_efs_anc <- read_rds(
        file.path(input_dir,
                  glue::glue("cox_{group_df$subtype[i]}_EFS_additive_terms_subtype_resection_{pc}.RDS")
        ))
    }else{
      survival_efs_anc <- read_rds(
        file.path(input_dir,
                  glue::glue("cox_{group_df$subtype[i]}_EFS_additive_terms_subtype_{pc}.RDS")
        ))
    }
    
    efs_df <- broom::tidy(survival_efs_anc)

    # extract EFS HR and p-value
    pc_res_list[[pc]][pc_res_list[[pc]]$group == group_df$group[i],]$EFS_HR <- round(exp(efs_df$estimate[efs_df$term == pc]), 2)
    pc_res_list[[pc]][pc_res_list[[pc]]$group == group_df$group[i],]$EFS_p <- efs_df$p.value[efs_df$term == pc]
    
  }

}

# Create OS HR matrix

merged_res <- Reduce(rbind, pc_res_list)

os_mat <- merged_res %>%
  dplyr::select(group, OS_HR) %>%
  dplyr::mutate(PC = rep(c("PC1", "PC2", "PC3", "PC4", "PC5"), each = nrow(group_df))) %>%
  spread(PC, OS_HR) %>%
  column_to_rownames("group")

os_p_mat <- merged_res %>%
  dplyr::select(group, OS_p) %>%
  dplyr::mutate(PC = rep(c("PC1", "PC2", "PC3", "PC4", "PC5"), each = nrow(group_df))) %>%
  spread(PC, OS_p) %>%
  column_to_rownames("group")

os_p_mat <- as.matrix(os_p_mat)

os_p_mat <- ifelse(os_p_mat >= 0.01, round(os_p_mat, 2),
                   "<0.01")

efs_mat <- merged_res %>%
  dplyr::select(group, EFS_HR) %>%
  dplyr::mutate(PC = rep(c("PC1", "PC2", "PC3", "PC4", "PC5"), each = nrow(group_df))) %>%
  spread(PC, EFS_HR) %>%
  column_to_rownames("group")

efs_p_mat <- merged_res %>%
  dplyr::select(group, EFS_p) %>%
  dplyr::mutate(PC = rep(c("PC1", "PC2", "PC3", "PC4", "PC5"), each = nrow(group_df))) %>%
  spread(PC, EFS_p) %>%
  column_to_rownames("group")

efs_p_mat <- as.matrix(efs_p_mat)

efs_p_mat <- ifelse(efs_p_mat >= 0.01, round(efs_p_mat, 2),
                   "<0.01")

# merge HR and p-valules in matrix
for (i in 1:nrow(efs_mat)){
  
  for (j in 1:ncol(efs_mat)){
    
    efs_p_mat[i,j] <- glue::glue("{efs_mat[i,j]} ({efs_p_mat[i,j]})")
    
    os_mat[i,j] <- ifelse(os_mat[i,j] > 4, NA, os_mat[i,j])
    os_p_mat[i,j] <- ifelse(os_mat[i,j] > 4, "", os_p_mat[i,j])
    
    os_p_mat[i,j] <- glue::glue("{os_mat[i,j]} ({os_p_mat[i,j]})")
    
  }
}

# plot OS and EFS HRs
col_fun = colorRamp2(c(0, 0.5, 1), c("green", "white", "orangered"))

pdf(file.path(plots_dir, "EFS_PC_HR_matrix.pdf"),
    width = 9, height = 7)

Heatmap(efs_mat,
        name = "EFS Hazard Ratio",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "black", lwd = 2),
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", efs_p_mat[i, j]), x, y, gp = gpar(fontsize = 12))
        })

dev.off()

pdf(file.path(plots_dir, "OS_PC_HR_matrix.pdf"),
    width = 9, height = 7)

Heatmap(os_mat,
        name = "OS Hazard Ratio",
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        rect_gp = gpar(col = "black", lwd = 2),
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", os_p_mat[i, j]), x, y, gp = gpar(fontsize = 12))
        })

dev.off()
