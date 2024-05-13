# R Corbett 2023
#
# Generate Kaplan-Meier survival curves and forest plots for PBTA ancestry cohort

# Load packages
library(tidyverse)
library(survival)
library(patchwork)
library(colorblindr)

# Define directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "survival")
input_dir <- file.path(analysis_dir, "results")
plots_dir <- file.path(analysis_dir, "plots")

# Load survival functions
source(file.path(analysis_dir, "util", "survival_models_additive.R"))

ancestry_file <- file.path(input_dir, "merged_ancestry_histology_survival.tsv")

# Define df of histology names and folder names
hist_df <- data.frame(group = c("ATRT", "Choroid plexus tumor",
                                "Craniopharyngioma", "HGG", 
                                 "Ependymoma", "Germ cell tumor", "GNG", "LGG", 
                                 "Medulloblastoma", "Meningioma", 
                                 "Mesenchymal tumor", "Neurofibroma Plexiform",
                                 "Oligodendroglioma", "Schwannoma"),
                       hist = c("ATRT", "CPT", "CRANIO", "HGG", "EPN", "GERM",
                                "GNG", "LGG", "MB", "MNG", "MES", "NFP",
                                  "OLIGO", "SWN"))

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
  dplyr::rename(group = mol_sub_group) %>%
  dplyr::mutate(subtype = glue::glue("{hist}_{subtype}")) %>%
  dplyr::select(group, hist, subtype) %>%
  bind_rows(hist_df) %>%
  dplyr::mutate(subtype = case_when(
    is.na(subtype) ~ hist,
    TRUE ~ subtype
  )) %>%
  arrange(group)
  
# define variables
variables <- c("predicted_ancestry", "race", "ethnicity", "EUR_status")

pdf(NULL)

# Loop through variables and histologies/subtypes to plot KM survival curves and hazard ratios
for (var in variables){
  
  for (i in 1:nrow(group_df)){
    
    # Define histology specific results and plots directories
    input_dir <- file.path(analysis_dir, "results", group_df$hist[i])
    
    plots_dir <- file.path(analysis_dir, "plots", group_df$hist[i])
    
    if (!dir.exists(plots_dir)) {
      dir.create(plots_dir)
      
    }
    
    # Load KM model files
    km_os_result <- read_rds(
      file.path(input_dir,
                glue::glue("logrank_{group_df$subtype[i]}_OS_{var}.RDS")
      )
    )
    
    km_efs_result <- read_rds(
      file.path(input_dir,
                glue::glue("logrank_{group_df$subtype[i]}_EFS_{var}.RDS")
      )
    )
    
    # Define output files
    km_os_output_pdf <- file.path(plots_dir, glue::glue("km_{group_df$subtype[i]}_OS_{var}.pdf"))
    km_efs_output_pdf <- file.path(plots_dir, glue::glue("km_{group_df$subtype[i]}_EFS_{var}.pdf"))
    
    # Generate KM plots
    km_os_plot <- plotKM(model = km_os_result,
                         variable = var,
                         combined = F, 
                         title = group_df$group[i])
    
    ggsave(km_os_output_pdf, km_os_plot,
           width = 10, height = 7, units = "in",
           device = "pdf")
    
    
    km_efs_plot <- plotKM(model = km_efs_result,
                          variable = var,
                          combined = F, 
                          title = group_df$group[i])
    
    ggsave(km_efs_output_pdf, km_efs_plot,
           width = 10, height = 7, units = "in",
           device = "pdf")
    
    # Load cox proportional hazards models and plot
    if (grepl("Low-grade glioma|LGG", group_df$group[i])){
      os_survival_result <- read_rds(
        file.path(input_dir,
                  glue::glue("cox_{group_df$subtype[i]}_OS_additive_terms_subtype_resection_{var}.RDS")
        ))
    }else{
      os_survival_result <- read_rds(
        file.path(input_dir,
                  glue::glue("cox_{group_df$subtype[i]}_OS_additive_terms_subtype_{var}.RDS")
        ))
    }
    
    if (sum(is.na(broom::tidy(os_survival_result)$estimate)) != nrow(broom::tidy(os_survival_result))) {
    
      forest_pdf <- file.path(plots_dir, 
                              glue::glue("forest_{group_df$subtype[i]}_OS_subtype_{var}.pdf"))
      
      
      forest_plot <- try(plotForest(os_survival_result), silent = TRUE)
      
      ggsave(forest_pdf, forest_plot, width = 8, height = 3)
    
    }
    
    
    if (grepl("Low-grade glioma|LGG", group_df$group[i])){
      efs_survival_result <- read_rds(
        file.path(input_dir,
                  glue::glue("cox_{group_df$subtype[i]}_EFS_additive_terms_subtype_resection_{var}.RDS")
        ))
    }else{
      efs_survival_result <- read_rds(
        file.path(input_dir,
                  glue::glue("cox_{group_df$subtype[i]}_EFS_additive_terms_subtype_{var}.RDS")
        ))
    }
    
    if (sum(is.na(broom::tidy(efs_survival_result)$estimate)) != nrow(broom::tidy(efs_survival_result))) {
      
      forest_pdf <- file.path(plots_dir, 
                              glue::glue("forest_{group_df$subtype[i]}_EFS_subtype_{var}.pdf"))
      
      
      forest_plot <- plotForest(efs_survival_result)
      
      ggsave(forest_pdf, forest_plot, width = 8, height = 3)
      
    }
    
  }
  
}
