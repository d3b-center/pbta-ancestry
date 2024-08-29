# R Corbett 2023
#
# Generate survival summary table for PBTA Ancestry project

# Load packages
library(tidyverse)
library(survival)
library(patchwork)
library(colorblindr)

# Establish directory paths
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "survival")
input_dir <- file.path(analysis_dir, "results")

# Load survival functions
source(file.path(analysis_dir, "util", "survival_models_additive.R"))

# set file paths
ancestry_file <- file.path(input_dir, "merged_ancestry_histology_survival.tsv")

# Define df of histology names and folder names
hist_key <- data.frame(plot_group = c("Atypical Teratoid Rhabdoid Tumor", "Choroid plexus tumor", 
                                      "Craniopharyngioma", "Other high-grade glioma", 
                                      "Ependymoma", "Germ cell tumor", "Mixed neuronal-glial tumor", 
                                      "Low-grade glioma", "Medulloblastoma", "Meningioma", 
                                      "Mesenchymal tumor", "Neurofibroma plexiform",
                                      "Oligodendroglioma", "Schwannoma"),
                       hist = c("ATRT", "CPT", "CRANIO", "HGG", "EPN", "GERM", 
                                "GNG", "LGG", "MB", "MNG", "MES", "NFP", "OLIGO", "SWN"))

# Load ancestry file
ancestry <- read_tsv(ancestry_file)

# Generate summary tables for both EFS and OS in a for-loop
for (surv in c("EFS", "OS")){
  
  # define correct status column
  status <- ifelse(surv == "EFS", "EFS_status", "OS_status")
  
  # filter plot_groups based on # events
  hist_df <- ancestry %>%
    dplyr::count(plot_group, !!sym(status)) %>%
    filter(!!sym(status) %in% c("DECEASED", "EVENT"),
           n >=10,
           # remove histologies with no survival models (DMG K28 was run separate from DIPG/DMG)
           !plot_group %in% c("Other tumor", "DIPG or DMG",
                              "Non-neoplastic tumor", "Other CNS embryonal tumor")) %>%
    left_join(hist_key) %>%
    dplyr::rename(group = plot_group) %>%
    dplyr::select(group, hist)
  
  # Create subtype df of molecular subgroup names and folder names
  subtype_df <- ancestry %>%
    dplyr::count(mol_sub_group, plot_group, !!sym(status)) %>%
    filter(!!sym(status) %in% c("DECEASED", "EVENT"),
           n >=10,
           !grepl("classified", mol_sub_group),
           !is.na(mol_sub_group),
           # no ATRT subtype models run due to low overall N
           !grepl("ATRT", mol_sub_group)) %>%
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
    arrange(group) %>%
    # No GNG/GNT WT models run due to low overall N
    dplyr::filter(group != "GNG/GNT, WT")
  
  # Create empty df to store survival model summary stats
  survival_anc_stats <- data.frame(row.names = c(group_df$group),
                                   chisq_superpop = rep(0, length(group_df$group)),
                                   p_superpop = rep(0, length(group_df$group)),
                                   HR_nonEUR = rep(0, length(group_df$group)),
                                   HR_p_nonEUR = rep(0, length(group_df$group)))
  
  # Loop through histologies and molecular subtypes and extract chi-squared and p-values for effect of predicted ancestry on survival 
  for (i in 1:nrow(group_df)){
    
    # set histology-specific input directory
    input_dir <- file.path(analysis_dir, "results", group_df$hist[i])
    
    # set file survival model file paths
    if (grepl("Low-grade glioma|LGG|Other high-grade glioma|HGG|ATRT|SWN|MB", group_df$hist[i])){
      
      survival_anc <- read_rds(
        file.path(input_dir,
                  glue::glue("cox_{group_df$subtype[i]}_{surv}_additive_terms_subtype_resection_predicted_ancestry.RDS")
        ))
      
      survival_eur <- read_rds(
        file.path(input_dir,
                  glue::glue("cox_{group_df$subtype[i]}_{surv}_additive_terms_subtype_resection_EUR_status.RDS")
        ))
      survival_eur <- broom::tidy(survival_eur)
      
    }else{
      
      survival_anc <- read_rds(
        file.path(input_dir,
                  glue::glue("cox_{group_df$subtype[i]}_{surv}_additive_terms_subtype_predicted_ancestry.RDS")
        ))
      
      survival_eur <- read_rds(
        file.path(input_dir,
                  glue::glue("cox_{group_df$subtype[i]}_{surv}_additive_terms_subtype_EUR_status.RDS")
        ))
      survival_eur <- broom::tidy(survival_eur)
      
    }
    
    # extract chi-squared a p-values from predicted ancestry coxph models
    survival_anc_stats[group_df$group[i],]$chisq_superpop <- round(anova(survival_anc)["predicted_ancestry", 2], 1)
    survival_anc_stats[group_df$group[i],]$p_superpop <- round(anova(survival_anc)["predicted_ancestry", 4], 2)
    
    # extract HRs and p-values from EUR status coxph models
    survival_anc_stats[group_df$group[i],]$HR_nonEUR <- survival_eur %>%
      dplyr::filter(term == "EUR_statusnon-EUR") %>%
      dplyr::mutate(HR = round(exp(estimate), 2)) %>%
      pull(HR)
    
    survival_anc_stats[group_df$group[i],]$HR_p_nonEUR <- survival_eur %>%
      dplyr::filter(term == "EUR_statusnon-EUR") %>%
      dplyr::mutate(p.value = round(p.value, 2)) %>%
      pull(p.value)
    
  }
    
  # calculaate FDRs
  survival_anc_stats <- survival_anc_stats %>%
    dplyr::mutate(fdr_superpop = round(p.adjust(p_superpop, method = "BH"), 2),
                  HR_fdr_nonEUR = round(p.adjust(HR_p_nonEUR, method = "BH"), 2))
  
  # define survival variable
  surv_var <- ifelse(surv == "EFS", "EFS_years", "OS_years")

  # Calculate median survival by predicted ancestry and plot group
  median_surv_group <- ancestry %>%
    dplyr::filter(plot_group %in% rownames(survival_anc_stats)) %>%
    group_by(plot_group, predicted_ancestry) %>%
    summarize(median_surv_years = round(median(!!sym(surv_var), na.rm = T), 1)) %>%
    spread(predicted_ancestry, median_surv_years) %>%
    rename(group = plot_group)

  # Calculate median survival by predicted ancestry and molecular subgroup
  median_surv_subgroup <- ancestry %>%
    dplyr::filter(mol_sub_group %in% rownames(survival_anc_stats)) %>%
    group_by(mol_sub_group, predicted_ancestry) %>%
    summarize(median_surv_years = round(median(!!sym(surv_var), na.rm = T), 1)) %>%
    spread(predicted_ancestry, median_surv_years) %>%
    rename(group = mol_sub_group)
  
  # Merge median survival data
  median_surv <- median_surv_group %>%
    bind_rows(median_surv_subgroup) %>%
    dplyr::rename(median_AFR = "AFR",
                  median_AMR = "AMR",
                  median_EAS = "EAS",
                  median_EUR = "EUR",
                  median_SAS = "SAS")
  
  # Calculate survival standard deviation by predicted ancestry within histologies and molecular subgroups
  sd_surv_group <- ancestry %>%
    dplyr::filter(plot_group %in% rownames(survival_anc_stats)) %>%
    group_by(plot_group, predicted_ancestry) %>%
    summarize(sd_surv_years = round(sd(!!sym(surv_var), na.rm = T), 1)) %>%
    spread(predicted_ancestry, sd_surv_years) %>%
    rename(group = plot_group)
  
  sd_surv_subgroup <- ancestry %>%
    dplyr::filter(mol_sub_group %in% rownames(survival_anc_stats)) %>%
    group_by(mol_sub_group, predicted_ancestry) %>%
    summarize(sd_surv_years = round(sd(!!sym(surv_var), na.rm = T), 1)) %>%
    spread(predicted_ancestry, sd_surv_years) %>%
    rename(group = mol_sub_group)
  
  # merge standard deviations
  sd_surv <- sd_surv_group %>%
    bind_rows(sd_surv_subgroup) %>%
    filter(group %in% rownames(survival_anc_stats)) %>%
    dplyr::rename(sd_AFR = AFR,
                  sd_AMR = AMR,
                  sd_EAS = EAS,
                  sd_EUR = EUR,
                  sd_SAS = SAS)
  
   # merge medians and sds 
  surv_merged <- median_surv %>%
    ungroup() %>%
    dplyr::mutate(median_AFR = glue::glue("{median_AFR} ({sd_surv$sd_AFR})"),
                  median_AMR = glue::glue("{median_AMR} ({sd_surv$sd_AMR})"),
                  median_EAS = glue::glue("{median_EAS} ({sd_surv$sd_EAS})"),
                  median_EUR = glue::glue("{median_EUR} ({sd_surv$sd_EUR})"),
                  median_SAS = glue::glue("{median_SAS} ({sd_surv$sd_SAS})")
    ) %>%
    dplyr::select(group, median_AFR, median_AMR, median_EAS, median_EUR, median_SAS)
  
  # merged all survival summary stats data
  survival_anc_stats <- survival_anc_stats %>%
    rownames_to_column("group") %>%
    # add median and sd survival
    left_join(surv_merged) %>%
    dplyr::arrange(group) %>%
    dplyr::select(group, 
                  median_AFR, median_AMR, median_EAS, 
                  median_EUR, median_SAS,
                  chisq_superpop, p_superpop, fdr_superpop,
                  HR_nonEUR, HR_p_nonEUR, HR_fdr_nonEUR) %>%
    # replace NAs with "--"
    mutate_all(funs(str_replace_all(., "NA", "--")))
  
  # Save table to output
  write_tsv(survival_anc_stats, 
            file.path(analysis_dir, "results", 
                      glue::glue("median-{surv}-by-ancestry-cancer-group.tsv"))) 
  
}

sessionInfo()