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

ancestry_file <- file.path(input_dir, "merged_ancestry_histology_survival.tsv")

# Define df of histology names and folder names
hist_df <- data.frame(group = c("Atypical Teratoid Rhabdoid Tumor", "Choroid plexus tumor", 
                                "Craniopharyngioma", "High-grade glioma", 
                                "Ependymoma", "Germ cell tumor", "Mixed neuronal-glial tumor", 
                                "Low-grade glioma", "Medulloblastoma", "Meningioma", 
                                "Mesenchymal tumor", "Neurofibroma plexiform",
                                "Oligodendroglioma", "Schwannoma"),
                      hist = c("ATRT", "CPT", "CRANIO", "HGG", "EPN", "GERM", 
                               "GNG", "LGG", "MB", "MNG", "MES", "NFP", "OLIGO", "SWN"))

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

# Create empty df to store survival model summary stats
survival_anc_stats <- data.frame(row.names = c(group_df$group),
                                 OS_chisq = rep(0, length(group_df$group)),
                                 OS_p = rep(0, length(group_df$group)),
                                 EFS_chisq = rep(0, length(group_df$group)),
                                 EFS_p = rep(0, length(group_df$group)))

# Loop through histologies and molecular subtypes and calculate chi-squared and p-values for effect of predicted ancestry on overall and event-free survival 
for (i in 1:nrow(group_df)){
  
  input_dir <- file.path(analysis_dir, "results", group_df$hist[i])
  
  if (grepl("Low-grade glioma|LGG|Other high-grade glioma|HGG|ATRT|Schwannoma", group_df$group[i])){
    survival_os_anc <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{group_df$subtype[i]}_OS_additive_terms_subtype_resection_predicted_ancestry.RDS")
      ))
  }else{
    survival_os_anc <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{group_df$subtype[i]}_OS_additive_terms_subtype_predicted_ancestry.RDS")
      ))
  }
  
  survival_anc_stats[group_df$group[i],]$OS_chisq <- round(anova(survival_os_anc)["predicted_ancestry", 2], 1)
  survival_anc_stats[group_df$group[i],]$OS_p <- round(anova(survival_os_anc)["predicted_ancestry", 4], 2)
  
  
  if (grepl("Low-grade glioma|LGG|Other high-grade glioma|HGG|ATRT|Schwannoma", group_df$group[i])){
    survival_efs_anc <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{group_df$subtype[i]}_EFS_additive_terms_subtype_resection_predicted_ancestry.RDS")
      ))
  }else{
    survival_efs_anc <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{group_df$subtype[i]}_EFS_additive_terms_subtype_predicted_ancestry.RDS")
      ))
  }
  
  survival_anc_stats[group_df$group[i],]$EFS_chisq <- round(anova(survival_efs_anc)["predicted_ancestry", 2], 1)
  survival_anc_stats[group_df$group[i],]$EFS_p <- round(anova(survival_efs_anc)["predicted_ancestry", 4], 2)
  
}


ancestry <- ancestry %>%
  mutate(plot_group = case_when(
    plot_group %in% c("DIPG or DMG", "Other high-grade glioma") ~ "High-grade glioma",
    TRUE ~ plot_group
  )) %>%
  mutate(predicted_ancestry = factor(predicted_ancestry,
                                     c("AFR", "AMR", "EAS", "EUR", "SAS")))

# Calculate median OS by predicted ancestry and histology
median_os_group <- ancestry %>%
  group_by(plot_group, predicted_ancestry) %>%
  summarize(median_OS_years = round(median(OS_years, na.rm = T), 1)) %>%
  spread(predicted_ancestry, median_OS_years) %>%
  rename(group = plot_group)

# Calculate median OS by predicted ancestry and molecular subgroup
median_os_subgroup <- ancestry %>%
  group_by(mol_sub_group, predicted_ancestry) %>%
  summarize(median_OS_years = round(median(OS_years, na.rm = T), 1)) %>%
  spread(predicted_ancestry, median_OS_years) %>%
  rename(group = mol_sub_group)

# Merge median OS data
median_os <- median_os_group %>%
  bind_rows(median_os_subgroup) %>%
  filter(group %in% rownames(survival_anc_stats)) %>%
  dplyr::rename(AFR_OS_median = AFR,
                AMR_OS_median = AMR,
                EAS_OS_median = EAS,
                EUR_OS_median = EUR,
                SAS_OS_median = SAS)

# Calculate OS standard deviation by predicted ancestry within histologies and molecular subgroups
sd_os_group <- ancestry %>%
  group_by(plot_group, predicted_ancestry) %>%
  summarize(sd_OS_years = round(sd(OS_years, na.rm = T), 1)) %>%
  spread(predicted_ancestry, sd_OS_years) %>%
  rename(group = plot_group)

sd_os_subgroup <- ancestry %>%
  group_by(mol_sub_group, predicted_ancestry) %>%
  summarize(sd_OS_years = round(sd(OS_years, na.rm = T), 1)) %>%
  spread(predicted_ancestry, sd_OS_years) %>%
  rename(group = mol_sub_group)

# merge standard deviations
sd_os <- sd_os_group %>%
  bind_rows(sd_os_subgroup) %>%
  filter(group %in% rownames(survival_anc_stats)) %>%
  dplyr::rename(AFR_OS_sd = AFR,
                AMR_OS_sd = AMR,
                EAS_OS_sd = EAS,
                EUR_OS_sd = EUR,
                SAS_OS_sd = SAS)

# Merge median and standard deviation data
os_merged <- median_os %>%
  ungroup() %>%
  dplyr::mutate(OS_AFR = glue::glue("{AFR_OS_median} ({sd_os$AFR_OS_sd})"),
                OS_AMR = glue::glue("{AMR_OS_median} ({sd_os$AMR_OS_sd})"),
                OS_EAS = glue::glue("{EAS_OS_median} ({sd_os$EAS_OS_sd})"),
                OS_EUR = glue::glue("{EUR_OS_median} ({sd_os$EUR_OS_sd})"),
                OS_SAS = glue::glue("{SAS_OS_median} ({sd_os$SAS_OS_sd})")
                ) %>%
  dplyr::select(group, OS_AFR, OS_AMR, OS_EAS, OS_EUR, OS_SAS)

# Repeat for EFS
median_EFS_group <- ancestry %>%
  group_by(plot_group, predicted_ancestry) %>%
  summarize(median_EFS_years = round(median(EFS_years, na.rm = T), 1)) %>%
  spread(predicted_ancestry, median_EFS_years) %>%
  rename(group = plot_group)

median_EFS_subgroup <- ancestry %>%
  group_by(mol_sub_group, predicted_ancestry) %>%
  summarize(median_EFS_years = round(median(EFS_years, na.rm = T), 1)) %>%
  spread(predicted_ancestry, median_EFS_years) %>%
  rename(group = mol_sub_group)

median_EFS <- median_EFS_group %>%
  bind_rows(median_EFS_subgroup) %>%
  filter(group %in% rownames(survival_anc_stats)) %>%
  dplyr::rename(AFR_EFS_median = AFR,
                AMR_EFS_median = AMR,
                EAS_EFS_median = EAS,
                EUR_EFS_median = EUR,
                SAS_EFS_median = SAS)

sd_EFS_group <- ancestry %>%
  group_by(plot_group, predicted_ancestry) %>%
  summarize(sd_EFS_years = round(sd(EFS_years, na.rm = T), 1)) %>%
  spread(predicted_ancestry, sd_EFS_years) %>%
  rename(group = plot_group)

sd_EFS_subgroup <- ancestry %>%
  group_by(mol_sub_group, predicted_ancestry) %>%
  summarize(sd_EFS_years = round(sd(EFS_years, na.rm = T), 1)) %>%
  spread(predicted_ancestry, sd_EFS_years) %>%
  rename(group = mol_sub_group)

sd_EFS <- sd_EFS_group %>%
  bind_rows(sd_EFS_subgroup) %>%
  filter(group %in% rownames(survival_anc_stats)) %>%
  dplyr::rename(AFR_EFS_sd = AFR,
                AMR_EFS_sd = AMR,
                EAS_EFS_sd = EAS,
                EUR_EFS_sd = EUR,
                SAS_EFS_sd = SAS)

EFS_merged <- median_EFS %>%
  ungroup() %>%
  dplyr::mutate(EFS_AFR = glue::glue("{AFR_EFS_median} ({sd_EFS$AFR_EFS_sd})"),
                EFS_AMR = glue::glue("{AMR_EFS_median} ({sd_EFS$AMR_EFS_sd})"),
                EFS_EAS = glue::glue("{EAS_EFS_median} ({sd_EFS$EAS_EFS_sd})"),
                EFS_EUR = glue::glue("{EUR_EFS_median} ({sd_EFS$EUR_EFS_sd})"),
                EFS_SAS = glue::glue("{SAS_EFS_median} ({sd_EFS$SAS_EFS_sd})")
  ) %>%
  dplyr::select(group, EFS_AFR, EFS_AMR, EFS_EAS, EFS_EUR, EFS_SAS)

# Merge chi-squared test, median and standard deviation survival data
survival_anc_stats <- survival_anc_stats %>%
  rownames_to_column("group") %>%
  left_join(os_merged) %>%
  left_join(EFS_merged) %>%
  dplyr::arrange(group) %>%
  dplyr::select(group, 
                OS_AFR, OS_AMR, OS_EAS, OS_EUR, OS_SAS,
                OS_chisq, OS_p,
                EFS_AFR, EFS_AMR, EFS_EAS, EFS_EUR, EFS_SAS,
                EFS_chisq, EFS_p) %>%
  mutate_all(funs(str_replace_all(., "NA", "--")))

# Save table to output
write_tsv(survival_anc_stats, 
          file.path(analysis_dir, "results", "median-survival-by-ancestry-cancer-group.tsv"))  

