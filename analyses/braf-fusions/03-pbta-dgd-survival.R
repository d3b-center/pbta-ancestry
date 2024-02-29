
library(data.table)
library(tidyverse)
library(survival)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "braf-fusions")
results_dir <- file.path(analysis_dir, "results")
input_dir <- file.path(analysis_dir, "input")
plot_dir <- file.path(analysis_dir, "plots")

source(file.path(root_dir, "analyses", "survival", "util", "survival_models.R"))


fusion_file <- file.path(results_dir, "braf-fusions-exon-annotation.tsv")

fusion_dgd_file <- file.path("/Users/corbettr/Desktop/OpenPedCan-analysis/data/v13/fusion-dgd.tsv.gz")

hist_file <- file.path("/Users/corbettr/Downloads/histologies.tsv")


hist <- read_tsv(hist_file)

hist_braf <- hist %>%
  dplyr::filter(sample_type == "Tumor") %>%
  dplyr::filter(grepl("KIAA1549-BRAF", molecular_subtype) & grepl("LGG", molecular_subtype)) %>%
  distinct(Kids_First_Participant_ID, .keep_all = T)

braf_fusions <- read_tsv(fusion_file) %>%
  dplyr::filter(Kids_First_Participant_ID %in% braf_hist$Kids_First_Participant_ID) %>%
  dplyr::filter(Fusion_Type == "in-frame",
                grepl("duplication", annots)) %>%
  mutate(keep = case_when(FusionName == "KIAA1549--BRAF" & DomainRetainedGene1B %in% c("Partial", "Yes") ~ "yes",
                          FusionName == "BRAF--KIAA1549" & DomainRetainedGene1A %in% c("Partial", "Yes") ~ "yes",
                          TRUE ~ "no"))



braf_fusions_dgd <- read_tsv(fusion_dgd_file) %>%
  dplyr::filter(FusionName == "KIAA1549--BRAF" & !grepl("Tier 3|Tier 4", variant_tier)) %>%
  dplyr::rename(Kids_First_Biospecimen_ID = Sample) %>%
  dplyr::select(Kids_First_Biospecimen_ID, FusionName, `5_prime_region`, `3_prime_region`) %>%
  dplyr::mutate(`5_prime_region` = gsub("^Exon |^exon ", "", `5_prime_region`),
         `3_prime_region` = gsub("^Exon |^exon ", "", `3_prime_region`),                        
         breakpoint_exons = paste(`5_prime_region`, `3_prime_region`, sep = ":"),
         breakpoint_type = case_when(breakpoint_exons %in% c("16:9", "15:9", "16:11", "18:10") ~ "common",
                                     TRUE ~ "rare/novel"),
         keep = "yes") %>%
  left_join(hist[,c("Kids_First_Participant_ID", "Kids_First_Biospecimen_ID")])


common_fusions_pbta <- braf_fusions %>%
  dplyr::filter(breakpoint_exons %in% c("16:9", "15:9", "16:11", "18:10") & keep == "yes") %>%
  distinct(Kids_First_Participant_ID, breakpoint_exons, .keep_all = T) %>%
  group_by(Kids_First_Participant_ID) %>%
  summarise(breakpoint_exons_common = str_c(breakpoint_exons, collapse = ";"))

common_fusions <- braf_fusions_dgd %>%
  dplyr::filter(breakpoint_exons %in% c("16:9", "15:9", "16:11", "18:10") & keep == "yes") %>%
  distinct(Kids_First_Participant_ID, breakpoint_exons, .keep_all = T) %>%
  group_by(Kids_First_Participant_ID) %>%
  summarise(breakpoint_exons_common = str_c(breakpoint_exons, collapse = ";")) %>%
  bind_rows(common_fusions_pbta) %>%
  distinct(Kids_First_Participant_ID, breakpoint_exons_common, .keep_all = T) %>%
  group_by(Kids_First_Participant_ID) %>%
  summarise(breakpoint_exons_common = str_c(breakpoint_exons_common, collapse = ";"))
  

rare_novel_fusions_pbta <- braf_fusions %>%
  dplyr::filter(!breakpoint_exons %in% c("16:9", "15:9", "16:11", "18:10") & keep == "yes") %>%
  distinct(Kids_First_Participant_ID, breakpoint_exons, .keep_all = T) %>%
  group_by(Kids_First_Participant_ID) %>%
  summarise(breakpoint_exons_rare = str_c(breakpoint_exons, collapse = ";"))

rare_novel_fusions <- braf_fusions_dgd %>%
  dplyr::filter(!breakpoint_exons %in% c("16:9", "15:9", "16:11", "18:10") & keep == "yes") %>%
  distinct(Kids_First_Participant_ID, breakpoint_exons, .keep_all = T) %>%
  group_by(Kids_First_Participant_ID) %>%
  summarise(breakpoint_exons_rare = str_c(breakpoint_exons, collapse = ";")) %>%
  bind_rows(rare_novel_fusions_pbta) %>%
  distinct(Kids_First_Participant_ID, breakpoint_exons_rare, .keep_all = T) %>%
  group_by(Kids_First_Participant_ID) %>%
  summarise(breakpoint_exons_rare = str_c(breakpoint_exons_rare, collapse = ";"))


hist_braf <- hist_braf %>%
  left_join(common_fusions) %>%
  left_join(rare_novel_fusions) %>%
  dplyr::filter(!is.na(breakpoint_exons_common) | !is.na(breakpoint_exons_rare))


hist_braf <- hist_braf %>%
  dplyr::mutate(breakpoint_group = case_when(
    !is.na(hist_braf$breakpoint_exons_common) & !is.na(hist_braf$breakpoint_exons_rare) ~ NA_character_,
    breakpoint_exons_common == "16:9" ~ "16:9",
    breakpoint_exons_common == "15:9" ~ "15:9",
    breakpoint_exons_common == "16:11" ~ "16:11",
    !is.na(breakpoint_exons_rare) ~ "rare",
    TRUE ~ NA_character_
  )) %>%
  dplyr::mutate(breakpoint_group = fct_relevel(breakpoint_group,
                                               c("16:11", "16:9", "15:9", "rare")
                                               )) %>%
  dplyr::mutate(breakpoint_type = case_when(
    !is.na(hist_braf$breakpoint_exons_common) & !is.na(hist_braf$breakpoint_exons_rare) ~ NA_character_,
    !is.na(breakpoint_exons_rare) ~ "rare",
    !is.na(breakpoint_exons_common) ~ "common",
    TRUE ~ NA_character_
  )) %>%
  dplyr::mutate(EFS_status = case_when(
    !is.na(EFS_event_type) & EFS_event_type != "Not Applicable" ~ "EVENT",
    TRUE ~ "NO EVENT"
  )) %>%
  dplyr::mutate(EFS_years = as.numeric(EFS_days)/365.25) %>%
  dplyr::mutate(is_cerebellum = case_when(
    grepl("Cerebellum", primary_site) ~ "Yes",
    TRUE ~ "No"
  ))


kap_efs <- survival_analysis(
  metadata  = hist_braf,
  ind_var = "breakpoint_group",
  test = "kap.meier",
  metadata_sample_col = "Kids_First_Biospecimen_ID", 
  days_col = "EFS_days",
  status_col = "EFS_status"
)

km_plot <- plotKM(model = kap_efs,
                  variable = "breakpoint_group",
                  combined = F, 
                  title = "LGG, BRAF fusion")

km_plot

ggsave(file.path(plot_dir, "opc_km_efs_breakpoint_group.pdf"),
       width = 9, height = 5)

kap_efs <- survival_analysis(
  metadata  = hist_braf,
  ind_var = "breakpoint_type",
  test = "kap.meier",
  metadata_sample_col = "Kids_First_Biospecimen_ID", 
  days_col = "EFS_days",
  status_col = "EFS_status"
)

km_plot <- plotKM(model = kap_efs,
                  variable = "breakpoint_type",
                  combined = F, 
                  title = "LGG, BRAF fusion")

km_plot

ggsave(file.path(plot_dir, "opc_km_efs_breakpoint_type.pdf"),
       width = 9, height = 5)


braf_model_os <- fit_save_model(hist_braf[hist_braf$extent_of_tumor_resection != "Not Reported",],
                                terms = "extent_of_tumor_resection+breakpoint_group",
                                file.path(results_dir, "coxph_opc_braf_efs.RDS"),
                                "multivariate",
                                years_col = "EFS_years", status_col = "EFS_status"
)

plotForest(readRDS(file.path(results_dir, "coxph_opc_braf_efs.RDS")))

ggsave(file.path(plot_dir, "opc_forest_efs_resection_breakpoint_group.pdf"),
       width = 8, height = 3)



braf_model_os <- fit_save_model(hist_braf[hist_braf$extent_of_tumor_resection != "Not Reported",],
                                terms = "extent_of_tumor_resection+breakpoint_type",
                                file.path(results_dir, "coxph_opc_braf_efs.RDS"),
                                "multivariate",
                                years_col = "EFS_years", status_col = "EFS_status"
)

plotForest(readRDS(file.path(results_dir, "coxph_opc_braf_efs.RDS")))

ggsave(file.path(plot_dir, "opc_forest_efs_resection_breakpoint_type.pdf"),
       width = 8, height = 3)
