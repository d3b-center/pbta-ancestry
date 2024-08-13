# Assess subpopulation distribution among genetic ancestry superpopulations in PBTA cohort

library(tidyverse)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set file paths
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "subpopulation-distribution")
results_dir <- file.path(analysis_dir, "results")
input_dir <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")

# call plot theme 
source(file.path(root_dir, "figures", "theme.R"))

# Set file paths

subpopulation_file <- file.path(input_dir, "PBTA_population.somalier-ancestry.tsv")

hist_file <- file.path(data_dir, "histologies.tsv")

ancestry_file <- file.path(root_dir, "analyses",
                           "add-histologies", "results",
                           "merged_ancestry_histology_data.tsv")

# Wrangle data

ancestry <- read_tsv(ancestry_file)

opc_hist <- read_tsv(hist_file) 

subpop_df <- read_tsv(subpopulation_file) %>%
  dplyr::filter(is.na(given_ancestry)) %>%
  # sample IDs are OPC aliquot IDs
  dplyr::rename(aliquot_id = `#sample_id`) %>%
  dplyr::mutate(aliquot_id = str_replace(aliquot_id, "_D1", "")) %>%
  # Get corresponding BS IDs from OPC hist
  dplyr::mutate(Kids_First_Biospecimen_ID = case_when(
    grepl("BS", aliquot_id) ~ aliquot_id,
    TRUE ~ opc_hist$Kids_First_Biospecimen_ID[match(aliquot_id, opc_hist$aliquot_id)]
  ))

# Create subpopulation-superpopulation key data frame

afr_subpopulations <- c("ASW", "ACB", "ESN", "GWD", "LWK", "MSL", "YRI")
amr_subpopulations <- c("CLM", "MXL", "PUR", "PEL")
eas_subpopulations <- c("CDX", "CHB", "CHS", "JPT", "KHV")
eur_subpopulations <- c("CEU", "GBR", "FIN", "IBS", "TSI")
sas_subpopulations <- c("BEB", "GIH", "ITU", "PJL", "STU")

subpop_key <- data.frame(subpopulation = c(afr_subpopulations,
                                             amr_subpopulations,
                                             eas_subpopulations,
                                             eur_subpopulations,
                                             sas_subpopulations),
                         superpopulation = c(rep("AFR", length(afr_subpopulations)),
                                             rep("AMR", length(amr_subpopulations)),
                                             rep("EAS", length(eas_subpopulations)),
                                             rep("EUR", length(eur_subpopulations)),
                                             rep("SAS", length(sas_subpopulations))))

# get list of PNOC patients

pnoc_pts <- opc_hist %>%
  dplyr::filter(sub_cohort == "PNOC") %>%
  pull(Kids_First_Biospecimen_ID)

# Loop through superpopulations and plot subpopulation composition among patients

populations <- c("AFR", "AMR", "EAS", "EUR", "SAS")

for (pop in populations){

  # extract patients assigned to given superpopulation
  superpop_pts <- ancestry %>%
    dplyr::filter(predicted_ancestry == pop) %>%
    pull(Kids_First_Biospecimen_ID)
  
  # pull all subpopulations associated with superpopulation
  subpopulations <- subpop_key %>%
    dplyr::filter(superpopulation == pop) %>%
    pull(subpopulation)

  # Create filtered df of subpopulation probabilities, concatenating probabilities of subpopulations associated with other superpopulations
  pop_subpop_df <- subpop_df %>%
    dplyr::filter(Kids_First_Biospecimen_ID %in% superpop_pts) %>%
    dplyr::select(-aliquot_id, -given_ancestry, -predicted_ancestry,
                  -PC1, -PC2, -PC3, -PC4, -PC5) %>%
    tidyr::gather(-Kids_First_Biospecimen_ID, key = "subpopulation", value = "probability") %>%
    dplyr::mutate(subpopulation = str_replace(subpopulation, "_prob", "")) %>%
    dplyr::mutate(probability = as.double(probability)) %>%
    dplyr::mutate(probability = case_when(
      is.na(probability) ~ 0,
      TRUE ~ probability
    )) %>%
    # Define all non-superpopulation subpopulations as "other" to be merged 
    dplyr::mutate(subpopulation = case_when(
      !subpopulation %in% subpopulations ~ glue::glue("non-{pop}"),
      TRUE ~ subpopulation
    )) %>%
    group_by(Kids_First_Biospecimen_ID, subpopulation) %>%
    # Sum non-superpopulation probabilities 
    summarise(probability = sum(probability)) %>%
    dplyr::mutate(cohort = case_when(
      Kids_First_Biospecimen_ID %in% pnoc_pts ~ "PNOC",
      TRUE ~ "CBTN"
    ))
  
  # Define color palette 
  color_palette <- c(brewer.pal(n = length(subpopulations), name = "Dark2"),
                     "#999999")

  if (pop == "AMR"){
    
    order_pts <- pop_subpop_df %>%
      dplyr::filter(subpopulation %in% c("MXL", "PEL")) %>%
      group_by(Kids_First_Biospecimen_ID) %>%
      summarize(probability = sum(probability)) %>%
      arrange(probability) %>%
      pull(Kids_First_Biospecimen_ID)
    
    pop_subpop_df <- pop_subpop_df %>%
      dplyr::mutate(Kids_First_Biospecimen_ID = factor(Kids_First_Biospecimen_ID,
                                                       levels = order_pts)) %>%
      arrange(Kids_First_Biospecimen_ID) %>%
      dplyr::mutate(subpopulation = factor(subpopulation,
                                           levels = c("MXL", "CLM", "PEL", "PUR", "non-AMR")))
    
  }
  
  if (pop == "AFR"){
    
    order_pts <- pop_subpop_df %>%
      dplyr::filter(subpopulation %in% c("ASW")) %>%
      group_by(Kids_First_Biospecimen_ID) %>%
      summarize(probability = sum(probability)) %>%
      arrange(probability) %>%
      pull(Kids_First_Biospecimen_ID)
    
    pop_subpop_df <- pop_subpop_df %>%
      dplyr::mutate(Kids_First_Biospecimen_ID = factor(Kids_First_Biospecimen_ID,
                                                       levels = order_pts)) %>%
      arrange(Kids_First_Biospecimen_ID) %>%
      dplyr::mutate(subpopulation = factor(subpopulation,
                                           levels = c("ASW", "ACB", "ESN", 
                                                      "GWD", "LWK", "MSL", 
                                                      "YRI", "non-AFR")))
    
  }
  
  if (pop == "EAS"){
    
    order_pts <- pop_subpop_df %>%
      dplyr::filter(subpopulation %in% c("CDX")) %>%
      group_by(Kids_First_Biospecimen_ID) %>%
      summarize(probability = sum(probability)) %>%
      arrange(probability) %>%
      pull(Kids_First_Biospecimen_ID)
    
    pop_subpop_df <- pop_subpop_df %>%
      dplyr::mutate(Kids_First_Biospecimen_ID = factor(Kids_First_Biospecimen_ID,
                                                       levels = order_pts)) %>%
      arrange(Kids_First_Biospecimen_ID) %>%
      dplyr::mutate(subpopulation = factor(subpopulation,
                                           levels = c("CDX", "CHB", "CHS", 
                                                      "JPT", "KHV", "non-EAS")))
    
  }
  
  if (pop == "EUR"){
    
    order_pts <- pop_subpop_df %>%
      dplyr::filter(subpopulation %in% c("CEU")) %>%
      group_by(Kids_First_Biospecimen_ID) %>%
      summarize(probability = sum(probability)) %>%
      arrange(probability) %>%
      pull(Kids_First_Biospecimen_ID)
    
    pop_subpop_df <- pop_subpop_df %>%
      dplyr::mutate(Kids_First_Biospecimen_ID = factor(Kids_First_Biospecimen_ID,
                                                       levels = order_pts)) %>%
      arrange(Kids_First_Biospecimen_ID) %>%
      dplyr::mutate(subpopulation = factor(subpopulation,
                                           levels = c("CEU", "GBR", "FIN",
                                                      "IBS", "TSI", "non-EUR")))
    
  }
  
  if (pop == "SAS"){
    
    order_pts <- pop_subpop_df %>%
      dplyr::filter(subpopulation %in% c("BEB")) %>%
      group_by(Kids_First_Biospecimen_ID) %>%
      summarize(probability = sum(probability)) %>%
      arrange(probability) %>%
      pull(Kids_First_Biospecimen_ID)
    
    pop_subpop_df <- pop_subpop_df %>%
      dplyr::mutate(Kids_First_Biospecimen_ID = factor(Kids_First_Biospecimen_ID,
                                                       levels = order_pts)) %>%
      arrange(Kids_First_Biospecimen_ID) %>%
      dplyr::mutate(subpopulation = factor(subpopulation,
                                           levels = c("BEB", "GIH", "ITU", 
                                                      "PJL", "STU", "non-SAS")))
    
  }
  

  cbtn_subpop_plot <- pop_subpop_df %>%
    dplyr::filter(cohort == "CBTN") %>%
    ggplot(aes(x = Kids_First_Biospecimen_ID, y = probability,
               fill = subpopulation)) +
    geom_col(size = 0.05, width = 1,
             show.legend = FALSE) +
    scale_fill_manual(values = color_palette,
                      breaks = c(subpopulations, glue::glue("non-{pop}"))) +
    labs(x = "Patient", fill = "Subpopulation", title = "CBTN") +
    theme_Publication() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  pnoc_subpop_plot <- pop_subpop_df %>%
    dplyr::filter(cohort == "PNOC") %>%
    ggplot(aes(x = Kids_First_Biospecimen_ID, y = probability,
               fill = subpopulation)) +
    geom_col(size = 0.25, width = 1) +
    scale_fill_manual(values = color_palette,
                      breaks = c(subpopulations, glue::glue("non-{pop}"))) +
    labs(x = "Patient", fill = "Subpopulation", y = NULL, title = "PNOC") +
    theme_Publication() +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  subpop_plots <- ggarrange(cbtn_subpop_plot,
                            pnoc_subpop_plot,
                            widths = c(0.75, 0.25))

  ggsave(file.path(plots_dir,
                   glue::glue("{pop}_subpopulation_probablilties.pdf")),
           width = 15, height = 4)
  
}

# save subpopulation data with BS IDs 
subpop_df %>%
  dplyr::filter(Kids_First_Biospecimen_ID %in% ancestry$Kids_First_Biospecimen_ID) %>%
  left_join(ancestry %>% dplyr::select(Kids_First_Biospecimen_ID, predicted_ancestry)) %>%
  write_tsv(file.path(results_dir, "pbta-somalier-subpopulations.tsv"))
