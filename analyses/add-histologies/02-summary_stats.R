# R Corbett 2023
#
# Generate summary statistics for PBTA ancestry project

# load libraries and set directories
library(data.table)
library(tidyverse)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)
library(colorblindr)
library(ggpubr)
library(cowplot)

root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Set file paths
data_dir <- file.path(root_dir, "data")
analysis_dir <- file.path(root_dir, "analyses", "add-histologies")
results_dir <- file.path(analysis_dir, "results")
input_dir <- file.path(analysis_dir, "input")
plots_dir <- file.path(analysis_dir, "plots")

# call plot theme 
source(file.path(root_dir, "figures", "theme.R"))
source(file.path(analysis_dir, "util", "heatmap_function.R"))

# set file paths
ancestry_file <- file.path(results_dir, "merged_ancestry_histology_data.tsv")

# wrangle data
ancestry <- read_tsv(ancestry_file)

# Assess major ancestry probability across patients
ancestry <- ancestry %>%
  mutate(major_perc = case_when(
    predicted_ancestry == "EUR" ~ EUR_prob,
    predicted_ancestry == "AFR" ~ AFR_prob,
    predicted_ancestry == "AMR" ~ AMR_prob,
    predicted_ancestry == "EAS" ~ EAS_prob,
    predicted_ancestry == "SAS" ~ SAS_prob
  ))

# define colorblind-friendly palette
okabe_palette <- colorblindr::palette_OkabeIto[c(1:3,5:6)]
names(okabe_palette) <- c("AFR", "AMR", "EAS", "EUR", "SAS")

# Plot major predicted ancestry
pdf(file.path(plots_dir, "major_predicted_ancestry_hist.pdf"),
    width = 5, height = 3)

ggplot(data=ancestry, aes(x=major_perc,fill = predicted_ancestry)) + 
  geom_histogram(color = "black") + 
  theme_minimal() + 
  labs(x = "major ancestry probability", y = "No. Patients",
       fill = NULL) +
  scale_fill_manual(values = okabe_palette) + 
  geom_vline(xintercept = 0.9, linetype = "dashed")

dev.off()

# generate heatmap of ancestry probability among patients with major anc prob < 0.9
low_major <- ancestry %>%
  filter(major_perc < 0.9) %>%
  dplyr::select(Kids_First_Biospecimen_ID, 
         EAS_prob, AFR_prob, AMR_prob, SAS_prob, EUR_prob) %>%
  dplyr::rename("AFR" = "AFR_prob",
                "AMR" = "AMR_prob",
                "EAS" = "EAS_prob",
                "EUR" = "EUR_prob",
                "SAS" = "SAS_prob") %>%
  column_to_rownames("Kids_First_Biospecimen_ID")

colors <- colorRampPalette(brewer.pal(9, "Reds"))(100)

pdf(file.path(plots_dir, "low_major_ancestry_heatmap.pdf"),
    height = 6, width = 3)

pheatmap(low_major,
         name = "probability",
         col = colors,
         treeheight_row = 20,
         treeheight_col = 10,
         show_rownames =  F)

dev.off()

# Plot first four Somalier principal components 
pdf(file.path(plots_dir, "ancestry-pcs.pdf"),
    height = 4, width = 10)

pc12 <- ancestry %>%
  ggplot(aes(x = PC1, y = PC2, fill = predicted_ancestry)) +
  geom_point(size=2, shape=23,
             show.legend = FALSE) +
  scale_fill_manual(values = okabe_palette) +
  theme_Publication()


pc34 <- ancestry %>%
  ggplot(aes(x = PC3, y = PC4, fill = predicted_ancestry)) +
  geom_point(size=2, shape=23) +
  labs(fill = "predicted ancestry") +
  scale_fill_manual(values = okabe_palette,
                    labels=c(glue::glue("AFR (n={sum(ancestry$predicted_ancestry == 'AFR')})"),
                             glue::glue("AMR (n={sum(ancestry$predicted_ancestry == 'AMR')})"),
                             glue::glue("EAS (n={sum(ancestry$predicted_ancestry == 'EAS')})"),
                             glue::glue("EUR (n={sum(ancestry$predicted_ancestry == 'EUR')})"),
                             glue::glue("SAS (n={sum(ancestry$predicted_ancestry == 'SAS')})")
                    )) +
  theme_Publication()

ggarrange(pc12, pc34,
          widths = c(1.2,2))

dev.off()

# consolidate unreported/unknown reported races and ethnicities 
ancestry <- ancestry %>%
  dplyr::mutate(race = case_when(
    race %in% c("Not Reported", "Other", "Reported Unknown", "Not Reported/Unknown", "Unknown") | is.na(race) ~ "Race Unknown",
    race == "White/Caucasian" ~ "White",
    race == "American Indian or Alaska Native" ~ "AI/AN",
    race == "Native Hawaiian or Other Pacific Islander" ~ "NHPI",
    race == "Black or African American" ~ "Black/Afr. Am.",
    race == "More Than One Race" ~ ">1 Race",
    TRUE ~ race
  )) %>%
  mutate(ethnicity = case_when(
    grepl("Reported|Available|Unavailable", ethnicity) | is.na(ethnicity) ~ "Ethnicity Unknown",
    grepl("Non-Hispanic|Not Hispanic", ethnicity) ~ "Not Hispanic/Latino",
    TRUE ~ "Hispanic/Latino"
  ))

# calculate reported race sums across cohort 
race_total <- ancestry %>%
  count(race) %>%
  rename("total" = "n") %>%
  mutate(race = fct_relevel(race,
                            c("AI/AN", "Asian", "Black/Afr. Am.", "NHPI", 
                                  "White", ">1 Race", "Race Unknown"))) %>%
  arrange(desc(race))

# plot ancestry count by reported race
pdf(file.path(plots_dir, "predicted_ancestry_counts_by_race.pdf"),
    width = 7, height = 4)

ancestry %>%
  count(predicted_ancestry, race) %>%
  mutate(race = fct_relevel(race,
                       rev(c("AI/AN", "Asian", "Black/Afr. Am.", "NHPI", 
                                      "White", ">1 Race", "Race Unknown")))) %>%
  arrange(race) %>%
  ggplot(aes(x = race, y = n, fill = predicted_ancestry)) +
  geom_bar(position="stack", stat="identity", color = "black") +
  labs(y = "No. Patients", title = "reported race", 
       fill = "predicted ancestry", x = NULL) +   scale_x_discrete(labels=paste0(race_total$race, " (n=", race_total$total, ")")) +
  scale_fill_manual(values = okabe_palette) + 
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=12)) +
  coord_flip() +
  theme_Publication()

dev.off()

# plot ancestry percent by reported race
pdf(file.path(plots_dir, "predicted_ancestry_percent_by_race.pdf"),
    width = 7, height = 4)

ancestry %>%
  count(predicted_ancestry, race) %>%
  left_join(race_total, by = "race") %>%
  mutate(race = fct_relevel(race,
                       rev(c("AI/AN", "Asian", "Black/Afr. Am.", "NHPI", 
                                  "White", ">1 Race", "Race Unknown")))) %>%
  arrange(race) %>%
  mutate(perc = n/total) %>%
  ggplot(aes(x = race, y = perc, fill = predicted_ancestry)) +
  geom_bar(position="stack", stat="identity", color = "black") +
  labs(y = "fraction", title = "reported race", 
       fill = "predicted ancestry", x = NULL) + 
  scale_x_discrete(labels=paste0(race_total$race, " (n=", race_total$total, ")")) +
  scale_fill_manual(values = okabe_palette) + 
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=12)) +
  coord_flip() +
  theme_Publication()

dev.off()

# calculate reported ethnicity sums across cohort
ethnicity_total <- ancestry %>%
  count(ethnicity) %>%
  rename("total" = "n") %>%
  arrange(desc(ethnicity))

# plot ancestry count by reported ethnicity
pdf(file.path(plots_dir, "predicted_ancestry_counts_by_ethnicity.pdf"),
    width = 7, height = 4)

ancestry %>%
  count(predicted_ancestry, ethnicity) %>%
  mutate(race = factor(ethnicity,
                       rev(unique(ethnicity)[order(unique(ethnicity))]))) %>%
  ggplot(aes(x = ethnicity, y = n, fill = predicted_ancestry)) +
  geom_bar(position="stack", stat="identity", color = "black") +
  labs(y = "No. patients", title = "reported ethnicity", 
       fill = "predicted ancestry", x = NULL) +   scale_x_discrete(labels=paste0(ethnicity_total$ethnicity, " (n=", ethnicity_total$total, ")")) +
  scale_fill_manual(values = okabe_palette) + 
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=12)) +
  coord_flip() +
  theme_Publication()

dev.off()

# plot ancestry percent by reported ethnicity
pdf(file.path(plots_dir, "predicted_ancestry_percent_by_ethnicity.pdf"),
    width = 7, height = 4)

ancestry %>%
  count(predicted_ancestry, ethnicity) %>%
  left_join(ethnicity_total, by = "ethnicity") %>%
  mutate(ethnicity = factor(ethnicity,
                            rev(unique(ethnicity)[order(unique(ethnicity))]))) %>%
  mutate(perc = n/total) %>%
  ggplot(aes(x = ethnicity, y = perc, fill = predicted_ancestry)) +
  geom_bar(position="stack", stat="identity", color = "black") +
  labs(y = "fraction", title = "reported ethnicity", 
       fill = "predicted ancestry", x = NULL) + 
  scale_x_discrete(labels=paste0(ethnicity_total$ethnicity, " (n=", ethnicity_total$total, ")")) +
  scale_fill_manual(values = okabe_palette) + 
  theme(axis.text.x = element_text(angle = 60, hjust=1),
        plot.title = element_text(hjust = 0.5),
        text = element_text(size=12)) +
  coord_flip() +
  theme_Publication()

dev.off()

# plot race-ancestry enrichment
pdf(file.path(plots_dir, "race_ancestry_ct_enr_heatmap.pdf"),
    height = 3.5, width = 6)

race_ht <- plot_enr(ancestry, "race", "predicted_ancestry",
                 var1_names = c("AI/AN", "Asian", "Black/Afr. Am.",
                                "NHPI", "White", ">1 Race",
                                "Race Unknown"),
                 var2_names = c("AFR", "AMR", "EAS", "EUR", "SAS"),
                 padjust = TRUE)

draw(race_ht)

invisible(dev.off())

# plot ethnicity-ancestry enrichment
pdf(file.path(plots_dir, "ethnicity_ancestry_ct_enr_heatmap.pdf"),
    height = 2, width = 6)

ethn_ht <- plot_enr(ancestry, "ethnicity", "predicted_ancestry",
                    var1_names = c("Hispanic/Latino", "Not Hispanic/Latino",
                                   "Ethnicity Unknown"),
                    var2_names = c("AFR", "AMR", "EAS", "EUR", "SAS"),
                    padjust = TRUE)

draw(ethn_ht)

invisible(dev.off())

# plot enrichment results
pdf(file.path(plots_dir, "plot_group_ancestry_ct_enr_heatmap.pdf"),
    height = 8, width = 6)

# Plot plot_group - ancestry enrichment
group_ht <- plot_enr(ancestry[!is.na(ancestry$plot_group),], "plot_group", "predicted_ancestry",
                    var1_names = sort(unique(ancestry$plot_group[!is.na(ancestry$plot_group)])),
                    var2_names = c("AFR", "AMR", "EAS", "EUR", "SAS"),
                    padjust = TRUE)

draw(group_ht)

invisible(dev.off())

# plot extent of tumor resection heatmap in pLGG

ancestry <- ancestry %>%
  dplyr::mutate(extent_of_tumor_resection = case_when(
    extent_of_tumor_resection %in% c("Unavailable", "Not Reported") ~ "Unavailable/Unreported",
    grepl("Gross/Near total resection", extent_of_tumor_resection) ~ "Gross/Near total resection",
    TRUE ~ extent_of_tumor_resection
  )) %>%
  dplyr::mutate(extent_of_tumor_resection = fct_relevel(extent_of_tumor_resection,
                                                   c("Gross/Near total resection", "Partial resection", "Biopsy only",
                                                     "Unavailable/Unreported")))

group_df <- data.frame(plot_group = c("Atypical Teratoid Rhabdoid Tumor",
                                      "Craniopharyngioma",
                                      "Ependymoma", "Mixed neuronal-glial tumor",
                                      "Low-grade glioma", "Medulloblastoma",
                                      "Mesenchymal tumor", "Other high-grade glioma",
                                      "Schwannoma"),
                       abbreviation = c("atrt", "cpg", "epn", "gnt", "lgg",
                                        "mb", "mes", "hgg", "swn"))

pdf(NULL)

for (i in 1:nrow(group_df)){
  
  group <- group_df$plot_group[i]
  abbrev <- group_df$abbreviation[i]
  
  
  group_anc_resection <- ancestry %>%
    dplyr::filter(plot_group == group,
                  !is.na(extent_of_tumor_resection),
                  !grepl("Unavailable", extent_of_tumor_resection))

  pdf(file.path(plots_dir, glue::glue("{abbrev}_tumor_resection_by_predicted_ancestry.pdf")),
      width = 6, height = 3)

  resection_ht <- plot_enr(group_anc_resection, "extent_of_tumor_resection", "predicted_ancestry",
                          var1_names = sort(unique(group_anc_resection$extent_of_tumor_resection)),
                          var2_names = sort(unique(group_anc_resection$predicted_ancestry)),
                          padjust = TRUE)
  
  draw(resection_ht)
  
  invisible(dev.off())
  
  
  group_anc_region <- ancestry %>%
    dplyr::filter(plot_group == group,
                  !is.na(CNS_region),
                  !CNS_region %in% c("Mixed", "Other"))
  
  pdf(file.path(plots_dir, glue::glue("{abbrev}_CNS_region_by_predicted_ancestry.pdf")),
      width = 5, height = 3.5)
  
  region_ht <- plot_enr(group_anc_region, "CNS_region", "predicted_ancestry",
                           var1_names = sort(unique(group_anc_region$CNS_region)),
                           var2_names = sort(unique(group_anc_region$predicted_ancestry)),
                           padjust = TRUE)
  
  draw(region_ht)
  
  invisible(dev.off())
  
}

# plot LGG molecular subtype by predicted ancestry
pdf(file.path(plots_dir, "lgg_subtype_by_predicted_ancestry.pdf"),
    width = 5, height = 14)

lgg_subtype_ht <- plot_enr(ancestry[ancestry$plot_group == "Low-grade glioma" & !is.na(ancestry$molecular_subtype),], "molecular_subtype", "predicted_ancestry",
                          var1_names = sort(unique(ancestry$molecular_subtype[ancestry$plot_group == "Low-grade glioma" & !is.na(ancestry$molecular_subtype)])),
                          var2_names = c("AFR", "AMR", "EAS", "EUR", "SAS"),
                          padjust = FALSE)

draw(lgg_subtype_ht)

dev.off()

# Plot frequency of MB SHH subtypes among predicted ancestries
mb_shh_subtypes <- read_tsv(file.path(input_dir, "mb_shh_molecular_subtypes.tsv")) %>%
  dplyr::rename(molecular_subtype_mb_shh = molecular_subtype)

mb_shh <- ancestry %>%
  left_join(mb_shh_subtypes %>% dplyr::select(Kids_First_Biospecimen_ID, molecular_subtype_mb_shh),
            by = c("Kids_First_Biospecimen_ID_tumor" = "Kids_First_Biospecimen_ID")) %>%
  dplyr::filter(!is.na(molecular_subtype_mb_shh),
                molecular_subtype_mb_shh != "MB, SHH")

pdf(file.path(plots_dir, "mb_shh_subtype_by_predicted_ancestry.pdf"),
    width = 5, height = 3)

mbshh_subtype_ht <- plot_enr(mb_shh, "molecular_subtype_mb_shh", "predicted_ancestry",
                           var1_names = c("MB, SHH alpha", "MB, SHH beta", "MB, SHH gamma", "MB, SHH delta"),
                           var2_names = c("AFR", "AMR", "EAS", "EUR"),
                           padjust = TRUE)

draw(mbshh_subtype_ht)

dev.off()
