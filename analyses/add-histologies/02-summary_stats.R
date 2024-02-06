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
library(ggalluvial)
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
  column_to_rownames("Kids_First_Biospecimen_ID")

colors <- colorRampPalette(brewer.pal(9, "Reds"))(100)

pdf(file.path(plots_dir, "low_major_ancestry_heatmap.pdf"),
    height = 6, width = 3)

pheatmap(low_major,
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
  scale_fill_manual(values = okabe_palette,
                    labels=c("AFR (n=155)", "AMR (n=224)", "EAS (n=67)",
                             "EUR (n=996)", "SAS (n=43)")) +
  theme_Publication()


pc34 <- ancestry %>%
  ggplot(aes(x = PC3, y = PC4, fill = predicted_ancestry)) +
  geom_point(size=2, shape=23) +
  labs(fill = "predicted ancestry") +
  scale_fill_manual(values = okabe_palette,
                    labels=c("AFR (n=155)", "AMR (n=224)", "EAS (n=67)",
                             "EUR (n=996)", "SAS (n=43)")) +
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
  #filter(Kids_First_Biospecimen_ID %in% survival$Kids_First_Biospecimen_ID) %>%
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


race_enr <- matrix(0, length(unique(ancestry$predicted_ancestry)),
                   length(unique(ancestry$race[!is.na(ancestry$race)])),
                   dimnames = list(c("AFR", "AMR", "EAS", "EUR", "SAS"),
                                   c("AI/AN", "Asian", "Black/Afr. Am.", "NHPI", 
                                     "White", ">1 Race", "Race Unknown")))
race_pval <- race_enr

# loop through cancer groups to calculate enrichment 
for (i in 1:nrow(race_enr)){
  no_ancestry <- sum(grepl(rownames(race_enr)[i], ancestry$predicted_ancestry))
  for (j in 1:ncol(race_enr)){
    no_race <- sum(ancestry$race == colnames(race_enr)[j] & !is.na(ancestry$race))
    no_anc_race <- sum(ancestry$predicted_ancestry == rownames(race_enr)[i] & ancestry$race == colnames(race_enr)[j])
    race_enr[i,j] <- (no_anc_race/no_race)/(no_ancestry/nrow(ancestry))
    race_pval[i,j] <- phyper(no_anc_race, no_ancestry, nrow(ancestry) - no_ancestry, no_race, lower.tail = F)
  }
}

race_enr <- t(race_enr)
race_pval <- t(race_pval)

sig_mat <- ifelse(race_pval < 0.05 & race_enr > 1, "*", "")

fill_mat <- matrix(glue::glue("{round(race_enr, 2)}{sig_mat}"), 
                   nrow(race_enr), ncol(race_enr))

count_mat <- table(ancestry$race, ancestry$predicted_ancestry)
count_mat <- count_mat[rownames(race_enr), colnames(race_enr)]
race_ct_enr_mat <- matrix(glue::glue("{count_mat}\n({fill_mat})"),
                     nrow(count_mat), ncol(count_mat))

col_fun = colorRamp2(c(0, 20), c("white", "orangered"))

# plot enrichment results
pdf(file.path(plots_dir, "race_ancestry_ct_enr_heatmap.pdf"),
    height = 3.5, width = 6)

race_ht <- Heatmap(race_enr,
              name = "Odds ratio",
              cluster_rows = F,
              cluster_columns = F,
              rect_gp = gpar(col = "black", lwd = 2),
              col = col_fun,
              cell_fun = function(j, i, x, y, width, height, fill) {
                grid.text(sprintf("%s", race_ct_enr_mat[i, j]), x, y, gp = gpar(fontsize = 12))
              })

draw(race_ht)

invisible(dev.off())



ethn_enr <- matrix(0, length(unique(ancestry$predicted_ancestry)),
                   length(unique(ancestry$ethnicity[!is.na(ancestry$ethnicity)])),
                   dimnames = list(c("AFR", "AMR", "EAS", "EUR", "SAS"),
                                   c("Hispanic/Latino", "Not Hispanic/Latino",
                                     "Ethnicity Unknown")))
ethn_pval <- ethn_enr

# loop through cancer groups to calculate enrichment 
for (i in 1:nrow(ethn_enr)){
  no_ancestry <- sum(grepl(rownames(ethn_enr)[i], ancestry$predicted_ancestry))
  for (j in 1:ncol(ethn_enr)){
    no_ethn <- sum(ancestry$ethnicity == colnames(ethn_enr)[j] & !is.na(ancestry$ethnicity))
    no_anc_ethn <- sum(ancestry$predicted_ancestry == rownames(ethn_enr)[i] & ancestry$ethnicity == colnames(ethn_enr)[j])
    ethn_enr[i,j] <- (no_anc_ethn/no_ethn)/(no_ancestry/nrow(ancestry))
    ethn_pval[i,j] <- phyper(no_anc_ethn, no_ancestry, nrow(ancestry) - no_ancestry, no_ethn, lower.tail = F)
  }
}

ethn_enr <- t(ethn_enr)
ethn_pval <- t(ethn_pval)

sig_mat <- ifelse(ethn_pval < 0.05 & ethn_enr > 1, "*", "")

fill_mat <- matrix(glue::glue("{round(ethn_enr, 2)}{sig_mat}"), 
                   nrow(ethn_enr), ncol(ethn_enr))

count_mat <- table(ancestry$ethnicity, ancestry$predicted_ancestry)
count_mat <- count_mat[rownames(ethn_enr), colnames(ethn_enr)]
ethn_ct_enr_mat <- matrix(glue::glue("{count_mat}\n({fill_mat})"),
                          nrow(count_mat), ncol(count_mat))

col_fun = colorRamp2(c(0, 6), c("white", "orangered"))

# plot enrichment results
pdf(file.path(plots_dir, "ethnicity_ancestry_ct_enr_heatmap.pdf"),
    height = 2, width = 6)

ethnicity_ht <- Heatmap(ethn_enr,
                   name = "Odds ratio",
                   cluster_rows = F,
                   cluster_columns = F,
                   rect_gp = gpar(col = "black", lwd = 2),
                   col = col_fun,
                   cell_fun = function(j, i, x, y, width, height, fill) {
                     grid.text(sprintf("%s", ethn_ct_enr_mat[i, j]), x, y, gp = gpar(fontsize = 12))
                   })

draw(ethnicity_ht)

invisible(dev.off())




# create heatmap of plot group group count by predicted ancestry
pdf(file.path(plots_dir, "plot_group_by_ancestry.pdf"),
    height = 6, width = 5)

count_mat <- table(ancestry$plot_group, ancestry$predicted_ancestry)

col_fun = colorRamp2(c(0, 300), c("white", "orangered"))

Heatmap(count_mat,
        name = "Count",
        cluster_rows = F,
        cluster_columns = F,
        rect_gp = gpar(col = "black", lwd = 2),
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", count_mat[i, j]), x, y, gp = gpar(fontsize = 12))
        })

dev.off()


# calculate enrichment and pvalues of ancestries within plot groups
enr <- matrix(0, length(unique(ancestry$predicted_ancestry)),
              length(unique(ancestry$plot_group[!is.na(ancestry$plot_group)])),
              dimnames = list(c("AFR", "AMR", "EAS", "EUR", "SAS"),
                              c(unique(ancestry$plot_group[!is.na(ancestry$plot_group)]))))
pval <- enr

ancestry <- ancestry[!is.na(ancestry$plot_group),]

# loop through cancer groups to calculate enrichment 
for (i in 1:nrow(enr)){
  no_ancestry <- sum(grepl(rownames(enr)[i], ancestry$predicted_ancestry))
  for (j in 1:ncol(enr)){
    no_plotGroup <- sum(ancestry$plot_group == colnames(enr)[j] & !is.na(ancestry$plot_group))
    no_anc_plotGroup <- sum(ancestry$predicted_ancestry == rownames(enr)[i] & ancestry$plot_group == colnames(enr)[j])
    enr[i,j] <- (no_anc_plotGroup/no_plotGroup)/(no_ancestry/nrow(ancestry))
    pval[i,j] <- phyper(no_anc_plotGroup, no_ancestry, nrow(ancestry) - no_ancestry, no_plotGroup, lower.tail = F)
  }
}

# calcualte FDRs
fdr <- t(apply(pval, 1, function(x) p.adjust(x, "fdr")))

enr <- t(enr[,order(colnames(enr))])
fdr <- t(fdr[,order(colnames(fdr))])

sig_mat <- ifelse(fdr < 0.05 & enr > 1, "*", "")

fill_mat <- matrix(glue::glue("{round(enr, 1)}{sig_mat}"), 
                   nrow(enr), ncol(enr))

# Plot enrichment heatmap

col_fun = colorRamp2(c(0, 6), c("white", "orangered"))

# plot enrichment results
pdf(file.path(plots_dir, "plot_group_ancestry_enrichment_heatmap.pdf"),
    height = 6, width = 6)

Heatmap(enr,
         name = "Odds ratio",
         cluster_rows = F,
         cluster_columns = F,
         rect_gp = gpar(col = "black", lwd = 2),
         col = col_fun,
         cell_fun = function(j, i, x, y, width, height, fill) {
           grid.text(sprintf("%s", fill_mat[i, j]), x, y, gp = gpar(fontsize = 12))
         })

dev.off()

# create heatmap of plot group count and enrichment by ancestry
ct_enr_mat <- matrix(glue::glue("{count_mat}\n({fill_mat})"),
                     nrow(count_mat), ncol(count_mat))

pdf(file.path(plots_dir, "plot_group_ancestry_ct_enr_heatmap.pdf"),
    height = 8, width = 6)

Heatmap(enr,
        name = "Odds ratio",
        cluster_rows = F,
        cluster_columns = F,
        rect_gp = gpar(col = "black", lwd = 2),
        col = col_fun,
        cell_fun = function(j, i, x, y, width, height, fill) {
          grid.text(sprintf("%s", ct_enr_mat[i, j]), x, y, gp = gpar(fontsize = 12))
        })

dev.off()

# Plot cancer group distribution by ancestry, subsetting for unknown/unreported race patients
pdf(file.path(plots_dir, "plot_group_by_ancestry_unk_race.pdf"),
    height = 4, width = 4)

ancestry %>%
  filter(!is.na(plot_group) & race %in% c("Not Reported/Unknown", "More Than One Race")) %>%
  mutate(predicted_ancestry = factor(predicted_ancestry, levels = unique(predicted_ancestry)),
         plot_group = factor(plot_group, rev(unique(plot_group)[order(unique(plot_group))]))) %>%
  count(predicted_ancestry, plot_group, name = "count", .drop = FALSE) %>%
  mutate(plot_group = factor(plot_group)) %>%
  ggplot(aes(x = predicted_ancestry, y = factor(plot_group), fill = count)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1, show.legend = FALSE) +
  scale_fill_gradient(low="white", high="orangered") +
  geom_text(aes(label = count), color = "black", size = 3) +
  labs(x = NULL, y = NULL) +
  theme_minimal()

dev.off()

# Get tumor resection type sums for LGG cohort
lgg_anc_total <- ancestry %>%
  filter(plot_group == "Low-grade glioma") %>%
  count(predicted_ancestry) %>%
  rename("total" = "n") %>%
  arrange(desc(predicted_ancestry))

# plot extent of tumor resection heatmap by ancestry (LGG only)
pdf(file.path(plots_dir, "lgg_tumor_resection_by_predicted_ancestry.pdf"),
    width = 5, height = 3)

ancestry %>%
  filter(plot_group == "Low-grade glioma") %>%
  dplyr::mutate(extent_of_tumor_resection = case_when(
    extent_of_tumor_resection %in% c("Unavailable", "Not Applicable", "Not Reported") ~ "Unavailable/Unreported/Not Applicable",
    grepl("Gross/Near total resection", extent_of_tumor_resection) ~ "Gross/Near total resection",
    TRUE ~ extent_of_tumor_resection
  )) %>%
  dplyr::mutate(extent_of_tumor_resection = factor(extent_of_tumor_resection,
                                                   rev(c("Gross/Near total resection", "Partial resection",
                                                         "Biopsy only", "Unavailable/Unreported/Not Applicable")
                                                   ))) %>%
  count(predicted_ancestry, extent_of_tumor_resection, .drop = FALSE) %>%
  left_join(lgg_anc_total, by = "predicted_ancestry") %>%
  mutate(perc = round(n/total*100, 2)) %>%
  ggplot(aes(x = predicted_ancestry, y = factor(extent_of_tumor_resection), fill = perc)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1, show.legend = FALSE) +
  scale_fill_gradient(low="white", high="orangered") +
  geom_text(aes(label = perc), color = "black", size = 3) +
  labs(x = NULL, y = NULL) +
  theme_minimal()

dev.off()

# Plot tumor location by predicted ancestry (LGG only)
pdf(file.path(plots_dir, "lgg_tumor_location_by_predicted_ancestry.pdf"),
    width = 3.5, height = 3)

ancestry %>%
  filter(plot_group == "Low-grade glioma") %>%
  dplyr::mutate(CNS_region = factor(CNS_region, levels = unique(CNS_region)[rev(order(unique(CNS_region)))])) %>%
  count(predicted_ancestry, CNS_region, .drop = FALSE) %>%
  left_join(lgg_anc_total, by = "predicted_ancestry") %>%
  mutate(perc = round(n/total*100, 2)) %>%
  ggplot(aes(x = predicted_ancestry, y = factor(CNS_region), fill = perc)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1, show.legend = FALSE) +
  scale_fill_gradient(low="white", high="orangered") +
  geom_text(aes(label = perc), color = "black", size = 3) +
  labs(x = NULL, y = "CNS region") +
  theme_minimal()

dev.off()


# plot LGG molecular subtype by predicted ancestry
pdf(file.path(plots_dir, "lgg_subtype_by_predicted_ancestry.pdf"),
    width = 5, height = 6)

ancestry %>%
  filter(plot_group == "Low-grade glioma") %>%
  dplyr::mutate(molecular_subtype = factor(molecular_subtype,
                                           levels = rev(unique(molecular_subtype[order(molecular_subtype)])))) %>%
  count(predicted_ancestry, molecular_subtype, .drop = FALSE) %>%
  left_join(lgg_anc_total, by = "predicted_ancestry") %>%
  mutate(perc = round(n/total*100, 2)) %>%
  ggplot(aes(x = predicted_ancestry, y = molecular_subtype, fill = n)) +
  geom_tile(color = "black",
            lwd = 1,
            linetype = 1, show.legend = FALSE) +
  scale_fill_gradient(low="white", high="orangered") +
  geom_text(aes(label = n), color = "black", size = 3) +
  labs(x = NULL, y = NULL) +
  theme_minimal()

dev.off()

