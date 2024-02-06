# R Corbett 2023
#
# Generate alluvial plot for PBTA ancestry project

# load libraries and set directories
library(data.table)
library(tidyverse)
# library(ComplexHeatmap)
# library(RColorBrewer)
# library(circlize)
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
    grepl("Reported|Available|Unavailable", ethnicity) | is.na(ethnicity) ~ "Unknown",
    grepl("Non-Hispanic|Not Hispanic", ethnicity) ~ "Not Hispanic/Latino",
    TRUE ~ "Hispanic/Latino"
  ))


# Create data frame for alluvial plots, including only predicted_ancestry, race, and ethnicity columns
alluvial_df <- as.data.frame(table(ancestry$predicted_ancestry, ancestry$race, ancestry$ethnicity)) %>% 
  # format for alluvial plot
  dplyr::rename(predicted_ancestry = Var1, 
                race = Var2,
                ethnicity = Var3) %>%
  to_lodes_form(axes = 1:3) %>% 
  dplyr::rename(Group = stratum) %>% 
  mutate(Group = factor(Group, levels = c("EAS", "SAS", "AFR", "AMR", "EUR",
                                          "Asian", "Black/Afr. Am.", "AI/AN", 
                                          "NHPI", "White", ">1 Race",
                                          "Race Unknown", 
                                          "Hispanic/Latino",
                                          "Not Hispanic/Latino",
                                          "Unknown"))) %>% 
  mutate(race = case_when(Group %in% c("Asian", "Black/Afr. Am.", "AI/AN", 
                                       "NHPI", "White", ">1 Race",
                                       "Race Unknown") ~ Group, 
                          TRUE ~ NA), 
         predicted_ancestry = case_when(Group %in% c("EAS", "SAS", "AFR", "EUR", "AMR") ~ Group, 
                                        TRUE ~ NA),
         ethnicity = case_when(Group %in% c("Hispanic/Latino",
                                            "Not Hispanic/Latino",
                                            "Unknown") ~ Group,
                               TRUE ~ NA))

# Generate alluvial plot
p1 <- ggplot(alluvial_df, aes(y = Freq, stratum = Group, alluvium = alluvium, x = x, fill = Group)) + 
  geom_alluvium(show.legend = F) + 
  geom_stratum(show.legend = F) +
  scale_fill_manual(values = c("Asian" = "#009E73", "SAS" = "#D55E00", "EAS" = "#009E73",
                               "White" = "#0072B2", "EUR" = "#0072B2",
                               "Black/Afr. Am." = "#E69F00", "AFR" = "#E69F00",
                               "NHPI" = "skyblue", "AMR" = "#56B4E9",
                               "AI/AN" = "#56B4E9",
                               "Race Unknown" = "grey", ">1 Race" = "brown",
                               "Hispanic/Latino" = "#56B4E9",
                               "Not Hispanic/Latino" = "#0072B2",
                               "Unknown" = "grey")) +
  xlab("") + 
  ylab("Number of Patients") +
  scale_x_discrete(labels = c("predicted ancestry", "reported race", "reported ethnicity")) + 
  theme_Publication()

# Create separate ancestry, race, and ethnicity dfs for legend generation
race_df <- data.frame(race = c("Asian", "Black/Afr. Am.", "AI/AN", 
                               "NHPI", "White", ">1 Race",
                               "Race Unknown"), 
                      value = 1)
ancestry_df <- data.frame(predicted_ancestry = c("EAS", "SAS", "AFR", "EUR", "AMR"), 
                          value = 1)
ethnicity_df <- data.frame(ethnicity = c("Hispanic/Latino",
                                         "Not Hispanic/Latino",
                                         "Unknown"), 
                           value = 1)

# plot race legend
lgd_race <- ggplot(race_df, aes(x = value, y = factor(race, levels = c("Race Unknown", 
                                                                       ">1 Race", "White", 
                                                                       "AI/AN",
                                                                       "NHPI", 
                                                                       "Black/Afr. Am.", "Asian")), 
                                fill = race)) + 
  geom_tile(show.legend = T, color = "black",
            lwd = 0.5, linetype = 1) + 
  scale_fill_manual(values = c("Asian" = "#009E73",  
                               "White" = "#0072B2", 
                               "Black/Afr. Am." = "#E69F00", 
                               "AI/AN" = "#56B4E9",
                               "NHPI" = "skyblue", 
                               ">1 Race" = "brown",
                               "Race Unknown" = "grey"),
                    breaks = c("Asian",  
                               "Black/Afr. Am.", 
                               "AI/AN",
                               "NHPI", 
                               "White", 
                               ">1 Race",
                               "Race Unknown")) +
  labs(fill = "Reported Race") +
  theme_Publication()
leg_race <- get_legend(lgd_race)

# Plot ancestry legend
lgd_anc <- ggplot(ancestry_df, aes(x = value, y = factor(predicted_ancestry, 
                                                         levels = c("SAS","EAS","EUR","AMR","AFR")), 
                                   fill = predicted_ancestry)) + 
  geom_tile(show.legend = T, col = "black",
            lwd = 0.5, linetype = 1) + 
  #xlim() + 
  scale_fill_manual(values = c("SAS" = "#D55E00", 
                               "EAS" = "#009E73",
                               "EUR" = "#0072B2", 
                               "AFR" = "#E69F00", 
                               "AMR" = "#56B4E9"), breaks = c("EAS", "SAS", "AFR", "AMR", "EUR")) +
  labs(fill = "Predicted Ancestry") +
  theme_Publication()
leg_anc <- get_legend(lgd_anc)

# Plot ethnicity legend
lgd_ethn <- ggplot(ethnicity_df, aes(x = value, y = factor(ethnicity, 
                                                           levels = c("Hispanic/Latino",
                                                                      "Not Hispanic/Latino",
                                                                      "Unknown")), 
                                     fill = ethnicity)) + 
  geom_tile(show.legend = T, col = "black",
            lwd = 0.5, linetype = 1) + 
  #xlim() + 
  scale_fill_manual(values = c("Hispanic/Latino" = "#56B4E9",
                               "Not Hispanic/Latino" = "#0072B2",
                               "Unknown" = "grey")) +
  labs(fill = "Reported Ethnicity") +
  theme_Publication()
leg_ethn <- get_legend(lgd_ethn)

leg <- plot_grid(leg_anc, leg_race, leg_ethn,
                 nrow = 3,
                 align = "v",
                 rel_heights = c(1.5, 0.4, 1.5))


final_p <- plot_grid(p1, leg,
                     nrow = 1,
                     align = "none",
                     axis = "t",
                     rel_widths = c(1,0.65))

pdf(file.path(plots_dir, "ancestry-race-ethnicity-alluvial.pdf"),
    width = 10, height = 6)

final_p

dev.off()