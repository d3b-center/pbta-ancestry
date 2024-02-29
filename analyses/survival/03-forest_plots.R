# S. Spielman for ALSF CCDL & Jo Lynne Rokita for D3b, 2022
#
# Makes a pdf panel of forest plot of survival analysis on HGG samples 
#  with molecular subtype as predictors


library(survival) # needed to parse model output
library(tidyverse)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))

# Declare output directory
plots_dir <- file.path(root_dir, "analyses", "survival", "plots")
if (!dir.exists(plots_dir)) {
  dir.create(plots_dir, recursive = TRUE)
}

# Input directory
input_dir <- file.path(root_dir, "analyses", "survival", "results")

lgg_forest_pdf <- file.path(plots_dir, "forest_lgg_ancestry.pdf")
hgg_forest_pdf <- file.path(plots_dir, "forest_hgg_ancestry.pdf")
hgg_forest_pc_pdf <- file.path(plots_dir, "forest_hgg_pc.pdf")
hgg_forest_eur_pdf <- file.path(plots_dir, "forest_hgg_eur.pdf")


hgg_wt_forest_pdf <- file.path(plots_dir, "forest_hgg_wt_ancestry.pdf")
hgg_wt_forest_pc_pdf <- file.path(plots_dir, "forest_hgg_wt_pc.pdf")
hgg_wt_forest_eur_pdf <- file.path(plots_dir, "forest_hgg_wt_eur.pdf")


hgg_k28_forest_pdf <- file.path(plots_dir, "forest_hgg_k28_ancestry.pdf")
hgg_k28_forest_pc_pdf <- file.path(plots_dir, "forest_hgg_k28_pc.pdf")
hgg_k28_forest_eur_pdf <- file.path(plots_dir, "forest_hgg_k28_eur.pdf")


epn_forest_pdf <- file.path(plots_dir, "forest_epn_ancestry.pdf")
mb_forest_pdf <- file.path(plots_dir, "forest_mb_ancestry.pdf")
mb4_forest_pdf <- file.path(plots_dir, "forest_mb4_ancestry.pdf")

atrt_forest_pdf <- file.path(plots_dir, "forest_atrt_ancestry.pdf")
atrt_forest_pc_pdf <- file.path(plots_dir, "forest_atrt_pc.pdf")

## Input and output files -------------------------------------
survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_lgg_additive_terms_resection_subtype_ancestry.RDS"
  )
)


## Make plot --------------------------------------------------
ref_term_anc <- "predicted_ancestryEUR"
ref_term_res <- "extent_of_tumor_resectionBiopsy only"
ref_term_type <- "mol_sub_groupWT"


# Set up ordering and labels for y-axis
term_order <- c("predicted_ancestrySAS",
                    "predicted_ancestryEAS",
                    "predicted_ancestryAMR",
                    "predicted_ancestryAFR",
                    ref_term_anc,
                    "mol_sub_groupOther alteration",
                    "mol_sub_groupBRAF fusion",
                    ref_term_type,
                    "extent_of_tumor_resectionPartial resection",
                    "extent_of_tumor_resectionGross/Near total resection",
                    ref_term_res)

term_labels <- c("SAS",
                     "EAS",
                     "AMR",
                     "AFR",
                     "EUR",
                     "Other alteration",
                     "BRAF fusion",
                     "WT",
                     "Partial resection",
                     "Gross/near total resection",
                     "Biopsy only")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add reference
  add_row(term = ref_term_anc, 
          estimate = 0) %>% 
  add_row(term = ref_term_type, 
          estimate = 0) %>%
  add_row(term = ref_term_res, 
          estimate = 0) %>%
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )


# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "EUR (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0, size = 3) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)



# Export plot
ggsave(lgg_forest_pdf, forest_panels, width = 10, height = 3)







survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_hgg_additive_terms_subtype_ancestry.RDS"
  )
)


ref_term_anc <- "predicted_ancestryEUR"
ref_term_type <- "mol_sub_groupH3 WT"


# Set up ordering and labels for y-axis
term_order <- c("predicted_ancestrySAS", "predicted_ancestryEAS",
                    "predicted_ancestryAMR", "predicted_ancestryAFR",
                    ref_term_anc,
                    "mol_sub_groupIDH", "mol_sub_groupH3 K28",
                    "mol_sub_groupH3 G35", ref_term_type)

term_labels <- c("SAS", "EAS", "AMR", "AFR", "EUR",
                     "IDH", "H3 K28", "H3 G35", "WT")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add reference
  add_row(term = ref_term_anc, 
          estimate = 0) %>% 
  add_row(term = ref_term_type, 
          estimate = 0) %>%
  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )

# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "EUR (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)

# Export plot
ggsave(hgg_forest_pdf, forest_panels, width = 10, height = 3)




survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_hgg_additive_terms_subtype_pc.RDS"
  )
)


ref_term_type <- "mol_sub_groupH3 WT"


# Set up ordering and labels for y-axis
term_order <- c("PC5",
                "mol_sub_groupIDH", "mol_sub_groupH3 K28",
                "mol_sub_groupH3 G35", ref_term_type)

term_labels <- c("PC5",
                 "IDH", "H3 K28", "H3 G35", "WT")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add reference

  add_row(term = ref_term_type, 
          estimate = 0) %>%
  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )

# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "WT", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)

# Export plot
ggsave(hgg_forest_pc_pdf, forest_panels, width = 10, height = 3)



survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_hgg_additive_terms_subtype_EURstatus.RDS"
  )
)

ref_term_ancestry <- "EUR_statusEUR"
ref_term_type <- "mol_sub_groupH3 WT"


# Set up ordering and labels for y-axis
term_order <- c("EUR_statusnon-EUR", ref_term_ancestry,
                "mol_sub_groupIDH", "mol_sub_groupH3 K28",
                "mol_sub_groupH3 G35", ref_term_type)

term_labels <- c("non-EUR", "EUR",
                 "IDH", "H3 K28", "H3 G35", "WT")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add reference
  add_row(term = ref_term_ancestry, 
          estimate = 0) %>%
  add_row(term = ref_term_type, 
          estimate = 0) %>%
  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )

# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "WT" | term == "EUR", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)

# Export plot
ggsave(hgg_forest_eur_pdf, forest_panels, width = 10, height = 3)




survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_hgg_wt_additive_terms_ancestry.RDS"
  )
)


ref_term_anc <- "predicted_ancestryEUR"

# Set up ordering and labels for y-axis
term_order <- c("predicted_ancestryEAS",
                "predicted_ancestryAMR",
                "predicted_ancestryAFR",
                ref_term_anc)

term_labels <- c("EAS",
                 "AMR",
                 "AFR",
                 "EUR")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add reference
  add_row(term = ref_term_anc, 
          estimate = 0) %>% 
  filter(!is.na(estimate)) %>%

  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )


# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "EUR (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)



# Export plot
ggsave(hgg_wt_forest_pdf, forest_panels, width = 10, height = 3)





survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_hgg_wt_additive_terms_pc.RDS"
  )
)



# Set up ordering and labels for y-axis
term_order <- c("PC5",
                "PC4",
                "PC3")

term_labels <- c("PC5",
                 "PC4",
                 "PC3")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%

  filter(!is.na(estimate)) %>%
  
  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )


# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "EUR (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)



# Export plot
ggsave(hgg_wt_forest_pc_pdf, forest_panels, width = 10, height = 3)




survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_hgg_wt_additive_terms_EURstatus.RDS"
  )
)

ref_term_anc = "EUR_statusEUR"

# Set up ordering and labels for y-axis
term_order <- c("EUR_statusnon-EUR",
                "EUR_statusEUR")

term_labels <- c("non-EUR",
                 "EUR")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  add_row(term = ref_term_anc, 
          estimate = 0) %>% 
  filter(!is.na(estimate)) %>%
  
  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )


# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "EUR (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)



# Export plot
ggsave(hgg_wt_forest_eur_pdf, forest_panels, width = 10, height = 3)




survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_hgg_k28_additive_terms_ancestry.RDS"
  )
)


ref_term_anc <- "predicted_ancestryEUR"


# Set up ordering and labels for y-axis
term_order <- c("predicted_ancestrySAS",
                "predicted_ancestryEAS",
                "predicted_ancestryAMR",
                "predicted_ancestryAFR",
                ref_term_anc)

term_labels <- c("SAS",
                 "EAS",
                 "AMR",
                 "AFR",
                 "EUR")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add reference
  add_row(term = ref_term_anc, 
          estimate = 0) %>% 
  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )


# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "EUR (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)



# Export plot
ggsave(hgg_k28_forest_pdf, forest_panels, width = 10, height = 3)



survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_hgg_k28_additive_terms_pc.RDS"
  )
)


# Set up ordering and labels for y-axis
term_order <- c("PC4",
                "PC3",
                "PC2")

term_labels <- c("PC4",
                 "PC3",
                 "PC2")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add reference
  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )


# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "EUR (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)



# Export plot
ggsave(hgg_k28_forest_pc_pdf, forest_panels, width = 10, height = 3)





survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_hgg_k28_additive_terms_EURstatus.RDS"
  )
)


ref_term_anc <- "EUR_statusEUR"


# Set up ordering and labels for y-axis
term_order <- c("EUR_statusnon-EUR",
                ref_term_anc)

term_labels <- c("non-EUR",
                 "EUR")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add reference
  add_row(term = ref_term_anc, 
          estimate = 0) %>% 
  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )


# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "EUR (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)



# Export plot
ggsave(hgg_k28_forest_eur_pdf, forest_panels, width = 10, height = 3)







survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_mb_additive_terms_subtype_ancestry.RDS"
  )
)


ref_term_anc <- "predicted_ancestryEUR"
ref_term_type <- "OPC_methyl_subtypeMB, Group3"


# Set up ordering and labels for y-axis
term_order <- c("predicted_ancestrySAS",
 #               "predicted_ancestryEAS",
                "predicted_ancestryAMR",
                "predicted_ancestryAFR",
                ref_term_anc,
 #               "OPC_methyl_subtypeMB, WNT",
                "OPC_methyl_subtypeMB, SHH",
                "OPC_methyl_subtypeMB, Group4",
                "OPC_methyl_subtypeMB, Group3")

term_labels <- c("SAS",
 #                "EAS",
                 "AMR",
                 "AFR",
                 "EUR",
 #                "MB, WNT",
                 "MB, SHH",
                 "MB, Group4",
                 "MB, Group3")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  dplyr::select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add reference
  add_row(term = ref_term_anc, 
          estimate = 0) %>% 
  add_row(term = ref_term_type, 
          estimate = 0) %>% 
  filter(!is.na(estimate)) %>%
  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )


# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  dplyr::select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "EUR (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)



# Export plot
ggsave(mb_forest_pdf, forest_panels, width = 10, height = 3)





survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_mb4_additive_terms_ancestry.RDS"
  )
)


ref_term_anc <- "predicted_ancestryEUR"


# Set up ordering and labels for y-axis
term_order <- c("predicted_ancestrySAS",
                "predicted_ancestryEAS",
                "predicted_ancestryAMR",
                "predicted_ancestryAFR",
                ref_term_anc)

term_labels <- c("SAS",
                 "EAS",
                 "AMR",
                 "AFR",
                 "EUR")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  dplyr::select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add reference
  add_row(term = ref_term_anc, 
          estimate = 0) %>% 
  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )


# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  dplyr::select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "EUR (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)



# Export plot
ggsave(mb4_forest_pdf, forest_panels, width = 10, height = 3)





survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_atrt_additive_terms_subtype_ancestry.RDS"
  )
)


ref_term_anc <- "predicted_ancestryAFR"
#ref_term_type <- "dkfz_v12_methylation_subclassATRT_MYC"



# Set up ordering and labels for y-axis
term_order <- c("predicted_ancestryEUR",
                "predicted_ancestryAMR",
                ref_term_anc)
#                "dkfz_v12_methylation_subclassATRT_TYR",
#                "dkfz_v12_methylation_subclassATRT_SHH",
#                ref_term_type)

term_labels <- c("EUR",
                 "AMR",
                 "AFR")
 #                "ATRT_TYR",
#                 "ATRT_SHH",
#                 "ATRT_MYC")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  dplyr::select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add reference
  add_row(term = ref_term_anc, 
          estimate = 0) %>% 
#  add_row(term = ref_term_type, 
#          estimate = 0) %>% 
  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )


# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  dplyr::select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "EUR (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)



# Export plot
ggsave(atrt_forest_pdf, forest_panels, width = 10, height = 3)











survival_result <- read_rds(
  file.path(
    input_dir,
    "cox_atrt_additive_terms_subtype_pc.RDS"
  )
)


ref_term_type <- "dkfz_v12_methylation_subclassATRT_MYC"



# Set up ordering and labels for y-axis
term_order <- c("PC2",
                "PC1",
                "dkfz_v12_methylation_subclassATRT_TYR",
                "dkfz_v12_methylation_subclassATRT_SHH",
                ref_term_type)

term_labels <- c("PC1",
                 "PC2",
                 "ATRT_TYR",
                 "ATRT_SHH",
                 "ATRT_MYC")


# Get n and event info from glance output
survival_n <- broom::glance(survival_result) %>%
  select(n, nevent)

# Convert survival model result to data frame, and exponentiate estimates/CIs to get HRs
survival_df <- broom::tidy(survival_result) %>%
  # Add reference
  filter(term != 'cohortCBTN2') %>%
  add_row(term = ref_term_type, 
          estimate = 0) %>% 
  
  mutate(
    conf.low = exp(estimate-std.error),
    conf.high = exp(estimate+std.error),
    estimate = exp(estimate),
    # significance indicator column for filling points.
    # Note T/F these are strings for type compatibility with "REF"
    significant = case_when(p.value <= 0.05 ~ "TRUE", 
                            p.value > 0.05 ~ "FALSE", 
                            is.na(p.value) ~ "REF"),
    # y-axis factor re-labeling
    term = factor(term, 
                  levels = term_order,
                  labels = term_labels)
  )


# Forest plot of the model
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 1 rows containing missing values (geom_errorbarh).
forest_plot <- ggplot(survival_df) +
  aes(x = estimate, 
      y = term,
      fill = significant
  ) + 
  # add CI first so line doesn't cover open point
  geom_errorbarh(
    aes(
      xmin = conf.low,
      xmax = conf.high,
    ),
    height = 0.15,
    size = 0.65
  ) + 
  geom_point(
    size = 3.5,
    shape = 23
  ) +
  # Point fill based on sigificance
  scale_fill_manual(
    values = c("FALSE" = "white", 
               "TRUE" = "black",
               "REF" = "gray"),
    guide = FALSE # turn off legend
  ) + 
  # Vertical guiding line at 1
  geom_vline(
    xintercept = 1, 
    linetype = 3
  ) +
  labs(
    x = "Hazard Ratio ± 95% CI",
    y = "",
    subtitle = glue::glue("N = {survival_n$n} with {survival_n$nevent} events")
  ) + 
  # log-scale the x-axis
  scale_x_log10() +
  ggpubr::theme_pubr() + 
  theme(
    plot.subtitle = element_text(face = "bold")
  ) +
  # grid makes it easier to follow lines
  cowplot::background_grid()


# Accompanying panel with sample sizes, P-values, etc.

# prepare data for panel
# note this warning is OK and EXPECTED because there is no CI for the reference group: 
#    Removed 2 rows containing missing values (geom_text). 
survival_df_spread <- survival_df %>%
  mutate(
    # Clean pvalues into labels. 
    p_string = if_else(
      p.value >= 0.001, 
      paste0("P = ",round(p.value,3)),
      "P < 0.001"
    ),
    # round to 2 digits and create single string with "hr (low-high)"
    conf.low = signif(conf.low, 2),
    conf.high = signif(conf.high, 2),
    estimate = signif(estimate, 2),
    hr_ci = glue::glue("{estimate} ({conf.low} - {conf.high})")
  ) %>%
  select(term, hr_ci, p_string) %>%
  # this throws a warning but it's ok
  # format tibble for plotting
  gather(hr_ci:p_string, key = "name", value = "value") %>%
  #remove values for reference
  mutate(value = ifelse(term == "EUR (ref)", NA, value))

labels_panel <- ggplot(survival_df_spread) +
  aes(x = name, y = term, label = value) + 
  geom_text(hjust = 0) +
  labs(
    # hack!
    subtitle = paste0("                    ",
                      "HR (95% CI)             P-value")
  ) +
  ggpubr::theme_pubr() + 
  # remove axes.
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    # -26 is as low as we can go before plot starts to get coverd
    plot.margin = margin(6, 0, 36, -25, unit = "pt"),
    plot.subtitle = element_text(face = "bold")
  ) 

forest_panels <- cowplot::plot_grid(forest_plot, labels_panel, nrow = 1, rel_widths = c(1,0.5), scale = 0.95)



# Export plot
ggsave(atrt_forest_pc_pdf, forest_panels, width = 10, height = 3)

