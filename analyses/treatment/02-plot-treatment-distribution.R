
res <- read_tsv("/Users/corbettr/Documents/dei_ancestry/treatment_frequencies.tsv") %>%
  mutate(group = factor(group, c("non-EUR", "EUR"))) %>%
  mutate(fraction = glue::glue("{n}/{total_n}")) %>%
  mutate(p_text = case_when(
    !is.na(pvalue) ~ glue::glue("p = {round(pvalue, 3)}"),
    TRUE ~ ""
  ))


# Create CPG Odds Ratio plot 

enr_plot <- res %>% 
  ggplot(aes(x = factor(group), y = OddsRatio, label = p_text)) +
  geom_point(size = 3, color = "#00A087FF",
             show.legend = FALSE) + 
  geom_errorbar(aes(ymin = lower_CI, ymax = upper_CI), width = 0.2, 
                show.legend = FALSE, color = "#00A087FF") +
  labs(y = "Odds Ratio (95% CI)", x = NULL) + 
  geom_text(y = 0, x = 2, hjust = 0, size = 4, fontface = 2) +
  coord_flip() +
  geom_hline(yintercept = 1, linetype = "dashed") + 
  facet_wrap(~category, nrow = 3, scale = "fixed") +
  expand_limits(y=0) +
  theme_Publication() +
  theme(plot.margin = unit(c(2,0.5,1,0.25), "lines"))



# Create % patients with CPG PLP plot separated by source call, and include fractions as text 

perc_plot <- res %>%
  ggplot(aes(x = perc, y = factor(group), label = fraction, fill = group)) +
  geom_bar(stat = "identity", color = "black",
           show.legend = FALSE) + 
  geom_text(x = 67, hjust = 0, size = 4, fontface = 2) +
  labs(x = "% Patients", y = NULL, fill = NULL) + 
  guides(fill = guide_legend(nrow = 1)) +
  facet_wrap(~category, nrow = 3, scale = "fixed") +
  expand_limits(x=2) +
  coord_cartesian(clip = 'off') +
  theme_Publication() +
  theme(plot.margin = unit(c(2,4,1,1), "lines"),
        legend.position = c(0.5, 1.07))




tiff(file.path("/Users/corbettr/Documents/dei_ancestry/treatment-enrichment-barplot.tiff"),
     width = 7, height = 5, units = "in", res = 300)

ggarrange(enr_plot, perc_plot,
          nrow = 1, widths = c(1.5,2))

dev.off()
