# R Corbett 2023
#
# Generate Kaplan-Meier survival curves for DEI ancestry project

library(tidyverse)
library(survival)
library(patchwork)
library(colorblindr)

# Establish base dir
root_dir <- rprojroot::find_root(rprojroot::has_dir(".git"))
analysis_dir <- file.path(root_dir, "analyses", "survival")

source(file.path(analysis_dir, "util", "survival_models.R"))

# Declare output directory
plots_dir <- file.path(analysis_dir, "plots")

# Input directory
input_dir <- file.path(analysis_dir, "results")


groups <- c("DIPG or DMG|Other high-grade glioma",
            "Mixed neuronal-glial tumor|Low-grade glioma"
)


file_names <- c("hgg", "lgg")

names(file_names) <- groups

titles <- c("High-grade glioma", "Low-grade glioma")
names(titles) <- groups


for (group in groups){
  
  km_os_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{file_names[group]}_OS_predicted_ancestry_raceUK.RDS")
    )
  )
  
  km_efs_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{file_names[group]}_EFS_predicted_ancestry_raceUK.RDS")
    )
  )
  
  km_output_pdf <- file.path(plots_dir, glue::glue("km_{file_names[group]}_OS_EFS_predicted_ancestry_raceUK.pdf"))
  
  km_plot <- plotKM(model = list(km_os_result, km_efs_result),
                    variable = "predicted_ancestry",
                    combined = T, 
                    title = titles[group])
  
  ggsave(km_output_pdf, km_plot,
         width = 10, height = 7, units = "in",
         device = "pdf")
  
  
}



for (group in groups[2]){
  
  km_os_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{file_names[group]}_OS_EUR_nonEUR_raceUK.RDS")
    )
  )
  
  km_efs_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{file_names[group]}_EFS_EUR_nonEUR_raceUK.RDS")
    )
  )
  
  km_output_pdf <- file.path(plots_dir, glue::glue("km_{file_names[group]}_OS_EFS_EUR_nonEUR_raceUK.pdf"))
  
  km_plot <- plotKM(model = list(km_os_result, km_efs_result),
                    variable = "EUR_status",
                    combined = T, 
                    title = titles[group])
  
  ggsave(km_output_pdf, km_plot,
         width = 10, height = 7, units = "in",
         device = "pdf")
  
  
}


for (group in groups){
  
  if (grepl("Low-grade glioma|LGG|Mixed neuronal-glial tumor|GNG", group)){
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_OS_additive_terms_subtype_resection_predicted_ancestry_raceUK.RDS")
      ))
  }else{
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_OS_additive_terms_subtype_predicted_ancestry_raceUK.RDS")
      ))
  }
  
  forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[group]}_OS_subtype_predicted_ancestry_raceUK.pdf"))
  
  
  forest_plot <- plotForest(survival_result)
  
  ggsave(forest_pdf, forest_plot, width = 8, height = 3)
  


  if (grepl("Low-grade glioma|LGG|Mixed neuronal-glial tumor|GNG", group)){
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_OS_additive_terms_subtype_resection_EUR_nonEUR_raceUK.RDS")
      ))
  }else{
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_OS_additive_terms_subtype_EUR_nonEUR_raceUK.RDS")
      ))
  }

  forest_pdf <- file.path(plots_dir,
                          glue::glue("forest_{file_names[group]}_OS_subtype_EUR_nonEUR_raceUK.pdf"))


  forest_plot <- plotForest(survival_result)

  ggsave(forest_pdf, forest_plot, width = 8, height = 3)
  
}




for (group in groups[2]){
  
  if (grepl("Low-grade glioma|LGG|Mixed neuronal-glial tumor|GNG", group)){
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_EFS_additive_terms_subtype_resection_predicted_ancestry_raceUK.RDS")
      ))
  }else{
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_EFS_additive_terms_subtype_predicted_ancestry_raceUK.RDS")
      ))
  }
  
  forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[group]}_EFS_subtype_predicted_ancestry_raceUK.pdf"))
  
  
  forest_plot <- plotForest(survival_result)
  
  ggsave(forest_pdf, forest_plot, width = 8, height = 3)
  
  
  

  if (grepl("Low-grade glioma|LGG|Mixed neuronal-glial tumor|GNG", group)){
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_EFS_additive_terms_subtype_resection_EUR_nonEUR_raceUK.RDS")
      ))
  }else{
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_EFS_additive_terms_subtype_EUR_nonEUR_raceUK.RDS")
      ))
  }

  forest_pdf <- file.path(plots_dir,
                          glue::glue("forest_{file_names[group]}_EFS_subtype_EUR_nonEUR_raceUK.pdf"))


  forest_plot <- plotForest(survival_result)

  ggsave(forest_pdf, forest_plot, width = 8, height = 3)
  
}

