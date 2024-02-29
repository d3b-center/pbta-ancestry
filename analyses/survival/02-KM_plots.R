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


groups <- c("Atypical Teratoid Rhabdoid Tumor",
            "Craniopharyngioma",
            "DIPG or DMG|Other high-grade glioma",
            "Ependymoma",
            "Mixed neuronal-glial tumor",
            "Low-grade glioma",
            "Medulloblastoma",
            "Meningioma",
            "Neurofibroma plexiform",
            "HGG, H3 wildtype",
            "HGG, H3 K28",
            "LGG, wildtype",
            "LGG, BRAF mutant",
            "GNG, wildtype", 
            "GNG, BRAF mutant",
            "MB, Group4"
            )


file_names <- c("atrt", "cranio", "hgg", "epn", "gng", "lgg", "mb", "mng", "nfp",
                "hgg_wildtype", "hgg_DMG", "lgg_wildtype", "lgg_-BRAF",
                "gng_wildtype", "gng_-BRAF", "mb_Group4")

names(file_names) <- groups

titles <- c("ATRT", "Craniopharyngioma", "HGG", "Ependymoma", "GNG", "LGG", "Medulloblastoma", "Meningioma", "Neurofibroma Plexiform",
            "HGG, H3 wildtype", "HGG, H3 K28", "LGG, wildtype", "LGG, BRAF fusion", "GNG, wildtype", "GNG, BRAF fusion",
            "MB, Group4")
names(titles) <- groups


for (group in groups){
  
  km_os_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{file_names[group]}_OS_predicted_ancestry.RDS")
    )
  )
  
  km_pfs_result <- read_rds(
    file.path(input_dir,
              glue::glue("logrank_{file_names[group]}_PFS_predicted_ancestry.RDS")
    )
  )
  
  km_output_pdf <- file.path(plots_dir, glue::glue("km_{file_names[group]}_OS_PFS_predicted_ancestry.pdf"))
  
  km_plot <- plotKM(model = list(km_os_result, km_pfs_result),
                    variable = "predicted_ancestry",
                    combined = T, 
                    title = titles[group])
  
  ggsave(km_output_pdf, km_plot,
         width = 10, height = 7, units = "in",
         device = "pdf")
  
  
}



groups <- groups[!grepl("HGG,|LGG,|GNG,|MB,", groups)]

for (group in groups){
  
  if (group == "Low-grade glioma"){
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_OS_additive_terms_subtype_resection_predicted_ancestry.RDS")
      ))
  }else{
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_OS_additive_terms_subtype_predicted_ancestry.RDS")
      ))
  }
  
  forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[group]}_OS_subtype_predicted_ancestry.pdf"))
  
  
  forest_plot <- plotForest(survival_result)
  
  ggsave(forest_pdf, forest_plot, width = 8, height = 3)
  
}




for (group in groups){
  
  if (group == "Low-grade glioma"){
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_PFS_additive_terms_subtype_resection_predicted_ancestry.RDS")
      ))
  }else{
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[group]}_PFS_additive_terms_subtype_predicted_ancestry.RDS")
      ))
  }
  
  forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[group]}_PFS_subtype_predicted_ancestry.pdf"))
  
  
  forest_plot <- plotForest(survival_result)
  
  ggsave(forest_pdf, forest_plot, width = 8, height = 3)
  
}


subtypes <- c("HGG, H3 wildtype",
              "HGG, H3 K28",
              "LGG, wildtype",
              "LGG, BRAF mutant",
              "GNG, wildtype", 
              "GNG, BRAF mutant",
              "MB, Group4")

file_names <- c("hgg_wildtype", "hgg_DMG", "lgg_wildtype", "lgg_-BRAF",
                "gng_wildtype", "gng_-BRAF", "mb_Group4")

names(file_names) <- subtypes




for (subtype in subtypes){
  
  if (grepl("LGG|GNG", subtype)){
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[subtype]}_OS_additive_terms_subtype_resection_predicted_ancestry.RDS")
      ))
  }else{
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[subtype]}_OS_additive_terms_subtype_predicted_ancestry.RDS")
      ))
  }
  
  forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[subtype]}_OS_predicted_ancestry.pdf"))
  
  
  forest_plot <- plotForest(survival_result)
  
  ggsave(forest_pdf, forest_plot, width = 8, height = 3)
  
}




for (subtype in subtypes){
  
  if (grepl("LGG|GNG", subtype)){
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[subtype]}_PFS_additive_terms_subtype_resection_predicted_ancestry.RDS")
      ))
  }else{
    survival_result <- read_rds(
      file.path(input_dir,
                glue::glue("cox_{file_names[subtype]}_PFS_additive_terms_subtype_predicted_ancestry.RDS")
      ))
  }
  
  forest_pdf <- file.path(plots_dir, 
                          glue::glue("forest_{file_names[subtype]}_PFS_predicted_ancestry.pdf"))
  
  
  forest_plot <- plotForest(survival_result)
  
  ggsave(forest_pdf, forest_plot, width = 8, height = 3)
  
}
