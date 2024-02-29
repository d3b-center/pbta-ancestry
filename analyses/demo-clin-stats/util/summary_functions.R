# Functions for generating summary statistics 
#
# Ryan Corbett
#
# 2023

# Function to stratify patient cohort by supplied variable, and generate numnber of frequency of patients in each group

summarize_freq <- function(hist, var){
  
  freq <- hist %>%
    count(!!rlang::ensym(var)) %>%
    dplyr::mutate(`Percent Cohort` = round(n/nrow(hist)*100, 1)) %>%
    dplyr::rename("Number Cohort" = "n") %>%
    dplyr::mutate(`Number Cohort` = glue::glue("{`Number Cohort`} ({`Percent Cohort`}%)")) %>%
    dplyr::select(-`Percent Cohort`) %>%
    column_to_rownames(var) %>%
    t() %>%
    as.data.frame()
  
  return(freq)
  
}


# Function to stratify patient cohort by two supplied variables, and generate number of patients in each group

summarize_count <- function(hist, var1, var2){
  
  summary <- hist %>%
    count(!!rlang::ensym(var1), !!rlang::ensym(var2)) %>%
    spread(!!rlang::ensym(var2), n) 
  
  return(summary)
}
