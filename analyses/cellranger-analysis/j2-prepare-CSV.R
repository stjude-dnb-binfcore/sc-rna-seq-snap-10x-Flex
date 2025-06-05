# Convert excelt o TSV file

# Load library
suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
  library(yaml)
})


#################################################################################
# load config file
configFile <- paste0("../../project_parameters.Config.yaml")
if (!file.exists(configFile)){
  cat("\n Error: configuration file not found:", configFile)
  stop("Exit...")}

# read `yaml` file defining the `params` of the project and strategy analysis
yaml <- read_yaml(configFile)

#################################################################################
# Set up directories and paths to file Inputs/Outputs
root_dir <- yaml$root_dir
metadata_dir <- yaml$metadata_dir
metadata_file <- yaml$metadata_file
cellranger_parameters <- yaml$cellranger_parameters

analysis_dir <- file.path(root_dir, "analyses", "cellranger-analysis") 
#results_dir <- file.path(analysis_dir, "results", "03_cellranger_count_summary", as.character(glue::glue("{cellranger_parameters}")), as.character(glue::glue("multi_run_{cellranger_parameters}"))) 

# Input data
input_dir <- yaml$data_dir
#input_dir 
results_dir <- file.path(input_dir, "outs", "per_sample_outs") 


#################################################################################
# Read metadata file and define `sample_name`
project_metadata_file <- file.path(metadata_dir, metadata_file) # metadata input file

# Read metadata file and define `sample_name`
project_metadata <- read.csv(project_metadata_file, sep = "\t", header = TRUE)
sample_name <- unique(as.character(project_metadata$ID))
sample_name <- sort(sample_name, decreasing = FALSE)
print(sample_name)
#################################################################################

##################################################################################################################################################################
# Process each sample
for (i in seq_along(sample_name)) {
  sample_input_dir <- file.path(input_dir, "outs", "per_sample_outs", sample_name[i]) 
  #print(sample_input_dir)
  input_file <- c(dir(path = sample_input_dir,  pattern = "metrics_summary.csv", full.names = TRUE, recursive = TRUE))
  print(input_file)
 
  if (length(input_file) == 0) {
    warning(glue::glue("No metrics_summary.csv found for sample {sample_name[i]}"))
    next
  } else if (length(input_file) > 1) {
    warning(glue::glue("Multiple metrics_summary.csv files found for sample {sample_name[i]}; using the first one."))
    input_file <- input_file[1]
  }
  
  # No numeric conversion is done — we’re preserving original text like "33,635 (81.37%)".
  
  processed_input_df <- read.csv(input_file, header = TRUE, stringsAsFactors = FALSE) %>%
    select(Metric.Name, Metric.Value) %>%
    rename(Metric_Name = dplyr::matches("Metric.?Name"),
           Metric_Value = dplyr::matches("Metric.?Value")) %>%
    pivot_wider(names_from = Metric_Name, values_from = Metric_Value, values_fn = dplyr::first) %>% # or `list` if you want to preserve duplicates
    write_csv(file.path(results_dir, sample_name[i], paste0("metrics_summary_updated.csv")))   # save data 
  
  head(processed_input_df)
}

##################################################################################################################################################################
