#!/usr/bin/env Rscript

# Load necessary libraries
library(dplyr)
library(readr)
library(ggplot2)
library(fs) # For creating directories

process_mgmt_data <- function(mgmt_bed, save_path) {
  
  # Ensure the destination directory exists
  if (!dir_exists(dirname(save_path))) {
    dir_create(dirname(save_path))
  }
  
  # Step 1: Read the BED file (assumed to be tab-separated)
 MGMT <- read_delim(file=mgmt_bed, col_names = TRUE, delim = "\t") 
  # Step 2: Filter the data for modBase == "m" or "5mC"
  MGMT <- MGMT %>% filter(between(Start, 129466683, 129467448)) %>% filter(modBase %in% c("m", "5mC"))
  
  # Step 3: (Optional) Create a ggplot (comment out if not needed)
  # ggplot(MGMT, aes(x = Start, y = Methylation)) + geom_point()

  # Step 4: Filter for the Pyro kit range and calculate the mean methylation
  Pyro <- MGMT %>% filter(Start > 129467253 & End < 129467274) %>% summarise(mean_methylation_pyro = mean(Methylation, na.rm = TRUE))
  
  # Step 5: Calculate the overall mean methylation
  FULL <- MGMT %>% summarise(mean_methylation_full = mean(Methylation, na.rm = TRUE))
  
  # Step 6: Combine Pyro and FULL into a single data frame
  results <- bind_cols(Pyro, FULL) %>% mutate(sample_id = basename(mgmt_bed))
  
  # Step 7: Add classification based on methylation values
  results <- results %>%
    mutate(Classification_by_Pyro = case_when(
      mean_methylation_pyro < 10 ~ "Unmethylated",
      mean_methylation_pyro > 30 ~ "Methylated",
      mean_methylation_pyro > 10 & mean_methylation_pyro <= 22 ~ "Grey zone, Unmethylated",
      mean_methylation_pyro > 22 & mean_methylation_pyro <= 30 ~ "Grey zone, Methylated"
    )) %>%
    mutate(Classification_by_Full = case_when(
      mean_methylation_full < 25 ~ "Unmethylated",
      mean_methylation_full > 30 ~ "Methylated",
      mean_methylation_full > 24 & mean_methylation_full <= 28 ~ "Grey zone, Unmethylated",
      mean_methylation_full > 28 & mean_methylation_full <= 31 ~ "Grey zone, Methylated"
    ))

  # Step 8: Write the combined results to a CSV file in the destination folder
  write_csv(results, save_path)
  
  # Return the results as a list
  return(results)
}

# Get the file paths from command-line arguments
args <- commandArgs(trailingOnly = TRUE)
mgmt_bed <- args[1]  # First argument: mgmt_bed file path
save_path <- args[2] # Second argument: output path for MGMT_results.csv

# Call the function with the provided paths
process_mgmt_data(mgmt_bed, save_path)

