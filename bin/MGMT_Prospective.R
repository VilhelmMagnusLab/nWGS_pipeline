process_mgmt_data <- function(sample_id) {
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(fs) # For creating directories
  
  # File paths
###  bed_file <- file.path("~/P24_pilot/Prospective/Methylation/MGMT_promoter_methylation/", paste0(sample_id, ".MGMT.bed"))
###  save_path <- file.path("~/P24_pilot/Prospective/Methylation/MGMT_results/", paste0(sample_id, ".MGMT_results.csv"))
  save(.MGMT_results.csv, file ="${MGMT_results.csv}")
  bed_file <- read.csv("${mgmt_bed}", sep = "\\t", header=TRUE) 
  # Ensure the destination directory exists
  if (!dir_exists(dirname(save_path))) {
    dir_create(dirname(save_path))
  }
  
  # Step 1: Read the CSV file
  MGMT <- read_table(bed_file)
  
  # Step 2: Filter the data for modBase == "m"
  
  MGMT <- MGMT %>% filter(between(Start, 129466683, 129467448)) %>% filter(modBase %in% c("m", "5mC"))
  
  # Step 3: Create a ggplot (optional, can be commented if not needed)
  # ggplot(MGMT, aes(x = Start_hg38, y = methylated_frequency)) + geom_point()
  
  # Step 4: Filter for the specific range and calculate the mean Methylation
  Pyro <- MGMT %>% filter(Start > 129467253 & End < 129467274) %>% summarise(mean_methylation_pyro = mean(Methylation, na.rm = TRUE))
  
  # Step 5: Calculate the overall mean Methylation
  FULL <- MGMT %>% summarise(mean_methylation_full = mean(Methylation, na.rm = TRUE))
  
  # Step 6: Combine Pyro and FULL into a single data frame
  results <- bind_cols(Pyro, FULL) %>% mutate(sample_id = sample_id)
  
  results <- results %>% mutate(Classification_by_Pyro = case_when(
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
  
  # Step 7: Write the combined results to a CSV file in the destination folder
  write_csv(results, save_path)
  
  # Return the results as a list
  return(results)
}

# Get sample_id from command line arguments
sample_id <- commandArgs(trailingOnly = TRUE)

# # Call the function with the provided sample_id
process_mgmt_data(sample_id)
