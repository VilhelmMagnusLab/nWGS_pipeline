#!/usr/bin/env Rscript
# Corrected R script
merge_variant_caller_output <- function(Merged_file, Pileup_file, Somatic_file, output_file, occgenes) {
  
  library(dplyr)
  library(tidyr)
  
  # Load files
  Merged <- read.delim(Merged_file, header = TRUE, colClasses = c("character"))
  Merged <- Merged %>%
    select(7,1:6,8:16,19) %>%
    separate(Otherinfo13, c("GT","Merged_GQ","Depth","AD","AF"), sep = ":")
  Merged$callerM <- "Merged"
  
  Pileup <- read.delim(Pileup_file, header = TRUE, colClasses = c("character"))
  Pileup <- Pileup %>%
    select(7,1:6,8:16,19) %>%
    separate(Otherinfo13, c("GT","Pileup_GQ","Depth","AD","AF"), sep = ":")
  Pileup$callerP <- "Pileup"
  
  # Check if somatic file has data
  Somatic <- read.delim(Somatic_file, header = TRUE, colClasses = c("character"))
  has_somatic <- TRUE
  
  if (nrow(Somatic) == 0) {
    print("No somatic mutations reported")
    has_somatic <- FALSE
  } else {
    Somatic <- Somatic %>%
      select(7,1:6,9:16,18) %>%
    #  select(7,1:6,9:16,18) %>%
      separate(Otherinfo10, c("GT","ClairS_GQ","ClairS_Depth","ClairS_AF","q","w","e","r","t"), sep = ":") %>%
      select(1:15,17:19)
    Somatic$callerS <- "ClairS_TO"
  }
  occgenes <- readRDS(occgenes)
  print(occgenes)
  # Merge files
  All_calls <- full_join(Pileup, Merged) %>%
    unite(Variant_caller, callerP, callerM, sep = ", ", na.rm = TRUE)
  
  # Only join with Somatic data if it exists
  if (has_somatic) {
    All_calls <- All_calls %>%
      full_join(Somatic) %>%
      unite(Variant_caller, Variant_caller, callerS, sep = ", ", na.rm = TRUE)
  }
  
  # Unite GQ values based on availability of somatic data
  if (has_somatic) {
    All_calls <- All_calls %>%
      unite(GQ, Pileup_GQ, Merged_GQ, sep = ",", na.rm = TRUE) %>%
      unite(GQ, GQ, ClairS_GQ, sep = ",", na.rm = TRUE)
  } else {
    # If no somatic data, do not unite ClairS_GQ
    All_calls <- All_calls %>%
      unite(GQ, Pileup_GQ, Merged_GQ, sep = ",", na.rm = TRUE)
  }
  
  # All_calls <- All_calls %>%
  #   separate(COSMIC100, c("cosmic100_ID", "y"), sep = ";") %>%
  #   mutate(cosmic100_ID = gsub("ID=", "", cosmic100_ID)) %>%
  #   select(25,1:7,9,10,15,16,20:22,19,23,24)
  
  if ("COSMIC100" %in% colnames(All_calls)) {
    All_calls <- All_calls %>%
      separate(COSMIC100, c("cosmic100_ID", "y"), sep = ";", remove = FALSE, fill = "right") %>%
      mutate(cosmic100_ID = gsub("ID=", "", cosmic100_ID))
  } else {
    # Print a message or log that the 'COSMIC100' column is missing
    print("COSMIC100 column is missing, skipping separation and mutation steps for this column.")
  }
  
  # Continue with the rest of the pipeline
  All_calls <- All_calls %>%
    select(any_of(c(any_of(c(25,1:7,9,10,15,16,20:22,19,23,24)))))
  All_calls <- All_calls %>% filter(Gene.refGene %in% occgenes)
  print(All_calls$Gene.refGene)
  
  # Write output
  write.table(All_calls, file = output_file, row.names = FALSE, sep = "\t")
}

# Parsing arguments passed from the Nextflow process
args <- commandArgs(trailingOnly = TRUE)
Merged_file <- args[1]
Pileup_file <- args[2]
Somatic_file <- args[3]
output_file <- args[4]
occgenes <- args[5]


# Call the function with provided arguments
merge_variant_caller_output(Merged_file, Pileup_file, Somatic_file, output_file, occgenes)
