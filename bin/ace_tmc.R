#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) {
    stop("Usage: ace_tmc.R <input_bam_folder> <output_dir> <sample_id>")
}

input_bam_folder <- args[1]
output_dir <- args[2]
sample_id <- args[3]

# Load required libraries
suppressPackageStartupMessages({
    library(ACE)
    library(QDNAseq)
    library(QDNAseq.hg38)
})

# Create output directory
dir.create(output_dir, showWarnings = FALSE)

# Run ACE analysis
runACE(
    input = input_bam_folder,
    outputdir = output_dir,
    genome = 'hg38',
    filetype = 'bam',
    ploidies = 2,
    imagetype = 'png',
    autopick = TRUE,
    method = 'RMSE',
    binsizes = 1000
)

# Read the fitpicker file and extract threshold
fitpicker_file <- file.path(output_dir, 
                           "1000kbp",
                           "2N", 
                           "fitpicker_2N.tsv")

# Ensure the file exists before reading
if (!file.exists(fitpicker_file)) {
    stop(paste("Fitpicker file not found:", fitpicker_file))
}

fit_data <- read.table(fitpicker_file, header=TRUE, sep="\t")
threshold_value <- fit_data$likely_fit[1]

# Write threshold to file
threshold_output <- file.path(output_dir, "threshold_value.txt")
cat(sprintf("%.4f", threshold_value), file=threshold_output) 
