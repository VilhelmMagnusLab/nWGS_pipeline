#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6) {
    stop("Usage: plot_genomic_regions.R <gviz_data> <sample_id> <bam_file> <output_egfr_coverage> <output_idh1_coverage> <output_tertp_coverage>")
}

gviz_data_path <- args[1]
sample_id <- args[2]
bam_file <- args[3]
egfr_output <- args[4]
idh1_output <- args[5]
tertp_output <- args[6]

# Load required libraries
suppressPackageStartupMessages({
    library(Gviz)
    library(GenomicRanges)
    library(BSgenome)
    library(BSgenome.Hsapiens.UCSC.hg38)
})

# Load Gviz object
gvizobject <- load(gviz_data_path)

# Verify BAM file exists
if (!file.exists(bam_file)) {
    stop(sprintf("BAM file not found: %s", bam_file))
}

# Create offline IdeogramTrack function
createOfflineIdeogram <- function(chromosome) {
    # Create a basic ideogram without UCSC connection
    itrack <- IdeogramTrack(
        chromosome = chromosome,
        genome = "hg38",
        bands = NULL,
        showId = TRUE,
        size = 2
    )
    return(itrack)
}

tryCatch({
    # Plot EGFR coverage
    pdf(egfr_output, width=8, height=4)  # Reduce PDF dimensions
    itrack7 <- createOfflineIdeogram("7")
    Sample_track <- AlignmentsTrack(
        bam_file, 
        name = "EGFR",
        type = "pileup",      
        size = 0.3,           
        col.coverage = "darkblue"  
    )
    
    # Narrow the genomic region to focus on EGFR
    ht <- HighlightTrack(
        trackList = Sample_track,
        start = 55142193,     # Single start position
        end = 55154167,       # Single end position
        chromosome = "7",
        col = "red",
        fill = "#FFE5E5"
    )
    
    plotTracks(
        list(itrack7, EGFR_annot, ht), 
        from = 55142000,      
        to = 55155000,        
        chromosome = "7",
        sizes = c(0.2, 0.3, 0.5)
    )
    dev.off()

    # Plot IDH1 p.R132
    pdf(idh1_output, width=8, height=4)
    itrack2 <- createOfflineIdeogram("2")
    sTrack <- SequenceTrack(Hsapiens, chromosome = "2")
    Sample_track <- AlignmentsTrack(bam_file, name = "IDH1 p.132")
    ht <- HighlightTrack(
        trackList = list(sTrack, Sample_track),
        start = 208248388,
        end = 208248390,    # Add small range around position
        chromosome = "2",
        col = "red",
        fill = "#FFE5E5"
    )
    plotTracks(
        list(itrack2, ht), 
        chromosome = "2",
        from = 208248370, 
        to = 208248405
    )
    dev.off()

    # Plot TERTp
    pdf(tertp_output, width=8, height=4)
    itrack5 <- createOfflineIdeogram("5")
    sTrack <- SequenceTrack(Hsapiens, chromosome = "5")
    Sample_track <- AlignmentsTrack(bam_file, name = "TERTp")
    ht <- HighlightTrack(
        trackList = list(sTrack, Sample_track),
        start = 1295113,
        end = 1295135,
        chromosome = "5",
        col = "red",
        fill = "#FFE5E5"
    )
    plotTracks(
        list(itrack5, ht), 
        chromosome = "5",
        from = 1295103, 
        to = 1295145
    )
    dev.off()
}, error = function(e) {
    message("Error creating plots: ", e$message)
    quit(status = 1)
}) 