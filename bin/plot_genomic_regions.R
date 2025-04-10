#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 7) {
    stop("Usage: plot_genomic_regions.R <gviz_data> <sample_id> <bam_file> <output_egfr_coverage> <output_idh1_coverage> <output_tertp_coverage> <cytoband_file>")
}

gviz_data_path <- args[1]
sample_id <- args[2]
bam_file <- args[3]
egfr_output <- args[4]
idh1_output <- args[5]
tertp_output <- args[6]
cytoband_file <- args[7]

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

# Function to create custom ideogram with local band data
createCustomIdeogram <- function(chromosome) {
    # Read the local cytoband file from passed argument
    bands_df <- read.table(
        cytoband_file,  # Use passed file path
        sep="\t",
        col.names=c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
        stringsAsFactors=FALSE
    )
    
    # Filter for the requested chromosome
    bands_df <- bands_df[bands_df$chrom == paste0("chr", chromosome), ]
    
    # Create ideogram track with full band information
    itrack <- IdeogramTrack(
        genome = "hg38",
        chromosome = chromosome,
        bands = bands_df,
        showId = TRUE,
        showBandId = TRUE,
        size = 2
    )
    return(itrack)
}

# Plot EGFR coverage
pdf(egfr_output, width=10, height=6)
itrack <- createCustomIdeogram("7")
gtrack <- GenomeAxisTrack()
Sample_track <- AlignmentsTrack(bam_file, name = "EGFR")
ht <- HighlightTrack(
    trackList = Sample_track,
    start = c(55142193, 55154167), 
    width = 10,
    chromosome = "7"
)
plotTracks(
    list(itrack, gtrack, EGFR_annot, ht), 
    from = 55019017, 
    to = 55211628,
    chromosome = "7",
    #main = paste("EGFR Coverage -", sample_id),
    cex = 0.9,
    cex.mismatch = 0.5
)
dev.off()

# Plot IDH1 p.R132
pdf(idh1_output, width=10, height=6)
itrack <- createCustomIdeogram("2")
gtrack <- GenomeAxisTrack()
sTrack <- SequenceTrack(Hsapiens, chromosome = "2")
Sample_track <- AlignmentsTrack(bam_file, name = "IDH1 p.132")
ht <- HighlightTrack(
    trackList = list(sTrack, Sample_track),
    start = c(208248388), 
    width = 2,
    chromosome = "2"
)
plotTracks(
    list(itrack, gtrack, ht), 
    chromosome = "2",
    from = 208248370, 
    to = 208248405,
    #main = paste("IDH1 p.R132 -", sample_id),
    cex = 0.9,
    cex.mismatch = 0.5
)
dev.off()

# Plot TERTp
pdf(tertp_output, width=10, height=6)
itrack <- createCustomIdeogram("5")
gtrack <- GenomeAxisTrack()
sTrack <- SequenceTrack(Hsapiens, chromosome = "5")
Sample_track <- AlignmentsTrack(bam_file, name = "TERTp")
ht <- HighlightTrack(
    trackList = list(sTrack, Sample_track),
    start = c(1295113, 1295135), 
    width = 0,
    chromosome = "5",
    name = "TERTp"
)
plotTracks(
    list(itrack, gtrack, ht), 
    chromosome = "5",
    from = 1295103, 
    to = 1295145,
    #main = paste("TERTp -", sample_id),
    cex = 0.9,
    cex.mismatch = 0.5
)
dev.off() 