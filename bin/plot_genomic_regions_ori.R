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

# Plot EGFR coverage
pdf(egfr_output, width=10, height=6)
itrack7 <- IdeogramTrack(genome = "hg38", chromosome = 7)
Sample_track <- AlignmentsTrack(bam_file, name = "EGFR")
ht <- HighlightTrack(trackList = Sample_track,
                    start = c(55142193, 55154167), width = 10,
                    chromosome = 7)
plotTracks(list(itrack7, EGFR_annot, ht), 
          from = 55019017, to = 55211628,
          main = paste("EGFR Coverage -", sample_id))
dev.off()

# Plot IDH1 p.R132
pdf(idh1_output, width=10, height=6)
itrack2 <- IdeogramTrack(genome = "hg38", chromosome = 2)
sTrack <- SequenceTrack(Hsapiens, chromosome = "chr2", genome = "hg38")
Sample_track <- AlignmentsTrack(bam_file, name = "IDH1 p.132")
ht <- HighlightTrack(trackList = list(sTrack, Sample_track),
                    start = c(208248388), width = 2,
                    chromosome = 2)
plotTracks(list(itrack2, ht), 
          from = 208248370, to = 208248405, 
          cex = 0.9, cex.mismatch = 0.5,
          main = paste("IDH1 p.R132 -", sample_id))
dev.off()

# Plot TERTp
pdf(tertp_output, width=10, height=6)
itrack5 <- IdeogramTrack(genome = "hg38", chromosome = 5)
sTrack <- SequenceTrack(Hsapiens, chromosome = "chr5", genome = "hg38")
Sample_track <- AlignmentsTrack(bam_file, name = "TERTp")
ht <- HighlightTrack(trackList = list(sTrack, Sample_track),
                    start = c(1295113, 1295135), width = 0,
                    chromosome = 5, name = "TERTp")
plotTracks(list(itrack5, ht), 
          from = 1295103, to = 1295145, 
          cex = 0.9, cex.mismatch = 0.5,
          main = paste("TERTp -", sample_id))
dev.off() 