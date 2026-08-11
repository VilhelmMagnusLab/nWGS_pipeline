#!/usr/bin/env Rscript

# Parse command line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 8) {
    stop("Usage: plot_genomic_regions_v2.R <gviz_data> <sample_id> <bam_file> <output_egfr_coverage> <output_idh1_coverage> <output_idh2_coverage> <output_tertp_coverage> <cytoband_file>")
}

gviz_data_path <- args[1]
sample_id      <- args[2]
bam_file       <- args[3]
egfr_output    <- args[4]
idh1_output    <- args[5]
idh2_output    <- args[6]
tertp_output   <- args[7]
cytoband_file  <- args[8]

suppressPackageStartupMessages({
    library(Gviz)
    library(GenomicRanges)
    library(Rsamtools)
    library(GenomicAlignments)
})

# Load BSgenome for reference sequence display in IDH1/IDH2/TERTp plots.
# showMismatch=FALSE on AlignmentsTrack prevents sequenceLayer from running
# (avoiding the off-limits crash), while SequenceTrack still renders the
# reference base row that shows the amino acid context.
has_bsgenome <- tryCatch({
    suppressPackageStartupMessages(library(BSgenome.Hsapiens.UCSC.hg38))
    TRUE
}, error = function(e) {
    message("BSgenome.Hsapiens.UCSC.hg38 not available — sequence row will be omitted: ",
            conditionMessage(e))
    FALSE
})

gvizobject <- load(gviz_data_path)

if (!file.exists(bam_file)) {
    stop(sprintf("BAM file not found: %s", bam_file))
}


createCustomIdeogram <- function(chromosome) {
    bands_df <- read.table(
        cytoband_file,
        sep = "\t",
        col.names = c("chrom", "chromStart", "chromEnd", "name", "gieStain"),
        stringsAsFactors = FALSE
    )
    bands_df <- bands_df[bands_df$chrom == paste0("chr", chromosome), ]
    IdeogramTrack(
        genome     = "hg38",
        chromosome = chromosome,
        bands      = bands_df,
        showId     = TRUE,
        showBandId = TRUE,
        size       = 2
    )
}

# ── EGFR ─────────────────────────────────────────────────────────────────────
# AlignmentsTrack triggers sequenceLayer even with type="coverage" because Gviz
# auto-detects the hg38 genome from the BAM header and loads BSgenome.
# Fix: compute coverage with Rsamtools and plot with DataTrack, which never
# calls sequenceLayer. Rsamtools is a dependency of Gviz so it should be
# present; fall back to gene-annotation-only plot if not.
EGFR_FROM <- 55019017
EGFR_TO   <- 55211628

pdf(egfr_output, width = 10, height = 6)
tryCatch({
    itrack <- createCustomIdeogram("7")
    gtrack <- GenomeAxisTrack()

    message("EGFR_annot class: ", paste(class(EGFR_annot), collapse = ", "))
    annot_is_seq <- inherits(EGFR_annot, "SequenceTrack")
    egfr_annot_tracks <- if (annot_is_seq) list() else list(EGFR_annot)

    has_rsamtools <- requireNamespace("Rsamtools", quietly = TRUE) &&
                     requireNamespace("GenomicAlignments", quietly = TRUE)

    if (has_rsamtools) {
        # Compute per-base coverage; try UCSC (chr7) then Ensembl (7) naming
        bam_obj  <- Rsamtools::BamFile(bam_file)
        egfr_cov <- NULL
        for (chr_name in c("chr7", "7")) {
            egfr_cov <- tryCatch({
                gr  <- GRanges(chr_name, IRanges(EGFR_FROM, EGFR_TO))
                cov <- GenomicAlignments::coverage(bam_obj,
                           param = Rsamtools::ScanBamParam(which = gr))[[chr_name]]
                as.numeric(cov[EGFR_FROM:EGFR_TO])
            }, error = function(e) NULL)
            if (!is.null(egfr_cov)) { message("BAM chr name: ", chr_name); break }
        }
        if (is.null(egfr_cov)) stop("Could not read BAM coverage for EGFR region")

        cov_pos   <- EGFR_FROM:EGFR_TO
        cov_track <- DataTrack(
            data       = egfr_cov,
            start      = cov_pos,
            end        = cov_pos,
            chromosome = "chr7",
            genome     = "hg38",
            name       = "EGFR",
            type       = "h",
            col        = "steelblue",
            fill       = "steelblue"
        )
        ht <- HighlightTrack(
            trackList  = cov_track,
            start      = c(55142193, 55154167),
            width      = 10,
            chromosome = "chr7"
        )
        plotTracks(
            c(list(itrack, gtrack), egfr_annot_tracks, list(ht)),
            from       = EGFR_FROM,
            to         = EGFR_TO,
            chromosome = "chr7",
            cex        = 0.9
        )
    } else {
        # Rsamtools not available — show gene annotation only (no coverage)
        message("Rsamtools not available — plotting EGFR annotation only")
        plotTracks(
            c(list(itrack, gtrack), egfr_annot_tracks),
            from       = EGFR_FROM,
            to         = EGFR_TO,
            chromosome = "7",
            cex        = 0.9
        )
    }
}, error = function(e) {
    message("EGFR plot failed: ", conditionMessage(e))
    plot.new()
    text(0.5, 0.5, paste("EGFR plot error:\n", conditionMessage(e)), cex = 0.8)
})
dev.off()

# ── IDH1 p.R132 ──────────────────────────────────────────────────────────────
# Chromosome 2, reverse-strand gene.
# Strategy: SequenceTrack shows the reference base row (amino acid context);
# showMismatch=FALSE prevents sequenceLayer from being called so ONT reads
# that extend far beyond the 35 bp window cannot trigger the off-limits crash.
# Falls back to pileup-only if SequenceTrack still causes an error.
pdf(idh1_output, width = 10, height = 6)
tryCatch({
    itrack <- createCustomIdeogram("2")
    gtrack <- GenomeAxisTrack()

    idh1_plotted <- FALSE
    if (has_bsgenome) {
        tryCatch({
            sTrack       <- SequenceTrack(Hsapiens, chromosome = "2")
            Sample_track <- AlignmentsTrack(bam_file, name = "IDH1 p.R132",
                                            reverseStacking = TRUE,
                                            showMismatch    = FALSE)
            ht <- HighlightTrack(
                trackList  = list(sTrack, Sample_track),
                start      = c(208248387),
                width      = 2,
                chromosome = "2"
            )
            plotTracks(
                list(itrack, gtrack, ht),
                chromosome = "2",
                from       = 208248370,
                to         = 208248405,
                type       = "pileup",
                cex        = 0.9
            )
            idh1_plotted <- TRUE
        }, error = function(e) {
            message("IDH1 SequenceTrack attempt failed, using pileup fallback: ",
                    conditionMessage(e))
        })
    }
    if (!idh1_plotted) {
        Sample_track <- AlignmentsTrack(bam_file, name = "IDH1 p.R132",
                                        reverseStacking = TRUE)
        ht <- HighlightTrack(
            trackList  = Sample_track,
            start      = c(208248387),
            width      = 2,
            chromosome = "2"
        )
        plotTracks(
            list(itrack, gtrack, ht),
            chromosome = "2",
            from       = 208248370,
            to         = 208248405,
            type       = "pileup",
            cex        = 0.9
        )
    }
}, error = function(e) {
    message("IDH1 plot failed: ", conditionMessage(e))
    plot.new()
    text(0.5, 0.5, paste("IDH1 plot error:\n", conditionMessage(e)), cex = 0.8)
})
dev.off()

# ── IDH2 p.R172 ──────────────────────────────────────────────────────────────
# Chromosome 15, reverse-strand gene.
pdf(idh2_output, width = 10, height = 6)
tryCatch({
    itrack <- createCustomIdeogram("15")
    gtrack <- GenomeAxisTrack()

    idh2_plotted <- FALSE
    if (has_bsgenome) {
        tryCatch({
            sTrack       <- SequenceTrack(Hsapiens, chromosome = "15")
            Sample_track <- AlignmentsTrack(bam_file, name = "IDH2 p.R172",
                                            reverseStacking = TRUE,
                                            showMismatch    = FALSE)
            ht <- HighlightTrack(
                trackList  = list(sTrack, Sample_track),
                start      = c(90088605),
                width      = 2,
                chromosome = "15"
            )
            plotTracks(
                list(itrack, gtrack, ht),
                chromosome = "15",
                from       = 90088587,
                to         = 90088622,
                type       = "pileup",
                cex        = 0.9
            )
            idh2_plotted <- TRUE
        }, error = function(e) {
            message("IDH2 SequenceTrack attempt failed, using pileup fallback: ",
                    conditionMessage(e))
        })
    }
    if (!idh2_plotted) {
        Sample_track <- AlignmentsTrack(bam_file, name = "IDH2 p.R172",
                                        reverseStacking = TRUE)
        ht <- HighlightTrack(
            trackList  = Sample_track,
            start      = c(90088605),
            width      = 2,
            chromosome = "15"
        )
        plotTracks(
            list(itrack, gtrack, ht),
            chromosome = "15",
            from       = 90088587,
            to         = 90088622,
            type       = "pileup",
            cex        = 0.9
        )
    }
}, error = function(e) {
    message("IDH2 plot failed: ", conditionMessage(e))
    plot.new()
    text(0.5, 0.5, paste("IDH2 plot error:\n", conditionMessage(e)), cex = 0.8)
})
dev.off()

# ── TERTp ─────────────────────────────────────────────────────────────────────
# Chromosome 5, reverse-strand promoter region.
pdf(tertp_output, width = 10, height = 6)
tryCatch({
    itrack <- createCustomIdeogram("5")
    gtrack <- GenomeAxisTrack()

    tertp_plotted <- FALSE
    if (has_bsgenome) {
        tryCatch({
            sTrack       <- SequenceTrack(Hsapiens, chromosome = "5")
            Sample_track <- AlignmentsTrack(bam_file, name = "TERTp",
                                            reverseStacking = TRUE,
                                            showMismatch    = FALSE)
            ht <- HighlightTrack(
                trackList  = list(sTrack, Sample_track),
                start      = c(1295113, 1295135),
                width      = 0,
                chromosome = "5",
                name       = "TERTp"
            )
            plotTracks(
                list(itrack, gtrack, ht),
                chromosome = "5",
                from       = 1295103,
                to         = 1295145,
                type       = "pileup",
                cex        = 0.9
            )
            tertp_plotted <- TRUE
        }, error = function(e) {
            message("TERTp SequenceTrack attempt failed, using pileup fallback: ",
                    conditionMessage(e))
        })
    }
    if (!tertp_plotted) {
        Sample_track <- AlignmentsTrack(bam_file, name = "TERTp",
                                        reverseStacking = TRUE)
        ht <- HighlightTrack(
            trackList  = Sample_track,
            start      = c(1295113, 1295135),
            width      = 0,
            chromosome = "5",
            name       = "TERTp"
        )
        plotTracks(
            list(itrack, gtrack, ht),
            chromosome = "5",
            from       = 1295103,
            to         = 1295145,
            type       = "pileup",
            cex        = 0.9
        )
    }
}, error = function(e) {
    message("TERTp plot failed: ", conditionMessage(e))
    plot.new()
    text(0.5, 0.5, paste("TERTp plot error:\n", conditionMessage(e)), cex = 0.8)
})
dev.off()
