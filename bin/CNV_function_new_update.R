#!/usr/bin/env Rscript

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 7){
  stop("Usage: CNV_function_new.R <calls_file> <annot_file> <segs_file> <out_new_cnv_plot> <out_chr9_plot> <out_chr7_plot> <sample_id>")
}

calls_file     <- args[1]
annot_file     <- args[2]
segs_file      <- args[3]
out_new_cnv_plot <- args[4]
out_chr9_plot    <- args[5]
out_chr7_plot    <- args[6]
sample_id        <- args[7]

# Load required libraries
library(ggplot2)
library(caTools)
library(dplyr)
library(ggrepel)
library(tidyr)
library(ggpp)

plot_CNV_data <- function(calls_file, annot_file, segs_file, out_new_cnv_plot, out_chr9_plot, out_chr7_plot, sample_id) {
  # Read Calls data
  Calls <- read.delim(calls_file, skip = 1, header = FALSE)
  names(Calls) <- c("Chr", "Start", "End", "Bin", "Log2ratio", "Strand")
  #print(Calls)
  
  # Read segmentation data
  Segs <- read.delim(segs_file, skip = 1, header = FALSE)
  names(Segs) <- c("Chr", "Start", "End", "Bin", "Seg", "Strand")
  
  Calls <- left_join(Calls, Segs)
  Calls$Chr <- factor(Calls$Chr, levels = c(1:22, "X", "Y"))
  
  chrom9 <- Calls %>% filter(Chr == 9)
  
  # Read Annot data
  Annot <- read.delim(annot_file, header = FALSE)
  #print(Annot)
  Annot <- Annot %>% 
    separate(V8, c(NA, "Gene"), sep = "=") %>% 
    separate(V6, c(NA, "Score"), sep = "=") %>% 
    separate(V7, c(NA, "LOG2CNT"), sep = "=", convert = TRUE) %>%
    select(c(1,2,6,7,8)) %>%
    filter(Score != "-1") %>%
    filter(Score != "1") %>%
    filter(Gene != "CRLF2")
  names(Annot) <- c("Chr", "Start", "Score", "LOG2CNT", "Gene")
  Annot$Start <- Annot$Start - 1
  Annot$Chr <- as.factor(gsub("chr", "", Annot$Chr))
  
  Annot9 <- Annot %>% filter(Gene == "CDKN2A")
  
  # Merge Calls and Annot
  Calls <- left_join(Calls, Annot)
  
  lim <- 3
  offset <- 0.1
  Calls$Log2ratio_Capped <- ifelse(abs(Calls$Log2ratio) > lim, sign(Calls$Log2ratio)*lim + sign(Calls$Log2ratio)*offset, Calls$Log2ratio)
  Calls$Segs_Capped <- ifelse(abs(Calls$Seg) > lim, sign(Calls$Seg)*lim + sign(Calls$Seg)*offset, Calls$Seg)
  Calls$flag <- abs(Calls$Log2ratio) > lim
  Calls$flagSegs <- abs(Calls$Segs) > lim
  
  # Plot CNV data for all chromosomes
  p <- ggplot(Calls) +
    geom_point(data = subset(Calls, !flag), aes(x = Start, y = Log2ratio), colour = "grey", size = 1) +
    geom_point(data = subset(Calls, flag), aes(x = Start, y = Log2ratio_Capped), colour = "red", shape = 17, size = 1) +
    facet_grid(~ Chr, scales = "free_x", space = "free_x") +
    geom_line(data = subset(Calls, !flagSegs), aes(x = Start, y = Seg), colour = "#e41a1c") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor.y = element_line(colour = "black", size = 0.5),
          panel.spacing.x = unit(0, "lines"),
          strip.text.x = element_text(size = 12),
          strip.background = element_blank(),
          panel.background = element_rect(colour = "grey", size = 0.4)) +
    scale_y_continuous(minor_breaks = 0)
  
  p2 <- p + geom_label_s(
    aes(x = Start, y = Log2ratio_Capped, label = Gene, fontface = 2),
    size = 4,
    color = "red",
    box.padding = unit(0.25, "lines")
  ) +
    ggtitle(sample_id)
  
  ggsave(filename = out_new_cnv_plot, plot = p2, width = 18, height = 5)
  
  ################## Chrom 9 Plot
  
  chrom9 <- left_join(chrom9, Annot9)
  chrom9[1,10] <- ""
  chrom9$Gene[chrom9$Start == 21950000] <- "CDKN2A" 
  chrom9$Gene[chrom9$Start == 22000000] <- "CDKN2B" 
  
  p9 <- ggplot(chrom9, aes(x = Start, y = Log2ratio, label = Gene)) +
    geom_point(colour = "grey", size = 1) +
    geom_line(aes(y = Seg), colour = "#e41a1c", linewidth = 1) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor.y = element_line(colour = "black", size = 0.5),
          panel.spacing.x = unit(0, "lines"),
          strip.text.x = element_text(size = 12),
          strip.background = element_blank(),
          panel.background = element_rect(colour = "grey", size = 0.4)) +
    scale_y_continuous(minor_breaks = 0)
  
  p92 <- p9 + geom_label_s(nudge_x = 5000000, 
                           nudge_y = -0.1,
                           point.padding = 0.4) +
    ggtitle(paste0(sample_id, " chromosome 9"))
  
  ggsave(filename = out_chr9_plot, plot = p92, width = 18, height = 5)
  
  ####################### Chrom 7 Plot
  
  chrom7 <- Calls %>% filter(Chr == 7)
  
  p7 <- ggplot(chrom7, aes(x = Start, y = Log2ratio, label = Gene)) +
    geom_point(colour = "grey", size = 1) +
    geom_line(aes(y = Seg), colour = "#e41a1c", linewidth = 1) +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor.y = element_line(colour = "black", size = 0.5),
          panel.spacing.x = unit(0, "lines"),
          strip.text.x = element_text(size = 12),
          strip.background = element_blank(),
          panel.background = element_rect(colour = "grey", size = 0.4)) +
    scale_y_continuous(minor_breaks = 0)
  
  p72 <- p7 + geom_label_s(nudge_x = 5000000,
                           nudge_y = -0.2,
                           point.padding = 0.4) +
    ggtitle(paste0(sample_id, " chromosome 7"))
  
  ggsave(filename = out_chr7_plot, plot = p72, width = 18, height = 5)
}

# Call the function with the parsed arguments
plot_CNV_data(calls_file, annot_file, segs_file, out_new_cnv_plot, out_chr9_plot, out_chr7_plot, sample_id)
