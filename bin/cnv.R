#!/usr/bin/env Rscript

plot_CNV_data <- function(calls_file, annot_file, CNV_plot_file, sample_id) {
  library(ggplot2)
  library(caTools)
  library(dplyr)
  library(ggrepel)
  library(tidyr)

  # Read Calls data
  Calls <- read.delim(file=calls_file, skip = 1)
  names(Calls) <- c("Chr", "Start", "End", "Bin", "Log2ratio", "Strand")
  Calls$Chr <- factor(Calls$Chr, levels = c(1:22,"X","Y"))
  
  # Read Annot data
  Annot <- read.delim(annot_file, header = FALSE)
  Annot <- Annot %>% 
    separate(V8, c(NA, "Gene"), sep = "=") %>% 
    separate(V6, c(NA, "Score"), sep = "=") %>% 
    separate(V7, c(NA, "LOG2CNT"), sep = "=", convert = TRUE) %>%
    select(c(1,2,6,7,8)) %>%
    filter(Score != "-1") %>%
    filter(Score != "1")
  names(Annot) <- c("Chr","Start", "Score", "LOG2CNT","Gene")
  Annot$Start <- Annot$Start - 1
  Annot$Chr <- as.factor(gsub("chr", "", Annot$Chr))
  
  # Merge Calls and Annot
  Calls <- left_join(Calls, Annot)
  
  # Calculate moving average
  Calls <- Calls %>% 
    group_by(Chr) %>% 
    arrange(Start) %>%
    mutate(movavg = zoo::rollmean(Log2ratio, k=100, fill = NA, align = "right", na.pad = TRUE))
  
  # Define y-axis limits
  y_min <- -2.5
  y_max <- 2.5
  
  # Create a new column to flag out-of-bounds points
  Calls <- Calls %>%
    mutate(
      Log2ratio_clipped = ifelse(Log2ratio > y_max, y_max, ifelse(Log2ratio < y_min, y_min, Log2ratio)),
      OutOfBounds = Log2ratio > y_max | Log2ratio < y_min
    )

  # Separate data for out-of-bounds points
  out_of_bounds_points <- Calls %>% filter(OutOfBounds)

  # Plot data
  p <- ggplot(Calls, aes(x = Start, y = Log2ratio_clipped)) +
    geom_point(colour = "grey", size = 1) +
    geom_point(data = out_of_bounds_points, aes(y = ifelse(Log2ratio > y_max, y_max, y_min)), color = "red", size = 2) +
    facet_grid(~ Chr, scales = "free_x", space = "free_x") +
    geom_line(aes(y = movavg), colour = "#e41a1c") +
    theme_classic() +
    theme(axis.line = element_blank(),
          axis.title.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.grid.minor.y = element_line(colour = "black", size = 0.5),
          panel.spacing.x = unit(0, "lines"),
          strip.text.x = element_text(size = 7),
          strip.background = element_blank(),
          panel.background = element_rect(colour = "grey", size = 0.4)) +
    scale_y_continuous(minor_breaks = 0, limits = c(y_min, y_max)) +
    geom_label_repel(
      aes(x = Start, y = LOG2CNT, label = Gene, fontface = 2, nudge_y = LOG2CNT + 0.5),
      size = 3,
      color = "red",
      text.padding = 2,
      box.padding = unit(0.25, 'lines')
    ) +
    ggtitle(sample_id)
  
  # Save the plot to the specified file
  ggsave(filename = CNV_plot_file, plot = p, width = 18, height = 5)
}

# Get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# Ensure there are exactly 3 arguments: calls_file, annot_file, and CNV_plot_file
if (length(args) != 4) {
  stop("Please provide 3 arguments: calls_file, annot_file, and CNV_plot_file")
}

# Assign arguments to variables
calls_file <- args[1]
annot_file <- args[2]
CNV_plot_file <- args[3]
sample_id <- args[4]

# Call the function with the provided arguments
plot_CNV_data(calls_file, annot_file, CNV_plot_file, sample_id)

