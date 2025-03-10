#!/usr/bin/env Rscript

plot_CNV_data <- function(calls_file, annot_file, CNV_plot_file, html_plot, sample_id) {
  library(ggplot2)
  library(caTools)
  library(dplyr)
  library(ggrepel)
  library(tidyr)
  library(htmltools)
  #install.packages("zoo")
  #library(zoo)  # for rollmean

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
  
  # Calculate dynamic y-axis limits based on Log2ratio
  y_min <- -2.5 #min(Calls$Log2ratio, na.rm = TRUE) - 1  # add a bit of padding
  y_max <- max(Calls$Log2ratio, na.rm = TRUE) + 1  # add a bit of padding

  # Plot data
  p <- ggplot(Calls, aes(x=Start, y=Log2ratio)) +
    geom_point(colour = "grey", size = 1) +
    facet_grid(~ Chr, scales = "free_x", space = "free_x") +
    geom_line(aes(y = movavg), colour="#e41a1c") +
    theme_classic() +
    theme(axis.line=element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.grid.minor.y=element_line(colour="black",size=0.5),
          panel.spacing.x=unit(0, "lines"),
          strip.text.x = element_text(size = 7),
          strip.background = element_blank(),
          panel.background = element_rect(colour = "grey", size=0.4)) +
    scale_y_continuous(minor_breaks = 0, limits = c(y_min, y_max)) +
    geom_label_repel(
      aes(x=Start, y=LOG2CNT, label=Gene, fontface=2, nudge_y = LOG2CNT + .5),
      size= 3,
      color="red",
      text.padding = 2,
      box.padding = unit(.25, 'lines')
    ) +
    ggtitle(sample_id)
  
  # Save the plot to the specified file
  ggsave(filename = CNV_plot_file, plot = p, width = 18, height = 5)
  html_plot <- as.tags(p) 
  save_html(html_plot, file="html_plot") 
}

# Get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# Ensure there are exactly 3 arguments: calls_file, annot_file, and CNV_plot_file
if (length(args) != 5) {
  stop("Please provide 3 arguments: calls_file, annot_file, and CNV_plot_file")
}

# Assign arguments to variables
calls_file <- args[1]
annot_file <- args[2]
CNV_plot_file <- args[3]
html_plot <- args[4]
sample_id <- args[5]

# Call the function with the provided arguments
plot_CNV_data(calls_file, annot_file, CNV_plot_file, html_plot, sample_id)

