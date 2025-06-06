#!/usr/bin/env Rscript

plot_CNV_data <- function(calls_file, annot_file, CNV_plot_file, html_plot, sample_id) {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(plotly)  # For interactive HTML plots
  library(htmltools)

  # Read Calls data
  Calls <- read.delim(file = calls_file, skip = 1)
  names(Calls) <- c("Chr", "Start", "End", "Bin", "Log2ratio", "Strand")
  Calls$Chr <- factor(Calls$Chr, levels = c(1:22, "X", "Y"))

  # Read Annot data and process
  Annot <- read.delim(annot_file, header = FALSE)
  Annot <- Annot %>%
    separate(V8, c(NA, "Gene"), sep = "=") %>%
    separate(V6, c(NA, "Score"), sep = "=") %>%
    separate(V7, c(NA, "LOG2CNT"), sep = "=", convert = TRUE) %>%
    select(c(1, 2, 6, 7, 8)) %>%
    filter(Score != "-1") %>%
    filter(Score != "1")
  names(Annot) <- c("Chr", "Start", "Score", "LOG2CNT", "Gene")
  Annot$Start <- Annot$Start - 1
  Annot$Chr <- as.factor(gsub("chr", "", Annot$Chr))

  # Merge Calls and Annot
  Calls <- left_join(Calls, Annot)

  # Calculate moving average
  Calls <- Calls %>%
    group_by(Chr) %>%
    arrange(Start) %>%
    mutate(movavg = zoo::rollmean(Log2ratio, k = 100, fill = NA, align = "right", na.pad = TRUE))

  # Set y-axis limits based on Log2ratio
  y_min <- -2.5
  y_max <- max(Calls$Log2ratio, na.rm = TRUE) + 1

  # Plot data with rectangular border for each chromosome facet
  p <- ggplot(Calls, aes(x = Start, y = Log2ratio)) +
    geom_point(colour = "grey", size = 1) +
    facet_wrap(~ Chr, nrow = 1, scales = "free_x") +
    geom_line(aes(y = movavg), colour = "#e41a1c") +
    theme_classic() +
    theme(
      axis.line = element_blank(),
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.minor.y = element_line(colour = "black", size = 0.5),
      panel.spacing = unit(0.1, "lines"),
      strip.text.x = element_text(size = 7),
      strip.background = element_blank(),
      panel.border = element_rect(colour = "black", fill = NA, size = 0.5),  # Add border to each panel
      panel.background = element_rect(colour = "black", size = 0.4)
    ) +
    scale_y_continuous(minor_breaks = 0, limits = c(y_min, y_max)) +
    ggtitle(sample_id)

  # Save as PNG for the PDF version
  ggsave(filename = CNV_plot_file, plot = p, width = 18, height = 5)

  # Convert to plotly for interactive HTML
  p_plotly <- ggplotly(p)
  
  # Wrap the plot in a box with CSS styling
  html_content <- tags$div(
    style = "border: 2px solid #333; padding: 20px; margin: 20px; box-shadow: 5px 5px 10px #888;",
    p_plotly
  )
  
  # Save the HTML content
  save_html(html_content, file = html_plot)
}

# Get arguments from the command line
args <- commandArgs(trailingOnly = TRUE)

# Ensure there are exactly 5 arguments
if (length(args) != 5) {
  stop("Please provide 5 arguments: calls_file, annot_file, CNV_plot_file, html_plot, and sample_id")
}

# Assign arguments to variables
calls_file <- args[1]
annot_file <- args[2]
CNV_plot_file <- args[3]
html_plot <- args[4]
sample_id <- args[5]

# Call the function with the provided arguments
plot_CNV_data(calls_file, annot_file, CNV_plot_file, html_plot, sample_id)
