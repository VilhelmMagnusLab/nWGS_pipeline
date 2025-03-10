
plot_CNV_data <- function(sample_id) {
  library(ggplot2)
  library(caTools)
  library(dplyr)
  library(ggrepel)
  library(tidyr)
  
  # File paths
  Calls_file <- paste0("~/P24_pilot/Prospective/OCC/CNV_call_files/", sample_id, "_bins.bed")
  Annot_file <- paste0("~/P24_pilot/Prospective/OCC/CNV_call_files/", sample_id, ".annotated.cnv.csv")
  
  # Read Calls data
  Calls <- read.delim(file=Calls_file, skip = 1)
  names(Calls) <- c("Chr", "Start", "End", "Bin", "Log2ratio", "Strand")
  Calls$Chr <- factor(file=Calls$Chr, levels = c(1:22,"X","Y"))
  
  # Read Annot data
  Annot <- read.delim(Annot_file, header = FALSE)
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
  
  # Plot data
p <-  ggplot(Calls, aes(x=Start, y=Log2ratio)) +
    geom_point(colour = "grey", size = 1)+
    facet_grid(~ Chr, scales = "free_x", space = "free_x") +
    geom_line(aes(y = movavg), colour="#e41a1c")+
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
    scale_y_continuous(minor_breaks = 0, limits = c(-5,7))+
    geom_label_repel(
      aes(x=Start, y=LOG2CNT, label=Gene, fontface=2, nudge_y = LOG2CNT + .5),
      size= 3,
      color="red",
      text.padding = 2,
      box.padding = unit(.25, 'lines')
    )+
    ggtitle(sample_id)
  
  ggsave(filename = paste0("~/P24_pilot/Prospective/OCC/CNV_plots/", sample_id, "_CNV_plot.pdf"), plot = p, width = 18, height = 5)
  #dev.off()
 
}

# Get sample_id from command line arguments
sample_id <- commandArgs(trailingOnly = TRUE)

# Call the function with the provided sample_id
plot_CNV_data(sample_id)
