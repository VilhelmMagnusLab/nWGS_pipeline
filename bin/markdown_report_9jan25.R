#!/usr/bin/env Rscript

library(argparse)
library(dplyr)
library(knitr)
library(rmarkdown)
library(pdftools)
library(magick) # Ensure this library is loaded for image handling
library(ggplot2)

# Function to convert PDF to images
convert_pdf_to_images <- function(pdf_path, output_dir, max_width = 1600, max_height = 400) {
  if (!file.exists(pdf_path)) {
    stop("PDF file not found: ", pdf_path)
  }

  cat("Converting PDF to images...\n")
  images <- pdf_convert(pdf_path, format = "png", dpi = 300)
  image_paths <- vector("list", length(images))

  for (i in seq_along(images)) {
    image_path <- file.path(output_dir, paste0("pdf_image_", i, ".png"))
    img <- image_read(images[i])
    img_resized <- image_scale(img, paste0(max_width, "x", max_height))
    image_write(img_resized, image_path)
    image_paths[[i]] <- image_path
  }

  return(unlist(image_paths))
}

# Function to generate the markdown report
generate_markdown_report <- function(sample_id, html_cnv, csv_cnv, csv_methyl_path, csv_merge_annotation, html_svanna, csv_annotsv_path, html_terp, svanna_html, annotsv_html, output_md) {
  today_date <- Sys.Date()
  cat("Generating Markdown report...\n")

  markdown_content <- paste0(
    "<div style='text-align: center; font-size: 34px;'>",
    "<u><strong>", sample_id, ": Initial WGS Report - P24 Nanopore Sequencing</strong></u>",
    "</div>",
    "<span style='font-size: 18px;'><strong>Sample ID</strong></span>: <span style='font-size: 18px; color: red;'><strong>", sample_id, "</strong></span>",
    "<span style='font-size: 18px;'><strong>Whole genome coverage</strong></span>: <span style='font-size: 18px;color: red;'><strong> 47.7x </strong></span>",
    "<span style='font-size: 18px;'><strong>MGMT promoter methylation</strong></span>: <span style='font-size: 18px;color: red;'><strong> Unmethylated, low confidence (grey-zone) </strong></span>",
    "<hr style='border: 1px solid black;'/>"
  )

  # Add CNV section
  markdown_content <- paste0(markdown_content, "<h2><span style='color:blue;'>Copy Number Variation</span></h2>\n", html_cnv)

  # Add CNV table
  if (file.exists(csv_cnv)) {
    df_cnv <- read.csv(csv_cnv, header = TRUE)
    markdown_content <- paste0(
      markdown_content,
      "<h2><span style='color:blue;'>Copy Number Variation Filter Table</span></h2>\n",
      knitr::kable(head(df_cnv, 50), format = "markdown"), "\n\n"
    )
  } else {
    warning("CNV CSV file not found: ", csv_cnv)
  }

  # Add methylation section
  if (file.exists(csv_methyl_path)) {
    df_methyl <- read.csv(csv_methyl_path, header = TRUE)
    markdown_content <- paste0(
      markdown_content,
      "<h2><span style='color:blue;'>Methylation</span></h2>\n",
      knitr::kable(head(df_methyl, 50), format = "markdown"), "\n\n"
    )
  } else {
    warning("Methylation CSV file not found: ", csv_methyl_path)
  }

  # Add annotation section
  if (file.exists(csv_merge_annotation)) {
    df_annotation <- read.csv(csv_merge_annotation, sep = "\t", header = TRUE)
    markdown_content <- paste0(
      markdown_content,
      "<h2><span style='color:blue;'>SNV calling and annotation</span></h2>\n",
      knitr::kable(head(df_annotation, 50), format = "markdown"), "\n\n"
    )
  } else {
    warning("Annotation CSV file not found: ", csv_merge_annotation)
  }

  # Add AnnotSV section
  if (file.exists(csv_annotsv_path)) {
    df_annotsv <- read.csv(csv_annotsv_path, sep = "\t", header = TRUE)
    markdown_content <- paste0(
      markdown_content,
      "<h2><span style='color:blue;'>AnnotSV Structure Variants</span></h2>\n",
      knitr::kable(head(df_annotsv, 50), format = "markdown"), "\n\n"
    )
  } else {
    warning("AnnotSV CSV file not found: ", csv_annotsv_path)
  }

  # Save to markdown file
  writeLines(markdown_content, output_md)
}

# Convert Markdown to HTML
convert_markdown_to_html <- function(md_path, output_html) {
  rmarkdown::render(md_path, output_file = output_html, output_format = "html_document", self_contained = TRUE)
}

# Convert HTML to PDF
convert_html_to_pdf <- function(html_path, output_pdf) {
  # Using wkhtmltopdf as an example (install via your package manager if missing)
  system(paste("wkhtmltopdf", html_path, output_pdf))
}

# Main function
main <- function() {
  parser <- ArgumentParser(description = "Generate a report with PDF, CSV, PNG, and HTML inputs.")
  parser$add_argument("--sample_id", type = "character", required = TRUE, help = "Sample ID")
  parser$add_argument("--html_cnv", type = "character", required = TRUE, help = "Path to the HTML file for CNV")
  parser$add_argument("--cnv_csv", type = "character", required = TRUE, help = "Path to the CNV CSV file")
  parser$add_argument("--met_csv", type = "character", required = TRUE, help = "Path to the Methylation CSV file")
  parser$add_argument("--merge_annotation", type = "character", required = TRUE, help = "Path to the merged annotation CSV file")
  parser$add_argument("--html_svanna", type = "character", required = TRUE, help = "Path to the HTML file for Svanna")
  parser$add_argument("--csv_annotsv", type = "character", required = TRUE, help = "Path to the AnnotSV CSV file")
  parser$add_argument("--html_terp", type = "character", required = TRUE, help = "Path to the HTML file for TERTp mutations")
  parser$add_argument("--svanna_html", type = "character", required = TRUE, help = "Path to the Svanna HTML file")
  parser$add_argument("--annotsv_html", type = "character", required = TRUE, help = "Path to the AnnotSV HTML file")
  parser$add_argument("--output_dir", type = "character", required = TRUE, help = "Output directory for reports")

  args <- parser$parse_args()

  if (!dir.exists(args$output_dir)) {
    dir.create(args$output_dir, recursive = TRUE)
  }

  output_md <- file.path(args$output_dir, "report.md")
  generate_markdown_report(
    args$sample_id, readLines(args$html_cnv), args$cnv_csv, args$met_csv,
    args$merge_annotation, readLines(args$html_svanna), args$csv_annotsv,
    readLines(args$html_terp), readLines(args$svanna_html),
    readLines(args$annotsv_html), output_md
  )

  output_html <- file.path(args$output_dir, "report.html")
  convert_markdown_to_html(output_md, output_html)

  output_pdf <- file.path(args$output_dir, "report.pdf")
  convert_html_to_pdf(output_html, output_pdf)

  cat("Report generated:\n", output_md, "\n", output_html, "\n", output_pdf, "\n")
}

if (!interactive()) {
  main()
}
