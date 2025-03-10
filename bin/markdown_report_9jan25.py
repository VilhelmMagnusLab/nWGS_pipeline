#!/usr/bin/env python3

import argparse
import pandas as pd
from pdf2image import convert_from_path
import markdown2
from bs4 import BeautifulSoup
import os
from tabulate import tabulate
from PIL import Image
from datetime import datetime
import sys
sys.setrecursionlimit(10000)  # Increase recursion limit

def convert_pdf_to_images(pdf_path, output_dir, max_width=1600, max_height=400):
    """Convert a PDF file to images, resize them, and save them."""
    images = convert_from_path(pdf_path)
    image_paths = []
    for i, image in enumerate(images):
        image_path = os.path.join(output_dir, f"pdf_image_{i + 1}.png")
        resized_image = image.resize((max_width, max_height), Image.Resampling.LANCZOS)
        resized_image.save(image_path, "PNG")
        image_paths.append(image_path)
    return image_paths

def generate_markdown_report(sample_id, html_cnv, csv_cnv, csv_methyl_path, csv_merge_annotation, html_svanna, csv_annotsv_path, html_terp, svanna_html, annotsv_html, output_md):
    """Generate a Markdown report combining various inputs."""
    today_date = datetime.now().strftime("%Y-%m-%d")
    markdown_content = f"""
<div style="text-align: center; font-size: 34px;">
    <u><strong> {sample_id}: Initial WGS Report - P24 Nanopore Sequencing</strong></u>
</div>
<span style="font-size: 18px;"><strong>Sample ID</strong></span>: <span style="font-size: 18px; color: red;"><strong> {sample_id}</strong></span>
<span style="font-size: 18px;"><strong>Whole genome coverage</strong></span>: <span style="font-size: 18px;color: red;"><strong> 47.7x </strong></span>
<span style="font-size: 18px;"><strong>MGMT promoter methylation</strong></span>: <span style="font-size: 18px;color: red;"><strong> Unmethylated, low confidence (grey-zone) </strong></span>
<hr style="border: 1px solid black;"/>
"""
    

    #nanodx

    markdown_content += '<h2><span style="color:blue;">NanoDX</span></h2>\n'
    df = pd.read_csv(nanodx, delimiter="\t")
    df2 = pd.read_csv(nanodict, delimiter="\t")
    

    # Copy Number Variation Section
    markdown_content += '<h2><span style="color:blue;">Copy Number Variation</span></h2>\n'
    markdown_content += f"<div>{html_cnv}</div>\n"

    # Read CNV CSV file and truncate if necessary
    markdown_content += '<h2><span style="color:blue;">Copy Number Variation Filter Table</span></h2>\n'
    df = pd.read_csv(csv_cnv, delimiter=",")
    markdown_content += df.head(50).to_markdown(index=False) + "\n\n"  # Limit rows to avoid large output

    # Descriptive text for CNV
    max_row = df.iloc[df['LOG2CNT'].idxmax()]
    min_row = df.iloc[df['LOG2CNT'].idxmin()]
    sorted_df = df.sort_values(by="LOG2CNT")
    second_min_row = sorted_df.iloc[1]
    description_text = f"""
<span style="font-size: 20px;">
    Gain on <strong style="color:red;">{max_row['Chrom']}</strong>, 
    partial loss on <strong style="color:red;">{min_row['Chrom']}q</strong> 
    and <strong style="color:red;">{second_min_row['Chrom']}q</strong>
</span>
"""
    markdown_content += f"<p>{description_text}</p>\n"

    # Methylation Section
    markdown_content += '<h2><span style="color:blue;">Methylation</span></h2>\n'
    df = pd.read_csv(csv_methyl_path)
    markdown_content += df.head(50).to_markdown(index=False) + "\n\n"  # Limit rows

    # SNV Calling and Annotation Section
    markdown_content += '<h2><span style="color:blue;">SNV calling and annotation</span></h2>\n'
    df = pd.read_csv(csv_merge_annotation, delimiter='\t', quotechar='"')
    markdown_content += df.head(50).to_markdown(index=False) + "\n\n"  # Limit rows

    # Structure Variants Report
    markdown_content += '<h2><span style="color:blue;">Svanna Structure Variants</span></h2>\n'
    markdown_content += f"<div>{html_svanna}</div>\n"

    # AnnotSV Structure Variants
    df = pd.read_csv(csv_annotsv_path, delimiter='\t')
    html_table = df.head(50).to_html(index=False)  # Limit rows
    markdown_content += '<h2><span style="color:blue;">AnnotSV Structure Variants</span></h2>\n'
    markdown_content += f'<div>{html_table}</div>\n\n'

    # Include HTML content for TERTp mutations
    markdown_content += '<h2><span style="color:blue;">TERTp mutations</span></h2>\n'
    markdown_content += f"<div>{html_terp}</div>\n"

    # Logo and Signatures
    markdown_content += f"""
<hr style="border: 1px solid black; margin-top: 40px;">
<div style="margin-top: 20px; font-size: 18px;">
    <span style="font-size: 18px; color: red;"><strong> Vilhelm Magnus Laboratory for Neurosurgical Research</strong></span> <br>
    <strong>Date:</strong> {today_date}
</div>
"""

    # Save to markdown file
    with open(output_md, 'w') as f:
        f.write(markdown_content)

def convert_markdown_to_html(md_path, output_html):
    """Convert Markdown file to HTML using markdown2."""
    with open(md_path, 'r') as f:
        markdown_text = f.read()
    html = markdown2.markdown(markdown_text, extras=["tables"])
    with open(output_html, 'w') as f:
        f.write(html)

def convert_html_to_pdf(html_path, output_pdf):
    """Convert HTML to PDF using weasyprint."""
    from weasyprint import HTML
    try:
        HTML(html_path).write_pdf(output_pdf)
    except RecursionError as e:
        print(f"Error: {e}. Check the HTML structure or increase recursion limit.")
        raise

def main():
    parser = argparse.ArgumentParser(description="Generate a report with PDF, CSV, PNG, and HTML inputs.")
    parser.add_argument("--sample_id", type=str, required=True, help="Sample ID")
    parser.add_argument("--html_cnv", type=str, required=True, help="Path to the HTML file for CNV")
    parser.add_argument("--cnv_csv", type=str, required=True, help="Path to the CNV CSV file")
    parser.add_argument("--met_csv", type=str, required=True, help="Path to the Methylation CSV file")
    parser.add_argument("--merge_annotation", type=str, required=True, help="Path to the merged annotation CSV file")
    parser.add_argument("--html_svanna", type=str, required=True, help="Path to the HTML file for Svanna")
    parser.add_argument("--csv_annotsv", type=str, required=True, help="Path to the AnnotSV CSV file")
    parser.add_argument("--html_terp", type=str, required=True, help="Path to the HTML file for TERTp mutations")
    parser.add_argument("--svanna_html", type=str, required=True, help="Path to the Svanna HTML file")
    parser.add_argument("--annotsv_html", type=str, required=True, help="Path to the AnnotSV HTML file")
    parser.add_argument("--output_dir", type=str, required=True, help="Output directory for reports")

    args = parser.parse_args()

    # Read HTML content
    with open(args.html_cnv, 'r') as f:
        html_cnv = f.read()
    with open(args.html_svanna, 'r') as f:
        html_svanna = f.read()
    with open(args.html_terp, 'r') as f:
        html_terp = f.read()
    with open(args.svanna_html, 'r') as f:
        svanna_html = f.read()
    with open(args.annotsv_html, 'r') as f:
        annotsv_html = f.read()

    # Generate the Markdown report
    output_md = os.path.join(args.output_dir, "report.md")
    generate_markdown_report(args.sample_id, html_cnv, args.cnv_csv, args.met_csv, args.merge_annotation, html_svanna, args.csv_annotsv, html_terp, svanna_html, annotsv_html, output_md)

    # Convert the Markdown to HTML
    output_html = os.path.join(args.output_dir, "report.html")
    convert_markdown_to_html(output_md, output_html)

    # Convert the HTML to PDF
    output_pdf = os.path.join(args.output_dir, "report.pdf")
    convert_html_to_pdf(output_html, output_pdf)

    print(f"Report generated: {output_md}, {output_html}, {output_pdf}")

if __name__ == "__main__":
    main()
