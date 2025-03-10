#!/usr/bin/env python3

import argparse
import pandas as pd
from pdf2image import convert_from_path
import markdown2
import os
from tabulate import tabulate
from PIL import Image
import shutil

def copy_to_assets(file_path, assets_dir):
    """Copy a file to the assets directory and return its relative path."""
    if os.path.exists(file_path):
        destination = os.path.join(assets_dir, os.path.basename(file_path))
        shutil.copy(file_path, destination)
        return destination  # Return full path within assets for reading
    else:
        raise FileNotFoundError(f"The file {file_path} does not exist.")

def generate_markdown_report(sample_id, html_cnv_path, csv_cnv_path, csv_methyl_path, csv_merge_annotation_path, html_svanna_path, csv_annovafusion_path, html_content_path, svanna_html_path, annotsv_html_path, output_md, report_dir):
    """Generate a Markdown report combining PDF images, a CSV table, PNG images, and HTML content."""
    
    # Create assets directory within the report folder
    assets_dir = os.path.join(report_dir, "assets")
    os.makedirs(assets_dir, exist_ok=True)  # Ensure assets directory is created
    
    # Copy all input files to the assets directory for portability
    html_cnv = copy_to_assets(html_cnv_path, assets_dir)
    csv_cnv = copy_to_assets(csv_cnv_path, assets_dir)
    csv_methyl = copy_to_assets(csv_methyl_path, assets_dir)
    csv_merge_annotation = copy_to_assets(csv_merge_annotation_path, assets_dir)
    html_svanna = copy_to_assets(html_svanna_path, assets_dir)
    csv_annovafusion = copy_to_assets(csv_annovafusion_path, assets_dir)
    svanna_html = copy_to_assets(svanna_html_path, assets_dir)
    annotsv_html = copy_to_assets(annotsv_html_path, assets_dir)
    html_content = copy_to_assets(html_content_path, assets_dir)

    # Generate the Markdown content
    markdown_content = f"""
<div style="text-align: center; font-size: 34px; margin-top: 20px;">
    <u><strong>{sample_id}: Initial WGS Report - P24 Nanopore Sequencing</strong></u>
</div>

<span style="font-size: 18px;"><strong>Sample ID:</strong> <span style="color: red;">{sample_id}</span></span><br>
<span style="font-size: 18px;"><strong>Whole genome coverage:</strong> <span style="color: red;">47.7x</span></span><br>
<span style="font-size: 18px;"><strong>MGMT promoter methylation:</strong> <span style="color: red;">Unmethylated, low confidence (grey-zone)</span></span><br>

<hr style="border: 1px solid black; margin-top: 20px;"/>
"""

    # Copy Number Variation Section
 ###   markdown_content += '<h2 style="color:blue; margin-top: 20px;">Copy Number Variation</h2>\n'
 ###   with open(html_cnv, 'r') as f:
 ###       html_cnv_content = f.read()
 ###   markdown_content += f"<div>{html_cnv_content}</div>\n\n"


    markdown_content += """
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Copy Number Variation Report</title>
    <script src="./assets/lib/htmlwidgets-1.5.4/htmlwidgets.js"></script>
    <script src="./assets/lib/plotly-binding-4.10.4/plotly.js"></script>
    <link href="./assets/lib/plotly-htmlwidgets-css-2.11.1/plotly-htmlwidgets.css" rel="stylesheet">
</head>
<body>
"""

# Add the section title and embed `html_cnv_content` directly in the HTML
    markdown_content += '<h2 style="color:blue; margin-top: 20px;">Copy Number Variation</h2>\n'
    with open(html_cnv, 'r') as f:
        html_cnv_content = f.read()

# Embed the content in a <div> for proper display
    markdown_content += f"""
<div style="border: 1px solid #ddd; padding: 10px; margin-top: 20px;">
    {html_cnv_content}
</div>
"""










    # Copy Number Variation Filter Table
    df_cnv = pd.read_csv(csv_cnv)
    markdown_content += '<h3 style="color:blue; margin-top: 20px;">Copy Number Variation Filter Table</h3>\n'
    markdown_content += df_cnv.to_markdown(index=False, tablefmt="pipe") + "\n\n"

    # Additional descriptive text based on data in CNV table
    max_row = df_cnv.iloc[df_cnv['LOG2CNT'].idxmax()]
    min_row = df_cnv.iloc[df_cnv['LOG2CNT'].idxmin()]
    sorted_df = df_cnv.sort_values(by="LOG2CNT")
    second_min_row = sorted_df.iloc[1]
    description_text = f"""
<span style="font-size: 16px; color: black;">
    Gain on <strong style="color:red;">{max_row['Chrom']}</strong>, 
    partial loss on <strong style="color:red;">{min_row['Chrom']}q</strong>, 
    and <strong style="color:red;">{second_min_row['Chrom']}q</strong>.
</span>
"""
    markdown_content += f"<p>{description_text}</p>\n"

    # Methylation Section
    df_methyl = pd.read_csv(csv_methyl)
    markdown_content += '<h2 style="color:blue; margin-top: 20px;">Methylation</h2>\n'
    markdown_content += df_methyl.to_markdown(index=False, tablefmt="pipe") + "\n\n"

    # SNV Calling and Annotation Section
    df_annotation = pd.read_csv(csv_merge_annotation, delimiter='\t', quotechar='"')
    def add_cosmic_link(value):
        if pd.notna(value) and "ID=" in value:
            cosmic_id = value.split("ID=")[-1].strip()
            return f'<a href="https://cancer.sanger.ac.uk/cosmic/search?q={cosmic_id}" target="_blank">{value}</a>'
        return value
    if 'COSMIC100' in df_annotation.columns:
        df_annotation['COSMIC100'] = df_annotation['COSMIC100'].apply(add_cosmic_link)
    markdown_content += '<h2 style="color:blue; margin-top: 20px;">SNV Calling and Annotation</h2>\n'
    markdown_content += df_annotation.to_markdown(index=False, tablefmt="pipe") + "\n\n"


    ###
    markdown_content += '<h2><span style="color:blue;">Svanna Structure Variants</span></h2>\n'
    with open(html_svanna, 'r') as f:
        html_svanna_content = f.read()
    markdown_content += f"<div>{html_svanna_content}</div>\n\n"



    # AnnotSV Structure Variants Section
  ###  df_annotsv = pd.read_csv(csv_annovafusion)
  ###  markdown_content += '<h2 style="color:blue; margin-top: 20px;">AnnotSV Structure Variants</h2>\n'
  ###  markdown_content += df_annotsv.to_markdown(index=False, tablefmt="pipe") + "\n\n"

    
    df = pd.read_csv(csv_annovafusion,  delimiter='\t')
    html_table = df.to_html(index=False)
    markdown_content += '<h2><span style="color:blue;">AnnotSV Structure Variants</span></h2>\n'
    markdown_content += f'<div>{html_table}</div>\n\n'
    markdown_content += '<h2><span ;"> Structure Variants Further results</span></h2>\n'
    

    # Links to Full HTML Reports
    markdown_content += f"""
<h3 style="margin-top: 20px;">Structure Variants Further Results</h3>
<p><span style="font-size: 15px;">
    The AnnotSV full result can be found in: <a href="assets/{os.path.basename(annotsv_html)}" target="_blank">AnnotSV Results</a><br>
    The Svanna full result can be found in: <a href="assets/{os.path.basename(svanna_html)}" target="_blank">Svanna Results</a>
</span></p>
"""

    # Additional HTML Content
    markdown_content += '<h2 style="color:blue; margin-top: 20px;">TERTp mutations</h2>\n'
    with open(html_content, 'r') as f:
        extra_content = f.read()
    markdown_content += f"<div>{extra_content}</div>\n"

    # Save the markdown content to a file
    with open(output_md, 'w') as f:
        f.write(markdown_content)

def convert_markdown_to_html(md_path, output_html):
    """Convert Markdown file to HTML using markdown2."""
    with open(md_path, 'r') as f:
        markdown_text = f.read()
    html = markdown2.markdown(markdown_text, extras=["tables"])
    with open(output_html, 'w') as f:
        f.write(html)

def main():
    parser = argparse.ArgumentParser(description="Generate a report with all dependencies in a portable folder.")
    parser.add_argument("--sample_id", type=str, help="Sample ID for the report", required=True)
    parser.add_argument("--html_cnv", type=str, help="Path to the HTML file for CNV", required=True)
    parser.add_argument("--cnv_csv", type=str, help="Path to the CSV file for CNV", required=True)
    parser.add_argument("--met_csv", type=str, help="Path to the CSV file for methylation", required=True)
    parser.add_argument("--merge_annotation", type=str, help="Path to the CSV file for SNV annotation", required=True)
    parser.add_argument("--html_svanna", type=str, help="Path to the Svanna HTML file", required=True)
    parser.add_argument("--csv_annotsv", type=str, help="Path to the AnnotSV CSV file", required=True)
    parser.add_argument("--html", type=str, help="Path to additional HTML content", required=True)
    parser.add_argument("--svanna_html", type=str, help="Path to the full Svanna HTML report", required=True)
    parser.add_argument("--annotsv_html", type=str, help="Path to the full AnnotSV HTML report", required=True)
    parser.add_argument("--output_dir", type=str, help="Output directory for reports", required=True)
    args = parser.parse_args()

    # Create the report directory and assets subdirectory
    report_dir = os.path.join(args.output_dir, "report")
    assets_dir = os.path.join(report_dir, "assets")
    os.makedirs(report_dir, exist_ok=True)
    os.makedirs(assets_dir, exist_ok=True)
    
    # Generate the Markdown report
    output_md = os.path.join(report_dir, "report.md")
    generate_markdown_report(args.sample_id, args.html_cnv, args.cnv_csv, args.met_csv, args.merge_annotation, args.html_svanna, args.csv_annotsv, args.html, args.svanna_html, args.annotsv_html, output_md, report_dir)

    # Convert the Markdown to HTML
    output_html = os.path.join(report_dir, "report.html")
    convert_markdown_to_html(output_md, output_html)

    print(f"Report generated at {output_md} and {output_html}")

if __name__ == "__main__":
    main()
