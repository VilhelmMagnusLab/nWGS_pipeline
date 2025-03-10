
#!/usr/bin/env python3

import argparse
import pandas as pd
from pdf2image import convert_from_path
import markdown2
import os

def convert_pdf_to_images(pdf_path, output_dir):
    """Convert a PDF file to images and save them."""
    images = convert_from_path(pdf_path)
    image_paths = []
    for i, image in enumerate(images):
        image_path = os.path.join(output_dir, f"pdf_image_{i + 1}.png")
        image.save(image_path, "PNG")
        image_paths.append(image_path)
    return image_paths

def generate_markdown_report(pdf_images, csv_table_path, png_images, html_content, output_md):
    """Generate a Markdown report combining PDF images, a CSV table, PNG images, and HTML content."""
    markdown_content = """
<div style="text-align: center; font-size: 34px;">
    <u><strong>Initial WGS Report - P24 Nanopore Sequencing</strong></u>
</div>

<span style="font-size: 18px;"><strong>Sample ID</strong></span>: T24-146

<span style="font-size: 18px;"><strong>Whole genome coverage</strong></span>: 47.7x

<span style="font-size: 18px;"><strong>MGMT promoter methylation</strong></span>: Unmethylated, low confidence (grey-zone)

<!-- Bold horizontal line to separate this section -->
<hr style="border: 1px solid black;"/>
"""


    # Embed PDF images
    markdown_content += '<h2><span style="color:blue;">Copy Number Variation</span></h2>\n'
    for img in pdf_images:
        markdown_content += f"![PDF Image]({img})\n\n"

    # Embed CSV table as Markdown table
    markdown_content += '<h2><span style="color:blue;">SNV calling and annotation</span></h2>\n'
    df = pd.read_csv(csv_table_path)
    markdown_content += df.to_markdown(index=False) + "\n\n"

    # Embed PNG images
    markdown_content += '<h2><span style="color:blue;">Annotation Circos plot</span></h2>\n'
    for img in png_images:
        markdown_content += f"![PNG Image]({img})\n\n"

    # Include HTML content
    markdown_content += '<h2><span style="color:blue;">TERTp mutations</span></h2>\n'
    markdown_content += f"<div>{html_content}</div>\n"

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
    HTML(html_path).write_pdf(output_pdf)

def main():
    parser = argparse.ArgumentParser(description="Generate a report with PDF, CSV, PNG, and HTML inputs.")
    parser.add_argument("--pdf", type=str, help="Path to the PDF file containing an image", required=True)
    parser.add_argument("--csv", type=str, help="Path to the CSV file", required=True)
    parser.add_argument("--pngs", nargs='+', help="Paths to the PNG images", required=True)
    parser.add_argument("--html", type=str, help="Path to the HTML file", required=True)
    parser.add_argument("--output_dir", type=str, help="Output directory for reports", required=True)
    args = parser.parse_args()

    # Read the HTML content
    with open(args.html, 'r') as f:
        html_content = f.read()

    # Convert the PDF file to images
    pdf_images = convert_pdf_to_images(args.pdf, args.output_dir)

    # Generate the Markdown report
    output_md = os.path.join(args.output_dir, "report.md")
    generate_markdown_report(pdf_images, args.csv, args.pngs, html_content, output_md)

    # Convert the Markdown to HTML
    output_html = os.path.join(args.output_dir, "report.html")
    convert_markdown_to_html(output_md, output_html)

    # Convert the HTML to PDF
    output_pdf = os.path.join(args.output_dir, "report.pdf")
    convert_html_to_pdf(output_html, output_pdf)

    print(f"Report generated: {output_md}, {output_html}, {output_pdf}")

if __name__ == "__main__":
    main()
