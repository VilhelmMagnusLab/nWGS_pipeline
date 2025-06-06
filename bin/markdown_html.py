#!/usr/bin/env python3

import argparse
import html2text
from PyPDF2 import PdfReader

def generate_markdown(sample_id, html_report, output_file):
    # Read the HTML report content
    with open(html_report, 'r') as f:
        html_content = f.read()

    # Convert the HTML content to Markdown
 ###   markdown_converter = html2text.HTML2Text()
 ###   markdown_converter.ignore_links = False
 ##   markdown_content = markdown_converter.handle(html_content)

    # Create the final Markdown report

def extract_text_from_pdf(pdf_path):
    with open(pdf_path, 'rb') as f:
        reader = PdfReader(pdf_path)
        text = ""
        for page in reader.pages:
            text += page.extract_text()
    
        return text

def generate_markdown(sample_id, pdf_report, output_file):
    # Extract text from the PDF report
    pdf_text = extract_text_from_pdf(pdf_report)

    markdown_content = f"""
# **Initial WGS report, P24 nanopore sequencing**

## Report for **{sample_id}**

## Whole genome coverage: 47.7x

## *MGMT promoter methylation*: Unmethylated, low confidence (grey-zone)

# *Initial nanoDx methylation analysis*:

# *CNV plot*:
{pdf_text}

---

Generated by Python Markdown Generator.
    """
#{html_content}
    # Write the markdown content to the output file
    with open(output_file, 'w') as f:
        f.write(markdown_content)

def main():
    # Initialize the argument parser
    parser = argparse.ArgumentParser(description="Generate a Markdown report from sample data.")

    # Add arguments for sample_id, HTML report, and output file
    parser.add_argument("--s", type=str, help='Sample ID for the report', required=True)
    parser.add_argument("--f", type=str, help='HTML report to be summarized', required=True)
    parser.add_argument("--out", type=str, help='Path to the output markdown file', required=True)

    # Parse the command-line arguments
    args = parser.parse_args()

    # Call the function to generate markdown using parsed arguments
    generate_markdown(args.s, args.f, args.out)

if __name__ == "__main__":
    main()
