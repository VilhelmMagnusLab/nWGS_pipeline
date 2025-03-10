
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

###def convert_pdf_to_images(pdf_path, output_dir):
#    """Convert a PDF file to images and save them."""
##    images = convert_from_path(pdf_path)
#    image_paths = []
#    for i, image in enumerate(images):
#        image_path = os.path.join(output_dir, f"pdf_image_{i + 1}.png")
#        image.save(image_path, "PNG")
#        image_paths.append(image_path)
#    return image_paths


def convert_pdf_to_images(pdf_path, output_dir, max_width=1600, max_height=400):
    """Convert a PDF file to images, resize them, and save them."""
    images = convert_from_path(pdf_path)
    image_paths = []
    for i, image in enumerate(images):
        image_path = os.path.join(output_dir, f"pdf_image_{i + 1}.png")

        # Resize the image using Image.Resampling.LANCZOS
        resized_image = image.resize((max_width, max_height), Image.Resampling.LANCZOS)
        resized_image.save(image_path, "PNG")
        image_paths.append(image_path)
    return image_paths




##def generate_markdown_report(sample_id, pdf_images, csv_cnv_path, csv_methyl_path, csv_table_path, html_svanna, csv_annovafusion_path, html_annotsv, html_content, output_md):
###def generate_markdown_report(sample_id, html_cnv, pdf_images, csv_cnv, csv_methyl_path, csv_merge_annotation, html_svanna, csv_annovafusion_path, html_content, svanna_html, annotsv_html, output_md):
def generate_markdown_report(sample_id, html_cnv, csv_cnv, csv_methyl_path, csv_merge_annotation, html_svanna, csv_annovafusion_path, html_content, svanna_html, annotsv_html, output_md):

    today_date = datetime.now().strftime("%Y-%m-%d")

    """Generate a Markdown report combining PDF images, a CSV table, PNG images, and HTML content."""
    markdown_content = f"""
<div style="text-align: center; font-size: 34px;">
    <u><strong> {sample_id}:Initial WGS Report - P24 Nanopore Sequencing</strong></u>
</div>

<span style="font-size: 18px;"><strong>Sample ID</strong></span>: <span style="font-size: 18px; color: red;"><strong> {sample_id}</strong></span> 


<span style="font-size: 18px;"><strong>Whole genome coverage</strong></span>: <span style="font-size: 18px;color: red;"><strong> 47.7x </strong></span>

<span style="font-size: 18px;"><strong>MGMT promoter methylation</strong></span>:  <span style="font-size: 18px;color: red;"><strong> Unmethylated, low confidence (grey-zone) </strong></span> 

<!-- Bold horizontal line to separate this section -->
<hr style="border: 1px solid black;"/>
"""

    #Copy number variation

    #markdwon cnv html file

    markdown_content += '<h2><span style="color:blue;">Copy Number Variation</span></h2>\n'
    markdown_content += f"<div>{html_cnv}</div>\n"

    # Embed PDF images 
 ###   markdown_content += '<h2><span style="color:blue;">Copy Number Variation</span></h2>\n'
 ###   for img in pdf_images:
 ###      markdown_content += f"![PDF Image]({img})\n\n"
    
     
    #read cnv csv file
    markdown_content += '<h2><span style="color:blue;">Copy Number Variation Fliter Table</span></h2>\n'
    df = pd.read_csv(csv_cnv, delimiter=",")
    print(df)
    
    markdown_content += df.to_markdown(index=False) + "\n\n"

    max_row = df.iloc[df['LOG2CNT'].idxmax()]
    min_row = df.iloc[df['LOG2CNT'].idxmin()]

# Sort the DataFrame by `LOG2CNT` to get the second-lowest
    sorted_df = df.sort_values(by="LOG2CNT")
    second_min_row = sorted_df.iloc[1]

# Format the descriptive text with identified chromosome values
    #description_text = f"Gain on ({max_row['Chrom']}), partial loss on ({min_row['Chrom']})q and chr({second_min_row['Chrom']})q"
    description_text = f"""
<span style="font-size: 20px;">
    Gain on <strong style="color:red;">{max_row['Chrom']}</strong>, 
    partial loss on <strong style="color:red;">{min_row['Chrom']}q</strong> 
    and <strong style="color:red;">{second_min_row['Chrom']}q</strong>
</span>
"""
# Add the descriptive text after the table in markdown content
    markdown_content += f"<p>{description_text}</p>\n"
    #Methylation

    markdown_content += '<h2><span style="color:blue;">Methylation</span></h2>\n'
    df = pd.read_csv(csv_methyl_path)
    markdown_content += df.to_markdown(index=False) + "\n\n"

    # Embed CSV table as Markdown table
###    markdown_content += '<h2><span style="color:blue;">SNV calling and annotation</span></h2>\n'
###    df = pd.read_csv(csv_merge_annotation, delimiter='\t', quotechar='"')
###    print(df)
###    markdown_content += df.to_markdown(index=False) + "\n\n"

    markdown_content += '<h2><span style="color:blue;">SNV calling and annotation</span></h2>\n'
    df = pd.read_csv(csv_merge_annotation, delimiter='\t', quotechar='"')
# Check for `ID=` in each row of the `COSMIC100` column and modify it
    def add_cosmic_link(value):
        if pd.notna(value) and "ID=" in value:
        # Extract the ID from the value, assuming the format is "ID=xxxxx"
            cosmic_id = value.split("ID=")[-1].strip()  # Get the ID part after "ID="
            return f'<a href="https://cancer.sanger.ac.uk/cosmic/search?q={cosmic_id}" target="_blank">{value}</a>'
        return value
# Apply the function to the 'COSMIC100' column
    df['COSMIC100'] = df['COSMIC100'].apply(add_cosmic_link)

# Convert the DataFrame to markdown with hyperlinks
    markdown_content += df.to_markdown(index=False) + "\n\n"

    # Structure Variants report

    markdown_content += '<h2><span style="color:blue;">Svanna Structure Variants</span></h2>\n'
    markdown_content += f"<div>{html_svanna}</div>\n"


    #read annotsv filter

    df = pd.read_csv(csv_annovafusion_path, delimiter='\t')
    html_table = df.to_html(index=False)

    markdown_content += '<h2><span style="color:blue;">AnnotSV Structure Variants</span></h2>\n'
    markdown_content += f'<div>{html_table}</div>\n\n'
    markdown_content += '<h2><span ;"> Structure Variants Further results</span></h2>\n'



    markdown_content += f"""
<span style="font-size: 15px;">
    The AnnotSV full result can be found in: <a href="{annotsv_html}" target="_blank">AnnotSV Results</a><br>
    The Svanna full result can be found in: <a href="{svanna_html}" target="_blank">Svanna Results</a>
</span>
"""

    
    # Include HTML content
    markdown_content += '<h2><span style="color:blue;">TERTp mutations</span></h2>\n'
    markdown_content += f"<div>{html_content}</div>\n"


    ###Logo and signatures
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
    HTML(html_path).write_pdf(output_pdf)



def main():
    parser = argparse.ArgumentParser(description="Generate a report with PDF, CSV, CSV, PNG, and HTML inputs.")
    parser.add_argument("--sample_id", type=str, help="Path to the PDF file containing an image", required=True)
    parser.add_argument("--html_cnv", type=str, help="Path to the HTML file", required=True)
###    parser.add_argument("--pdf", type=str, help="Path to the PDF file containing an image", required=True)
    parser.add_argument("--cnv_csv", type=str, help="Path to the CSV file", required=True)
    parser.add_argument("--met_csv", type=str, help="Path to the CSV file", required=True)
    parser.add_argument("--merge_annotation", type=str, help="Path to the CSV file", required=True)
 #   parser.add_argument("--pngs", nargs='+', help="Paths to the PNG images", required=True)
    parser.add_argument("--html_svanna", type=str, help="Path to the HTML file", required=True)
    parser.add_argument("--csv_annotsv", type=str, help="Path to the HTML file", required=True)
 ###   parser.add_argument("--html_annotsv", type=str, help="Path to the HTML file", required=True)
    parser.add_argument("--html", type=str, help="Path to the HTML file", required=True)
    parser.add_argument("--svanna_html", type=str, help="Path to the HTML file", required=True)
    parser.add_argument("--annotsv_html", type=str, help="Path to the HTML file", required=True)
    parser.add_argument("--output_dir", type=str, help="Output directory for reports", required=True)
    
    
    args = parser.parse_args()

    
    #read cnv html file

    with open(args.html_cnv, 'r') as f:
        html_cnv = f.read()
        
    # Read the HTML content circos plot

    with open(args.html_svanna, 'r') as f:
        html_svanna = f.read()

   #### with open(args.html_annotsv, 'r') as f:
   ###     html_annotsv = f.read()

    with open(args.html, 'r') as f:
        html_content = f.read()


    with open(args.svanna_html, 'r') as f:
        svanna_html = f.read()

    with open(args.annotsv_html, 'r') as f:
        annotsv_html = f.read()

    # Convert the PDF file to images
    ###pdf_images = convert_from_path(args.pdf, dpi=100, output_folder=args.output_dir)
    ###pdf_images = convert_pdf_to_images(pdf_images, args.output_dir)
###    pdf_images = convert_pdf_to_images(args.pdf, args.output_dir)

    # Generate the Markdown report
    output_md = os.path.join(args.output_dir, "report.md")
###    generate_markdown_report(args.sample_id, html_cnv, pdf_images, args.cnv_csv, args.met_csv, args.merge_annotation, html_svanna, args.csv_annotsv, html_content, args.svanna_html, args.annotsv_html, output_md)
    generate_markdown_report(args.sample_id, html_cnv, args.cnv_csv, args.met_csv, args.merge_annotation, html_svanna, args.csv_annotsv, html_content, args.svanna_html, args.annotsv_html, output_md)

    # Convert the Markdown to HTML
    output_html = os.path.join(args.output_dir, "report.html")
    convert_markdown_to_html(output_md, output_html)

    # Convert the HTML to PDF
    output_pdf = os.path.join(args.output_dir, "report.pdf")
    convert_html_to_pdf(output_html, output_pdf)

    print(f"Report generated: {output_md}, {output_html}, {output_pdf}")

if __name__ == "__main__":
    main()
