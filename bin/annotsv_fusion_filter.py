#!/usr/bin/env python
import csv
import re
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description="Extract matching rows and specific columns from File 1 based on gene patterns in File 2.")
parser.add_argument('file1', type=str, help='Path to File 1 (annotated variants file)')
parser.add_argument('file2', type=str, help='Path to File 2 (gene names file)')
parser.add_argument('output_file', type=str, help='Path to output file for saving matching rows')

# Parse arguments
args = parser.parse_args()

# File paths
file1_path = args.file1
file2_path = args.file2
output_file_path = args.output_file

# Load gene names from File 2
with open(file2_path, 'r') as gene_file:
    gene_reader = csv.reader(gene_file)
    next(gene_reader, None)  # Skip the header if present
    gene_names = {row[0].strip() for row in gene_reader if row}  # Set of gene names for fast lookup

# Columns to extract from File 1 (zero-based indices for specific columns)
columns_to_extract = [1, 2, 3, 4, 5, 6, 8, 10, 11, 13, 14, 15, 16, 114]

# Read File 1 and find matching rows
matching_rows = []
with open(file1_path, 'r') as file1:
    reader = csv.reader(file1, delimiter='\t', quoting=csv.QUOTE_NONE)
    headers = next(reader)

    # Find the index of the "AnnotSV_ranking_criteria" column
    try:
        criteria_column_index = headers.index("AnnotSV_ranking_criteria")
    except ValueError:
        raise Exception("The header 'AnnotSV_ranking_criteria' was not found in the file.")

    # Add headers for extracted columns to output
    extracted_headers = [headers[i] for i in columns_to_extract]
    matching_rows.append(extracted_headers)

    # Check each row for matching patterns
    for row in reader:
        if len(row) > criteria_column_index:
            criteria_column = row[criteria_column_index]

            # Extract and retain only matching gene patterns in the criteria column
            matched_genes = []
            for gene_name in gene_names:
                pattern = rf'\b{gene_name}/\w+\b'
                matches = re.findall(pattern, criteria_column)
                if matches:
                    matched_genes.extend(matches)

            # If there are matched genes, add them to the extracted data
            if matched_genes:
                extracted_data = [row[i] for i in columns_to_extract]
                # Replace the last item in extracted_data with the matched genes, joined by commas
                extracted_data[-1] = ','.join(matched_genes)
                matching_rows.append(extracted_data)

# Save the matching rows and extracted columns to a new file
with open(output_file_path, 'w', newline='') as output_file:
    writer = csv.writer(output_file, delimiter='\t')
    writer.writerows(matching_rows)

