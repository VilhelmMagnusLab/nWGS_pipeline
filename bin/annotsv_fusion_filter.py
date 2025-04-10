#!/usr/bin/env python3
import csv
import sys
import os

# Increase CSV field size limit
maxInt = sys.maxsize
while True:
    try:
        csv.field_size_limit(maxInt)
        break
    except OverflowError:
        maxInt = int(maxInt/10)

def main():
    if len(sys.argv) != 4:
        print("Usage: annotsv_fusion_filter.py <input_file> <fusion_genes_list> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    fusion_genes_list = sys.argv[2]
    output_file = sys.argv[3]

    # Read fusion genes list
    with open(fusion_genes_list, 'r') as f:
        fusion_genes = set(line.strip() for line in f)

    # Process AnnotSV file
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        # Write header
        header = next(reader)
        writer.writerow(header)

        # Process rows
        for row in reader:
            try:
                gene_name = row[header.index('Gene_name')]
                if any(gene in fusion_genes for gene in gene_name.split('/')):
                    writer.writerow(row)
            except (ValueError, IndexError) as e:
                print(f"Warning: Error processing row: {e}", file=sys.stderr)
                continue

if __name__ == "__main__":
    main()

