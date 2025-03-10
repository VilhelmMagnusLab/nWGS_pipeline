#!/usr/bin/env python3

import csv
import sys

def generate_methylation_file(input_file, output_file):
    print(f"Processing input file: {input_file}")
    print(f"Writing to output file: {output_file}")
    
    try:
        with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
            reader = csv.DictReader(infile, delimiter='\t')
            writer = csv.writer(outfile, delimiter='\t')
            
            # Write header
            writer.writerow(["index", "probe_id", "methylation_call"])
            
            # Process rows
            for index, row in enumerate(reader):
                probe_id = row.get("IlmnID")
                methylation_frequency = row.get("methylated_frequency")
                
                try:
                    methylation_frequency = float(methylation_frequency) if methylation_frequency != 'NA' else 0.0
                    methylation_call = -1.0 if methylation_frequency < 60 else 1.0
                except (ValueError, TypeError):
                    print(f"Warning: Invalid methylation value in row {index}")
                    methylation_call = -1.0
                
                writer.writerow([index, probe_id, methylation_call])

    except FileNotFoundError as e:
        print(f"Error: Could not find file - {e}")
        sys.exit(1)
    except Exception as e:
        print(f"Error processing files: {e}")
        sys.exit(1)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_file> <output_file>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    generate_methylation_file(input_file, output_file)