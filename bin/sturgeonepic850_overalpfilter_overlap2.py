#!/usr/bin/env python3
import csv
import re
import sys

def filter_methylation_data(file1, file2, output_file):
    # Read the valid IDs from file1
    valid_ids = set()
    with open(file1, 'r') as f1:
        for line in f1:
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                valid_ids.add(parts[3])  # cgxxxx pattern is in column 4

    # Compile regex pattern for cgxxxx IDs
    pattern = re.compile(r'^cg\d+$')

    # Open file2 and filter based on criteria
    with open(file2, 'r') as f2, open(output_file, 'w', newline='') as out:
    # Preprocess the input file to replace spaces with tabs
        lines = f2.readlines()
        processed_lines = [re.sub(r'\s+', '\t', line.strip()) for line in lines]  # Replace multiple spaces/tabs with a single tab

    # Use csv.DictReader on the processed lines
        reader = csv.DictReader(processed_lines, delimiter='\t')
        writer = csv.writer(out, delimiter='\t')

        # Write the header for the output file
        writer.writerow(["index", "probe_id", "methylation_call"])

        output_index = 0  # Reset index counter for the output file
        for row in reader:
            # Filter out rows with "h" in the "modBase" column
            if 'h' in row.get("modBase", "") or 'h' in row.get("ModBase", ""):
                continue  # Skip this row

            # Check if column 4 matches the pattern and if the probe_id is in valid_ids
            if len(row) >= 6 and pattern.match(row["IlmnID"]) and row["IlmnID"] in valid_ids:
                probe_id = row.get("IlmnID")
                methylation_frequency = row.get("methylated_frequency")

                try:
                    # Convert methylation_frequency to a float
                    methylation_frequency = float(methylation_frequency) if methylation_frequency != 'NA' else 0.0
                    methylation_call = 0 if methylation_frequency < 60 else 1
                except (ValueError, TypeError):
                    print(f"Warning: Invalid methylation value")
                    methylation_call = -1.0

                # Write the filtered and processed row to the output file
                writer.writerow([output_index, probe_id, methylation_call])
                output_index += 1  # Increment the output index

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python script.py <file1> <file2> <output_file>")
        sys.exit(1)

    file1 = sys.argv[1]
    file2 = sys.argv[2]
    output_file = sys.argv[3]

    filter_methylation_data(file1, file2, output_file)