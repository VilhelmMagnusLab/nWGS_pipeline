#!/usr/bin/env python3
import sys
import csv
import re

def main(file1_path, file2_path, output_path, methylation_output_path, filtered_output_path):
    # Step 1: Compile the regex pattern for "cg" followed by exactly eight digits
    pattern = re.compile(r"^cg\d{8}$")

    # Step 2: Read the third column from file1 and store entries that match the "cgxxxxxxxx" pattern
    pattern_elements = set()
    with open(file1_path, 'r') as file1:
        reader = csv.reader(file1, delimiter='\t')  # Assuming tab-separated files
        for row in reader:
            if len(row) >= 4 and pattern.match(row[3]):  # Column 4 (index 3) matches "cgxxxxxxxx"
                pattern_elements.add(row[3])

    # Step 3: Filter rows from file2 where the seventh column matches the "cgxxxxxxxx" pattern from file1
    matching_rows = []
    with open(file2_path, 'r') as file2:
        reader = csv.reader(file2, delimiter='\t')  # Assuming tab-separated files
        header = next(reader)  # Read the header line
        for row in reader:
            if len(row) >= 7 and pattern.match(row[6]) and row[6] in pattern_elements:  # Column 7 (index 6) matches "cgxxxxxxxx" and is in patterns
                matching_rows.append(row)  # Store the entire row

    # Step 4: Write the header and matching rows to the output file
    with open(output_path, 'w', newline='') as output_file:
        writer = csv.writer(output_file, delimiter='\t')
        writer.writerow(header)  # Write the header
        writer.writerows(matching_rows)  # Write matching rows

    # Step 5: Generate the methylation output file
    generate_methylation_file(output_path, methylation_output_path)

    # Step 6: Filter the methylation output file to remove rows with methylation_call -1
    filter_methylation_file(methylation_output_path, filtered_output_path)

def generate_methylation_file(input_file, output_file):
    # Read the input file and create the methylation output file with required criteria
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.DictReader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        # Write the header for the methylation output file
        writer.writerow(["index", "probe_id", "methylation_call"])

        # Process each row according to the methylation call criteria
        for index, row in enumerate(reader):
            probe_id = row["Illumina_ID"]
            methylation_frequency = float(row["Methylation_frequency"])

            # Apply the methylation call criteria
            if methylation_frequency < 60:
                methylation_call = -1.0
            else:
                methylation_call = 1.0

            # Write the row to the methylation output file
            writer.writerow([index, probe_id, methylation_call])

def filter_methylation_file(input_file, output_file):
    # Filter the methylation output file to only include rows where methylation_call is 1
    with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
        reader = csv.reader(infile, delimiter='\t')
        writer = csv.writer(outfile, delimiter='\t')

        # Write the header to the filtered output file
        header = next(reader)
        writer.writerow(header)

        # Write only rows with methylation_call == 1
        new_index=0
        for row in reader:
            if row[2] == "1.0":  # Column 3 (index 2) is methylation_call
                row[0] =new_index
                writer.writerow(row)
                new_index +=1

# Run the script with file arguments
if __name__ == "__main__":
    if len(sys.argv) != 6:
        print("Usage: python script.py <bedfile> <450K_hg19bedfile> <output> <methylation_output> <filtered_output>")
        sys.exit(1)

    file1_path = sys.argv[1]
    file2_path = sys.argv[2]
    output_path = sys.argv[3]
    methylation_output_path = sys.argv[4]
    filtered_output_path = sys.argv[5]
    main(file1_path, file2_path, output_path, methylation_output_path, filtered_output_path)

