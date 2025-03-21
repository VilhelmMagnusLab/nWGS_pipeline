#!/usr/bin/env python3
import math
import os
import pandas as pd  # Import pandas for handling data files
import argparse  # Import argparse for command line arguments

def classify_copy_number_event(log2_ratio, f_T, C_N=2):
    """
    Classify the copy number event (loss, gain, or no change) based on the log2 ratio.

    Parameters:
        log2_ratio (float): The observed log2 ratio.
        f_T (float): The tumor fraction (e.g., 0.8 for 80% tumor).
        C_N (int, optional): The normal copy number (default is 2 for diploid).

    Returns:
        tuple: (C_T, event_type, change_value)
            C_T (float): Calculated tumor copy number.
            event_type (str): 'Loss', 'Gain', 'No Change', or 'Indeterminated'.
            change_value (float or str): Difference in copy number from C_N, or 'Indeterminated'.
    """
    # Calculate the normal fraction
    f_N = 1 - f_T

    try:
        # Calculate tumor copy number (C_T)
        C_T = (C_N * (2**log2_ratio - f_N)) / f_T
    except ZeroDivisionError:
        C_T = float('nan')

    if C_T < 0 or math.isnan(C_T):
        return "0", "Loss", "2"

    C_T_rounded = round(C_T)  # Round to the nearest integer

    # Determine the event type and the change value
    if C_T_rounded < C_N:
        event_type = "Loss"
        change_value = C_N - C_T_rounded  # Amount of loss
    elif C_T_rounded > C_N:
        event_type = "Gain"
        change_value = C_T_rounded - C_N  # Amount of gain
    else:
        event_type = "No Change"
        change_value = 0

    return C_T_rounded, event_type, change_value

def match_and_create_file_old(file1, file2, matched_file):
    """
    Match rows between file1 and file2 based on columns 1 and 2, and create a new file with matched rows.

    Args:
        file1 (str): Path to the first file.
        file2 (str): Path to the second file.
        matched_file (str): Path to the output file with matched rows.
    """
    df1 = pd.read_csv(file1, sep="\t", header=None, skiprows=1)  # File 1
    #print(df1)
    df2 = pd.read_csv(file2, sep="\t", header=None)  # File 2
    #print(df2)

    # Perform matching based on columns 1 and 2
  #  matched = df1.merge(df2, left_on=[1, 2], right_on=[1, 2])
    matched = df1.merge(df2, left_on=[1, 2], right_on=[1, 2], how='inner', suffixes=('', '_file2'))
    matched = matched.iloc[:, :5]
    #print(matched)

    # Save matched rows to a new file
    matched.to_csv(matched_file, sep="\t", index=False, header=False)

def filter_by_interval(file1, file2, output_file):
    """
    Retain all rows from file1 where the value in the second column of file2 falls within the intervals defined by file1.

    Args:
        file1 (str): Path to the first file.
        file2 (str): Path to the second file.
        output_file (str): Path to the output file with filtered rows.
    """
    df1 = pd.read_csv(file1, sep="\t", header=None, skiprows=1)  # Read file1
    df2 = pd.read_csv(file2, sep="\t", header=None)  # Read file2

    # Extract the second column from file2
    values_to_check = df2[1].values  # Assuming the second column is at index 1

    # Create a mask to filter rows in df1
    mask = df1.apply(lambda row: any(row[1] <= value <= row[2] for value in values_to_check), axis=1)

    # Filter df1 based on the mask
    filtered_df = df1[mask]

    # Save the filtered rows to a new file
    filtered_df.to_csv(output_file, sep="\t", index=False, header=False)

def process_file(input_file, output_file, f_T):
    if not os.path.exists(input_file):
        print(f"Error: The file {input_file} does not exist.")
        return

    if os.path.getsize(input_file) == 0:
        print(f"Error: The file {input_file} is empty.")
        return
    # Read the input file
    df = pd.read_csv(input_file, sep="\t", header=None)

    # Prepare lists to store results
    tumor_copies = []
    event_types = []
    change_values = []

    # Iterate through each row in the DataFrame
    for index, row in df.iterrows():
        log2_ratio = row[4]  # Get the log2 ratio from the fifth column

        # Classify the event
        C_T_rounded, event_type, change_value = classify_copy_number_event(log2_ratio, f_T)

        # Append results to lists
        tumor_copies.append(C_T_rounded)
        event_types.append(event_type)
        change_values.append(change_value)

    # Add new columns to the DataFrame
    df['Tumor Copy'] = tumor_copies
    df['Event.Type'] = event_types
    df['Change.Value'] = change_values

    # Create a new DataFrame with the desired columns and headers
    output_df = pd.DataFrame({
        'Chrom': df[0],
        'Start': df[1],
        'End': df[2],
        'Log ratio': df[4],
        'Tumor Copy': tumor_copies,
        'Event.Type': event_types,
        'Change.Value': change_values
    })

    # Replace negative values in 'Tumor Copy' with "Indeterminated" and update corresponding columns
   # output_df.loc[output_df['Tumor Copy'] == "Indeterminated", ['Event Type', 'Change Value']] = "Indeterminated"

    # Write the updated DataFrame to the output file with headers
    output_df.to_csv(output_file, sep="\t", index=False, header=True)  # Save as tab-separated values with headers

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Process copy number events from a file.')
    parser.add_argument('file1', type=str, help='Path to the seg bed file')
    parser.add_argument('file2', type=str, help='Path to the occ fusion file')
    #parser.add_argument('matched_file', type=str, help='Path to the intermediate matched file')
    parser.add_argument('output_file', type=str, help='Path to the tumor copy output file')
    parser.add_argument('filtered_output_file', type=str, help='Path to the bins seg filtered output file')
    parser.add_argument('f_T', type=float, help='Tumor fraction (e.g., 0.8 for 80% tumor)')

    args = parser.parse_args()

    # Match rows and create a new file
    #match_and_create_file(args.file1, args.file2, args.matched_file)

    # Process the matched file
   # process_file(args.matched_file, args.output_file, args.f_T)
    filter_by_interval(args.file1, args.file2, args.filtered_output_file)
    process_file(args.filtered_output_file, args.output_file, args.f_T)
    ##filter_by_interval(args.file1, args.file2, args.filtered_output_file)

