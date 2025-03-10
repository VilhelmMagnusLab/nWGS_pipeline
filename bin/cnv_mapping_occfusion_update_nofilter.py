#!/usr/bin/env python3
import math
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

def process_file(input_file, output_file, f_T):
    # Read the input file
   # df = pd.read_csv(input_file, sep="\t", header=None)  # Assuming tab-separated 
    df = pd.read_csv(input_file, sep="\t", header=None, skiprows=1) 

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
    df['Event Type'] = event_types
    df['Change Value'] = change_values

    # Create a new DataFrame with the desired columns and headers
    output_df = pd.DataFrame({
        'Chrom': df[0],
        'Start': df[1],
        'End': df[2],
        'Log ratio': df[4],
        'Tumor Copy': tumor_copies,
        'Event Type': event_types,
        'Change Value': change_values
    })

    # Replace negative values in 'Tumor Copy' with "Indeterminated" and update corresponding columns
    output_df.loc[output_df['Tumor Copy'] == "Indeterminated", ['Event Type', 'Change Value']] = "Indeterminated"

    # Write the updated DataFrame to the output file with headers
    output_df.to_csv(output_file, sep="\t", index=False, header=True)  # Save as tab-separated values with headers

if __name__ == "__main__":
    # Set up argument parsing
    parser = argparse.ArgumentParser(description='Process copy number events from a file.')
    parser.add_argument('input_file', type=str, help='Path to the bed input file')
    parser.add_argument('output_file', type=str, help='Path to the output file')
    parser.add_argument('f_T', type=float, help='Tumor fraction (e.g., 0.8 for 80% tumor)')

    args = parser.parse_args()

    # Process the files
    process_file(args.input_file, args.output_file, args.f_T)
