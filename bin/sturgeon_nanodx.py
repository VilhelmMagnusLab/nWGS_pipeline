import os
import pandas as pd
import sys

def process_files(input_folder, output_folder):
    """
    Process all *.EpicSelect_sturgeon_EpicSelect.bed files in the input folder.
    Replace methylation_call column values:
        - 1 -> 1.0
        - 0 -> -1.0
    Write the modified data to new files in the output folder.

    Args:
        input_folder (str): Path to the folder containing the .EpicSelect_sturgeon_EpicSelect.bed files.
        output_folder (str): Path to the folder where modified files will be saved.
    """

    # Ensure the input folder exists
    if not os.path.isdir(input_folder):
        print(f"Error: Input folder '{input_folder}' does not exist.")
        sys.exit(1)

    # Ensure the output folder exists or create it
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
        print(f"Created output folder: {output_folder}")

    # Process all files ending with '.EpicSelect_sturgeon_EpicSelect.bed'
    for filename in os.listdir(input_folder):
        if filename.endswith(".EpicSelect_sturgeon_EpicSelect.bed"):
            input_file_path = os.path.join(input_folder, filename)
            output_file_path = os.path.join(output_folder, filename.replace(".EpicSelect_sturgeon_EpicSelect.bed", "_updated.bed"))

            try:
                # Load file into a DataFrame
                df = pd.read_csv(input_file_path, sep="\t", header=0)

                # Check if 'methylation_call' exists
                if 'methylation_call' not in df.columns:
                    print(f"Warning: 'methylation_call' column not found in {filename}. Skipping...")
                    continue

                # Replace 1 -> 1.0 and 0 -> -1.0
                df['methylation_call'] = df['methylation_call'].replace({1: 1.0, 0: -1.0})

                # Write the modified DataFrame to the output file
                df.to_csv(output_file_path, sep="\t", index=False)
                print(f"Processed and saved: {output_file_path}")

            except Exception as e:
                print(f"Error processing file {filename}: {e}")

if __name__ == "__main__":
    # Check for correct number of arguments
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_folder> <output_folder>")
        sys.exit(1)

    # Get input and output folder paths from command-line arguments
    input_folder = sys.argv[1]
    output_folder = sys.argv[2]

    # Call the function
    process_files(input_folder, output_folder)

