#!/usr/bin/env python3
import os
import shutil
import sys

def copy_files_from_root(root_folder, destination_folder):
    """
    Copy all files from all folders within the root folder to a destination folder.
    
    Args:
        root_folder (str): Path to the root folder containing subfolders.
        destination_folder (str): Path to the destination folder where files will be copied.
    """
    # Ensure destination folder exists
    os.makedirs(destination_folder, exist_ok=True)

    # Get all subfolders in the root folder
    subfolders = [os.path.join(root_folder, folder) for folder in os.listdir(root_folder) if os.path.isdir(os.path.join(root_folder, folder))]

    for folder in subfolders:
        for root, _, files in os.walk(folder):
            for file in files:
                source_file = os.path.join(root, file)
                dest_file = os.path.join(destination_folder, file)

                try:
                    # Copy file to destination folder
                    shutil.copy2(source_file, dest_file)
                    print(f"Copied: {source_file} -> {dest_file}")
                except Exception as e:
                    print(f"Error copying {source_file} to {dest_file}: {e}")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python script.py <root_folder> <destination_folder>")
        sys.exit(1)

    root_folder = sys.argv[1]
    destination_folder = sys.argv[2]

    copy_files_from_root(root_folder, destination_folder)

