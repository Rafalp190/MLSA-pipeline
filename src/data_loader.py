import pandas as pd
import os
from datetime import datetime
from Bio import SeqIO

# Define paths (assuming you're running from the root of the project)
input_excel = '/mnt/data/input.xlsx'  # Path to the Excel file (mounted device)
seq_dir = '/mnt/data/seq/'            # Path to the directory containing the sequences
log_dir = '/mnt/output/log/'            # Log directory

# Create the log directory if it doesn't exist
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

# Load the Excel file
try:
    df = pd.read_excel(input_excel)
    print(f"Successfully loaded the Excel file: {input_excel}")
except Exception as e:
    print(f"Error loading Excel file: {e}")
    exit(1)

# Skip the first column (Organism names) and process only the loci columns
loci_columns = df.columns[1:]  # Skip the first column, which is organism names

# Calculate number of organisms (rows) and loci (columns)
num_organisms = df.shape[0]
num_loci = len(loci_columns)

# Create a terminal report
print(f"Number of organisms loaded: {num_organisms}")
print(f"Number of loci per organism: {num_loci}")

# Function to validate a FASTA or FASTQ file
def validate_sequence_file(file_path):
    if file_path.endswith(".fasta") or file_path.endswith(".fa"):
        file_format = "fasta"
    elif file_path.endswith(".fastq"):
        file_format = "fastq"
    else:
        return False, "Invalid file extension"
    
    try:
        # Test reading the file using Biopython
        with open(file_path, "r") as handle:
            records = list(SeqIO.parse(handle, file_format))
            if not records:
                return False, "No records found"
        return True, "Valid"
    except Exception as e:
        return False, f"Error: {str(e)}"

# Check if each file in the matrix exists and is valid
missing_files = []
invalid_files = []

for organism in df.index:
    for locus in loci_columns:
        file_path = os.path.join(seq_dir, df.loc[organism, locus])
        if not os.path.isfile(file_path):
            missing_files.append(file_path)
        else:
            # Validate the sequence file
            is_valid, message = validate_sequence_file(file_path)
            if not is_valid:
                invalid_files.append((file_path, message))

# Print any missing or invalid files
if missing_files:
    print("Missing files:")
    for missing_file in missing_files:
        print(f" - {missing_file}")
else:
    print("All files are present.")

if invalid_files:
    print("Invalid files:")
    for invalid_file, message in invalid_files:
        print(f" - {invalid_file}: {message}")
else:
    print("All sequence files are valid.")

# Generate the log file
current_date = datetime.now().strftime('%Y-%m-%d-%H-%M')
log_file = os.path.join(log_dir, f"data_load.log")

# Write to log file
with open(log_file, 'w') as log:
    log.write(f"Log generated on {current_date}\n")
    log.write(f"Number of organisms: {num_organisms}\n")
    log.write(f"Number of loci per organism: {num_loci}\n")
    
    if missing_files:
        log.write("Missing files:\n")
        for missing_file in missing_files:
            log.write(f" - {missing_file}\n")
    else:
        log.write("All sequence files are present.\n")
    
    if invalid_files:
        log.write("Invalid files:\n")
        for invalid_file, message in invalid_files:
            log.write(f" - {invalid_file}: {message}\n")
    else:
        log.write("All sequence files are valid.\n")

print(f"Log file created: {log_file}")
