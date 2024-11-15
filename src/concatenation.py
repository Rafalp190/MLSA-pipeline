import os
import subprocess
import sys
from Bio import AlignIO
from pymsaviz import MsaViz as msaviz

# Function to generate an alignment image to use as output later
def generate_alignment_image():
    # Load the sequence alignment
    alignment = AlignIO.read("/mnt/output/concatenated/concatenated_sequences.fasta", "fasta")

    # Create and save alignment visualization using pyMSAviz
    msa_viz = msaviz(alignment, wrap_length=200, show_grid=True, show_consensus=True)
    msa_viz.savefig("/mnt/output/alignment_plot.png")

# Load configuration from Snakemake
alignment_files = snakemake.input
output_file = snakemake.output.concatenated
concatenation_tool = snakemake.params.concatenation_tool
amas_params = snakemake.params.amas_params  # Parameters for AMAS
seqkit_params = snakemake.params.seqkit_params  # Parameters for SeqKit
threads = snakemake.threads
log_file = snakemake.log[0]

try:
    if concatenation_tool == "amas":
        # Use AMAS via subprocess with the correct CLI format
        partition_file = output_file.replace(".fasta", "_partition.txt")
        cmd = [
            "python3", "/app/src/AMAS.py", "concat", "-f", "fasta", "-d", "dna", "-i"
        ] + list(alignment_files) + ["-t", output_file, "-p", partition_file, "-c", str(threads)] 

        # Add optional parameters from amas_params
        if amas_params.get("check_align", False):
            cmd.append("-e")

        print(f"Running command: {' '.join(cmd)}")
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        with open(log_file, "w") as log:
            log.write(f"Command: {' '.join(cmd)}\n")
            log.write("Standard Output:\n" + result.stdout + "\n")
            log.write("Standard Error:\n" + result.stderr + "\n")

        if result.returncode != 0:
            raise RuntimeError(f"AMAS concatenation failed with exit code {result.returncode}. See log for details.")

    elif concatenation_tool == "seqkit":
        # Use SeqKit for concatenation
        cmd = ["seqkit", "concat"] + list(alignment_files) + ["-o", output_file]

        if seqkit_params.get("full", False):
            cmd.append("--full")
        
        print(f"Running command: {' '.join(cmd)}")
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        with open(log_file, "w") as log:
            log.write(f"Command: {' '.join(cmd)}\n")
            log.write("Standard Output:\n" + result.stdout + "\n")
            log.write("Standard Error:\n" + result.stderr + "\n")

        if result.returncode != 0:
            raise RuntimeError(f"SeqKit concatenation failed with exit code {result.returncode}. See log for details.")

    else:
        raise ValueError(f"Unknown concatenation tool: {concatenation_tool}")
    
    # Generate alignment visualization after successful concatenation
    generate_alignment_image()

except Exception as e:
    with open(log_file, "a") as log:
        log.write(f"Error: {str(e)}\n")
    sys.exit(1)
