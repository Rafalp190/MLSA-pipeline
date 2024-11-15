import os
import subprocess
import sys
import pandas as pd

# Load configuration from Snakemake
excel = snakemake.input.excel
seq_dir = snakemake.input.seq_dir
output_file = snakemake.output.aligned
aligner = snakemake.params.aligner
log_file = snakemake.log[0]
threads = snakemake.threads
muscle_params = snakemake.params.muscle_params  # Dynamically pass muscle parameters
clustalw_params = snakemake.params.clustalw_params  # Dynamically pass ClustalW parameters

# Reverse sanitization for dataframe access
safe_locus = snakemake.wildcards.locus
original_locus = safe_locus.replace("_", " ")

# Load the input matrix
df = pd.read_excel(excel)

# Create a temporary FASTA file for the sequences
temp_fasta = f"/mnt/output/temp_{safe_locus}.fasta"
with open(temp_fasta, "w") as fasta:
    for organism, seq_file in zip(df.iloc[:, 0], df[original_locus]):
        seq_path = os.path.join(seq_dir, seq_file)
        if not os.path.isfile(seq_path):
            print(f"Sequence file missing: {seq_path}", file=sys.stderr)
            sys.exit(1)

        # Ensure each sequence has a proper header and no extra line breaks
        with open(seq_path, "r") as seq:
            lines = seq.readlines()
            fasta.write(f">{organism}\n")
            fasta.write("".join(line.strip() for line in lines if not line.startswith(">")) + "\n")

# Construct MUSCLE or ClustalW command
if aligner == "muscle":
    # Start with the basic MUSCLE command
    cmd = [
        "muscle", "-align", temp_fasta, "-output", output_file, "--threads", str(threads)
    ]

    # Add additional MUSCLE parameters if specified
    if "perturb" in muscle_params:
        cmd += ["-perturb", str(muscle_params["perturb"])]
    if "perm" in muscle_params:
        cmd += ["-perm", muscle_params["perm"]]
    if "stratified" in muscle_params and muscle_params["stratified"]:
        cmd += ["-stratified"]
    if "diversified" in muscle_params and muscle_params["diversified"]:
        cmd += ["-diversified"]
    if "replicates" in muscle_params:
        cmd += ["-replicates", str(muscle_params["replicates"])]
    if "consiters" in muscle_params:
        cmd += ["-consiters", str(muscle_params["consiters"])]
    if "refineiters" in muscle_params:
        cmd += ["-refineiters", str(muscle_params["refineiters"])]

elif aligner == "clustalw":
    # Start with the basic ClustalW command
    cmd = [
        "clustalw", f"-INFILE={temp_fasta}", f"-OUTFILE={output_file}", "-OUTPUT=FASTA"
    ]

    # Add additional ClustalW parameters if specified
    if "gapopen" in clustalw_params:
        cmd += ["-GAPOPEN={}".format(clustalw_params["gapopen"])]
    if "gapextend" in clustalw_params:
        cmd += ["-GAPEXT={}".format(clustalw_params["gapextend"])]
    if "numiter" in clustalw_params:
        cmd += ["-NUMITER={}".format(clustalw_params["numiter"])]
    if "num_threads" in clustalw_params:
        cmd += ["-NUMTHREADS={}".format(str(threads))]
    if "matrix" in clustalw_params:
        cmd += ["-MATRIX={}".format(clustalw_params["matrix"])]
else:
    raise ValueError(f"Unknown aligner: {aligner}")

# Run the alignment
try:
    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    with open(log_file, "w") as log:
        log.write(f"Command: {' '.join(cmd)}\n")
        log.write("Standard Output:\n" + result.stdout + "\n")
        log.write("Standard Error:\n" + result.stderr + "\n")

    if result.returncode != 0:
        raise RuntimeError(f"Alignment failed with exit code {result.returncode}. See log for details.")

    # Handle ClustalW .dnd file after command execution
    if aligner == "clustalw":
        temp_dnd = f"/mnt/output/temp_{safe_locus}.dnd"
        if os.path.exists(temp_dnd):
            with open(temp_dnd, "r") as dnd_file:
                with open(log_file, "a") as log:
                    log.write("\nGuide Tree (.dnd file):\n")
                    log.write(dnd_file.read())
            os.remove(temp_dnd)

except Exception as e:
    with open(log_file, "a") as log:
        log.write(f"Error: {str(e)}\n")
    sys.exit(1)

finally:
    if os.path.exists(temp_fasta):
        os.remove(temp_fasta)
