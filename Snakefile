from datetime import datetime
import os

# Ensure the log directory exists
log_dir = "/mnt/data/log/snakemake/"
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

# Calculate the current date
current_date = datetime.now().strftime('%Y-%m-%d-%H-%M')

# Top-level rule that defines the final target output (without wildcards)
rule all:
    input:
        f"/mnt/data/log/data_load_{current_date}.txt"

# Rule to run the data loader
rule data_loader:
    input:
        excel="/mnt/data/input.xlsx",
        seq_dir="/mnt/data/seq/"
    output:
        log=f"/mnt/data/log/data_load_{current_date}.txt"
    log:
        f"/mnt/data/log/snakemake/snakemake_{current_date}.log"  # Place log before script
    params:
        date=current_date
    script:
        "/app/src/data_loader.py"
