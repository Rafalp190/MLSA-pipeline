from datetime import datetime
import os
import pandas as pd

# Ensure the log directory exists
log_dir = "/mnt/output/log/"
if not os.path.exists(log_dir):
    os.makedirs(log_dir)

# Calculate the current date
current_date = datetime.now().strftime('%Y-%m-%d-%H-%M')

def get_loci():
    # Load the Excel file to extract loci (columns)
    df = pd.read_excel("/mnt/data/input.xlsx")
    # Return sanitized names for Snakemake but maintain the original names for mapping
    return [locus.replace(" ", "_") for locus in df.columns[1:]]

configfile: "/mnt/config.yaml"

wildcard_constraints:
    locus="[a-zA-Z0-9_]+"

# Top-level rule that defines the final target output (without wildcards)
rule all:
    threads: int(config["cores"])
    input:
        "/mnt/output/output_inference.nwk",  # Final pipeline output
        "/mnt/output/dashboard_ready.log"   # Final dashboard step

# Rule to run the data loader
rule data_loader:
    input:
        excel="/mnt/data/input.xlsx",
        seq_dir="/mnt/data/seq/"
    output:
        log=f"/mnt/output/log/data_load.log"
    params:
        date=current_date
    script:
        "/app/src/data_loader.py"

rule alignment:
    input:
        excel="/mnt/data/input.xlsx",
        data_load_log="/mnt/output/log/data_load.log",
        seq_dir="/mnt/data/seq/"
    output:
        aligned="/mnt/output/aligned/{locus}.fasta"
    params:
        aligner=config["aligner"],
        muscle_params=config.get("muscle_params", {}),  # Default to an empty dict if not provided
        clustalw_params=config.get("clustalw_params", {})  # Default to an empty dict if not provided
    log:
        "/mnt/output/log/alignment_{locus}.log"
    threads:
        config.get("cores", 1)
    script:
        "/app/src/alignment.py"

rule concatenation:
    input:
         expand("/mnt/output/aligned/{locus}.fasta", locus=get_loci())
    output:
        concatenated="/mnt/output/concatenated/concatenated_sequences.fasta",
        align_image="/mnt/output/alignment_plot.png"
    params:
        concatenation_tool=config["concatenator"],
        amas_params=config.get("amas_params", {}),  # Default to an empty dict if not provided
        seqkit_params=config.get("seqkit_params", {})  # Default to an empty dict if not provided
    log:
        "/mnt/output/log/concatenation.log"
    threads:
        config.get("cores", 1)
    script:
        "/app/src/concatenation.py"

rule inference:
    input:
        concatenated="/mnt/output/concatenated/concatenated_sequences.fasta"
    output:
        tree="/mnt/output/output_inference.nwk"
    params:
        inference_tool=config["inferer"],
        mr_bayes_params=config.get("mr_bayes_params", {}),
        dendropy_params=config.get("dendropy_params", {}),
        biophylo_params=config.get("biophylo_params", {})
    log:
        "/mnt/output/log/inference.log"
    threads:
        config.get("cores", 1)
    script:
        "/app/src/inference.py"

# Rule to run the Shiny dashboard
rule dashboard:
    input:
        inference_log="/mnt/output/output_inference.nwk"  # Ensure inference step is complete
    output:
        touch("/mnt/output/dashboard_ready.log")
    shell:
        """
        shiny run --reload /app/src/app.py --host 0.0.0.0 --port 8000
        touch {output}
        """