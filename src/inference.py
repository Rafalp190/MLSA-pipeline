import os
import sys
from Bio import Phylo, AlignIO, SeqIO
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceCalculator
from dendropy import DnaCharacterMatrix
import subprocess
import glob
import re

concatenated_file = snakemake.input.concatenated
output_tree = snakemake.output.tree
inference_tool = snakemake.params.inference_tool
log_file = snakemake.log[0]
threads = snakemake.threads

biophylo_params = snakemake.params.get("biophylo_params", {})
dendropy_params = snakemake.params.get("dendropy_params", {})
mr_bayes_params = snakemake.params.get("mr_bayes_params", {})

try:
    if inference_tool == "biophylo":
        # Use BioPhylo for Neighbor Joining Tree
        distance_metric = biophylo_params.get("distance_metric", "identity")
        outgroup = biophylo_params.get("outgroup")

        alignment = AlignIO.read(concatenated_file, "fasta")
        calculator = DistanceCalculator(distance_metric)
        dm = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor(calculator)
        tree = constructor.nj(dm)

        if outgroup:
            tree.root_with_outgroup(outgroup)

        Phylo.write(tree, output_tree, "newick")

    elif inference_tool == "dendropy":
        # Use BioPhylo for Neighbor Joining Tree
        distance_metric = biophylo_params.get("distance_metric", "identity")
        outgroup = biophylo_params.get("outgroup")

        alignment = AlignIO.read(concatenated_file, "fasta")
        calculator = DistanceCalculator(distance_metric)
        dm = calculator.get_distance(alignment)
        constructor = DistanceTreeConstructor(calculator)
        tree = constructor.nj(dm)

        if outgroup:
            tree.root_with_outgroup(outgroup)

        Phylo.write(tree, output_tree, "newick")

    elif inference_tool == "mrbayes":
        # Check the number of taxa
        alignment = AlignIO.read(concatenated_file, "fasta")
        num_taxa = len(alignment)
        if num_taxa < 4:
            raise RuntimeError(f"MrBayes requires at least 4 taxa, but only {num_taxa} were provided.")

        # Convert FASTA to Nexus format
        nexus_file = concatenated_file.replace(".fasta", ".nex")
        SeqIO.convert(concatenated_file, "fasta", nexus_file, "nexus", "DNA")

        # Create the MrBayes command file
        mrbayes_input = f"{os.path.splitext(output_tree)[0]}.mb"
        with open(mrbayes_input, "w") as mb_input:
            mb_input.write(f"begin mrbayes;\n")
            mb_input.write(f"   set autoclose=yes;\n")
            mb_input.write(f"   execute {nexus_file};\n")
            mb_input.write(f"   mcmc;\n")
            mb_input.write(f"   sump;\n")
            mb_input.write(f"   sumt;\n")
            mb_input.write(f"end;\n")

        # Execute MrBayes using the command file
        cmd = ["mb", mrbayes_input]
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

        # Check for errors in MrBayes execution
        if result.returncode != 0:
            raise RuntimeError(f"MrBayes failed with error: {result.stderr}")
        # Convert Nexus tree output to Newick format
        latest_tree_file = "/mnt/output/concatenated/concatenated_sequences.nex.con.tre"
        if not os.path.exists(latest_tree_file):
            raise RuntimeError("Tree file (.con.tre) not generated by MrBayes.")
        
        cleaned_tree_file = latest_tree_file.replace(".con.tre", ".cleaned.nex")
        with open(latest_tree_file, "r") as infile, open(cleaned_tree_file, "w") as outfile:
            for line in infile:
                # Remove metadata annotations within brackets
                clean_line = re.sub(r"\[.*?\]", "", line)
                outfile.write(clean_line)

        # Read the cleaned Nexus tree and write it as Newick format
        tree = Phylo.read(cleaned_tree_file, "nexus")
        Phylo.write(tree, output_tree, "newick")

    else:
        raise ValueError(f"Unknown inference tool: {inference_tool}")
    with open(log_file, "w") as log:
        log.write("Inference completed successfully.\n")

except Exception as e:
    with open(log_file, "a") as log:
        log.write(f"Error: {str(e)}\n")
    sys.exit(1)
