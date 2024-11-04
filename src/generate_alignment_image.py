from Bio import Phylo, AlignIO
from pymsaviz import MsaViz as msaviz
from io import StringIO
import plotly.figure_factory as ff
import plotly.graph_objects as go
import numpy as np

def generate_alignment_image():
        # Load the sequence alignment
        alignment = AlignIO.read("./mnt/output/alignment.fasta", "fasta")

        # Create and save alignment visualization using pyMSAviz
        msa_viz = msaviz(alignment, wrap_length=200, show_grid=True, show_consensus=True)
        msa_viz.savefig("./mnt/output/alignment_plot.png")