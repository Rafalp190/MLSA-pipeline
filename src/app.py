from shiny import App, render, ui
from Bio import Phylo, AlignIO
from pymsaviz import MsaViz
from io import StringIO
import plotly.figure_factory as ff
import plotly.graph_objects as go
import numpy as np

# Define UI
app_ui = ui.page_fluid(
    ui.h1("MLSA Pipeline Visualizer"),
    ui.p("Check out the results of the multi locus sequence analysis you just performed"),
    ui.navset_card_pill(
        ui.nav_panel("Dendrogram",
                    ui.card(
                        ui.card_header("Dendrogram Plot"),
                        ui.layout_sidebar(
                            ui.sidebar("", 
                                    ui.download_button("download_nwk", "Download Newick File") ,position="right", background="#f8f8f8"),
                                    ui.output_ui("dendrogram_plot")  # Use output_ui for HTML-based Plotly rendering
                     ))
   
        ),
        ui.nav_panel("Alignment plot",
                    ui.card(
                        ui.layout_columns(
                            ui.download_button("download_align", "Download Alignment image"),
                            ui.download_button("download_align_FA", "Download Alignment FASTA")
                        ),
                        ui.div(
                            ui.output_image("alignment_plot"),
                            style="max-height: 500px; overflow: scroll;"  # CSS to limit image height
                        )
                    ))))


# Server logic
def server(input, output, session):
    
    def load_tree():
        try:
            # Read the Newick file into a Phylo tree
            with open("/mnt/output/output_inference.nwk", "r") as f:
                tree = Phylo.read(f, "newick")
            return tree
        except FileNotFoundError:
            return None

    def newick_to_distance_matrix(tree):
        """
        Convert Newick tree to a distance matrix for use with scipy's linkage.
        This function calculates pairwise distances between all terminal nodes.
        """
        # Get terminal clades (leaf nodes)
        terminals = tree.get_terminals()
        # Replace underscores with spaces in the labels for display
        labels = [f"<i>{term.name.replace('_', ' ')}</i>" for term in terminals]  
        n = len(terminals)

        # Create an empty distance matrix
        distance_matrix = np.zeros((n, n))

        # Fill the distance matrix with pairwise distances
        for i in range(n):
            for j in range(i + 1, n):
                distance = tree.distance(terminals[i], terminals[j])
                distance_matrix[i, j] = distance
                distance_matrix[j, i] = distance

        return distance_matrix, labels

    def get_plotly_dendrogram(tree):
        # Convert the Newick tree to a distance matrix
        distance_matrix, labels = newick_to_distance_matrix(tree)

        # Generate a dendrogram with Plotly's figure factory directly from the distance matrix
        fig = ff.create_dendrogram(distance_matrix, labels=labels, color_threshold=1.5)
        
        # Update layout for italicized axis titles and labels
        fig.update_layout(
            title="Dendrogram from Newick File",
            xaxis=dict(title="<i>Clade</i>", showticklabels=True),
            yaxis=dict(title="<i>Distance</i>"),
            height=600,
            width=1100
        )
        # Add branch lengths as text annotations
        dendrogram_data = fig['data']
        for trace in dendrogram_data:
            if trace['mode'] == 'lines':
                x = trace['x']
                y = trace['y']
                for i in range(0, len(x) - 1, 2):
                    midpoint_x = ((x[i] + x[i + 1]) / 2)+1
                    midpoint_y = (y[i] + y[i + 1]) / 2
                    length = abs(y[i] - y[i + 1])
                    fig.add_trace(go.Scatter(
                        x=[midpoint_x],
                        y=[midpoint_y],
                        mode='text',
                        text=[f"{length:.2f}"],
                        textposition="middle right",
                        showlegend=False
                    ))

        # Update font style for labels (italic)
        fig.update_traces(textfont=dict(size=10, color="black", family="Times New Roman"))
        
        return fig

    # Render the dendrogram plot as HTML
    @output
    @render.ui
    def dendrogram_plot():
        tree_data = load_tree()
        if tree_data is None:
            return ui.div("File not found: /mnt/output/output_inference.nwk", style="color:red; font-size:20px;")
        
        try:
            fig = get_plotly_dendrogram(tree_data)
        except ValueError as e:
            return ui.div(f"Error: {e}", style="color:red; font-size:20px;")

        return ui.HTML(fig.to_html(full_html=False))  # Render Plotly figure as HTML
    
    # Render the alignment visualization as an image
    @output
    @render.image
    def alignment_plot():
        # Generate the alignment image file
        return {"src": "/mnt/output/alignment_plot.png", "alt": "Sequence Alignment"}

    # Create a downloadable file for the Newick data
    @output
    @render.download
    def download_nwk():
        # Directly specify the path of the Newick file for download
        return "/mnt/output/output_inference.nwk"
    
    @output
    @render.download
    def download_align():
        # Directly specify the path of the Newick file for download
        return "/mnt/output/alignment_plot.png"
    
        # Create a downloadable file for the Newick data
    @output
    @render.download
    def download_align_FA():
        # Directly specify the path of the Newick file for download
        return "/mnt/output/concatenated/concatenated_sequences.fasta"
# Run the app
app = App(app_ui, server)
