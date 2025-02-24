from Bio import SeqIO
import plotly.graph_objects as go
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from Bio import AlignIO
from Bio.Align import substitution_matrices

def parse_fasta(fasta_files):
    """Parses FASTA files and returns a list of sequence dictionaries."""
    sequences = []
    for fasta_file in fasta_files:
        try:  # Handle potential file errors
            for record in SeqIO.parse(fasta_file, "fasta"):
                sequences.append({
                    'id': record.id,
                    'sequence': str(record.seq),
                    'length': len(record.seq)
                })
        except Exception as e:
            print(f"Error parsing {fasta_file}: {e}")
    return sequences

def create_separate_ideograms(sequences, output_dir="ideograms"):
    """Creates and saves separate ideogram plots for each sequence."""

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)  # Create the directory if it doesn't exist

    for seq in sequences:
        fig = go.Figure()

        # Draw ideogram as a horizontal line
        fig.add_trace(go.Scatter(
            x=[0, seq['length']],
            y=[0, 0],
            mode='lines',
            line=dict(color='black', width=4),
            name=seq['id'],
            hoverinfo='text',
            hovertext=f"ID: {seq['id']}<br>Length: {seq['length']} bp"
        ))

        # Example marker (midpoint) - You can customize this
        fig.add_trace(go.Scatter(
            x=[seq['length'] // 2],
            y=[0],
            mode='markers',
            marker=dict(color='red', size=10),
            hoverinfo='text',
            hovertext=f"Feature at {seq['length'] // 2} bp"
        ))

        fig.update_layout(
            title=f"Ideogram for {seq['id']}",
            xaxis_title="Base Pair Position",
            yaxis=dict(showticklabels=False),  # Hide y-axis labels
            plot_bgcolor='white',
            height=150,
            showlegend=False
        )

        # Save the figure to a file
        filename = os.path.join(output_dir, f"{seq['id']}.html")  # Use .html for interactive plots
        fig.write_html(filename)  # Or fig.write_image() for static images
        print(f"Ideogram saved to: {filename}")

def create_sequence_alignment_graph(alignments, output_file="alignment_graph.png"):
    """Creates a sequence alignment graph using matplotlib."""
    fig, ax = plt.subplots(figsize=(10, len(alignments) * 0.5))  # Adjust figure size

    for i, alignment in enumerate(alignments):
        seq1, seq2 = alignment.split('\n')  # Assuming two sequences per alignment
        seq1 = seq1.replace("-","")
        seq2 = seq2.replace("-","")
        for j, (base1, base2) in enumerate(zip(seq1, seq2)):
            color = 'green' if base1 == base2 else 'red'  # Match/mismatch coloring
            ax.plot([j, j], [i, i + 0.5], color=color, linewidth=2)  # Vertical lines

    ax.set_yticks(np.arange(0.25, len(alignments) * 0.5 + 0.25, 0.5))
    ax.set_yticklabels([f"Seq {i+1}" for i in range(len(alignments)//2)])
    ax.set_xlabel("Base Position")
    ax.set_title("Sequence Alignment")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

def create_probability_distribution_graph(data, output_file="probability_distribution.png"):
    """Creates a probability distribution graph (example using histogram)."""
    plt.figure(figsize=(8, 6))
    plt.hist(data, bins=20, density=True, alpha=0.7, color='skyblue')  # Example: histogram
    plt.title("Probability Distribution (Example)")
    plt.xlabel("Data Value")
    plt.ylabel("Probability Density")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

def create_mixture_deconvolution_plot(data, output_file="mixture_deconvolution.png"):
    """Creates a mixture deconvolution plot (example using stacked bar chart)."""
    # Example data (replace with your actual data)
    categories = ['Component A', 'Component B', 'Component C']
    values1 = [30, 40, 30]
    values2 = [10, 50, 40]
    values3 = [20, 25, 55]

    plt.figure(figsize=(8, 6))
    plt.bar(categories, values1, color='skyblue', label='Sample 1')
    plt.bar(categories, values2, bottom=values1, color='lightcoral', label='Sample 2')
    plt.bar(categories, values3, bottom=np.array(values1) + np.array(values2), color='lightgreen', label='Sample 3')

    plt.title("Mixture Deconvolution (Example)")
    plt.xlabel("Components")
    plt.ylabel("Proportion")
    plt.legend()
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

def create_phylogenetic_tree(tree_data, output_file="phylogenetic_tree.png"): # Placeholder
    """Creates a phylogenetic tree (requires a tree library like DendroPy or ete3)."""
    print("Phylogenetic tree generation is complex and requires a dedicated library.  This is a placeholder.")
    # You would use a library like DendroPy or ete3 here to create and save the tree.
    # Example (conceptual - replace with actual library usage):
    # from dendropy import Tree
    # tree = Tree.get(data=tree_data, schema="newick")  # Assuming Newick format
    # tree.plot(output_file)


def create_heatmap(data, output_file="heatmap.png"):
    """Creates a heatmap using matplotlib."""
    plt.figure(figsize=(8, 6))
    plt.imshow(data, cmap='viridis', aspect='auto')  # Example heatmap
    plt.colorbar(label="Value")
    plt.title("Heatmap")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

def create_allele_frequency_histogram(allele_frequencies, output_file="allele_frequency_histogram.png"):
    """Creates a histogram of allele frequencies."""
    plt.figure(figsize=(8, 6))
    plt.hist(allele_frequencies, bins=20, color='skyblue')
    plt.title("Allele Frequency Histogram")
    plt.xlabel("Allele Frequency")
    plt.ylabel("Count")
    plt.tight_layout()
    plt.savefig(output_file)
    plt.show()

def process_dna_visualization(fasta_dir="dataset"):
    """Processes FASTA files and creates various visualizations."""
    try:
        fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith(".fasta")]
        sequences = parse_fasta(fasta_files)
        create_separate_ideograms(sequences)

        # Example usage of other graph functions (replace with your data and logic):
        # Example: Alignment graph (requires alignment data)
        # alignments = [ # Placeholder. You need to generate alignment objects.
        #     "ATGC-GTA\nAT-C-GTA",
        #     "G-T-C\nGA-TC"
        # ]

        # create_sequence_alignment_graph(alignments)

        # Example: Probability distribution
        # data = np.random.normal(0, 1, 100)  # Example data
        # create_probability_distribution_graph(data)

        # ... (similarly for other graph types)

    except FileNotFoundError:
        print(f"Error: Directory '{fasta_dir}' not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")



# Example usage:
# process_dna_visualization()
