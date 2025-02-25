import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
import seaborn as sns
from Bio import SeqIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceTreeConstructor, DistanceMatrix
from Bio.Align import MultipleSeqAlignment, PairwiseAligner
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler

# Ensure SVG output
plt.switch_backend('Agg')

# Create the "graphs" folder if it doesn't exist
if not os.path.exists("graphs"):
    os.makedirs("graphs")

def load_sequences(folder_path):
    """Load DNA sequences from all FASTA files in a directory."""
    sequences = {}
    for filename in os.listdir(folder_path):
        if filename.endswith(('.fasta', '.fa', '.fna')):
            filepath = os.path.join(folder_path, filename)
            try:
                records = list(SeqIO.parse(filepath, "fasta"))
                if records:
                    sequences[filename] = records
            except Exception as e:
                print(f"Error parsing {filename}: {e}")
    return sequences


def save_svg(fig, filename):
    """Save a matplotlib figure as SVG into the 'graphs' folder."""
    filepath = os.path.join("graphs", filename)  # Add "graphs" folder to the path
    fig.savefig(filepath, format='svg', bbox_inches='tight')

def plot_mixture_deconvolution():
    """Simulates a mixture deconvolution plot."""
    fig, ax = plt.subplots(figsize=(10, 4))
    x = np.linspace(100, 1000, 100)
    colors = ['blue', 'green', 'red']
    contributors = ['Contributor 1', 'Contributor 2', 'Contributor 3']

    for i in range(3):
        # Simulate peaks for each contributor
        peak1 = np.exp(-0.005 * (x - np.random.randint(200, 400)) ** 2) * np.random.uniform(5, 10)
        peak2 = np.exp(-0.008 * (x - np.random.randint(500, 700)) ** 2) * np.random.uniform(3, 8)
        peak3 = np.exp(-0.01 * (x - np.random.randint(700, 900)) ** 2) * np.random.uniform(6, 12)
        total_peak = peak1 + peak2 + peak3  # Combine peaks

        ax.plot(x, total_peak, color='black', alpha=0.5, label='Mixture')  # Plot the mixture
        ax.plot(x, peak1, color=colors[i], linestyle='--', label=contributors[i])  # Plot individual contributions

    ax.set_xlabel("Fragment Size (bp)")
    ax.set_ylabel("Fluorescence Intensity")
    ax.set_title("Simulated Mixture Deconvolution")
    ax.legend()
    save_svg(fig, "mixture_deconvolution.svg")


def plot_electropherogram():
    """Simulates an electropherogram (DNA peak intensity plot)."""
    fig, ax = plt.subplots(figsize=(10, 4))
    x = np.linspace(100, 1000, 100)
    colors = ['purple', 'pink', 'black', 'red']
    bases = ['A', 'T', 'G', 'C']

    for i in range(4):
        peak = np.exp(-0.01 * (x - np.random.randint(200, 800)) ** 2) * np.random.uniform(5, 10)
        ax.plot(x, peak, color=colors[i], label=bases[i])

    ax.set_xlabel("Fragment Size (bp)")
    ax.set_ylabel("Fluorescence Intensity")
    ax.set_title("Simulated Electropherogram")
    ax.legend()
    save_svg(fig, "electropherogram.svg")


def plot_phylogenetic_tree(sequences):
    """Constructs a simple phylogenetic tree from DNA sequences."""
    names = list(sequences.keys())
    if not names:
        print("No sequences to build phylogenetic tree.")
        return

    matrix = []
    for i in range(len(names)):
        row = []
        for j in range(i + 1):  # Include the diagonal
            row.append(np.random.rand())
        matrix.append(row)

    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(DistanceMatrix(names, matrix))
    fig = plt.figure(figsize=(10, 5))
    ax = fig.add_subplot(111)
    Phylo.draw(tree, do_show=False, axes=ax)
    ax.set_title("Phylogenetic Tree")
    save_svg(fig, "phylogenetic_tree.svg")


def plot_similarity_heatmap(sequences):
    """Generates a heatmap of sequence similarity."""
    names = list(sequences.keys())
    if not names:
        print("No sequences to build similarity heatmap.")
        return

    similarity_matrix = np.zeros((len(names), len(names)))

    for i, key1 in enumerate(names):
        seq1 = str(sequences[key1][0].seq)
        for j, key2 in enumerate(names):
            seq2 = str(sequences[key2][0].seq)
            matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
            similarity_matrix[i, j] = matches / min(len(seq1), len(seq2))

    fig, ax = plt.subplots(figsize=(8, 6))
    sns.heatmap(similarity_matrix, annot=True, cmap='PiYG', xticklabels=names, yticklabels=names)
    ax.set_title("Sequence Similarity Heatmap")
    save_svg(fig, "similarity_heatmap.svg")


def plot_histogram_of_allele_frequencies():
    """Simulates a histogram of allele frequencies."""
    fig, ax = plt.subplots(figsize=(8, 5))
    allele_counts = np.random.poisson(20, 50)
    ax.hist(allele_counts, bins=15, color='pink', edgecolor='black', alpha=0.7)
    ax.set_xlabel("Allele Count")
    ax.set_ylabel("Frequency")
    ax.set_title("Histogram of Allele Frequencies")
    save_svg(fig, "allele_frequencies.svg")


def plot_probability_distribution():
    """Simulates probability distribution for DNA sequence classification."""
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.linspace(0, 1, 100)
    probs = [np.exp(-((x - mu) ** 2) / 0.02) for mu in [0.2, 0.5, 0.8]]

    for prob, color, label in zip(probs, ['pink', 'purple', 'red'], ['Group A', 'Group B', 'Group C']):
        ax.plot(x, prob, color=color, label=label)

    ax.set_xlabel("Probability")
    ax.set_ylabel("Density")
    ax.set_title("DNA Sequence Probability Graph")
    ax.legend()
    save_svg(fig, "probability_distribution.svg")


def plot_dot_plot(seq1, seq2):
    """Creates a dot plot comparison of two sequences."""
    fig, ax = plt.subplots(figsize=(8, 8))
    seq1, seq2 = seq1.upper(), seq2.upper()
    dot_matrix = np.zeros((len(seq1), len(seq2)))

    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if seq1[i] == seq2[j]:
                dot_matrix[i, j] = 1

    ax.imshow(dot_matrix, cmap='Blues', aspect='auto', origin='lower')
    ax.set_xlabel("Sequence 2")
    ax.set_ylabel("Sequence 1")
    ax.set_title("Dot Plot Comparison")
    save_svg(fig, "dot_plot.svg")

def plot_ideograms_combined(sequences):
    """Creates a combined ideogram plot for all suspects."""
    num_suspects = len(sequences)
    if num_suspects == 0:
        print("No suspects to plot ideograms.")
        return

    chromosome_length = 200  # Fixed chromosome length
    ideogram_height = 0.5  # Height of each ideogram (adjust as needed)
    gap = 0.2  # Gap between ideograms

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.set_xlim(0, chromosome_length)
    ax.set_ylim(-1, num_suspects * (ideogram_height + gap))
    ax.axis('off')

    suspect_names = list(sequences.keys())

    for i, suspect_name in enumerate(suspect_names):
        num_bands = np.random.randint(5, 15)
        bands = []
        for _ in range(num_bands):
            start = np.random.randint(0, chromosome_length - 20)
            end = start + np.random.randint(10, 30)
            color = np.random.choice(['lightgray', 'darkgray', 'gray'])
            bands.append((start, end, color))

        y_center = i * (ideogram_height + gap)
        for start, end, color in bands:
            ax.add_patch(patches.Rectangle((start, y_center), end - start, ideogram_height, color=color))

        ax.add_patch(patches.Rectangle((0, y_center), chromosome_length, ideogram_height, fill=False, edgecolor='black')) #Outline

        ax.text(-10, y_center + ideogram_height / 2, suspect_name, va='center', ha='right')  # Add suspect names to the left

    save_svg(fig, "combined_ideograms.svg")

# Example usage:
sequences = load_sequences("suspects")
plot_mixture_deconvolution()
plot_electropherogram()
if sequences:
    plot_phylogenetic_tree(sequences)
    plot_similarity_heatmap(sequences)
    plot_ideograms_combined(sequences)  # Call the combined ideogram plotting function
plot_histogram_of_allele_frequencies()
plot_probability_distribution()
plot_dot_plot("AGCTTAGCTA", "AGCTAGCTAA")