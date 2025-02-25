from analyze_results import analyze_results
from align_sequences import align_sequences
from load_sequences import load_fasta_sequences
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

REPORT_FILE = "alignment_report.txt"
GRAPH_FILE = "graphs/dna_similarity_graph.png"
MATCH_THRESHOLD = 0.9

def generate_report(results, output_file=REPORT_FILE):
    """Generates a report summarizing alignment results."""
    with open(output_file, "w") as f:
        f.write("Crime Scene DNA Comparison Report\n")
        f.write("=================================\n\n")
        for result in results:
            f.write(f"Suspect File: {result['filename']}\n")
            f.write(f"Similarity Score: {result['similarity_score']:.2f}\n")
            f.write("---------------------------------\n")
    print(f"Report generated: {output_file}")


def generate_seaborn_graph(results, output_file=GRAPH_FILE):
    """Generates a bar graph visualizing similarity scores."""
    data = {
        'Suspect': [result['filename'].replace('.fasta', '') for result in results],
        'Similarity (%)': [result['similarity_score'] * 100 for result in results],
        'Match': ['MATCH' if result['similarity_score'] >= MATCH_THRESHOLD else 'NO MATCH' for result in results]
    }
    df = pd.DataFrame(data)

    color_palette = {'MATCH': 'green', 'NO MATCH': 'red'}

    plt.figure(figsize=(10, 6))
    sns.barplot(x='Suspect', y='Similarity (%)', hue='Match', data=df, palette=color_palette, dodge=False)
    plt.axhline(y=MATCH_THRESHOLD * 100, color='blue', linestyle='--', label=f"{MATCH_THRESHOLD * 100}% Match Threshold")

    plt.title("DNA Similarity Scores with Crime Scene")
    plt.xlabel("Suspects")
    plt.ylabel("Similarity (%)")
    plt.legend(title="Match Status")
    plt.tight_layout()

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file)
    plt.show()

def process_dna_analysis(data_directory="suspects", crime_scene_file="crime_scene_dna.fasta"):
    """
    Processes DNA analysis by aligning crime scene DNA with suspects' DNA from FASTA files in a directory, generates a report, and creates a graph.

    Args:
        data_directory (str): The directory containing the suspects' DNA FASTA files (defaults to "suspects/").
        crime_scene_file (str): The filename of the crime scene DNA FASTA file.
    """
    try:
        fasta_filenames = [f for f in os.listdir(data_directory) if f.endswith(".fasta")]
        sequences = load_fasta_sequences(data_directory, fasta_filenames + [crime_scene_file])

        crime_scene_seq = sequences.get(crime_scene_file)
        if not crime_scene_seq:
            raise ValueError(f"Error: Crime scene sequence file '{crime_scene_file}' not found.")

        results = []
        for filename, suspect_seq in sequences.items():
            if filename != crime_scene_file:
                alignments = align_sequences(crime_scene_seq, suspect_seq)
                result = analyze_results(alignments, filename)
                results.append(result)

                print(f"\nAlignment between Crime Scene and {filename}:")
                print(format_alignment(*alignments[0]))
                print(f"Similarity Score: {result['similarity_score']:.2f}")

        generate_report(results)
        generate_seaborn_graph(results)

    except FileNotFoundError:
        print(f"Error: Directory '{data_directory}' not found.")
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Example usage:
process_dna_analysis()