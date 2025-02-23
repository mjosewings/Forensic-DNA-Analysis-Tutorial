from analyze_results import analyze_results
from align_sequences import align_sequences
from load_sequences import load_fasta_sequences
from Bio.pairwise2 import format_alignment
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os


def generate_report(results, output_file="alignment_report.txt"):
    """
    Generate a report summarizing the alignment results.

    Parameters:
    - results (list): List of dictionaries containing filenames and similarity scores.
    - output_file (str): Path to save the report.
    """
    with open(output_file, "w") as f:
        f.write("Crime Scene DNA Comparison Report\n")
        f.write("=================================\n\n")
        for result in results:
            f.write(f"Suspect File: {result['filename']}\n")
            f.write(f"Similarity Score: {result['similarity_score']:.2f}\n")
            f.write("---------------------------------\n")
    print(f"Report generated: {output_file}")


def generate_seaborn_graph(results, output_file="graphs/dna_similarity_graph.png"):
    """
    Generate a bar graph visualizing the similarity scores using Seaborn.

    Parameters:
    - results (list): List of dictionaries containing filenames and similarity scores.
    - output_file (str): Path to save the graph.
    """
    # Prepare data for the DataFrame
    data = {
        'Suspect': [result['filename'].replace('.fasta', '') for result in results],
        'similarity_score': [result['similarity_score'] * 100 for result in results],  # Convert to percentage
        'Match': ['MATCH' if result['similarity_score'] >= 0.9 else 'NO MATCH' for result in results]
    }

    df = pd.DataFrame(data)

    # Assign colors based on match status
    color_palette = {'MATCH': 'green', 'NO MATCH': 'red'}

    # Plotting with Seaborn
    plt.figure(figsize=(10, 6))
    sns.barplot(x='Suspect', y='similarity_score', hue='Match', data=df, palette=color_palette, dodge=False)
    plt.axhline(y=90, color='blue', linestyle='--', label="90% Match Threshold")

    plt.title("DNA Similarity Scores with Crime Scene")
    plt.xlabel("Suspects")
    plt.ylabel("Similarity Score (%)")
    plt.legend(title="Match Status")
    plt.tight_layout()

    plt.savefig(output_file)
    plt.show()


def process_dna_analysis(crime_scene_file="crime_scene_dna.fasta", suspects_dir="suspects"):
    """
    Processes DNA analysis by aligning crime scene DNA with suspects' DNA from FASTA files in a directory.

    Args:
        crime_scene_file (str): The filename of the crime scene DNA FASTA file.
        suspects_dir (str): The directory containing the suspects' DNA FASTA files.
    """
    try:
        fasta_filenames = [f for f in os.listdir(suspects_dir) if f.endswith(".fasta")]
        sequences = load_fasta_sequences(suspects_dir, fasta_filenames + [crime_scene_file])  # Include crime scene file

        crime_scene_seq = sequences.get(crime_scene_file)
        if not crime_scene_seq:
            raise ValueError(f"Crime scene sequence file '{crime_scene_file}' not found.")

        results = []
        for filename, suspect_seq in sequences.items():
            if filename != crime_scene_file:  # Compare with all suspects
                alignments = align_sequences(crime_scene_seq, suspect_seq)
                result = analyze_results(alignments, filename)
                results.append(result)

                print(f"\nAlignment between Crime Scene and {filename}:")
                print(format_alignment(*alignments[0]))
                print(f"Similarity Score: {result['similarity_score']:.2f}")

        generate_report(results)
        generate_seaborn_graph(results)

    except FileNotFoundError:
        print(f"Error: Directory '{suspects_dir}' not found.")
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Example of how to use the function:
# process_dna_analysis() # Uses default filenames and directory
# process_dna_analysis(crime_scene_file="my_crime_scene.fasta", suspects_dir="dna_profiles") # Specify files