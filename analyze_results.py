from align_sequences import align_sequences
from load_sequences import load_fasta_sequences
from Bio.pairwise2 import format_alignment
import os

def analyze_results(alignments, filename):
    """
    Analyze alignment results and determine similarity score.

    Args:
        alignments (list): List of alignments.
        filename (str): Filename of the suspect's DNA sequence.

    Returns:
        dict: Filename and similarity score.
    """
    if not alignments: #Handle the case where there are no alignments
        return {"filename": filename, "similarity_score": 0.0}

    best_alignment = alignments[0]
    similarity_score = best_alignment[2] / best_alignment[4]  # Matches / Alignment length
    return {"filename": filename, "similarity_score": similarity_score}

def process_and_analyze(data_directory="suspects/", crime_scene_file="crime_scene_dna.fasta"):
    """
    Loads sequences from the 'suspects' folder, aligns them to the crime scene, and analyzes the results.

    Args:
        data_directory (str): The directory containing the suspect FASTA files.
        crime_scene_file (str): The filename of the crime scene DNA FASTA file.
    """
    try:
        fasta_filenames = [f for f in os.listdir(data_directory) if f.endswith(".fasta")]
        sequences = load_fasta_sequences(data_directory, fasta_filenames + [crime_scene_file])

        crime_scene_seq = sequences.get(crime_scene_file)
        if not crime_scene_seq:
            raise ValueError(f"Error: Crime scene sequence '{crime_scene_file}' not found!")

        results = []
        for filename, suspect_seq in sequences.items():
            if filename != crime_scene_file:
                alignments = align_sequences(crime_scene_seq, suspect_seq)
                result = analyze_results(alignments, filename)
                results.append(result)

                print(f"\nAlignment between Crime Scene and {filename}:")
                if alignments:
                    print(format_alignment(*alignments[0]))
                else:
                    print ("No Alignments found")
                print(f"Similarity Score: {result['similarity_score']:.2f}")

    except FileNotFoundError:
        print(f"Error: Directory '{data_directory}' not found.")
    except ValueError as e:
        print(e)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Example usage (will run when the script is executed):
process_and_analyze()