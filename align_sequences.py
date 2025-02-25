from Bio import pairwise2
from Bio.pairwise2 import format_alignment
from load_sequences import load_fasta_sequences
import os

def align_sequences(seq1, seq2):
    """
    Performs global sequence alignment using the Needleman-Wunsch algorithm.

    Args:
        seq1 (str): The first DNA sequence.
        seq2 (str): The second DNA sequence.

    Returns:
        list: A list of top alignments.
    """
    alignments = pairwise2.align.globalxx(seq1, seq2)
    return alignments

def process_and_print_alignments(data_directory="suspects", crime_scene_file="crime_scene_dna.fasta"):
    """
    Loads FASTA sequences from the 'suspects' directory, aligns them to a crime scene sequence, and prints the alignments.

    Args:
        data_directory (str): The directory containing the FASTA files (defaults to "suspects/").
        crime_scene_file (str): The filename of the crime scene DNA FASTA file.
    """
    try:
        fasta_filenames = [f for f in os.listdir(data_directory) if f.endswith(".fasta")]

        sequences = load_fasta_sequences(data_directory, fasta_filenames + [crime_scene_file]) # Ensure crime scene is loaded.

        crime_scene_seq = sequences.get(crime_scene_file)
        if not crime_scene_seq:
            raise ValueError(f"Error: Crime scene sequence '{crime_scene_file}' not found!")

        for filename, suspect_seq in sequences.items():
            if filename != crime_scene_file:
                alignments = align_sequences(crime_scene_seq, suspect_seq)

                print(f"\nAlignment between Crime Scene and {filename}:")
                if alignments:  # Check if alignments were found
                    print(format_alignment(*alignments[0]))
                else:
                    print("No alignments found.")

    except FileNotFoundError:
        print(f"Error: Directory '{data_directory}' not found.")
    except ValueError as e:
        print(e)
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Example usage:
process_and_print_alignments()
