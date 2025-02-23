from align_sequences import align_sequences
from load_sequences import load_fasta_sequences
from Bio.pairwise2 import format_alignment
import os

def analyze_results(alignments, filename):
    """
    Analyze alignment results and determine similarity score.

    Parameters:
    - alignments (list): List of alignments.
    - filename (str): Filename of the suspect's DNA sequence.

    Returns:
    - dict: Filename and similarity score.
    """
    if not alignments: # Handle cases where there are no alignments
        return {"filename": filename, "similarity_score": 0.0}  # Or another appropriate default

    best_alignment = alignments[0]
    similarity_score = best_alignment[2] / best_alignment[4] if best_alignment[4] else 0.0 # Avoid division by zero
    return {"filename": filename, "similarity_score": similarity_score}


def process_dna_alignment(crime_scene_file="crime_scene_dna.fasta", suspects_dir="dataset"):
    """
    Processes DNA alignment by comparing a crime scene sequence to suspects' sequences.

    Args:
        crime_scene_file (str): Filename of the crime scene FASTA file.
        suspects_dir (str): Directory containing the suspects' FASTA files.
    """
    try:
        fasta_files = [f for f in os.listdir(suspects_dir) if f.endswith(".fasta")]
        sequences = load_fasta_sequences(suspects_dir, fasta_files + [crime_scene_file])

        crime_scene_seq = sequences.get(crime_scene_file)
        if not crime_scene_seq:
            raise ValueError(f"Crime scene sequence file '{crime_scene_file}' not found.")

        results = []
        for filename, suspect_seq in sequences.items():
            if filename != crime_scene_file:
                alignments = align_sequences(crime_scene_seq, suspect_seq)
                result = analyze_results(alignments, filename)
                results.append(result)

                print(f"\nAlignment between Crime Scene and {filename}:")
                if alignments: # Print only if there are alignments
                    print(format_alignment(*alignments[0]))
                else:
                    print("No alignments found.")
                print(f"Similarity Score: {result['similarity_score']:.2f}")

    except FileNotFoundError:
        print(f"Error: Directory '{suspects_dir}' not found.")
    except ValueError as e:
        print(f"Error: {e}")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")

# Example usage (remove or comment out when using as a module):
# if __name__ == "__main__":
#     process_dna_alignment() # Uses default file names and directory
#     # process_dna_alignment(crime_scene_file="my_crime_scene.fasta", suspects_dir="my_suspects")