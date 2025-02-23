import os
from Bio import SeqIO

def load_fasta_sequences(data_dir, fasta_files):
    """
    Load multiple FASTA files and return a dictionary of sequences.

    Parameters:
    - data_dir (str): Directory containing the FASTA files.
    - fasta_files (list): List of FASTA filenames to load.

    Returns:
    - dict: Keys are filenames, values are sequence strings.
    """
    sequences = {}

    for file in fasta_files:
        path = os.path.join(data_dir, file)

        try:
            with open(path, "r") as f:
                for record in SeqIO.parse(f, "fasta"):
                    sequences[file] = str(record.seq)  # Store sequence as string
        except FileNotFoundError:
            print(f"Warning: {file} not found in {data_dir}.")
        except Exception as e: # Catch other potential Bio.SeqIO parse errors
            print(f"Error parsing {file}: {e}")

    return sequences


def process_fasta_loading(data_dir="dataset", fasta_files=None):
    """Loads FASTA sequences from a directory and prints summaries.

    Args:
        data_dir (str): The directory containing the FASTA files (default: "dataset").
        fasta_files (list, optional): A list of specific FASTA files to load.
                                     If None, all .fasta files in the directory are loaded.
    """
    try:
        if fasta_files is None:  # Load all fasta files in the directory if file list not given
            fasta_files = [f for f in os.listdir(data_dir) if f.endswith(".fasta")]

        sequences = load_fasta_sequences(data_dir, fasta_files)

        for filename, sequence in sequences.items():
            print(f"{filename}: {len(sequence)} bases loaded.")

    except FileNotFoundError:
        print(f"Error: Directory '{data_dir}' not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


# Example usage (remove or comment out when using as a module):
# if __name__ == "__main__":
#     process_fasta_loading()  # Uses default directory and loads all .fasta files
#     # Or specify the directory and a list of files:
#     # process_fasta_loading(data_dir="my_data", fasta_files=["seq1.fasta", "seq2.fasta"])