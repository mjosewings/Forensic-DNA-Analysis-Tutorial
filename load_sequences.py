import os
from Bio import SeqIO

def load_fasta_sequences(data_dir, fasta_files):
    """
    This function reads FASTA files and extracts the DNA sequences.

    It takes the directory where the files are located and a list of filenames as input.
    It returns a dictionary where the keys are the filenames and the values are the DNA sequences.
    """
    sequences = {}  # Create an empty dictionary to store the sequences

    for file in fasta_files:  # Loop through each filename
        path = os.path.join(data_dir, file)  # Create the full path to the file

        try:  # Try to open and read the file
            with open(path, "r") as f:  # Open the file for reading
                for record in SeqIO.parse(f, "fasta"):  # Use Bio.SeqIO to parse the FASTA file
                    sequences[file] = str(record.seq)  # Extract the sequence and store it in the dictionary
        except FileNotFoundError:  # If the file is not found
            print(f"Warning: {file} not found in {data_dir}.")
        except Exception as e:  # If there's another error during parsing
            print(f"Error parsing {file}: {e}")

    return sequences  # Return the dictionary of sequences


def process_fasta_loading(data_dir="dataset", fasta_files=None):
    """
    This function loads FASTA sequences and prints a summary.

    It takes the directory containing the FASTA files and an optional list of filenames.
    If no filenames are provided, it loads all .fasta files in the directory.
    It then prints the number of bases (characters) in each sequence.
    """
    try:
        if fasta_files is None:  # If no list of files is given
            fasta_files = [f for f in os.listdir(data_dir) if f.endswith(".fasta")]  # Find all .fasta files

        sequences = load_fasta_sequences(data_dir, fasta_files)  # Load the sequences

        for filename, sequence in sequences.items():  # Loop through the sequences
            print(f"{filename}: {len(sequence)} bases loaded.")  # Print the filename and sequence length

    except FileNotFoundError:  # If the directory is not found
        print(f"Error: Directory '{data_dir}' not found.")
    except Exception as e:  # If another error occurs
        print(f"An unexpected error occurred: {e}")


# Example usage (now always runs when the script is executed):
process_fasta_loading()  # Uses default directory "dataset" and loads all .fasta files in it.  Make sure you have a folder called 'dataset' with .fasta files in it.
# To use a different directory and specific files, uncomment and modify the lines below:
# process_fasta_loading(data_dir="my_data", fasta_files=["seq1.fasta", "seq2.fasta"])  # Replace with your directory and filenames. Make sure the files are in the directory.
