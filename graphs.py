from Bio import SeqIO
import plotly.graph_objects as go
import os

def parse_fasta(fasta_files):
    """Parses FASTA files and returns a list of sequence dictionaries."""
    sequences = []
    for fasta_file in fasta_files:
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append({
                'id': record.id,
                'sequence': str(record.seq),
                'length': len(record.seq)
            })
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


def process_fasta_visualization(fasta_dir="dataset"):
    """Processes FASTA files in a directory and creates ideograms."""
    try:
        fasta_files = [os.path.join(fasta_dir, f) for f in os.listdir(fasta_dir) if f.endswith(".fasta")]
        sequences = parse_fasta(fasta_files)
        create_separate_ideograms(sequences)

    except FileNotFoundError:
        print(f"Error: Directory '{fasta_dir}' not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


# Example usage:
# process_fasta_visualization()  # Uses default directory "dataset"
# process_fasta_visualization(fasta_dir="my_fasta_files") # Specify directory