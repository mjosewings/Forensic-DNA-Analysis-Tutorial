import random

def generate_dna_sequence(length=100):
    """Generates a random DNA sequence of a given length."""
    return ''.join(random.choices('ATCG', k=length))

def create_fasta_entry(identifier, sequence):
    """Formats a DNA sequence as a FASTA entry."""
    return f'>{identifier}\n{sequence}'

def generate_crime_scene_sample():
    """Generates a simulated crime scene DNA sample."""
    return create_fasta_entry('Crime_Scene_Sample', generate_dna_sequence(200))

def generate_suspect_samples(num_suspects=3):
    """Generates simulated DNA samples for a given number of suspects."""
    return [create_fasta_entry(f'Suspect_{i+1}', generate_dna_sequence(200)) for i in range(num_suspects)]

crime_scene_dna = generate_crime_scene_sample()
suspect_dna_samples = generate_suspect_samples()

with open('suspects/crime_scene_dna.fasta', 'w') as file:
    file.write(crime_scene_dna + '\n')
    file.write('\n'.join(suspect_dna_samples))

print("Crime scene and suspect DNA samples have been generated and saved in 'crime_scene_dna.fasta'")
