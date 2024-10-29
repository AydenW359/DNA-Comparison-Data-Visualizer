import random
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


# Generate DNA sequences
def generate_dna(sample_size, length):
    return [''.join(random.choices('ATCG', k=length)) for _ in range(sample_size)]


# DNA to RNA Conversion
def dna_to_rna(dna_sequence):
    return dna_sequence.replace('T', 'U')


# Codon table for RNA to Amino Acids translation
codon_table = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'M', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': 'Stop', 'UAG': 'Stop', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': 'Stop', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}


# RNA to Amino Acids Conversion
def rna_to_amino_acids(rna_sequence):
    return ''.join([codon_table.get(rna_sequence[i:i + 3], '-') for i in range(0, len(rna_sequence), 3)])


# Function to find subsequence matches in a sequence
def find_subsequence_positions(sequence, subsequence):
    positions = []
    sub_len = len(subsequence)
    for i in range(len(sequence) - sub_len + 1):
        if sequence[i:i + sub_len] == subsequence:
            positions.append(i)
    return positions


# Function to generate Circos-like circular plot with connection lines
def plot_circos_comparison_with_connections(samples, subsequences, seq_type):
    num_samples = len(samples)
    angles = np.linspace(0, 2 * np.pi, num_samples, endpoint=False).tolist()
    radii = 10  # Static radii for the circular layout

    fig, ax = plt.subplots(figsize=(8, 8), subplot_kw={'projection': 'polar'})

    # Plot each sample as an arc around the circle
    for i, (sample, angle) in enumerate(zip(samples, angles)):
        # Get positions of each subsequence in the current sample
        for subseq in subsequences:
            positions = find_subsequence_positions(sample, subseq)
            # Plot each subsequence match as a point on the arc
            for pos in positions:
                ax.plot(angle, radii, 'o', label=f'{subseq} in Sample {i + 1}')

    # Draw connection lines between matching subsequences in different samples
    for i, sample1 in enumerate(samples):
        for j, sample2 in enumerate(samples[i + 1:], start=i + 1):
            angle1 = angles[i]
            angle2 = angles[j]
            for subseq in subsequences:
                # Find matching subsequences in both samples
                positions1 = find_subsequence_positions(sample1, subseq)
                positions2 = find_subsequence_positions(sample2, subseq)

                # Draw lines connecting the matching subsequences
                for pos1 in positions1:
                    for pos2 in positions2:
                        # Add a connecting line between positions in the two arcs
                        ax.plot([angle1, angle2], [radii, radii], color='blue', lw=1, alpha=0.7)

    # Customize the circular plot
    ax.set_title(f"Circos-like Plot with Connections: {seq_type} Subsequence Comparison", pad=20)
    ax.set_xticks(angles)
    ax.set_xticklabels([f"Sample {i + 1}" for i in range(num_samples)], fontsize=10)
    ax.set_yticklabels([])  # Hide radial ticks

    plt.show()


# Example DNA sequences
dna_samples = generate_dna(5, 100)  # 5 DNA samples, each of length 100
rna_samples = [dna_to_rna(sample) for sample in dna_samples]
amino_acid_samples = [rna_to_amino_acids(sample) for sample in rna_samples]

# Define subsequences of interest
subsequences_dna = ['ATG', 'CGT']  # Example DNA subsequences
subsequences_rna = ['AUG', 'CGU']  # Corresponding RNA subsequences
subsequences_amino_acids = ['M', 'R']  # Corresponding amino acid subsequences

# Plot Circos-like comparisons for DNA, RNA, and Amino Acids with connections
plot_circos_comparison_with_connections(dna_samples, subsequences_dna, "DNA")
plot_circos_comparison_with_connections(rna_samples, subsequences_rna, "RNA")
plot_circos_comparison_with_connections(amino_acid_samples, subsequences_amino_acids, "Amino Acids")
