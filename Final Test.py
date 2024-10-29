import random
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
import matplotlib.patches as mpatches

# DNA Sequence Generation
def generate_dna(sample_size, length):
    return [''.join(random.choices('ATCG', k=length)) for _ in range(sample_size)]

# DNA to RNA Conversion
def dna_to_rna(dna_sequence):
    return dna_sequence.replace('T', 'U')

# RNA to Amino Acids Translation
codon_table = {
    'UUU': 'F', 'UUC': 'F', 'UUA': 'L', 'UUG': 'L', 'CUU': 'L', 'CUC': 'L', 'CUA': 'L', 'CUG': 'L',
    'AUU': 'I', 'AUC': 'I', 'AUA': 'I', 'AUG': 'START', 'GUU': 'V', 'GUC': 'V', 'GUA': 'V', 'GUG': 'V',
    'UCU': 'S', 'UCC': 'S', 'UCA': 'S', 'UCG': 'S', 'CCU': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACU': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCU': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'UAU': 'Y', 'UAC': 'Y', 'UAA': 'STOP', 'UAG': 'STOP', 'CAU': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAU': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAU': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'UGU': 'C', 'UGC': 'C', 'UGA': 'STOP', 'UGG': 'W', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGU': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGU': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'
}

def rna_to_amino_acids(rna_sequence):
    return ''.join([codon_table.get(rna_sequence[i:i + 3], '-') for i in range(0, len(rna_sequence), 3)])

# Similarity Calculation
def calculate_similarity(sequence1, sequence2):
    if not sequence1 or not sequence2:  # Check for empty sequences
        return 0
    return sum(1 for a, b in zip(sequence1, sequence2) if a == b) / len(sequence1)

# Function to find a different color
def get_different_color(color, colors, index):
    for i in range(len(colors)):
        if i != index:  # Ensure it's not the same sample
            return colors[i]  # Use the actual color value
    return color  # Fallback, shouldn't happen if samples > 1

# Circular Plot for Comparisons
def plot_circular_comparison(ax, sequences, sequence_type, colors, searched_strings, num_samples):
    theta = np.linspace(0, 2 * np.pi, num_samples, endpoint=False)  # Change this line

    # Calculate a dynamic font size based on the number of samples
    base_font_size = 12  # Base font size for a small number of samples
    font_size = max(6, base_font_size - (num_samples // 2))  # Scale font size down

    # Outer Ring: Match Frequencies (placed at bottom for better visibility)
    for i in range(num_samples):
        for j in range(num_samples):
            if i != j:  # Prevent comparing the same sample
                similarity = calculate_similarity(sequences[i], sequences[j])
                ax.bar(theta[i], similarity, width=(2 * np.pi / num_samples), bottom=2.75,
                       color=colors[i], alpha=0.7)

    # Middle Ring: Most Similar Sample Highlight
    #'''
    for i in range(num_samples):
        similarities = [calculate_similarity(sequences[i], sequences[j]) for j in range(num_samples)]
        most_similar_index = np.argmax(similarities)  # Index of the most similar sample
        most_similar_value = similarities[most_similar_index]  # Value of that similarity

        middle_color = get_different_color(colors[most_similar_index], colors, most_similar_index)
        ax.bar(theta[i], most_similar_value, width=(2 * np.pi / num_samples), bottom=1.5,
               color=colors[most_similar_index], alpha=0.8)
    #'''

    # Inner Ring: Searched String Frequencies
    for i, sample in enumerate(sequences):
        freq = sum(sample.count(s) for s in searched_strings) / (len(sample) / 10)
        ax.bar(theta[i], freq, width=(2 * np.pi / num_samples), bottom=0.5,
               color=colors[i], alpha=0.9)

    # Set labels and title
    if num_samples <= 10:
        ax.set_xticks(theta)
        ax.set_xticklabels([f"Sample {i + 1}" for i in range(num_samples)])
        ax.set_yticklabels([])  # Remove radial ticks
        ax.set_yticks([])
    else:
        ax.set_xticks(theta)
        ax.set_xticklabels([f"Sample {i + 1}" for i in range(num_samples)], fontsize=font_size)
        ax.set_yticklabels([])  # Remove radial ticks
        ax.set_yticks([])

    ax.set_title(f"{sequence_type} Sequence Comparison", va='bottom')

    # Create a legend for this subplot
    if num_samples <= 10:
        patches = [mpatches.Patch(color=colors[i], label=f"Sample {i + 1}") for i in range(num_samples)]
        ax.legend(handles=patches, loc='upper right', title="Samples")

# Customize
# Generate DNA, RNA, and Amino Acid Samples
num_samples = 10  # Change this to test with different numbers of samples
dna_samples = generate_dna(num_samples, 1000)  # DNA samples, each of length 100
rna_samples = [dna_to_rna(sample) for sample in dna_samples]
amino_acid_samples = [rna_to_amino_acids(sample) for sample in rna_samples]
searched_strings = ['TTC', 'UAUU', 'SCY', 'STOP', 'FS']  # Ensure these are RNA codons
'''
    String pre sets
    'TGA', 'TAG', 'TAA', 'ATG', 'AUG', 'UGA', 'UAG','UAA', 'STOP', 'START'
    'ATGG', 'AAUGGA', 'UUAT', 'CCCCC', 'CA', 'FF', 'SC'
    'TTC', 'UAUU', 'SCY', 'STOP', 'FS'
    'UAG', 'UCG'
'''


# Plot the comparison for DNA, RNA, and Amino Acids
fig, axs = plt.subplots(1, 3, subplot_kw={'projection': 'polar'}, figsize=(20, 8))

# Update colormap retrieval
colors = plt.colormaps['tab10']  # Get a colormap with enough distinct colors
color_list = [colors(i) for i in range(colors.N)]  # Create a list of colors


# If more samples are requested than available colors, cycle through the colors
if num_samples > len(color_list):
    color_list = [color_list[i % len(color_list)] for i in range(num_samples)]

# Plot DNA
plot_circular_comparison(axs[0], dna_samples, "DNA", color_list, searched_strings, num_samples)

# Plot RNA
plot_circular_comparison(axs[1], rna_samples, "RNA", color_list, searched_strings, num_samples)

# Plot Amino Acids
plot_circular_comparison(axs[2], amino_acid_samples, "Amino Acids", color_list, searched_strings, num_samples)

# Display the plot
plt.tight_layout()
plt.show()
