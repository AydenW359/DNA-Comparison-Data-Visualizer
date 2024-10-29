# DNA-Comparison-Data-Visualizer

This project generates DNA samples, converts them to RNA, translates them to amino acid sequences, and visualizes the similarity and searched patterns across these sequences in circular plots.

## Features
- **DNA Sample Generation**: Creates a specified number of DNA samples with randomized sequences of ATCG nucleotides.
- **DNA to RNA Conversion**: Transforms each DNA sample into RNA sequences by substituting thymine (T) with uracil (U).
- **RNA to Amino Acid Translation**: Translates RNA sequences into amino acid sequences based on the standard genetic codon table.
- **Circular Comparison Plot**: Visualizes comparisons between samples in three layers:
  - **Outer Layer**: Similarity scores between samples.
  - **Middle Layer**: Color of sample.
  - **Inner Layer**: Frequency of specific searched-for strings within each sample.

## Getting Started
### Prerequisites
- **Python 3.x**
- **Required Libraries**: `random`, `matplotlib`, `numpy`

### Installation
1. Clone this repository.
2. Install the required libraries:
   ```bash
   pip install matplotlib numpy
