# TaqMan PCR Primer and Probe Design Tool

A Python command-line tool that replicates key functionality from Primer Express, specifically for TaqMan-based PCR primer and probe design.

## Features

- Accepts DNA sequence input in FASTA format
- Cleans and parses the FASTA input
- Generates three primer pairs (forward and reverse)
- For each primer pair:
  - Selects a valid TaqMan probe following specific constraints
  - Calculates GC content
  - Simulates NCBI BLAST requests
  - Generates GC content plots

## Installation

1. Clone this repository:

```bash
git clone <repository-url>
cd <repository-directory>
```

2. Create a virtual environment and activate it:

```bash
python -m venv venv
source venv/bin/activate  # On Windows, use: venv\Scripts\activate
```

3. Install the required dependencies:

```bash
pip install -r requirements.txt
```

## Usage

The tool can be used in two ways:

### 1. Using a FASTA file:

```bash
python primer_express.py -f input.fasta -o results.json -p plots_directory
```

### 2. Providing a FASTA sequence directly:

```bash
python primer_express.py -s ">Sequence1
ATGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCT
AGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
```

### Command-line options:

- `-f, --file`: Input FASTA file path
- `-s, --sequence`: Input DNA sequence string in FASTA format
- `-o, --output`: Output JSON file path (default: output.json)
- `-p, --plot-dir`: Directory to save GC plots (default: plots)

## Output

The tool generates:

1. A JSON file with the primer pairs, TaqMan probes, and GC content data
2. GC content plots for each amplicon region and the entire sequence

Example JSON output:

```json
{
  "primers": [
    {
      "forward": "ATGCGT...",
      "reverse": "TACGAT...",
      "forward_tm": 60.3,
      "reverse_tm": 59.9,
      "forward_gc": 50.0,
      "reverse_gc": 48.0,
      "penalty": 1.23,
      "probe": {
        "sequence": "ACGTAGCTAGCTAGCTA",
        "gc": 52.6,
        "tm": 69.0
      },
      "blast_forward": "Submitted to NCBI BLAST",
      "blast_reverse": "Submitted to NCBI BLAST",
      "gc_plot": "plots/gc_plot_pair_1.png"
    }
  ],
  "gc_plot": "plots/gc_plot_overall.png"
}
```

## Requirements

- Python 3.6+
- Biopython
- primer3-py
- matplotlib
- requests

## License

[MIT License](LICENSE) 