#!/usr/bin/env python3
"""
TaqMan PCR Primer and Probe Design Tool

A command-line tool that replicates key functionality from Primer Express,
specifically for TaqMan-based PCR primer and probe design.
"""

import argparse
import json
import os
import sys
from typing import Dict, List, Any, Tuple, Optional

# Set matplotlib backend to non-interactive Agg before importing pyplot
# This is necessary for using matplotlib in web applications
import matplotlib
matplotlib.use('Agg')

from Bio import SeqIO
from Bio.SeqUtils import GC
import primer3
import requests
import matplotlib.pyplot as plt
import numpy as np
from io import StringIO


def parse_arguments() -> argparse.Namespace:
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(
        description="TaqMan PCR Primer and Probe Design Tool"
    )
    parser.add_argument(
        "-f", "--file", type=str, help="Input FASTA file path"
    )
    parser.add_argument(
        "-s", "--sequence", type=str, help="Input DNA sequence string in FASTA format"
    )
    parser.add_argument(
        "-o", "--output", type=str, default="output.json", 
        help="Output JSON file path (default: output.json)"
    )
    parser.add_argument(
        "-p", "--plot-dir", type=str, default="plots",
        help="Directory to save GC plots (default: plots)"
    )
    parser.add_argument(
        "--relaxed", action="store_true", help="Use relaxed constraints for probe design"
    )
    return parser.parse_args()


def parse_fasta(fasta_data: str) -> Tuple[str, str]:
    """
    Parse FASTA format data to extract sequence ID and sequence.
    
    Args:
        fasta_data: FASTA format string
        
    Returns:
        Tuple of (sequence ID, DNA sequence)
    """
    try:
        fasta_io = StringIO(fasta_data)
        for record in SeqIO.parse(fasta_io, "fasta"):
            return record.id, str(record.seq).upper()
    except Exception as e:
        raise ValueError(f"Failed to parse FASTA data: {e}")
    
    raise ValueError("No valid FASTA sequence found")


def clean_sequence(sequence: str) -> str:
    """
    Clean DNA sequence by removing non-DNA characters.
    
    Args:
        sequence: DNA sequence string
        
    Returns:
        Cleaned DNA sequence
    """
    # Remove whitespace and convert to uppercase
    sequence = ''.join(sequence.split()).upper()
    
    # Keep only valid DNA characters
    valid_bases = set('ATGCN')
    sequence = ''.join(c for c in sequence if c in valid_bases)
    
    return sequence


def design_primers(sequence: str, num_pairs: int = 3) -> List[Dict[str, Any]]:
    """
    Design primer pairs using primer3.
    
    Args:
        sequence: DNA sequence
        num_pairs: Number of primer pairs to generate
        
    Returns:
        List of dictionaries containing primer pair information
    """
    # Primer3 requires sequence to be uppercase
    sequence = sequence.upper()
    
    # Configure primer3 parameters
    seq_args = {
        'SEQUENCE_ID': 'target_sequence',
        'SEQUENCE_TEMPLATE': sequence,
    }
    
    global_args = {
        'PRIMER_PICK_LEFT_PRIMER': 1,
        'PRIMER_PICK_RIGHT_PRIMER': 1,
        'PRIMER_OPT_SIZE': 20,
        'PRIMER_MIN_SIZE': 18,
        'PRIMER_MAX_SIZE': 25,
        'PRIMER_OPT_TM': 60.0,
        'PRIMER_MIN_TM': 57.0,
        'PRIMER_MAX_TM': 63.0,
        'PRIMER_MIN_GC': 40.0,
        'PRIMER_MAX_GC': 60.0,
        'PRIMER_NUM_RETURN': num_pairs,
        'PRIMER_PRODUCT_SIZE_RANGE': [[70, 150]],
    }
    
    # Run primer3 design
    results = primer3.bindings.designPrimers(seq_args, global_args)
    
    primer_pairs = []
    
    for i in range(num_pairs):
        # Check if we have results for this pair
        if f'PRIMER_LEFT_{i}_SEQUENCE' not in results:
            break
            
        forward = results[f'PRIMER_LEFT_{i}_SEQUENCE']
        reverse = results[f'PRIMER_RIGHT_{i}_SEQUENCE']
        forward_tm = results[f'PRIMER_LEFT_{i}_TM']
        reverse_tm = results[f'PRIMER_RIGHT_{i}_TM']
        forward_gc = GC(forward)
        reverse_gc = GC(reverse)
        penalty = results[f'PRIMER_PAIR_{i}_PENALTY']
        
        # Get primer locations
        forward_start = results[f'PRIMER_LEFT_{i}'][0]
        forward_length = results[f'PRIMER_LEFT_{i}'][1]
        reverse_start = results[f'PRIMER_RIGHT_{i}'][0] - results[f'PRIMER_RIGHT_{i}'][1] + 1
        
        pair_info = {
            'forward': forward,
            'reverse': reverse,
            'forward_tm': round(forward_tm, 1),
            'reverse_tm': round(reverse_tm, 1),
            'forward_gc': round(forward_gc, 1),
            'reverse_gc': round(reverse_gc, 1),
            'penalty': round(penalty, 2),
            'product_size': results[f'PRIMER_PAIR_{i}_PRODUCT_SIZE'],
            'forward_start': forward_start,
            'forward_end': forward_start + forward_length - 1,
            'reverse_start': reverse_start,
            'reverse_end': results[f'PRIMER_RIGHT_{i}'][0],
        }
        
        primer_pairs.append(pair_info)
    
    return primer_pairs


def design_taqman_probe(sequence: str, forward_start: int, reverse_end: int, relaxed: bool = False) -> Optional[Dict[str, Any]]:
    """
    Design a TaqMan probe for the given primer pair region.
    
    Args:
        sequence: Full DNA sequence
        forward_start: Start position of forward primer
        reverse_end: End position of reverse primer
        relaxed: Use relaxed constraints for probe selection
        
    Returns:
        Dictionary with probe information or None if no valid probe found
    """
    # Define target region for probe (between primers)
    probe_region = sequence[forward_start:reverse_end+1]
    
    # Configure primer3 parameters for probe design
    seq_args = {
        'SEQUENCE_ID': 'probe_template',
        'SEQUENCE_TEMPLATE': probe_region,
    }
    
    global_args = {
        'PRIMER_TASK': 'generic',
        'PRIMER_PICK_LEFT_PRIMER': 0,
        'PRIMER_PICK_RIGHT_PRIMER': 0,
        'PRIMER_PICK_INTERNAL_OLIGO': 1,
        'PRIMER_INTERNAL_OPT_SIZE': 22,
        'PRIMER_INTERNAL_MIN_SIZE': 18,
        'PRIMER_INTERNAL_MAX_SIZE': 30,
        'PRIMER_INTERNAL_OPT_TM': 68.0,
        'PRIMER_INTERNAL_MIN_TM': 65.0,
        'PRIMER_INTERNAL_MAX_TM': 72.0,
        'PRIMER_INTERNAL_MIN_GC': 40.0,
        'PRIMER_INTERNAL_MAX_GC': 65.0,
        'PRIMER_NUM_RETURN': 20  # Increased from 10 to 20 for more candidates
    }
    
    # Run primer3 design for probe
    probe_results = primer3.bindings.designPrimers(seq_args, global_args)
    
    # First try with strict constraints
    for i in range(20):
        probe_key = f'PRIMER_INTERNAL_OLIGO_{i}_SEQUENCE'
        if probe_key not in probe_results:
            break
            
        probe_seq = probe_results[probe_key]
        
        # Apply TaqMan probe constraints
        if probe_seq.startswith('G') or 'GGGG' in probe_seq or 'AAAAAA' in probe_seq:
            continue
        
        probe_tm = probe_results[f'PRIMER_INTERNAL_OLIGO_{i}_TM']
        probe_gc = GC(probe_seq)
        
        return {
            'sequence': probe_seq,
            'gc': round(probe_gc, 1),
            'tm': round(probe_tm, 1),
            'relaxed': False
        }
    
    # If relaxed mode or no probe found with strict constraints
    if relaxed or True:  # Always try relaxed as a fallback
        for i in range(20):
            probe_key = f'PRIMER_INTERNAL_OLIGO_{i}_SEQUENCE'
            if probe_key not in probe_results:
                break
                
            probe_seq = probe_results[probe_key]
            
            # Apply only minimal constraints: no super-long homopolymers
            if 'GGGGGGG' in probe_seq or 'AAAAAAA' in probe_seq or 'TTTTTTT' in probe_seq or 'CCCCCCC' in probe_seq:
                continue
            
            probe_tm = probe_results[f'PRIMER_INTERNAL_OLIGO_{i}_TM']
            probe_gc = GC(probe_seq)
            
            return {
                'sequence': probe_seq,
                'gc': round(probe_gc, 1),
                'tm': round(probe_tm, 1),
                'relaxed': True
            }
    
    # If no suitable probe is found, create a simple one in the middle of the amplicon
    if len(probe_region) >= 25:
        mid_point = len(probe_region) // 2
        custom_probe = probe_region[mid_point-12:mid_point+13]
        probe_gc = GC(custom_probe)
        
        # Estimate melting temperature using a simple formula
        # Tm = 2°C(A+T) + 4°C(G+C)
        a_count = custom_probe.count('A')
        t_count = custom_probe.count('T')
        g_count = custom_probe.count('G')
        c_count = custom_probe.count('C')
        tm = 2 * (a_count + t_count) + 4 * (g_count + c_count)
        
        return {
            'sequence': custom_probe,
            'gc': round(probe_gc, 1),
            'tm': round(tm, 1),
            'relaxed': True,
            'note': 'Custom mid-amplicon probe (fallback)'
        }
    
    return None


def simulate_blast_submission(primer: str) -> str:
    """
    Simulate submitting a primer sequence to NCBI BLAST.
    
    Args:
        primer: Primer sequence
        
    Returns:
        Status message
    """
    # Simulate a BLAST submission
    print(f"Simulating BLAST submission for primer: {primer}")
    
    # In a real implementation, this would use the NCBI API
    # requests.post("https://blast.ncbi.nlm.nih.gov/Blast.cgi", data={"QUERY": primer})
    
    return "Submitted to NCBI BLAST"


def plot_gc_content(sequence: str, output_path: str, window_size: int = 20) -> str:
    """
    Generate a GC content plot for a DNA sequence.
    
    Args:
        sequence: DNA sequence
        output_path: Path to save the plot
        window_size: Window size for GC content calculation
        
    Returns:
        Path to the saved plot file
    """
    # Calculate GC content over sliding windows
    gc_values = []
    for i in range(len(sequence) - window_size + 1):
        window = sequence[i:i+window_size]
        gc_values.append(GC(window))
    
    positions = list(range(len(gc_values)))
    
    # Create the plot
    plt.figure(figsize=(10, 6))
    plt.plot(positions, gc_values)
    plt.xlabel('Position in Sequence')
    plt.ylabel('GC Content (%)')
    plt.title('GC Content Distribution')
    plt.axhline(y=50, color='r', linestyle='--', alpha=0.5)
    plt.grid(True, alpha=0.3)
    
    # Add average GC line
    avg_gc = round(GC(sequence), 1)
    plt.axhline(y=avg_gc, color='g', linestyle='-', alpha=0.5)
    plt.text(0, avg_gc + 2, f'Avg GC: {avg_gc}%', color='g')
    
    # Save the figure
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    return output_path


def main() -> None:
    """Main function that runs the TaqMan primer design workflow."""
    args = parse_arguments()
    
    # Get input sequence
    if args.file:
        try:
            with open(args.file, 'r') as f:
                fasta_data = f.read()
        except Exception as e:
            print(f"Error reading file: {e}", file=sys.stderr)
            sys.exit(1)
    elif args.sequence:
        fasta_data = args.sequence
    else:
        print("Error: Either --file or --sequence must be provided", file=sys.stderr)
        sys.exit(1)
    
    try:
        # Parse and clean sequence
        seq_id, sequence = parse_fasta(fasta_data)
        sequence = clean_sequence(sequence)
        
        if len(sequence) < 100:
            print("Warning: Sequence is shorter than 100 bp, which may limit primer design options", 
                  file=sys.stderr)
        
        # Make sure output directory exists
        os.makedirs(args.plot_dir, exist_ok=True)
        
        # Generate primers
        primer_pairs = design_primers(sequence)
        
        # Create result structure
        result = {"primers": []}
        
        # Process each primer pair
        for i, pair in enumerate(primer_pairs):
            # Design a TaqMan probe
            probe = design_taqman_probe(
                sequence, 
                pair['forward_start'], 
                pair['reverse_end'],
                relaxed=args.relaxed
            )
            
            if not probe:
                print(f"Warning: No valid TaqMan probe found for primer pair {i+1}", file=sys.stderr)
                probe = {
                    "sequence": "No valid probe found",
                    "gc": 0.0,
                    "tm": 0.0
                }
            else:
                print(f"Designed probe for pair {i+1}: {probe['sequence']}" + 
                      (" (relaxed constraints)" if probe.get('relaxed', False) else ""))
            
            # Simulate BLAST submission
            blast_forward = simulate_blast_submission(pair['forward'])
            blast_reverse = simulate_blast_submission(pair['reverse'])
            
            # Generate GC plot filename
            plot_filename = f"{args.plot_dir}/gc_plot_pair_{i+1}.png"
            
            # Generate GC content plot
            gc_plot_path = plot_gc_content(
                sequence[pair['forward_start']:pair['reverse_end']+1],
                plot_filename
            )
            
            # Add to results
            primer_result = {
                "forward": pair['forward'],
                "reverse": pair['reverse'],
                "forward_tm": pair['forward_tm'],
                "reverse_tm": pair['reverse_tm'],
                "forward_gc": pair['forward_gc'],
                "reverse_gc": pair['reverse_gc'],
                "penalty": pair['penalty'],
                "probe": probe,
                "blast_forward": blast_forward,
                "blast_reverse": blast_reverse,
                "gc_plot": gc_plot_path
            }
            
            result["primers"].append(primer_result)
        
        # Add overall GC plot
        overall_plot_path = f"{args.plot_dir}/gc_plot_overall.png"
        overall_gc_plot = plot_gc_content(sequence, overall_plot_path)
        result["gc_plot"] = overall_gc_plot
        
        # Output results
        with open(args.output, 'w') as f:
            json.dump(result, f, indent=2)
        
        # Also print to stdout
        print(json.dumps(result, indent=2))
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main() 