#!/usr/bin/env python3
"""
Test script for the TaqMan PCR Primer and Probe Design Tool.
"""

import os
import subprocess
import json
import sys


def run_test():
    """Run a basic test of the primer_express.py script."""
    print("Testing TaqMan PCR Primer and Probe Design Tool...")
    
    # Create test output directory
    test_dir = "test_output"
    os.makedirs(test_dir, exist_ok=True)
    
    # Run the tool with example.fasta
    cmd = [
        sys.executable, 
        "primer_express.py", 
        "-f", "example.fasta",
        "-o", f"{test_dir}/test_results.json",
        "-p", f"{test_dir}/plots"
    ]
    
    print(f"Running command: {' '.join(cmd)}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    
    if result.returncode != 0:
        print("Test failed! Error output:")
        print(result.stderr)
        return False
    
    # Check if output file exists
    output_file = f"{test_dir}/test_results.json"
    if not os.path.exists(output_file):
        print(f"Test failed! Output file {output_file} not found.")
        return False
    
    # Check output JSON structure
    try:
        with open(output_file, 'r') as f:
            data = json.load(f)
        
        # Verify structure
        if "primers" not in data:
            print("Test failed! Output JSON missing 'primers' key.")
            return False
        
        if "gc_plot" not in data:
            print("Test failed! Output JSON missing 'gc_plot' key.")
            return False
        
        # Check if we have primers
        if not data["primers"]:
            print("Test failed! No primers were generated.")
            return False
        
        # Check a primer entry
        primer = data["primers"][0]
        required_keys = ["forward", "reverse", "forward_tm", "reverse_tm", 
                         "forward_gc", "reverse_gc", "probe", "blast_forward", 
                         "blast_reverse", "gc_plot"]
        
        for key in required_keys:
            if key not in primer:
                print(f"Test failed! Primer data missing '{key}' key.")
                return False
        
        # Check if GC plot file exists
        gc_plot_file = data["gc_plot"]
        if not os.path.exists(gc_plot_file):
            print(f"Test failed! GC plot file {gc_plot_file} not found.")
            return False
        
        print("Test passed! Output looks good.")
        return True
        
    except Exception as e:
        print(f"Test failed! Exception: {e}")
        return False


if __name__ == "__main__":
    success = run_test()
    if not success:
        sys.exit(1) 