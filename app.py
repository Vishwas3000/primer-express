#!/usr/bin/env python3
"""
Web interface for TaqMan PCR Primer and Probe Design Tool
"""

import os
import json
import tempfile
import uuid
from flask import Flask, render_template, request, jsonify, send_file, url_for
from io import StringIO

# Import functions from primer_express
from primer_express import (
    parse_fasta,
    clean_sequence,
    design_primers,
    design_taqman_probe,
    simulate_blast_submission,
    plot_gc_content
)

app = Flask(__name__)
app.config['UPLOAD_FOLDER'] = 'static/uploads'
app.config['PLOT_DIR'] = 'static/plots'

# Ensure upload and plot directories exist
os.makedirs(app.config['UPLOAD_FOLDER'], exist_ok=True)
os.makedirs(app.config['PLOT_DIR'], exist_ok=True)


@app.route('/')
def index():
    """Render the main page"""
    return render_template('index.html')


@app.route('/design', methods=['POST'])
def design():
    """Handle primer design requests"""
    try:
        # Get FASTA data from form or file
        fasta_data = request.form.get('sequence', '')
        fasta_file = request.files.get('file')
        
        if fasta_file and fasta_file.filename:
            fasta_data = fasta_file.read().decode('utf-8')
        
        if not fasta_data:
            return jsonify({'error': 'No sequence provided'}), 400
        
        # Parse and clean sequence
        seq_id, sequence = parse_fasta(fasta_data)
        sequence = clean_sequence(sequence)
        
        if len(sequence) < 100:
            return jsonify({'error': 'Sequence is too short (< 100 bp)'}), 400
        
        # Generate a unique ID for this design job
        job_id = str(uuid.uuid4())
        plot_dir = os.path.join(app.config['PLOT_DIR'], job_id)
        os.makedirs(plot_dir, exist_ok=True)
        
        # Generate primers
        primer_pairs = design_primers(sequence)
        
        # Create result structure
        result = {"primers": [], "job_id": job_id}
        
        # Process each primer pair
        for i, pair in enumerate(primer_pairs):
            # Design a TaqMan probe
            probe = design_taqman_probe(
                sequence, 
                pair['forward_start'], 
                pair['reverse_end']
            )
            
            if not probe:
                probe = {
                    "sequence": "No valid probe found",
                    "gc": 0.0,
                    "tm": 0.0
                }
            
            # Simulate BLAST submission
            blast_forward = simulate_blast_submission(pair['forward'])
            blast_reverse = simulate_blast_submission(pair['reverse'])
            
            # Generate GC plot filename
            plot_filename = f"{plot_dir}/gc_plot_pair_{i+1}.png"
            rel_path = os.path.relpath(plot_filename, 'static')
            
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
                "gc_plot": rel_path
            }
            
            result["primers"].append(primer_result)
        
        # Add overall GC plot
        overall_plot_path = f"{plot_dir}/gc_plot_overall.png"
        overall_gc_plot = plot_gc_content(sequence, overall_plot_path)
        result["gc_plot"] = os.path.relpath(overall_gc_plot, 'static')
        
        # Save results to file
        results_file = os.path.join(plot_dir, "results.json")
        with open(results_file, 'w') as f:
            json.dump(result, f, indent=2)
        
        return jsonify(result)
        
    except Exception as e:
        return jsonify({'error': str(e)}), 500


@app.route('/download/<job_id>')
def download_results(job_id):
    """Download results as JSON file"""
    results_file = os.path.join(app.config['PLOT_DIR'], job_id, "results.json")
    if not os.path.exists(results_file):
        return jsonify({'error': 'Results not found'}), 404
    
    return send_file(results_file, as_attachment=True, download_name='primer_results.json')


if __name__ == '__main__':
    app.run(debug=True) 