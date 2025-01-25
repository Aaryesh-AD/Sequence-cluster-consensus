#!/usr/bin/env python3

"""
This script reads a FASTA file, merges the sequences, and writes them to a new FASTA file.

"""
import os
import warnings
import tempfile

warnings.filterwarnings("ignore", category=DeprecationWarning)

def read_fasta(file_path):
    with open(file_path, 'r') as file:
        fasta_data = file.read().split(">")[1:]
    sequences = {}
    for lines in fasta_data:
        header_key = f"Sample_{len(sequences)}"
        seq_value = lines.split("\n", 1)[1].replace("\n", "").strip()
        sequences[header_key] = seq_value
    return sequences

def write_fasta(sequences, output_file):
    try:
        with open(output_file, 'w') as file:
            for idx, seq in enumerate(sequences.values(), start=1):
                file.write(f">Sample_{idx}\n{seq}\n")
    except Exception as e:
        raise RuntimeError(f"Error writing to the output file: {e}")

def merge_fasta(input_fasta, output_fasta=None):
    try:
        sequences = read_fasta(input_fasta)
    except Exception as e:
        raise RuntimeError(f"Error reading the FASTA file: {e}")
    
    # If no output file is provided, create a temporary file
    if output_fasta is None:
        temp_file = tempfile.NamedTemporaryFile(delete=False, suffix=".fasta")
        output_fasta = temp_file.name
        temp_file.close()  # Close the temp file to allow writing later

    write_fasta(sequences, output_fasta)
    return output_fasta
