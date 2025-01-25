#!/usr/bin/env python3

"""
This script performs clustering and generates consensus sequences from multiple sequence alignments using Biopython.

"""

import os
import glob
from Bio.Align import MultipleSeqAlignment
from Bio.Align import AlignInfo
from Bio import AlignIO
from Bio import BiopythonWarning

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)


def perform_msa(fasta_files, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for fasta_file in fasta_files:
        cluster_name = os.path.splitext(os.path.basename(fasta_file))[0]
        output_file = os.path.join(output_dir, f"{cluster_name}_aligned.aln")

        # Read the sequences
        records = list(AlignIO.read(fasta_file, "fasta"))
        if len(records) < 2:
            print(f"Skipping {fasta_file}: Less than two sequences available.")
            continue

        # Create a dummy alignment (identity alignment)
        alignment = MultipleSeqAlignment(records)

        # Save the alignment to a file
        AlignIO.write(alignment, output_file, "clustal")
        print(f"Alignment completed for {fasta_file}. Output saved to {output_file}")


def consensus_sequences(msa_dir, output_file, threshold=0.7, ambiguous_char="N"):
    msa_files = glob.glob(os.path.join(msa_dir, "*.aln"))
    if not msa_files:
        print(f"No MSA files found in {msa_dir}.")
        return

    os.makedirs(os.path.dirname(output_file), exist_ok=True)

    with open(output_file, "w") as output_handle:
        for msa_file in msa_files:
            try:
                alignment = AlignIO.read(msa_file, "clustal")
                if not alignment:
                    print(f"Skipping empty or invalid file: {msa_file}")
                    continue

                # Generate consensus sequence
                summary_info = AlignInfo.SummaryInfo(alignment)
                consensus = summary_info.dumb_consensus(threshold=threshold, ambiguous=ambiguous_char)

                cluster_name = os.path.splitext(os.path.basename(msa_file))[0]
                output_handle.write(f">{cluster_name}_consensus\n{consensus}\n")
                print(f"Consensus sequence generated for {msa_file} and saved.")
            except Exception as e:
                print(f"Error processing {msa_file}: {e}")

    print(f"All consensus sequences saved to {output_file}.")



