#! usr/bin/env python3

'''
This script performs clustering using Clustal Omega and generates consensus sequences from the clustered alignments.

'''

import os
import glob
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio import AlignIO
from Bio.Align.AlignInfo import SummaryInfo

def clustal_cluster(input_dir, output_dir, clustalo_path="clustalo"):
    os.makedirs(output_dir, exist_ok=True)

    fasta_files = glob.glob(os.path.join(input_dir, "*.fasta"))
    if not fasta_files:
        print(f"No FASTA files found in {input_dir}.")
        return

    for fasta_file in fasta_files:
        cluster_name = os.path.splitext(os.path.basename(fasta_file))[0]
        output_file = os.path.join(output_dir, f"{cluster_name}_aligned.aln")

        clustalomega_cline = ClustalOmegaCommandline(
            cmd=clustalo_path,
            infile=fasta_file,
            outfile=output_file,
            outfmt="clustal",
            verbose=True,
            auto=True
        )

        print(f"Running Clustal Omega for {fasta_file}...")
        try:
            stdout, stderr = clustalomega_cline()
            print(f"Alignment completed for {fasta_file}. Output saved to {output_file}")
        except Exception as e:
            print(f"Error processing {fasta_file}: {e}")


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

                summary_info = SummaryInfo(alignment)
                consensus = summary_info.dumb_consensus(threshold=threshold, ambiguous=ambiguous_char)

                cluster_name = os.path.splitext(os.path.basename(msa_file))[0]
                output_handle.write(f">{cluster_name}_consensus\n{consensus}\n")
                print(f"Consensus sequence generated for {msa_file} and saved.")
            except Exception as e:
                print(f"Error processing {msa_file}: {e}")

    print(f"Consensus sequences saved to {output_file}.")
