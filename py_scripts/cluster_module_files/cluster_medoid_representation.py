#! usr/bin/env python3

'''
This script processes clustering results to identify and save representative sequences based on finding the medoids for each cluster.
'''

import os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool, cpu_count
import warnings


def calculate_medoid(distance_matrix):
    """Find the medoid index for the cluster."""
    total_distances = np.sum(distance_matrix, axis=1)
    medoid_index = np.argmin(total_distances)
    return medoid_index


def hamming_distance(seq1, seq2):
    """Calculate Hamming distance between two sequences."""
    return sum(a != b for a, b in zip(seq1, seq2))


def generate_distance_matrix(sequences):
    """Generate a symmetric Hamming distance matrix."""
    num_sequences = len(sequences)
    distance_matrix = np.zeros((num_sequences, num_sequences), dtype=int)

    # Compute only the upper triangle
    for i in range(num_sequences):
        for j in range(i + 1, num_sequences):
            dist = hamming_distance(sequences[i], sequences[j])
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist  # Mirror the upper triangle
    return distance_matrix


def process_cluster(args):
    """Process a single cluster and find its medoid."""
    cluster_path, cluster_name = args
    cluster_sequences = [str(record.seq) for record in SeqIO.parse(cluster_path, "fasta")]

    if len(cluster_sequences) < 2:
        print(f"Skipping {cluster_name}: Less than two sequences.")
        return None

    # Calculate distance matrix and find medoid
    distance_matrix = generate_distance_matrix(cluster_sequences)
    medoid_index = calculate_medoid(distance_matrix)
    medoid_sequence = cluster_sequences[medoid_index]

    # Return medoid record
    return SeqRecord(
        Seq(medoid_sequence),
        id=f"{cluster_name}_medoid",
        description=f"Medoid sequence for {cluster_name}",
    )


def process_medoid_for_methods(cluster_dir, output_dir, sequences):
    os.makedirs(output_dir, exist_ok=True)
    clustering_methods = [d for d in os.listdir(cluster_dir) if os.path.isdir(os.path.join(cluster_dir, d))]

    for method in clustering_methods:
        method_dir = os.path.join(cluster_dir, method)
        medoid_output_file = os.path.join(output_dir, f"{method}_medoid_sequences.fasta")

        print(f"Processing medoids for method: {method}")
        cluster_files = [
            (os.path.join(method_dir, cluster_file), os.path.splitext(cluster_file)[0])
            for cluster_file in os.listdir(method_dir)
            if cluster_file.endswith(".fasta") and cluster_file != "cluster_-1.fasta"
        ]

        # Parallelize cluster processing
        with Pool(cpu_count()) as pool:
            medoid_records = pool.map(process_cluster, cluster_files)

        # Write medoid sequences to file
        medoid_records = [record for record in medoid_records if record is not None]
        with open(medoid_output_file, "w") as output_handle:
            SeqIO.write(medoid_records, output_handle, "fasta")

        print(f"Medoid sequences saved to {medoid_output_file}")