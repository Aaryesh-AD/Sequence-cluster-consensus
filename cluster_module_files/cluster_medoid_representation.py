#! usr/bin/env python3

'''
This script processes clustering results to identify and save representative sequences based on finding the medoids for each cluster.

'''

import os
import numpy as np
from sklearn.metrics import pairwise_distances


def compute_medoid(cluster_data):
    distances = pairwise_distances(cluster_data)
    medoid_index = np.argmin(distances.sum(axis=0))
    return medoid_index


def medoid_seq_with_data(features, labels, sequences):
    unique_labels = np.unique(labels)
    medoids = {}
    medoid_sequences = {}

    for label in unique_labels:
        if label == -1:  # Skip noise points
            continue

        cluster_indices = np.where(labels == label)[0]
        cluster_points = features[cluster_indices]

        # Compute medoid
        medoid_index_within_cluster = compute_medoid(cluster_points)
        medoid_index = cluster_indices[medoid_index_within_cluster]

        medoids[label] = features[medoid_index]
        medoid_sequences[label] = sequences[medoid_index]

    return medoids, medoid_sequences


def save_medoid_sequences(medoid_sequences, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, "w") as f:
        for cluster_id, medoid_sequence in medoid_sequences.items():
            f.write(f">{cluster_id}\n{medoid_sequence}\n")
    print(f"Medoid sequences saved to {output_file}")


def process_medoid_for_methods(root_dir, output_dir, clustering_methods):
    os.makedirs(output_dir, exist_ok=True)

    for method in clustering_methods:
        method_dir = os.path.join(root_dir, method)
        if not os.path.exists(method_dir):
            print(f"Skipping {method}: Directory not found.")
            continue
        print(f"Processing medoid sequences for {method}...")

        try:
            features = np.load(os.path.join(method_dir, "features.npy"))
            labels = np.load(os.path.join(method_dir, "labels.npy"))
            sequences = np.load(os.path.join(method_dir, "sequences.npy"), allow_pickle=True)
        except Exception as e:
            print(f"Error loading data for {method}: {e}")
            continue
        medoids, medoid_sequences = medoid_seq_with_data(features, labels, sequences)

        # Save medoid sequences
        output_file = os.path.join(output_dir, f"representative_cluster_sequence_{method}.fasta")
        save_medoid_sequences(medoid_sequences, output_file)

    print("All representative medoid sequences generated and saved.")
