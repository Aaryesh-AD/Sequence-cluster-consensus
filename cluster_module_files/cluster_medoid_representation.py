# Compute the medoid
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
        if label == -1:
            continue

        cluster_indices = np.where(labels == label)[0]
        cluster_points = features[cluster_indices]
        

        medoid_index_within_cluster = compute_medoid(cluster_points)
        medoid_index = cluster_indices[medoid_index_within_cluster]
        
        medoids[label] = features[medoid_index]
        medoid_sequences[label] = sequences[medoid_index]
    
    return medoids, medoid_sequences

def save_medoid_sequences(medoid_sequences, output_file):
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as f:
        for cluster_id, medoid_sequence in medoid_sequences.items():
            f.write(f">{cluster_id}\n{medoid_sequence}\n")
    print(f"Medoid sequences saved to {output_file}")