#! /usr/bin/env python3

'''
This script saves clustered sequences into separate FASTA files for each clustering method.

'''
import os
import shutil
from tempfile import NamedTemporaryFile

def save_cluster_sequences(sequences, labels_dict, output_dir, original_headers):
    for method_name, labels in labels_dict.items():
        # Create method-specific directory
        method_dir = os.path.join(output_dir, method_name)
        os.makedirs(method_dir, exist_ok=True)

        # Group sequences by cluster
        clusters = {}
        for seq, label, header in zip(sequences, labels, original_headers):
            if label not in clusters:
                clusters[label] = []
            clusters[label].append((seq, header))

        # Save sequences for each cluster
        for label, cluster_seqs in clusters.items():
            try:
                with NamedTemporaryFile(delete=False, mode="w", dir="/tmp") as temp_file:
                    for seq, header in cluster_seqs:
                        temp_file.write(f"{header}\n{seq}\n")
                    temp_file_path = temp_file.name
                
                # Move the temp file to the target directory
                cluster_file = os.path.join(method_dir, f"cluster_{label}.fasta")
                shutil.move(temp_file_path, cluster_file)

            except Exception as e:
                print(f"Error writing cluster {label} for {method_name}: {e}")


