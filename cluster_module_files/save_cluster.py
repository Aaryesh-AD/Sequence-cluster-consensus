import os

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
            cluster_file = os.path.join(method_dir, f"cluster_{label}.fasta")
            with open(cluster_file, 'w') as f:
                for seq, header in cluster_seqs:
                    f.write(f"{header}\n{seq}\n")

        print(f"Saved {len(clusters)} clusters for {method_name} in {method_dir}")
