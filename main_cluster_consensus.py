import os
import argparse
from cluster_module_files.fasta_processor import merge_fasta
from cluster_module_files.feature_extractor import load_sequences
from cluster_module_files.sequence_clustering import cluster_sequences, generate_labels_dict
from cluster_module_files.save_cluster import save_cluster_sequences
from cluster_module_files.plot_clusters import visualize_all_clusters
from cluster_module_files.cluster_consensus import clustal_cluster, consensus_sequences
from cluster_module_files.cluster_medoid_representation import process_medoid_for_methods


def argument():
    """
    Parses command-line arguments for clustering, medoid, and consensus sequence generation.
    """
    parser = argparse.ArgumentParser(description="Cluster sequences and generate representative sequences.")
    parser.add_argument("-i", "--input_fasta", type=str, required=True, help="Input FASTA file containing amino acid sequences.")
    parser.add_argument("-o", "--output_dir", type=str, required=True, help="Output directory to save results.")
    parser.add_argument("--consensus_threshold", type=float, default=0.7, help="Threshold for consensus sequence generation.")
    parser.add_argument("--ambiguous_char", type=str, default="N", help="Ambiguous character for consensus sequence generation.")

    return parser.parse_args()


def main():
    """
    Main function to orchestrate clustering, medoid, and consensus sequence generation.
    """
    args = argument()
    os.makedirs(args.output_dir, exist_ok=True)

    # Step 1: Merge sequences
    merged_fasta = merge_fasta(args.input_fasta, os.path.join(args.output_dir, "merged_sequences.fasta"))
    sequences = load_sequences(merged_fasta)

    # Step 2: Perform clustering
    best_labels, reduced_features = cluster_sequences(sequences)
    labels_dict = generate_labels_dict(reduced_features, sequences, n_clusters=10, eps=0.5, min_samples=5)

    # Step 3: Visualize and save cluster plots
    visualize_all_clusters(reduced_features, labels_dict, os.path.join(args.output_dir, "plots"))

    # Step 4: Save clustered sequences
    original_headers = [f">Sample_{i}" for i in range(len(sequences))]
    save_cluster_sequences(sequences, labels_dict, os.path.join(args.output_dir, "clusters"), original_headers)

    # Step 5: Generate consensus sequences for all clustering methods
    consensus_output_dir = os.path.join(args.output_dir, "consensus_sequences")
    os.makedirs(consensus_output_dir, exist_ok=True)

    clustering_methods = ["K-means", "DBSCAN", "Hierarchical", "Birch", "Agglomerative", "OPTICS", "HDBSCAN"]

    for method in clustering_methods:
        # Ensure cluster subfolder exists
        cluster_method_dir = os.path.join(args.output_dir, "clusters", method)
        if not os.path.exists(cluster_method_dir):
            print(f"Skipping {method}: No cluster folder found.")
            continue

        # Create method-specific alignment and consensus output paths
        method_alignment_dir = os.path.join(consensus_output_dir, f"{method}_alignments")
        method_consensus_file = os.path.join(consensus_output_dir, f"consensus_cluster_sequence_{method}.fasta")

        print(f"\nProcessing consensus for clustering method: {method}")
        # Align sequences within each cluster and generate consensus
        clustal_cluster(cluster_method_dir, method_alignment_dir)
        consensus_sequences(
            method_alignment_dir,
            method_consensus_file,
            threshold=args.consensus_threshold,
            ambiguous_char=args.ambiguous_char
        )

    # Step 6: Generate representative medoid sequences for all clustering methods
    representative_output_dir = os.path.join(args.output_dir, "representative_sequences")
    process_medoid_for_methods(os.path.join(args.output_dir, "clusters"), representative_output_dir, clustering_methods)

    print("All clustering, consensus, and medoid sequences processed and saved.")


if __name__ == "__main__":
    main()
