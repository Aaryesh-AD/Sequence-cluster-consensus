'''
Argparse: Fasta input, output directory setup (merged_fasta, cluster files, plots, representative sequences, consensus sequences)

parameters: {K-mers, [lambda, w (from PAAC)] from feature_extractor.py},

n clusters: {K-means, hierarchical, Agglomerative, Birch}, top 10 clusters by size, silhouette score,
eps for DBSCAN, xi, min_samples for OPTICS and HDBSCAN, min_cluster_size and min_samples for HDBSCAN

plot figure size

consensus sequences: threshold, ambiguous character
'''