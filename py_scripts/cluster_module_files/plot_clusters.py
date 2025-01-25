#! /usr/bin/env python3

'''
This script visualizes clustering results using PCA and dendrograms.

'''

import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
import os

import warnings
from sklearn.utils._testing import ignore_warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

def plot_clusters(features, labels, title, output_dir):
    pca = PCA(n_components=2)
    pca_features = pca.fit_transform(features)

    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    colors = plt.cm.tab10(np.linspace(0, 1, n_clusters))

    plt.figure(figsize=(10, 8))
    for label, color in zip(unique_labels, colors):
        if label == -1:  # Noise in DBSCAN or HDBSCAN
            plt.scatter(
                pca_features[labels == label, 0],
                pca_features[labels == label, 1],
                c="gray",
                marker="x",
                label="Noise"
            )
        else:
            plt.scatter(
                pca_features[labels == label, 0],
                pca_features[labels == label, 1],
                c=[color],
                label=f"Cluster {label}"
            )

    plt.title(title)
    plt.xlabel("PCA Component 1")
    plt.ylabel("PCA Component 2")
    plt.legend(loc="best", bbox_to_anchor=(1.05, 1), title="Clusters")
    plt.grid(alpha=0.3)
    plt.tight_layout()

    # Save the plot
    plot_path = os.path.join(output_dir, f"{title.replace(' ', '_')}.png")
    plt.savefig(plot_path)
    plt.close()



def visualize_all_clusters(features, labels_dict, output_dir):
    os.makedirs(output_dir, exist_ok=True)

    for method, labels in labels_dict.items():
        # Create PCA-based cluster plot
        plot_clusters(features, labels, f"{method} Clustering", output_dir)

