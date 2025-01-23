import matplotlib.pyplot as plt
import numpy as np
from sklearn.decomposition import PCA
from scipy.cluster.hierarchy import dendrogram

def plot_clusters(features, labels, title):
    
    pca = PCA(n_components=2)
    pca_features = pca.fit_transform(features)

    unique_labels = np.unique(labels)
    n_clusters = len(unique_labels)
    colors = plt.cm.tab10(np.linspace(0, 1, n_clusters))

    plt.figure(figsize=(10, 8))
    for label, color in zip(unique_labels, colors):
        if label == -1: 
            plt.scatter(
                pca_features[labels == label, 0],
                pca_features[labels == label, 1],
                c='gray',
                marker='x',
                label='Noise'
            )
        else:
            plt.scatter(
                pca_features[labels == label, 0],
                pca_features[labels == label, 1],
                c=[color],
                label=f'Cluster {label}'
            )

    plt.title(title)
    plt.xlabel('PCA Component 1')
    plt.ylabel('PCA Component 2')
    plt.legend(loc='best', bbox_to_anchor=(1.05, 1), title="Clusters")
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.show()

def plot_dendrogram(linkage_matrix, title):
    plt.figure(figsize=(10, 8))
    dendrogram(linkage_matrix, truncate_mode='lastp', p=30, show_leaf_counts=True, leaf_rotation=45, leaf_font_size=12)
    plt.title(title)
    plt.xlabel('Cluster Size')
    plt.ylabel('Distance')
    plt.show()

def visualize_all_clusters(features, labels_dict):
    for method, labels in labels_dict.items():
        plot_clusters(features, labels, f"{method} Clustering")

