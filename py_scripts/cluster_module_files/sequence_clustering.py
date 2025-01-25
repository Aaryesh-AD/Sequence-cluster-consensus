#! /usr/bin/env python3

'''
This script performs clustering on sequences using various clustering algorithms and evaluates their performance.

'''

from .feature_extractor import combine_features
import numpy as np
import hdbscan
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans, DBSCAN, Birch, AgglomerativeClustering, OPTICS
from sklearn.metrics import silhouette_score
from scipy.cluster.hierarchy import linkage, fcluster

import warnings
from sklearn.utils._testing import ignore_warnings
warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=DeprecationWarning)

def kmeans_clustering(features, n_clusters):
    kmeans = KMeans(n_clusters=n_clusters, random_state=69)
    labels = kmeans.fit_predict(features)
    return labels

def dbscan_clustering(features, eps, min_samples):
    dbscan = DBSCAN(eps=eps, min_samples=min_samples)
    labels = dbscan.fit_predict(features)
    return labels

def hierarchical_clustering(features, n_clusters):
    linkage_matrix = linkage(features, method='ward')
    labels = fcluster(linkage_matrix, n_clusters, criterion='maxclust')
    return labels

def birch_clustering(features, n_clusters):
    birch = Birch(n_clusters=n_clusters)
    labels = birch.fit_predict(features)
    return labels

def agglomerative_clustering(features, n_clusters):
    agglomerative = AgglomerativeClustering(n_clusters=n_clusters)
    labels = agglomerative.fit_predict(features)
    return labels

def optics_clustering(features, min_samples, xi, min_cluster_size):
    optics = OPTICS(min_samples=min_samples, xi=xi, min_cluster_size=min_cluster_size)
    labels = optics.fit_predict(features)
    return labels

def hdbscan_clustering(features, min_cluster_size, min_samples):
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=min_samples)
    labels = clusterer.fit_predict(features)
    
    unique_labels, counts = np.unique(labels[labels != -1], return_counts=True)
    
    sorted_clusters = sorted(zip(unique_labels, counts), key=lambda x: x[1], reverse=True)
    
    top_clusters = sorted_clusters[:min(10, len(sorted_clusters))]
    top_cluster_labels = [label for label, _ in top_clusters]
    
    new_labels = np.full_like(labels, -1)
    for new_label, (old_label, _) in enumerate(top_clusters):
        new_labels[labels == old_label] = new_label
    
    return new_labels

# Evaluation function
def evaluate_clustering(features, labels):
    if len(set(labels)) > 1 and -1 not in set(labels):
        return silhouette_score(features, labels)
    else:
        return -1 

def cluster_sequences(sequences, k=3, lambda_=30, w=0.05, n_clusters=10, eps = 0.05, min_samples=5):
    features = combine_features(sequences, k, lambda_, w)
    scaler = StandardScaler()
    normalized_features = scaler.fit_transform(features)
    pca = PCA(n_components=min(50, len(normalized_features[0])))
    reduced_features = pca.fit_transform(normalized_features)
    
    kmeans_labels = kmeans_clustering(reduced_features, n_clusters)
    dbscan_labels = dbscan_clustering(reduced_features, eps, min_samples)
    hierarchical_labels = hierarchical_clustering(reduced_features, n_clusters)
    birch_labels = birch_clustering(reduced_features, n_clusters)
    agglomerative_labels = agglomerative_clustering(reduced_features, n_clusters)
    optics_labels = optics_clustering(reduced_features, min_samples, xi=eps, min_cluster_size=int(eps * len(sequences)))
    hdbscan_labels = hdbscan_clustering(reduced_features, min_cluster_size=5, min_samples = 5)
    

    kmeans_score = evaluate_clustering(reduced_features, kmeans_labels)
    dbscan_score = evaluate_clustering(reduced_features, dbscan_labels)
    hierarchical_score = evaluate_clustering(reduced_features, hierarchical_labels)
    birch_score = evaluate_clustering(reduced_features, birch_labels)
    agglomerative_score = evaluate_clustering(reduced_features, agglomerative_labels)
    optics_score = evaluate_clustering(reduced_features, optics_labels)
    hdbscan_score = evaluate_clustering(reduced_features, hdbscan_labels)
    
    print(f"K-means Silhouette Score: {kmeans_score:.3f}")
    print(f"DBSCAN Silhouette Score: {dbscan_score:.3f}")
    print(f"Hierarchical Silhouette Score: {hierarchical_score:.3f}")
    print(f"Birch Silhouette Score: {birch_score:.3f}")
    print(f"Agglomerative Silhouette Score: {agglomerative_score:.3f}")
    print(f"OPTICS Silhouette Score: {optics_score:.3f}")
    print(f"HDBSCAN Silhouette Score: {hdbscan_score:.3f}")
    
    all_labels = [kmeans_labels, dbscan_labels, hierarchical_labels, birch_labels, agglomerative_labels, optics_labels, hdbscan_labels]
    best_labels = max(all_labels, key=lambda x: evaluate_clustering(reduced_features, x))
    return best_labels, reduced_features

def generate_labels_dict(reduced_features, sequences, n_clusters=10, eps=0.5, min_samples=5):
    return {
        "K-means": KMeans(n_clusters=n_clusters, random_state=42).fit_predict(reduced_features),
        "DBSCAN": DBSCAN(eps=eps, min_samples=min_samples).fit_predict(reduced_features),
        "Hierarchical": fcluster(linkage(reduced_features, method='ward'), n_clusters, criterion='maxclust'),
        "Birch": Birch(n_clusters=n_clusters).fit_predict(reduced_features),
        "Agglomerative": AgglomerativeClustering(n_clusters=n_clusters).fit_predict(reduced_features),
        "OPTICS": OPTICS(min_samples=min_samples, xi=0.05, min_cluster_size=int(0.05 * len(sequences))).fit_predict(reduced_features),
        "HDBSCAN": hdbscan.HDBSCAN(min_cluster_size=min_samples, min_samples=min_samples).fit_predict(reduced_features),
    }
    