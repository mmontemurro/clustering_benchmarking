#!/usr/bin/env python
import os, sys
import argparse
import pandas as pd
from sklearn.cluster import AffinityPropagation
from sklearn.cluster import AgglomerativeClustering
from sklearn.cluster import DBSCAN
import hdbscan
from sklearn.base import ClusterMixin
from sklearn.mixture import GaussianMixture
from sklearn.cluster import KMeans
from kneed import KneeLocator
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors

# Wrapper function for Gaussian Mixture Modelling
class GMClusters(GaussianMixture, ClusterMixin):

    def __init__(self, n_clusters=2, **kwargs):
        kwargs["n_components"] = n_clusters
        super(GMClusters, self).__init__(**kwargs)

    def fit(self, X):
        super(GMClusters, self).fit(X)
        self.labels_ = self.predict(X)
        return self 

def APN(X, clusters, clust_method, **kwargs):
    # Average Proportion of Non-overlap (APN)
    #model = clust_method(**kwargs).fit(X)
    #labels = model.labels_
    
    return apn 

def AD(X, clusters, clust_method, **kwargs):
    # Average Distance (AD)
    return ad 

def ADM(X, clusters, clust_method, **kwargs):
    # Average Distance between Means (ADM)
    return adm

def FOM(X, clusters, clust_method, **kwargs):
    # Figure of Merit (FOM)
    return fom 

def main():
    parser = argparse.ArgumentParser(description="Cluster Stability Estimator.")

    parser.add_argument("-c", "--clusters_file", help="Clusters file", required=True)
    parser.add_argument("-f", "--file", help="Raw data file.", required=True, type=str)
    parser.add_argument("-s", "--stats_file", help="Stats file", required=True)
    parser.add_argument("-m", "--method", 
            choices=["affinity", "agglomerative_average", "agglomerative_complete", "agglomerative_single", 
            "agglomerative_ward", "birch", "dbscan", "gaussian", "hdbscan", "kmeans"], 
            help="Clustering method", required=True)
    parser.add_argument('-p', '--preproc', action='store_true', help="To be specified if the input file requires preprocessing")
    parser.add_argument("o", "--out_prefix", help="Output prefix", required=True)

    args = parser.parse_args()

    clusters = pd.read_csv(args.clusters_file, sep="\t")
    clusters.columns = ['cell', 'cluster'] #rename columns, otherwise the cell column is named 'Unnamed: 0'

    if preproc:
        matrix = pd.read_csv(args.file, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END']).transpose()
    else:
        matrix = pd.read_csv(args.file, index_col=0, sep="\t")

    stats = pd.read_csv(filepath, index_col=0, sep="\t", header=None).transpose()
    #to fix an error
    stats = stats.loc[:,~stats.columns.duplicated()]

    outprefix = args.out_prefix

    method = args.method

    k = stats["n_clusters_unfiltered"].values[0]

    if "DBSCAN_eps" in stats.columns.values:
        eps = stats["DBSCAN_eps"]
    if "HDBSCAN_min_cluster_size" in stats.columns.values:
        min_cluster_size = stats["HDBSCAN_min_cluster_size"]
    
    if method == "affinity":
        apn = APN(matrix, clusters, AffinityPropagation)
        ad = AD(matrix, clusters, AffinityPropagation)
        adm = ADM(matrix, clusters, AffinityPropagation)
        fom = FOM(matrix, clusters, AffinityPropagation)
    elif method == "agglomerative_average":
        apn = APN(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l1", linkage="average")
        ad = AD(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l1", linkage="average")
        adm = ADM(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l1", linkage="average")
        fom = FOM(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l1", linkage="average")
    elif method == "agglomerative_complete":
        apn = APN(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l1", linkage="complete")
        ad = AD(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l1", linkage="complete")
        adm = ADM(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l1", linkage="complete")
        fom = FOM(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l1", linkage="complete")
    elif method == "agglomerative_single":
        apn = APN(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l1", linkage="single")
        ad = AD(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l1", linkage="single")
        adm = ADM(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l1", linkage="single")
        fom = FOM(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l1", linkage="single")
    elif method == "agglomerative_ward":
        apn = APN(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l2", linkage="ward")
        ad = AD(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l2", linkage="ward")
        adm = ADM(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l2", linkage="ward")
        fom = FOM(matrix, clusters, AgglomerativeClustering, n_clusters=k, affinity="l2", linkage="ward")
    elif method == "birch":
        apn = APN(matrix, clusters, Birch, branching_factor=20, n_clusters=int(k), threshold=0.5, compute_labels=True)
        ad = AD(matrix, clusters, Birch, branching_factor=20, n_clusters=int(k), threshold=0.5, compute_labels=True)
        adm = ADM(matrix, clusters, Birch, branching_factor=20, n_clusters=int(k), threshold=0.5, compute_labels=True)
        fom = FOM(matrix, clusters, Birch, branching_factor=20, n_clusters=int(k), threshold=0.5, compute_labels=True)
    elif method == "dbscan":
        apn = APN(matrix, clusters, DBSCAN, eps=eps)
        ad = AD(matrix, clusters, DBSCAN, eps=eps)
        adm = ADM(matrix, clusters, DBSCAN, eps=eps)
        fom = FOM(matrix, clusters, DBSCAN, eps=eps)
    elif method == "hdbscan":
        apn = APN(matrix, clusters, hdbscan.HDBSCAN, min_cluster_size=min_cluster_size, min_samples=1, cluster_selection_epsilon=0.5)
        ad = AD(matrix, clusters, hdbscan.HDBSCAN, min_cluster_size=min_cluster_size, min_samples=1, cluster_selection_epsilon=0.5)
        adm = ADM(matrix, clusters, hdbscan.HDBSCAN, min_cluster_size=min_cluster_size, min_samples=1, cluster_selection_epsilon=0.5)
        fom = FOM(matrix, clusters, hdbscan.HDBSCAN, min_cluster_size=min_cluster_size, min_samples=1, cluster_selection_epsilon=0.5)
    elif method == "gaussian":
        apn = APN(matrix, clusters, GMClusters, n_components=k)
        ad = AD(matrix, clusters, GMClusters, n_components=k)
        adm = ADM(matrix, clusters, GMClusters, n_components=k)
        fom = FOM(matrix, clusters, GMClusters, n_components=k)
    elif method == "kmeans":
        apn = APN(matrix, clusters, KMeans, n_clusters=k)
        ad = AD(matrix, clusters, KMeans, n_clusters=k)
        adm = ADM(matrix, clusters, KMeans, n_clusters=k)
        fom = FOM(matrix, clusters, KMeans, n_clusters=k)

    with open(args.stats, 'a') as f:
        f.write("\n")
        f.write("APN\t" + str(apn) + "\n")
        f.write("AD\t" + str(ad)+ "\n")
        f.write("ADM\t" + str(adm)+ "\n")
        f.write("FOM\t" + str(fom)+ "\n")

if __name__ == "__main__":
    sys.exit(main())