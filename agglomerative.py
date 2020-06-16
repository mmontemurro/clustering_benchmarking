#!/usr/bin/envs python

import time
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from yellowbrick.cluster import KElbowVisualizer

import argparse, os, sys

def locate_elbow(matrix, outprefix, linkage):
    if linkage == "ward":
        affinity="euclidean"
    else:
        affinity="l1"
    model = AgglomerativeClustering(affinity=affinity, linkage=linkage)
    visualizer = KElbowVisualizer(model, k=(2,50))
    visualizer.fit(matrix.values)
    visualizer.finalize()
    plt.savefig(outprefix+"_agglomerative_" + linkage + "_elbow.png")
    plt.clf()
    return visualizer.elbow_value_


def cluster_and_output(k, matrix, outprefix, linkage):
    start_time = time.time()
    if linkage == "ward":
        affinity="euclidean"
    else:
        affinity="l1"
    model = AgglomerativeClustering(n_clusters=k, affinity=affinity, linkage=linkage).fit(matrix.values)
    end_time = time.time()
    clusters = pd.DataFrame(model.labels_, index=matrix.index, columns=["cluster"]) 
    clusters.to_csv(outprefix+"_agglomerative_" + linkage + "_clusters.csv", sep="\t")
    with open(outprefix+"_agglomerative_" + linkage + "_performance.csv", 'w+') as f:
        f.write("computation_time\t"+str(end_time - start_time))
    plt.scatter(matrix.values[:,0],matrix.values[:,1], c=model.labels_, cmap='rainbow')
    plt.suptitle("Agglomerative clustering result (linkage = {})".format(linkage))
    plt.savefig(outprefix+"_agglomerative_" + linkage + "_scatter.png")
    plt.clf()
    

def main():
    parser = argparse.ArgumentParser(description='Agglomerative clustering')
    parser.add_argument('-f', '--file', type=str, required=True,
                    help='Principal components matrix.')
    parser.add_argument('-p', '--preproc', action='store_true', 
                    help="To be specified if the input file requires preprocessing")
    parser.add_argument('-o', '--outprefix', type=str, required=True,
                    help='Output prefix')
    parser.add_argument('-l', '--linkage', type=str, required=True,
                    help='Linkage method')

    args = parser.parse_args()

    if args.preproc:
        matrix_input = pd.read_csv(args.file, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END']).transpose()
    else:
        matrix_input = pd.read_csv(args.file, index_col=0, sep='\t')    

    outprefix = args.outprefix
    linkage = args.linkage

    k = locate_elbow(matrix_input, outprefix, linkage)
    cluster_and_output(k, matrix_input, outprefix, linkage)

if __name__ == "__main__":
    sys.exit(main())