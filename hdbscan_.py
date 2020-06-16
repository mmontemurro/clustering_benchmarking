#!/usr/bin/env python

import time
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import hdbscan

import argparse, os, sys, operator

def cluster_size(matrix):
    model_dict = {}
    for min_size in range(3,50+1):
            model = hdbscan.HDBSCAN(min_cluster_size=min_size, min_samples=1, cluster_selection_epsilon=0.5).fit(matrix.values) 
            count = 0
            for p in model.probabilities_:
                if p < 0.05:
                    count = count +1
            model_dict[min_size] = count
    min_cluster_size=min(model_dict.items(), key=operator.itemgetter(1))[0]
    return min_cluster_size

def cluster_and_output(min_cluster_size, matrix, outprefix):
    start_time = time.time()
    #cluster_selection_method='eom', cluster_selection_epsilon=0.5, metric='l1'
    model = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=1, cluster_selection_epsilon=0.5).fit(matrix.values)
    end_time = time.time()
    clusters = pd.DataFrame(model.labels_, index=matrix.index, columns=["cluster"])
    clusters.to_csv(outprefix+"_hdbscan_clusters.csv", sep="\t")
    with open(outprefix+"_hdbscan_performance.csv", 'w+') as f:
        f.write("computation_time\t"+str(end_time - start_time))
    color_palette = sns.color_palette("rainbow", len(np.unique(model.labels_)))
    cluster_colors = [color_palette[x] if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in model.labels_]
    plt.scatter(matrix.values[:,0],matrix.values[:,1], c=cluster_colors)
    plt.suptitle("HDBSCAN clustering result")
    plt.savefig(outprefix+"_hdbscan_scatter.png")
    plt.clf()
    

def main():
    parser = argparse.ArgumentParser(description='HDBSCAN clustering')
    parser.add_argument('-f', '--file', type=str, required=True,
                    help='Principal components matrix.')
    parser.add_argument('-p', '--preproc', action='store_true', 
                    help="To be specified if the input file requires preprocessing")
    parser.add_argument('-o', '--outprefix', type=str, required=True,
                    help='Output prefix')

    args = parser.parse_args()
    if args.preproc:
        matrix_input = pd.read_csv(input_file, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END']).transpose()
    else:
        matrix_input = pd.read_csv(args.file, index_col=0, sep='\t')
    outprefix = args.outprefix

    min_cluster_size = cluster_size(matrix_input)
    cluster_and_output(min_cluster_size, matrix_input, outprefix)

if __name__ == "__main__":
    sys.exit(main())


