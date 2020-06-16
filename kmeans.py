import time
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from yellowbrick.cluster import KElbowVisualizer

import argparse, os, sys

def calinski_harabasz(matrix, outprefix):
    d = {}
    for k in range(2, 21):
        model = KMeans(n_clusters=k).fit(matrix.values)
        labels = model.labels_
        d[k] = metrics.calinski_harabasz_score(matrix.values, labels)
    k = max(d.items(), key=operator.itemgetter(1))[0]
    lists = sorted(d.items())
    x, y = zip(*lists) 
    plt.plot(x, y, 'bx-')
    plt.xlabel('k')
    plt.ylabel('calinski_harabaz_score')
    plt.axvline(x=k, color='k', linestyle='--', label="optK")
    plt.gca().grid(True)
    plt.legend(fontsize=18)
    plt.savefig(outprefix+"_kmeans_calinski_harabasz.png")
    plt.clf()
    return k

def locate_elbow(matrix, outprefix):
    model = KMeans()
    visualizer = KElbowVisualizer(model, k=(2,20))
    visualizer.fit(matrix.values)
    visualizer.finalize()
    plt.savefig(outprefix+"_kmeans_elbow.png")
    plt.clf()
    return visualizer.elbow_value_


def cluster_and_output(k, matrix, outprefix):
    start_time = time.time()
    model = KMeans(n_clusters=k).fit(matrix.values)
    end_time = time.time()
    clusters = pd.DataFrame(model.labels_, index=matrix.index, columns=["cluster"])
    clusters.to_csv(outprefix+"_kmeans_clusters.csv", sep="\t")
    with open(outprefix+"_kmeans_performance.csv", 'w+') as f:
        f.write("computation_time\t"+str(end_time - start_time))
    plt.scatter(matrix.values[:,0],matrix.values[:,1], c=model.labels_, cmap='rainbow')
    plt.suptitle("K-Means clustering result")
    plt.savefig(outprefix+"_kmeans_scatter.png")
    plt.clf()
    

def main():
    parser = argparse.ArgumentParser(description='K-Means clustering')
    parser.add_argument('-f', '--file', type=str, required=True,
                    help='Principal components matrix.')
    parser.add_argument('-p', '--preproc', action='store_true', 
                    help="To be specified if the input file requires preprocessing")
    parser.add_argument('-o', '--outprefix', type=str, required=True,
                    help='Output prefix')

    args = parser.parse_args()
    if args.preproc:
        matrix_input = pd.read_csv(args.file, sep="\t", usecols = lambda column : column not in ['CHR', 'START', 'END']).transpose()
    else:
        matrix_input = pd.read_csv(args.file, index_col=0, sep='\t')
    outprefix = args.outprefix

    k = locate_elbow(matrix_input, outprefix)
    cluster_and_output(k, matrix_input, outprefix)

if __name__ == "__main__":
    sys.exit(main())


