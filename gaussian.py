import time
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.base import ClusterMixin
from sklearn.mixture import GaussianMixture
from yellowbrick.cluster import KElbowVisualizer


#gmm = mixture.GaussianMixture(n_components=4, covariance_type="full").fit(pca_result)
#labels = gmm.predict(pca_result)

import argparse, os, sys

class GMClusters(GaussianMixture, ClusterMixin):

    def __init__(self, n_clusters=2, **kwargs):
        kwargs["n_components"] = n_clusters
        super(GMClusters, self).__init__(**kwargs)

    def fit(self, X):
        super(GMClusters, self).fit(X)
        self.labels_ = self.predict(X)
        return self 

def locate_elbow(matrix, outprefix):
    model = GMClusters()
    visualizer = KElbowVisualizer(model, k=(2,20), force_model=True)
    visualizer.fit(matrix.values)
    visualizer.finalize()
    plt.savefig(outprefix+"_gaussian_elbow.png")
    plt.clf()
    return visualizer.elbow_value_


def cluster_and_output(k, matrix, outprefix):
    start_time = time.time()
    model = GaussianMixture(n_components=k).fit(matrix.values)
    end_time = time.time()
    labels_ = model.predict(matrix.values)
    clusters = pd.DataFrame(labels_, index=matrix.index, columns=["cluster"])
    clusters.to_csv(outprefix+"_gaussian_clusters.csv", sep="\t")
    with open(outprefix+"_gaussian_performance.csv", 'w+') as f:
        f.write("computation_time\t"+str(end_time - start_time))
    plt.scatter(matrix.values[:,0],matrix.values[:,1], c=labels_, cmap='rainbow')
    plt.suptitle("Gaussian Mixture Modelling result")
    plt.savefig(outprefix+"_gaussian_scatter.png")
    plt.clf()
    

def main():
    parser = argparse.ArgumentParser(description='Gaussian Mixture Modelling clustering')
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