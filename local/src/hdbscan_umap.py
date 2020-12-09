#!/usr/bin/env python

import sys
import operator
import argparse
import numpy as np
from math import sqrt
import pandas as pd
import seaborn as sns
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import umap
import hdbscan
import matplotlib.colors as colors
from statistics import mean
from itertools import combinations
from sklearn.metrics import davies_bouldin_score
from sklearn.feature_selection import VarianceThreshold

def heatmap(cnvs, boundaries, method, metric, outdir, clusters, verbose, sample=None):
    divnorm = colors.DivergingNorm(vmin=0, vcenter=2, vmax=12)

    chr_limits = boundaries.index[boundaries['END'].isin(boundaries.groupby('CHROM', sort=False)['END'].max().values)].tolist()
    chr_boundaries = np.append(0, chr_limits)
    chr_list = boundaries['CHROM'].unique().tolist()
    chrN_list = []

    for x in chr_list:
        x = x[3:] #remove 'chr' for readability
        chrN_list.append(x)

    #compute the position where chromosome labels will be placed on the plots
    start = 0
    pos_list = []
    for end in chr_limits:
        pos_list.append((start+end)/2)
        start = end+1

    if clusters:
        yticklabels = True
    else:
        yticklabels = False

    cbar_kws={"ticks":np.arange(0,7,1)}
    h = sns.clustermap(cnvs, method=method, metric=metric, col_cluster=False, yticklabels = yticklabels,  cmap='RdBu_r', vmin=0, vmax=6,center=2, cbar_kws=cbar_kws)
    Z = h.dendrogram_row.linkage
        


    """
    if clusters == False:
        h.cax.set_position([0.05, .2, .03, .45])
        c, coph_dist = cophenet(Z,pdist(cnvs, metric))
        textstr = '(cophenet coeff= ' + str(c)+')'
        print_msg("cophenet coefficient: " + str(c), 1, verbose)
    """

    ax = h.ax_heatmap
    #place vertical lines to identify chromosomes
    for pos in chr_limits:
        ax.axvline(x=pos, color='black')

    #place chromosome ticks at the right position
    ax.xaxis.set_major_locator(ticker.FixedLocator((pos_list)))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter((chrN_list)))
    ax.tick_params(axis='x', rotation=0, labelsize=14)

    ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_boundaries))
    ax.tick_params(axis='x', length=20, which='minor')

    ax.set_xlabel("Chromosomes", fontweight='bold', fontsize=14)
    
    if clusters:
        ax.set_ylabel("Clusters", fontweight='bold', fontsize=14)
    else:
        ax.set_ylabel("cells", fontsize=14, fontweight='bold')

    plt.gcf().set_size_inches(37, 21)

    if clusters:
        plt.gcf().suptitle("Clusters mean CNV heatmap", fontsize=16, fontweight='bold')
        plt.savefig(outdir+".dendr.heatmap.png")
    else:
        plt.gcf().suptitle("CNV heatmap", fontsize=16, fontweight='bold')
        plt.savefig(outdir+".dendr.heatmap.png")
    plt.clf()
    return Z



def main():

    parser = argparse.ArgumentParser(description="hdbscan")


    parser.add_argument("input", metavar='SegCopy', action='store',
            help='cnvs file.', type=str)

    #parser.add_argument("input2", metavar='SegStats', action='store',
    #        help='cell info file.', type=str)

    parser.add_argument("prefix", metavar='out/path/prefix', action='store',
            help='Output directory prefix.', type=str)

    parser.add_argument("--algorithm", metavar='algorithm', action='store', choices=['best', 'generic', 'prims_kdtree', 'prims_balltree', 'boruvka_kdtree', 'boruvka_balltree'], default='best',
            help='Number of principal components.', type=str)

    parser.add_argument("--min_cluster_size", metavar='N', action='store',
            help='Min cluster size (hdbscan parameter).', type=str)


    args=parser.parse_args()

    input_file = args.input
    #input2 = args.input2
    prefix = args.prefix
    algorithm = args.algorithm
    min_cluster_size = None

    if args.min_cluster_size:
        min_cluster_size = min_cluster_size

    df = pd.read_csv(input_file, sep="\t")

    cnvs = df.drop(['CHR', 'START', 'END'], axis=1).transpose()
    boundaries = df[['CHR', 'START', 'END']].copy()
    chr_limits = boundaries.index[boundaries['END'].isin(boundaries.groupby('CHR', sort=False)['END'].max().values)].tolist()
    chr_boundaries = np.append(0, chr_limits)
    chr_list = boundaries['CHR'].unique().tolist()
    chrN_list = []

    for x in chr_list:
        x = x[3:] #remove 'chr' for readability
        chrN_list.append(x)

    #compute the position where chromosome labels will be placed on the plots
    start = 0
    pos_list = []
    for end in chr_limits:
        pos_list.append((start+end)/2)
        start = end+1


    # Cells with high MAD are filtered out (noisy)
    
    #df_stats = pd.read_csv(input2, sep="\t")
    #threshold = df_stats['Disp'].quantile([0.9]).values[0]
    #noisy_cells = df_stats[df_stats['Disp'] >= threshold].index.values
    #noisy = df_stats[df_stats['Disp'] > threshold].index.values

    #sns.distplot(df_stats['Disp'], rug=True)
    #plt.axvline(x=threshold, color='r', linestyle='--', label="90th percentile")
    #plt.legend()
    #plt.savefig(prefix + ".noisy.distr.png")
    #plt.clf()

    #print("N cells: {} - N noisy cells {}".format(len(cnvs), len(noisy_cells)))
    #cnvs = cnvs.loc[~cnvs.index.isin(noisy_cells)]
    
    # cell id dictionary
    cell_id_dict = dict(zip(list(range(len(cnvs))), cnvs.index ))
    cell_ids = cnvs.index
        
    #heatmap
    #heatmap(cnvs, boundaries, 'single', 'euclidean', prefix, False, True)

    #umap
    N=len(cnvs)
    selector = VarianceThreshold() #filters out features with 0 variance from cell to cell (not informative for clustering)
    cnvs_filtered = pd.DataFrame(data=selector.fit_transform(cnvs), index=cnvs.index)
    
    standard_embedding = umap.UMAP( 
                                    #n_neighbors=15,
                                    #min_dist=0.0,

                                    random_state=42
            
                                ).fit_transform(cnvs_filtered.values)
    

    #hdbscan
    
    
    if min_cluster_size == None:
        clusterer_dict = {}
        for min_size in range(3,50+1):
            clusterer = hdbscan.HDBSCAN(min_cluster_size=min_size,  min_samples=1, cluster_selection_method='eom').fit(standard_embedding) 
            count = 0
            for p in clusterer.probabilities_:
                if p < 0.05:
                    count = count +1
            clusterer_dict[min_size] = count

        min_cluster_size=min(clusterer_dict.items(), key=operator.itemgetter(1))[0]
        print("optimal min cluster size: " + str(min_cluster_size))
    


    #clusterer = hdbscan.HDBSCAN(metric='hamming').fit(pca_result)
    clusterer = hdbscan.HDBSCAN(min_cluster_size=min_cluster_size, min_samples=1, cluster_selection_epsilon=0.5).fit(standard_embedding)
    
    clustered = (clusterer.labels_>=0)
    unclustered = (clusterer.labels_ < 0)
    clustered_perc = np.sum(clustered) / cnvs.values.shape[0]
    print("Percentage of clustered data: " + str(clustered_perc))



    #outliers detection
    threshold = pd.Series(clusterer.outlier_scores_).quantile(0.90)
    
    sns.distplot(clusterer.outlier_scores_[np.isfinite(clusterer.outlier_scores_)], rug=True)
    plt.axvline(x=threshold, color='r', linestyle='--', label="90th percentile")
    plt.legend()
    plt.savefig(prefix+".outlier.distr.png")
    plt.clf()
    

    outliers = np.where(clusterer.outlier_scores_ > threshold)[0]
    plt.scatter(standard_embedding[: , 0], standard_embedding[: , 1], s=50, linewidth=0, c='gray', alpha=0.25)
    plt.scatter(standard_embedding[outliers , 0], standard_embedding[outliers , 1], s=50, linewidth=0, c='red', alpha=0.5)
    plt.gcf().set_size_inches(12,12)
    plt.savefig(prefix+".outliers.scatter.png")
    plt.clf()


    # get names of outlier cells
    outlier_cells = cell_ids[outliers]

    # get index of non-outliers
    non_outliers = [i for i in range(len(cnvs)) if i not in outliers]
    
    #remove outliers from dataframe
    print("N cells: {} - N outlier cells {}".format(len(cnvs), len(outlier_cells)))
    cnvs = cnvs.loc[~cnvs.index.isin(outlier_cells)]

    #remove  labels of outliers
    labels = [clusterer.labels_[i] for i in range(len(clusterer.labels_)) if i not in outliers]

    #remove outliers from umap res
    outliers_embedding = standard_embedding[outliers]
    standard_embedding = standard_embedding[non_outliers]
    
    # hdbscan results

    color_palette = sns.color_palette("hls", len(np.unique(labels)))
    cluster_colors = [color_palette[x] if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in labels]
    colors_dict = {v: k for v, k in enumerate(color_palette)}
    markers = [plt.Line2D([0,0],[0,0],color=color, marker='o', linestyle='') for color in colors_dict.values() ]
    #print(standard_embedding)
    plt.scatter(standard_embedding[: , 0], standard_embedding[: , 1],
                linewidth=0,
                c=cluster_colors,
                #s=0.9,
                alpha=0.5)
    plt.legend(markers, colors_dict.keys(), numpoints=1, fontsize='large')
    #plt.scatter(outliers_embedding[:, 0], outliers_embedding[: , 1], s=50, linewidth=0, c='gray', alpha=0.5)
    plt.gcf().suptitle('hdbscan clusters', fontsize=16, fontweight='bold' )
    plt.gcf().set_size_inches(12,12)
    plt.savefig(prefix+'.hdbscan.clusters.png')
    plt.clf()

    

    #hdbscab+cnv profiles
    cnvs['cluster'] = labels

    cnvs = cnvs.sort_values(by='cluster')
    
    res_df = cnvs[["cluster"]]
    res_df.index.rename("cellid", inplace=True)
    res_df.to_csv(prefix+".clusters.csv")

    
    sorted_labels = cnvs['cluster'].values
    #color_palette = sns.color_palette("Spectral", len(np.unique(cnvs['cluster'][cnvs['cluster'] >= 0].values)))
    cluster_colors = [color_palette[x] if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in sorted_labels]
    #print(row_colors)
    
    cnvs = cnvs.drop(['cluster'], axis=1)
    cbar_kws={"ticks":np.arange(0,7,1)}
    
    sns.set(font_scale=2)
    h = sns.clustermap(cnvs,
            row_cluster=False,
            col_cluster=False,
            yticklabels = False,
            row_colors=cluster_colors,
            cmap='RdBu_r',
            vmin=0, vmax=6,
            center = 2,
            #norm=divnorm,
            cbar_kws=cbar_kws)

    h.cax.set_position([0.2, .2, .03, .45])
    for label in np.unique(sorted_labels):
        h.ax_col_dendrogram.bar(0, 0, color=color_palette[label],label=label, linewidth=0)
    h.ax_col_dendrogram.legend(loc="center",  ncol=len(sorted_labels))

    ax = h.ax_heatmap
    #place vertical lines to identify chromosomes
    for pos in chr_limits:
            ax.axvline(x=pos, color='black')

    #place chromosome ticks at the right position
    ax.xaxis.set_major_locator(ticker.FixedLocator(pos_list))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter((chrN_list)))
    ax.tick_params(axis='x', rotation=0, labelsize=14, which='major')

    ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_boundaries))
    ax.tick_params(axis='x', length=20, which='minor')

    ax.set_xlabel("chromosomes", fontsize=16, fontweight='bold')
    ax.set_ylabel("cells", fontsize=16, fontweight='bold')
       
    #plt.gcf().suptitle('Cnv profiles', fontsize=16, fontweight='bold' )
    plt.gcf().set_size_inches(37, 21)
    plt.savefig(prefix+'.clusters.heatmap.png', bbox_inches = 'tight')
    plt.clf()


    davies_bouldin = davies_bouldin_score(cnvs, labels)
    print(davies_bouldin)

if __name__ == "__main__":
    sys.exit(main())
