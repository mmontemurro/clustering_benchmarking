#!/usr/bin/envs python

import time
import umap
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib import cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from sklearn.cluster import AffinityPropagation
from sklearn.metrics import davies_bouldin_score
import argparse, os, sys
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectFromModel
from sklearn.metrics import pairwise_distances

def cluster_and_output(matrix, outprefix, boundaries, chr_limits, pos_list, chrN_list, chr_boundaries):
    model = AffinityPropagation().fit(matrix.values)
    clusters_df = pd.DataFrame({'signature':model.cluster_centers_.tolist()}) #uns["clusters"]
    clusters_df.index.name = 'label'
    clusters = pd.DataFrame(model.labels_, index=matrix.index, columns=["cluster"]) #obs["cluster"]

    """
    Simpson's Index of Diversity 1 - D

    The value of this index also ranges between 0 and 1, but now, the greater the value, the greater the sample diversity. 
    """
    
    #TODO filter out singletons?
    clusters_df["n_items"] = clusters.groupby("cluster").size()
    clusters_df = clusters_df[clusters_df["n_items"] > 1]


    clusters_df.to_csv(outprefix+"_cluster_signatures.tsv", sep="\t")

    #filter out clusters with < 14 items
    clusters = clusters.groupby("cluster").filter(lambda x: len(x) >= 10)

    clusters.to_csv(outprefix+"_affinity_clusters.csv", sep="\t")

    # get index of non outliers
    non_outliers = clusters.index


    #remove outliers from dataframe
    print("N cells: {} - N outlier cells {}".format(len(matrix), len(matrix) - len(non_outliers)))
    matrix = matrix.loc[matrix.index.isin(non_outliers)]    
        
    labels = clusters["cluster"]
    
    clusters_dict = {k: v for v, k in enumerate(np.unique(labels))}
    
    color_palette = sns.color_palette("hls", len(clusters_dict))

    cluster_colors = [color_palette[clusters_dict[x]] if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in labels]    
    
    colors_dict = {v: k for v, k in enumerate(color_palette)}
    markers = [plt.Line2D([0,0],[0,0],color=color, marker='o', linestyle='') for color in colors_dict.values() ]

    standard_embedding = umap.UMAP(random_state=42).fit_transform(matrix.values)
    plt.scatter(standard_embedding[:,0],standard_embedding[:,1], 
                linewidth=0,
                c=cluster_colors,
                #s=0.9,
                alpha=0.5)

    plt.legend(markers, colors_dict.keys(), numpoints=1, fontsize='large')

    plt.suptitle("Affinity propagation clusters", fontsize=16, fontweight='bold')
    plt.gcf().set_size_inches(12,12)
    plt.savefig(outprefix+"_affinity_scatter.png")
    plt.clf()
    
    # Use a random forest to compute relevance of features for clustering
    rf = RandomForestClassifier(n_estimators = 100)
    rf.fit(matrix, labels)
    sel = SelectFromModel(rf, threshold=0.001)
    sel.fit(matrix.values, labels)
    selected_feat= matrix.columns[(sel.get_support())]
    boundaries.loc[selected_feat].to_csv(outprefix+'_important_features.csv')

    plt.style.use('fivethirtyeight')
    importances = list(rf.feature_importances_)
    x_values = list(range(len(importances)))
    plt.bar(x_values, importances, orientation = 'vertical')
    #plt.xticks(x_values, feature_list, rotation='vertical')
    plt.gca().xaxis.set_major_locator(ticker.FixedLocator(pos_list))
    plt.gca().xaxis.set_major_formatter(ticker.FixedFormatter((chrN_list)))
    plt.gca().tick_params(axis='x', rotation=0, labelsize=16, which='major')

    plt.gca().xaxis.set_minor_locator(ticker.FixedLocator(chr_boundaries))
    plt.gca().tick_params(axis='x', length=20, which='minor')

    plt.ylabel('Importance')
    plt.xlabel('Variable')
    plt.title('Variable Importances')
    plt.gcf().set_size_inches(30, 20)
    plt.savefig(outprefix+"_feature_importances.png")
    plt.clf()

    print(sel.threshold_)

    matrix['cluster'] = labels

    matrix = matrix.sort_values(by='cluster')
    res_df = matrix[["cluster"]]
    res_df.index.rename("cellid", inplace=True)
    res_df.to_csv(outprefix+"_clusters.csv")


    sorted_labels = matrix['cluster'].values
    #color_palette = sns.color_palette("Spectral", len(np.unique(cnvs['cluster'][cnvs['cluster'] >= 0].values)))
    cluster_colors = [color_palette[clusters_dict[x]] if x >= 0
                      else (0.5, 0.5, 0.5)
                      for x in sorted_labels]
    sns.set(font_scale=2)
    matrix = matrix.drop(['cluster'], axis=1)
    cbar_kws={"ticks":np.arange(0,7,1)}
    h = sns.clustermap(matrix,
            row_cluster=False,
            col_cluster=False,
            yticklabels = False,
            row_colors=cluster_colors,
            cmap='RdBu_r',
            vmin=0, vmax=6,
            center = 2,
            #norm=divnorm,
            cbar_kws=cbar_kws)

    for label in np.unique(list(colors_dict.keys())):
        h.ax_col_dendrogram.bar(0, 0, color=color_palette[label],label=label, linewidth=0)
    h.cax.set_position([0.2, .2, .03, .45])
    h.ax_col_dendrogram.legend(loc="center",  ncol=len(np.unique(list(colors_dict.keys()))))
    ax = h.ax_heatmap
    #place vertical lines to identify chromosomes
    for pos in chr_limits:
        ax.axvline(x=pos, color='black')
    
    
    #place chromosome ticks at the right position
    ax.xaxis.set_major_locator(ticker.FixedLocator(pos_list))
    ax.xaxis.set_major_formatter(ticker.FixedFormatter((chrN_list)))
    ax.tick_params(axis='x', rotation=0, labelsize=16, which='major')

    ax.xaxis.set_minor_locator(ticker.FixedLocator(chr_boundaries))
    ax.tick_params(axis='x', length=20, which='minor')

    ax.set_xlabel("chromosomes", fontsize=16)
    ax.set_ylabel("cells", fontsize=16)

    #plt.gcf().suptitle('Cnv profiles', fontsize=16, fontweight='bold' )
    plt.gcf().set_size_inches(37, 21)
    plt.savefig(outprefix+'_clusters_heatmap.png', bbox_inches = 'tight')
    plt.clf()
    
    davies_bouldin = davies_bouldin_score(matrix, labels)
    print(davies_bouldin)

def main():
    parser = argparse.ArgumentParser(description='Affinity propagation clustering')
    parser.add_argument('-f', '--file', type=str, required=True,
                    help='input matrix.')
    #parser.add_argument("-s", "--statsfile", type=str, required=True,
    #                help='Stats file.')
    parser.add_argument('-o', '--outprefix', type=str, required=True,
                    help='Output prefix')

    args = parser.parse_args()

    outprefix = args.outprefix
    df = pd.read_csv(args.file, sep="\t")
    matrix_input = df.drop(['CHR', 'START', 'END'], axis=1).transpose()
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


    #stats_f = pd.read_csv(args.statsfile, sep="\t")
    #threshold = stats_f['Disp'].quantile([0.9]).values[0]
    #noisy_cells = stats_f[stats_f['Disp'] >= threshold].index.values
    #noisy = df_stats[df_stats['Disp'] > threshold].index.values

    #sns.distplot(stats_f['Disp'], rug=True)
    #plt.axvline(x=threshold, color='r', linestyle='--', label="90th percentile")
    #plt.legend()
    #plt.savefig(outprefix + "_noisy_distr.png")
    plt.clf()

    #print("N cells: {} - N noisy cells {}".format(len(matrix_input), len(noisy_cells)))
    #matrix_input = matrix_input.loc[~matrix_input.index.isin(noisy_cells)]


    #k = locate_elbow(matrix_input, outprefix)

    # cell id dictionary
    #cell_id_dict = dict(zip(list(range(len(matrix_input))), matrix_input.index ))
    #cell_ids = matrix_input.index

    cluster_and_output(matrix_input, outprefix, boundaries, chr_limits, pos_list, chrN_list, chr_boundaries)
    
    



if __name__ == "__main__":
    sys.exit(main())
