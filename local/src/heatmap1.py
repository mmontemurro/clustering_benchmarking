import sys, os
import argparse
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib 
import matplotlib.cm as cm
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import matplotlib.colors as colors


def heatmap(cnvs, boundaries, method, metric, outprefix):
    divnorm = colors.DivergingNorm(vmin=0, vcenter=2, vmax=12)  
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
    cbar_kws={"ticks":np.arange(0,13,1)}
    h = sns.clustermap(cnvs, method=method, metric=metric, row_cluster=False, col_cluster=False, yticklabels = False,  cmap='RdBu_r', vmin=0, vmax=12,norm=divnorm, cbar_kws=cbar_kws)
    h.cax.set_position([0.05, .2, .03, .45])
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
    ax.set_xlabel("Chromosomes", fontweight='bold', fontsize=18)
    ax.set_ylabel("cells", fontsize=18, fontweight='bold')
    plt.gcf().set_size_inches(37, 21)   
    plt.gcf().suptitle("CNV heatmap", fontsize=24, fontweight='bold')
    plt.savefig(outprefix+"_heatmap.png")
    plt.clf()
    

def main():

    parser = argparse.ArgumentParser(description="Draw heatmap.")

    parser.add_argument("input", metavar="SegCopy", action="store", type=str, help="CNV file")
    parser.add_argument("outdir", metavar="out_dir_path", action="store", type=str, help="Output directory path")
    parser.add_argument("-o", "--outprefix", metavar="out_prefix", action="store", type=str, help="Output prefix")

    args = parser.parse_args()

    segcopy = args.input
    outdir = args.outdir
    outprefix = args.outprefix if args.outprefix else ""

    df = pd.read_csv(segcopy, sep="\t")

    cnvs = df.drop(['CHR', 'START', 'END'], axis=1).transpose()
    boundaries = df[['CHR', 'START', 'END']].copy()

    heatmap(cnvs, boundaries, 'complete', 'cityblock', os.path.join(outdir, outprefix))

if __name__ == "__main__":
    sys.exit(main())