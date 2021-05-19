#!/usr/bin/env python3

import os, sys
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns 
from scipy.stats import pearsonr, gaussian_kde
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

def plot_dist(shscores, outpath):
    sns.set(font_scale=1.5)
    sns.violinplot(data=shscores)
    plt.savefig(os.path.join(outpath, "sh_scores_violin.png"))

def compute_correlations_and_plot3D(shscores, thetas, rec_ps, outpath):
    c1, p1 = pearsonr(shscores, thetas)
    print("Pearson correlation (mean CNV size): c = {}, p={}".format(c1, p1))
    plt.scatter(shscores, thetas)
    plt.suptitle('SHscore vs mean CNV size', fontsize=24)
    plt.title("Pearson correlation coef = " + str(round(c1, 3)) + ", pvalue = " + str(round(p1, 3)), fontsize=18)
    plt.xlabel("SHscore", fontsize=20)
    plt.ylabel("Mean CNV size (expected value)", fontsize=20)
    plt.gcf().set_size_inches(7,7)
    plt.savefig(os.path.join(outpath, "shscores_theta_pearson.png"))
    plt.cla()
    c2, p2 = pearsonr(shscores, rec_ps)
    print("Pearson correlation (mean gained copies): c = {}, p={}".format(c2, p2))
    plt.scatter(shscores, rec_ps)
    plt.suptitle('SHscore vs mean gained copy number', fontsize=24)
    plt.title("Pearson correlation coef = " + str(round(c2, 3)) + ", pvalue = " + str(round(p2, 3)), fontsize=18)
    plt.xlabel("SHscore", fontsize=20)
    plt.ylabel("Mean gained copy numbers (expected value)", fontsize=20)
    plt.gcf().set_size_inches(7,7)
    plt.savefig(os.path.join(outpath, "shscores_p_pearson.png"))
    plt.cla()

    x = thetas
    y = rec_ps
    z = shscores

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    xyz = np.vstack([x,y,z])
    k = gaussian_kde(xyz)(xyz)
    plot = ax.scatter(x, y, z, c=k, s=5)
    #rotate xy plane
    ax.azim = -60
    ax.dist = 10
    ax.elev = 10
    ax.set_xlabel("Mean CNA size", linespacing=3.2)
    ax.set_ylabel("Mean gained copy numbers", linespacing=3.2)
    ax.set_zlabel("SHscore", linespacing=3.2)
    cax = fig.add_axes([0.15, 0.03, 0.03, 0.80])
    fig.colorbar(plot, cax=cax)
    plt.savefig(os.path.join(outpath, "shscores_pearson.png"))
    plt.close("all")

def main():
    parser = argparse.ArgumentParser(description="Compute Pearson correlation coefficient between sh_score and params of simulator (theta, p)")
    parser.add_argument("-i", "--input", required=True, help="Input files", nargs="*")
    parser.add_argument("-o", "--outdir", required=True, help="Outdir path")

    args = parser.parse_args()
    table = pd.DataFrame(columns=["supersample", "theta", "1/p", "sh_score"])

    for f in args.input:
        df = pd.read_csv(f, sep="\t")
        table = table.append(df, ignore_index=True)

    table.to_csv(os.path.join(args.outdir, "shscores.tsv"), sep="\t", index=False)
    plot_dist(table.sh_score.values, args.outdir)
    compute_correlations_and_plot3D(table.sh_score.values, table.theta.values, table["1/p"].values, args.outdir)

if __name__ == "__main__":
    sys.exit(main())