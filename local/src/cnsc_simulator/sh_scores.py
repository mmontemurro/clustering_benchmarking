#!/usr/bin/env python3

import os, sys
import argparse
import pandas as pd
import numpy as np
import phylics
import itertools
import matplotlib.pyplot as plt
from tree_metrics import distanceBetweenNodes
import multiprocessing
from joblib import Parallel, delayed
from scipy.stats import pearsonr
import random
import seaborn as sns

def dist_plot(shscores, outpath):
    sns.set(font_scale=1.5)
    sns.violinplot(data=shscores)
    #textstr = 'Mann-Whitney-U test pval=%.3E' % (p1, )
    #props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    #plt.gca().text(0.05, 0.95, textstr, transform=plt.gca().transAxes, fontsize=10,
    #    verticalalignment='top', bbox=props)
    plt.savefig(os.path.join(outpath, "sh_scores_violin.png"))


def sh_score(ex1, ex2):
    ss = phylics.MultiSample.from_list(ex1, ex2)
    s = ss.SHscore(n_jobs=10)
    return s['score'].values[0]


def main():
    parser = argparse.ArgumentParser(description="Generates combination of trees, computes SHscore and phylo distance between MRCAs and computes correlaiton between the two measurements")
    parser.add_argument("-d", "--directories", required=True, help="List of paths to the dirs containing both the CNV matrix and the tree structure of subtrees", nargs="+")
    parser.add_argument("-c", "--cnv_suffix", required=True, help="Suffix to reach cnv matrix from dir path")
    parser.add_argument("-o", "--outpath", required=True, help="Output directory path")
    parser.add_argument("-j", "--n_jobs", default=1, type=int, help="Number of parallel jobs to launch")

    args = parser.parse_args()

    n_subtrees = len(args.directories)
    
    # compute all possibile combinations of 2 out of n_subtrees
    pairs = list(itertools.combinations(range(n_subtrees), 2))
    #print("Number of combinations = " + str(len(pairs)))
    #pairs = tqdm(pairs)
    #results = Parallel(n_jobs=args.n_jobs)(delayed(shscore_mrca)(args.directories[pair[0]], args.directories[pair[1]], args.cnv_suffix, args.npy_suffix, copy(tree)) for pair in pairs)
    
    df = pd.DataFrame(columns=['sample1', 'sample2', 'sh_score'])                                          
    for i, pair in enumerate(pairs):
        
        # compute the full path of the cnv matrix and of the npy structure of both datasets
        cnv_path1 = os.path.join(args.directories[pair[0]],  args.cnv_suffix)
        cnv_path2 = os.path.join(args.directories[pair[1]],  args.cnv_suffix)

        #retrieve sample names from the last part of the directory path
        sample1 = os.path.basename(os.path.normpath(args.directories[pair[0]]))
        sample2 = os.path.basename(os.path.normpath(args.directories[pair[1]]))

        # compute the SHscore between the two datasets
        cnvs1 = phylics.Sample.from_file(cnv_path1, sample_name=sample1)
        cnvs2 = phylics.Sample.from_file(cnv_path2, sample_name=sample2)
        sh = sh_score(cnvs1, cnvs2)
        print("Iteration n " + str(i) + ': sh_score = ' + str(sh))
        df = df.append({'sample1': sample1, 'sample2':sample2, 'sh_score':sh}, ignore_index=True)
    
    #df = pd.DataFrame(results)
    df.to_csv(os.path.join(args.outpath, "shscores.tsv"), sep="\t", index=False)
    print(df.describe()) 
    third_quart = df.quantile(0.75)
    first_quart = df.quantile(0.25)
    df1 = third_quart - first_quart
    print(df1)

    dist_plot(df["sh_score"].values, args.outpath)

if __name__ == "__main__":
    sys.exit(main())

