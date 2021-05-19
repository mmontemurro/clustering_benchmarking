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
from copy import deepcopy
from tqdm import tqdm
import random
from scipy.stats import pearsonr

def correlate_and_plot(shscores, mrca_dists, outpath):
    coef, p = pearsonr(shscores, mrca_dists)
    plt.scatter(shscores, mrca_dists)
    plt.suptitle('SHscore vs MRCA distance', fontsize=24)
    plt.title("Pearson correlation coefficient = {:.3}, pval = {:.3g}".format(coef, p), fontsize=18)
    plt.xlabel("SHscore", fontsize=20)
    plt.ylabel("MRCA distance", fontsize=20)
    plt.gcf().set_size_inches(7,7)
    plt.savefig(os.path.join(outpath, "shscores_mrcadist_pearson.png"))
    plt.cla()
    """
    x_min = np.min(x)
    x_max = np.max(x)
    y_min = np.min(y)
    y_max = np.max(y)
    x_bins = np.linspace(x_min,x_max,50)
    y_bins = np.linspace(y_min,y_max,20)
    data , x_e, y_e = np.histogram2d( x, y, bins = [x_bins, y_bins], density = True )
    z = interpn( ( 0.5*(x_e[1:] + x_e[:-1]) , 0.5*(y_e[1:]+y_e[:-1]) ) , data , np.vstack([x,y]).T , method = "splinef2d", bounds_error = False)
    z[np.where(np.isnan(z))] = 0.0
    idx = z.argsort()
    x, y, z = x[idx], y[idx], z[idx]
    plt.scatter( x, y, c=z)
    norm = Normalize(vmin = np.min(z), vmax = np.max(z))
    plt.colorbar(cm.ScalarMappable(norm = norm), ax=plt.gca())
    plt.colorbar()
    plt.suptitle('SHscore vs MRCA distance', fontsize=24)
    plt.title("Pearson correlation coefficient = {:.3}, pval = {:.3g}".format(coef, p), fontsize=18)
    plt.xlabel("SHscore", fontsize=20)
    plt.ylabel("MRCA distance", fontsize=20)
    plt.gcf().set_size_inches(7,7)
    plt.savefig("shscores_mrcadist_pearson3.png")
    plt.close("all")
    """
    return coef, p

def mrca_distance(tree, subtree1, subtree2):
    root1 = subtree1[0].id 
    root2 = subtree2[1].id
    d = distanceBetweenNodes(tree, root1, root2)
    return d

def sh_score(ex1, ex2):
    ss = phylics.MultiSample.from_list(ex1, ex2)
    s = ss.SHscore(n_jobs=10)
    return s['score'].values[0]

def shscore_mrca(dir1, dir2, cnv_suffix, npy_suffix, tree):
    # compute the full path of the cnv matrix and of the npy structure of both datasets
    cnv_path1 = os.path.join(dir1, cnv_suffix)
    npy_path1 = os.path.join(dir2, npy_suffix)
    print(cnv_path1)
    print(npy_path1)

    cnv_path2 = os.path.join(dir2,  cnv_suffix)
    npy_path2 = os.path.join(dir2,  npy_suffix)

    #retrieve sample names from the last part of the directory path
    sample1 = os.path.basename(os.path.normpath(dir1))
    sample2 = os.path.basename(os.path.normpath(dir2))

    tree1 = np.load(npy_path1, allow_pickle=True)
    tree2 = np.load(npy_path2, allow_pickle=True)

    # compute the phylogenetic distance between the MRCA of the two datasets
    mrca_dist = mrca_distance(tree, tree1, tree2)

    # compute the SHscore between the two datasets
    cnvs1 = phylics.Sample.from_file(cnv_path1, sample_name=sample1)
    cnvs2 = phylics.Sample.from_file(cnv_path2, sample_name=sample2)
    sh = sh_score(cnvs1, cnvs2)
    #print('mrca_distance = ' + str(mrca_dist) + ' - sh_score = ' + str(sh))
    res = {'sample1': sample1, 'sample2':sample2, 'mrca_distance':mrca_dist, 'sh_score':sh}
    return res

def main():
    parser = argparse.ArgumentParser(description="Generates combination of trees, computes SHscore and phylo distance between MRCAs and computes correlaiton between the two measurements")
    parser.add_argument("-t", "--tree", required=True, help="Path to the npy of the parental tree")
    parser.add_argument("-d", "--directories", required=True, help="List of paths to the dirs containing both the CNV matrix and the tree structure of subtrees", nargs="+")
    parser.add_argument("-c", "--cnv_suffix", required=True, help="Suffix to reach cnv matrix from dir path")
    parser.add_argument("-n", "--npy_suffix", required=True, help="Suffix to reach subtree npy from dir path")
    parser.add_argument("-o", "--outpath", required=True, help="Output directory path")
    parser.add_argument("-j", "--n_jobs", default=1, type=int, help="Number of parallel jobs to launch")

    args = parser.parse_args()

    tree = np.load(args.tree, allow_pickle=True)
    n_subtrees = len(args.directories)
    
    # compute all possibile combinations of 2 out of n_subtrees
    pairs = list(itertools.combinations(range(n_subtrees), 2))
    print("Number of combinations = " + str(len(pairs)))
    #pairs = tqdm(pairs)
    #results = Parallel(n_jobs=args.n_jobs)(delayed(shscore_mrca)(args.directories[pair[0]], args.directories[pair[1]], args.cnv_suffix, args.npy_suffix, tree) for pair in pairs)
    subset = random.sample(pairs, 500)
    df = pd.DataFrame()                                          
    for i, pair in enumerate(subset):
        
        # compute the full path of the cnv matrix and of the npy structure of both datasets
        cnv_path1 = os.path.join(args.directories[pair[0]],  args.cnv_suffix)
        npy_path1 = os.path.join(args.directories[pair[0]],  args.npy_suffix)

        cnv_path2 = os.path.join(args.directories[pair[1]],  args.cnv_suffix)
        npy_path2 = os.path.join(args.directories[pair[1]],  args.npy_suffix)

        #retrieve sample names from the last part of the directory path
        sample1 = os.path.basename(os.path.normpath(args.directories[pair[0]]))
        sample2 = os.path.basename(os.path.normpath(args.directories[pair[1]]))

        tree1 = np.load(npy_path1, allow_pickle=True)
        tree2 = np.load(npy_path2, allow_pickle=True)

        # compute the phylogenetic distance between the MRCA of the two datasets
        mrca_dist = mrca_distance(tree, tree1, tree2)

        # compute the SHscore between the two datasets
        cnvs1 = phylics.Sample.from_file(cnv_path1, sample_name=sample1)
        cnvs2 = phylics.Sample.from_file(cnv_path2, sample_name=sample2)
        sh = sh_score(cnvs1, cnvs2)
        print("Iteraion n " + str(i) + ': mrca_distance = ' + str(mrca_dist) + ' - sh_score = ' + str(sh))
        df = df.append({'sample1': sample1, 'sample2':sample2, 'mrca_distance':mrca_dist, 'sh_score':sh}, ignore_index=True)
    
    #f = pd.DataFrame(results)
    df.to_csv(os.path.join(args.outpath, "shscores_mrcadist.tsv"), sep="\t", index=False)
    c, p = correlate_and_plot(df["sh_score"].values, df["mrca_distance"].values, args.outpath)
    print("Pearson correlation coefficient = " + str(c) + ", pvalue = " + str(p))

if __name__ == "__main__":
    sys.exit(main())

