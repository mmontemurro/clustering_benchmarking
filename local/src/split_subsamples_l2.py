import os, sys, pickle, random
import numpy as np 
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description="Split cell tree in top second level cell lineages")
parser.add_argument("-i", "--input_file", required=True, help="Subtrees dictionary (pickle format)")
parser.add_argument("-s", "--stats_file", required=True, help="Subtrees statistics")
parser.add_argument("-o", "--out_prefix", required=True, help="Outpath prefix")
#parser.add_argument("-m", "--max_iterations", required=True, type=int, help="Maximum number of iterations to find subsamples and avoiding infine loops")
#parser.add_argument("-s", "--n_subtrees", required=True, type=int, help="Number of second level cell lineages per clade to save")

args = parser.parse_args()

# read data structure containing subtrees data
f = open(args.input_file, "rb")
subtrees = pickle.load(f)
f.close()

#read the file containing the statistics on subtrees (n_leaves)
df = pd.read_csv(args.stats_file, sep="\t")
first_quartile = df.n_leaves.quantile(0.25)
third_quartile = df.n_leaves.quantile(0.75)

#keep on trees with first_quartile <= tree.n_leaves <= third_quartile (100 trees)
df = df[(df.n_leaves >= first_quartile ) & (df.n_leaves <= third_quartile) ]
keep_trees = df.subtree_id.values 

for count, k in enumerate(keep_trees):
    tree = subtrees[k]["tree"]
    leaves = subtrees[k]["tree"][0].leaves
    outdir = args.out_prefix + str(count) 
    np.save(os.path.join(outdir, "from_first_step.tree.npy"), tree)
    np.save(os.path.join(outdir, "from_first_step.leaf_index.npy"), leaves)
        