import os, sys, pickle, random
import numpy as np 
import argparse

parser = argparse.ArgumentParser(description="Split cell tree in top level cell lineages")

parser.add_argument("-i", "--input_file", required=True, help="Subtrees dictionary (pickle format)")
parser.add_argument("-n", "--n_cells", required=True, type=int, help="Number of final leaves to be subsampled")
parser.add_argument("-o", "--out_prefix", required=True, help="Oupath prefix")

args = parser.parse_args()

f = open(args.input_file, "rb")
p = pickle.load(f)
f.close()

for k in p.keys():
    tree = p[k]["tree"]
    leaves = p[k]["leaves"]
    leaves_ = random.sample(leaves, args.n_cells)
    np.save(os.path.join(args.out_prefix + str(k), "from_first_step.subtree.npy"), tree)
    np.save(os.path.join(args.out_prefix + str(k), "from_first_step.leaf_index.npy"), leaves_)