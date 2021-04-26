import os, sys, pickle, random
import numpy as np 
import argparse

parser = argparse.ArgumentParser(description="Split cell tree in top second level cell lineages")
parser.add_argument("-i", "--input_file", required=True, help="Subtrees dictionary (pickle format)")
parser.add_argument("-o", "--out_prefix", required=True, help="Outpath prefix")
parser.add_argument("-m", "--max_iterations", required=True, type=int, help="Maximum number of iterations to find subsamples and avoiding infine loops")
parser.add_argument("-s", "--n_subtrees", required=True, type=int, help="Number of second level cell lineages per clade to save")

args = parser.parse_args()

f = open(args.input_file, "rb")
p = pickle.load(f)
f.close()

subtrees = {}
n_leaves = []
for k in p.keys():
    subtrees_l1 = p[k]["subtree_l1"]
    if subtrees_l1 not in subtrees.keys():
        subtrees[subtrees_l1] = []
    subtrees[subtrees_l1].append(p[k])
    n_leaves.append(p[k]["n_leaves"])
        
first_quartile = np.quantile(n_leaves, 0.25)
third_quartile = np.quantile(n_leaves, 0.75)
for idx, k in enumerate(subtrees.keys()):
#iterate on the top level clades
    print("Clade_n = " + str(idx))
    #max_iterations = args.max_iterations #this variable avoid infinite loops, in case the desired number of subsamples is not found
    for i in range(args.n_subtrees):
        #max_iterations -= 1
        if len(subtrees[k]) == 1:
            tree = subtrees[k][0]["tree"]
            leaves = subtrees[k][0]["leaves"]
            n_leaves = subtrees[k][0]["n_leaves"]
        else:
            # sample a random number in the interval [0, (number of subtrees for this clade - 1)]
            n = random.randint(0, len(subtrees[k])-1)
            #if the current tree has a number of leaves which stands 
            # in between the second or the third quartile, take it 
            st = subtrees[k][n]
            """
            if st["n_leaves"] > first_quartile and  st["n_leaves"] < third_quartile:
            """    
            tree = st["tree"]
            leaves = st["leaves"]
            n_leaves = st["n_leaves"]
            """
            else:
                i -= 1
                continue
            """
        outdir = args.out_prefix + str(idx) + "_" + str(i)
        np.save(os.path.join(outdir, "from_first_step.subtree.npy"), tree)
        np.save(os.path.join(outdir, "from_first_step.leaf_index.npy"), leaves)
        print("\tSaving subtree: n_leaves " + str(n_leaves))
        