#!/usr/bin/env python
import os
import sys
import argparse
import numpy as np
import pandas as pd 

def extract_clone_cells(segcopy, leaves, clones, c):
    cols = ["CHR", "START", "END"]
    for idx, cc in enumerate(clones):
        if cc == c:
            cols.append("cell" + str(leaves[idx]))
    return segcopy[cols]
    

def main():
    parser = argparse.ArgumentParser(description="Split ground truth by clonal subtree")
    parser.add_argument("-s", "--segcopy", type=str)
    parser.add_argument("-n1", "--npy_leaves", type=str)
    parser.add_argument("-n2", "--npy_subtrees", type=str)
    parser.add_argument("-o", "--outpath", type=str)

    args = parser.parse_args()

    segcopy = pd.read_csv(args.segcopy, sep="\t")
    leaves = np.load(args.npy_leaves)
    clones = np.load(args.npy_subtrees)

    clones_id = np.unique(clones)
    for c in clones_id:
        segcopy_part = extract_clone_cells(segcopy, leaves, clones, c)
        segcopy_part.to_csv(os.path.join(args.outpath, "SegCopy_" + str(c)), sep="\t", index=False)


if __name__ == "__main__":
    sys.exit(main())