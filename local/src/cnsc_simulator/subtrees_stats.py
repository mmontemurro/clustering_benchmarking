#!/usr/bin/env python3

import sys, os
import pickle
import argparse
import pandas as pd
import gen_tree


def subtrees_stats(trees_d, level="l1"):
    df = pd.DataFrame()
    if level == "l1":
        for k in trees_d.keys():
            root_id = trees_d[k][0].id
            n_leaves = len(trees_d[k][0].leaves)
            df = df.append({'subtree_id':k, 'root_id':root_id, 'n_leaves':n_leaves}, ignore_index=True)
    else:
        for k in trees_d.keys():
            subtree_l1 = trees_d[k]["subtree_l1"]
            root_id = trees_d[k]["tree"][0].id
            n_leaves = len(trees_d[k]["tree"][0].leaves)
            df = df.append({'subtree_id':k, 'subtree_l1':subtree_l1, 'root_id':root_id, 'n_leaves':n_leaves}, ignore_index=True)
    return df

def main():
    parser = argparse.ArgumentParser(description="Compute stats on subtrees")
    parser.add_argument("-t", "--subtrees", required=True, help="Subtrees pickle")
    parser.add_argument("-o", "--outpath", required=True, help="Output directory path")
    parser.add_argument("-l", "--level", required=True, help="Level of the subtree (l1/l2)")

    args = parser.parse_args()

    f = open(args.subtrees, "rb")
    st = pickle.load(f)

    df = subtrees_stats(st, args.level)
    print(df)

    df.to_csv(os.path.join(args.outpath, "subtrees_metric.tsv"), index=False, sep="\t")

if __name__ == "__main__":
    sys.exit(main())            

