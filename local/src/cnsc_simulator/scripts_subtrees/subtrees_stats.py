#!/usr/bin/env python3

import sys, os
from ..gen_tree import gen_tree
import pickle
import argparse
import pandas as pd

def subtrees_stats(trees_d, level="l1"):
    df = pd.DataFrame()
    if level == "l1":
        for k in trees_d.keys():
            root_id = trees_d[k][0].id
            n_leaves = len(trees_d[k][0].leaves)
            df.append({'subtree_id':k, 'root_id':root_id, 'n_leaves':n_leaves})
    else:
        for k in trees_d.keys():
            subtree_l1 = trees_d[k]["subtree_l1"]
            root_id = trees_d[k]["tree"][0].id
            n_leaves = len(trees_d[k]["tree"][0].leaves)
            df.append({'subtree_id':k, 'subtree_l1':subtree_l1, 'root_id':root_id, 'n_leaves':n_leaves})
            

