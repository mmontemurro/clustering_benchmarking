import os
import natsort as ns
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "From intersectBed to ginkgo beds")

parser.add_argument("-i", "-input", required=True, metavar="input/dir/path/cell.bed", action="store", help="Input file path", type=str)
parser.add_argument("-o", "--output", required=True, metavar="input/dir/path/cell.bed", action="store", help="Output file path", type=str)

args = parser.parse_args()

beds = args.beds_dir
outdir = args.out_dir

wm = lambda x: round(np.average(x, weights=df.loc[x.index, "len"]))

df = pd.read_csv(args.input, header=None, sep ="\t")
df[7][df[7]== '.'] = 2
df[5][df[5] == -1] = df[1]
df[6][df[6] == -1] = df[2]
df['len'] = df[6] - df[5]
df[7] = df[7].astype(int)
res = df.groupby([0,1,2]).agg(cn_weighted_mean=(7, wm))
res = res.reset_index()
#natural sort by chromosome (chr1, chr2, ..., chr10, chr11, ..., chr20)
res[0] = pd.Categorical(res[0], ordered=True, categories= ns.natsorted(res[0].unique()))
res = res.sort_values([0,1,2])
res.to_csv(args.output, sep = "\t", header=None, index=False)
    
