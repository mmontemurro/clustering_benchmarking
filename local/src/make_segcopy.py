import os
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser(description = "From intersectBed to ginkgo beds")

parser.add_argument("-b", "--beds", required=True, metavar="[file1.bed, file2.bed...]", action="store", help="Input bed files", nargs="+")
parser.add_argument("-o", "--output", required=True, metavar="output/dir/path/SegCopy", action="store", help="Output directory file", type=str)

args = parser.parse_args()

beds = args.beds
segcopy = pd.DataFrame(columns=["CHR", "START", "END"])
for bed in beds:
        df = pd.read_csv(bed, sep="\t", header=None)
        cellid = os.path.splitext(os.path.basename(bed))[0]
        if len(segcopy) == 0: #first cell
            segcopy["CHR"] = df[0]
            segcopy["START"] = df[1]
            segcopy["END"] = df[2]
        
        segcopy[cellid] = df[3]

segcopy.to_csv(args.output, sep="\t", index=False)
