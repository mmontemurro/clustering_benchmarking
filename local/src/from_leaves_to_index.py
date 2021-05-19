#!/usr/bin/env python3

import sys, os 
import argparse
import numpy as np 

def main():
    parser = argparse.ArgumentParser(description="Extracts tree leaf indices")

    parser.add_argument("-i", "--input", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    nodes = np.load(args.input, allow_pickle=True)
    leaf_index = []
    for n in nodes:
        leaf_index.append(n.id)

    np.save(args.output, leaf_index)

if __name__ == "__main__":
    sys.exit(main())