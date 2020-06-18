#!/usr/bin/envs python
import os, sys, glob
import argparse
import pandas as pd 

def main():
    parser = argparse.ArgumentParser(description="Merge clustering statistics.")
    parser.add_argument("-d", "--data_folder", help="Data folder path.", required=True)
    parser.add_argument("-s","--simulations", help="Simulation number.", required=True, type=int)
    parser.add_argument("-c", "--clusterings", nargs="+", help="Clustering algorithm list.", required=True)
    parser.add_argument("-o", "--out_folder", help="Output folder.", required=True)

    args = parser.parse_args()

    data_folder = args.data_folder
    simulations = args.simulations
    clusterings = args.clusterings
    out_folder = args.out_folder


    samples = []
    for i in range(1, simulations + 1):
        samples.append("Sample" + str(i))

    cophenetic_scores = pd.DataFrame(columns=clusterings, index=samples)
    computation_times = pd.DataFrame(columns=clusterings, index=samples)
    clusters_n = pd.DataFrame(columns=clusterings, index=samples)
    for c in clusterings:
        scores = pd.Series(index=samples)
        times = pd.Series(index=samples)
        clusters = pd.Series(index=samples)
        prefix = os.path.join(data_folder, c)
        for s in samples:
            filepath = os.path.join(prefix, s + "_" + c + "_performance.csv")
            stats = pd.read_csv(filepath, index_col=0, sep="\t", header=None).transpose()
            #to fix an error
            stats = stats.loc[:,~stats.columns.duplicated()]
            scores.loc[s] = stats["cophenetic_affinity_score"].values[0]
            times.loc[s] = stats["computation_time"].values[0]      
            clusters.loc[s] = stats["n_clusters_unfiltered"].values[0] 
        cophenetic_scores[c] = scores
        computation_times[c] = times
        clusters_n[c] = clusters

    cophenetic_scores.to_csv(os.path.join(out_folder, "cophenetic_affinity_scores.csv"), sep="\t")
    computation_times.to_csv(os.path.join(out_folder, "computation_times.csv"), sep="\t")
    clusters_n.to_csv(os.path.join(out_folder, "cluster_number.csv"), sep="\t")

if __name__ == "__main__":
    sys.exit(main())