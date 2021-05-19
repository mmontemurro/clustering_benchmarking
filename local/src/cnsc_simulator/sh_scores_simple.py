import os, sys
import argparse
import pandas as pd
import phylics 

def main():
    parser = argparse.ArgumentParser(description="Compute SHscore on set of samples")
    parser.add_argument("-i", "--input", required=True, help="CNV matrices", nargs="+")
    parser.add_argument("-o", "--output", required=True, help="Ouput file")
    parser.add_argument("-c", "--config", required=True, help="Configuration file")
    parser.add_argument("-s", "--supersample", help="Supersample name")

    args = parser.parse_args()
    sample = args.supersample
    cnvs = []
    for i, cnv_path in enumerate(args.input):
        cnvs_ = phylics.Sample.from_file(cnv_path, sample_name=sample +"_"+ str(i+1))
        cnvs.append(cnvs_)
    ss = phylics.MultiSample.from_list(cnvs)
    s = ss.SHscore(n_jobs=10)

    c = pd.read_csv(args.config, sep="\t")
    df = pd.DataFrame({"supersample":args.supersample, "theta":str(c["exp_theta"].values[0]), "1/p":str(1/c["amp_p"].values[0]), "sh_score":s['score'].values[0]}, index=[0])
    df.to_csv(args.output, sep="\t", index=False)


if __name__ == "__main__":
    sys.exit(main())