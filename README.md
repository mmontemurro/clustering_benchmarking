# Clustering benchmarking

This repository contains the code used to perform a first formal and systematic performance evaluation of nine well-known clustering methods on scCNA data, on simulated sequencing data representing known phylogenies of cancer cells.

It contains a compiled version of the simulator presented by Fan et al. [1] and freely available at https://bitbucket.org/xianfan/cnsc_simulator/src/master/. It also contains a customized version of TreeCluster, a tool for cluster extraction from phylogenetic trees, presented by Balaban et al.[2] and freely available at https://github.com/niemasd/TreeCluster.


## Requirements
* Python >= 3.6.7
* Perl >= 5.16.3
* R >= 3.6.2

## Python dependencies:
* argparse
* json
* pprint
* gzip
* time
* math
* glob
* natsort
* operator
* itertools
* pandas
* numpy
* subprocess
* multiprocessing
* joblib
* yellowbrick
* sklearn
* scipy
* umap
* hdbscan
* kneed
* matplotlib
* seaborn
* ete3
* biopython
* treeswift
* niemads
* anytree
* graphviz
 
## Reference genome
We have used the GrCh38 reference genome which can be found at following link:
  https://drive.google.com/file/d/1Q1e0QTltVQYgMK0POFZRPLeSwj0BsUwW/view?usp=sharing
  
## Steps to reproduce our analysis
  
	#download reference genome and unzip
  	wget --load-cookies /tmp/cookies.txt "https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://docs.google.com/uc?export=download&id=1Q1e0QTltVQYgMK0POFZRPLeSwj0BsUwW' -O- local/share/data/genome.fa.gz | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1\n/p')&id=1Q1e0QTltVQYgMK0POFZRPLeSwj0BsUwW" -O genome.fa.gz && rm -rf /tmp/cookies.txt 
  	unzip local/share/data/genome.fa.gz -d local/share/data

  	#generate synthetic data
  	snakemake -s local/share/snakerule/Snakefile all

  	#run clustering and evaluation pipeline on non-reduced data
  	snakemake -s local/share/snakerule/Snakefile_clusterings.smk all

  	#run clustering and evaluation pipeline on PCA-reduced data
  	snakemake -s local/share/snakerule/Snakefile_clusterings_pca.smk all

  	#run clustering and evaluation pipeline on UMAP-reduced data
  	snakemake -s local/share/snakerule/Snakefile_clusterings_umap.smk all

  	#run merge all
  	snakemake -s local/share/snakerule/Snakefile_merge.smk merge_all

---
## References
  1. Fan, Xian, et al. "Benchmarking tools for copy number aberration detection from single-cell DNA sequencing data." bioRxiv (2019): 696179.
  2. Balaban, Metin, et al. "TreeCluster: Clustering biological sequences using phylogenetic trees." PloS one 14.8 (2019): e0221068.
---
