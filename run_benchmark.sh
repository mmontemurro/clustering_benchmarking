#!/usr/bin/bash

#gdrive URL for reference genome = https://drive.google.com/file/d/1Q1e0QTltVQYgMK0POFZRPLeSwj0BsUwW/view?usp=sharing

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