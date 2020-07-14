include: "../snakemake/conf.sk"

rule merge_all:
    input: DATASET_DIR
    output: DATASET_DIR+"/cross_merge.csv"
    params: prgm = MERGE_ALL
    shell:
        """
            python {params.prgm} -i {input} -o {output}
        """