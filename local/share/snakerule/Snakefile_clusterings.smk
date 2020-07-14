include: "../snakemake/conf.sk"

def get_stats(wildcards):
    import glob, os
    stats = glob.glob(os.path.join(DATASET_DIR, "{wildcards.exp}/clusterings/*/", "*_performance.csv"))
    return stats

def get_ext_val_done(wildcards):
    import glob, os
    pl = glob.glob(os.path.join(DATASET_DIR, "snake_workflow", "{wildcards.exp}-{wildcards.tool}-{wildcards.sample}--external-validation.done"))
    return pl

_data_folder_ = DATASET_DIR+"/{exp}/data"
_out_folder_ = DATASET_DIR+"/{exp}/clusterings/"
_plh_dir_ = DATASET_DIR+"/snake_workflow"


rule run_tools_nopca:
    output:
        touch(_plh_dir_+'/clustering.done')
    input:
        expand(_out_folder_ + '/affinity/{sample}_affinity_clusters.csv', sample=SAMPLES,
                exp=EXP),
        expand(_out_folder_ + '/agglomerative_{linkage}/{sample}_agglomerative_{linkage}_clusters.csv', sample=SAMPLES, linkage=LINKAGE,
                exp=EXP),
        expand(_out_folder_ + '/birch/{sample}_birch_clusters.csv', sample=SAMPLES,
                exp=EXP),
        expand(_out_folder_ + '/kmeans/{sample}_kmeans_clusters.csv', sample=SAMPLES,
                exp=EXP),
        expand(_out_folder_ + '/dbscan/{sample}_dbscan_clusters.csv', sample=SAMPLES,
                exp=EXP),
        expand(_out_folder_ + '/hdbscan/{sample}_hdbscan_clusters.csv', sample=SAMPLES,
                exp=EXP)


# ---------------------------------
# -------- Run Clusterings --------
# ---------------------------------

rule run_affinity_nopca:
    output:
        out_cluster = _out_folder_ + '/affinity/{sample}_affinity_clusters.csv'
    input:
        sim_file = _data_folder_ + '/{sample}/gt_merged/SegCopy'
    params:
        prgm = TOOLS['affinity'],
        out_folder = _out_folder_ + '/affinity', out_prefix=_out_folder_ + '/affinity/{sample}'
    shell:
        '''
            if [ ! -d {params.out_folder} ]; then
                mkdir -p {params.out_folder}
            fi

            python3 {params.prgm} -f {input.sim_file} -o {params.out_prefix} --preproc
        '''

rule run_agglomerative_nopca:
    output:
        out_cluster = _out_folder_ + '/agglomerative_{linkage}/{sample}_agglomerative_{linkage}_clusters.csv'
    input:
        sim_file = _data_folder_ + '/{sample}/gt_merged/SegCopy'
    params:
        prgm = TOOLS['agglomerative'],
        out_folder = _out_folder_ + '/agglomerative_{linkage}', out_prefix=_out_folder_ + '/agglomerative_{linkage}/{sample}',
        n_clusters = '/{sample}/gt_merged/tree_clusters_stats.csv'
    shell:
        '''
            if [ ! -d {params.out_folder} ]; then
                mkdir -p {params.out_folder}
            fi

            python3 {params.prgm} -f {input.sim_file} -o {params.out_prefix} -l {wildcards.linkage} --preproc
        '''

rule run_birch_nopca:
    output:
        out_cluster = _out_folder_ + '/birch/{sample}_birch_clusters.csv'
    input:
        sim_file = _data_folder_ + '/{sample}/gt_merged/SegCopy'
    params:
        prgm = TOOLS['birch'],
        out_folder = _out_folder_ + '/birch', out_prefix=_out_folder_ + '/birch/{sample}',
        n_clusters = '/{sample}/gt_merged/tree_clusters_stats.csv'
    shell:
        '''
            if [ ! -d {params.out_folder} ]; then
                mkdir -p {params.out_folder}
            fi

            python3 {params.prgm} -f {input.sim_file} -o {params.out_prefix} --preproc
        '''

rule run_kmeans_nopca:
    output:
        out_cluster = _out_folder_ + '/kmeans/{sample}_kmeans_clusters.csv'
    input:
        sim_file = _data_folder_ + '/{sample}/gt_merged/SegCopy'
    params:
        prgm = TOOLS['kmeans'],
        out_folder = _out_folder_ + '/kmeans', out_prefix=_out_folder_ + '/kmeans/{sample}',
        n_clusters = '/{sample}/gt_merged/tree_clusters_stats.csv'
    shell:
        '''
            if [ ! -d {params.out_folder} ]; then
                mkdir -p {params.out_folder}
            fi

            python3 {params.prgm} -f {input.sim_file} -o {params.out_prefix} --preproc
        '''

rule run_dbscan_nopca:
    output:
        out_cluster = _out_folder_ + '/dbscan/{sample}_dbscan_clusters.csv'
    input:
        sim_file = _data_folder_ + '/{sample}/gt_merged/SegCopy'
    params:
        prgm = TOOLS['dbscan'],
        out_folder = _out_folder_ + '/dbscan', out_prefix=_out_folder_ + '/dbscan/{sample}'
    shell:
        '''
            if [ ! -d {params.out_folder} ]; then
                mkdir -p {params.out_folder}
            fi

            python3 {params.prgm} -f {input.sim_file} -o {params.out_prefix} --preproc
        '''
rule run_hdbscan_nopca:
    output:
        out_cluster = _out_folder_ + '/hdbscan/{sample}_hdbscan_clusters.csv'
    input:
        sim_file = _data_folder_ + '/{sample}/gt_merged/SegCopy'
    params:
        prgm = TOOLS['hdbscan'],
        out_folder = _out_folder_ + '/hdbscan', out_prefix=_out_folder_ + '/hdbscan/{sample}'
    shell:
        '''
            if [ ! -d {params.out_folder} ]; then
                mkdir -p {params.out_folder}
            fi

            python3 {params.prgm} -f {input.sim_file} -o {params.out_prefix} --preproc
        '''

rule stability_score:
    input: prev_step = _plh_dir_+'/clustering.done',
        stats = get_stats,
        file = _data_folder_ + "/{sample}/gt_merged/SegCopy",      
    output:
        touch(_plh_dir_+"/{exp}-{tool}-{sample}-stability.done")
    params: 
        clusters = _out_folder_ + '/{tool}/{sample}_{tool}_clusters.csv',
        prgm = STABILITY,
        n_iter = 100 
    shell:
        """
            python3 {params.prgm} -f {input.file} -c {params.clusters} -s {input.stats} -m {wildcards.tool} --preproc --n_iter {params.n_iter} 
        """   

rule ext_validation:
    input:
        prev_step= _plh_dir_+"/{exp}-{tool}-{sample}-stability.done",
        true_labels = _data_folder_+"/{sample}/gt_merged/tree_clusters.csv"
    output:
        touch(_plh_dir_+"/{exp}-{tool}-{sample}-external-validation.done")
    params:
        prgm =  EXT_VALIDATION,
        pred_labels = _out_folder_+"/{tool}/{sample}_{tool}_clusters.csv",
        out_prefix = _out_folder_+"/{tool}/{sample}_{tool}"
    shell:
        """
            python {params.prgm} -t {input.true_labels} -p {params.pred_labels} -o {params.out_prefix}
        """

rule merge_metrics:
    input: prev_step=get_ext_val_done
    output:
            apn_stats = _out_folder_ + "/apn_scores.csv",
            time_stats = _out_folder_ + "/computation_times.csv",
            out_file = _out_folder_ + "/validation_indices_heatmap.png"
    params:
            data_folder = _out_folder_,
            out_folder = _out_folder_, _sim_=50, 
            clusterings = CLUSTERINGS_STR,
            prgm = MERGE_METRICS
    shell:
        """
            python3 {params.prgm} -d {params.data_folder} -s {params._sim_} -c {params.clusterings} -o {params.out_folder}
        """

rule apn_stats:
    input:
        file = _out_folder_ + "/apn_scores.csv"
    output:
        out_file = _out_folder_ + "/apn_stats.csv"
    params:
        prgm = APN_STATS,
        out_folder = _out_folder_
    shell:
        """
            python3 {params.prgm} -f {input.file} -o {params.out_folder}
        """


rule time_stats:
    input:
        file = _out_folder_ + "/computation_times.csv"
    output:
        out_file = _out_folder_ + "/comp_time_stats.csv"
    params:
        prgm = TIME_STATS,
        out_folder = _out_folder_
    shell:
        """
            python3 {params.prgm} -f {input.file} -o {params.out_folder}
        """


rule all:
    input: expand(_out_folder_ + "/validation_indices_heatmap.png", exp=EXP),
            expand(_out_folder_ + "/apn_stats.csv", exp=EXP),
            expand(_out_folder_ + "/comp_time_stats.csv", exp=EXP)
    output: touch(_plh_dir_+"/clustering_evaluation.done")