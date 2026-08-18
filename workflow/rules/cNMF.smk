import os

# Define the list of worker indices (e.g., ["0", "1", "2", ...]) based on config cores
WORKERS = [str(i) for i in range(config["cNMF"]["cNMF_factorize_worker"]["cores"])]


rule cNMF_prepare:
    """
    Prepare step which normalizes the count matrix and prepares the factorization step (https://github.com/dylkot/cNMF).

    """
    input:
        matrix=f"results/{config['analysis_name']}/preprocess/mtx/matrix.mtx",
        barcodes=f"results/{config['analysis_name']}/preprocess/mtx/barcodes.tsv",
        genes=f"results/{config['analysis_name']}/preprocess/mtx/genes.tsv",
    output:
        done=f"results/{config['analysis_name']}/cNMF/.done_prepare.txt",
    log:
        f"logs/{config['analysis_name']}/cNMF/cNMF_prepare.log",
    conda:
        "../envs/cNMF.yaml"
    threads: config["cNMF"]["cNMF_prepare"]["cores"]
    resources:
        mem_mb=config["cNMF"]["cNMF_prepare"]["mem_mb"],
        time=config["cNMF"]["cNMF_prepare"]["time"],
        queue=config["queues"]["cpu"],
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.done),
        analysis_name=config["analysis_name"],
        max_nmf_iter=config["cNMF"]["cNMF_prepare"]["max_nmf_iter"],
        n_iter=config["cNMF"]["cNMF_prepare"]["n_iter"],
        k_vals=" ".join(
            map(
                str,
                range(
                    config["cNMF"]["cNMF_prepare"]["k_min"],
                    config["cNMF"]["cNMF_prepare"]["k_max"],
                    config["cNMF"]["cNMF_prepare"]["k_step"],
                ),
            )
        ),
    message:
        "Normalize count matrix and prepare the factorization step"
    shell:
        """
        cnmf prepare \
            --output-dir {params.out_dir} \
            --name {params.analysis_name} \
            -c {input.matrix} \
            --max-nmf-iter {params.max_nmf_iter} \
            -k {params.k_vals} \
            --n-iter {params.n_iter}

        touch {output.done}
        """


rule cNMF_factorize_worker:
    """
    Performs the actual matrix factorization.
    Submits an independent 1-core job for a single worker.
    """
    input:
        f"results/{config['analysis_name']}/cNMF/.done_prepare.txt",
    output:
        done=f"results/{config['analysis_name']}/cNMF/.done_factorize_worker_{{worker}}.txt",
    log:
        f"logs/{config['analysis_name']}/cNMF/cNMF_factorize_worker_{{worker}}.log",
    conda:
        "../envs/cNMF.yaml"
    threads: 1
    resources:
        mem_mb=config["cNMF"]["cNMF_factorize_worker"]["mem_mb"],
        time=config["cNMF"]["cNMF_factorize_worker"]["time"],
        queue=config["queues"]["cpu"],
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.done),
        analysis_name=config["analysis_name"],
        total_workers=config["cNMF"]["cNMF_factorize_worker"]["cores"],
    message:
        "Performs the factorization of the expression matrix"
    shell:
        """
        cnmf factorize \
            --output-dir {params.out_dir} \
            --name {params.analysis_name} \
            --worker-index {wildcards.worker} \
            --total-workers {params.total_workers}

        touch {output.done}
        """


rule cNMF_factorize:
    """
    Tells Snakemake to wait until all individual worker files exist.
    """
    input:
        expand(
            f"results/{config['analysis_name']}/cNMF/.done_factorize_worker_{{worker}}.txt",
            worker=WORKERS,
        ),
    output:
        done=f"results/{config['analysis_name']}/cNMF/.done_factorize.txt",
    log:
        f"logs/{config["analysis_name"]}/cNMF/cNMF_factorize.log",
    localrule: True
    conda:
        "../envs/cNMF.yaml"
    shell:
        """
        touch {output.done}
        """


rule cNMF_combine:
    """
    Combine the results from different values of K.
    """
    input:
        f"results/{config['analysis_name']}/cNMF/.done_factorize.txt",
    output:
        done=f"results/{config['analysis_name']}/cNMF/.done_combine.txt",
    log:
        f"logs/{config['analysis_name']}/cNMF/cNMF_combine.log",
    conda:
        "../envs/cNMF.yaml"
    threads: config["cNMF"]["cNMF_combine"]["cores"]
    resources:
        mem_mb=config["cNMF"]["cNMF_combine"]["mem_mb"],
        time=config["cNMF"]["cNMF_combine"]["time"],
        queue=config["queues"]["cpu"],
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.done),
        analysis_name=config["analysis_name"],
    message:
        "Combine the results from different values of K."
    shell:
        """
        cnmf combine \
            --output-dir {params.out_dir} \
            --name {params.analysis_name}

        touch {output.done}
        """


rule cNMF_k_selection_plot:
    """
    Make a plot estimating the trade-off between higher values of K and stability and error.
    (Just for debuggig purpose, the actual k will be selected automatically)
    """
    input:
        done=f"results/{config['analysis_name']}/cNMF/.done_combine.txt",
    output:
        k_plot_data=f"results/{config["analysis_name"]}/cNMF/{config["analysis_name"]}/{config["analysis_name"]}.k_selection_stats.df.npz",
    log:
        f"logs/{config["analysis_name"]}/cNMF/cNMF_k_selection_plot.log",
    conda:
        "../envs/cNMF.yaml"
    threads: config["cNMF"]["cNMF_k_selection_plot"]["cores"]
    resources:
        mem_mb=config["cNMF"]["cNMF_k_selection_plot"]["mem_mb"],
        time=config["cNMF"]["cNMF_k_selection_plot"]["time"],
        queue=config["queues"]["cpu"],
    params:
        out_dir=lambda wildcards, input: os.path.dirname(input.done),
        analysis_name=config["analysis_name"],
    message:
        "Make a plot estimating the trade-off between higher values of K and stability and error."
    shell:
        """
        cnmf k_selection_plot \
            --output-dir {params.out_dir} \
            --name {params.analysis_name}
        touch {output.k_plot_data}
        """


rule extract_best_k:
    """
    On the basis of the statistics and plot generated by the previous rule,
    select the value of k (number of gene programs to consider) with the
    best tradeoff between stability and error
    #TODO is this the best way to choose k? check the script
    """
    input:
        k_selection_stats=f"results/{config["analysis_name"]}/cNMF/{config["analysis_name"]}/{config["analysis_name"]}.k_selection_stats.df.npz",
    output:
        best_k=f"results/{config["analysis_name"]}/cNMF/{config["analysis_name"]}/best_k.txt",
    log:
        f"logs/{config["analysis_name"]}/cNMF/extract_best_k.log",
    conda:
        "../envs/cNMF.yaml"
    threads: config["cNMF"]["extract_best_k"]["cores"]
    resources:
        mem_mb=config["cNMF"]["extract_best_k"]["mem_mb"],
        time=config["cNMF"]["extract_best_k"]["time"],
        queue=config["queues"]["cpu"],
    message:
        "On the basis of the statistics and plot generated by the previous rule select the value of k."
    script:
        "../scripts/cNMF/extract_best_k.R"


rule cNMF_consensus:
    """
    Generate program usage tables for the selected k
    """
    input:
        best_k=f"results/{config["analysis_name"]}/cNMF/{config["analysis_name"]}/best_k.txt",
    output:
        done=f"results/{config['analysis_name']}/cNMF/.done_consensus.txt",
    log:
        f"logs/{config["analysis_name"]}/cNMF/cNMF_consensus.log",
    conda:
        "../envs/cNMF.yaml"
    threads: config["cNMF"]["cNMF_consensus"]["cores"]
    resources:
        mem_mb=config["cNMF"]["cNMF_consensus"]["mem_mb"],
        time=config["cNMF"]["cNMF_consensus"]["time"],
        queue=config["queues"]["cpu"],
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output.done),
        k=lambda wildcards, input: open(input.best_k, "r").read().strip(),
        ldt=config["cNMF"]["cNMF_consensus"]["local_density_threshold"],
        analysis_name=config["analysis_name"],
    message:
        "Generate program usage tables for the selected k."
    shell:
        """
        cnmf consensus \
            --output-dir {params.output_dir} \
            --name {params.analysis_name} \
            --components {params.k} \
            --local-density-threshold {params.ldt} \
            --show-clustering #TODO in README.md warn that adjusting ldt may be needed
        touch {output.done}
        """


rule cNMF_rename_programs:
    """
    Run Over Representation Analysis on gene program's markers to label them
    """
    input:
        best_k=f"results/{config["analysis_name"]}/cNMF/{config["analysis_name"]}/best_k.txt",
        done=f"results/{config['analysis_name']}/cNMF/.done_consensus.txt",
        cytopus_list=f"results/{config["analysis_name"]}/spectra/cytopus_Gene_sets.json",
        metadata=f"results/{config["analysis_name"]}/DE_analysis/metadata.csv",
    output:
        labeled_cell_scores=f"results/{config["analysis_name"]}/cNMF/{config["analysis_name"]}/{config["analysis_name"]}.labeled_cell_score.csv",
    log:
        f"logs/{config["analysis_name"]}/cNMF/cNMF_rename_programs.log",
    conda:
        "../envs/cNMF.yaml"
    threads: config["cNMF"]["cNMF_rename_programs"]["cores"]
    resources:
        mem_mb=config["cNMF"]["cNMF_rename_programs"]["mem_mb"],
        time=config["cNMF"]["cNMF_rename_programs"]["time"],
        queue=config["queues"]["cpu"],
    params:
        output_dir=lambda wildcards, output: os.path.dirname(output.labeled_cell_scores),
        k=lambda wildcards, input: open(input.best_k, "r").read().strip(),
        ldt=config["cNMF"]["cNMF_consensus"]["local_density_threshold"],
        analysis_name=config["analysis_name"],
    message:
        "Run Over Representation Analysis on gene program's markers to label them."
    script:
        "../scripts/cNMF/rename_programs.R"
