import os


rule cNMF_prepare:
    """
    Prepare step which normalizes the count matrix and prepares the factorization step (https://github.com/dylkot/cNMF).
    The script only requires a output directory not a filename so a fake output is used instead and the output dir is
    passed as a parameter.
    """
    input:
        matrix=f"results/{config['analysis_name']}/preprocess/mtx/matrix.mtx",
        barcodes=f"results/{config['analysis_name']}/preprocess/mtx/barcodes.tsv",
        genes=f"results/{config['analysis_name']}/preprocess/mtx/genes.tsv",
    output:
        done=f"results/{config['analysis_name']}/cNMF/cNMF_prepare/.done.txt",
    log:
        f"results/{config['analysis_name']}/cNMF/cNMF_prepare.log",
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


rule cNMF_factorize:
    input:
        f"results/{config['analysis_name']}/cNMF/cNMF_prepare/.done.txt",
    output:
        done=f"results/{config['analysis_name']}/cNMF/cNMF_factorize/.done.txt",
    log:
        f"logs/{config['analysis_name']}/cNMF/cNMF_factorize.log",
    conda:
        "../envs/cNMF.yaml"
    threads: config["cNMF"]["cNMF_factorize"]["cores"]
    resources:
        mem_mb=config["cNMF"]["cNMF_factorize"]["mem_mb"],
        time=config["cNMF"]["cNMF_factorize"]["time"],
        queue=config["queues"]["cpu"],
    params:
        out_dir=lambda wildcards, output: os.path.dirname(output.done),
        analysis_name=config["analysis_name"],
    shell:
        """
        for core in $(seq 0 1 $(({threads} - 1))); do
            cnmf factorize \
                --output-dir {params.out_dir} \
                --name {params.analysis_name} \
                --worker-index $core \
                --total-workers {threads} &
        done
        wait
        touch {output.done}
        """
