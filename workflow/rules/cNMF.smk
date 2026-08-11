rule cNMF_prepare:
    """
    Prepare step which normalizes the count matrix and prepares the factorization step (https://github.com/dylkot/cNMF).
    The script only requires a output directory not a filename so a fake output is used instead and the output dir is
    passed as a parameter.
    """
    input:
        matrix=f"results/{config["analysis_name"]}/preprocess/mtx/matrix.mtx",
        barcodes=f"results/{config["analysis_name"]}/preprocess/mtx/barcodes.tsv",
        genes=f"results/{config["analysis_name"]}/preprocess/mtx/genes.tsv",
    output:
        f"results/{config["analysis_name"]}/cNMF/cNMF_prepare/.done.txt",
    log:
        f"results/{config["analysis_name"]}/cNMF/cNMF_prepare.log",
    conda:
        "../envs/cNMF.yaml"
    threads: config["cNMF"]["cNMF_prepare"]["cores"]
    resources:
        mem_mb=config["cNMF"]["cNMF_prepare"]["mem_mb"],
        time=config["cNMF"]["cNMF_prepare"]["time"],
        queue=config["queues"]["cpu"],
    params:
        out_dir=f"results/{config["analysis_name"]}/cNMF/cNMF_prepare/",
    shell:
        "cnmf prepare \
            --output-dir {params.out_dir} \
            --name " + config["analysis_name"] + " \
            -c {input.matrix} \
            --max-nmf-iter " + str(config["cNMF"]["cNMF_prepare"]["max_nmf_iter"]) + " \
            -k " + str(list(range(config["cNMF"]["cNMF_prepare"]["k_min"], config["cNMF"]["cNMF_prepare"]["k_max"], config["cNMF"]["cNMF_prepare"]["k_step"]))) + " \
            --n-iter " + str(config["cNMF"]["cNMF_prepare"]["n_iter"])
