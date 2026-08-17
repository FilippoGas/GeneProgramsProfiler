# Create logdir for PBS inside the repo's results directory
PBS_LOGDIR="logs/PBS/"
if [ ! -e $PBS_LOGDIR ]; then
    mkdir -p $PBS_LOGDIR
fi

snakemake \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "qsub -e $PBS_LOGDIR -o $PBS_LOGDIR -l select=1:ncpus={threads}:mem={resources.mem_mb}mb -l walltime={resources.time} -q {resources.queue}" \
    --cluster-generic-status-cmd "python workflow/scripts/cluster_status/PBS_status.py" \
    --jobs 50 \
    --keep-going \
    --sdm conda \
    --rerun-incomplete \
    --latency-wait 120
