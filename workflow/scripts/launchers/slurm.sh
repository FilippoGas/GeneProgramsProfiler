# Create logdir for Slurm inside the repo's results directory
SLURM_LOGDIR="logs/SLURM/"
if [ ! -e $SLURM_LOGDIR ]; then
    mkdir -p $SLURM_LOGDIR
fi

snakemake \
    --executor cluster-generic \
    --cluster-generic-submit-cmd "sbatch --output=$SLURM_LOGDIR --error=$SLURM_LOGDIR --cpus-per-task={threads} --mem={resources.mem_mb}M --time={resources.time} --partition={resources.queue}" \
    --cluster-generic-status-cmd "python workflow/scripts/cluster_status/slurm_status.py" \
    --jobs 100 \
    --keep-going \
    --sdm conda \
    --rerun-incomplete \
    --latency-wait 600
