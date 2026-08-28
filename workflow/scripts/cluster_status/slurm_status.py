"""
Script: slurm_status.py
Author: Filippo Gastaldello
Date: 28/08/26
Description:
    Periodically checks job status to notice snakemake in the case a job fails
"""
import subprocess
import sys

job_id = sys.argv[1]

try:
    output = subprocess.check_output(
        f"sacct --job {job_id} --format=State,ExitCode --noheader -X",
        shell=True,
        stderr=subprocess.STDOUT,
    ).decode().strip()

    if not output:
        # Accounting may not have recorded the job yet; check the queue instead
        output = subprocess.check_output(
            f"squeue --job {job_id} --noheader",
            shell=True,
            stderr=subprocess.STDOUT,
        ).decode()

    if not output:
        # No accounting or queue record -> treat as running (as in PBS_status.py)
        print("running")
    elif "COMPLETED" in output and "0:0" in output:
        print("success")
    elif any(s in output for s in ["RUNNING", "PENDING", "CONFIGURING", "COMPLETING"]):
        print("running")
    else:
        print("failed")

except subprocess.CalledProcessError:
    print("success")
