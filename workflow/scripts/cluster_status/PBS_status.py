"""
Script: PBS_status.py
Author: Filippo Gastaldello
Description:
    Checks PBS job status to avoid Snakemake from hanging indefinetely if errors occur
"""

import subprocess
import sys

# Snakemake passes the cluster job ID as the first argument
job_id = sys.argv[1]

try:
    # Query PBS Pro for the full job state
    output = subprocess.check_output(f"qstat -f {job_id}", shell=True, stderr=subprocess.STDOUT).decode()
    
    # Check the standard PBS Pro state codes
    if "job_state = Q" in output or "job_state = R" in output or "job_state = H" in output or "job_state = W" in output:
        print("running")
    elif "job_state = F" in output or "job_state = E" in output:
        # F means finished. We must check the exit status to know if it succeeded or crashed.
        if "Exit_status = 0" in output:
            print("success")
        else:
            print("failed")
    else:
        print("failed")

except subprocess.CalledProcessError:
    # If qstat throws an error (e.g., the job was purged from the queue history entirely),
    # we report 'failed' so Snakemake knows to drop the job and stop hanging.
    print("failed")