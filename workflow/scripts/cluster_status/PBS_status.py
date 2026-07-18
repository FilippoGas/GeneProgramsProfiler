"""
Script: PBS_status.py
Author: Filippo Gastaldello
Date: 18/07/26
Description:
    Periodically checks job status to notice snakemake in the case a job fails
"""
import subprocess
import sys

job_id = sys.argv[1]

try:
    output = subprocess.check_output(f"qstat -x -f {job_id}", shell=True, stderr=subprocess.STDOUT).decode()
    
    if "job_state = Q" in output or "job_state = R" in output or "job_state = H" in output or "job_state = W" in output or "job_state = M" in output:
        print("running")
    elif "job_state = F" in output or "job_state = E" in output:
        if "Exit_status = 0" in output:
            print("success")
        else:
            print("failed")
    else:
        print("running")

except subprocess.CalledProcessError:
    print("success")