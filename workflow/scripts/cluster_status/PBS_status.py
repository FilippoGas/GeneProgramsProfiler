"""
Script: PBS_status.py
Author: Filippo Gastaldello
Date: 18/07/26
Description:
    Periodically checks job status to notice snakemake in the case a job fails.
    Success is only reported when the job is genuinely finished with exit status 0.
    A job that cannot be found by qstat is reported as "running" (never "success"),
    otherwise a live queue lookup race would cause a false success and a
    spurious MissingOutputException in snakemake.
"""
import subprocess
import sys

job_id = sys.argv[1]


def qstat(job_id, flag):
    try:
        return subprocess.check_output(
            f"qstat {flag} -f {job_id}", shell=True, stderr=subprocess.STDOUT
        ).decode()
    except subprocess.CalledProcessError:
        return None


def job_state_in(output, states):
    return any(f"job_state = {s}" in output for s in states)


# 1. Live queue: a job currently queued/running/held is clearly still running.
live = qstat(job_id, "")
if live is not None:
    if job_state_in(live, ["Q", "R", "H", "W", "M"]):
        print("running")
        sys.exit(0)

# 2. Finished-job history (qstat -x): decide success/failed from the real exit status.
history = qstat(job_id, "-x")
if history is not None:
    if job_state_in(history, ["F", "E"]):
        if "Exit_status = 0" in history or "exit_status = 0" in history:
            print("success")
        else:
            print("failed")
    else:
        print("running")
    sys.exit(0)

# 3. Job not found anywhere (still in a transition or record purged): be conservative.
print("running")