#!/usr/bin/env python3
#
# Ed Mountjoy
#
# Script adapted from: https://github.com/slowkow/snakefiles/blob/master/bsub.py
#
# bsub.py
#
# This script checks a Snakemake job"s properties (threads, resources) and chooses
# an appropriate LSF queue that meets the requirements. It also automatically
# chooses the queue that is least busy unless you already specified a queue.
#
# Usage
# -----
#
# Add "threads" and "resources" to your resource-intensive rules:
#
#     rule my_rule:
#         input: ...
#         output ...
#         threads: 4
#         resources:
#             mem=8000,                    # megabytes
#             runtime=35,                  # minutes
#             queue="my_favorite_queue"   # queue name
#
# Invoke snakemake with the path to bsub.py:
#
#     snakemake --jobs 999 --cluster "path/to/bsub.py -o bsub.stdout"
#
# Consider adding bsub.py to a folder in your $PATH, so you can do:
#
#     snakemake --jobs 999 --cluster "bsub.py -o bsub.stdout"

import os
import sys
import json
import argparse
import time

from subprocess import check_output

from snakemake.utils import read_job_properties

def main():

    # Parse command line
    parser = argparse.ArgumentParser()
    parser.add_argument("jobscript")
    args = parser.parse_args()

    # Parse the job properties
    job_properties = read_job_properties(args.jobscript)

    # By default, we use 1 thread.
    threads = job_properties.get("threads", 1)

    # Get defualt mem, runtimes and output files from cluster.json
    mem = int(job_properties["cluster"]["mem"])
    runtime = int(job_properties["cluster"]["runtime"])
    stdout = job_properties["cluster"]["output"]
    stderr = job_properties["cluster"]["error"]
    jobname = job_properties["cluster"]["name"]

    # If the rule has specified resources, replace with those
    mem = int(job_properties["resources"].get("mem", mem))
    runtime = int(job_properties["resources"].get("runtime", runtime))

    # Make log file directories
    os.makedirs(os.path.dirname(stdout), exist_ok=True)
    os.makedirs(os.path.dirname(stderr), exist_ok=True)

    # Let the user specify the queue.
    queue = job_properties["resources"].get("queue", None)

    # Otherwise, choose an appropriate queue based on required resources.
    if not queue:
        queue = get_queue(threads, mem, runtime)

    # If we fail to find a queue, exit with an error.
    if not queue:
        msg = "No valid queue! job_properties:\n"
        js = json.dumps(job_properties, indent=4, sort_keys=True)
        sys.stderr.write(msg + js)
        sys.exit(1)

    # Submit the job to the queue.
    run_bsub(queue, threads, mem, runtime, args.jobscript, jobname, stdout, stderr)
    time.sleep(1)

def run_bsub(queue, threads, mem, runtime, script, jobname, stdout, stderr):
    cmd = "bsub -J {j} -q {q} -n {t}".format(j=jobname, q=queue, t=threads)
    if mem:
        cmd += ' -R "select[mem>{m}] rusage[mem={m}] span[hosts=1]" -M{m}'.format(m=mem) # "resources" : "\"select[mem>2000] rusage[mem=2000] span[hosts=1]\"",
    if runtime:
        cmd += " -W {}".format(runtime)
    if stdout:
        cmd += " -o {}".format(stdout)
    if stderr:
        cmd += " -e {}".format(stderr)
    cmd += " {s}".format(s=script)
    print(cmd)
    return os.system(cmd)

def get_queue(threads, mem, runtime):
    # All the Sanger farm queues.
    queues = ["small", "normal", "long", "basement", "hugemem", "teramem"]
    # Find valid queues for this job"s requirements.
    retval = []
    # The other queues are all ok if we leave runtime=0.
    if threads == 24 and mem <= 256000 and runtime <= 30:
        retval.append("small")
    if threads <= 24 and mem <= 256000 and runtime <= 60 * 12:
        retval.append("normal")
    if threads <= 24 and mem <= 256000 and runtime <= 60 * 24 * 2:
        retval.append("long")
    if threads <= 24 and mem <= 256000 and runtime <= 60 * 24 * 30:
        retval.append("basement")
    if threads <= 24 and 196000 < mem < 727500 and runtime <= 60 * 24 * 15:
        retval.append("hugemem")
    if threads <= 24 and 727500 < mem < 2.9e6 and runtime <= 60 * 24 * 15:
        retval.append("teramem")
    # Make sure we have at least one valid queue.
    if not len(retval):
        return None

    # # Get the number of currently running jobs on each queue.
    # lines = check_output("bqueues").split(b"\n")[1:-1]
    # lines = [line.decode("utf-8").split() for line in lines]
    # njobs = {x[0]: int(x[7]) for x in lines}
    # # Among valid queues, choose the one with fewest running jobs.
    # return min(retval, key=lambda j: njobs[j])

    # Return the first of the suitable queues
    return retval[0]

if __name__ == "__main__":
    main()
