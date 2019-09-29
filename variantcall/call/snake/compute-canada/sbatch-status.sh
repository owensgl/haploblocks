#!/bin/bash

#
# Prints the status of a job launched on slurm/sbatch
#

JOBID="$1"

exec 3>>./sbatch-status.log

function log_result ()
{
    echo "[$(date -Iseconds)] $JOBID $@" >&3
}

# http://snakemake.readthedocs.io/en/latest/executable.html?highlight=cluster-status
# return one of {'success','failed','running'}

set -eo pipefail

# BF  BOOT_FAIL       Job  terminated  due to launch failure, typically due to a hardware failure (e.g. unable to boot the node or block and
#                     the job can not be requeued).
# CA  CANCELLED       Job was explicitly cancelled by the user or system administrator.  The job may or may not have been initiated.
# CD  COMPLETED       Job has terminated all processes on all nodes with an exit code of zero.
# CF  CONFIGURING     Job has been allocated resources, but are waiting for them to become ready for use (e.g. booting).
# CG  COMPLETING      Job is in the process of completing. Some processes on some nodes may still be active.
# DL  DEADLINE        Job missed its deadline.
# F   FAILED          Job terminated with non-zero exit code or other failure condition.
# NF  NODE_FAIL       Job terminated due to failure of one or more allocated nodes.
# PD  PENDING         Job is awaiting resource allocation. Note for a job to be selected in this state it must have  "EligibleTime"  in  the
#                     requested  time  interval or different from "Unknown". The "EligibleTime" is displayed by the "scontrol show job" com‐
#                     mand.  For example jobs submitted with the "--hold" option will have "EligibleTime=Unknown" as they are pending indef‐
#                     initely.
# PR  PREEMPTED       Job terminated due to preemption.
# R   RUNNING         Job currently has an allocation.
# RS  RESIZING        Job is about to change size.
# S   SUSPENDED       Job has an allocation, but execution has been suspended.
# TO  TIMEOUT         Job terminated upon reaching its time limit.

# UNDOCUMENTED
# OUT_OF_ME+          Job ran over its memory allocation

function parse_state ()
{
    local state="$1"
    case "${state^^}" in
	BOOT_FAIL*|CANCELLED*|DEADLINE*|FAILED*|NODE_FAIL*|PREEMPTED*|TIMEOUT*|OUT_OF_ME*)
	    log_result "$state -> failed"
	    echo "failed"
	    ;;
	COMPLETED*)
	    echo "success"
	    ;;
	CONFIGURING*|COMPLETING*|PENDING*|SUSPENDED*|RUNNING*|RESIZING*)
	    log_result "$state -> running"
	    echo "running"
	    ;;
	*)
	    log_result "$state (unknown) -> failed"
	    echo "Unknown state: $state" >&2
	    echo "failed"
	    exit 1
    esac
}

try_seconds=3000
try_sleep=5

# round up
tries=$(( (try_seconds+try_sleep-1) / try_sleep))

while [[ $tries -gt 0 ]]; do
    # FIXME compute-canada/sbatch-status.sh: Cannot send after transport endpoint shutdown
    while read JOBi STATEi; do
	if [[ "$JOBi" == "$JOBID" ]]; then
	    parse_state "$STATEi"
	    exit 0
	fi
    done < <(sacct -n -j "$JOBID" --format=JobID,State 2>/dev/null | head -n 1)

    # no such job -- could be freshly started
    tries=$((tries - 1))
    if [[ tries -gt 0 ]]; then
	log_result "no such job. (${tries} left)"
	sleep ${try_sleep}
    fi
done
log_result "no such job -> failed"

# FIXME just in case we timed out and it's still running, we cancel
#       the job.  snakemake's behavior is to resubmit failed jobs, and
#       things get worse if there are two copies of the same job
#       running. (esp. around logging and cleaning up outputs on
#       failed results)
scancel "$jobid" || :
echo "failed"
exit 0
