#!/bin/bash

set -euo pipefail
HOURS=24
LONGNAMES=60

ARGS=()
while [[ "$#" -gt 0 ]]; do
    case "$1" in
	--hours|-h)
	    HOURS="$2"
	    shift 2
	    ;;
	--xlong)
	    LONGNAMES=180
	    shift
	    ;;
	*)
	    ARGS+=( "$1" )
	    shift
	    ;;
    esac
done

START=$(date -d "$HOURS hours ago" +"%Y-%m-%d"T%H:%M:%S)

sacct --format "jobid%-10,elapsed,state%15,jobname%${LONGNAMES},maxrss,reqtres%25,timelimit,nodelist" -S "$START" "${ARGS[@]:+${ARGS[@]}}"
