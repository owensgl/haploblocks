#!/bin/bash

#
# Converts a fasta file to a BED file containing the entirety
# of the genome defined in the fasta file.
#
# output is sorted in the order listed inside fasta
#

set -euo pipefail

FASTA="$1"
HERE=$(cd "$(dirname "$0")" && pwd)

while read CHR START END CHROMLEN REST; do
    if [[ "$START" -ne 1 ]]; then
	echo "window size too small. bump" >&2
	exit 1
    fi
    printf "${CHR}\t0\t${CHROMLEN}\n"
done < <(grep ">" "$FASTA" | "$HERE"/window-creator-all.py --size 999999999 -)

