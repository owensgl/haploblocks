#!/bin/bash

function usage () {
    echo \
"  Usage: $(basename $0) REF [--size N]

    Divides the genome in the reference file into contig windows bounded
    by size N (bp). Prints regions to stdout. Windows will not straddle
    two chromosomes.

    This is a wrapper script around window-creator.py

    Arguments:

       REF       reference sequence file pathname (.fa, .fa.gz, .bz2), or
                 the special value '-' to read from stdin.

       --size N  the max size of each contig windows, in basepairs.
"
}

set -euo pipefail

HERE=$(cd "$(dirname "$0")" && pwd)
CONTIG_CMD="${HERE}"/window-creator-all.py
REF=""

CONTIG_ARGS=( - )

while [[ "$#" -gt 0 ]]; do
    case "$1" in
	--size)
	    shift
	    CONTIG_ARGS+=(--size "$1")
	    ;;
	--help)
	    shift
	    usage
	    exit 0
	    ;;
	-)
	    shift
	    REF="$1"
	    ;;
	-*)
	    shift
	    usage >&2
	    exit 1
	    ;;
	*)
	    REF="$1"
	    ;;
    esac
    shift
done

if [[ -z "$REF" ]]; then
    echo "you must specify a reference" >&2
    exit 1
fi

case "$REF" in
    *.fai)
	CONTIG_ARGS+=(--faidx)
	;;
esac

case "$REF" in
    *.fa|*.fasta)
        egrep "^>" "$REF"
        ;;
    *.fa.gz)
        zcat "$REF" | egrep "^>"
        ;;
    *.bz2)
        bzcat "$REF" | egrep "^>"
        ;;
    *.fai)
	cat "$REF"
	;;
    -)
	egrep "^>"
	;;
    *)
        echo "Unsupported reference file format $REF" >&2
	exit 1
        ;;
esac |  "${CONTIG_CMD}" "${CONTIG_ARGS[@]}"

