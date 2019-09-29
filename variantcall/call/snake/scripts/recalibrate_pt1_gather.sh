#!/bin/bash

#
# Combines the partial models produced by VariantRecalibrator into a single
# model, with GatherTranches.
#
# https://software.broadinstitute.org/gatk/documentation/article.php?id=2805
#
# This needs to be run from a GATK singularity container.
#
# It produces two tranche files file usable as the input to recalibrate_pt2.sh
#
# usage:
#
#  recalibrate_pt1_gather.sh [--tmpdir TMPDIR] OUTSNP OUTINDEL INMODELTGZ+
#
#  If MODELTGZ starts with "@", then the parameter is assumed to be a
#  filename containing a list of input file names.
set -eo pipefail

CFGTMP=/tmp
INPUTS=()
OUTPUTS=()

MAX_MEM_MB=$((16*1024))

# shut up perl.
#
# the LANG environment variable is propagated from the login node into
# the compute job. but a different set of locales are installed.
export LANG=C

while [[ $# -gt 0 ]]; do
    arg="$1"
    shift
    case "$arg" in
	--tmpdir)
	    CFGTMP="$1"
	    shift
	    ;;
	--max-mem-mb)
	    MAX_MEM_MB="$1"
	    shift
	    ;;
	@*)
	    # read filenames from file
	    infile="${arg##@}"
	    while read inputfile rest; do
		INPUTS+=( "$inputfile" )
	    done < "$infile"
	    ;;
	*)
	    if [[ "${#OUTPUTS[@]}" -lt 2 ]]; then
		OUTPUTS+=( "$arg" )
	    else
		INPUTS+=( "$arg" )
	    fi
	    ;;
    esac
done

# We're running two of the same tool in parallel.
# And we keep some extra.
MAX_MEM_MB=$((MAX_MEM_MB / 2 - 2048))

if [[ ${MAX_MEM_MB} -lt 1024 ]]; then
    echo "Provide more ram to --max-mem-mb <AMOUNT_MB>" >&2
    exit 1
fi

JAVA_OPTIONS="-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx${MAX_MEM_MB}m"

set -x

OUTPUT_SNP="${OUTPUTS[0]}"
OUTPUT_INDEL="${OUTPUTS[1]}"

if [[ -z "${OUTPUT_SNP}" || -z "${OUTPUT_INDEL}" ]]; then
    echo "Missing output file(s)" >&2
    exit 1
fi

if [[ "${#INPUTS[@]}" -lt 1 ]]; then
    echo "Missing input file(s)" >&2
    exit 1
fi

XTMP="${SLURM_TMPDIR:-$CFGTMP}" # use local disk if on cluster
mkdir -p "$XTMP" "$XTMP"/gatk
workdir=$(mktemp -d -p $XTMP tmp.recalibrate_pt1_gather.XXXXXX);
trap 'rm --one-file-system -r ${workdir}' EXIT;

count=0
argsfile_snps="${workdir}/arguments_snps.txt"
argsfile_indels="${workdir}/arguments_indels.txt"

for modeltgz in "${INPUTS[@]}"; do
    partdir="${workdir}/model_part_${count}"
    echo "--input ${partdir}/snp.tranches"   >> "${argsfile_snps}"
    echo "--input ${partdir}/indel.tranches" >> "${argsfile_indels}"
    echo "mkdir -p ${partdir} && tar -C ${partdir} -xzf ${modeltgz}"
    let count=count+1
done | "${CONDA_PREFIX}/bin/parallel" -j2 --line-buffer

snp_tranches="snp.tranches"
indel_tranches="indel.tranches"

outputs=(
    "${snp_tranches}"
    "${indel_tranches}"
)

function gather_tranches ()
{
    local logprefix="$1"
    local LINE
    shift

    set -o pipefail
    (
	set -e
	HOME="${workdir}" /gatk/gatk GatherTranches \
	    --java-options "$JAVA_OPTIONS" \
	    --TMP_DIR "${XTMP}/gatk" "$@"
	echo "GatherTranches Done"
    ) |& ( while read LINE; do printf "%s%s\n" "$logprefix" "$LINE"; done )
}

export JAVA_OPTIONS
export workdir
export XTMP
export -f gather_tranches

(
    # STEP 1 - gather tranches for snps
    echo gather_tranches SNP_ --arguments_file "${argsfile_snps}" \
	 --mode SNP \
	 --truth-sensitivity-tranche 100.0 \
	 --truth-sensitivity-tranche  99.0 \
	 --truth-sensitivity-tranche  90.0 \
	 --truth-sensitivity-tranche  70.0 \
	 --truth-sensitivity-tranche  50.0 \
	 --output        "$workdir/${snp_tranches}"

    # STEP 2 - build model for indels
    echo gather_tranches IND_ --arguments_file "${argsfile_indels}" \
	 --mode INDEL \
	 --truth-sensitivity-tranche 100.0 \
	 --truth-sensitivity-tranche  99.0 \
	 --truth-sensitivity-tranche  90.0 \
	 --truth-sensitivity-tranche  70.0 \
	 --truth-sensitivity-tranche  50.0 \
	 --output        "$workdir/${indel_tranches}"

) | "${CONDA_PREFIX}/bin/parallel" -j2 --line-buffer

mv -vf "${workdir}/${snp_tranches}"   "${OUTPUT_SNP}"
mv -vf "${workdir}/${indel_tranches}" "${OUTPUT_INDEL}"
