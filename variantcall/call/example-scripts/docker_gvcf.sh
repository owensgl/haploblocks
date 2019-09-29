#!/bin/bash -ex

HERE=$(cd "$(dirname "$0")" && pwd)
REPODIR=$(cd "${HERE}"/../ && pwd)
REFERENCE="${DATADIR:-/data/ref_genome/}"
#INFILE="${INFILE:-/data/aligned/DBGBGBS_GB334.bam}"
INFILE="${INFILE:-/data/aligned/merged_ANN1002.bam}"
OUTDIR="${OUTDIR:-/data/gvcf/}"
BEDFILE="${BEDFILE:-HanXRQr1.0-20151230_allTEs_ubc.non-repetitive-regions.2017.sorted.bed}"


readonly DOCKERIMAGE="${DOCKERIMAGE:-broadinstitute/gatk:4.0.1.2}"

if [[ ! -f "$INFILE" || ! -f "$INFILE".bai ]]; then
    echo "missing bam file: $INFILE and $INFILE.bai" >&2
    exit 1
fi

VOLUMEOPTS=(
    -v "${REPODIR}/":/mnt/genomics:ro
    -v "${REFERENCE}/":/mnt/data:ro
    -v "${INFILE}":/mnt/aligned/"$(basename "$INFILE")":ro
    -v "${INFILE}.bai":/mnt/aligned/"$(basename "$INFILE")".bai:ro
    -v "${OUTDIR}/":/mnt/gvcf
)

[[ -x "$REPODIR"/bin/vc ]] || {
    echo "missing vc binary. make sure you run $REPODIR/docker-build.sh" >&2
    exit 1
}

[[ -d "$OUTDIR" ]] || {
    mkdir -p "$OUTDIR"
}

set -o pipefail

tmplog=tmp-$RANDOM-$$
mkfifo "$tmplog"
trap 'rm -f "$tmplog"' EXIT
tee -a "docker_gvcf.log" < "$tmplog" &
tailpid=$$
(
    rm "$tmplog"
    docker run --rm "${VOLUMEOPTS[@]}" "${DOCKERIMAGE}" /mnt/genomics/bin/vc -n 24 -g \
	   -r file:/mnt/data/HanXRQr1.0-20151230.fa \
	   -b file:/mnt/data/"${BEDFILE}" \
           -i file:/mnt/aligned/"$(basename "$INFILE")" \
	   -o file:/mnt/gvcf/ \
	   -w /mnt/gvcf/ \
	   -minbp 0 \
	   -nsegments 32 \
	   -gatk4
) 1>"$tmplog" 2>&1 
