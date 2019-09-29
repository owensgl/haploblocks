#!/bin/bash -ex

HERE=$(cd "$(dirname "$0")" && pwd)
REPODIR=$(cd "${HERE}"/../ && pwd)
REFERENCE="${REFERENCE:-/data/ref_genome/HanXRQr1.0-20151230.fa}"
INDIR="${INDIR:-/data/gendb/}"
OUTDIR="${OUTDIR:-/data/jointcalled/}"
MAPFILE="${MAPFILE:-${HERE}/gvcf.sample.map}"

readonly DOCKERIMAGE="${DOCKERIMAGE:-broadinstitute/gatk:4.beta.6}"

#line number selection in reference file
CHRLINE=1
INTERVAL=$(grep ">" "$REFERENCE" | cut -f 1 -d " " | sed s/\>//g | sed -n "${CHRLINE}"p)

VOLUMEOPTS=(
    -v "$(dirname "${REFERENCE}")":/mnt/reference:ro
    -v "${REPODIR}/":/mnt/genomics:ro
    -v "${INDIR}/":/mnt/gendb:ro
    -v "${OUTDIR}/":/mnt/output
)

[[ -d "$OUTDIR" ]] || {
    mkdir -p "$OUTDIR"
}

set -o pipefail

# Tool documentation
# https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_genomicsdb_GenomicsDBImport.php#--intervals

tmplog=tmp-$RANDOM-$$
mkfifo "$tmplog"
trap 'rm -f "$tmplog"' EXIT
tee -a "$(basename "$0" .sh).log" < "$tmplog" &
tailpid=$$
(
    rm "$tmplog"
    docker run --rm "${VOLUMEOPTS[@]}" "${DOCKERIMAGE}" /gatk/gatk-launch GenotypeGVCFs \
	   -R /mnt/reference/"$(basename "${REFERENCE}")" \
	   -V gendb:///mnt/gendb/"db-${INTERVAL}" \
	   -L "$INTERVAL" \
	   -O /mnt/output/contig-"${INTERVAL}".tmp.vcf.gz \
	   --secondsBetweenProgressUpdates 5 --verbosity DEBUG
) 1>"$tmplog" 2>&1 
