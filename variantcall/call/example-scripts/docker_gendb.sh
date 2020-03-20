#!/bin/bash -ex

HERE=$(cd "$(dirname "$0")" && pwd)
REPODIR=$(cd "${HERE}"/../ && pwd)
REFERENCE="${REFERENCE:-/data/ref_genome/HanXRQr1.0-20151230.fa}"
OUTDIR="${OUTDIR:-/data/gendb/}"
GVCFDIR="${GVCFDIR:-/data/gvcf}"
MAPFILE="${MAPFILE:-${HERE}/gvcf.sample.map}"

readonly DOCKERIMAGE="${DOCKERIMAGE:-broadinstitute/gatk:4.beta.6}"

#line number selection in reference file
CHRLINE=1
INTERVAL=$(grep ">" "$REFERENCE" | cut -f 1 -d " " | sed s/\>//g | sed -n "${CHRLINE}"p)

VOLUMEOPTS=(
    -v "${REPODIR}/":/mnt/genomics:ro
    -v "${MAPFILE}":/mnt/mapfile:ro
    -v "${OUTDIR}/":/mnt/gendb
    -v "${GVCFDIR}/":/mnt/gvcf:ro
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
    echo '
while read sample file rest; do
[[ -e "$file" ]] || { echo "sample $sample ($file) cannot be found" >&2; exit 1; }
  echo "SAMPLE $sample SIZE $(du -sh "$file")"
done < /mnt/mapfile
' | docker run --rm -i "${VOLUMEOPTS[@]}" "${DOCKERIMAGE}" /bin/bash -x

    docker run --rm "${VOLUMEOPTS[@]}" "${DOCKERIMAGE}" /gatk/gatk-launch GenomicsDBImport \
	   --genomicsDBWorkspace /mnt/gendb/"db-${INTERVAL}" \
	   -L "$INTERVAL" \
	   --batchSize 50 \
	   --validateSampleNameMap true\
           --sampleNameMap /mnt/mapfile

) 1>"$tmplog" 2>&1 
