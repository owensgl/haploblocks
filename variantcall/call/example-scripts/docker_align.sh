#!/bin/bash -ex

HERE=$(cd "$(dirname "$0")" && pwd)
REPODIR=$(cd "${HERE}"/../ && pwd)
DATADIR="${DATADIR:-/data/}"
OUTDIR=/data/aligned/

readonly DOCKERIMAGE=genomics/analytics

VOLUMEOPTS=(
    -v "${REPODIR}/":/mnt/genomics:ro
    -v "${DATADIR}/":/mnt/data:ro
    -v "${OUTDIR}/":/mnt/aligned
)

[[ -x "$REPODIR"/bin/align ]] || {
    echo "missing align binary. make sure you run $REPODIR/docker-build.sh" >&2
    exit 1
}

[[ -d "$OUTDIR" ]] || {
    mkdir -p "$OUTDIR"
}

set -o pipefail

tmplog=tmp-$RANDOM-$$
mkfifo "$tmplog"
trap 'rm -f "$tmplog"' EXIT
tee -a "docker_align.log" < "$tmplog" &
tailpid=$$
(
    rm "$tmplog"
    docker run --rm "${VOLUMEOPTS[@]}" "${DOCKERIMAGE}" /mnt/genomics/bin/align \
	   -cas /mnt/data/cas \
    	   -r file:/mnt/data/ref_genome/HanXRQr1.0-20151230.fa \
    	   -i file:/mnt/genomics/example-scripts/align.job \
    	   -x file:/mnt/data/ref_genome/PE-Marco.fa \
	   -sb \
           -m \
    	   -w /mnt/aligned \
    	   -o file:/mnt/aligned \
    	   -lossy -d 1 | tee -a docker_align.log
) 1>"$tmplog" 2>&1 
