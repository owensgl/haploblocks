#!/bin/bash

set -Eexo pipefail

#
# Builds and places the go programs under ./bin, using a docker container
#

BUILDERIMAGE=rieseberglab/go-builder:3
OUTPUTIMAGE=genomicsbuild:latest
KEEP="${KEEP:-y}"

HERE=$(cd "$(dirname "$0")" && pwd)
REPODIR="$(dirname "$(dirname "$HERE")")"

function _has_docker_image ()
{
    # check if image exists
    local img="$1"
    docker image inspect "$img" >/dev/null 2>&1 || {
	return $?
    }
    return 0;
}

function write_dockerfile ()
{
    echo "
FROM $BUILDERIMAGE

# get the source in the container
COPY . \$GOPATH/src/genomics/
WORKDIR \$GOPATH/src/genomics/

RUN make dep
# to remove symbols: -ldflags=”-w -s”
# FIXME fully static: --ldflags '-linkmode external -extldflags "-static"'
RUN make XBUILDFLAGS=-a CGO_ENABLED=0 GOOS=linux GOARCH=amd64 build
RUN make GOBIN=/build/bin install
"
}

function _run_container ()
{
    local img=$1
    shift
    (
	tartmp=""
	trap 'if [[ -d "${tartmp}" ]]; then rm --one-file-system -r -- "${tartmp}"; fi' EXIT;
	cd "$HERE"

	tartmp=$(mktemp --tmpdir -d docker-build.XXXXXXXX)
	write_dockerfile > "${tartmp}/Dockerfile"
	git rev-parse HEAD | tee "${tartmp}/BUILDNO"

	# send contents to Docker builder context
	cp -r . "${tartmp}"
	docker build -t "$OUTPUTIMAGE" "${tartmp}"

	:
	: docker image "$OUTPUTIMAGE" contains build output in /build/bin/
	:
    )
    : extract compiled files out of build archive
    docker run --rm "$OUTPUTIMAGE" tar c -C /build/bin -vf - . | tar -xf - -C "$REPODIR/bin"

    if [[ "$KEEP" == "n" ]]; then
	docker rmi "$OUTPUTIMAGE"
    fi
}

if ! _has_docker_image "$BUILDERIMAGE"; then
    cd docker/builder2 && docker build -t "$BUILDERIMAGE" .
fi

set -o pipefail
echo "Using current directory as sources and output (./bin  and ./vendor)..."
_run_container "$BUILDERIMAGE" "$@" 2>&1 | tee docker-build.log
