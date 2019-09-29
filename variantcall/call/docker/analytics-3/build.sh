#!/bin/bash

set -exo pipefail

HERE=$(readlink -f $(dirname "$0"))
FROMIMAGE="rieseberglab/analytics:2"
BUILDIMAGE="genomicsbuild:latest"  # where the output is
TOIMAGE="rieseberglab/analytics:3"

function _has_docker_image ()
{
    # check if image exists
    local img="$1"
    docker image inspect "$img" >/dev/null 2>&1 || {
	return $?
    }
    return 0;
}


if ! _has_docker_image "$FROMIMAGE"; then
    docker pull "$FROMIMAGE"
fi

(
    tartmp=$(mktemp --tmpdir -d docker-build.XXXXXXXX)
    trap 'if [[ -d "${tartmp}" ]]; then rm --one-file-system -r -- "${tartmp}"; fi' EXIT;

    cp Dockerfile "$tartmp"/
    mkdir -p "$tartmp"/fs "$tartmp"/fs/usr/local "$tartmp"/fs/etc

    (
	echo "$TOIMAGE"
	echo -n "git "
	git rev-parse HEAD
    ) | tee "$tartmp"/fs/etc/analytics-revision.txt
    # GNU tar
    docker run --rm "$FROMIMAGE" tar -C / -czvf - /home | tar -C "${tartmp}/fs" -xzf -
    # busybox tar
    docker run --rm "$BUILDIMAGE" tar c -C /build/ -zvf - ./bin | tar -C "${tartmp}/fs/usr/local" -xzf -
    docker build -t "$TOIMAGE"-tmp "${tartmp}"
    name=tmp-build-$RANDOM-$RANDOM
    (
	trap 'docker rm "$name" || :; docker rmi "$TOIMAGE"-tmp || :;' EXIT;

	docker run --name "$name" --entrypoint _INVALID_ "$TOIMAGE"-tmp || :
	docker export "$name" | docker import - "$TOIMAGE"
    )
)
