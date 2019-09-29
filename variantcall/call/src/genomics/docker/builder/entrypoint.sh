#!/bin/bash -ex

# Volumes:
# /build/src should contain the sources
# /build/bin will contain the compiled output
# /go/       will contain go tools

function usage ()
{
    echo "Usage $(basename "$0") [--build] [--dep]

    --build      Build the tools

    --dep        Download/Prepare third party deps
                 (do this before the first build)
"
}

ACTIONS=()

while [[ $# -gt 0 ]]; do
    arg="$1"
    shift
    case "$arg" in
	--build)
	    ACTIONS+=( build )
	    ;;
	--dep)
	    ACTIONS+=( dep )
	    ;;
	--help|-h)
	    usage
	    exit 0;
	    ;;
	*)
	    echo "unexpected arg: $arg" >&2
	    exit 1
	    ;;
    esac
done

[[ 0 -eq "${#ACTIONS[@]}" ]] && {
    echo "missing flag" >&2
    exit 1
}

if [[ -z "$GOPATH" ]]; then
    export GOPATH="/go:/usr/local/go:/build"
fi

export PATH="/usr/lib/go-1.9/bin:/usr/local/go/bin/:/build/bin:/go/bin:$PATH"

for ACTION in "${ACTIONS[@]}"; do
    case "$ACTION" in
	build)
	    ( cd /build/src/genomics
	      make GOBIN=/build/bin install
	    )
	    ;;
	dep)
	    ( cd /build/src/genomics
	      make dep
	    )
	    ;;
    esac
done

