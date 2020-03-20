#!/bin/bash -ex

# Run this script to setup a local host that will serve the role of
# the genomics job controller.

ALLSTEPS=(
    install_debs
    install_conda
    check_operational
    permissions
    install_repo
    docker_images
    sample_genome
)


DEBDEPS=(
    git
    curl
    docker-ce

    # debugging
    nano
    vim
    emacs
)

GITHOST=github.com
GENOMICSREPO=git@"${GITHOST}":rieseberglab/snp-calling.git
RUSER=ubuntu
REPODIR="$HOME/src/genomics"
ANALYTICSIMAGE=genomics/analytics
GATK3IMAGE=broadinstitute/gatk3:3.8-0
#broadinstitute/gatk:4.beta.3 Created on 2017-07-26T15:55:37.228533128Z
GATK43IMAGE=broadinstitute/gatk:4.beta.3
GATK46IMAGE=broadinstitute/gatk:4.beta.6

#reference genome
SAMPLEGENOME="https://www.cs.ubc.ca/~jslegare/genomics-files/ref_genome.tgz"

#
MINICONDA="https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh"


POSARGS=()
STEPS=()

function usage ()
{
    echo "Usage: $(basename "$0") [--host HOST] [--step STEP]*

    STEP is one of: ${ALLSTEPS[@]}
"
}

function _has_step ()
{
    ( set +x
      local s="$1"
      for s2 in "${ALLSTEPS[@]}"; do
	  if [[ "$s" == "$s2" ]]; then
	      return 0
	  fi
      done
      return 1
    )
}

function _has_docker_image ()
{
    # check if image exists
    local img="$1"
    docker image inspect "$img" >/dev/null 2>&1 || {
	return $?
    }
    return 0;
}

function _docker_apt ()
{
    # make sure docker packages are coming from docker.com

    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo apt-key add -
    sudo add-apt-repository "deb [arch=amd64] https://download.docker.com/linux/ubuntu $(lsb_release -cs) stable"
    sudo apt-get update -y
    # make sure it's coming from docker.com
    apt-cache policy docker-ce | grep -q https://download.docker.com
    
}

function install_debs ()
{
    _docker_apt
    sudo apt-get update -y && \
	sudo apt-get install -y "${DEBDEPS[@]}"
}

function check_operational ()
{
    sudo systemctl status docker | grep -q "active (running)"
}

function install_conda ()
{
    # to get the environment file
    #install_repo

    if [[ ! -d ~/miniconda2/bin ]]; then
	curl "$MINICONDA" > conda.sh
	bash conda.sh -b
	rm conda.sh
    fi

    (
	export PATH=~/miniconda2/bin/:"$PATH"
	conda info --envs | grep -q "genomics-py36" 2>/dev/null >&2 || {
	    yes | conda env create -f "$REPODIR/snake/environment.yaml" -n genomics-py36
	}
    )
}

function permissions ()
{
    # configure services
    if egrep -q ^docker: /etc/group; then
	if ! ( egrep ^docker: /etc/group | grep -q "$RUSER" ); then
	    sudo usermod -aG docker "${RUSER}"
	    # might require logging out and in
	fi
    fi
}

function install_repo ()
{
    touch ~/.ssh/known_hosts
    chmod 644 ~/.ssh/known_hosts
    ssh-keygen -R "${GITHOST}"
    ssh-keyscan -H "${GITHOST}" >> ~/.ssh/known_hosts

    if [[ ! -d "$(dirname "$REPODIR")" ]]; then
	mkdir -p "$(dirname "$REPODIR")"
    fi

    # install genomics repo
    if [[ ! -d "$REPODIR" ]]; then
	rm --one-file-system -rf "$REPODIR.tmp"
	git clone "$GENOMICSREPO" "$REPODIR.tmp"
	mv "$REPODIR.tmp" "$REPODIR"
    fi
}

function docker_images ()
{
    local img
    
    # install docker images
    # FIXME setup an external docker registry instead.
    [[ -d "$REPODIR" ]] # prereq install_repo

    if ! _has_docker_image "$ANALYTICSIMAGE"; then
	docker import --message "imported by setup-controller.sh" \
	       "$REPODIR"/docker/analytics_container.tar.bz2 \
	       "$ANALYTICSIMAGE"
    fi

    for img in "${GATK3IMAGE}" "${GATK43IMAGE}" "${GATK46IMAGE}"; do
	_has_docker_image "$img" || docker pull "$img"
    done
}

function sample_genome ()
{
    [[ -d "/data/" ]] || {
	echo "Expected data directory on remote host" >&2
	return 1;
    }

    if [[ -d "/data/ref_genome" ]]; then return 0; fi

    local tmp=`mktemp -d -p /data/ sample-genome.dl.XXXXXX`
    (
	set -o pipefail
	trap "rm -rf --one-file-system $tmp" EXIT
	cd "$tmp"
	curl "${SAMPLEGENOME}" | tee >(tar -xzvf -) >(md5sum > ./sample.md5) > /dev/null
	expected=$(curl "${SAMPLEGENOME}.md5" | cut -f 1 -d' ')
	match=$(cut -f1 -d' ' ./sample.md5)
	if [[ "$match" != "$expected" ]]; then
	    return 1;
	fi
	mv ref_genome /data/
    )
}

function remote ()
{
    # SSH on remote host
    ssh -A "${RUSER}@${RHOST}" "$@"
}

while [[ "$#" -gt 0 ]]; do
    arg="$1"
    shift
    case "$arg" in
	--help|-h)
	    ( set +x
	      usage
	    )
	    exit 0
	    ;;
	--)
	    POSARGS+=("$@")
	    shift $#
	    ;;
	--step)
	    [[ $# -gt 0 ]] || { echo "Missing stepname" >&2; exit 1; }
	    STEPS+=("$1")
	    shift
	    ;;
	--host)
	    [[ $# -gt 0 ]] || { echo "Missing hostname" >&2; exit 1; }
	    RHOST="$1"
	    shift
	    ;;
	-*)
	    echo "invalid flag: $arg" >&2
	    exit 1
	    ;;
	*)
	    POSARGS+=("$arg")
	    ;;
    esac
done

if [[ ${#STEPS[@]} -lt 1 ]]; then STEPS=("${ALLSTEPS[@]}"); fi

function _cleanup ()
{
    local f
    for f in "$@"; do
	[[ -n "$f" ]] && ssh "${RUSER}@${RHOST}" rm --one-file-system "$f" || :
    done
}

if [[ -n "$RHOST" ]]; then
    # run remotely
    TEMPBIN=$(ssh "${RUSER}@${RHOST}" mktemp "$(basename "$0")".XXXXXX.sh)
    (
	trap "cleanup $TEMPBIN" EXIT
	REMOTE_ARGS=()
	for s in "${STEPS}"; do
	    REMOTE_ARGS+=(--step "$s")
	done
	remote tee "${TEMPBIN}" > /dev/null < "$0"
	remote bash -ex "${TEMPBIN}" "${REMOTE_ARGS[@]}"
    ) || {
	err="$?"
	exit $err
    }
else
    # local run
    for step in "${STEPS[@]}"; do
	if ! _has_step "$step"; then
	    echo "Invalid step: '$step'" >&2
	    echo "valid choices: ${ALLSTEPS[@]}" >&2
	    exit 1
	else
	    "$step"
	fi
    done
fi




