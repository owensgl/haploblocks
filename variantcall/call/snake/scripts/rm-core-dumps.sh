#!/bin/bash

# remove core dumps from given folders
set -eu

if [[ "$#" -lt 1 ]]; then
    FOLDERS=( . )
else
    FOLDERS=( "$@" )
fi

for folder in "${FOLDERS[@]}"; do
    find "$folder" -xdev -regextype egrep -regex '.*/core[.][0-9]+' -type f -print0 | \
	xargs -0 --no-run-if-empty rm -v --one-file-system -- || :
done
