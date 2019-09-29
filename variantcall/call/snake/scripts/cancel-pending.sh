#!/bin/bash

set -euxo pipefail

ACCOUNT="${ACCOUNT:-rrg-rieseber-ac_cpu}"

scancel -t PENDING -f -u "$(whoami)" -A "${ACCOUNT}" -vvv
exit $?

# while read jobid state rest; do
#     case "${state}" in
# 	PENDING)
# 	    scancel "$jobid"
# 	    sleep 1
# 	    ;;
# 	*)
# 	    ;;
#     esac
# done < <(sacct -b -n | grep -v '[.]')
