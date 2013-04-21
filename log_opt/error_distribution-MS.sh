#!/bin/bash

#
# Error distribution with default model parameters
#
SCRIPT=""
for c in $( cat meritve/lte/cells_ms.txt ); do
    SCRIPT="${SCRIPT},${c} A0=38 A1=32 A2=-12.0 A3=0.1"
done
SCRIPT="$( echo "${SCRIPT}" | sed -e 's/^,//g' )"
SCRIPT="./log_opt/error_distribution.sh ${SCRIPT}"
${SCRIPT}
cp /tmp/error_distribution.eps /tmp/error_distribution-MS-def_params.eps

#
# Calculate the optimal parameters for each cell
#
for c in $( cat meritve/lte/cells_ms.txt ); do
    ./run_lte.sh - ${c} -mp
    PARAMS="$( grep -A4 'found optimal' /tmp/worker.log | grep 'A' | cut -f2,3 | tr '\t' '=' | tr '\n' ' ' )"
    echo "${c} ${PARAMS}" > /tmp/.${c}
done

#
# Error distribution with optimal model parameters
#
SCRIPT=""
for c in $( cat meritve/lte/cells_ms.txt ); do
    SCRIPT="${SCRIPT},$( cat /tmp/.${c} )"
done
SCRIPT="$( echo "${SCRIPT}" | sed -e 's/^,//g' )"
SCRIPT="./log_opt/error_distribution.sh ${SCRIPT}"
${SCRIPT}
cp /tmp/error_distribution.eps /tmp/error_distribution-MS-opt_params.eps
