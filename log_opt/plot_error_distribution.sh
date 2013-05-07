#!/bin/bash

REGION=$1

if [ -n "${REGION}" ]; then
    #
    # Error distribution with default model parameters
    #
    SCRIPT=""
    for c in $( cat meritve/lte/cells_${REGION}.txt ); do
        SCRIPT="${SCRIPT},${c} A0=38 A1=32 A2=-12.0 A3=0.1"
    done
    SCRIPT="$( echo "${SCRIPT}" | sed -e 's/^,//g' )"
    SCRIPT="./log_opt/error_distribution.sh ${SCRIPT}"
    ${SCRIPT}
    cp /tmp/error_distribution.eps ./log_opt/error_distribution-${REGION}-def_params.eps

    #
    # Calculate the optimal parameters for each cell
    #
    for c in $( cat meritve/lte/cells_${REGION}.txt ); do
        rm -f /tmp/worker*.log
        ./run_lte.sh - ${c} -m p=1
        PARAMS="$( grep -A4 'found optimal' /tmp/worker*.log | grep 'A' | cut -f2,3 | tr '\t' '=' | tr '\n' ' ' )"
        echo "${c} ${PARAMS}" > /tmp/.${c}
    done

    #
    # Error distribution with optimal model parameters
    #
    SCRIPT=""
    for c in $( cat meritve/lte/cells_${REGION}.txt ); do
        SCRIPT="${SCRIPT},$( cat /tmp/.${c} )"
    done
    SCRIPT="$( echo "${SCRIPT}" | sed -e 's/^,//g' )"
    SCRIPT="./log_opt/error_distribution.sh ${SCRIPT}"
    ${SCRIPT}
    cp /tmp/error_distribution.eps ./log_opt/error_distribution-${REGION}-opt_params.eps
else
    echo "Usage: $0 [region]"
    echo "Generates error-distribution plots (with default and fitted parameters) for the given region.-"
    echo "[region]	Region is one of MS,LJ,HI"
    echo
fi
