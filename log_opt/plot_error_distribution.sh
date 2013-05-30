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
    echo "${SCRIPT}"
    cp /tmp/error_distribution.eps ./log_opt/error_distribution-${REGION}-def_params.eps

    #
    # Error distribution with optimal model parameters
    #
    SCRIPT=""
    for c in $( cat meritve/lte/cells_${REGION}.txt ); do
        A0="$( cat log_opt/${REGION}_opt_params.dat | grep "${c}" | tr -s ' ' | cut -d' ' -f2 )"
        A1="$( cat log_opt/${REGION}_opt_params.dat | grep "${c}" | tr -s ' ' | cut -d' ' -f3 )"
        A2="$( cat log_opt/${REGION}_opt_params.dat | grep "${c}" | tr -s ' ' | cut -d' ' -f4 )"
        A3="$( cat log_opt/${REGION}_opt_params.dat | grep "${c}" | tr -s ' ' | cut -d' ' -f5 )"
        SCRIPT="${SCRIPT},${c} A0=${A0} A1=${A1} A2=${A2} A3=${A3}"
    done
    SCRIPT="$( echo "${SCRIPT}" | sed -e 's/^,//g' )"
    SCRIPT="./log_opt/error_distribution.sh ${SCRIPT}"
    echo "${SCRIPT}"
    cp /tmp/error_distribution.eps ./log_opt/error_distribution-${REGION}-opt_params.eps
else
    echo "Usage: $0 [region]"
    echo "Generates error-distribution plots (with default and fitted parameters) for the given region.-"
    echo "[region]	Region is one of MS,LJ,HI"
    echo
fi
