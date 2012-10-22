#!/bin/bash


NTX=$1
LOG_FILES=$@
BASE_CASE=''

if [ $# -gt 2 ]; then
    for LOG in ${@:2}; do
        if [ -z "${BASE_CASE}" ]; then
            BASE_CASE="$( echo ${LOG} | grep 'NP\.1\_' - )"
            echo "${BASE_CASE}"
        fi
    done
    BASE_CASE=$( cat ${BASE_CASE} | grep 'Process took' | tr -d '[:alpha:]' )
    for LOG in ${@:2}; do
        NP=$( grep 'Number of processors' ${LOG} | tr -dc '[:digit:]' )
        TIME=$( cat ${LOG} | grep 'Process took' | tr -d '[:alpha:]' )
        SPEEDUP=$( echo "(${BASE_CASE}/${TIME})/${NP}" | bc -l )
        echo "${NTX}	${NP}	${SPEEDUP}"
    done
else
    echo "Usage:"
    echo -e "\t $0 [number of transmitters] [log files ...]"
fi
