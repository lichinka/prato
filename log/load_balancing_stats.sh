#!/bin/bash

NTX=$1
LOG_FILES=$@

if [ $# -gt 1 ]; then
    for LOG in ${@:2}; do
        NP=$( grep -i 'number of processors' ${LOG} | cut -f2 )
        FASTEST_WORKER=$( ./load_balancing.py ${NP} ${LOG} | grep fastest | cut -f2 )
        MASTER_TOTAL_TIME=$( tail ${LOG} | grep -i 'process took' | tr -d '[:alpha:]' )
        MASTER_PROC_TIME="$( grep '^TIME' ${LOG} | grep -v '(' | cut -f3 | tr -d '[:alpha:]' | tr '\n' '+' ) 0.0"
        MASTER_PROC_TIME=$( echo "${MASTER_TOTAL_TIME} - ( ${MASTER_PROC_TIME} )" | bc -l )
        LB=$( echo "${FASTEST_WORKER} / ${MASTER_PROC_TIME}" | bc -l )
        echo "${NTX}	${NP}	${LB}"
    done
else
    echo "Usage:"
    echo -e "\t $0 [number of transmitters] [log files ...]"
fi

