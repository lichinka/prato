#!/bin/bash


NTX=$1
LOG_FILES=$@

if [ $# -gt 2 ]; then
    for LOG in ${@:2}; do
        NP=$( grep 'Number of processors' ${LOG} | tr -dc '[:digit:]' )
        TIME=$( cat ${LOG} | grep 'Process took' | tr -d '[:alpha:]' | tr -d '[:blank:]' )
        echo "${NTX}	${NP}	${TIME}"
    done
else
    echo "Usage:"
    echo -e "\t $0 [number of transmitters] [log files ...]"
fi
