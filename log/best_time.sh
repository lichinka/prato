#!/bin/bash

DIR=$1
NP=$2
NRUNS=$3
BEST_TIME=99999

if [ -n "${DIR}" ] && [ -n "${NP}" ] && [ -n "${NRUNS}" ]; then
    INPUT="weak_scaling_${NP}"
    for i in $(seq -w 1 ${NRUNS}); do 
        F="${INPUT}.${i}.txt"
        TOTAL_TIME=$(tail -n1 ${DIR}/${F} | tr -d '[:alpha:]' | tr -d '[:blank:]')
        if [ ${BEST_TIME} -gt ${TOTAL_TIME} ]; then
            BEST_TIME=${TOTAL_TIME}
        fi
    done
    echo "Best time is ${BEST_TIME}"
else
    echo "Usage:"
    echo -e "\t $0 [dir containing log files] [number of processes] [number of runs]"
fi
