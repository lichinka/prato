#!/bin/bash

LOG_FILES=$@
BEST_TIME=99999

if [ $# -gt 1 ]; then
    for LOG in ${@:1}; do
        TIME=$( grep "Process took" ${LOG} | tr -d '[:alpha:]' | tr -d '[:blank:]' )
        if [ ${BEST_TIME} -gt ${TIME} ]; then
            BEST_FILE=${LOG}
            BEST_TIME=${TIME}
        fi
    done
    echo "${BEST_FILE}	${BEST_TIME}"
else
    echo -e "Usage:\t $0 [log files ...]"
    echo "Extracts the best running time from a group of simulation results."
fi
