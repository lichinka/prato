#!/bin/bash

LOG_FILES=$@

if [ $# -gt 0 ]; then
    echo "# NP	Spawning	Read_and_distribute	Processing_loop	Result_aggregation	Total_time"
    for LOG in ${@:1}; do
        NP=$( grep "Number of processors" ${LOG} | tr -dc '[:digit:]' )
        SPAWN_TIME=$( grep "process spawning" ${LOG} | cut -f3 | tr -d '[:alpha:]' )
        COMMON_TIME=$( grep "Common data" ${LOG} | cut -f3 | tr -d '[:alpha:]' )
        PROCESSING_TIME=$( grep "Calculation" ${LOG} | cut -f3 | sed -e 's/ sec/+/g' | tr -d '\n' | sed -e 's/+$/\n/g' | bc -l )
        PROCESSING_TIME=$( echo "${PROCESSING_TIME} / ${NP}" | bc -l )
        TOTAL_TIME=$( cat ${LOG} | grep "Process took" | tr -d '[:alpha:]' | tr -d '[:blank:]' )
        DIFF_TIME=$( echo "${TOTAL_TIME}-${PROCESSING_TIME}-${COMMON_TIME}-${SPAWN_TIME}" | bc -l )
        echo "${NP}	${SPAWN_TIME}	${COMMON_TIME}	${PROCESSING_TIME}	${DIFF_TIME}	${TOTAL_TIME}"
    done
else
    echo -e "Usage:\t $0 [log files ...]"
    echo "Extracts partial times from strong-scaling simulation results."
fi
