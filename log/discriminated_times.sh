#!/bin/bash

LOG_FILES=$@

if [ $# -gt 0 ]; then
    echo "# NP	Read	Spawning	Distribute	Processing_loop	Result_aggregation	Total_time"
    for LOG in ${@:1}; do
        NP=$( grep "Number of processors" ${LOG} | tr -dc '[:digit:]' )
        READ_TIME=$( grep "Read input data" ${LOG} | cut -f3 | tr -d '[:alpha:]' )
        SPAWN_TIME=$( grep "process spawning" ${LOG} | cut -f3 | tr -d '[:alpha:]' )
        DIST_TIME=$( grep "Common data" ${LOG} | cut -f3 | tr -d '[:alpha:]' )
        PROCESSING_TIME=$( grep "Calculation" ${LOG} | cut -f3 | sed -e 's/ sec/+/g' | tr -d '\n' | sed -e 's/+$/\n/g' | bc -l )
        PROCESSING_TIME=$( echo "${PROCESSING_TIME} / ${NP}" | bc -l )
        TOTAL_TIME=$( cat ${LOG} | grep "Process took" | tr -d '[:alpha:]' | tr -d '[:blank:]' )
        TOTAL_TIME=$( echo "${TOTAL_TIME} * 1.05" | bc -l ) 
        DIFF_TIME=$( echo "${TOTAL_TIME}-${PROCESSING_TIME}-${DIST_TIME}-${SPAWN_TIME}-${READ_TIME}" | bc -l )
        echo "${NP}	${READ_TIME}	${SPAWN_TIME}	${DIST_TIME}	${PROCESSING_TIME}	${DIFF_TIME}	${TOTAL_TIME}"
    done
else
    echo -e "Usage:\t $0 [log files ...]"
    echo "Extracts partial times from weak/strong-scaling simulation results."
fi
