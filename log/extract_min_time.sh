#!/bin/bash

LOG=$1
WHAT=$2

if [ -n "${WHAT}" ] && [ -n "${LOG}" ]; then
    TIMES="$(cat ${LOG} | grep ${WHAT} | cut -f3 | cut -d' ' -f1 | sort -n)"
    LINES_IN_LOG=$(wc -l ${LOG} | cut -d' ' -f1)
    LINES_IN_LOG=$(( ${LINES_IN_LOG} / 3 ))
    NTIMES=$(echo "${TIMES}" | wc -l)
    if [ ${LINES_IN_LOG} -eq ${NTIMES} ]; then
        MIN_TIME=$(echo "${TIMES}" | head -n1)
        echo $MIN_TIME
    else
        echo "Lines in log file ${LINES_IN_LOG}"
        echo "Lines after grepping ${NTIMES}"
        echo "ERROR Give me more keywords to uniquely identify the target line"
    fi
else
    echo "Usage"
    echo "  $0 [log file to plot] [keywords contained in lines to be plotted]"
    echo
fi
