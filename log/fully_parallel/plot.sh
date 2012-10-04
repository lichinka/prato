#!/bin/bash

WHAT=$1
OUTFILE=/tmp/.${WHAT}.dat

if [ -n "${WHAT}" ]; then
    rm -f ${OUTFILE}
    for i in 001 002 004 008 016 032 064 100 128; do
        T=$(./extract_min_time.sh out_${i}.txt ${WHAT})
        echo "${i}	${T}" >> ${OUTFILE}
    done 
else
    echo "Usage"
    echo "  $0 [keywords contained in lines to be plotted]"
    echo
fi
