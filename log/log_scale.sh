#!/bin/bash

BEG=$1
END=$2

if [ -n "${BEG}" ] && [ -n "${END}" ]; then
    CURR=${BEG}
    while [ ${CURR} -le ${END} ]; do
        echo "${CURR}"
        CURR=$( echo "${CURR}*2" | bc -l )
    done
else
    echo "Usage:"
    echo -e "\t $0 [begin of LOG() series] [end of LOG() series]"
fi

