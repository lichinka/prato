#!/bin/bash

HFILE=$1
TX=$2
OUT=$3
INI="$(pwd)/parameters.ini"
CMD="mpirun --mca btl tcp,sm,self --hostfile ${HFILE} -n 1 r.coverage -p"

if [ -n "${HFILE}" ] && [ -n "${TX}" ]; then
    CMD="${CMD} ini_file=${INI} tx_ini_sections=${TX}"
    echo ${CMD} >&2 && ${CMD}
    FINISHED=$?
    echo "Parallel process took $SECONDS sec"
    if [ ${FINISHED} -eq 0 ]; then
        if [ -z "${OUT}" ]; then
            OUT="final_coverage"
        fi
        ./aggregate.sh ${TX} ${INI} ${OUT}
        echo "Total time $SECONDS sec"
    fi
else
    echo -e "Usage:\t$0 [hostfile to use] [transmitter's section name in the INI file] [output raster name]"
    echo "Runs the radio coverage calculation prediction over MPI, generating a raster map."
    echo
fi
