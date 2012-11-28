#!/bin/bash

HFILE=$1
TX=$2
OUT=$3
INI="$(pwd)/parameters.ini"
CMD="mpirun --mca btl tcp,sm,self --hostfile ${HFILE} -n 1 r.coverage -p"

if [ -n "${HFILE}" ] && [ -n "${TX}" ]; then
    CMD="${CMD} ini_file=${INI} tx_ini_sections=${TX}"
    if [ -n "${OUT}" ]; then
        CMD="${CMD} output_raster=${OUT}"
    fi
    echo ${CMD} >&2 && ${CMD}
    echo "Parallel process took $SECONDS sec"
    ./aggregate.sh ${TX} ${INI}
    echo "Total time $SECONDS sec"
else
    echo -e "Usage:\t$0 [hostfile to use] [transmitter's section name in the INI file] [output raster name]"
    echo "Runs the radio coverage calculation prediction over MPI."
    echo
fi
