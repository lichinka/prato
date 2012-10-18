#!/bin/bash

HFILE=$1
NTX=$2
TX=CGBRIA,CGBRIB,CGBRIC,CSOSTA,CSOSTB,CSOSTC,CTOPOLA,CVELCEA,CVELCEB,SBANOVA
TX_LEN=10

if [ -n ${HFILE} ] && [ -n "${NTX}" ]; then
    for i in $(seq 1 10);
    do
        TX="${TX},${TX}"
    done
    TX=$( echo "${TX}" | cut -d',' -f 1-${NTX} )
    echo "./run_mpi.sh ${HFILE} ${TX}"
else
    echo -e "Usage:\t$0 [hostfile to use] [number of transmitters]"
    echo "Creates a command for processing [number of transmitters] over MPI."
fi
