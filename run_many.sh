#!/bin/bash

NPROC=$1
TX=CGBRIA,CGBRIB,CGBRIC,CSOSTA,CSOSTB,CSOSTC,CTOPOLA,CVELCEA,CVELCEB,SBANOVA
TX_LEN=10

if [ -n "${NPROC}" ]; then
    for i in $(seq 1 10);
    do
        TX=${TX},${TX}
    done
    TX=$(echo "${TX}" | cut -d',' -f 1-${NPROC})
    echo "./run_mpi.sh ${TX}"
else
    echo "Usage"
    echo "  $0 [number of processes]"
    echo
fi
