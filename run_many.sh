#!/bin/bash

NTX=$1
TX=CGBRIA,CGBRIB,CGBRIC,CSOSTA,CSOSTB,CSOSTC,CTOPOLA,CVELCEA,CVELCEB,SBANOVA
TX_LEN=10

if [ -n "${NTX}" ]; then
    for i in $(seq 1 10);
    do
        TX="${TX},${TX}"
    done
    TX=$( echo "${TX}" | cut -d',' -f 1-${NTX} )
    echo "${TX}"
else
    echo -e "Usage:\t$0 [number of transmitters]"
    echo "Returns a comma-separated list of [number of transmitters]."
fi
