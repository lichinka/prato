#!/bin/bash

NTX=$1
TX="$( cat meritve/lte/cells_??.txt )"

if [ -n "${NTX}" ]; then
    rm -f /tmp/.cells
    for i in $(seq 1 100);
    do
        echo "${TX}" >> /tmp/.cells
    done
    TX="$( cat /tmp/.cells | shuf -n${NTX} | tr '\n' ',' | sed -e 's/,$//' )"
    echo "${TX}"
else
    echo -e "Usage:\t$0 [number of transmitters]"
    echo "Returns a comma-separated list of [number of transmitters]."
fi
