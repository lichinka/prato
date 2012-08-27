#!/bin/bash

TX=$1
OUT=$2
CMD="./r.coverage"

if [ -n "${TX}" ]; then
    CMD="${CMD} ini_file=$(pwd)/parameters.ini tx_ini_section=${TX}"
    if [ -n "${OUT}" ]; then
        CMD="${CMD} output_raster=${OUT}"
    fi
    ${CMD}
else
    echo "Usage"
    echo "  $0 [transmitter's section name in the INI file] [output raster name]"
    echo
fi
