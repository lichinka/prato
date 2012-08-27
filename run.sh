#!/bin/bash

TX=$1
OUT=$2

if [ -n "${TX}" ] && [ -n "${OUT}" ]; then
    ./r.coverage ini_file=$(pwd)/parameters.ini tx_ini_section=${TX} output_raster=${OUT}
else
    echo "Usage"
    echo "  $0 [transmitter's section name in the INI file] [output raster name]"
    echo
fi
