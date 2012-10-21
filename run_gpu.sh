#!/bin/bash

TX=$1
OUT=$2
CMD="./r.coverage"
PSQL_SERVER=localhost

if [ -n "${TX}" ]; then
    CMD="${CMD} -g ini_file=$(pwd)/parameters.ini tx_ini_section=${TX}"
    if [ -n "${OUT}" ]; then
        CMD="${CMD} output_raster=${OUT}"
    fi
    echo ${CMD} >&2 
    #${CMD} | psql -h ${PSQL_SERVER} -U garufa grass_backend
    ${CMD}
else
    echo "Usage"
    echo "  $0 [transmitter's section name in the INI file] [output raster name]"
    echo
fi
