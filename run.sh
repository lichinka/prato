#!/bin/bash

TX=$1
OUT=$2
CMD="./r.coverage"

if [ -n "${TX}" ]; then
    CMD="${CMD} ini_file=$(pwd)/parameters.ini tx_ini_section=${TX}"
    if [ -n "${OUT}" ]; then
        CMD="${CMD} output_raster=${OUT}"
    fi
    #
    # dynamically create a table for the command output
    #
    CREATE_SQL=$(cat create.sql | sed -e "s/@table/coverage_$$/g" -)
    echo "${CREATE_SQL}"
    psql -h 192.168.1.160 -U garufa grass_backend -c "${CREATE_SQL}"
    COPY_SQL=$(cat copy.sql | sed -e "s/@table/coverage_$$/g" -)
    echo "${COPY_SQL}"
    echo ${CMD} >&2 && ${CMD} | psql -h 192.168.1.160 -U garufa grass_backend -c "${COPY_SQL}"
else
    echo "Usage"
    echo "  $0 [transmitter's section name in the INI file] [output raster name]"
    echo
fi
