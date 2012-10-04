#!/bin/bash

TX=$1
OUT=$2
#CMD="mpirun -x LD_LIBRARY_PATH --mca btl_tcp_endpoint_cache 122880 --hostfile hostfile.local -n 1 r.coverage -p"
CMD="mpirun --mca btl tcp,sm,self --hostfile hostfile.local -n 1 r.coverage -p"

if [ -n "${TX}" ]; then
    CMD="${CMD} ini_file=$(pwd)/parameters.ini tx_ini_sections=${TX}"
    if [ -n "${OUT}" ]; then
        CMD="${CMD} output_raster=${OUT}"
    fi
    echo ${CMD} >&2 && ${CMD}
    echo "Process took $SECONDS sec"
else
    echo "Usage"
    echo "  $0 [transmitter's section name in the INI file] [output raster name]"
    echo
fi
