#!/bin/bash

OUTPUT=$1
CELLS=$2
OPTIONS=${@:3}

if [ -n "${OUTPUT}" ] && [ -n "${CELLS}" ] && [ -n "${OPTIONS}" ]; then
    #
    # create as many workers as there are cells
    #
    NP="$( echo "${CELLS}" | tr ',' '\n' | wc -l )"
    #
    # start MPI jobs for LTE network
    #
    mpirun --mca btl_tcp_if_include 192.168.1.0/24 --mca btl tcp,sm,self --host localhost -np 1 r.coverage ${OPTIONS} ini_file=./parameters_lte.ini tx_ini_sections=${CELLS} : --hostfile hostfile.ninestein -np ${NP} ${HOME}/etc/dr/tun_par/prato/src/run_worker.sh ${OUTPUT}
else
    echo "Usage: $0 [output] [comma-separated cell list] [module opts ...]"
    echo "Runs the [r.coverage] tool module for the LTE network with the given parameters.-"
    echo "[output]	May be one of -db, -rast or - for standard output."
    echo "[cell list]	A comma-separated list of cell names."
    echo "[module opts]	Any combination of flags and options for the 'r.coverage' GRASS module. Help follows."
    ./r.coverage --help
    exit -1
fi