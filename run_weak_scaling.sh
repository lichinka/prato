#!/bin/bash

NP=$1
OUTPUT=$2
CELLS=$3
OPTIONS=${@:4}
INI=./parameters_lte.ini

if [ -n "${NP}" ] && [ -n "${OUTPUT}" ] && [ -n "${CELLS}" ] && [ -n "${OPTIONS}" ]; then
    #
    # start MPI jobs for LTE network
    #
    mpirun --mca btl_tcp_if_include 192.168.1.0/24 --mca btl tcp,sm,self --host localhost -np 1 r.coverage ${OPTIONS} ini_file=${INI} tx_ini_sections=${CELLS} : --hostfile hostfile.ninestein -np ${NP} ${HOME}/etc/dr/tun_par/prato/src/run_worker.sh ${OUTPUT}
    #
    # aggregate the partial prediction only if output is database
    #
    if [ "${OUTPUT}" = "-db" ]; then
        AGG="$( { time ./aggregate.sh ${CELLS} ${INI} final_coverage; } 2>&1 )"
	echo "${AGG}"
    fi
else
    echo "Usage: $0 [worker] [output] [comma-separated cell list] [module opts ...]"
    echo "Runs the [r.coverage] tool module for the LTE network with the given parameters.-"
    echo "[workers]	Number of worker processes to start."
    echo "[output]	May be one of -db, -rast or - for standard output."
    echo "[cell list]	A comma-separated list of cell names."
    echo "[module opts]	Any combination of flags and options for the 'r.coverage' GRASS module. Help follows."
    ./r.coverage --help
    exit -1
fi
