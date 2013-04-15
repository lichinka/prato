#!/bin/bash

OUTPUT=$1
OPTIONS=$2
CELLS=$3

if [ -n "${OUTPUT}" ] && [ -n "${OPTIONS}" ] && [ -n "${CELLS}" ]; then
    if [ "${OUTPUT}" = "-" ]; then
        OUTPUT=" "
    fi
    #
    # create as many workers as there are cells
    #
    NP="$( echo "${CELLS}" | tr ',' '\n' | wc -l )"
    #
    # start MPI jobs for LTE network
    #
    mpirun --mca btl tcp,sm,self -np 1 -host localhost r.coverage ${OPTIONS} ini_file=./parameters_lte.ini tx_ini_sections=${CELLS} :  -np ${NP} --hostfile hostfile.local ./run_worker.sh ${OUTPUT}
else
    echo "Usage: $0 [output] [module flags] [comma-separated cell list]"
    echo "Runs the [r.coverage] tool module for the LTE network with the given parameters.-"
    echo "[output]	May be one of -db, -rast or - for standard output."
    ./r.coverage --help
    exit -1
fi
