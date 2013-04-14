#!/bin/bash

OPTIONS=$1
CELLS=$2

if [ -n "${OPTIONS}" ] && [ -n "${CELLS}" ]; then
    #
    # create as many workers as there are cells
    #
    NP="$( echo "${CELLS}" | tr ',' '\n' | wc -l )"
    #
    # start MPI jobs for LTE network
    #
    mpirun --mca btl tcp,sm,self -np 1 -host localhost r.coverage $1 ini_file=./parameters_lte.ini tx_ini_sections=${CELLS} :  -np ${NP} --hostfile hostfile.local ./run_worker.sh
else
    echo "Usage: $0 [module flags] [comma-separated cell list]"
    echo "Runs the [r.coverage] tool module for the LTE network with the given parameters.-"
    ./r.coverage --help
    exit -1
fi
