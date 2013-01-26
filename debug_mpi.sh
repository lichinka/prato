#!/bin/bash

if [ -n "$1" ]; then
    gdb --pid=$1 --command=commands_for_mpi_debug.gdb
else
    echo "Attaches the debugger to a running MPI session."
    echo "Usage: $0 [pid]"
    echo "pid	PID of the MPI process to attach to"
    echo
    exit -1
fi
