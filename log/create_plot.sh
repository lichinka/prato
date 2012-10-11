#!/bin/bash

NTX=$@

#set term postscript eps enhanced; 
#set output "weak_scaling_${NTX}_Tx/weak_scaling_plot.eps";
if [ $# -gt 1 ]; then
    CMD=$( cat <<EOF
set xlabel "Number of processes";
set xtics 2;
set grid xtics;
set xrange [1:132];
set log x; 
set ylabel "Wall clock time (sec)";
set yrange [90:];
set log y;
EOF
    )
    PLOT="plot "
    for n in ${NTX}; do
        PLOT="${PLOT} \"weak_scaling_${n}_Tx/best_times.log\" with lines lw 2 title \"${n} Tx per process\", "
    done
    PLOT="$( echo "${PLOT}" | sed -e 's/, $//g' )"
    CMD="$( echo -e "${CMD}\n${PLOT}" )"
    echo "${CMD}"
    gnuplot -p -e "${CMD}"
else
    echo "Usage:"
    echo -e "\t $0 [number of transmitters per process] ..."
fi
