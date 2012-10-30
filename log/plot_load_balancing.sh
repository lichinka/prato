#!/bin/bash

NTX=$( ./log_scale.sh 64 4096 )
NP=$( ./log_scale.sh 1 128 )

#
# Returns the files with the best execution time,
# for all processors, out of a set of runs for fixed NTX.-
#
# $1    Number of transmitters.-
#
function get_files {
    ntx=$1
    for np in ${NP}; do
        F=$( find strong_scaling -iname "NTX.${ntx}_NP.${np}_?.txt" )
        if [ -n "${F}" ]; then
            ./best_time.sh ${F};  
        fi
    done | cut -f1
}



######################
# Load balancing 
#
######################
echo -n "Calculating for "
for ntx in ${NTX}; do 
    echo -n "${ntx}...	"
    ./load_balancing_stats.sh ${ntx} $( get_files ${ntx} ) | sort -n -k2 > /tmp/${ntx}.dat
done
echo "done"

CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "strong_scaling/strong_scaling-load_balancing_plot.eps";
set title ".:. Load balancing .:.";
set key bottom left;
set xlabel "Number of cores";
set xtics 2;
set grid xtics;
set xrange [1:140];
set log x; 
set ylabel "Load balancing factor";
set ytics 0.05;
set yrange [0.70:1];
set grid ytics;
unset log y;
EOF
)
PLOT="plot "
for ntx in ${NTX}; do
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using 2:3 with lines lw 2 title \"${ntx} Transmitters\", "
done
PLOT="$( echo "${PLOT}" | sed -e 's/, $//g' )"
CMD="$( echo -e "${CMD}\n${PLOT}" )"
echo "${CMD}"
gnuplot -p -e "${CMD}"



################
# Time plot
#
################
for ntx in 64 128; do 
    ./load_balancing.py ${ntx} strong_scaling/NTX.${ntx}_NP.${ntx}_2.txt | grep -i 'total running time' | tr -d '[:alpha:]' | sed -e 's/[->()]//g' | sed -e 's/^[ ]*//g' > /tmp/${ntx}.dat
    CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "strong_scaling/per_worker_times-${ntx}_NTX.eps";
set title ".:. Processing times per worker .:.";
set key upper right;
set xlabel "Worker process ID";
set xtics 2;
set xrange [0:${ntx}];
unset log x; 
set ylabel "Wall clock time (sec)";
set ytics 1;
set yrange [0:14];
unset log y;
EOF
)
    PLOT="plot \"/tmp/${ntx}.dat\" with lines lw 3 title \"${ntx} Transmitters per core\"; "
    CMD="$( echo -e "${CMD}\n${PLOT}" )"
    echo "${CMD}"
    gnuplot -p -e "${CMD}"
done
