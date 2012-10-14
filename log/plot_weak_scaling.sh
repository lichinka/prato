#!/bin/bash

NTX=$( ./log_scale.sh 5 80 )

for ntx in ${NTX}; do
    rm -f /tmp/${ntx}.dat
    touch /tmp/${ntx}.dat
    for np in  $( ./log_scale.sh 1 128 ); do
        #PROBLEM_SIZE=$( echo "${ntx}*${np}" | bc -l )
        BEST_TIME=$( ./best_time.sh weak_scaling_${ntx}_Tx/weak_scaling_${np}.*.txt )
        echo "${np}	${BEST_TIME}" >> /tmp/${ntx}.dat
    done
done

################
# Time plot
#
################
#CMD=$( cat <<EOF
#set term postscript eps enhanced; 
#set output "weak_scaling-time_plot.eps";
#set title ".:. Weak scalability .:.";
#set key bottom right;
#set xlabel "Number of cores";
#set xtics 2;
#set grid xtics;
#set xrange [1:140];
#set log x; 
#set ylabel "Wall clock time (sec)";
#set format y "10^%T";
#set ytics 10;
#set yrange [40:1100];
#set log y;
#EOF
#)
CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "weak_scaling-time_plot.eps";
set title ".:. Weak scalability .:.";
set key top right;
set style data histograms;
set style histogram rowstacked;
set boxwidth 1 relative;
set style fill solid 1.0 border -1;
set xlabel "Number of cores";
set xtics 2;
set grid xtics;
set xrange [1:140];
set nolog x; 
set ylabel "Wall clock time (sec)";
set format y "10^%T";
set ytics 10;
set yrange [40:1100];
set nolog y;
EOF
)
PLOT="plot "
for ntx in ${NTX}; do
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" title \"${ntx} Transmitters per core\", "
done
PLOT="$( echo "${PLOT}" | sed -e 's/, $//g' )"
CMD="$( echo -e "${CMD}\n${PLOT}" )"
echo "${CMD}"
gnuplot -p -e "${CMD}"
