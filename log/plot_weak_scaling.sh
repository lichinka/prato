#!/bin/bash

NTX=$( ./log_scale.sh 5 80 )
NP=$( ./log_scale.sh 1 128 )

for ntx in ${NTX}; do
    rm -f /tmp/${ntx}.dat
    touch /tmp/${ntx}.dat
    for np in  ${NP}; do
        BEST_TIME=$( ./best_time.sh weak_scaling_${ntx}_Tx/weak_scaling_${np}.*.txt | cut -f2 )
        echo "${np}	${BEST_TIME}" >> /tmp/${ntx}.dat
    done
done

################
# Time plot
#
################
CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "weak_scaling-time_plot.eps";
set title ".:. Weak scalability .:.";
set key bottom right;
set xlabel "Number of cores";
set xtics 2;
set grid xtics;
set xrange [1:140];
set log x; 
set ylabel "Wall clock time (sec)";
set format y "10^%T";
set ytics 10;
set yrange [40:1100];
set log y;
EOF
)
PLOT="plot "
for ntx in ${NTX}; do
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" with lines lw 3 title \"${ntx} Transmitters per core\", "
done
PLOT="$( echo "${PLOT}" | sed -e 's/, $//g' )"
CMD="$( echo -e "${CMD}\n${PLOT}" )"
echo "${CMD}"
gnuplot -p -e "${CMD}"



#################################
# Relative processing time
#
#################################
for ntx in ${NTX}; do 
    ./discriminated_times.sh $( for np in ${NP}; do ./best_time.sh weak_scaling_${ntx}_Tx/weak_scaling_${np}.*.txt; done | cut -f1 ) | sort -n -k1 > /tmp/${ntx}.dat
done

for ntx in 5 20 80; do
    CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "weak_scaling_relative_time_plot_${ntx}.eps";
set title ".:. Weak scalability -- ${ntx} Transmitters .:.";
set style data histograms;
set style histogram rowstacked;
set boxwidth 1 relative;
set style fill pattern 3 border;
set xlabel "Number of cores";
set nolog x; 
set ylabel "Relative processing time";
set ytics 0.1;
set yrange [0:1.20];
set nolog y;
EOF
    )
    PLOT="plot "
    PLOT="$( echo "${PLOT}" | sed -e 's/, $//g' )"
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using (\$2/\$6) title \"Dynamic worker-process spawning\", "
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using (\$3/\$6):xticlabels(1) title \"Input data broadcasting\", "
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using (\$4/\$6):xticlabels(1) title \"Processing loop\", "
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using (\$5/\$6):xticlabels(1) title \"Create final coverage prediction\"; "
    CMD="$( echo -e "${CMD}\n${PLOT}" )"
    echo "${CMD}"
    gnuplot -p -e "${CMD}"
done
