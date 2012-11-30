#!/bin/bash

NTX=$( ./log_scale.sh 5 80 )
NP=$( ./log_scale.sh 1 128 )


#
# find the fastest run, i.e. the one with the minimum total time
#
for ntx in ${NTX}; do
    rm -f /tmp/${ntx}.dat
    touch /tmp/${ntx}.dat
    for np in  ${NP}; do
        BEST_TIME=$( ./best_time.sh weak_scaling/NTX.${ntx}_NP.${np}_*.txt | cut -f2 )
        echo "${np}	${BEST_TIME}" >> /tmp/${ntx}.dat
    done
done

################
# Time plot
#
################
CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "weak_scaling/weak_scaling-time_plot.eps";
set title ".:. Weak scalability .:.";
set key bottom right;
set xlabel "Number of cores";
set xtics 2;
set grid xtics;
set xrange [1:140];
set log x; 
set ylabel "Wall-clock time [s]";
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
    ./discriminated_times.sh $( for np in ${NP}; do ./best_time.sh weak_scaling/NTX.${ntx}_NP.${np}_*.txt; done | cut -f1 ) | sort -n -k1 > /tmp/${ntx}.dat
done

for ntx in 5 10 20 80; do
    CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "weak_scaling/weak_scaling_relative_time_plot_${ntx}.eps";
set title ".:. Weak scalability -- ${ntx} transmitters per core .:.";
set style data histograms;
set style histogram rowstacked;
set boxwidth 1 relative;
set style fill pattern 2 border;
set xlabel "Number of cores";
set nolog x; 
set ylabel "Relative processing time";
set ytics (0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0);
set yrange [0:1.25];
set nolog y;
EOF
    )
    PLOT="plot "
    PLOT="$( echo "${PLOT}" | sed -e 's/, $//g' )"
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using (\$2/\$7) title \"Read input data\", "
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using (\$3/\$7):xticlabels(1) title \"Dynamic worker-process spawning\", "
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using (\$4/\$7):xticlabels(1) title \"Input data broadcasting\", "
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using (\$5/\$7):xticlabels(1) title \"Processing loop\", "
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using (\$6/\$7):xticlabels(1) title \"Create final coverage prediction\"; "
    CMD="$( echo -e "${CMD}\n${PLOT}" )"
    echo "${CMD}"
    gnuplot -p -e "${CMD}"
done
