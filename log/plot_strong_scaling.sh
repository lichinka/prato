#!/bin/bash

NTX=$(./log_scale.sh 64 4096)
NP=$( ./log_scale.sh 1 128 )

#
# Returns the files with the best execution time,
# out of a set of runs. 
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
# Speed-up 
#
######################
for ntx in ${NTX}; do 
    ./speedup_stats.sh ${ntx} $( get_files ${ntx} ) | sort -n -k2 > /tmp/${ntx}.dat
done

CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "strong_scaling/strong_scaling-speedup_plot.eps";
set title ".:. Strong scalability .:.";
set key top left;
set xlabel "Number of cores";
set xtics 2;
set grid xtics;
set xrange [1:140];
set log x; 
set ylabel "Speed-up";
set ytics 2;
set yrange [1:150];
set grid ytics;
set log y;
EOF
)
PLOT="plot "
for ntx in ${NTX}; do
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using 2:3 with lines lw 2 title \"${ntx} Transmitters\", "
done
PLOT="${PLOT} x with lines lw 2 title \"Perfect speed-up\""
CMD="$( echo -e "${CMD}\n${PLOT}" )"
echo "${CMD}"
gnuplot -p -e "${CMD}"



#################################
# Relative processing time
#
#################################
for ntx in ${NTX}; do 
    ./discriminated_times.sh $( get_files ${ntx} ) | sort -n -k1 > /tmp/${ntx}.dat
done

for ntx in ${NTX}; do
    CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "strong_scaling/strong_scaling-relative_time_plot_${ntx}.eps";
set title ".:. Strong scalability -- ${ntx} transmitters .:.";
set style data histograms;
set style histogram rowstacked;
set boxwidth 1 relative;
set style fill pattern 3 border;
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



######################
# Wall clock time
#
######################
for ntx in ${NTX}; do 
    ./time_stats.sh ${ntx} $( get_files ${ntx} ) | sort -n -k2 > /tmp/${ntx}.dat
done

CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "strong_scaling/strong_scaling-time_plot.eps";
set title ".:. Strong scalability .:.";
set xlabel "Number of cores";
set xtics 2;
set grid xtics;
set xrange [1:140];
set log x; 
set ylabel "Wall-clock time [s]";
set format y "10^%T";
set ytics 10;
set yrange [60:50000];
set log y;
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



######################
# Efficiency
#
######################
for ntx in ${NTX}; do 
    ./efficiency_stats.sh ${ntx} $( get_files ${ntx} ) | sort -n -k2 > /tmp/${ntx}.dat
done

CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "strong_scaling/strong_scaling-efficiency_plot.eps";
set title ".:. Strong scalability .:.";
set key bottom left;
set xlabel "Number of cores";
set xtics 2;
set grid xtics;
set xrange [1:140];
set log x; 
set ylabel "Parallel efficiency";
set ytics 0.1;
set yrange [0:1];
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

