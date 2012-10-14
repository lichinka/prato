#!/bin/bash

NTX=4096


#################################
# Discriminated wall clock time
#
#################################
for ntx in $( ./log_scale.sh 64 ${NTX} ); do 
    ./discriminated_times.sh strong_scaling/NTX.${ntx}_NP.*.txt | sort -n -k1 > /tmp/${ntx}.dat
done

for ntx in 64 256 1024 4096; do
    CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "strong_scaling/discriminated_time_plot_${ntx}.eps";
set title ".:. Strong scalability -- ${ntx} Transmitters .:.";
set style fill pattern 1 border;
set xlabel "Number of cores";
set xtics 2;
set grid xtics;
set xrange [1:140];
set log x; 
set ylabel "Wall clock time (sec)";
set format y "10^%T";
set ytics 10;
set yrange [1:50000];
set log y;
EOF
    )
    PLOT="plot "
    PLOT="$( echo "${PLOT}" | sed -e 's/, $//g' )"
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using 1:(\$2+\$3+\$4+\$5) with filledcurve y1=0 title \"Result aggregation\", "
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using 1:(\$2+\$3+\$4) with filledcurve y1=0 title \"Processing loop\", "
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using 1:(\$2+\$3) with filledcurve y1=0 title \"Read and Broadcast common data\", "
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using 1:(\$2) with filledcurve y1=0 title \"Dynamic worker-process spawning\";"
    CMD="$( echo -e "${CMD}\n${PLOT}" )"
    echo "${CMD}"
    gnuplot -p -e "${CMD}"
done



######################
# Wall clock time
#
######################
for ntx in $( ./log_scale.sh 64 ${NTX} ); do 
    ./time_stats.sh ${ntx} strong_scaling/NTX.${ntx}_NP.*.txt | sort -n -k2 > /tmp/${ntx}.dat
done

CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "strong_scaling/time_plot.eps";
set title ".:. Strong scalability .:.";
set xlabel "Number of cores";
set xtics 2;
set grid xtics;
set xrange [1:140];
set log x; 
set ylabel "Wall clock time (sec)";
set format y "10^%T";
set ytics 10;
set yrange [60:30000];
set log y;
EOF
)
PLOT="plot "
for ntx in $( ./log_scale.sh 64 ${NTX} ); do
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using 2:3 with lines lw 2 title \"${ntx} Transmitters\", "
done
PLOT="$( echo "${PLOT}" | sed -e 's/, $//g' )"
CMD="$( echo -e "${CMD}\n${PLOT}" )"
echo "${CMD}"
gnuplot -p -e "${CMD}"



######################
# Speed-up 
#
######################
for ntx in $( ./log_scale.sh 64 ${NTX} ); do 
    ./speedup_stats.sh ${ntx} strong_scaling/NTX.${ntx}_NP.*.txt | sort -n -k2 > /tmp/${ntx}.dat
done

CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "strong_scaling/speedup_plot.eps";
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
for ntx in $( ./log_scale.sh 64 ${NTX} ); do
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using 2:3 with lines lw 2 title \"${ntx} Transmitters\", "
done
PLOT="${PLOT} x with lines lw 2 title \"Perfect speed-up\""
CMD="$( echo -e "${CMD}\n${PLOT}" )"
echo "${CMD}"
gnuplot -p -e "${CMD}"


######################
# Efficiency
#
######################
for ntx in $( ./log_scale.sh 64 ${NTX} ); do 
    ./efficiency_stats.sh ${ntx} strong_scaling/NTX.${ntx}_NP.*.txt | sort -n -k2 > /tmp/${ntx}.dat
done

CMD=$( cat <<EOF
set term postscript eps enhanced; 
set output "strong_scaling/efficiency_plot.eps";
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
for ntx in $( ./log_scale.sh 64 ${NTX} ); do
    PLOT="${PLOT} \"/tmp/${ntx}.dat\" using 2:3 with lines lw 2 title \"${ntx} Transmitters\", "
done
PLOT="$( echo "${PLOT}" | sed -e 's/, $//g' )"
CMD="$( echo -e "${CMD}\n${PLOT}" )"
echo "${CMD}"
gnuplot -p -e "${CMD}"
