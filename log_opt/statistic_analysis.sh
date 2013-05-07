#!/bin/bash

REGION=$1
OUT=/tmp/.plot.dat

if [ -n "${REGION}" ]; then
    #
    # create a raster to contain all field measurements
    #
    CMD="null()"
    for c in $( cat meritve/lte/cells_${REGION}.txt ); do
        CMD="if(!isnull(field_meas_${c}),field_meas_${c},${CMD})"
    done
    CMD="$( echo "${CMD}" | sed -e 's/,$//g' )"
    CMD="r.mapcalc field_meas_${REGION}=\'${CMD}\'"
    CMD="$( echo "${CMD}" | tr -d '\\' )"
    echo "${CMD}" > /tmp/.run
    #
    # get clutter categories for the points with field measurements
    #
    echo "r.mapcalc temp=\'if(!isnull(clut_cat_25_no_road) && !isnull(field_meas_${REGION}),clut_cat_25_no_road,null())\'" | tr -d '\\' >> /tmp/.run
    echo "r.report -ni map=temp units=p" >> /tmp/.run
    echo "-------------------------------------------------"
    echo "Field-measurement proportion per clutter category"
    echo "-------------------------------------------------"
    sh /tmp/.run

    #
    # there are 12 different clutter categories
    #
    echo "-------------------------------------------------"
    echo "Clutter-category loss: solution stability"
    echo "-------------------------------------------------"
    echo "#x	box_min	whisker_min	whisker_high	box_high	mean" > ${OUT}
    for i in $( seq 0 11 ); do
        tail -n15 log_opt/${REGION}.0??  | grep "p${i}" | cut -d'=' -f2 | tr -d ' ' > /tmp/.data.dat
        /usr/bin/env python3 - <<EOF >> ${OUT}
import numpy as np 

v       = np.genfromtxt ('/tmp/.data.dat')
min     = np.min (v)
max     = np.max (v)
mean    = np.mean (v)
std_dev = np.std (v)

if max < 40:
    print ('%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f' % (${i}, mean-std_dev, min, max, mean+std_dev, mean))
EOF
    done
    #
    # create the box plot
    #
    PLOT_CMD=$( cat <<EOF
set term postscript eps enhanced font "Helvetica,20";
set output "log_opt/boxplot-${REGION}.eps";
set xlabel "Clutter category";
set xrange [-1:12];
set xtics (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11);
set ylabel "Signal loss due to clutter (dB)";
set yrange [3:25];
set ytics (0, 5, 10, 15, 20, 25, 30, 35, 40);
set grid ytics;
set boxwidth 0.5 absolute;
plot "${OUT}" using 1:2:3:4:5 with candlestick lt -1 lw 2 notitle whiskerbars,
     "${OUT}" using 1:6:6:6:6 with candlestick lt -1 lw 2 notitle;
EOF
)
    echo "${PLOT_CMD}"
    gnuplot -p -e "${PLOT_CMD}"
else
    echo "Usage: $0 [region]"
    echo "Generates the statistical analysis of the optimization results for the given region.-"
    echo "[region]	Region is one of MS,LJ,HI"
    echo
fi
