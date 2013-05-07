#!/bin/bash

ALL_PARAMS=$@
TMP_FILE=/tmp/.diff.xyz

if [ -n "${ALL_PARAMS}" ]; then
    rm -f ${TMP_FILE}
    ALL_PARAMS="$( echo "${ALL_PARAMS}" | tr ' ' '|' )"
    #
    # each prediction has a cell name and E/// parameters
    #
    for PARAMS in $( echo "${ALL_PARAMS}" | tr ',' ' ' ); do
        CELL="$( echo "${PARAMS}" | tr '|' '\n' | head -n1 )"
        ERIC_PARAMS="$( echo "${PARAMS}" | tr '|' '\n' | tail -n+2 - )"
        #
        # calculate cell prediction
        #
        ./run_lte.sh -rast ${CELL} -m ${ERIC_PARAMS}
        g.copy rast=temp,prediction_${CELL} --overwrite
        #
        # transmit power of the given cell
        #
        PWR="$( cat parameters_lte.ini | grep -A7 ${CELL} | grep power | cut -d'=' -f2 | cut -d';' -f1 )"
        #
        # calculate the difference between prediction and measurements
        #
        r.mapcalc diff_${CELL}="abs(${PWR}-prediction_${CELL}-field_meas_${CELL})"
        r.out.xyz input=diff_${CELL} output=- >> ${TMP_FILE}
        #
        # clean unused rasters
        #
        g.remove rast=prediction_${CELL},diff_${CELL}
    done
    #
    # calculate histogram ...
    #
    /usr/bin/env python3 ./log_opt/error_distribution.py ${TMP_FILE} > /tmp/.plot.dat
    #
    # take the total to display the proportions in percent
    #
    TOT="$( cat /tmp/.plot.dat | cut -f2 | tr '\n' '+' )0"
    TOT="$( echo "${TOT}" | bc -l )"
    #
    # ... and display it
    #
    PLOT_CMD=$( cat <<EOF
set term postscript eps enhanced font "Helvetica,20";
set output "/tmp/error_distribution.eps";
set xlabel "Absolute error (dB)";
set xtics ("[0,5)" 0, "[5,10)" 1, "[10,15)" 2, "[15,20)" 3, "[20,25)" 4, "[25,30)" 5, "[30,35)" 6, "[35,...)" 7);
set ylabel "Proportion (%)";
set yrange [0:40];
set grid ytics;
plot '/tmp/.plot.dat' using ((\$2/${TOT})*100) with boxes lt -1 fs pattern 2 notitle;
EOF
)
    echo "${PLOT_CMD}"
    gnuplot -p -e "${PLOT_CMD}"
else
    echo "Usage: $0 [cell] [E/// parameter values], [cell] [E/// parameter values], ..."
    echo "Generates the cumulative error distribution of the prediction against the field measurements for the given cells, e.g."
    echo "	$0 CELL1 A0=38 A1=32 A2=-12 A3=0.1,CELL2 A0=39 A1=33 A2=-13 A3=0.2"
    exit -1
fi
