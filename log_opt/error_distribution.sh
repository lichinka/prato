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
    done
    #
    # calculate histogram ...
    #
    /usr/bin/env python3 ./log_opt/error_distribution.py ${TMP_FILE} > /tmp/.plot.dat
    #
    # ... and disply it
    #
    PLOT_CMD=$( cat <<EOF
set term postscript eps enhanced font "Helvetica,20";
set output "/tmp/error_distribution.eps";
plot '/tmp/.plot.dat' with boxes;
EOF
)
    gnuplot -p -e "${PLOT_CMD}"
else
    echo "Usage: $0 [cell] [E/// parameter values], [cell] [E/// parameter values], ..."
    echo "Generates the error distribution of the prediction against the field measurements for the given cells.-"
    echo "e.g."
    echo "	$0 CELL1 A0=38 A1=32 A2=-12 A3=0.1,CELL2 A0=39 A1=33 A2=-13 A3=0.2"
    exit -1
fi
