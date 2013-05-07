#!/bin/bash

REGION=$1

if [ -n "${REGION}" ]; then
    #
    # all field measurements for this region
    #
    FIELD_MEAS=""
    for c in $( cat meritve/lte/cells_${REGION}.txt ); do
        FIELD_MEAS="${FIELD_MEAS} meritve/lte/${c}.xyz"
    done
    FIELD_MEAS="cat ${FIELD_MEAS}"
    #
    # the minimum and maximum eastern coordinates
    #
    MIN_E="$( ${FIELD_MEAS} | cut -d'|' -f1 | sort -n | head -n1 )"
    MAX_E="$( ${FIELD_MEAS} | cut -d'|' -f1 | sort -nr | head -n1 )"
    #
    # the minimum and maximum northern coordinates
    #
    MIN_N="$( ${FIELD_MEAS} | cut -d'|' -f2 | sort -n | head -n1 )"
    MAX_N="$( ${FIELD_MEAS} | cut -d'|' -f2 | sort -nr | head -n1 )"
    #
    # calculate the aggregated prediction for all cells within this region
    #
    #CELLS="$( cat meritve/lte/cells_${REGION}.txt| tr '\n' ',' | sed -e 's/.$//g' )"
    #./run_lte.sh -db ${CELLS} -m
    #
    # set the region 
    #
    g.region n=${MAX_N} e=${MAX_E} s=${MIN_N} w=${MIN_E}
    g.region -p
    exit -1

    echo "-------------------"
    echo "Clutter proportions"
    echo "-------------------"
    r.report -ni map=clut_cat_25_no_road@tun_par units=p
    
    #
    # calculate the surface covered by the prediction of all cells
    # combined, within this region (in pixels)
    #
    echo "-----------------------------"
    echo "Field-measurement proportions"
    echo "-----------------------------"
    PREDICTION="$( r.report -ni map=final_coverage units=c | grep 'TOTAL' | cut -d'|' -f3 )"
    MEASUREMENTS="$( r.report -ni map=field_meas_${REGION} units=c | grep TOTAL | cut -d'|' -f3 )"
    echo "Pixels covered by the prediction:     ${PREDICTION}"
    echo "Pixels covered by field measurements: ${MEASUREMENTS}"
    echo "Proportion:                           $( echo "(${MEASUREMENTS}/${PREDICTION})*100" | bc -l )%"
    #
    # reset the region
    #
    g.region rast=clut_cat_25_no_road
else
    echo "Usage: $0 [region]"
    echo "Displays the clutter-category and field-measurement proportions for the given region.-"
    echo "[region]	Region is one of MS,LJ,HI"
    echo
fi
