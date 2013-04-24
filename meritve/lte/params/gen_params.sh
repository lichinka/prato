#!/bin/bash

for c in $( sort ../cells_lj.txt ); do
    DATA="$( grep -h -B2 "info: \"${c}" *.params )"

    NAME="${c}"
    BEAM_DIR="$( echo "${DATA}" | egrep -o 'Azimut: [0-9]+' | egrep -o '[0-9]+' )"
    MECH="$( echo "${DATA}" | egrep -o 'Nagib: [0-9]+ \+ [0-9]+' | egrep -o '[0-9]+' | tail -n2 | head -n1 )"
    ELECT="$( echo "${DATA}" | egrep -o 'Nagib: [0-9]+ \+ [0-9]+' | egrep -o '[0-9]+' | tail -n1 )"
    HEIGHT="$( echo "${DATA}" | egrep -o 'Višina: [0-9]+' | egrep -o '[0-9]+' )"
    EAST="$( echo "${DATA}" | grep 'x:' | egrep -o '[0-9]+' )"
    NORTH="$( echo "${DATA}" | grep 'y:' | egrep -o '[0-9]+' )"
    ANT="$( echo "${DATA}" | egrep -o 'Tip:.+÷ K' | tr -d ' ' | egrep -o '[0-9]+' )"
    ANT="${ANT}/$( ls ../../../antene/*${ANT}* | grep 1855 | grep 0${ELECT}T )"
    PWR="$( echo "${DATA}" | egrep -o 'ant=[0-9]+' | egrep -o '[0-9]+' )"
    PWR="$( echo "${PWR}/10.0" | bc -l )"
    MEAS="field_meas_${NAME}@tun_par"

    echo "[${NAME}]"
    echo "cellName = ${NAME}"
    echo "beamDirection = ${BEAM_DIR}"
    echo "electricalTiltAngle = ${ELECT}"
    echo "mechanicalTiltAngle = ${MECH}"
    echo "heightAGL = ${HEIGHT}"
    echo "antennaFile = ${ANT}"
    echo "positionEast = ${EAST}"
    echo "positionNorth = ${NORTH}"
    echo "power = ${PWR}"
    echo "measurementsMap = ${MEAS}"
    echo
done

