#!/bin/bash

CELL=$1

if [ -n "${CELL}" ]; then
    for y in 2013; do
        if [ -n "${sql}" ]; then
            sql=" UNION ${sql}"
        fi
        sql="SELECT * FROM romes.romes_${y} R JOIN network.cell C ON R.id_cell = C.id_cell AND C.name='${CELL}'\n${sql}"
    done
    #
    # get the coordinates
    #
    sql="SELECT Y(GK.location) AS X, X(GK.location) AS Y, rsrp FROM (${sql}) AS A JOIN gis.coordinate_gk GK ON A.id_coordinate = GK.id_coordinate WHERE rsrp IS NOT NULL ORDER BY X,Y;"
    #
    # create raster map
    #
    RASTER="field_meas_${CELL}"
    echo "Creating raster [${RASTER}] ..."
    #echo -e "${sql}" | psql -h mg_joe.mobitel.si -U garufa -d lte | tr -d ' ' | egrep '^[0-9]' > ${CELL}.xyz 
    echo -e "${sql}" | psql -h mg_joe.mobitel.si -U garufa -d lte | tr -d ' ' | egrep '^[0-9]' | r.in.xyz input=- output=${RASTER} method=median type=FCELL
else
    echo "Usage: $0 [cell name]"
    echo "Generates SQL statement for retrieving RSRP measurements of the given cell.-"
    exit -1
fi
