#!/bin/bash

TX=$1
INI=$2
RAST=$3
#PSQL_SERVER=ninestein
PSQL_SERVER=localhost
#PSQL_USER=grassuser
PSQL_USER=garufa
PSQL_DB=grass

if [ -n "${TX}" ] && [ -n "${INI}" ] && [ -n "${RAST}" ]; then
    echo "*** INFO: Aggregating partial predictions for ${TX} ..."
    SQL="SELECT x, y, MAX(rscp) AS rscp FROM ("
    TX="$( echo "${TX}" | tr ',' ' ' )"
    for tx in ${TX}; do
        PWR="$( grep -A 9 ${tx} ${INI}  | grep power | cut -d'=' -f2 | cut -d';' -f1 | tr -d ' ' )"
        SQL="${SQL} SELECT east::integer-(east::integer % 25) AS x, north::integer-(north::integer % 25) AS y, (${PWR} - pl) AS rscp FROM pathloss_${tx} WHERE pl > 0 UNION"
    done 
    SQL="$( echo "${SQL}" | sed -e 's/UNION$//g' )"
    SQL="${SQL}) AS agg GROUP BY x, y "
    echo ${SQL}

    echo "*** INFO: Importing path-loss predictions from the database ..."
    echo -e "\t${SQL}" | psql -q -t -h ${PSQL_SERVER} -U ${PSQL_USER} -d ${PSQL_DB} | tr -d ' ' | gzip -c - > /tmp/.prediction.dat.gz
    gunzip -c /tmp/.prediction.dat.gz | v.in.ascii -t output=temp format=point -z z=3 --overwrite
    echo "*** INFO: Converting vector to raster map ..."
    v.to.rast input=temp type=point output=${RAST} use=z --overwrite
    #r.colors -n map=${RAST} color=elevation
else
    echo -e "Usage:\t$0 [comma-separated transmitter's section names] [INI file] [raster name]"
    echo "Aggregates individual path-loss predictions for the given transmitters into a raster map."
    echo
fi
