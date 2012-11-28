#!/bin/bash

TX=$1
INI=$2
PSQL_SERVER=localhost

if [ -n "${TX}" ] && [ -n "${INI}" ]; then
    echo "*** INFO: Aggregating partial predictions for ${TX} ..."
    SQL="SELECT east, north, rscp FROM ( SELECT east, north, MAX(rscp) AS rscp FROM ("
    TX="$( echo "${TX}" | tr ',' ' ' )"
    for tx in ${TX}; do
        PWR="$( grep -A 9 ${tx} ${INI}  | grep power | cut -d'=' -f2 | cut -d';' -f1 | tr -d ' ' )"
        SQL="${SQL} SELECT east, north, (${PWR} - pl) AS rscp FROM coverage_${tx} WHERE pl > 0 UNION"
    done
    SQL="$( echo "${SQL}" | sed -e 's/UNION$//g' )"
    SQL="${SQL}) AS agg GROUP BY agg.east, agg.north ) AS sub "
    SQL="${SQL} WHERE rscp < 0"
    echo "${SQL}" | psql -q -t -h ${PSQL_SERVER} -U garufa grass_backend | tr -d ' ' | v.in.ascii -t output=temp format=point -z z=3 --overwrite
    v.to.rast input=temp type=point output=temp use=z --overwrite
    r.resamp.rst input=temp ew_res=25 ns_res=25 elev=final_coverage --overwrite
    r.colors -n map=final_coverage color=elevation
else
    echo -e "Usage:\t$0 [comma-separated transmitter's section names] [INI file]"
    echo "Aggregates individual path-loss predictions for the given transmitters."
    echo
fi
