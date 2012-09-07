#!/bin/bash

RANK=$1
CMD="./worker/worker"
PSQL_SERVER=192.168.1.160

if [ -n "${RANK}" ]; then
    echo ${RANK}
    #
    # dynamically create a table for the command output
    #
    CREATE_SQL=$(cat create.sql | sed -e "s/@table/coverage_${RANK}/g" -)
    echo "${CREATE_SQL}"
    psql -h ${PSQL_SERVER} -U garufa grass_backend -c "${CREATE_SQL}"
    COPY_SQL=$(cat copy.sql | sed -e "s/@table/coverage_${RANK}/g" -)
    echo "${COPY_SQL}"
    echo ${CMD} >&2 && ${CMD} | psql -h ${PSQL_SERVER} -U garufa grass_backend -c "${COPY_SQL}"
else
    echo "Usage"
    echo "  $0 [worker's process ID or rank]"
    echo
fi

