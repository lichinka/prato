#!/bin/bash

RANK=$1
CMD="./worker/worker"
PSQL_SERVER=ninestein
PSQL_USER=grassuser
PSQL_DB=grass

if [ -n "${RANK}" ]; then
    #
    # redirect worker's output to the database server
    #
    ${CMD} | psql -q -h ${PSQL_SERVER} -U ${PSQL_USER} -d ${PSQL_DB} 
else
    echo "Usage"
    echo "  $0 [worker's process ID or rank]"
    echo
fi

