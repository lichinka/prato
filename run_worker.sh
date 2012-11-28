#!/bin/bash

RANK=$1
CMD="./worker/worker"
PSQL_SERVER=localhost

if [ -n "${RANK}" ]; then
    #
    # redirect worker's output to the database server
    #
    ${CMD} | psql -q -h ${PSQL_SERVER} -U garufa grass_backend 
else
    echo "Usage"
    echo "  $0 [worker's process ID or rank]"
    echo
fi

