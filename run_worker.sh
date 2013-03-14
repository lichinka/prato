#!/bin/bash

REDIR=$1
#CMD="valgrind --leak-check=yes ./worker/worker"
CMD="./worker/worker"
#PSQL_SERVER=ninestein
PSQL_SERVER=localhost
#PSQL_USER=grassuser
PSQL_USER=garufa
PSQL_DB=grass

if [ "${REDIR}" = "-db" ]; then
    #
    # redirect worker's output to the database server
    #
    ${CMD} | grep -v 'GPU' - | grep -v 'INFO' - | psql -q -h ${PSQL_SERVER} -U ${PSQL_USER} -d ${PSQL_DB}
else 
    if [ "${REDIR}" = "-rast" ]; then
        ${CMD} | grep '^[0-9]' - | v.in.ascii -t output=temp format=point -z z=3 --overwrite 
        v.to.rast input=temp type=point output=temp use=z --overwrite
    else
        if [ -z "${REDIR}" ]; then
            ${CMD} > /tmp/worker.log 2>&1
        else
            echo "Usage: $0 [-db] [-rast]"
            echo "Starts a worker process, writing its results to stdout."
            echo
            echo " -db		redirects output to the database server, stdout otherwise"
            echo " -rast	redirects output to a raster map, named 'temp'"
        fi
    fi
fi

