#!/bin/bash

REDIR=$1
CMD="./worker/worker"
PSQL_SERVER=ninestein
PSQL_USER=grassuser
PSQL_DB=grass

if [ "${REDIR}" = "-db" ]; then
    #
    # redirect worker's output to the database server
    #
    ${CMD} | psql -q -h ${PSQL_SERVER} -U ${PSQL_USER} -d ${PSQL_DB}
else
    if [ -z "${REDIR}" ]; then
        ${CMD} > /tmp/worker.log
        #${CMD} | grep '^[0-9]' - | v.in.ascii -t output=temp format=point -z z=3 --overwrite 
        #v.to.rast input=temp type=point output=temp use=z --overwrite
    else
        echo "Usage: $0 [-db]"
        echo "Starts a worker process."
        echo
        echo " -db	redirects output to the database server, stdout otherwise"
    fi
fi

