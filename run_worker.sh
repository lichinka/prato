#!/bin/bash

REDIR=$1
RUN_DIR="${HOME}/etc/dr/tun_par/prato/src"
#CMD="valgrind --leak-check=yes ${RUN_DIR}/worker/worker"
CMD="${RUN_DIR}/worker/worker"
PSQL_SERVER=k1
PSQL_USER=grassuser
PSQL_DB=grass

#
# update LD_LIBRARY_PATH only if necessary
#
if [ -z "$( echo ${LD_LIBRARY_PATH} | grep -f "${RUN_DIR}" )" ]; then
	export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${RUN_DIR}
fi

if [ "${REDIR}" = "-db" ]; then
    #
    # redirect worker's output to the database server
    #
    #${CMD} | grep -v 'GPU' - | grep -v 'INFO' - | grep -v 'TIME' - | psql -q -h ${PSQL_SERVER} -U ${PSQL_USER} -d ${PSQL_DB}
    #${CMD} | grep '|' - | gzip -c - | ssh ${PSQL_SERVER} 'cat - > /dev/null'
    ${CMD} | gzip -c - > /dev/null
else 
    if [ "${REDIR}" = "-rast" ]; then
        ${CMD} | grep '^[0-9]' - | v.in.ascii -t output=temp format=point -z z=3 --overwrite 
        v.to.rast input=temp type=point output=temp use=z --overwrite
    else
        if [ "${REDIR}" = "-" ]; then
            rm -f /tmp/worker.log
            #ln -s /tmp/worker.$$.log /tmp/worker.log
            ${CMD} > /tmp/worker.$$.log 2>&1
        else
            echo "Usage: $0 [-] [-db] [-rast]"
            echo "Starts a worker process, writing its results to stdout."
            echo
            echo " -		redirects output to stdout"
            echo " -db		redirects output to the database server"
            echo " -rast	redirects output to a raster map, named 'temp'"
        fi
    fi
fi
