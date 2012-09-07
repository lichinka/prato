#!/bin/bash

TABLE_COUNT=$1
PSQL_SERVER=localhost

if [ -n "${TABLE_COUNT}" ]; then
    TABLE_COUNT=$(( ${TABLE_COUNT} - 1 ))
    for i in $(seq 0 ${TABLE_COUNT}); do
        psql -h ${PSQL_SERVER} -U garufa grass_backend -c "DROP TABLE coverage_${i}"
    done
else
    echo "Usage"
    echo "  $0 [number of coverage tables to drop]"
    echo
fi
