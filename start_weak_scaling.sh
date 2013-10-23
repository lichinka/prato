#!/bin/bash

NP="64 32 16 8 4 2 1"
NTX="80 40 20 10 5"
TARGET_DIR="log_perf/dist_db"

echo "*** WARNING: This script will overwrite the log files for weak-scaling experiments in [${TARGET_DIR}]"
read

for np in ${NP}; do
    for ntx in ${NTX}; do
        ntx="$( echo "${np}*${ntx}" | bc -l )"
        for i in $( seq -w 001 020 ); do
	    OUT="weak-NTX_${ntx}-NP_${np}.${i}"
            ./run_scaling.sh ${np} -db $( ./run_rand.sh ${ntx} ) -m | grep -v '|' > /tmp/${OUT}
	    mv /tmp/${OUT} ${TARGET_DIR}/${OUT}
	    #
	    # free tmp space
            #
	    ${HOME}/remote_run.sh hostfile.ninestein 'rm -f /tmp/worker.*.log'
            #
            # wait for the cluster to cool down
            #
	    sleep 60s
        done
    done
done
