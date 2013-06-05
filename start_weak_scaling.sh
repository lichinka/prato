#!/bin/bash

NP="1 2 4 8 16 32 64 128"
NTX="5 10 20 40 80"

echo "*** WARNING: This script will overwrite the log files for weak-scaling experiments"
read

for np in ${NP}; do
    for ntx in ${NTX}; do
        for i in $( seq -w 001 020 ); do
            ntx="$( echo "${np}*${ntx}" | bc -l )"
            ./run_weak_scaling.sh ${np} - $( ./run_rand.sh ${ntx} ) -m | gzip -c - > log_perf/no_db/weak-NTX_${ntx}-NP_${np}.${i}
        done
	done
done

