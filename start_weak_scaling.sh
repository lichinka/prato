#!/bin/bash

NP="128 64 32 16 8 4 2 1"
NTX="80 40 20 10 5"

echo "*** WARNING: This script will overwrite the log files for weak-scaling experiments"
read

for np in ${NP}; do
    for ntx in ${NTX}; do
        ntx="$( echo "${np}*${ntx}" | bc -l )"
        for i in $( seq -w 001 020 ); do
            ./run_weak_scaling.sh ${np} - $( ./run_rand.sh ${ntx} ) -m | grep -v '|' > log_perf/no_db/weak-NTX_${ntx}-NP_${np}.${i}
        done
    done
done
