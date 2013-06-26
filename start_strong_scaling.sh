#!/bin/bash

NP="64 32 16 8 4 2 1"
#NTX="4096 2048 1024 512 256 128 64"
NTX="4096"

echo "*** WARNING: This script will overwrite the log files for strong-scaling experiments"
read

for ntx in ${NTX}; do
    for np in ${NP}; do
        for i in $( seq -w 001 005 ); do
            ./run_scaling.sh ${np} -db $( ./run_rand.sh ${ntx} ) -m | grep -v '|' > log_perf/with_db/strong-NTX_${ntx}-NP_${np}.${i}
        done
    done
done
