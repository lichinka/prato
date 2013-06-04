#!/bin/bash

NP="1 2 4 8 16"
NTX="5"

for i in $( seq -w 001 020 ); do
	for np in ${NP}; do
		ntx="$( echo "${np}*${NTX}" | bc -l )"
		./run_weak_scaling.sh ${np} - $( ./run_rand.sh ${ntx} ) -m | gzip -c - > log_perf/no_db/weak-NTX_${NTX}-NP_${np}.${i}
	done
done

