#!/bin/bash


for i in $(seq 0 255); do 
	echo -n "${i}     "
	SUM=$(grep 'result' log/weak_scaling_256.txt | cut -d '(' -f2 | grep "^${i})" - | cut -d')' -f2 | sed -e 's/sec/+/g' | tr -d '\n' | sed -e 's/+$/\n/g' | tr -d '[:blank:]' | bc)
	COUNT=$(grep 'result' log/weak_scaling_256.txt | cut -d '(' -f2 | grep "^${i})" - | wc -l)
	AVG=$(echo -e "scale=10\n ${SUM}/${COUNT}" | bc)
	echo "${SUM}	${COUNT}	${AVG}"
done
