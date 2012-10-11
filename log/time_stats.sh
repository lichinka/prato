#!/bin/bash


INPUT=$1
NP=$2

if [ -n "${INPUT}" ] && [ -n "${NP}" ]; then
    for i in $(seq 0 ${NP}); do 
        echo -n "${i}     "
        SUM=$(grep 'result' ${INPUT} | cut -d '(' -f2 | grep "^${i})" - | cut -d')' -f2 | sed -e 's/sec/+/g' | tr -d '\n' | sed -e 's/+$/\n/g' | tr -d '[:blank:]' | bc)
        COUNT=$(grep 'result' ${INPUT} | cut -d '(' -f2 | grep "^${i})" - | wc -l)
        AVG=$(echo -e "scale=10\n ${SUM}/${COUNT}" | bc)
        echo "${SUM}	${COUNT}	${AVG}"
    done
else
    echo "Usage:"
    echo -e "\t $0 [dir containing log files] [number of parallel processes]"
fi
