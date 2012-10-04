#!/bin/bash

NP=$1
NTX_NP=$2
NP_NODE=$3
HOSTFILE=hostfile.local

if [ -n "${NP}" ] && [ -n "${NTX_NP}" ] && [ -n "${NP_NODE}" ]; then
	echo -e "localhost	slots=1\n" > ${HOSTFILE}
	NODES=$(echo "${NP}/${NP_NODE}" | bc)
	for h in $(grep -v '^#' ~/hosts.20121003 | head -n ${NODES}); 
	do 
		echo "${h}        slots=${NP_NODE}" >> ${HOSTFILE}
	done
	TOTAL_NTX=$( echo "${NP} * ${NTX_NP}" | bc )
	echo "Number of processors:	${NP}"
	echo "Transmitters per processor:	${NTX_NP}"
	echo "Processes per node:	${NP_NODE}"
	$(./run_many.sh ${TOTAL_NTX})
else
	echo "Usage:"
	echo "	$0 [number of processors] [number of transmitters per processor] [number of processes per node]"
fi
