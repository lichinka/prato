#!/bin/bash

NTX=$1
NP=$2
NP_NODE=$3
HOSTFILE=hostfile.local

if [ -n "${NTX}" ] && [ -n "${NP}" ] && [ -n "${NP_NODE}" ]; then
	echo -e "localhost	slots=1\n" > ${HOSTFILE}
	NODES=$( echo "${NP}/${NP_NODE}" | bc )
	for h in $(grep -v '^#' ~/hosts.20121003 | head -n ${NODES}); 
	do 
		echo "${h}        slots=${NP_NODE}" >> ${HOSTFILE}
	done
	NTX_NP=$( echo "${NTX}/${NP}" | bc )
	TOTAL_NTX=$( echo "${NP} * ${NTX_NP}" | bc )
	echo "Number of processors:	${NP}"
	echo "Transmitters per processor:	${NTX_NP}"
	echo "Processes per node:	${NP_NODE}"
	$(./run_many.sh ${TOTAL_NTX})
else
	echo "Usage:"
	echo "	$0 [number of transmitters] [number of processors] [number of processes per node]"
fi
