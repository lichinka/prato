#!/bin/bash

HFILE=$1
NP=$2
NTX_NP=$3
NP_NODE=$4

if [ -n "${HFILE}" ] && [ -n "${NP}" ] && [ -n "${NTX_NP}" ] && [ -n "${NP_NODE}" ]; then
	echo -e "localhost	slots=1\n" > ${HFILE}
	if [ ${NP} -lt ${NP_NODE} ]; then
		NODES=1
		NP_NODE=${NP}
	else
		NODES=$(echo "${NP}/${NP_NODE}" | bc)
	fi
	for h in $(grep -v '^#' ~/hosts.20121003 | head -n ${NODES}); 
	do 
		echo "${h}        slots=${NP_NODE}" >> ${HFILE}
	done
	TOTAL_NTX=$( echo "${NP} * ${NTX_NP}" | bc )
	echo "Number of processors:	${NP}"
	echo "Transmitters per processor:	${NTX_NP}"
	echo "Processes per node:	${NP_NODE}"
	$( ./run_many.sh ${HFILE} ${TOTAL_NTX} )
else
	echo "Usage:"
	echo "	$0 [hostfile to use] [number of processors] [number of transmitters per processor] [number of processes per node]"
fi
