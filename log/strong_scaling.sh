#!/bin/bash

HFILE=$1
NTX=$2
NP=$3
NP_NODE=$4

if [ -n "${HFILE}" ] && [ -n "${NTX}" ] && [ -n "${NP}" ] && [ -n "${NP_NODE}" ]; then
	echo -e "localhost	slots=1\n" > ${HFILE}
	NODES=$( echo "${NP}/${NP_NODE}" | bc )
	for h in $(grep -v '^#' ~/hosts.20121003 | head -n ${NODES}); 
	do 
		echo "${h}        slots=${NP_NODE}" >> ${HFILE}
	done
	NTX_NP=$( echo "${NTX}/${NP}" | bc )
	TOTAL_NTX=$( echo "${NP} * ${NTX_NP}" | bc )
	echo "Number of processors:	${NP}"
	echo "Transmitters per processor:	${NTX_NP}"
	echo "Processes per node:	${NP_NODE}"
	$( ./run_many.sh ${HFILE} ${TOTAL_NTX} )
else
	echo -e "Usage:\t$0 [hostfile to use] [number of transmitters] [number of processors] [number of processes per node]"
	echo "Runs the strong-scaling simulations over MPI."
fi
