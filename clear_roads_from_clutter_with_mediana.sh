#!/bin/bash

INPUT=$1
ROAD_CLUT=$2
OUTPUT=$3

if [ -n "${INPUT}" ] && [ -n "${ROAD_CLUT}" ] && [ -n "${OUTPUT}" ]; then
    echo "Replacing road clutter with the mode of the neighboring ones ..."
    r.mapcalc "${OUTPUT}=if(${INPUT}==${ROAD_CLUT},\
mode(${INPUT}[-1,0],\
${INPUT}[-1,1],\
${INPUT}[0,1],\
${INPUT}[1,1],\
${INPUT}[1,0],\
${INPUT}[1,-1],\
${INPUT}[0,-1],\
${INPUT}[-1,-1]),\
${INPUT})"
else
    echo "Replaces the road clutter with the mode category (i.e. most frequent) of the neighboring pixels."
    echo "Usage:"
    echo "	$0 [input raster name] [road-clutter category number] [output raster name]"
    exit -1
fi

