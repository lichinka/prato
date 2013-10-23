#!/bin/bash

INPUT=$1
BS_X=$2
BS_Y=$3
ROAD_CLUT=$4
OUTPUT=$5

if [ -n "${INPUT}" ] && [ -n "${BS_X}" ] && [ -n "${BS_Y}" ] && [ -n "${ROAD_CLUT}" ] && [ -n "${OUTPUT}" ]; then
    echo "Moving coordinate origin to (${BS_X}, ${BS_Y}) ..."
    r.mapcalc "x=x()-${BS_X}"
    r.mapcalc "y=y()-${BS_Y}"
    echo "Calculating rotation angles (clock-wise) for all four quadrants (0 degrees point to the north) ..."
    r.mapcalc "angle=if(x*y>0, if(x>0, atan(x/y)+180, atan(x/y)), if(x>0, 360+atan(x/y), 180+atan(x/y)))"
    echo "Calculating rotation angles for points with partial coordinates equal to zero, i.e. (x, 0) and (0, y) ..."
    r.mapcalc "angle=if(isnull(angle), if(x>0, 270, 90), angle)"
    r.mapcalc "angle=if(angle==180 && y<0, float (0), float (angle))"
    echo "Replacing road clutter with the adjacent one, towards the transmitter ..."
    r.mapcalc "${OUTPUT}=if(angle>=0.0 && angle<22.5, ${INPUT}[-1,0],\
if(angle>=22.5 && angle<67.5,   ${INPUT}[-1,1],\
if(angle>=67.5 && angle<112.5,  ${INPUT}[0,1],\
if(angle>=112.5 && angle<157.5, ${INPUT}[1,1],\
if(angle>=157.5 && angle<202.5, ${INPUT}[1,0],\
if(angle>=202.5 && angle<247.5, ${INPUT}[1,-1],\
if(angle>=247.5 && angle<292.5, ${INPUT}[0,-1],\
if(angle>=292.5 && angle<337.5, ${INPUT}[-1,-1],\
if(angle>=337.5 && angle<360.0, ${INPUT}[-1,0],\
null())))))))))"
    echo "Cleaning up temporary raster files ..."
    g.remove rast=x
    g.remove rast=y
    g.remove rast=angle
    echo "done"
else
    echo "Replaces the road clutter with the nearest neighboring one, found in the direction towards the transmitter."
    echo "Usage:"
    echo "	$0 [input raster name] [transmitter X coordinate] [transmitter Y coordinate] [road-clutter category number] [output raster name]"
    exit -1
fi

