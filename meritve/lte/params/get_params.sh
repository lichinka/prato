#!/bin/bash

#for s in $( cat ../sites_lj.txt | tr -d ' ' ); do
for s in $( cat /tmp/sites.txt | tr -d ' ' ); do
    CELL="$( echo "${s}" | cut -d'|' -f1 | tr -d ' ' )*"
    X="$( echo "${s}" | cut -d'|' -f2 | tr -d ' ' )"
    Y="$( echo "${s}" | cut -d'|' -f3 | tr -d ' ' )"
    TGT="$( cat ./url.txt| sed -e "s/@@/${CELL}/g" )"
    TGT="$( echo "${TGT}" | sed -e "s/#x#/${X}/g" )"
    TGT="$( echo "${TGT}" | sed -e "s/#y#/${Y}/g" )"
    wget --no-check-certificate -O ../${CELL}.params --private-key=./key.pem --certificate=./cert.pem --ca-certificate=./ca_cert.pem ${TGT}
done
