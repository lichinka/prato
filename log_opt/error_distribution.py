#!/usr/bin/env python3
# coding: utf-8

import sys

from numpy import (genfromtxt, array, concatenate, histogram)


if __name__ == "__main__":
    if len (sys.argv) < 2:
        print ("Run 'error_distribution.sh' instead")
        sys.exit (-1)
    diff = array ([])
    for f in sys.argv[1:]:
        a = genfromtxt (f,
                        delimiter='|',
                        usecols=(2))
        diff = concatenate ( (diff, a) )
    hist, bins = histogram (diff,
                            bins=[0,5,10,15,20,25,30,35,255])
    for idx in range (len (hist)):
        print ('%d\t%d' % (bins[idx],
                           hist[idx]))
