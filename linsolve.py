#!/bin/env python3
# coding: utf-8

import scipy.sparse as sps

from numpy import genfromtxt, array, linalg
from scipy.sparse.linalg.dsolve import linsolve



#
# number of field measurements
#
field_meas_count = 242

#
# transmit power at the antenna
#
ptx = 47.0
tx_power = array ([ptx] * field_meas_count)

#
# logarithm of the distance between the base station and the field measurement
#
log_d = genfromtxt ('/tmp/worker.log',
                    delimiter='|',
                    skip_header=1,
                    usecols=(2))
#
# logarithm of the effective height difference between transmitter and receiver
#
log_HEBK = genfromtxt ('/tmp/worker.log',
                       delimiter='|',
                       skip_header=1,
                       usecols=(4))
#
# knjife-edge defraction factor
#
KDFR = genfromtxt ('/tmp/worker.log',
                   delimiter='|',
                   skip_header=1,
                   usecols=(5))
#
# - 3.2 * [log(11.75 hm)]^2 + 44.49 * log(Freq) - 4.78 * log(Freq)^2
#
f4 = genfromtxt ('/tmp/worker.log',
                 delimiter='|',
                 skip_header=1,
                 usecols=(10))
#
# loss due to clutter
#
clutter = genfromtxt ('/tmp/worker.log',
                   delimiter='|',
                   skip_header=1,
                   usecols=(11))
#
# loss due to the antenna diagram
#
antenna = genfromtxt ('/tmp/worker.log',
                  delimiter='|',
                  skip_header=1,
                  usecols=(14))
#
# received signal strength on the field
#
field_meas = genfromtxt ('/tmp/worker.log',
                         delimiter='|',
                         skip_header=1,
                         usecols=(13))
f3 = log_d * log_HEBK

#
# fill the arrays to solve the system
#
#   A x = b
#
A = array ([[field_meas_count,             log_d.sum( ),            log_HEBK.sum ( ),            f3.sum ( )],
            [   log_d.sum ( ),    (log_d*log_d).sum ( ),                  f3.sum ( ),    (log_d*f3).sum ( )],
            [log_HEBK.sum ( ),               f3.sum ( ), (log_HEBK*log_HEBK).sum ( ), (log_HEBK*f3).sum ( )],
            [      f3.sum ( ),       (log_d*f3).sum ( ),       (log_HEBK*f3).sum ( ),       (f3*f3).sum ( )]])

b = array ([(tx_power - clutter - KDFR - antenna - f4 - field_meas).sum ( ),
            (tx_power*log_d - clutter*log_d - KDFR*log_d - antenna*log_d - f4*log_d - field_meas*log_d).sum ( ),
            (tx_power*log_HEBK - clutter*log_HEBK - KDFR*log_HEBK - antenna*log_HEBK - f4*log_HEBK - field_meas*log_HEBK).sum ( ),
            (tx_power*f3 - clutter*f3 - KDFR*f3 - antenna*f3 - f4*f3 - field_meas*f3).sum ( )])
"""
#
# without clutter and antenna influence
#
b = array ([(tx_power - KDFR - f4 - field_meas).sum ( ),
            (tx_power*log_d - KDFR*log_d - f4*log_d - field_meas*log_d).sum ( ),
            (tx_power*log_HEBK - KDFR*log_HEBK - f4*log_HEBK - field_meas*log_HEBK).sum ( ),
            (tx_power*f3 - KDFR*f3 - f4*f3 - field_meas*f3).sum ( )])
"""

mtx_A = sps.lil_matrix (A)
mtx_A = mtx_A.tocsr ( )
x = linsolve.spsolve (mtx_A, b)
print ('Input matrix A:\n', A)
print ('Input matrix b:\n', b)
print ('Solution:\n', x)
print ('Rezidual:\n', linalg.norm (mtx_A * x - b))

#
# matrix singularity, represents the "weight" each parameter has in the
# solution vector
#
from numpy.linalg import svd
mtx_A = mtx_A.toarray ( )
U, S, Vh = svd (mtx_A)
print ('Singularity values:\n', S)

