#ifndef _PRATO_WORKER_H_
#define _PRATO_WORKER_H_

#include "worker/gpu.h"
#include "worker/coverage.h"



/**
 * Initializes the transmitter parameters structure and returns a pointer
 * to the parameter `tx_params`.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 * rcv_tx_params    a structure holding transmitter-specific configuration 
 *                  parameters as received from the master process.-
 *
 */
Tx_parameters *
init_tx_params (Parameters    *params,
                Tx_parameters *tx_params,
                Tx_parameters *rcv_tx_params);

/**
 * Deallocates all internal structures contained in the 
 * transmitter-parameters structure.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters.-
 *
 */
void 
free_tx_params (Parameters    *params,
                Tx_parameters *tx_params);

/**
 *
 * Starts a worker process.
 *
 * rank The rank of this worker process;
 * comm the MPI communicator to use.-
 *
 */
void 
worker (const int rank,
        MPI_Comm comm);

#endif
