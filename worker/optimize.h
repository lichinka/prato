#include <float.h>

#include "worker/coverage.h"


/**
 * The random number generator defined by the URAND macro should return
 * double-precision floating-point values, uniformly distributed over the 
 * interval [0.0, 1.0)
 *
 */
#define URAND  ((double)rand ( ) / ((double)RAND_MAX + 1.0))


/**
 * Optimizes the parameters of the prediction model to fit field measurements
 * within different radio zones.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration 
 *              parameters;
 *
 */
void 
optimize_on_worker (Parameters    *params,
                    Tx_parameters *tx_params);

/**
 * Optimizes the parameters of the prediction model to fit field measurements
 * within different radio zones using a deterministic (analytical) approach.
 * It also calculates the objective-function value using the solutions 
 * received from the master process, where the optimization of the 
 * clutter-category losses (using differential evolution (DE)) runs.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration 
 *              parameters;
 * comm         the communicator used to receive messages from the master;
 *
 */
void 
optimize_from_master (Parameters    *params,
                      Tx_parameters *tx_params,
                      MPI_Comm      *comm);
