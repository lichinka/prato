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
optimize (Parameters    *params,
          Tx_parameters *tx_params);
