#include "worker/coverage.h"



/**
 * Initializes the OpenCL environment that enabled calculation using
 * GPU hardware on the workers, if available.
 *
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 *
 */
void
init_gpu (Tx_parameters *tx_params);
