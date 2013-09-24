#include "worker/coverage.h"



/**
 * Deactivates the OpenCL environment, releasing all associated memory.
 *
 * tx_params    a structure holding transmitter-specific configuration
 *              parameters;
 *
 */
void
close_gpu (Tx_parameters *tx_params);

/**
 * Deactivates the OpenCL environment, releasing all associated memory.
 *
 * tx_params    a structure holding transmitter-specific configuration
 *              parameters;
 *
 */
void
release_gpu (Tx_parameters *tx_params);

/**
 * Initializes the OpenCL environment that enabled calculation using
 * GPU hardware on the workers, if available.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 * device_hint      sends a hint to the OpenCL backend about which device
 *                  id to select and use; useful for using multiple GPUs.-
 *
 */ 
void
init_gpu (Parameters    *params,
          Tx_parameters *tx_params,
          const int     device_hint);

/**
 * Defines a 2D execution range for a kernel, consisting of square tiles,
 * which size is based on the execution capabilities of the available GPU 
 * hardware, i.e. the number of concurrent threads it can handle.
 *
 * params   a structure holding configuration parameters which are 
 *          common to all transmitters;
 * global_sizes the two-dimensional global execution range (output parameter);
 * local_sizes  the two-dimensional local execution range (output parameter).-
 *
 */
void 
define_2D_range (const Parameters *params,
                 size_t *global_sizes,
                 size_t *local_sizes);
