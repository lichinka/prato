#ifndef _PERFORMANCE_MEASURE_H_
#define _PERFORMANCE_MEASURE_H_

#include <stdlib.h>
#include <unistd.h>



/**
 * Measures, displays and resets performance measures.-
 */
void measure_flops (const char *func_name, const char start_measuring);

/**
 * Measures memory access bandwidth.-
 */
void memory_access (const unsigned int access_count, const size_t access_size);

/**
 * A helper function to debug MPI processes.-
 */
void debug_mpi ( );

#endif
