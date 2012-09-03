#ifndef _PERFORMANCE_MEASURE_H_
#define _PERFORMANCE_MEASURE_H_

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include <papi.h>


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
static void debug_mpi ( )
{
    int i = 0;
    char hostname[256];
    gethostname(hostname, sizeof(hostname));
    printf("PID %d on %s ready for attach\n", getpid(), hostname);
    fflush(stdout);
    while (0 == i)
        sleep(5);
}

#endif
