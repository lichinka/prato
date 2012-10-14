#ifndef _PERFORMANCE_MEASURE_H_
#define _PERFORMANCE_MEASURE_H_

#define _METRIC_CHAR_BUFF_SIZE_    2048

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>

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
 * Measures the time passed between calls.
 *
 * message          A message to be printed when the measurement starts;
 *                  passing NULL to this parameter makes the measurement
 *                  stop and display the results.-
 *
 */
void measure_time (const char *message);

/**
 * Measures the time passed between calls.
 *
 * message          A message to be printed when the measurement starts;
 *                  passing NULL to this parameter marks the end of the 
 *                  measurement stop and displays the elapsed time;
 * clock_id         an ID to uniquely identify this clock, so that different
 *                  things can be measured at the same time.-
 *
 */
void measure_time_id (const char *message, const unsigned int clock_id);

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
