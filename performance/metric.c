#include <stdio.h>
#include <stdlib.h>
#include <papi.h>

#include "metric.h"



//
// Variables to count memory bandwidth
//
static unsigned int mem_access_count = 0;

//
// Variables to count the MFlops 
// 
static float last_real_time, last_proc_time;
static long long last_flpins;



/**
 * Measures memory access bandwidth.-
 */
void memory_access (const unsigned int access_count, const size_t access_size)
{
    mem_access_count += access_count * access_size;
}



/**
 * Function to test the correctness of the code performance counters.-
 */
static void test_fail (char *file, int line, char *call, int ret_value)
{
    if (ret_value != PAPI_OK)
    {
        fprintf (stderr, "PAPI returned error %d\n", ret_value);
        exit (ret_value);
    }
}



/**
 * Measures and displays performance measures using the PAPI library.-
 */
void measure_flops (const char *func_name, const char start_measuring)
{
    int ret_value;
    float real_time, proc_time, mflops;
    long long flpins;

    //
    // collect data from the performance counters
    //
    if ((ret_value = PAPI_flops (&real_time, &proc_time, &flpins, &mflops)) < PAPI_OK)
    {
        test_fail (__FILE__, __LINE__, "measure_flops", ret_value);
    }
    else
    {
        if (start_measuring == 1)
        {
            printf ("%s\n", func_name);
            last_real_time   = real_time;
            last_proc_time   = proc_time;
            last_flpins      = flpins;
            mem_access_count = 0;
        }
        else
        {
            float mem_bandwidth = mem_access_count / 1024.0 / 1024.0;
            mem_bandwidth /= real_time - last_real_time;

            printf ("\treal time:\t%f\n\tproc time:\t%f\n", real_time - last_real_time,
                                                            proc_time - last_proc_time);
            printf ("\ttotal flpins:\t%lld\n\tMFlops:\t\t%f\n", flpins - last_flpins,
                                                                mflops);
            printf ("\tmemory (MB/s)\t%.3f\n", mem_bandwidth);
        }
    }
}

