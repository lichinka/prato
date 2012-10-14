#include "performance/metric.h"

#define _NUMBER_OF_MULTIPLE_CLOCKS_ 2048


//
// Variables to count memory bandwidth
//
static unsigned int mem_access_count = 0;

//
// Variables to count the MFlops 
// 
static float last_real_time, last_proc_time;
static long long last_flpins;

//
// variables to measure time
//
static char time_message  [_NUMBER_OF_MULTIPLE_CLOCKS_ + 1][_METRIC_CHAR_BUFF_SIZE_];
static long last_clock [_NUMBER_OF_MULTIPLE_CLOCKS_ + 1];



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
    int errno;
    float real_time, proc_time, mflops;
    long long flpins;

    //
    // collect data from the performance counters
    //
    if ((errno = PAPI_flops (&real_time, &proc_time, &flpins, &mflops)) < PAPI_OK)
    {
        test_fail (__FILE__, __LINE__, "measure_flops", errno);
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
void measure_time_id (const char *message, const unsigned int clock_id)
{
    if (clock_id <= _NUMBER_OF_MULTIPLE_CLOCKS_)
    {
        if (message != NULL)
        {
            if (clock_id == _NUMBER_OF_MULTIPLE_CLOCKS_)
                snprintf (time_message[clock_id],
                          _METRIC_CHAR_BUFF_SIZE_,
                          "%s", 
                          message);
            else
                snprintf (time_message[clock_id],
                          _METRIC_CHAR_BUFF_SIZE_,
                          "%s (%d)", 
                          message,
                          clock_id);
            last_clock[clock_id] = PAPI_get_real_usec ( );
        }
        else
        {
            //
            // ignore uninitialized calls
            //
            if (last_clock[clock_id] != 0)
            {
                double elapsed_time = PAPI_get_real_usec ( ) - last_clock[clock_id];
                double current_time = elapsed_time + last_clock[clock_id];
                elapsed_time /= 1000000.0;
		current_time /= 1000000.0;
                fprintf (stdout, "TIME:%.8f:\t%s\t%.8f sec\n", current_time,
                                                               time_message[clock_id],
                                                               elapsed_time);
                last_clock[clock_id] = 0;
            }
        }
    }
    else
    {
        fprintf (stderr, 
                 "ERROR You can have a maximum of %d concurrent time measurements\n",
                 _NUMBER_OF_MULTIPLE_CLOCKS_);
        exit (1);
    }
}



/**
 * Measures the time passed between calls.
 *
 * message          A message to be printed when the measurement starts;
 *                  passing NULL to this parameter marks the end of the 
 *                  measurement stop and displays the elapsed time;
 *
 */
void measure_time (const char *message)
{
    measure_time_id (message, 
                     _NUMBER_OF_MULTIPLE_CLOCKS_);
}

