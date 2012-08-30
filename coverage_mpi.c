#include "coverage.h"
#include "performance/metric.h"




/**
 * The master process takes care of sending the input data to the workers 
 * and gathering the results from them.
 *
 * params       A structure holding all parameters needed for 
 *              calculation;
 * nworkers     the number of workers this master has available;
 * comm         the communicator used to communicate with the workers.-
 *
 */
static void master (const Parameters *params,
                    const int nworkers,
                    MPI_Comm *comm)
{
    int ncols = params->ncols;
    int rows_per_worker = params->nrows / nworkers;

    //
    // synchronize all processes before starting
    //
    MPI_Barrier (*comm);

    //
    // broadcast the row size and the number of rows per worker
    //
    MPI_Bcast (&ncols,
               1, 
               MPI_INT, 
               MPI_ROOT,
               *comm);
    MPI_Bcast (&rows_per_worker, 
               1, 
               MPI_INT, 
               MPI_ROOT,
               *comm);
    /*
    // start processing the subarea
    //
    process_subarea (params,
                     _COVERAGE_MASTER_RANK_,
                     comm,
                     params->ncols,
                     rows_per_worker);*/
}



/**
 * Calculates the area coverage using MPI to achieve parallelization.
 *
 * argc             Number of command line parameters;
 * argv             array containing command line parameters;
 * params           a structure holding all parameters needed for 
 *                  calculation;
 * eric_params      contains the four tunning parameters for the 
 *                  Ericsson 9999 model, set by the optimization 
 *                  algorithm;
 * eric_params_len  the number of parameters within the received vector,
 *                  four in this case (A0, A1, A2 and A3);
 * output_raster    the name of the output raster created;
 *                  no output is generated if this parameter is NULL.-
 *
 */
void coverage_mpi (int argc, 
                   char *argv[],
                   const Parameters *params,
                   const double *eric_params,
                   const int eric_params_len,
                   const char *output_raster)
{
    MPI_Comm everyone;                          // intercommunicator
    int rank, world_size, universe_size, flag;
    int *universe_size_ptr;
    char worker_program [] = "worker/worker";

    MPI_Init  (&argc, &argv); 
    MPI_Comm_size (MPI_COMM_WORLD, &world_size); 
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    //
    // only one master process should be started at once
    //
    if (world_size > 1)
        G_fatal_error ("You're supposed to start only one master process at a time");

    //
    // master takes care of starting the workers and 
    // dividing the work among them
    //
    MPI_Comm_get_attr (MPI_COMM_WORLD, 
                       MPI_UNIVERSE_SIZE,  
                       &universe_size_ptr, 
                       &flag);
    universe_size = *universe_size_ptr;
    if (universe_size == 1)
        G_fatal_error ("Allocate more slots so that the master will start workers");

    //
    // Now spawn the workers. Note that there is a run-time determination 
    // of what type of worker to spawn, and presumably this calculation must 
    //  be done at run time and cannot be calculated before starting 
    // the program. If everything is known when the application is  
    // first started, it is generally better to start them all at once 
    // in a single MPI_COMM_WORLD.  
    //
    int nworkers = universe_size - 1;
    int worker_errcodes [nworkers];

    printf ("About to spawn %d worker processes ...\n", nworkers);
    MPI_Comm_spawn (worker_program, 
                    MPI_ARGV_NULL, 
                    nworkers,
                    MPI_INFO_NULL, 
                    _COVERAGE_MASTER_RANK_,
                    MPI_COMM_SELF, 
                    &everyone,  
                    worker_errcodes); 

    int i;
    for (i = 0; i < nworkers; i ++)
    {
        if (worker_errcodes[i] != 0)
            G_fatal_error ("Worker %d exited with code %d\n", i, 
                                                              worker_errcodes[i]);
    }

    //
    // the number of workers should divide the number of rows,
    // so that we will be able to evenly separate the work
    //
    if ((params->nrows % universe_size) != 0)
        G_fatal_error ("The number of processes (%d) should divide the number of rows in the input raster", 
                       universe_size);

    master (params, nworkers, &everyone);

    MPI_Finalize ( ); 
} 

