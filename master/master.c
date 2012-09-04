#include "worker/coverage.h"
#include "performance/metric.h"

/**
 * The master process takes care of sending the input data to the workers 
 * and gathering the results from them.
 *
 * params           A structure holding all parameters needed for calculation;
 * nworkers         the number of workers this master has available;
 * comm             the communicator used to communicate with the workers.-
 *
 */
static void master (Parameters *params,
                    const int nworkers,
                    MPI_Comm *comm)
{
    int worker_rank;

    /*
     * distribute_coverage_calculation_by_tx (params,
     *                                        nworkers,
     *                                        comm);
     */
    //
    // broadcast the common parameters structure
    //
    MPI_Bcast (params,
               sizeof (Parameters),
               MPI_BYTE,
               MPI_ROOT,
               *comm);
    //
    // broadcast the DEM array content
    // 
    MPI_Bcast (&(params->m_dem[0][0]),
               params->nrows * params->ncols,
               MPI_DOUBLE,
               MPI_ROOT,
               *comm);
    //
    // broadcast the Clutter array content
    // 
    MPI_Bcast (&(params->m_clut[0][0]),
               params->nrows * params->ncols,
               MPI_DOUBLE,
               MPI_ROOT,
               *comm);
    //
    // send the transmitter-specific parameters to each of the workers
    //
    for (worker_rank = 0; worker_rank < params->ntx; worker_rank ++)
    {
        MPI_Send (&(params->tx_params[worker_rank]),
                   sizeof (Tx_parameters),
                   MPI_BYTE,
                   worker_rank,
                   1,
                   *comm);
    }
}



/**
 * Calculates the area coverage using MPI to achieve parallelization.
 *
 * argc             Number of command line parameters;
 * argv             array containing command line parameters;
 * params           a structure holding all parameters needed for calculation;
 * eric_params      contains the four tunning parameters for the 
 *                  Ericsson 9999 model, set by the optimization 
 *                  algorithm;
 * eric_params_len  the number of parameters within the received vector,
 *                  four in this case (A0, A1, A2 and A3);
 *
 */
void coverage_mpi (int argc, 
                   char *argv [],
                   Parameters *params,
                   const double *eric_params,
                   const int eric_params_len)
{
    int i;
    int rank, world_size, universe_size, flag;
    int *universe_size_ptr;

    MPI_Comm everyone;
    MPI_Init  (&argc, &argv); 
    MPI_Comm_size (MPI_COMM_WORLD, &world_size); 
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);

    //
    // only one master process should be started at once
    //
    if (world_size > 1)
    {
        fprintf (stderr, "ERROR You're supposed to start only one master process at a time");
        exit (1);
    }

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
    {
        fprintf (stderr, "ERROR Allocate more slots so that the master may start some workers\n");
        exit (1);
    }

    //
    // spawn the workers. Note that there is a run-time determination 
    // of what type of worker to spawn, and presumably this calculation must 
    // be done at run time and cannot be calculated before starting 
    // the program. If everything is known when the application is  
    // first started, it is generally better to start them all at once 
    // in a single MPI_COMM_WORLD.  
    //
    int nworkers = universe_size - 1;

    if (nworkers < params->ntx)
    {
        fprintf (stderr, "WARNING There are not enough slots to process all transmitters in parallel\n");
        params->ntx = nworkers;
    }
    else if (nworkers > params->ntx)
    {
        nworkers = params->ntx;
        fprintf (stderr, "WARNING Spawning %d processes only\n", nworkers);
    }
    //
    // prepare workers' commands and arguments
    //
    int buff_size = 128;
    char    *worker_command = "run_worker.sh";
    char    *worker_program  [nworkers];
    char    *worker_arg      [nworkers];
    char    *worker_args     [nworkers][2];
    char   **worker_argv     [nworkers];
    int      worker_maxproc  [nworkers];
    MPI_Info worker_info     [nworkers];
    int      worker_errcodes [nworkers];

    for (i = 0; i < nworkers; i ++)
    {
        // worker's arguments 
        worker_arg[i] = (char *) malloc (buff_size);
        snprintf (worker_arg[i], buff_size, "%d", i);
        worker_args[i][0] = worker_arg[i];
        worker_args[i][1] = NULL;
        worker_argv[i]    = &(worker_args[i][0]);

        // worker's command
        worker_program[i] = (char *) malloc (buff_size);
        snprintf (worker_program[i], buff_size, worker_command);

        worker_maxproc[i] = 1;
        worker_info[i] = MPI_INFO_NULL;
    }
    //
    // spawn worker processes
    //
    printf ("Spawning %d. worker processes ...\n", nworkers);
    MPI_Comm_spawn_multiple (nworkers,
                             worker_program,
                             &(worker_argv[0]),
                             worker_maxproc,
                             worker_info,
                             _COVERAGE_MASTER_RANK_,
                             MPI_COMM_SELF,
                             &everyone,
                             worker_errcodes);
    //
    // start the master process
    //
    master (params, 
            nworkers, 
            &everyone);
    //
    // deallocate memory
    //
    for (i = 0; i < nworkers; i ++)
    {
        free (worker_arg[i]);
        free (worker_program[i]);
    }
    MPI_Finalize ( );
} 

