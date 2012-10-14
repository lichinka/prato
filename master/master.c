#include "worker/coverage.h"
#include "performance/metric.h"



/**
 * Distribute common input data among all workers.
 *
 * params       A structure holding all parameters needed for calculation;
 * comm         the communicator used to communicate with the workers.-
 *
 */
static void distribute_common_data (Parameters *params,
                                    MPI_Comm *comm)
{
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
}



/**
 * Sends transmitter input data to the specified worker.
 *
 * tx_params    A structure holding transmitter-specific parameters;
 * comm         the communicator used to communicate with the worker;
 * worker_rank  the target worker.-
 *
 */
static void send_tx_data (Tx_parameters *tx_params,
                          MPI_Comm *comm,
                          int worker_rank)
{
    //
    // send the transmitter-specific parameters to each of the workers
    //
    MPI_Send (tx_params,
              sizeof (Tx_parameters),
              MPI_BYTE,
              worker_rank,
              1,
              *comm);
}



/**
 * Dynamically spawns worker processes to calculate the area coverage using MPI.
 * It returns the number of spawned workers, and a reference to the communicator
 * used to talk to them in the 'worker_comm' parameter.
 *
 * params           a structure holding all parameters needed for calculation;
 * worker_comm      output parameter: the communicator used to talk to the 
 *                  spawned workers.-
 *
 */
static int spawn_workers (Parameters *params,
                          MPI_Comm *worker_comm)
{
    int i;
    int rank, world_size, universe_size, flag;
    int *universe_size_ptr;

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
    }
    else if (nworkers > params->ntx)
    {
        nworkers = params->ntx;
    }
    //
    // prepare workers' commands and arguments
    //
    int buff_size = _CHAR_BUFFER_SIZE_;
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
    printf ("Spawning %d workers to process %d transmitters ...\n", nworkers,
                                                                    params->ntx);
    MPI_Comm_spawn_multiple (nworkers,
                             worker_program,
                             &(worker_argv[0]),
                             worker_maxproc,
                             worker_info,
                             _COVERAGE_MASTER_RANK_,
                             MPI_COMM_SELF,
                             worker_comm,
                             worker_errcodes);
    for (i = 0; i < nworkers; i ++)
        if (worker_errcodes[i] != 0)
            fprintf (stderr, "WARNING Worker %d. returned with exit code %d\n", i,
                                                                                worker_errcodes[i]);

    //
    // deallocate memory
    //
    for (i = 0; i < nworkers; i ++)
    {
        free (worker_arg[i]);
        free (worker_program[i]);
    }

    return nworkers;
} 



/**
 * Starts coverage calculation over MPI.
 *
 * argc             Number of command line parameters;
 * argv             array containing command line parameters;
 * params           a structure holding all parameters needed for calculation.-
 *
 */
void coverage_mpi (int argc, 
                   char *argv [],
                   Parameters *params)
{
    int nworkers = -1;
    char buff [_CHAR_BUFFER_SIZE_];
    MPI_Status status;

    measure_time ("Dynamic process spawning");
    MPI_Init  (&argc, &argv); 
    MPI_Comm worker_comm;

    nworkers = spawn_workers (params, &worker_comm);
    if (nworkers < 1)
    {
        fprintf (stderr, "ERROR Could not start any worker processes\n");
        exit (1);
    }
    measure_time (NULL);

    //
    // sync point: pass common input data to all workers
    //
    MPI_Barrier (worker_comm);
    measure_time ("Common data distribution");
    distribute_common_data (params,
                            &worker_comm);
    //
    // sync point: common data distribution finished
    //
    MPI_Barrier (worker_comm);
    measure_time (NULL);

    //
    // start processing loop
    //
    int worker_rank;
    int tx_count = params->ntx;
    int running_workers = nworkers;
    while (running_workers > 0)
    {   
        MPI_Recv (&buff, 
                  0,
                  MPI_BYTE,
                  MPI_ANY_SOURCE,
                  MPI_ANY_TAG,
                  worker_comm,
                  &status);
        worker_rank = status.MPI_SOURCE;
        //
        // coverage calculation and result dump finished
        //
        measure_time_id (NULL, worker_rank);

        switch (status.MPI_TAG)
        {   
            case (_WORKER_IS_IDLE_TAG_):
                //
                // send calculation data, if we still have transmitters
                //
                if (tx_count > 0)
                {
                    //
                    // tell the worker to keep working
                    //
                    MPI_Send (NULL,
                              0,
                              MPI_BYTE,
                              worker_rank,
                              _WORKER_KEEP_WORKING_TAG_,
                              worker_comm);
                    //
                    // starting transmitter-data send
                    //
                    measure_time_id ("Transmitter data send", 
                                     worker_rank);
                    send_tx_data (&(params->tx_params[-- tx_count]),
                                  &worker_comm,
                                  worker_rank);
                    //
                    // transmitter-data send finished,
                    // start coverage calculation
                    // 
                    measure_time_id (NULL, worker_rank);
                    measure_time_id ("Calculation and result dump",
                                     worker_rank);
                }
                else
                {
                    //
                    // send the shutdown order to this worker
                    //
                    MPI_Send (NULL,
                              0,
                              MPI_BYTE,
                              worker_rank,
                              _WORKER_SHUTDOWN_TAG_,
                              worker_comm);
                    //
                    // coverage calculation and result dump finished
                    //
                    measure_time_id (NULL, worker_rank);
                    running_workers --;
                }
                break;

            default:
                fprintf (stderr, 
                         "WARNING Unknown message from %d. worker\n", 
                         worker_rank);
        }   
    }

    MPI_Finalize ( );
}

