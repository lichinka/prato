#include "worker/coverage.h"
#include "performance/metric.h"




/**
 * Receives data distributed by the master process, that is common
 * to all transmitters.
 *
 * params   the parameters structure into which the data should be saved;
 * comm     the MPI communicator to use.-
 *
 */
static void receive_common_data (Parameters *params,
                                 MPI_Comm *comm)
{
    int i;

    //
    // receive the broadcasted common parameters structure
    //
    MPI_Bcast (params,
               sizeof (Parameters),
               MPI_BYTE,
               _COVERAGE_MASTER_RANK_,
               *comm);
    //
    // the worker doesn't parse the INI file content, it
    // receives the relevant parts from master
    //
    params->ini_file_content = NULL;
    params->ini_file_content_size = -1;

    //
    // allocate memory for the output Path-loss array
    //
    double *m_loss_data = (double *) calloc (params->nrows * params->ncols,
                                             sizeof (double));
    params->m_loss = (double **) calloc (params->nrows,
                                        sizeof (double *));
    for (i = 0; i < params->nrows; i ++)
        params->m_loss[i] = &(m_loss_data[i * params->ncols]);

    //
    // receive the broadcasted the DEM array content
    //
    double *m_dem_data = (double *) calloc (params->nrows * params->ncols,
                                            sizeof (double));
    params->m_dem = (double **) calloc (params->nrows,
                                       sizeof (double *));
    for (i = 0; i < params->nrows; i ++)
        params->m_dem[i] = &(m_dem_data[i * params->ncols]);

    MPI_Bcast (&(params->m_dem[0][0]),
               params->nrows * params->ncols,
               MPI_DOUBLE,
               _COVERAGE_MASTER_RANK_,
               *comm);
    //
    // receive the broadcasted the Clutter array content
    //
    double *m_clut_data = (double *) calloc (params->nrows * params->ncols,
                                             sizeof (double));
    params->m_clut = (double **) calloc (params->nrows,
                                         sizeof (double *));
    for (i = 0; i < params->nrows; i ++)
        params->m_clut[i] = &(m_clut_data[i * params->ncols]);

    MPI_Bcast (&(params->m_clut[0][0]),
               params->nrows * params->ncols,
               MPI_DOUBLE,
               _COVERAGE_MASTER_RANK_,
               *comm);
}



/**
 * Receives data distributed by the master process, that is specific to
 * one transmitter only.
 *
 * params   the parameters structure into which the data should be saved;
 * comm     the MPI communicator to use.-
 *
 */
static void receive_tx_data (Parameters *params,
                             MPI_Comm *comm)
{
    MPI_Status status;

    //
    // receive the transmitter-specific parameters
    //
    if (params->tx_params == NULL)
        params->tx_params = (Tx_parameters *) malloc (sizeof (Tx_parameters));

    MPI_Recv (params->tx_params,
              sizeof (Tx_parameters),
              MPI_BYTE,
              _COVERAGE_MASTER_RANK_,
              MPI_ANY_TAG,
              *comm,
              &status);
}



/**
 * Starts a worker process.
 *
 * rank The rank of this worker process;
 * comm the MPI communicator to use.-
 *
 */
void worker (const int rank,
             MPI_Comm *comm)
{
    int has_finished = 0;
    MPI_Request request;
    double ericsson_params [4] = {38.0, 32.0, -12.0, 0.1};
    Parameters *params = (Parameters *) malloc (sizeof (Parameters));

    //
    // sync point: receive common input data from master
    //
    MPI_Barrier (*comm);
    receive_common_data (params, comm);
    params->tx_params = NULL;
    //
    // sync point: common data passing finished, 
    // starting coverage processing
    //
    MPI_Barrier (*comm);

    //
    // this worker will receive the 'shutdown' order with this message tag
    //
    MPI_Irecv (NULL,
               0,
               MPI_BYTE,
               _COVERAGE_MASTER_RANK_,
               _WORKER_SHUTDOWN_TAG_,
               *comm,
               &request);
    //
    // start processing loop
    //
    while (!has_finished)
    {
        //
        // inform master this worker is ready to receive work
        //
        MPI_Send (NULL,
                  0,
                  MPI_BYTE,
                  _COVERAGE_MASTER_RANK_,
                  _WORKER_IS_IDLE_TAG_,
                  *comm);
        //
        // receive data for processing the next transmitter
        //
        receive_tx_data (params, comm);
        //
        // calculate coverage for the received transmitter
        //
        coverage (params,
                  params->tx_params,
                  ericsson_params,
                  4);
        //
        // starting result dump
        //
        output_to_stdout (params);
        //
        // check if this worker should shutdown
        //
        MPI_Test (&request, 
                  &has_finished,
                  MPI_STATUS_IGNORE);
    }

    //
    // deallocate memory before exiting
    //
    free (params->tx_params);
    free (&(params->m_clut[0][0]));
    free (params->m_clut);
    free (&(params->m_dem[0][0]));
    free (params->m_dem);
    free (&(params->m_loss[0][0]));
    free (params->m_loss);
    free (params);
}

