#include "worker/coverage.h"
#include "performance/metric.h"



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
    int i;
    MPI_Status status;

    //
    // synch point: input data to all workers
    //
    MPI_Barrier (*comm);

    //
    // receive the broadcasted common parameters structure
    //
    Parameters *params = (Parameters *) malloc (sizeof (Parameters));

    MPI_Bcast (params,
               sizeof (Parameters),
               MPI_BYTE,
               _COVERAGE_MASTER_RANK_,
               *comm);
    //
    // the worker doesn't use the INI file content
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
    //
    // receive the transmitter-specific parameters
    //
    params->tx_params = (Tx_parameters *) malloc (sizeof (Tx_parameters));
    MPI_Recv (params->tx_params,
              sizeof (Tx_parameters),
              MPI_BYTE,
              _COVERAGE_MASTER_RANK_,
              MPI_ANY_TAG,
              *comm,
              &status);
    //
    // sync point: data passing finished, starting coverage processing
    //
    MPI_Barrier (*comm);

    //
    // start processing coverage prediction for this transmitter
    //
    double ericsson_params [4] = {38.0, 32.0, -12.0, 0.1};
    coverage (params,
              params->tx_params,
              ericsson_params,
              4);
    output_to_stdout (params);

    //
    // sync point: coverage finished
    //
    MPI_Barrier (*comm);

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

