#include "worker/optimize.h"




/**
 * Receives data distributed by the master process, that is common
 * to all transmitters.
 *
 * params   a structure holding configuration parameters which are 
 *          common to all transmitters;
 * comm     the MPI communicator to use.-
 *
 */
static void 
receive_common_data (Parameters *params,
                     MPI_Comm comm)
{
    //
    // receive the broadcasted common parameters structure
    //
    MPI_Bcast (params,
               sizeof (Parameters),
               MPI_BYTE,
               _COVERAGE_MASTER_RANK_,
               comm);
    //
    // the worker doesn't use these matrices;
    // it uses a reduced version of them, which are received from master
    //
    params->m_dem  = NULL;
    params->m_clut = NULL;
    params->m_loss = NULL;

    //
    // the worker doesn't parse the INI file content, it
    // receives the relevant parts from master
    //
    params->ini_file_content = NULL;
    params->ini_file_content_size = -1;

    //
    // check flag indicating whether to run in optimization mode
    //
    if (params->use_opt)
        fprintf (stdout, 
                 "*** INFO: Optimization mode enabled\n");
    else
        fprintf (stdout, 
                 "*** INFO: Coverage prediction mode enabled\n");
}



/**
 * Receives data distributed by the master process, that is specific to
 * one transmitter only.
 *
 * params   a structure holding configuration parameters which are 
 *          common to all transmitters;
 * comm     the MPI communicator to use.-
 *
 */
static void 
receive_tx_data (Parameters *params,
                 MPI_Comm comm)
{
    Tx_parameters *tx_params;
    MPI_Status     status;

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
              comm,
              &status);

    if (status.MPI_ERROR)
        fprintf (stderr, 
                 "*** ERROR: Transmitter parameters incorrectly received\n");
    else
        tx_params = params->tx_params;

    //
    // allocate memory for the transmitter matrices;
    // they contain strictly the data needed for the calculation inside 
    // the user-defined calculation radius (params->radius)
    //
    int r;

    //
    // digital elevation model
    //
    double *m_dem_data = (double *) calloc (tx_params->nrows * tx_params->ncols, 
                                            sizeof (double));
    tx_params->m_dem = (double **) calloc (tx_params->nrows,
                                           sizeof (double *));
    for (r = 0; r < tx_params->nrows; r ++)
        tx_params->m_dem[r] = &(m_dem_data[r * tx_params->ncols]);
    //
    // receive DEM data from master
    //
    MPI_Recv (tx_params->m_dem[0],
              tx_params->nrows * tx_params->ncols,
              MPI_DOUBLE,
              _COVERAGE_MASTER_RANK_,
              MPI_ANY_TAG,
              comm,
              &status);
    if (status.MPI_ERROR)
        fprintf (stderr, 
                 "*** ERROR: Incorrect receive of DEM data\n");

    //
    // clutter data
    //
    double *m_clut_data = (double *) calloc (tx_params->nrows * tx_params->ncols,
                                             sizeof (double));
    tx_params->m_clut = (double **) calloc (tx_params->nrows,
                                            sizeof (double *));
    for (r = 0; r < tx_params->nrows; r ++)
        tx_params->m_clut[r] = &(m_clut_data[r * tx_params->ncols]);
    //
    // receive clutter data from master
    //
    MPI_Recv (tx_params->m_clut[0],
              tx_params->nrows * tx_params->ncols,
              MPI_DOUBLE,
              _COVERAGE_MASTER_RANK_,
              MPI_ANY_TAG,
              comm,
              &status);
    if (status.MPI_ERROR)
        fprintf (stderr, 
                 "*** ERROR: Incorrect receive of clutter data\n");

    //
    // radio zones - only allocation
    //
    char *m_radio_zone_data = (char *) calloc (tx_params->nrows * tx_params->ncols, 
                                               sizeof (char));
    tx_params->m_radio_zone = (char **) calloc (tx_params->nrows,
                                                sizeof (char *));
    for (r = 0; r < tx_params->nrows; r ++)
        tx_params->m_radio_zone[r] = &(m_radio_zone_data[r * tx_params->ncols]);

    //
    // antenna losses - only allocation
    //
    double *m_antenna_loss_data = (double *) calloc (tx_params->nrows * tx_params->ncols, 
                                                     sizeof (double));
    tx_params->m_antenna_loss = (double **) calloc (tx_params->nrows,
                                                    sizeof (double *));
    for (r = 0; r < tx_params->nrows; r ++)
        tx_params->m_antenna_loss[r] = &(m_antenna_loss_data[r * tx_params->ncols]);

    //
    // path loss - only allocation
    //
    double *m_loss_data = (double *) calloc (tx_params->nrows * tx_params->ncols, 
                                             sizeof (double));
    tx_params->m_loss = (double **) calloc (tx_params->nrows,
                                            sizeof (double *));
    for (r = 0; r < tx_params->nrows; r ++)
        tx_params->m_loss[r] = &(m_loss_data[r * tx_params->ncols]);

    //
    // antenna losses, radio zones and field-measurements matrices 
    // are only used in optimization mode
    //
    if (params->use_opt)
    {
        //
        // field measurements
        //
        double *m_field_meas_data = (double *) calloc (tx_params->nrows * tx_params->ncols,
                                                       sizeof (double));
        tx_params->m_field_meas = (double **) calloc (tx_params->nrows, 
                                                      sizeof (double *));
        for (r = 0; r < tx_params->nrows; r ++)
            tx_params->m_field_meas[r] = &(m_field_meas_data[r * tx_params->ncols]);
        //
        // receive field-measurement data from master
        //
        MPI_Recv (tx_params->m_field_meas[0],
                  tx_params->nrows * tx_params->ncols,
                  MPI_DOUBLE,
                  _COVERAGE_MASTER_RANK_,
                  MPI_ANY_TAG,
                  comm,
                  &status);
        if (status.MPI_ERROR)
            fprintf (stderr, 
                     "*** ERROR: Incorrect receive of field measurements data\n");
    }
}



/**
 * Starts a worker process.
 *
 * rank The rank of this worker process;
 * comm the MPI communicator to use.-
 *
 */
void worker (const int rank,
             MPI_Comm comm)
{
    int has_finished;
    MPI_Status status;
    double ericsson_params [4] = {38.0, 32.0, -12.0, 0.1};
    Parameters *params = (Parameters *) malloc (sizeof (Parameters));

    //
    // sync point: receive common input data from master
    //
    MPI_Barrier (comm);
    receive_common_data (params, comm);
    params->tx_params = NULL;

    //
    // sync point: common data distribution finished 
    //
    MPI_Barrier (comm);

    //
    // start coverage-processing loop
    //
    has_finished = 0;
    while (!has_finished)
    {
        //
        // inform master we are ready to receive work
        //
        MPI_Send (NULL,
                  0,
                  MPI_BYTE,
                  _COVERAGE_MASTER_RANK_,
                  _WORKER_IS_IDLE_TAG_,
                  comm);
        //
        // check if this worker should shutdown
        // or continue working
        //
        MPI_Recv (NULL,
                  0,
                  MPI_BYTE,
                  _COVERAGE_MASTER_RANK_,
                  MPI_ANY_TAG,
                  comm,
                  &status);
        if (status.MPI_TAG == _WORKER_SHUTDOWN_TAG_)
            has_finished = 1;
        else
        {
            //
            // receive input data for the next transmitter
            //
            receive_tx_data (params, comm);

            //
            // calculate coverage prediction or optimize parameters?
            // 
            if (params->use_opt)
            {
                optimize (params,
                          params->tx_params);
            }
            else
            {
                //
                // calculate coverage for the received transmitter
                //
                coverage (params,
                          params->tx_params,
                          ericsson_params,
                          4);
                //
                // start result dump
                //
                output_to_stdout (params,
                                  params->tx_params);
            }
        }
    }

    //
    // deallocate memory before exiting
    //
    if (params->use_opt)
    {
        free (&(params->tx_params->m_field_meas[0][0]));
        free (params->tx_params->m_field_meas);
    }
    free (&(params->tx_params->m_antenna_loss[0][0]));
    free (params->tx_params->m_antenna_loss);
    free (&(params->tx_params->m_radio_zone[0][0]));
    free (params->tx_params->m_radio_zone);
    free (&(params->tx_params->m_loss[0][0]));
    free (params->tx_params->m_loss);
    free (&(params->tx_params->m_dem[0][0]));
    free (params->tx_params->m_dem);
    free (&(params->tx_params->m_clut[0][0]));
    free (params->tx_params->m_clut);
    free (params->tx_params);
    free (params);
}

