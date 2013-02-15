#include "worker/optimize.h"




/**
 * Allocates a 2D matrix of the specified dimensions as a continous
 * chunck of memory. Despite this, the matrix can be referenced with
 * two indices as `matrix[i][j]`.
 * This function returns the target pointer to where the memory has 
 * been allocated.
 *
 * nrows    the number of rows of the matrix;
 * ncols    the number of columns in the matrix;
 * m_ptr    target pointer, where the address of the allocated memory
 *          is saved (output parameter).-
 *
 */
static double **
allocate_double_matrix (const int    nrows,
                        const int    ncols,
                        double     **m_ptr)
{
    int r;

    //
    // only allocate new memory if the target pointer is NULL
    //
    if (m_ptr == NULL)
     {
        double *m_ptr_data = (double *) calloc (nrows * ncols, 
                                                sizeof (double));
        m_ptr = (double **) calloc (nrows,
                                    sizeof (double *));
        for (r = 0; r < nrows; r ++)
            m_ptr[r] = &(m_ptr_data[r * ncols]);
    }
    return m_ptr;
}



/**
 * Allocates a 2D matrix of the specified dimensions as a continous
 * chunck of memory. Despite this, the matrix can be referenced with
 * two indices as `matrix[i][j]`.
 * This function returns the target pointer to where the memory has 
 * been allocated.
 *
 * nrows    the number of rows of the matrix;
 * ncols    the number of columns in the matrix;
 * m_ptr    target pointer, where the address of the allocated memory
 *          is saved (output parameter).-
 *
 */
static char **
allocate_char_matrix (const int    nrows,
                      const int    ncols,
                      char       **m_ptr)
{
    int r;

    //
    // only allocate new memory if the target pointer is NULL
    //
    if (m_ptr == NULL)
    {
        char *m_ptr_data = (char *) calloc (nrows * ncols, 
                                            sizeof (char));
        m_ptr = (char **) calloc (nrows,
                                  sizeof (char *));
        for (r = 0; r < nrows; r ++)
            m_ptr[r] = &(m_ptr_data[r * ncols]);
    }
    return m_ptr;
}



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
    // it uses a reduced version of them, which relevant data are 
    // received from master; 
    // they are kept in the `Tx_parameters` structure
    //
    params->m_dem     = NULL;
    params->m_clut    = NULL;
    params->m_loss    = NULL;

    //
    // mark the `Tx_parameters` structure as uninitialized
    //
    params->tx_params = NULL;

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
    //
    // check flag indicating whether to use the GPU if available
    //
    if (params->use_gpu)
        fprintf (stdout, 
                 "*** INFO: GPU hardware will be used on the workers, if available\n");
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
    char          uninitialized_pointers;
    Tx_parameters *tx_params;
    MPI_Status     status;

    //
    // allocate the `Tx_parameters` structure, if we haven't already,
    // and mark all its pointers as uninitialized
    //
    if (params->tx_params == NULL)
    {
        params->tx_params = (Tx_parameters *) malloc (sizeof (Tx_parameters));
        uninitialized_pointers = 1;
    }
    else
        uninitialized_pointers = 0;

    //
    // receive the transmitter-specific parameters
    //
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
    {
        //
        // an alias for less writing
        //
        tx_params = params->tx_params;

        //
        // after receiving the `Tx_parameters` structure, all contained
        // pointer are invalid, because they contain addresses mapped 
        // within the master process; make sure we clean this up
        //
        if (uninitialized_pointers)
        {
            tx_params->eric_params[0]     = 38.0;
            tx_params->eric_params[1]     = 32.0;
            tx_params->eric_params[2]     = -12.0;
            tx_params->eric_params[3]     = 0.1;
            tx_params->diagram            = NULL;
            tx_params->m_dem              = NULL;
            tx_params->m_dem_dev          = NULL;
            tx_params->m_clut             = NULL;
            tx_params->m_clut_dev         = NULL;
            tx_params->m_field_meas       = NULL;
            tx_params->m_field_meas_dev   = NULL;
            tx_params->m_loss             = NULL;
            tx_params->m_loss_dev         = NULL;
            tx_params->m_antenna_loss     = NULL;
            tx_params->m_antenna_loss_dev = NULL;
            tx_params->m_radio_zone       = NULL;
            tx_params->m_radio_zone_dev   = NULL;
            tx_params->m_obst_height      = NULL;
            tx_params->m_obst_dist        = NULL;
            tx_params->m_obst_offset      = NULL;
            tx_params->ocl_obj            = NULL;
        }
    }

    //
    // allocate memory for the transmitter matrices;
    // they contain strictly the data needed for calculation, i.e. 
    // within the user-defined calculation radius (params->radius)
    //

    //
    // digital elevation model
    //
    tx_params->m_dem = allocate_double_matrix (tx_params->nrows,
                                               tx_params->ncols,
                                               tx_params->m_dem);
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
    // transmitter height above sea level
    //
    tx_params->total_tx_height  = tx_params->m_dem[tx_params->tx_east_coord_idx]
                                                  [tx_params->tx_north_coord_idx];
    tx_params->total_tx_height += tx_params->antenna_height_AGL;

    //
    // clutter data
    //
    tx_params->m_clut = allocate_double_matrix (tx_params->nrows,
                                                tx_params->ncols,
                                                tx_params->m_clut);
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
    // radio zones 
    //
    tx_params->m_radio_zone = allocate_char_matrix (tx_params->nrows,
                                                    tx_params->ncols,
                                                    tx_params->m_radio_zone);
    //
    // antenna-introduced losses
    //
    tx_params->m_antenna_loss = allocate_double_matrix (tx_params->nrows,
                                                        tx_params->ncols,
                                                        tx_params->m_antenna_loss);
    //
    // path loss
    //
    tx_params->m_loss = allocate_double_matrix (tx_params->nrows,
                                                tx_params->ncols,
                                                tx_params->m_loss);
    //
    // heights of obstacles - line-of-sight
    //
    tx_params->m_obst_height = allocate_double_matrix (tx_params->nrows,
                                                       tx_params->ncols,
                                                       tx_params->m_obst_height);
    //
    // distances to obstacles - line-of-sight
    //
    tx_params->m_obst_dist = allocate_double_matrix (tx_params->nrows,
                                                     tx_params->ncols,
                                                     tx_params->m_obst_dist);
    //
    // offset distances to obstacles - line-of-sight
    //
    tx_params->m_obst_offset = allocate_double_matrix (tx_params->nrows,
                                                       tx_params->ncols,
                                                       tx_params->m_obst_offset);
    //
    // field-measurements matrix is only used in optimization mode
    //
    if (params->use_opt)
    {
        //
        // field measurements
        //
        tx_params->m_field_meas = allocate_double_matrix (tx_params->nrows,
                                                          tx_params->ncols,
                                                          tx_params->m_field_meas);
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

    Parameters *params = (Parameters *) malloc (sizeof (Parameters));

    //
    // sync point: receive common input data from master
    //
    MPI_Barrier (comm);
    receive_common_data (params, comm);

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
                          params->tx_params);
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
        if (params->use_gpu)
        {
            free (params->tx_params->m_field_meas_dev);
            free (params->tx_params->v_partial_sum_dev);
            free (params->tx_params->v_partial_sum);
        }
    }
    if (params->use_gpu)
    {
        free (params->tx_params->ocl_obj);
        free (params->tx_params->m_dem_dev);
        free (params->tx_params->m_clut_dev);
        free (params->tx_params->m_loss_dev);
        free (params->tx_params->m_radio_zone_dev);
        free (params->tx_params->m_antenna_loss_dev);
        free (params->tx_params->m_obst_height_dev);
        free (params->tx_params->m_obst_dist_dev);
    }
    free (&(params->tx_params->m_antenna_loss[0][0]));
    free (params->tx_params->m_antenna_loss);
    free (&(params->tx_params->m_radio_zone[0][0]));
    free (params->tx_params->m_radio_zone);
    free (&(params->tx_params->m_obst_height[0][0]));
    free (params->tx_params->m_obst_height);
    free (&(params->tx_params->m_obst_dist[0][0]));
    free (params->tx_params->m_obst_dist);
    free (&(params->tx_params->m_obst_offset[0][0]));
    free (params->tx_params->m_obst_offset);
    free (&(params->tx_params->m_loss[0][0]));
    free (params->tx_params->m_loss);
    free (&(params->tx_params->m_dem[0][0]));
    free (params->tx_params->m_dem);
    free (&(params->tx_params->m_clut[0][0]));
    free (params->tx_params->m_clut);
    free (params->tx_params);
    free (params);
}

