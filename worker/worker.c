#include "worker/worker.h"
#include "worker/optimize.h"



/**
 * Initializes the transmitter parameters structure and returns a pointer
 * to the parameter `tx_params`.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 * rcv_tx_params    a structure holding transmitter-specific configuration 
 *                  parameters as received from the master process.-
 *
 */
Tx_parameters * 
init_tx_params (Parameters    *params,
                Tx_parameters *tx_params,
                Tx_parameters *rcv_tx_params)
{
    if (tx_params == NULL)
    {
        //
        // allocate the local structure and initialize it
        //
        tx_params = (Tx_parameters *) malloc (sizeof (Tx_parameters));

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
        tx_params->m_obst_height_dev  = NULL;
        tx_params->m_obst_dist        = NULL;
        tx_params->m_obst_dist_dev    = NULL;
        tx_params->m_obst_offset      = NULL;
        tx_params->v_partial_sum      = NULL;
        tx_params->v_partial_sum_dev  = NULL;
        tx_params->v_clutter_loss_dev = NULL;
        tx_params->ocl_obj            = NULL;
    } 
    //
    // copy the received values to the local structure, keeping
    // all the pointers untouched
    //
    tx_params->beam_direction       = rcv_tx_params->beam_direction;
    tx_params->electrical_tilt      = rcv_tx_params->electrical_tilt;
    tx_params->mechanical_tilt      = rcv_tx_params->mechanical_tilt;
    tx_params->antenna_height_AGL   = rcv_tx_params->antenna_height_AGL;
    tx_params->total_tx_height      = rcv_tx_params->total_tx_height;
    tx_params->tx_east_coord        = rcv_tx_params->tx_east_coord;
    tx_params->tx_north_coord       = rcv_tx_params->tx_north_coord;
    tx_params->tx_east_coord_idx    = rcv_tx_params->tx_east_coord_idx;
    tx_params->tx_north_coord_idx   = rcv_tx_params->tx_north_coord_idx;
    tx_params->tx_power             = rcv_tx_params->tx_power;
    tx_params->nrows                = rcv_tx_params->nrows;
    tx_params->ncols                = rcv_tx_params->ncols;
    tx_params->map_north            = rcv_tx_params->map_north;
    tx_params->map_east             = rcv_tx_params->map_east;
    tx_params->map_south            = rcv_tx_params->map_south;
    tx_params->map_west             = rcv_tx_params->map_west;
    tx_params->map_north_idx        = rcv_tx_params->map_north_idx;
    tx_params->map_east_idx         = rcv_tx_params->map_east_idx;
    tx_params->map_south_idx        = rcv_tx_params->map_south_idx;
    tx_params->map_west_idx         = rcv_tx_params->map_west_idx;

    //
    // also copy the received character arrays
    //
    strncpy (tx_params->tx_name, 
             rcv_tx_params->tx_name,
             _CHAR_BUFFER_SIZE_);
    strncpy (tx_params->antenna_diagram_file,
             rcv_tx_params->antenna_diagram_file,
             _CHAR_BUFFER_SIZE_);
    strncpy (tx_params->field_meas_map,
             rcv_tx_params->field_meas_map,
             _CHAR_BUFFER_SIZE_);

    //
    // clear the previously used antenna data
    //
    if (tx_params->diagram != NULL)
    {
        free (tx_params->diagram->horizontal);
        free (tx_params->diagram->vertical);
        free (tx_params->diagram);
        tx_params->diagram = NULL;
    }
    //
    // clear the previously used GPU buffers
    //
    if (params->use_gpu)
        release_gpu (tx_params);

    //
    // allocate memory for the transmitter matrices only if needed;
    // they contain strictly the data needed for calculation, i.e. 
    // within the user-defined calculation radius (params->radius)
    //

    //
    // digital elevation model
    //
    tx_params->m_dem = prato_alloc_double_matrix (tx_params->nrows,
                                                  tx_params->ncols,
                                                  tx_params->m_dem);
    //
    // clutter data
    //
    tx_params->m_clut = prato_alloc_double_matrix (tx_params->nrows,
                                                   tx_params->ncols,
                                                   tx_params->m_clut);
    //
    // radio zones 
    //
    tx_params->m_radio_zone = prato_alloc_char_matrix (tx_params->nrows,
                                                       tx_params->ncols,
                                                       tx_params->m_radio_zone);
    //
    // antenna-introduced losses
    //
    tx_params->m_antenna_loss = prato_alloc_double_matrix (tx_params->nrows,
                                                           tx_params->ncols,
                                                           tx_params->m_antenna_loss);
    //
    // path loss
    //
    tx_params->m_loss = prato_alloc_double_matrix (tx_params->nrows,
                                                   tx_params->ncols,
                                                   tx_params->m_loss);
    //
    // heights of obstacles - line-of-sight
    //
    tx_params->m_obst_height = prato_alloc_double_matrix (tx_params->nrows,
                                                          tx_params->ncols,
                                                          tx_params->m_obst_height);
    //
    // distances to obstacles - line-of-sight
    //
    tx_params->m_obst_dist = prato_alloc_double_matrix (tx_params->nrows,
                                                        tx_params->ncols,
                                                        tx_params->m_obst_dist);
    //
    // offset distances to obstacles - line-of-sight
    //
    tx_params->m_obst_offset = prato_alloc_double_matrix (tx_params->nrows,
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
        tx_params->m_field_meas = prato_alloc_double_matrix (tx_params->nrows,
                                                             tx_params->ncols,
                                                             tx_params->m_field_meas);
    }
    return tx_params;
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
    MPI_Status     status;
    Tx_parameters *rcv_tx_params;

    //
    // the received structure will be saved here
    //
    rcv_tx_params = (Tx_parameters *) malloc (sizeof (Tx_parameters));

    //
    // receive the transmitter-specific parameters
    //
    MPI_Recv (rcv_tx_params,
              sizeof (Tx_parameters),
              MPI_BYTE,
              _COVERAGE_MASTER_RANK_,
              MPI_ANY_TAG,
              comm,
              &status);
    if (status.MPI_ERROR)
    {
        fprintf (stderr, 
                 "*** ERROR: Transmitter parameters incorrectly received\n");
        exit (-1);
    }

    //
    // initialize or clear the transmitter structure
    // before receiving the rest of the data
    //
    params->tx_params = init_tx_params (params,
                                        params->tx_params,
                                        rcv_tx_params);
    //
    // receive DEM data from master
    //
    MPI_Recv (params->tx_params->m_dem[0],
              params->tx_params->nrows * params->tx_params->ncols,
              MPI_DOUBLE,
              _COVERAGE_MASTER_RANK_,
              MPI_ANY_TAG,
              comm,
              &status);
    if (status.MPI_ERROR)
    {
        fprintf (stderr, 
                 "*** ERROR: Incorrect receive of DEM data\n");
        exit (-1);
    }
    //
    // transmitter height above sea level
    //
    params->tx_params->total_tx_height  = params->tx_params->m_dem[params->tx_params->tx_east_coord_idx]
                                                  [params->tx_params->tx_north_coord_idx];
    params->tx_params->total_tx_height += params->tx_params->antenna_height_AGL;

    //
    // receive clutter data from master
    //
    MPI_Recv (params->tx_params->m_clut[0],
              params->tx_params->nrows * params->tx_params->ncols,
              MPI_DOUBLE,
              _COVERAGE_MASTER_RANK_,
              MPI_ANY_TAG,
              comm,
              &status);
    if (status.MPI_ERROR)
    {
        fprintf (stderr, 
                 "*** ERROR: Incorrect receive of clutter data\n");
        exit (-1);
    }
    //
    // field-measurements matrix is only used in optimization mode
    //
    if (params->use_opt)
    {
        //
        // receive field-measurement data from master
        //
        MPI_Recv (params->tx_params->m_field_meas[0],
                  params->tx_params->nrows * params->tx_params->ncols,
                  MPI_DOUBLE,
                  _COVERAGE_MASTER_RANK_,
                  MPI_ANY_TAG,
                  comm,
                  &status);
        if (status.MPI_ERROR)
        {
            fprintf (stderr, 
                     "*** ERROR: Incorrect receive of field measurements data\n");
            exit (-1);
        }
    }
    free (rcv_tx_params);
}



/**
 * Deallocates all internal structures contained in the 
 * transmitter-parameters structure.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters.-
 *
 */
void 
free_tx_params (Parameters    *params,
                Tx_parameters *tx_params)
{
    if (params->use_opt)
    {
        if (tx_params->m_field_meas != NULL)
        {
            free (tx_params->m_field_meas[0]);
            free (tx_params->m_field_meas);
            tx_params->m_field_meas = NULL;
        }
    }
    if (tx_params->m_antenna_loss != NULL)
    {
        free (tx_params->m_antenna_loss[0]);
        free (tx_params->m_antenna_loss);
        tx_params->m_antenna_loss = NULL;
    }
    if (tx_params->m_radio_zone != NULL)
    {
        free (tx_params->m_radio_zone[0]);
        free (tx_params->m_radio_zone);
        tx_params->m_radio_zone = NULL;
    }
    if (tx_params->m_obst_height != NULL)
    {
        free (tx_params->m_obst_height[0]);
        free (tx_params->m_obst_height);
        tx_params->m_obst_height = NULL;
    }
    if (tx_params->m_obst_dist != NULL)
    {
        free (tx_params->m_obst_dist[0]);
        free (tx_params->m_obst_dist);
        tx_params->m_obst_dist = NULL;
    }
    if (tx_params->m_obst_offset != NULL)
    {
        free (tx_params->m_obst_offset[0]);
        free (tx_params->m_obst_offset);
        tx_params->m_obst_offset = NULL;
    }
    if (tx_params->m_loss != NULL)
    {
        free (tx_params->m_loss[0]);
        free (tx_params->m_loss);
        tx_params->m_loss = NULL;
    }
    if (tx_params->m_dem != NULL)
    {
        free (tx_params->m_dem[0]);
        free (tx_params->m_dem);
        tx_params->m_dem = NULL;
    }
    if (tx_params->m_clut != NULL)
    {
        free (tx_params->m_clut[0]);
        free (tx_params->m_clut);
        tx_params->m_clut = NULL;
    }
    if (tx_params->diagram != NULL)
    {
        free (tx_params->diagram->horizontal);
        free (tx_params->diagram->vertical);
        free (tx_params->diagram);
        tx_params->diagram = NULL;
    }
    if (params->use_gpu)
        close_gpu (tx_params);
}



/**
 * Starts a worker process.
 *
 * rank The rank of this worker process;
 * comm the MPI communicator to use.-
 *
 */
void 
worker (const int rank,
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
                optimize (params,
                          params->tx_params);
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
    free_tx_params (params,
                    params->tx_params);
    free (params->tx_params);
    free (params);
}

