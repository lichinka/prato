#include "master/master.h"



/**
 * Distribute common input data among all workers.
 *
 * params       A structure holding all parameters needed for calculation;
 * comm         the communicator used to communicate with the workers.-
 *
 */
static void 
distribute_common_data (Parameters *params,
                        MPI_Comm comm)
{
    //
    // broadcast the common parameters structure
    //
    MPI_Bcast (params,
               sizeof (Parameters),
               MPI_BYTE,
               _COVERAGE_MASTER_RANK_,
               comm);
}



/**
 * Sends transmitter input data to the specified worker.
 *
 *
 * params       A structure holding all parameters needed for calculation;
 * tx_params    a structure holding transmitter-specific parameters;
 * comm         the communicator used to communicate with the worker;
 * worker_rank  the target worker.-
 *
 */
static void 
send_tx_data (Parameters *params,
              Tx_parameters *tx_params,
              MPI_Comm comm,
              int worker_rank)
{
    //
    // MPI data type for the radius-calculation area of this transmitter
    //
    MPI_Datatype Radius_area;
    MPI_Type_vector (tx_params->nrows,
                     tx_params->ncols,
                     params->ncols,
                     MPI_DOUBLE,
                     &Radius_area);
    MPI_Type_commit (&Radius_area);

    //
    // send the transmitter-specific parameters to this worker
    //
    MPI_Send (tx_params,
              sizeof (Tx_parameters),
              MPI_BYTE,
              worker_rank,
              1,
              comm);
    //
    // send DEM data to worker
    //
    MPI_Send (&(params->m_dem[tx_params->map_north_idx][tx_params->map_west_idx]),
              1,
              Radius_area,
              worker_rank,
              1,
              comm);
    //
    // send clutter data to worker
    //
    MPI_Send (&(params->m_clut[tx_params->map_north_idx][tx_params->map_west_idx]),
              1,
              Radius_area,
              worker_rank,
              1,
              comm);
    //
    // also send the field measurements matrix if in optimization mode
    //
    if (params->use_opt)
    {
        MPI_Send (&(params->m_field_meas[tx_params->map_north_idx][tx_params->map_west_idx]),
                  1,
                  Radius_area,
                  worker_rank,
                  1,
                  comm);
    }
}



/**
 * Receives the path-loss results of one transmitter from a worker process.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * worker_rank  rank of the worker sending the data;
 * tx_name      name of the transmitter for which the results are being 
 *              received;
 * comm         the MPI communicator to use.-
 *
 */
static void 
receive_tx_results (Parameters *params,
                    const int  worker_rank,
                    const char *tx_name,
                    MPI_Comm   *comm)
{
    MPI_Status status;
    Tx_parameters *tx_params = params->tx_params;

    //
    // look for the correct transmitter
    //
    int tx_id;
    for (tx_id = 0; tx_id < params->ntx; tx_id ++)
    {
        if (strncmp (tx_name, tx_params[tx_id].tx_name, strlen (tx_name)) == 0)
            break;
    }
    //
    // the transmitter is always found
    // 
    assert (tx_id < params->ntx);

    //
    // the received results will be saved here
    //
    double *rcv_results = (double *) calloc (sizeof (double), 
                                             tx_params[tx_id].nrows *
                                             tx_params[tx_id].ncols);
    //
    // receive the transmitter-specific results
    //
    MPI_Recv (rcv_results,
              tx_params[tx_id].nrows * tx_params[tx_id].ncols,
              MPI_DOUBLE,
              worker_rank,
              _WORKER_SENDING_RESULT_TAG_,
              *comm,
              &status);
    if (status.MPI_ERROR)
    {
        fprintf (stderr, 
                 "*** ERROR: Transmitter parameters incorrectly received\n");
        fflush (stderr);
        exit (1);
    }
    else
    {
        //
        // save the received results on the corresponding
        // location of the transmitter
        //
        int r_work, c_work;

        /*
        // DEBUG only!
        //
        fprintf (stdout,
                 "-- Rank %d (%s) -- Map extent (N, W)\t(%.2f,%.2f)\n",
                 worker_rank,
                 tx_name,
                 tx_params[tx_id].map_north,
                 tx_params[tx_id].map_west);*/

        for (r_work = 0; r_work < tx_params[tx_id].nrows; r_work ++)
        {
            int r_mast  = r_work + tx_params[tx_id].map_north_idx;

            for (c_work = 0; c_work < tx_params[tx_id].ncols; c_work ++)
            {
                int c_mast  = c_work + tx_params[tx_id].map_west_idx;
                int elem_id = r_work * tx_params[tx_id].ncols + c_work;

                //
                // best server aggregation
                //
                if (!isnan (rcv_results[elem_id]))
                {
                    double signal_level = tx_params[tx_id].tx_power;
                    signal_level -= rcv_results[elem_id];

                    if (isnan (params->m_loss[r_mast][c_mast]))
                        params->m_loss[r_mast][c_mast] = signal_level;
                    else
                    {
                        if (signal_level > params->m_loss[r_mast][c_mast])
                            params->m_loss[r_mast][c_mast] = signal_level;
                    }
                }
                /*
                // DEBUG only!
                //
                float east_coord  = tx_params[tx_id].map_west + c_work * params->map_ew_res;
                float north_coord = tx_params[tx_id].map_north - r_work * params->map_ns_res;

                float rcv_signal = (float) rcv_results[elem_id];

                if ((!isnan (rcv_signal)) && (rcv_signal != params->fcell_null_value))
                    fprintf (stdout, "%.2f|%.2f|%.5f\n", east_coord,
                                                         north_coord,
                                                         rcv_signal);
                */
            }
        }
    }
    free (rcv_results);
}



/**
 * Initializes the coverage calculation by reading the configuration 
 * parameters in the INI file passed as argument.
 * This function uses the 'params' parameter to save its results.-
 *
 * ini_file     Pointer to the stream object containing the INI file;
 * tx_section   name of the section containing transmitter-specific parameters;
 * params       a structure holding configuration parameters which are common
 *              to all transmitters;
 * tx_params    the output parameter: a structure holding transmitter-specific
 *              configuration parameters needed for calculation.-
 *
 */
static void 
init_coverage_for_tx (FILE          *ini_file,
                      const char    *tx_section,
                      Parameters    *params,
                      Tx_parameters *tx_params)
{
    int errno;

    //
    // parse the INI file containing the transmitter's configuration values
    //
    rewind (ini_file);
    errno = ini_parse_file (ini_file, 
                            tx_params_handler, 
                            tx_params,
                            tx_section);
    if (errno < 0)
        G_fatal_error ("Can't parse INI memory buffer\n");

    //
    // calculate the subregion (within the area) where this transmitter is 
    // located, taking into account its location and the calculation radius
    //
    double radius_in_meters       = params->radius * 1000;
    int radius_in_pixels          = (int) (radius_in_meters / params->map_ew_res);
    tx_params->nrows              = 2 * radius_in_pixels;
    tx_params->ncols              = 2 * radius_in_pixels;
    tx_params->map_north          = tx_params->tx_north_coord + radius_in_meters;
    tx_params->map_east           = tx_params->tx_east_coord + radius_in_meters;
    tx_params->map_south          = tx_params->tx_north_coord - radius_in_meters;
    tx_params->map_west           = tx_params->tx_east_coord - radius_in_meters;
    tx_params->map_north_idx      = (int) ((params->map_north - tx_params->map_north) /
                                            params->map_ns_res);
    tx_params->map_east_idx       = tx_params->map_west_idx + tx_params->ncols;
    tx_params->map_south_idx      = tx_params->map_north_idx + tx_params->nrows;
    tx_params->map_west_idx       = (int) ((tx_params->map_west - params->map_west) / 
                                            params->map_ew_res);
    tx_params->tx_north_coord_idx = radius_in_pixels;
    tx_params->tx_east_coord_idx  = radius_in_pixels;

    //
    // initialize the pointers within the transmitter structure
    //
    tx_params->diagram            = NULL;
    tx_params->m_dem              = params->m_dem;
    tx_params->m_dem_dev          = NULL;
    tx_params->m_clut             = params->m_clut;
    tx_params->m_clut_dev         = NULL;
    tx_params->m_field_meas       = NULL;
    tx_params->m_field_meas_dev   = NULL;
    tx_params->field_meas_count   = 0;
    tx_params->m_loss             = params->m_loss;
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



/**
 * Splits the comma-separated list of INI sections received.
 * This function writes its output in the last parameter, returning the
 * number of sections found.
 *
 * tx_section_list  The comma-separated string to be splitted;
 * tx_sections      an array of strings containing each of the sections.-
 *                  
 */
static int 
split_sections (char *tx_section_list,
                char *tx_sections [])
{
    int ret_value = 0;
    char *list_ptr;

    list_ptr = &(tx_section_list[0]);

    while (list_ptr != NULL)
    {
        char *found;
        // look for a comma
        found = strchr (list_ptr, ',');
        // ignore possible white spaces
        while (isspace (*list_ptr))
            list_ptr ++;
        // section name found
        tx_sections[ret_value ++] = list_ptr;
        // move to the first character after the comma
        list_ptr = found + 1;
        // smash the comma found with the end-string character
        if (found != NULL)
            *found = '\0';
        else
            break;
    }
    return ret_value;
}



/**
 * Calculates the coverage prediction of several transmitters over MPI.
 *
 * params       a structure holding all parameters needed for calculation;
 * nworkers     the number of available workers within the group;
 * worker_comm  the MPI communicator, to which the workers are connected;
 *
 */
static void 
coverage_mpi (Parameters *params,
              const int   nworkers,
              MPI_Comm   *worker_comm)
{
    MPI_Status status;
    char buff [_CHAR_BUFFER_SIZE_];

    //
    // start optimization-processing loop
    //
    int worker_rank;
    int tx_count = params->ntx;
    int running_workers = nworkers;

    while (running_workers > 0)
    {   
        MPI_Recv (&buff, 
                  _CHAR_BUFFER_SIZE_,
                  MPI_BYTE,
                  MPI_ANY_SOURCE,
                  MPI_ANY_TAG,
                  *worker_comm,
                  &status);
        worker_rank = status.MPI_SOURCE;

#ifdef _PERFORMANCE_METRICS_
        //
        // previous result dump finished
        //
        measure_time_id (NULL, 
                         worker_rank);
        //
        // async messaging for heterogeneous-system support
        //
        measure_time_id ("Async message",
                         worker_rank);
#endif

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
                              *worker_comm);
#ifdef _PERFORMANCE_METRICS_
                    //
                    // async messaging finished
                    //
                    measure_time_id (NULL,
                                     worker_rank);
                    //
                    // start transmitter-data send
                    //
                    measure_time_id ("Transmitter data send", 
                                     worker_rank);
#endif
                    //
                    // starting transmitter-data send
                    //
                    send_tx_data (params,
                                  &(params->tx_params[-- tx_count]),
                                  *worker_comm,
                                  worker_rank);
#ifdef _PERFORMANCE_METRICS_
                    //
                    // transmitter-data send finished
                    //
                    measure_time_id (NULL, 
                                     worker_rank);
                    //
                    // start coverage calculation
                    // 
                    measure_time_id ("Coverage calculation",
                                     worker_rank);
#endif
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
                              *worker_comm);
                    running_workers --;
#ifdef _PERFORMANCE_METRICS_
                    //
                    // async messaging finished
                    //
                    measure_time_id (NULL,
                                     worker_rank);
                    //
                    // coverage calculation and result dump finished
                    //
                    measure_time_id (NULL, 
                                     worker_rank);
#endif
                }
                break;

            case (_WORKER_SENDING_RESULT_TAG_):
#ifdef _PERFORMANCE_METRICS_
                //
                // coverage calculation finished,
                // 
                //
                measure_time_id (NULL, 
                                 worker_rank);
                //
                // start result dump
                // 
                measure_time_id ("Result dump",
                                 worker_rank);
#endif
                receive_tx_results (params,
                                    worker_rank,
                                    buff,
                                    worker_comm);
                break;

            default:
                fprintf (stderr, 
                         "*** WARNING: Unknown message from %d. worker\n", 
                         worker_rank);
        }   
    }
    //
    // output the coverage data that has been already aggregated
    //
    int r, c;

    fprintf (stderr,
             "*** WARNING: the following data is valid only in the case when workers send the\n");
    fprintf (stderr,
             "\tintermediate results back to the master; otherwise external DB\n");
    fprintf (stderr,
             "\taggregation is needed for the coverage to be valid\n");

    for (r = 0; r < params->nrows; r ++)
    {
        for (c = 0; c < params->ncols; c ++)
        {
            float east_coord  = params->map_west + c * params->map_ew_res;
            float north_coord = params->map_north - r * params->map_ns_res;

            float rcv_signal = (float) params->m_loss[r][c];

            if ((!isnan (rcv_signal)) && (rcv_signal != params->fcell_null_value) && (rcv_signal != 0))
                fprintf (stdout, "%.2f|%.2f|%.5f\n", east_coord,
                                                     north_coord,
                                                     rcv_signal);
        }
    }
}



/**
 * Optimizes the clutter-category losses using differential evolution on
 * the master process and evaluating the objective function on the workers.
 *
 * params       a structure holding all parameters needed for calculation;
 * nworkers     the number of available workers within the group;
 * worker_comm  the MPI communicator, to which the workers are connected;
 *
 */
static void 
optimize_mpi (Parameters *params,
              const int   nworkers,
              MPI_Comm   *worker_comm)
{
    MPI_Status status;
    char buff [_CHAR_BUFFER_SIZE_];

    //
    // start coverage-processing loop
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
                  *worker_comm,
                  &status);
        worker_rank = status.MPI_SOURCE;
        switch (status.MPI_TAG)
        {   
            case (_WORKER_IS_IDLE_TAG_):
                //
                // send transmitter data, before starting the optimization loop
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
                              *worker_comm);
#ifdef _PERFORMANCE_METRICS_
                    measure_time_id ("Load field measurements", 
                                     worker_rank);
#endif
                    //
                    // allocate a matrix for loading the field measurements
                    //
                    if (params->m_field_meas == NULL)
                        params->m_field_meas = prato_alloc_double_matrix (params->nrows,
                                                                          params->ncols,
                                                                          params->m_field_meas);
                    //
                    // load the measurements from a raster map
                    //
                    load_field_measurements_from_map (params->tx_params[-- tx_count].field_meas_map,
                                                      params->m_field_meas);
#ifdef _PERFORMANCE_METRICS_
                    measure_time_id (NULL,
                                     worker_rank);
                    measure_time_id ("Transmitter data send", 
                                     worker_rank);
#endif
                    //
                    // starting transmitter-data send
                    //
                    send_tx_data (params,
                                  &(params->tx_params[tx_count]),
                                  *worker_comm,
                                  worker_rank);
#ifdef _PERFORMANCE_METRICS_
                    //
                    // transmitter-data send finished,
                    //
                    measure_time_id (NULL, 
                                     worker_rank);
#endif
                    //
                    // clear the field-measurement matrix to
                    // avoid data clashing with the next transmitter
                    //
                    if (params->m_field_meas != NULL)
                    {
                        free (params->m_field_meas[0]);
                        free (params->m_field_meas);
                        params->m_field_meas = NULL;
                    }
                    //
                    // reduce the number of transmitter-data still to be sent
                    //
                    running_workers --;
                }
                break;

            default:
                fprintf (stderr, 
                         "WARNING Unknown message from %d. worker\n", 
                         worker_rank);
        }   
    }
    fprintf (stdout,
             "*** INFO: All transmitter data sent. Starting optimization ...\n");
    optimize_on_master (params,
                        params->tx_params,
                        worker_comm);
}



/**
 * Initializes the coverage calculation by reading the configuration 
 * parameters in the [Common] section of the INI file passed as argument.
 * This function returns a pointer to the newly created parameters structure.
 *
 * ini_file         Pointer to the stream containing the configuration read
 *                  from the INI file;
 * tx_sections_list the names of the sections containing transmitter-specific
 *                  configuration;
 * params           the output parameter, into which everything is saved.-
 *
 */
void 
init_coverage (FILE       *ini_file,
               char       *tx_sections_list,
               Parameters *params)
{
    int i, errno;

    //
    // initialize matrix pointers inside the Parameters structure
    //
    params->m_dem        = NULL;
    params->m_clut       = NULL;
    params->m_loss       = NULL;
    params->m_field_meas = NULL;

    //
    // initialize all clutter categories to zero (0)
    //
    params->clutter_category_count = 0;
    for (i = 0; i < _CHAR_BUFFER_SIZE_; i ++)
        params->clutter_loss[i] = 0;

    //
    // parse the INI file containing the common configuration values
    //
    rewind (ini_file);
    errno = ini_parse_file (ini_file,
                            common_params_handler, 
                            params,
                            NULL);
    if (errno < 0)
        G_fatal_error ("Can't parse INI memory buffer\n");
    
    //
    // mapset names are kept here
    //
    char *mapset  = (char *) calloc (_CHAR_BUFFER_SIZE_, sizeof (char));
    char *mapset2 = (char *) calloc (_CHAR_BUFFER_SIZE_, sizeof (char));
    int  row, col;
    int  infd, infd2;                           // file descriptors

    //
    // returns NULL if the map was not found in any mapset
    //
    mapset = G_find_cell2 (params->dem_map, mapset);
    if (mapset == NULL)
        G_fatal_error("Raster map <%s> not found", params->dem_map);
   
    mapset2 = G_find_cell2 (params->clutter_map, mapset2);
    if (mapset2 == NULL)
        G_fatal_error("Raster map <%s> not found", params->clutter_map);
    //
    // returns a negative value if the map cannot be opened
    //
    infd = G_open_cell_old (params->dem_map, mapset);
    if (infd < 0)
        G_fatal_error("Unable to open raster map <%s>", params->dem_map);

    infd2 = G_open_cell_old (params->clutter_map, mapset2);
    if (infd2 < 0)
        G_fatal_error("Unable to open raster map <%s>", params->clutter_map);

    //
    // read metadata of each map, making sure they match
    //
    struct Cell_head *metadata = (struct Cell_head *) malloc (sizeof (struct Cell_head));
    struct Cell_head *window   = (struct Cell_head *) malloc (sizeof (struct Cell_head));

    //
    // DEM metadata
    //
    errno = G_get_cellhd (params->dem_map,
                          mapset,
                          metadata);
    if (errno == 0)
    {
        //
        // NULL value in DEM maps
        //
        G_set_f_null_value ((FCELL *) &(params->fcell_null_value), 1);
        if (metadata->format != sizeof (params->fcell_null_value) - 1)
        {
            fprintf (stderr, 
                     "*** WARNING: DEM map-cell format (%d) is not recognized.\n",
                     metadata->format);
            fprintf (stderr, 
                     "\tMake sure it is of type FCELL to avoid undefined behaviour!\n");
        }
        //
        // get map metadata
        //
        params->map_east = metadata->east;
        params->map_west = metadata->west;
        params->map_north = metadata->north;
        params->map_south = metadata->south;
        params->map_ew_res = metadata->ew_res;
        params->map_ns_res = metadata->ns_res;
    }
    else
        G_fatal_error ("Unable to open raster map <%s>", 
                       params->dem_map);
    //
    // clutter metadata
    //
    errno = G_get_cellhd (params->clutter_map,
                          mapset2,
                          metadata);
    if (errno == 0)
    {
        //
        // NULL value in clutter maps 
        //
        G_set_c_null_value ((CELL *) &(params->cell_null_value), 1);
        if (metadata->format != sizeof (params->cell_null_value) - 1)
        {
            fprintf (stderr, 
                     "*** WARNING: Clutter map-cell format (%d) is not recognized.\n",
                     metadata->format);
            fprintf (stderr, 
                     "\tMake sure it is of type CELL to avoid undefined behaviour!\n");
        }
        if (!(params->map_east >= metadata->east &&
              params->map_west <= metadata->west &&
              params->map_north >= metadata->north &&
              params->map_south <= metadata->south &&
              params->map_ew_res == metadata->ew_res &&
              params->map_ns_res == metadata->ns_res))
            G_fatal_error ("Map metadata of input maps do not match.");
    }
    else
        G_fatal_error ("Unable to open raster map <%s>", params->clutter_map);

    //
    // check that the current active window matches the loaded data
    //
	G_get_window (window);
    if (params->map_east != window->east ||
        params->map_west != window->west ||
        params->map_north != window->north ||
        params->map_south != window->south ||
        params->map_ew_res != window->ew_res ||
        params->map_ns_res != window->ns_res)
        G_fatal_error ("Loaded map metadata does not match with your current GRASS window.\nRun 'g.region -p' to check the settings.");

    //
    // number of rows and columns within the maps
    //
    params->nrows = metadata->rows;
    params->ncols = metadata->cols;

    //
    // allocate memory for the resulting path-loss matrix
    //
    if (params->m_loss == NULL)
    {
        params->m_loss = prato_alloc_double_matrix (params->nrows,
                                                    params->ncols,
                                                    params->m_loss);
        //
        // initialize all matrix elements to the NULL value
        //
        int r, c;
        for (r = 0; r < params->nrows; r ++)
            for (c = 0; c < params->ncols; c ++)
                params->m_loss[r][c] = params->fcell_null_value;
    }
    //
    // allocate the reading buffers for DEM and clutter data
    //
    void *inrast  = G_allocate_raster_buf (FCELL_TYPE);
    void *inrast2 = G_allocate_raster_buf (FCELL_TYPE);

    //
    // allocate memory to contain the both DEM and clutter maps
    //
    if ((params->m_dem == NULL) && (params->m_clut == NULL))
    {
        // DEM
        params->m_dem = prato_alloc_double_matrix (params->nrows,
                                                   params->ncols,
                                                   params->m_dem);
        // CLUTTER
        params->m_clut = prato_alloc_double_matrix (params->nrows,
                                                    params->ncols,
                                                    params->m_clut);
        //
        // read files (DEM and clutter) into their respective arrays
        //
        for (row = 0; row < params->nrows; row ++) 
        {	
            // read DEM map data
            if (G_get_raster_row (infd, inrast, row, FCELL_TYPE) < 0)
                G_fatal_error ("Unable to read raster map <%s> row %d", 
                               params->dem_map, 
                               row);

            // read Clutter map data
            if (G_get_raster_row (infd2, inrast2, row, FCELL_TYPE) < 0)
                G_fatal_error ("Unable to read raster map <%s> row %d",
                               params->clutter_map, 
                               row);

            // process the data
            for (col = 0; col < params->ncols; col ++) 
            { 
                FCELL f_in = ((FCELL *) inrast)[col];
                FCELL f_in2 = ((FCELL *)inrast2)[col];
                params->m_dem[row][col] = (double) f_in;
                params->m_clut[row][col] = (double) f_in2;
             }
        } 
    }
    //
    // create an array containing the transmitters to be processed
    //
    char *tx_sections [20 * _CHAR_BUFFER_SIZE_];
    params->ntx = split_sections (tx_sections_list,
                                  tx_sections);
    // 
    // allocate and parse the per-transmitter parameters structure
    //
    params->tx_params = (Tx_parameters *) calloc (params->ntx,
                                                  sizeof (Tx_parameters));
    for (i = 0; i < params->ntx; i ++)
    {
        init_coverage_for_tx (ini_file,
                              tx_sections[i],
                              params,
                              &(params->tx_params[i]));
        printf ("*** INFO: Tx [%s] initialized!\n",
                params->tx_params[i].tx_name);
    }
    //
    // free the allocated metadata
    //
    free (window);
    free (metadata);
    free (mapset);
    free (mapset2);
    //
    // free the read buffers for DEM and clutter data
    //
    G_free (inrast2);
    G_free (inrast);
    //
    // close raster maps
    // 
    G_close_cell (infd);
    G_close_cell (infd2);

#ifdef _DEBUG_INFO_
    int r, c;
    for (r = 0; r < params->tx_params->nrows; r ++)
    {
        for (c = 0; c < params->tx_params->ncols; c ++)
        {
            fprintf (stdout,
                     "%.f|%.f|%.f\n",
                     params->map_west + 12.5 + (c * params->map_ew_res),
                     params->map_north - 12.5 - (r * params->map_ns_res),
                     params->tx_params->m_clut[r][c]);
        }
    }
    fflush (stderr);
    exit (1);
#endif
}



/**
 * Initializes the MPI environment for master and workers.
 *
 * argc             Number of command line parameters;
 * argv             array containing command line parameters;
 * params           a structure holding all parameters needed for calculation.-
 *
 */
void 
init_mpi (int argc, 
          char *argv [],
          Parameters *params)
{
    int nworkers = -1;

    MPI_Comm worker_comm;
    MPI_Init  (&argc, &argv);

    worker_comm = MPI_COMM_WORLD;
    MPI_Comm_size (worker_comm, &nworkers);

    //
    // one is the master process
    //
    nworkers -= 1;

    //
    // in optimization mode, we may only process as many transmitters
    // as there are worker processes
    //
    if ((params->use_opt) && (params->ntx > nworkers))
    {
        fprintf (stderr, 
                 "*** WARNING Only the first %d transmitters will take part in the optimization\n",
                 nworkers);
        params->ntx = nworkers;
    }
    //
    // sync point: pass common input data to all workers
    //
    MPI_Barrier (worker_comm);
#ifdef _PERFORMANCE_METRICS_
    measure_time ("Common data distribution");
#endif
    distribute_common_data (params,
                            worker_comm);
    //
    // sync point: common data distribution finished
    //
    MPI_Barrier (worker_comm);
#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
#endif

    //
    // do we have to start coverage or optimization calculation?
    //
    if (params->use_master_opt)
        optimize_mpi (params,
                      nworkers,
                      &worker_comm);
    else 
        coverage_mpi (params,
                      nworkers,
                      &worker_comm);
    //
    // close the MPI environment
    //
    MPI_Finalize ( );
}

