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
        MPI_Send (&(tx_params->m_field_meas[tx_params->map_north_idx][tx_params->map_west_idx]),
                  1,
                  Radius_area,
                  worker_rank,
                  1,
                  comm);
    }
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
    // by default, a transmitter subregion is the whole input region;
    // this is recalculated on the workers when using the MPI implementation
    // to lower the memory consumption
    //
    tx_params->nrows              = params->nrows;
    tx_params->ncols              = params->ncols;
    tx_params->map_north          = params->map_north;
    tx_params->map_east           = params->map_east;
    tx_params->map_south          = params->map_south;
    tx_params->map_west           = params->map_west;
    tx_params->map_north_idx      = params->nrows - 1;
    tx_params->map_west_idx       = 0;
    tx_params->map_east_idx       = params->ncols - 1;
    tx_params->map_south_idx      = 0;
    tx_params->tx_north_coord_idx = (tx_params->map_north - tx_params->tx_north_coord) /
                                    params->map_ns_res;
    tx_params->tx_east_coord_idx  = (tx_params->map_east - tx_params->tx_east_coord) /
                                    params->map_ew_res;
    //
    // initialize the pointers within the transmitter structure
    //
    tx_params->m_dem            = params->m_dem;
    tx_params->m_clut           = params->m_clut;
    tx_params->m_field_meas     = NULL;
    tx_params->field_meas_count = 0;
    tx_params->m_loss           = params->m_loss;
    tx_params->m_antenna_loss   = NULL;
    tx_params->m_radio_zone     = NULL;

    //
    // initialize the rest of the data structures within this structure
    // 
    init_tx_params (params,
                    tx_params,
                    0);
    //
    // if in optimization mode, load the field measurements of this Tx
    //
    if (params->use_opt)
    {
        //
        // load the measurements from a raster map
        //
        load_field_measurements_from_map (tx_params->field_meas_map,
                                          tx_params->m_field_meas);
        tx_params->m_field_meas = params->m_field_meas;
    }
    else
        fprintf (stdout, 
                 "*** INFO: Not loading transmitter measurements in coverage prediction mode\n");
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
    // initialize the clutter categories to zero (0)
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
        // check the size of a map cell in bytes, so that later casts won't fail
        //
        if (metadata->format == sizeof (params->null_value) - 1)
        {
            G_set_f_null_value ((FCELL *) &(params->null_value), 1);
        }
	    else
        {
            fprintf (stdout, 
                     "*** WARNING: DEM map-cell format (%d) is not recognized.\n",
                     metadata->format);
            fprintf (stdout, 
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
        // check the size of a map cell in bytes, so that later casts won't fail
        //
        if (metadata->format == sizeof (params->null_value) - 1)
        {
            G_set_f_null_value ((FCELL *) &(params->null_value), 1);
        }
        else
        {
            fprintf (stdout, 
                     "*** WARNING: Clutter map-cell format (%d) is not recognized.\n",
                     metadata->format);
            fprintf (stdout, 
                     "\tMake sure it is of type FCELL to avoid undefined behaviour!\n");
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
        params->m_loss = prato_alloc_double_matrix (params->nrows,
                                                    params->ncols,
                                                    params->m_loss);
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
    char *tx_sections [10 * _CHAR_BUFFER_SIZE_];
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

    MPI_Comm worker_comm;
    MPI_Init  (&argc, &argv);

    worker_comm = MPI_COMM_WORLD;
    MPI_Comm_size (worker_comm, &nworkers);

    //
    // one is the master process
    //
    nworkers -= 1;

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

#ifdef _PERFORMANCE_METRICS_
        //
        // previous coverage calculation and result dump finished
        //
        measure_time_id (NULL, worker_rank);
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
                              worker_comm);
#ifdef _PERFORMANCE_METRICS_
                    measure_time_id ("Transmitter data send", 
                                     worker_rank);
#endif
                    //
                    // starting transmitter-data send
                    //
                    send_tx_data (params,
                                  &(params->tx_params[-- tx_count]),
                                  worker_comm,
                                  worker_rank);
#ifdef _PERFORMANCE_METRICS_
                    //
                    // transmitter-data send finished,
                    //
                    measure_time_id (NULL, 
                                     worker_rank);
                    //
                    // start coverage calculation
                    // 
                    measure_time_id ("Calculation and result dump",
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

