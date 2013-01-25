#include <performance_metric.h>
#include "measurement.h"
#include "worker/coverage.h"



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

	printf ("*** DEBUG: World size %d\n", world_size);
	printf ("*** DEBUG: Rank %d\n", world_size);

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

	printf ("*** DEBUG: Universe size %d\n", universe_size);

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
        fprintf (stderr, 
                 "*** WARNING There are not enough slots to process all transmitters in parallel\n");
    }
    else if (nworkers > params->ntx)
    {
        nworkers = params->ntx;
    }
    //
    // prepare workers' commands and arguments
    //
    int buff_size = _CHAR_BUFFER_SIZE_;
    //char    *worker_command = "run_worker.sh";
    char    *worker_command = "hostname";
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
		//
		// test dynamic spawning with the hostname command
		//
        //snprintf (worker_arg[i], buff_size, "%d", i);
        snprintf (worker_arg[i], buff_size, "-s");
        worker_args[i][0] = worker_arg[i];
        worker_args[i][1] = NULL;
        worker_argv[i]    = &(worker_args[i][0]);

        // worker's command
        worker_program[i] = (char *) malloc (buff_size);
        snprintf (worker_program[i], buff_size, worker_command);

        worker_maxproc[i] = 1;

		// 
		// manually set host info to test the problem on ninestein
		//
        worker_info[i] = MPI_INFO_NULL;
		//MPI_Info_create (&(worker_info[i]));
		//MPI_Info_set (worker_info[i], "host", "k1");

		printf ("*** DEBUG: worker_args[%d][0]\t%s\n", i, worker_args[i][0]);
		printf ("*** DEBUG: worker_argv[%d]\t%s\n", i, worker_args[i][0]);
		printf ("*** DEBUG: worker_program[%d]\t%s\n", i, worker_program[i]);
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
            fprintf (stderr, 
                     "*** WARNING Worker %d. returned with exit code %d\n", 
                     i,
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
    // returns NULL if the map was not found in any mapset
    //
    int infd, tr_row, tr_col;
    char *mapset = G_find_cell2 (params->dem_map, "");

    if (mapset == NULL)
        G_fatal_error ("DEM raster map <%s> not found", 
                       params->dem_map);
    //
    // returns file descriptor to the raster map
    //
    if ((infd = G_open_cell_old (params->dem_map, mapset)) < 0)
        G_fatal_error ("Unable to open DEM raster map <%s>",
                       params->dem_map);
    //
    // check if the specified transmitter location is inside the DEM map
    // 
    if (tx_params->tx_east_coord  < params->map_west || 
        tx_params->tx_east_coord  > params->map_east ||
        tx_params->tx_north_coord > params->map_north || 
        tx_params->tx_north_coord < params->map_south)
        G_fatal_error ("Specified TX coordinates are outside the DEM map.");
   
    //
    // map array coordinates for transmitter
    //
    tr_row = (params->map_north - tx_params->tx_north_coord) / params->map_ns_res;
    tr_col = (tx_params->tx_east_coord - params->map_west) / params->map_ew_res;

    //
    // allocate the reading buffer for DEM 
    //
    void *inrast = G_allocate_raster_buf (FCELL_TYPE);

    //
    // total height of transmitter 
    // 
    if (G_get_raster_row (infd, inrast, tr_row, FCELL_TYPE) < 0)
        G_fatal_error ("Unable to read raster map <%s> row %d", params->dem_map, 
                                                                tr_row);
    FCELL trans_elev = ((FCELL *) inrast)[tr_col];

    //
    // check if transmitter is on DEM
    //
    if (isnan ((double) trans_elev))							
        G_fatal_error ("Transmitter outside DEM raster map.");

    tx_params->total_tx_height = (double) trans_elev + tx_params->antenna_height_AGL;
  
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
    tx_params->map_north_idx      = (int) ceil ((params->map_north - tx_params->map_north) /
                                                 params->map_ns_res);
    tx_params->map_east_idx       = tx_params->map_west_idx + tx_params->ncols;
    tx_params->map_south_idx      = tx_params->map_north_idx + tx_params->nrows;
    tx_params->map_west_idx       = (int) floor ((tx_params->map_west - params->map_west) / 
                                                  params->map_ew_res);
    tx_params->tx_east_coord_idx  = radius_in_pixels - 1;
    tx_params->tx_north_coord_idx = radius_in_pixels - 1;

    //
    // if in optimization mode, load the field measurements of this Tx
    //
    if (params->use_opt)
    {
        //
        // allocate memory for the field-measurement matrix
        //
        if (tx_params->m_field_meas == NULL)
        {
            int i;
            double *m_field_meas_data = (double *) calloc (params->nrows * params->ncols,
                                                     sizeof (double));
            tx_params->m_field_meas = (double **) calloc (params->nrows, 
                                                 sizeof (double *));
            for (i = 0; i < params->nrows; i ++)
                tx_params->m_field_meas[i] = &(m_field_meas_data[i * params->ncols]);
        }
        else
            fprintf (stderr, 
                     "*** WARNING: Not initializing field measurements matrix\n");
        //
        // load the measurements from a raster map
        //
        load_field_measurements_from_map (tx_params->field_meas_map,
                                          tx_params->m_field_meas);
    }
    else
        fprintf (stdout, 
                 "*** INFO: Not loading transmitter measurements in coverage prediction mode\n");

    //
    // free the read buffer for DEM data
    //
    G_free (inrast);
    
    //
    // close raster maps
    //
    G_close_cell (infd);
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

    //
    // GPU-specific parameters are kept here
    //
    params->gpu_params = (GPU_parameters *) malloc (sizeof (GPU_parameters));
    params->gpu_params->ocl_obj = NULL;

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
    // set all evaluation parameters up
    //
    char *mapset, *mapset2;		    // mapset names
    int row, col;
    int infd, infd2;                // file descriptors

    //
    // returns NULL if the map was not found in any mapset
    //
    mapset = G_find_cell2 (params->dem_map, "");
    if (mapset == NULL)
        G_fatal_error("Raster map <%s> not found", params->dem_map);
   
    mapset2 = G_find_cell2 (params->clutter_map, "");
    if (mapset2 == NULL)
        G_fatal_error("Raster map <%s> not found", params->clutter_map);

    // G_open_cell_old - returns file destriptor (>0) 
    if ((infd = G_open_cell_old (params->dem_map, mapset)) < 0)
        G_fatal_error("Unable to open raster map <%s>", params->dem_map);

    if ((infd2 = G_open_cell_old (params->clutter_map, mapset2)) < 0)
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
        if (metadata->format == 0)
        {
            fprintf (stdout, 
                     "*** WARNING: DEM map-cell format (%d) is not recognized.\n",
                     metadata->format);
            fprintf (stdout, 
                     "\tMake sure it is of type FCELL to avoid undefined behaviour!\n");
        }
        else
        {
            if (metadata->format != sizeof (params->null_value) - 1)
                G_fatal_error ("DEM: the number of bytes per map-cell (%d) is not castable to `float`", 
                               metadata->format);
            else
                G_set_f_null_value ((FCELL *) &(params->null_value), 1);
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
        G_fatal_error("Unable to open raster map <%s>", params->dem_map);
       
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
        if (metadata->format == 0)
        {
            fprintf (stdout, 
                     "*** WARNING: Clutter map-cell format (%d) is not recognized.\n",
                     metadata->format);
            fprintf (stdout, 
                     "\tMake sure it is of type FCELL to avoid undefined behaviour!\n");
        }
        else
        {
            if (metadata->format != sizeof (params->null_value) - 1)
                G_fatal_error ("Clutter: the number of bytes per map-cell (%d) is not castable to `float`", 
                               metadata->format);
            else
                G_set_f_null_value ((FCELL *) &(params->null_value), 1);
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
        double *m_loss_data = (double *) calloc (params->nrows * params->ncols,
                                                 sizeof (double));
        params->m_loss = (double **) calloc (params->nrows, 
                                             sizeof (double *));
        for (i = 0; i < params->nrows; i ++)
            params->m_loss[i] = &(m_loss_data[i * params->ncols]);
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
        double *m_dem_data = (double *) calloc (params->nrows * params->ncols,
                                                sizeof (double));
        params->m_dem = (double **) calloc (params->nrows,
                                            sizeof (double *));
        for (i = 0; i < params->nrows; i ++)
            params->m_dem[i] = &(m_dem_data[i * params->ncols]);

        // CLUTTER
        double *m_clut_data = (double *) calloc (params->nrows * params->ncols,
                                                 sizeof (double));
        params->m_clut = (double **) calloc (params->nrows, 
                                             sizeof (double *));
        for (i = 0; i < params->nrows; i ++)
            params->m_clut[i] = &(m_clut_data[i * params->ncols]);

        //
        // read files (DEM and clutter) into their respective arrays
        //
        for (row = 0; row < params->nrows; row ++) 
        {	
            // read DEM map data
            if (G_get_raster_row (infd, inrast, row, FCELL_TYPE) < 0)
              G_fatal_error ("Unable to read raster map <%s> row %d", params->dem_map, row);

            // read Clutter map data
            if (G_get_raster_row (infd2, inrast2, row, FCELL_TYPE) < 0)
              G_fatal_error ("Unable to read raster map <%s> row %d", params->clutter_map, row);

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
    char *tx_sections [10240];
    params->ntx = split_sections (tx_sections_list,
                                  tx_sections);
    // 
    // allocate and parse the per-transmitter parameters structure
    //
    params->tx_params = (Tx_parameters *) calloc (params->ntx,
                                                  sizeof (Tx_parameters));
    for (i = 0; i < params->ntx; i ++)
    {
        params->tx_params[i].m_field_meas = NULL;
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
    measure_time ("Common data distribution");
    distribute_common_data (params,
                            worker_comm);
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
        // previous coverage calculation and result dump finished
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
                    send_tx_data (params,
                                  &(params->tx_params[-- tx_count]),
                                  worker_comm,
                                  worker_rank);
                    //
                    // transmitter-data send finished,
                    //
                    measure_time_id (NULL, worker_rank);
                    //
                    // start coverage calculation
                    // 
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

