#include "worker/coverage.h"
#include "worker/antenna.h"
#include "worker/eric.h"




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
static void init_coverage_for_tx (FILE          *ini_file,
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
    // set all evaluation parameters for this transmitter up
    //
    char *mapset;
    int infd, tr_row, tr_col;

    //
    // returns NULL if the map was not found in any mapset
    //
    mapset = G_find_cell2 (params->dem_map, "");
    if (mapset == NULL)
        G_fatal_error(_("Raster map <%s> not found"), params->dem_map);
  
    //
    // G_open_cell_old - returns file descriptor
    //
    if ((infd = G_open_cell_old (params->dem_map, mapset)) < 0)
        G_fatal_error(_("Unable to open raster map <%s>"), params->dem_map);

    //
    // check if the specified transmitter location is inside the DEM map
    // 
    if (tx_params->tx_east_coord  < params->map_west || 
        tx_params->tx_east_coord  > params->map_east ||
        tx_params->tx_north_coord > params->map_north || 
        tx_params->tx_north_coord < params->map_south)
        G_fatal_error (_("Specified BS coordinates are outside current region bounds."));
   
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
        G_fatal_error (_("Unable to read raster map <%s> row %d"), params->dem_map, 
                                                                   tr_row);
    FCELL trans_elev = ((FCELL *) inrast)[tr_col];

    //
    // check if transmitter is on DEM
    //
    if (isnan ((double) trans_elev))							
        G_fatal_error (_("Transmitter outside raster DEM map."));

    tx_params->total_tx_height = (double) trans_elev + tx_params->antenna_height_AGL;
   
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
static int split_sections (char *tx_section_list,
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
void init_coverage (FILE       *ini_file,
                    char       *tx_sections_list,
                    Parameters *params)
{
    int i, errno;

    //
    // initialize some pointers inside the Parameters structure
    //
    params->m_dem        = NULL;
    params->m_clut       = NULL;
    params->m_loss       = NULL;
    params->m_field_meas = NULL;

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
        G_fatal_error(_("Raster map <%s> not found"), params->dem_map);
   
    mapset2 = G_find_cell2 (params->clutter_map, "");
    if (mapset2 == NULL)
        G_fatal_error(_("Raster map <%s> not found"), params->clutter_map);

    // G_open_cell_old - returns file destriptor (>0) 
    if ((infd = G_open_cell_old (params->dem_map, mapset)) < 0)
        G_fatal_error(_("Unable to open raster map <%s>"), params->dem_map);

    if ((infd2 = G_open_cell_old (params->clutter_map, mapset2)) < 0)
        G_fatal_error(_("Unable to open raster map <%s>"), params->clutter_map);

    //
    // read metadata of each map, making sure they match
    //
    struct Cell_head *metadata = (struct Cell_head *) malloc (sizeof (struct Cell_head));
    struct Cell_head *window   = (struct Cell_head *) malloc (sizeof (struct Cell_head));
    errno = G_get_cellhd (params->dem_map,
                          mapset,
                          metadata);
    //
    // check the size of a map cell in bytes, so that later casts won't fail
    //
    if (metadata->format != sizeof (params->null_value) - 1)
        G_fatal_error ("The number of bytes per map-cell is castable to Float");
    else
        G_set_f_null_value ((FCELL *) &(params->null_value), 1);

    if (errno == 0)
    {
        params->map_east = metadata->east;
        params->map_west = metadata->west;
        params->map_north = metadata->north;
        params->map_south = metadata->south;
        params->map_ew_res = metadata->ew_res;
        params->map_ns_res = metadata->ns_res;
    }
    else
        G_fatal_error(_("Unable to open raster map <%s>"), params->dem_map);
        
    errno = G_get_cellhd (params->clutter_map,
                          mapset2,
                          metadata);
    if (errno == 0)
    {
        if (params->map_east != metadata->east ||
            params->map_west != metadata->west ||
            params->map_north != metadata->north ||
            params->map_south != metadata->south ||
            params->map_ew_res != metadata->ew_res ||
            params->map_ns_res != metadata->ns_res)
            G_fatal_error (_("Map metadata of input maps do not match."));
    }
    else
        G_fatal_error (_("Unable to open raster map <%s>"), params->clutter_map);
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
        G_fatal_error (_("Loaded map metadata does not match with your current GRASS window. Run 'g.region -p' to check the settings."));

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
              G_fatal_error (_("Unable to read raster map <%s> row %d"), params->dem_map, row);

            // read Clutter map data
            if (G_get_raster_row (infd2, inrast2, row, FCELL_TYPE) < 0)
              G_fatal_error (_("Unable to read raster map <%s> row %d"), params->clutter_map, row);

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
 * Calculates the coverage prediction for one transmitter, using the 
 * Ericsson 9999 model.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 *                  configuration parameters needed for calculation;
 * eric_params      contains the four tunning parameters for the Ericsson 9999
 *                  model;
 * eric_params_len  the number of parameters within the received vector, four 
 *                  in this case (A0, A1, A2 and A3);
 *
 */
void coverage (const Parameters     *params,
               const Tx_parameters  *tx_params,
               const double         *eric_params, 
               const unsigned int   eric_params_len)
{
    //
    // calculate the subregion (within the area) where this transmitter is 
    // located, taking into account its location and the calculation radius
    //
    double radius_in_meters     = params->radius * 1000;
    int radius_in_pixels        = (int) (radius_in_meters / params->map_ew_res);
    int diameter_in_pixels      = 2 * radius_in_pixels;
    double mini_west_border     = tx_params->tx_east_coord - radius_in_meters;
    double mini_east_border     = tx_params->tx_east_coord + radius_in_meters;
    double mini_south_border    = tx_params->tx_north_coord - radius_in_meters;
    double mini_north_border    = tx_params->tx_north_coord + radius_in_meters;
    int mini_west_border_idx    = (mini_west_border - params->map_west) / 
                                  params->map_ew_res;
    int mini_east_border_idx    = (mini_east_border - params->map_west) / 
                                   params->map_ew_res;
    int mini_south_border_idx   = (params->map_north - mini_south_border) / 
                                  params->map_ns_res;
    int mini_north_border_idx   = (params->map_north - mini_north_border) / 
                                  params->map_ns_res;
    int mini_nrows              = 2 * radius_in_pixels;
    int mini_ncols              = 2 * radius_in_pixels;
    int tx_east_coord_index     = radius_in_pixels - 1;
    int tx_north_coord_index    = radius_in_pixels - 1;

    //
    // allocate memory for the mini matrices, that contain strictly the data needed
    // for the calculation inside the transmition radius
    //
    int r, c, mini_r, mini_c;
    //
    // digital elevation model
    //
    double *mini_m_dem_data = (double *) calloc (diameter_in_pixels * diameter_in_pixels, 
                                                 sizeof (double));
    double **mini_m_dem = (double **) calloc (diameter_in_pixels,
                                              sizeof (double *));
    for (r = 0; r < diameter_in_pixels; r ++)
        mini_m_dem[r] = &(mini_m_dem_data[r * diameter_in_pixels]);
    //
    // clutter data
    //
    double *mini_m_clut_data = (double *) calloc (diameter_in_pixels * diameter_in_pixels, 
                                                  sizeof (double));
    double **mini_m_clut = (double **) calloc (diameter_in_pixels,
                                               sizeof (double *));
    for (r = 0; r < diameter_in_pixels; r ++)
        mini_m_clut[r] = &(mini_m_clut_data[r * diameter_in_pixels]);
    //
    // path loss
    //
    double *mini_m_loss_data = (double *) calloc (diameter_in_pixels * diameter_in_pixels, 
                                                  sizeof (double));
    double **mini_m_loss = (double **) calloc (diameter_in_pixels,
                                               sizeof (double *));
    for (r = 0; r < diameter_in_pixels; r ++)
        mini_m_loss[r] = &(mini_m_loss_data[r * diameter_in_pixels]);
  
    //
    // copy data from the big matrices to the mini ones
    //
    mini_r = 0;
    for (r = 0; r < params->nrows; r ++)
    {
        if ((r > mini_north_border_idx) && (r <= mini_south_border_idx))
        {
            mini_c = 0;
            for (c = 0; c < params->ncols; c ++)
            {
                if ((c > mini_west_border_idx) && (c <= mini_east_border_idx))
                {
                    mini_m_dem[mini_r][mini_c]  = params->m_dem[r][c];
                    mini_m_clut[mini_r][mini_c] = params->m_clut[r][c];
                    mini_m_loss[mini_r][mini_c] = params->m_loss[r][c];
                    mini_c ++;
                }
            }
            mini_r ++;
        }
    }

    //
    // initialize the structure variables for Ericsson 9999 model
    //
    struct StructEric IniEric = {tx_east_coord_index, 
                                 tx_north_coord_index,
                                 tx_params->antenna_height_AGL, 
                                 params->rx_height_AGL,
                                 mini_ncols,
                                 mini_nrows,
                                 params->map_ew_res,
                                 params->frequency, 
                                 (double) eric_params[0],
                                 (double) eric_params[1],
                                 (double) eric_params[2], 
                                 (double) eric_params[3], 
                                 1, 
                                 params->radius};
    //
    // execute the path-loss calculation on CPU or GPU?
    //
    if (params->use_gpu)
    {
        eric_pathloss_on_gpu (tx_params->tx_east_coord,
                              tx_params->tx_north_coord,
                              tx_east_coord_index,
                              tx_north_coord_index,
                              tx_params->antenna_height_AGL,
                              tx_params->total_tx_height,
                              tx_params->beam_direction,
                              tx_params->mechanical_tilt,
                              params->frequency,
                              params->radius,  
                              params->rx_height_AGL,
                              mini_nrows,       
                              mini_ncols,      
                              mini_west_border,
                              mini_north_border,
                              params->map_ew_res,  
                              params->map_ns_res,  
                              params->null_value,
                              mini_m_dem,          
                              mini_m_clut,
                              mini_m_loss);
    }
    else
    {
#ifdef _PERFORMANCE_METRICS_
        measure_time ("E/// on CPU");
#endif
        EricPathLossSub (mini_m_dem, 
                         mini_m_clut, 
                         mini_m_loss, 
                         &IniEric);
#ifdef _PERFORMANCE_METRICS_
        measure_time (NULL);
#endif
    }

    /*
    // DEBUG memory
    //
    for (r=0;r<mini_nrows;r++)
    {
        for (c=0;c<mini_ncols;c++)
        {
            printf ("%.5f\t", mini_m_loss[r][c]);
        }
        printf ("\n");
    }
    exit (1);
    */

    //
    // calculate the antenna influence, 
    // overwriting the isotrophic path-loss
    //
#ifdef _PERFORMANCE_METRICS_
    measure_flops ("Antenna", 1);
#endif
    calculate_antenna_influence (params->use_gpu,
                                 tx_params->tx_east_coord,
                                 tx_params->tx_north_coord,
                                 tx_params->antenna_height_AGL,
                                 tx_params->total_tx_height,
                                 tx_params->beam_direction,
                                 tx_params->mechanical_tilt,
                                 params->frequency,
                                 params->radius,  
                                 params->rx_height_AGL,
                                 mini_nrows,       
                                 mini_ncols,      
                                 mini_west_border,
                                 mini_north_border,
                                 params->map_ew_res,  
                                 params->map_ns_res,  
                                 params->null_value,
                                 params->antenna_diagram_dir,
                                 tx_params->antenna_diagram_file, 
                                 mini_m_dem,          
                                 mini_m_loss);
#ifdef _PERFORMANCE_METRICS_
    measure_flops ("Antenna", 0);
#endif

    //
    // copy the PathLoss mini matrix back into the big one;
    // we ignore the other ones, because they are used only as input
    //
    mini_r = 0;
    for (r = 0; r < params->nrows; r ++)
    {
        if ((r > mini_north_border_idx) && (r <= mini_south_border_idx))
        {
            mini_c = 0;
            for (c = 0; c < params->ncols; c ++)
            {
                if ((c > mini_west_border_idx) && (c <= mini_east_border_idx))
                {
                    params->m_loss[r][c] = mini_m_loss[mini_r][mini_c];
                    mini_c ++;
                }
                else
                {
                    //
                    // this point is outside the mini matrix (i.e. calculation 
                    // within radius), write null to it
                    //
                    params->m_loss[r][c] = params->null_value;
                }
            }
            mini_r ++;
        }
        else
        {
            //
            // this row is outside the mini matrix (i.e. calculation 
            // within radius), write nulls to it
            //
            for (c = 0; c < params->ncols; c ++)
            {
                params->m_loss[r][c] = params->null_value;
            }
        }
    }
    
    //
    // deallocate the mini matrices
    //
    free (mini_m_dem[0]);
    free (mini_m_dem);
    free (mini_m_clut[0]);
    free (mini_m_clut);
    free (mini_m_loss[0]);
    free (mini_m_loss);

    /*
     * This part is only used when using this module as the objective
     * funcion evaluation component of an optimization process
     *
    //
    // calculate the least square difference between the field
    // measurements and the predicted RSCP values
    //
    if (params->m_field_meas == NULL)
    {
        //
        // allocate memory for the matrix where the measurements
        // will be saved
        //
        int i;

        double *m_field_meas_data = (double *) G_calloc (params->nrows * params->ncols,
                                                         sizeof (double));
        params->m_field_meas = (double **) G_calloc (params->nrows, 
                                                     sizeof (double *));
        for (i = 0; i < params->nrows; i ++)
            params->m_field_meas[i] = &(m_field_meas_data[i * params->ncols]);

        //
        // load the field measurements from the DB
        //
        //load_field_measurements_from_db (params->tx_name,
        //                                 params->tx_east_coord,
        //                                 params->tx_north_coord,
        //                                 params->radius,
        //                                 params->map_west,
        //                                 params->map_north,
        //                                 params->map_ew_res,
        //                                 params->map_ns_res,
        //                                 params->m_field_meas);

        //
        // load the field measurements from a binary file,
        // useful when no database access is available
        //
        char filename [1024];
        strcpy (filename, params->antenna_diagram_dir);
        strcat (filename, "/");
        strcat (filename, params->tx_name);
        strcat (filename, "_field_measurements.bin");
        load_field_measurements_from_file  (filename,
                                            params->m_field_meas);
    }
    *
    */
}


/**
 * Displays the calculation result in the standard output.
 *
 * params           a structure holding configuration parameters which are common
 *                  to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 */
void output_to_stdout (const Parameters *params,
                       const Tx_parameters *tx_params)
{
    int r, c;
    
    // 
    // tell the DB server to get ready before sending the data
    //
    fprintf (stdout, 
             "CREATE TABLE IF NOT EXISTS coverage_%s (east int, north int, rscp float);\n",
             tx_params->tx_name);
    fprintf (stdout,
             "\\COPY coverage_%s (east, north, rscp) FROM STDIN WITH DELIMITER ' '\n",
             tx_params->tx_name);
    //
    // output the data
    //
    for (r = 0; r < params->nrows; r ++)
    {
        for (c = 0; c < params->ncols; c ++)
        {
            float rscp = (float) params->m_loss[r][c];
            if ((!isnan (rscp)) && (rscp != params->null_value))
            {
                float east_coord  = params->map_west + c * params->map_ew_res;
                float north_coord = params->map_north - r * params->map_ns_res;

                fprintf (stdout, "%.f %.f %.5f\n", east_coord,
                                                   north_coord,
                                                   rscp);
            }
        }
    }
}

