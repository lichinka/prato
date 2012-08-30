#include "coverage.h"




/**
 * Initializes the coverage calculation by reading the configuration 
 * parameters in the INI file passed as argument.
 * Returns a structure holding all parameters needed for calculation.
 *
 * ini_file     parameter INI file to be used.-
 * tx_section   the section of the INI file containing the
 *              transmitter's configuration parameters.-
 *
 */
static Parameters* coverage_init (const char *ini_file,
                                  const char *tx_section)
{
    //
    // allocate the parameters structure
    //
    Parameters *params = (Parameters *) malloc (sizeof (Parameters));

    //
    // parse the INI file containing the configuration values
    //
    int i = ini_parse (ini_file, 
                       params_handler, 
                       params,
                       tx_section);
    if (i < 0)
        G_fatal_error ("Can't parse '%s'\n", ini_file);

    //
    // set all evaluation parameters up
    //
    char *mapset, *mapset2;		    // mapset names
    int row, col;
    int tr_row, tr_col;
    int infd, infd2;                // file descriptors

    // returns NULL if the map was not found in any mapset
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
    // read metadata of each map, making sure they match
    struct Cell_head *metadata = (struct Cell_head *) malloc (sizeof (struct Cell_head));
    int errno = G_get_cellhd (params->dem_map,
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

    // number of rows and columns within the maps
    params->nrows = metadata->rows;
    params->ncols = metadata->cols;

    // free the allocated map metadata
    free (metadata);

    // check if the specified transmitter location is inside the window
    if (params->tx_east_coord < params->map_west || 
        params->tx_east_coord > params->map_east ||
        params->tx_north_coord > params->map_north || 
        params->tx_north_coord < params->map_south)
        G_fatal_error (_("Specified BS coordinates are outside current region bounds."));
    
    // map array coordinates for transmitter
    tr_row = (params->map_north - params->tx_north_coord) / params->map_ns_res;
    tr_col = (params->tx_east_coord - params->map_west) / params->map_ew_res;

    //
    // allocate the reading buffer for DEM and clutter data
    //
    void *inrast  = G_allocate_raster_buf (FCELL_TYPE);
    void *inrast2 = G_allocate_raster_buf (FCELL_TYPE);

    // total height of transmitter 
    if (G_get_raster_row (infd, inrast, tr_row, FCELL_TYPE) < 0)
        G_fatal_error (_("Unable to read raster map <%s> row %d"), params->dem_map, 
                                                                   tr_row);
    FCELL trans_elev = ((FCELL *) inrast)[tr_col];
    params->total_tx_height = (double)trans_elev + params->antenna_height_AGL;
    
    // check if transmitter is on DEM
    if (isnan ((double) trans_elev))							
        G_fatal_error(_("Transmitter outside raster DEM map."));

    //
    // allocate memory to contain the both DEM and clutter maps
    //
    // DEM
    double *m_rast_data = (double *) G_calloc (params->nrows * params->ncols,
                                               sizeof (double));
    params->m_rast = (double **) G_calloc (params->nrows,
                                           sizeof (double *));
    for (i = 0; i < params->nrows; i ++)
        params->m_rast[i] = &(m_rast_data[i * params->ncols]);

    // CLUTTER
    double *m_clut_data = (double *) G_calloc (params->nrows * params->ncols,
                                               sizeof (double));
    params->m_clut = (double **) G_calloc (params->nrows, 
                                           sizeof (double *));
    for (i = 0; i < params->nrows; i ++)
        params->m_clut[i] = &(m_clut_data[i * params->ncols]);

    //
    // read files (DEM and clutter) into their respective arrays
    //
    for (row = 0; row < params->nrows; row ++) 
    {	
        // read input map
        if (G_get_raster_row (infd, inrast, row, FCELL_TYPE) < 0)
          G_fatal_error(_("Unable to read raster map <%s> row %d"), params->dem_map, row);

        /* read input map 2 */	
        if (G_get_raster_row (infd2, inrast2, row, FCELL_TYPE) < 0)
          G_fatal_error(_("Unable to read raster map <%s> row %d"), params->clutter_map, row);

        /* process the data */
        for (col = 0; col < params->ncols; col ++) 
        { 
            FCELL f_in = ((FCELL *) inrast)[col];
            FCELL f_in2 = ((FCELL *)inrast2)[col];
            params->m_rast[row][col] = (double) f_in;
            params->m_clut[row][col] = (double) f_in2;
        }
    }
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

    //
    // allocate memory for the path-loss matrix
    //
    double *m_loss_data = (double *) G_calloc (params->nrows * params->ncols,
                                               sizeof (double));
    params->m_loss = (double **) G_calloc (params->nrows, 
                                           sizeof (double *));
    for (i = 0; i < params->nrows; i ++)
        params->m_loss[i] = &(m_loss_data[i * params->ncols]);
    
    //
    // return a pointer to the created parameter structure
    //
    return params;
}




/**
 * Entry point of the GRASS module.-
 */
int main (int argc, char *argv [])
{
    struct GModule *module;
    struct Option  *ini_file, *tx_ini_section, *output;
    struct Flag    *use_mpi;

    //
    // initialize the GIS environment
    //
    G_gisinit (argv[0]);

    //
    // initialize the module itself
    //
    module = G_define_module ( );
    module->keywords = _("raster, coverage, E/// 9999, directional antennas");
    module->description = _("Coverage module using the E/// 9999 model and antenna diagrams.-");
  
    //
    // define the module options ...
    //
    ini_file = G_define_option ( );
    ini_file->key = "ini_file";
    ini_file->type = TYPE_STRING;
    ini_file->required = YES;
    ini_file->label = _("Full path to the parameters INI file");	

    tx_ini_section = G_define_option ( );
    tx_ini_section->key = "tx_ini_section";
    tx_ini_section->type = TYPE_STRING;
    tx_ini_section->required = YES;
    tx_ini_section->label = _("Name of the INI-file section containing the transmitter configuration");

    output = G_define_option ( );
    output->key = "output_raster";
    output->type = TYPE_STRING;
    output->required = NO;
    output->label = _("Output raster name");

    use_mpi = G_define_flag ( );
    use_mpi->key = 'p';
    use_mpi->description = _("Whether to use the MPI implementation");

    //
    // ... and parse them
    //
    if (G_parser (argc, argv) < 0)
	    exit (EXIT_FAILURE);

    //
    // initialize the coverage calculation ...
    //
    Parameters *params = coverage_init (ini_file->answer,
                                        tx_ini_section->answer);

    //
    // ... and execute it
    //
    double ericsson_params [4] = {38.0, 32.0, -12.0, 0.1};

#ifdef _PERFORMANCE_METRICS_
    for (i = 0; i < 3; i ++)
#endif
    if (use_mpi->answer)
        coverage_mpi (argc,
                      argv,
                      params,
                      ericsson_params,
                      4,
                      output->answer);
    else
        coverage_serial (params,
                         ericsson_params,
                         4);
    //
    // do we have to write the raster output?
    //
    if (output->answer != NULL)
    {
        int outfd, row, col;
        void *outrast = G_allocate_raster_buf (FCELL_TYPE);    // output buffer
        double path_loss_num;
        FCELL  null_f_out;
        G_set_f_null_value (&null_f_out, 1);   

        // controlling, if we can write the raster 
        if ((outfd = G_open_raster_new (output->answer, FCELL_TYPE)) < 0)
            G_fatal_error(_("Unable to create raster map <%s>"), output->answer);

        for (row = 0; row < params->nrows; row++)
        {
            G_percent(row, params->nrows, 2);
            for (col = 0; col < params->ncols; col++) 
            {
                path_loss_num = params->m_loss[row][col];
                if (path_loss_num == 0)
                    ((FCELL *) outrast)[col] = null_f_out;
                else
                    ((FCELL *) outrast)[col] = (FCELL)path_loss_num;
            }
            // write raster row to output raster map
            if (G_put_raster_row (outfd, outrast, FCELL_TYPE) < 0)
                G_fatal_error (_("Failed writing raster map <%s>"), output->answer);
        }
        G_close_cell (outfd);
        G_free (outrast);
    }

    //
    // deallocate memory before exiting
    //
    free (&(params->m_loss[0][0]));
    free (params->m_loss);
    free (&(params->m_clut[0][0]));
    free (params->m_clut);
    free (&(params->m_rast[0][0]));
    free (params->m_rast);
    free (params);

    exit (EXIT_SUCCESS);
}

