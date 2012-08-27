#include "coverage.h"
#include "antenna.h"
#include "eric.h"



/**
 * Calculates the objective function value, based on the arrays received,
 * and using MPI.-
 *
 * eric_params      Contains the four tunning parameters for the 
 *                  Ericsson 9999 model, set by the optimization 
 *                  algorithm;
 * eric_params_len  the number of parameters within the received vector,
 *                  four in this case (A0, A1, A2 and A3);
 * tx_section_name  the section of the INI file containing the
 *                  transmitter's configuration parameters.-
 *
 *
static void coverage_mpi (const double *eric_params, 
                          const unsigned int eric_params_len,
                          const char *tx_section_name)
{
}
*/



/**
 * Calculates the coverage of the transmitte defined in the INI file,
 * under the section 'tx_section_name'.
 *
 * eric_params      Contains the four tunning parameters for the 
 *                  Ericsson 9999 model, set by the optimization 
 *                  algorithm;
 * eric_params_len  the number of parameters within the received vector,
 *                  four in this case (A0, A1, A2 and A3);
 * output_raster    the name of the output raster created;
 *                  no output is generated if this parameter is NULL.-
 *
 */
static void coverage_serial (const double *eric_params, 
                             const unsigned int eric_params_len,
                             const char *output_raster)
{
    
    //
    // calculate the path-loss over the area, using the Ericsson 9999 model
    //

    // define structure variables
    double BSAntHeight = params->antenna_height_AGL;
    double MSAntHeight = params->rx_height_AGL;
    int xN = params->nrows;
    int yN = params->ncols;
    double scale = params->map_ew_res;
    double A0 = (double) eric_params[0];
    double A1 = (double) eric_params[1];
    double A2 = (double) eric_params[2];
    double A3 = (double) eric_params[3];
    double ResDist = 1;
    //Patrik double BSyIndex = (params->tx_east_coord-params->map_west)/scale+0.5; 
    //Patrik double BSxIndex = (params->map_north-params->tx_north_coord)/scale+0.5;
    double BSyIndex = (params->tx_east_coord - params->map_west) / scale; 
    double BSxIndex = (params->map_north - params->tx_north_coord) / scale;
    double radi = params->radius;

    //G_message(_("2.Parameter BSxIndex: %f, BSyIndex: %f, BSAntHeight: %f, MSAntHeight: %f, xN: %d, yN: %d, scale: %f, freq: %f, A0: %f, A1: %f, A2: %f, A3: %f, ResDist: %f, radi: %f"),BSxIndex,BSyIndex, BSAntHeight, MSAntHeight,xN, yN, scale, freq, A0, A1, A2, A3, ResDist, radi);

    struct StructEric IniEric = {BSxIndex, BSyIndex, 
                                 BSAntHeight, MSAntHeight,
                                 xN, yN, 
                                 scale, params->tx_frequency, 
                                 A0, A1, A2, A3, 
                                 ResDist, radi};
#ifdef _PERFORMANCE_METRICS_
    measure_flops ("Ericsson", 1);
#endif
    EricPathLossSub (m_rast, m_clut, m_loss, &IniEric);
#ifdef _PERFORMANCE_METRICS_
    measure_flops ("Ericsson", 0);
#endif

    //
    // calculate the antenna influence, 
    // overwriting the isotrophic path-loss
    //
#ifdef _PERFORMANCE_METRICS_
    measure_flops ("Antenna", 1);
#endif
    calculate_antenna_influence (params->antenna_diagram_dir,
                                 params->antenna_diagram_file,
                                 params->tx_east_coord,
                                 params->tx_north_coord,
                                 params->total_tx_height,
                                 params->beam_direction,
                                 params->mechanical_tilt,
                                 params->radius,
                                 params->rx_height_AGL,
                                 params->nrows,
                                 params->ncols,
                                 params->map_west,
                                 params->map_north,
                                 params->map_ew_res,
                                 params->map_ns_res,
                                 m_rast,
                                 m_loss);
#ifdef _PERFORMANCE_METRICS_
    measure_flops ("Antenna", 0);
#endif

    /*
     * This part is only used when using this module as the objective
     * funcion evaluation component of an optimization process
     *
    //
    // calculate the least square difference between the field
    // measurements and the predicted RSCP values
    //
    if (m_field_meas == NULL)
    {
        //
        // allocate memory for the matrix where the measurements
        // will be saved
        //
        int i;
        m_field_meas = (double **) G_calloc (params->nrows, 
                                             sizeof (double *));
        for (i = 0; i < params->nrows; i ++)
            m_field_meas[i] = (double *) G_calloc (params->ncols,
                                                   sizeof (double));
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
        //                                 m_field_meas);

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
                                            m_field_meas);
    }
    *
    */

    //
    // do we have to write the raster output?
    //
    if (output_raster != NULL)
    {
        int outfd, row, col;
        void *outrast = G_allocate_raster_buf (FCELL_TYPE);    // output buffer
        double path_loss_num;
        FCELL  null_f_out;
        G_set_f_null_value (&null_f_out, 1);   

        // controlling, if we can write the raster 
        if ((outfd = G_open_raster_new (output_raster, FCELL_TYPE)) < 0)
            G_fatal_error(_("Unable to create raster map <%s>"), output_raster);

        for (row = 0; row < params->nrows; row++)
        {
            G_percent(row, params->nrows, 2);
            for (col = 0; col < params->ncols; col++) 
            {
                path_loss_num = m_loss[row][col];
                if (path_loss_num == 0)
                    ((FCELL *) outrast)[col] = null_f_out;
                else
                    ((FCELL *) outrast)[col] = (FCELL)path_loss_num;
            }
            // write raster row to output raster map
            if (G_put_raster_row(outfd, outrast, FCELL_TYPE) < 0)
                G_fatal_error(_("Failed writing raster map <%s>"), output_raster);
        }
        G_close_cell (outfd);
        G_free (outrast);
    }
}




/**
 * Initializes the coverage calculation by reading the configuration 
 * parameters in the INI file passed as argument.
 *
 * ini_file     Parameter INI file to be used.-
 * tx_section   the section of the INI file containing the
 *              transmitter's configuration parameters.-
 *
 */
static void coverage_init (const char *ini_file,
                           const char *tx_section)
{
    //
    // check that the parameters object has yet to be initialized
    //
    if (params == NULL)
    {
        //
        // allocate the parameters structure
        //
        params = (Parameters *) malloc (sizeof (Parameters));
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
        m_rast = (double **) G_calloc (params->nrows,
                                       sizeof (double *));
        for (i = 0; i < params->nrows; i ++)
            m_rast[i] = (double *) G_calloc (params->ncols, 
                                             sizeof (double));

        // CLUTTER
        m_clut = (double **) G_calloc (params->nrows, sizeof (double *));
        for (i = 0; i < params->nrows; i ++)
            m_clut[i] = (double *) G_calloc (params->ncols, 
                                             sizeof (double));

        //
        // read files (DEM and clutter) into their respective arrays
        //
        for (row = 0; row < params->nrows; row ++) 
        {	
            /* read input map */
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
                m_rast[row][col] = (double) f_in;
                m_clut[row][col] = (double) f_in2;
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
    }
    //
    // allocate memory for the path-loss matrix only once
    //
    if (m_loss == NULL)
    {
        int i;
        m_loss = (double **) G_calloc (params->nrows, 
                                       sizeof (double *));
        for (i = 0; i < params->nrows; i ++)
            m_loss[i] = (double *) G_calloc (params->ncols, 
                                             sizeof (double));
    }
}




/**
 * Entry point of the GRASS module.-
 */
int main (int argc, char *argv [])
{
    struct GModule *module;
    struct Option *ini_file, *tx_ini_section, *output;

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

    //
    // ... and parse them
    //
    if (G_parser (argc, argv) < 0)
	    exit (EXIT_FAILURE);

    //
    // initialize the coverage calculation ...
    //
    coverage_init (ini_file->answer,
                   tx_ini_section->answer);

    //
    // ... and execute it
    //
    double ericsson_params [4] = {38.0, 32.0, -12.0, 0.1}; 
    coverage_serial (ericsson_params,
                     4,
                     output->answer);
    //
    // deallocate memory before exiting
    //
    int i;
    for (i = 0; i < params->nrows; i ++)
    {
        free (m_loss[i]);
        free (m_clut[i]);
        free (m_rast[i]);
    }
    free (m_loss);
    free (m_clut);
    free (m_rast);
    free (params);

    exit (EXIT_SUCCESS);
}

