#include "worker/coverage.h"
#include "worker/antenna.h"
#include "worker/eric.h"



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
void 
coverage (const Parameters     *params,
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
    // radio zones
    //
    char *mini_m_radio_zone_data = (char *) calloc (diameter_in_pixels * diameter_in_pixels, 
                                                    sizeof (char));
    char **mini_m_radio_zone = (char **) calloc (diameter_in_pixels,
                                                 sizeof (char *));
    for (r = 0; r < diameter_in_pixels; r ++)
        mini_m_radio_zone[r] = &(mini_m_radio_zone_data[r * diameter_in_pixels]);
  
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
                    mini_m_dem[mini_r][mini_c]        = params->m_dem[r][c];
                    mini_m_clut[mini_r][mini_c]       = params->m_clut[r][c];
                    mini_m_loss[mini_r][mini_c]       = params->m_loss[r][c];
                    mini_m_radio_zone[mini_r][mini_c] = params->tx_params->m_radio_zone[r][c];
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
#ifdef _PERFORMANCE_METRICS_
        measure_time ("E/// on GPU");
#endif
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
                              params->gpu_params,
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
    }
#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
#endif

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
                                 params->main_zone_horiz,
                                 params->main_zone_vert,
                                 params->sec_zone_horiz,
                                 params->sec_zone_vert,
                                 params->gpu_params,
                                 mini_m_dem,          
                                 mini_m_loss,
                                 mini_m_radio_zone);
#ifdef _PERFORMANCE_METRICS_
    measure_flops ("Antenna", 0);
#endif

    //
    // copy the path-loss and radio-zones mini matrices back into the big ones;
    // we ignore the other ones, because they are only used as input
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
                    params->m_loss[r][c]                  = mini_m_loss[mini_r][mini_c];
                    params->tx_params->m_radio_zone[r][c] = mini_m_radio_zone[mini_r][mini_c];
                    mini_c ++;
                }
                else
                {
                    //
                    // this point is outside the calculation radius
                    //
                    params->m_loss[r][c]                   = params->null_value;
                    params->tx_params->m_radio_zone[r][c] &= _RADIO_ZONE_MODEL_DISTANCE_OFF_;
                }
            }
            mini_r ++;
        }
        else
        {
            //
            // this row is outside the calculation radius
            //
            for (c = 0; c < params->ncols; c ++)
            {
                params->m_loss[r][c]                   = params->null_value;
                params->tx_params->m_radio_zone[r][c] &= _RADIO_ZONE_MODEL_DISTANCE_OFF_;
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
    free (mini_m_radio_zone[0]);
    free (mini_m_radio_zone);

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
void 
output_to_stdout (const Parameters *params,
                  const Tx_parameters *tx_params)
{
    int r, c;
    
    // 
    // prepare the DB server before sending the data
    //
    fprintf (stdout, 
             "CREATE TABLE IF NOT EXISTS coverage_%s (east float, north float, pl float);\n",
             tx_params->tx_name);
    fprintf (stdout, 
             "TRUNCATE TABLE coverage_%s;\n",
             tx_params->tx_name);
    fprintf (stdout,
             "\\COPY coverage_%s (east, north, pl) FROM STDIN WITH DELIMITER '|'\n",
             tx_params->tx_name);
    //
    // output the data
    //
    for (r = 0; r < params->nrows; r ++)
    {
        for (c = 0; c < params->ncols; c ++)
        {
            if (params->use_opt)
            {
                char  rz = params->tx_params->m_radio_zone[r][c];
                float east_coord  = params->map_west + c * params->map_ew_res;
                float north_coord = params->map_north - r * params->map_ns_res;
                if (rz)
                    fprintf (stdout, "%.2f|%.2f|%02x\n", east_coord,
                                                         north_coord,
                                                         rz);
            }
            else
            {
                float pl = (float) params->m_loss[r][c];

                if ((!isnan (pl)) && (pl != params->null_value))
                {
                    float east_coord  = params->map_west + c * params->map_ew_res;
                    float north_coord = params->map_north - r * params->map_ns_res;

                    fprintf (stdout, "%.2f|%.2f|%.5f\n", east_coord,
                                                         north_coord,
                                                         pl);
                }
            }
        }
    }
    //
    // mark end-of transmitter data
    //
    fprintf (stdout, "\\.\n");
}

