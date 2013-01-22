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
    // initialize the structure variables for Ericsson 9999 model
    //
    struct StructEric IniEric = {tx_params->tx_east_coord_idx,
                                 tx_params->tx_north_coord_idx,
                                 tx_params->antenna_height_AGL, 
                                 params->rx_height_AGL,
                                 tx_params->ncols,
                                 tx_params->nrows,
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
                              tx_params->tx_east_coord_idx,
                              tx_params->tx_north_coord_idx,
                              tx_params->antenna_height_AGL,
                              tx_params->total_tx_height,
                              tx_params->beam_direction,
                              tx_params->mechanical_tilt,
                              params->frequency,
                              params->radius,  
                              params->rx_height_AGL,
                              tx_params->nrows,       
                              tx_params->ncols,      
                              tx_params->map_west,
                              tx_params->map_north,
                              params->map_ew_res,  
                              params->map_ns_res,  
                              params->null_value,
                              params->gpu_params,
                              tx_params->m_dem,          
                              tx_params->m_clut,
                              tx_params->m_loss);
    }
    else
    {
#ifdef _PERFORMANCE_METRICS_
        measure_time ("E/// on CPU");
#endif
        EricPathLossSub (tx_params->m_dem, 
                         tx_params->m_clut, 
                         tx_params->m_loss, 
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
                                 tx_params->nrows,       
                                 tx_params->ncols,      
                                 tx_params->map_west,
                                 tx_params->map_north,
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
                                 tx_params->m_dem,          
                                 tx_params->m_loss,
                                 tx_params->m_radio_zone,
                                 tx_params->m_antenna_loss);
#ifdef _PERFORMANCE_METRICS_
    measure_flops ("Antenna", 0);
#endif
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
    for (r = 0; r < tx_params->nrows; r ++)
    {
        for (c = 0; c < tx_params->ncols; c ++)
        {
            float east_coord  = tx_params->map_west + c * params->map_ew_res;
            float north_coord = tx_params->map_north - r * params->map_ns_res;

            float pl = (float) params->tx_params->m_loss[r][c];

            if ((!isnan (pl)) && (pl != params->null_value))
                fprintf (stdout, "%.2f|%.2f|%.5f\n", east_coord,
                                                     north_coord,
                                                     pl);
        }
    }
    //
    // mark end-of transmitter data
    //
    fprintf (stdout, "\\.\n");
}

