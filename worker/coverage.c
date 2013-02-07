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
coverage (Parameters           *params,
          Tx_parameters        *tx_params,
          const double         *eric_params, 
          const unsigned int   eric_params_len)
{
    //
    // execute the path-loss calculation on CPU or GPU?
    //
    if (params->use_gpu)
    {
#ifdef _PERFORMANCE_METRICS_
        measure_time ("E/// on GPU");
#endif
        eric_pathloss_on_gpu (params,
                              tx_params,
                              eric_params);
    }
    else
    {
#ifdef _PERFORMANCE_METRICS_
        measure_time ("E/// on CPU");
#endif
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
        EricPathLossSub (tx_params->m_obst_height,
                         tx_params->m_obst_dist,
                         tx_params->m_obst_offset,
                         tx_params->m_dem, 
                         tx_params->m_clut, 
                         tx_params->m_loss, 
                         tx_params->m_field_meas,
                         tx_params->m_antenna_loss,
                         tx_params->m_radio_zone,
                         &IniEric);
    }
#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
#endif

    //
    // calculate the antenna influence, 
    // overwriting the isotrophic path-loss
    //
#ifdef _PERFORMANCE_METRICS_
    measure_time ("Antenna influence");
#endif
    calculate_antenna_influence (params,
                                 tx_params);
#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
#endif

    /*
    // DEBUG memory
    //
    int r, c;
    for (r = 0; r < tx_params->nrows; r ++)
    {
        for (c = 0; c < tx_params->ncols; c++)
        {
            fprintf (stdout,
                     "%.5f\t", tx_params->m_antenna_loss[r][c]);
        }
        fprintf (stdout, "\n");
    }
    fflush (stdout);
    exit (-1);
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
    // bring the data from the GPU if it has been used
    //
    if (params->use_gpu)
    {
        size_t buff_size = tx_params->nrows * 
                           tx_params->ncols * 
                           sizeof (tx_params->m_loss[0][0]);
        read_buffer_blocking (tx_params->ocl_obj,
                              0,
                              tx_params->m_loss_dev,
                              buff_size,
                              tx_params->m_loss[0]);
    }

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

