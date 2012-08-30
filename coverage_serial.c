#include "coverage.h"
#include "antenna.h"
#include "eric.h"




/**
 * Calculates the area coverage in a serial fashion.
 *
 * params           A structure holding all parameters needed for 
 *                  calculation;
 * eric_params      contains the four tunning parameters for the 
 *                  Ericsson 9999 model, set by the optimization 
 *                  algorithm;
 * eric_params_len  the number of parameters within the received vector,
 *                  four in this case (A0, A1, A2 and A3);
 * output_raster    the name of the output raster created;
 *                  no output is generated if this parameter is NULL.-
 *
 */
void coverage_serial (const Parameters *params,
                      const double *eric_params, 
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
    EricPathLossSub (params->m_rast, params->m_clut, params->m_loss, &IniEric);
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
                                 params->m_rast,
                                 params->m_loss);
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
                path_loss_num = params->m_loss[row][col];
                if (path_loss_num == 0)
                    ((FCELL *) outrast)[col] = null_f_out;
                else
                    ((FCELL *) outrast)[col] = (FCELL)path_loss_num;
            }
            // write raster row to output raster map
            if (G_put_raster_row (outfd, outrast, FCELL_TYPE) < 0)
                G_fatal_error (_("Failed writing raster map <%s>"), output_raster);
        }
        G_close_cell (outfd);
        G_free (outrast);
    }
}

