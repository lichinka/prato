
/****************************************************************************
 *
 * MODULE:      r.ericsson
 * AUTHOR(S):   Andrej Vilhar, Tomaz Javornik, Andrej Hrovat Jozef Stefan Institute
 * UPDATE:		Patrik Ritosa, Telekom Slovenije                
 *
 * PURPOSE:     Calculates radio coverage from a single base station 
 *              according to Ericsson model
 *             
 *
 * COPYRIGHT:   (C) 2009 Jozef Stefan Institute
 *
 *
 *****************************************************************************/
/****************************************************************************
 *
 * MODULE:       do_profile function for r.ericsson
 * AUTHOR(S):    Patrik Ritosa, Telekom Slovenije                
 *
 * PURPOSE:      Calculates obstacle high and distance betwen BS and MS
 *             
 *
 * COPYRIGHT:
 *
 *
 *****************************************************************************/
/****************************************************************************
 *
 * MODULE:       E/// on GPU and parameter_fine_tuning for E///
 * AUTHOR(S):    Lucas Benedicic, Telekom Slovenije                
 *
 * PURPOSE:      Fine tunes the various parameters of the prediction model in 
 *               order to minimize the mean square error with field 
 *               measurements
 *
 * COPYRIGHT:    GNU GPL
 *
 *
 *****************************************************************************/

//
// From the GSL library
//
//      http://www.gnu.org/software/gsl/
//
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>

//
// From the performance metrics library
//  
//      https://github.com/lichinka/performance_metrics
//
#include <performance_metric.h>

#include "worker/eric.h"
#include "worker/antenna.h"




/**
 * Calculates the path loss (in dB) on a specific coordinate, 
 * using E/// 9999 path-loss formula.
 *			
 * Ericsson 9999 model update (Patrik Ritosa, May 2012):
 * 			- Spherical earth diffraction komponent (JDFR)
 * 			- NLOS komponent (KDFR) 
 * 			- do_profile calculation speed optimisation
 * 
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 * ix               northern-coordinate index within the area;
 * iy               eastern-coordinate index within the area;
   DistBS2MSKm      distance between MS and BS in km;
 * log10Zeff        logarithm of the effective antenna height (output parameter);
 * nlos             NLOS component of the path-loss formula (output parameter).-
 *
 */
static void
eric_pathloss_on_point (const Parameters    *params,
                        const Tx_parameters *tx_params,
                        const int           ix,
                        const int           iy,
                        const double        DistBS2MSKm,
                        double             *log10Zeff,
                        double             *nlos)
{
	int BSxIndex       = tx_params->tx_north_coord_idx;	// position of BS in pixels
	int BSyIndex       = tx_params->tx_east_coord_idx;	// position of BS in pixels
	double AntHeightBS = tx_params->antenna_height_AGL;	// antenna height of BS [m]
	double AntHeightMS = params->rx_height_AGL;	        // antenna height of MS [m]
	double scale       = params->map_ew_res;            // terrain resolution (in mts)
	double A0          = tx_params->eric_params[0];		// model parameters
	double A1          = tx_params->eric_params[1];		// model parameters
	double A2          = tx_params->eric_params[2];		// model parameters
	double A3          = tx_params->eric_params[3];	    // model parameters
	double freq        = params->frequency;             // carrier frequency
	double Lambda      = 300.0 / freq;			        //	wave lenght

	double ZoBS;					// BS and MS height about the sea level
	double ZObs2LOS = 0;
	double DistObs2BS = 0;
	double ZoTransBS,ZoTransMS;
    double tiltBS2MS;				// (ZoBS-ZoMS)/distBS2MSNorm	
	double PathLossFreq = 0;	    // path loss due to carrier frequency
	double PathLossTmp = 0;			// tmp path loss
    double Zeff;			        // Difference in X and Y direction
	double PathLossAntHeightMS;
    double DistBS2MSNorm;           // normalized distance between MS and BS
							
	double ElevAngCos, Hdot, Ddot, Ddotdot, PathLossDiff;
    double KDFR, Alfa, Fresnel, JDFR;

	// BS height above the sea level calculated from raster DEM file
	ZoBS = tx_params->m_dem[(int)BSxIndex][(int)BSyIndex];

	// BS transmitter height above the sea level
	ZoTransBS = ZoBS + AntHeightBS;

	PathLossFreq = 44.49*log10(freq) - 4.78*pow(log10(freq),2);	// Loss due to carrier frequency
									/*POPRAVJNEO (4.2.2010)*/	
	PathLossAntHeightMS = 3.2*pow(log10(11.75*AntHeightMS),2);

    //
    // Path Loss due to Hata model
    //
    // ZoMS = tx_params->m_dem[ix][iy];
    ZoTransMS = tx_params->m_dem[ix][iy]+AntHeightMS;  // ZoMS
    Zeff = ZoTransBS - ZoTransMS;		// ??
    DistBS2MSNorm = (DistBS2MSKm * 1000) / scale;

    //height correction due to earth sphere
    Zeff = Zeff + (DistBS2MSKm*DistBS2MSKm)/((6370 * 8000) / 3);

    if ((- AntHeightBS < Zeff) && (Zeff < AntHeightBS))
    {
        Zeff = AntHeightBS;		// Preventing Log10(Zeff) to go toward -inf
    }
    
    *log10Zeff = log10 (fabs (Zeff));

    double log10DistBS2MSKm = log10(DistBS2MSKm);			
    
    PathLossTmp = A0 + A1 * (log10DistBS2MSKm); 
    PathLossTmp = PathLossTmp + A2 * (*log10Zeff);
    PathLossTmp = PathLossTmp + A3 * (log10DistBS2MSKm) * (*log10Zeff);
    PathLossTmp = PathLossTmp - PathLossAntHeightMS + PathLossFreq;

    //
    // Calc position of the height and position of the highest obstacle
    // 
    tiltBS2MS = ZoTransBS - ZoTransMS; 	//STARO: tiltBS2MS = Zeff; Zeff je vmes lahko spremenjena /* Sprememba (4.2.2010)*/

    if (DistBS2MSNorm > 0) 
        tiltBS2MS = -tiltBS2MS/DistBS2MSNorm; 
    else
        tiltBS2MS = 0; 
    
    ZObs2LOS = tx_params->m_obst_height[ix][iy];
    DistObs2BS = tx_params->m_obst_dist[ix][iy];
    // Calc path loss due to NLOS conditions
/*Patrik/scale*/	ElevAngCos = cos(atan(tiltBS2MS/scale));

    Ddot = DistObs2BS; 
    if (ElevAngCos != 0) 
    {	
        Ddot = DistObs2BS/ElevAngCos;
    }
    Ddotdot = DistBS2MSNorm - Ddot;
    if (ElevAngCos != 0) {
        Ddotdot = (DistBS2MSNorm)/ElevAngCos - Ddot;
    }

// Obstacle height korrection due to earth sphere
    if (Ddot <= Ddotdot){
        ZObs2LOS = ZObs2LOS + (Ddot*scale/1000*Ddot*scale/1000)/(2 * 4/3 * 6370)*1000;
    }
    else{
        ZObs2LOS = ZObs2LOS + (Ddotdot*scale/1000*Ddotdot*scale/1000)/(2 * 4/3 * 6370)*1000;
    }

//Hight correction due to BS2MS line angle
    Hdot = ZObs2LOS*ElevAngCos;

    PathLossDiff = 0;
    KDFR = 0;

    if (Ddot > 0 && Ddotdot > 0) {
        Fresnel=sqrt((Lambda*Ddot*Ddotdot*scale)/(Ddot+Ddotdot)); // First Fresnel elipsoid radius

        PathLossDiff = Hdot/Fresnel;

// NLOS komponent calculation KDFR

        if (PathLossDiff < -0.49 ) {
            KDFR = 0; 
        }
        else if(-0.49 <= PathLossDiff && PathLossDiff < 0.5) {
            KDFR = 6 + 12.2 * PathLossDiff;
        }
        else if(0.5 <= PathLossDiff && PathLossDiff < 2) {
            KDFR = 11.61 * PathLossDiff - 2 * PathLossDiff *PathLossDiff + 6.8;
        }
        else if(2 <= PathLossDiff) {
            KDFR = 16 + 20 * log10(PathLossDiff);
        }

// Alfa correction factor
        
        if (Hdot > Fresnel/2)
        {
            Alfa = 1;
        }
        else if (Fresnel/4 <= Hdot && Hdot <= Fresnel/2)
        {
            Alfa = 4 * Hdot/Fresnel - 1;
        }
        else if (Hdot < Fresnel/4)
        {
            Alfa = 0;
        }
    }

//Spherical earth diffraction komponent JDFR
    if (Hdot > 0)
    {
        JDFR = 20 + 0.112 * pow(freq/(16/9),1/3) * (DistBS2MSKm - sqrt(12.73 * 4/3) * ( sqrt(ZoTransBS) + sqrt(ZoTransMS) )  );
        if (JDFR < 0)
            JDFR=0;
    }
    else
    {
        JDFR = 0;
    }
    //
    // NLOS path-loss component
    //
    *nlos = sqrt(pow(Alfa*KDFR,2) + pow(JDFR,2));

    PathLossTmp += *nlos;

    //
    // get the clutter loss, based on the category of this point
    //
    double clutter_loss = 0;
    int    clutter_category = (int) tx_params->m_clut[ix][iy];

    if (clutter_category == params->cell_null_value)
    {
        fprintf (stderr,
                 "*** WARNING: clutter category (%d) at (%.f, %.f) does not exist. Assuming 0 (zero).\n",
                 clutter_category,
                 tx_params->map_west  + (iy * params->map_ns_res),
                 tx_params->map_north - (ix * params->map_ew_res));
    }
    else
    {
        clutter_loss = params->clutter_loss[clutter_category];
    }
    tx_params->m_loss[ix][iy] = PathLossTmp + clutter_loss;

#ifdef _DEBUG_INFO_
    //
    // DEBUG: parameter dump for external approximation using least squares
    //
    if ((ix == 0) && (iy == 0))
    {
        // column titles: only at the beggining
        printf ("xi|yi|log(d)|HEBK|log(HEBK)|NLOS|A0|A1|A2|A3|-k0+k1|clut|path_loss|field_meas|antenna\n");
    }
    //if ((params->use_opt) && (! isnan (tx_params->m_field_meas[ix][iy])))
    if (! isnan (tx_params->m_clut[ix][iy]))
        if ((tx_params->m_radio_zone[ix][iy] & _RADIO_ZONE_MAIN_BEAM_ON_) > 0)
            printf ("%d|%d|%20.10f|%20.10f|%20.10f|%20.10f|%20.10f|%20.10f|%20.10f|%20.10f|%20.10f|%20.10f|%20.10f|%20.10f\n",
                    ix,
                    iy,
                    log10DistBS2MSKm,
                    fabs (Zeff),
                    *log10Zeff,
                    *nlos,
                    A0,
                    A1,
                    A2,
                    A3,
                    -PathLossAntHeightMS+PathLossFreq,
                    tx_params->m_clut[ix][iy],
                    tx_params->m_loss[ix][iy],
                    //tx_params->m_field_meas[ix][iy],
                    tx_params->m_antenna_loss[ix][iy]);
#endif
}



/**
 * Calculates the path loss using E/// 9999 model implementation on GPU.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters.-
 *
 */
void 
eric_pathloss_on_gpu (Parameters    *params,
                      Tx_parameters *tx_params)
{
    //
    // activate the compiled kernel
    //
    activate_kernel (tx_params->ocl_obj,
                     "eric_per_tx");
    //
    // send the clutter losses to the device
    //
    write_buffer_blocking (tx_params->ocl_obj,
                           0,
                           tx_params->v_clutter_loss_dev,
                           params->clutter_category_count * sizeof (params->clutter_loss[0]),
                           params->clutter_loss);
    //
    // set scalar kernel parameters 
    // 
    set_kernel_double_arg (tx_params->ocl_obj,
                           0,
                           &params->map_ew_res);
    set_kernel_double_arg (tx_params->ocl_obj,
                           1,
                           &tx_params->map_north);
    set_kernel_double_arg (tx_params->ocl_obj,
                           2,
                           &tx_params->map_west);
    set_kernel_double_arg (tx_params->ocl_obj,
                           5,
                           &params->rx_height_AGL);
    set_kernel_double_arg (tx_params->ocl_obj,
                           6,
                           &params->frequency);
    //
    // set prediction model parameters
    //
    cl_double4 model_params;
    model_params.s[0] = tx_params->eric_params[0];
    model_params.s[1] = tx_params->eric_params[1];
    model_params.s[2] = tx_params->eric_params[2];
    model_params.s[3] = tx_params->eric_params[3];
    set_kernel_value_arg (tx_params->ocl_obj,
                          7,
                          sizeof (cl_double4),
                          &model_params);
    //
    // set memory pointer parameters
    //
    set_kernel_mem_arg (tx_params->ocl_obj,
                        8,
                        tx_params->m_obst_height_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        9,
                        tx_params->m_obst_dist_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        10,
                        tx_params->m_dem_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        11,
                        tx_params->m_clut_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        12,
                        tx_params->v_clutter_loss_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        13,
                        tx_params->m_loss_dev);

    // reserve local memory on the device
    size_t lmem_size = _WORK_ITEMS_PER_DIMENSION_ *
                       _WORK_ITEMS_PER_DIMENSION_ *
                       sizeof (tx_params->m_loss[0][0]);
    set_local_mem (tx_params->ocl_obj,
                   14,
                   lmem_size);
    //
    // calculation radius in meters
    //
    double radius_in_meters = params->radius * 1000;

    //
    // calculation tile offset within the target area, given in pixel coordinates
    //
    cl_int2 tile_offset;
    tile_offset.s[0]  = (int) ((tx_params->tx_east_coord - tx_params->map_west) - radius_in_meters);
    tile_offset.s[0] /= params->map_ew_res;
    tile_offset.s[1]  = (int) ((tx_params->map_north - tx_params->tx_north_coord) - radius_in_meters);
    tile_offset.s[1] /= params->map_ns_res;

    //
    // transmitter data
    //
    cl_double4 tx_data;
    tx_data.s[0] = (double) tx_params->tx_east_coord;   // transmitter coordinate
    tx_data.s[1] = (double) tx_params->tx_north_coord;  // transmitter coordinate
    tx_data.s[2] = (double) tx_params->total_tx_height; // antenna height above sea level
    tx_data.s[3] = (double) tx_params->antenna_height_AGL; // antenna height above ground

    //
    // set kernel parameters, specific for this transmitter
    //
    set_kernel_value_arg (tx_params->ocl_obj,
                          3,
                          sizeof (cl_double4),
                          &tx_data);
    set_kernel_value_arg (tx_params->ocl_obj,
                          4,
                          sizeof (cl_int2),
                          &tile_offset);
    //
    // define a 2D execution range for the kernel ...
    //
    size_t global_sizes [2],
           local_sizes  [2];
    define_2D_range (params,
                     global_sizes,
                     local_sizes);
    //
    // ... and execute it
    //
    run_kernel_2D_blocking (tx_params->ocl_obj,
                            0,
                            NULL,
                            global_sizes,
                            local_sizes);
    //
    // no need to bring the path-loss matrix from the device to the host
    //
    /*
    size_t buff_size = tx_params->nrows * 
                       tx_params->ncols * 
                       sizeof (tx_params->m_loss[0][0]);
    read_buffer_blocking (tx_params->ocl_obj,
                          0,
                          tx_params->m_loss_dev,
                          buff_size,
                          tx_params->m_loss[0]);
     */
}



/**
 * Calculates path loss using the E/// 9999 prediction model on CPU.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters.-
 *
 */
void 
eric_pathloss_on_cpu (Parameters    *params,
                      Tx_parameters *tx_params)
{
	int BSxIndex       = tx_params->tx_north_coord_idx;	// position of BS in pixels
	int BSyIndex       = tx_params->tx_east_coord_idx;	// position of BS in pixels
	int xN             = tx_params->nrows;				// dimension of the input(Raster) and output (PathLoss)
	int yN             = tx_params->ncols;				// dimension of the input(Raster) and output (PathLoss)
	double scale       = params->map_ew_res;            // terrain resolution (in mts)
	double radi        = params->radius;			    // calculation radius around transmiter

	int     ix, iy;	
	int     DiffX, DiffY;
    double  nlos;
	double  log10Zeff;
    double  DistBS2MSKm;	// distance between MS and BS in km

#ifdef _PERFORMANCE_METRICS_
    measure_time ("E/// on CPU");
#endif
	for (ix = 0; ix < xN; ix++)
	{
		for (iy = 0; iy < yN; iy++)
		{
            //
            // distance from the current point to the transmitter
            //
			DiffX = (BSxIndex-ix); 
            DiffY = (BSyIndex-iy);
			DistBS2MSKm = (sqrt (DiffX*DiffX + DiffY*DiffY) * scale) / 1000;

			if (DistBS2MSKm < 0.01)
				DistBS2MSKm = 0.01;
                
            //
            // calculate path loss for points within the user-defined radius
            //
			if (DistBS2MSKm <= radi)
            {
                eric_pathloss_on_point (params,
                                        tx_params,
                                        ix,
                                        iy,
                                        DistBS2MSKm,
                                        &log10Zeff,
                                        &nlos);
            }
        }
    }
#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
#endif
}



/**
 * Fine tunes the A0...A3 E/// parameters to best fit a set of measurements 
 * within a given radio zone. The result is saved in the vector
 * `tx_params->eric_params`.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters.-
 *
 */
void
parameter_fine_tuning (Parameters    *params,
                       Tx_parameters *tx_params)
{
	int BSxIndex       = tx_params->tx_north_coord_idx;	// position of BS in pixels
	int BSyIndex       = tx_params->tx_east_coord_idx;	// position of BS in pixels
	double AntHeightMS = params->rx_height_AGL;	        // antenna height of MS [m]
	int xN             = tx_params->nrows;				// dimension of the input(Raster) and output (PathLoss)
	int yN             = tx_params->ncols;				// dimension of the input(Raster) and output (PathLoss)
	double scale       = params->map_ew_res;            // terrain resolution (in mts)
	double freq        = params->frequency;             // carrier frequency
	double radi        = params->radius;			    // calculation radius around transmiter

    double nlos;
	double log10Ha;
	double log10DistBS2MSKm;
	double PathLossFreq = 0;	                        // path loss due to carrier frequency
	int ix; int iy;	
	int DiffX, DiffY;
	double PathLossAntHeightMS;
    double DistBS2MSKm;		                            // distance between MS and BS in km

    //
    // define and initialize matrices for solving the linear system
    // of equations that will return the best-fiting parameters
    // for the prediction model (a0 and a1 only)
    //
    //  A x = b
    //
    const int Dim = 2;
    gsl_matrix *A = gsl_matrix_alloc (Dim, Dim);
    gsl_vector *b = gsl_vector_alloc (Dim);
    double A_data [Dim] [Dim];
    double b_data [Dim];

    gsl_matrix_set_zero (A);
    gsl_vector_set_zero (b);

    for (ix = 0; ix < Dim; ix ++)
    {
        for (iy = 0; iy < Dim; iy ++)
            A_data[ix][iy] = 0;
        b_data[ix] = 0;
    }

	PathLossFreq = 44.49*log10(freq) - 4.78*pow(log10(freq),2);	// Loss due to carrier frequency
									/*POPRAVJNEO (4.2.2010)*/	
	PathLossAntHeightMS = 3.2*pow(log10(11.75*AntHeightMS),2);

    //
    // these parameters (a2 and a3) are linear combination of the other two (a0 and a1)
    //
    tx_params->eric_params[2] = -12.0;
    tx_params->eric_params[3] = 0.1;

    //
    // reset the counter of valid field measurements
    //
    tx_params->field_meas_count = 0;

#ifdef _PERFORMANCE_METRICS_
    measure_time ("E/// parameter fine tuning on CPU");
#endif
    //
    // DEBUG
    //
    printf ("## pwr\trcv\tlog10Ha\tpl_Ha\tpl_freq\tclut\tnlos\tlog10D\tpl_ant\n");
	for (ix = 0; ix < xN; ix++)
	{
		for (iy = 0; iy < yN; iy++)
		{
            //
            // distance from the current point to the transmitter
            //
			DiffX = (BSxIndex-ix); 
            DiffY = (BSyIndex-iy);
			DistBS2MSKm = (sqrt (DiffX*DiffX + DiffY*DiffY) * scale) / 1000;

			if (DistBS2MSKm < 0.01)
				DistBS2MSKm = 0.01;
                
            //
            // calculate path loss for points within the user-defined radius
            //
			if (DistBS2MSKm <= radi)
		   	{    
                //
                // accumulate the elements for solving the linear system of equations,
                // only where there are valid field measurements
                //
                if (! isnan (tx_params->m_field_meas[ix][iy]))
                {
                    eric_pathloss_on_point (params,
                                            tx_params,
                                            ix,
                                            iy,
                                            DistBS2MSKm,
                                            &log10Ha,
                                            &nlos);
                    log10DistBS2MSKm = log10 (DistBS2MSKm);

                    //
                    // constant part of the Hata Open Area formula
                    //
                    double HOA_k = tx_params->eric_params[2] * log10Ha +
                                   tx_params->eric_params[3] * log10DistBS2MSKm * log10Ha -
                                   PathLossAntHeightMS +
                                   PathLossFreq;
                    //
                    // valid field measurement and prediction within the user-defined radio zone
                    //
                    tx_params->field_meas_count ++;

                    //
                    // A matrix - element accumulation
                    //
                    A_data[0][1] += log10DistBS2MSKm;

                    A_data[1][0] += log10DistBS2MSKm;
                    A_data[1][1] += log10DistBS2MSKm * log10DistBS2MSKm;

                    //
                    // b matrix - element accumulation
                    //
                    int clutter_category = (int) tx_params->m_clut[ix][iy];

                    b_data[0] += tx_params->tx_power - 
                                 tx_params->m_field_meas[ix][iy] -
                                 HOA_k -
                                 params->clutter_loss[clutter_category] - 
                                 tx_params->m_antenna_loss[ix][iy] -
                                 nlos;
                    b_data[1] += (tx_params->tx_power - 
                                 tx_params->m_field_meas[ix][iy] -
                                 HOA_k -
                                 params->clutter_loss[clutter_category] - 
                                 tx_params->m_antenna_loss[ix][iy] -
                                 nlos) * log10DistBS2MSKm;
                    //
                    // DEBUG
                    //
                    printf ("#%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\t%.5f\n", tx_params->tx_power,
                                                                                       tx_params->m_field_meas[ix][iy],
                                                                                       log10Ha,
                                                                                       PathLossAntHeightMS,
                                                                                       PathLossFreq,
                                                                                       params->clutter_loss[clutter_category],
                                                                                       nlos,
                                                                                       log10DistBS2MSKm,
                                                                                       tx_params->m_antenna_loss[ix][iy]);
                }
            }
		}
	}
    //
    // the only matrix element that does not need accumulation
    //
    A_data[0][0] = tx_params->field_meas_count;
    
    //
    // use the GSL library to solve the linear system of equations
    //
    int s;
	double det;
    gsl_vector      *x  = gsl_vector_alloc (Dim);
    gsl_vector      *r  = gsl_vector_alloc (Dim);
    gsl_permutation *p  = gsl_permutation_alloc (Dim);

    //
    // copy the accumulated sums into the matrices
    //
    for (ix = 0; ix < Dim; ix ++)
    {
        for (iy = 0; iy < Dim; iy ++)
            gsl_matrix_set (A, ix, iy, A_data[ix][iy]);
        gsl_vector_set (b, ix, b_data[ix]);
    }
	//
	// if the matrix is singular (i.e. its determinant is 0), the system is not solvable
	//
    gsl_linalg_LU_decomp (A, p, &s);
	det = gsl_linalg_LU_det (A, s);
	if (det != 0)
	{
		//
		// solve the linear system
		//
		gsl_linalg_LU_solve  (A, p, b, x);

		//
		// copy the original values back into A, i.e.
		// the values before the decomposition
		//
		for (ix = 0; ix < Dim; ix ++)
			for (iy = 0; iy < Dim; iy ++)
				gsl_matrix_set (A, ix, iy, A_data[ix][iy]);

		//
		// r = A * x - b
		//
		gsl_blas_dgemv (CblasNoTrans, 1.0, A, x, 0.0, r);
		gsl_vector_sub (r, b);

		// 
		// if the norm of the residual is not 0 (zero), then the system
		// does not have a unique solution
		//
		double residual_norm = gsl_blas_dnrm2 (r);

		//
		// save the solution in the parameters of the propagation model
		//
		fprintf (stdout,
				 "*** INFO: found optimal values for E/// (residual is %g)\n", 
				 residual_norm);
		for (ix = 0; ix < Dim; ix ++)
            if (! isnan (gsl_vector_get (x, ix)))
			    tx_params->eric_params[ix] = gsl_vector_get (x, ix);
	}
	else
	{
		fprintf (stderr,
				 "*** WARNING: Linear system for transmitter [%s] is not solvable. Using [%.2f, %.2f, %.2f, %.2f]\n",
				 tx_params->tx_name,
                 tx_params->eric_params[0],
                 tx_params->eric_params[1],
                 tx_params->eric_params[2],
                 tx_params->eric_params[3]);
	}

	//
	// display the parameter values for E///
	//
	for (ix = 0; ix < 4; ix ++)
		fprintf (stdout,
				 "\tA%d\t%.3f\n",
				 ix,
				 tx_params->eric_params[ix]);
	
#ifdef _DEBUG_INFO_
    //
    // DEBUG: dump the matrices
    //
    printf ("A = \n");
    for (ix = 0; ix < Dim; ix ++)
    {
        for (iy = 0; iy < Dim; iy ++)
            printf ("%.5f\t", gsl_matrix_get (A, ix, iy));
        printf ("\n");
    }

    printf ("b = \n");
    for (ix = 0; ix < Dim; ix ++)
        printf ("%.5f\t", gsl_vector_get (b, ix));

    printf ("\nx = \n");
    for (ix = 0; ix < Dim; ix ++)
        printf ("%.5f\t", gsl_vector_get (x, ix));
    printf ("\n");
#endif

    //
    // free allocated elements
    //
    gsl_permutation_free (p);
    gsl_vector_free (r);
    gsl_vector_free (x);
    gsl_vector_free (b);
    gsl_matrix_free (A);

#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
#endif
}

