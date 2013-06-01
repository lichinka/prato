#include "worker/coverage.h"
#include "worker/antenna.h"
#include "worker/eric.h"



/**
 * Calculates the terrain profile for line-of-sight.
 *
 * Obst_high    output 2D matrix of heights;
 * Obst_dist    output 2D matrix of distances;
 * Offset       output 2D matrix of offsets.-
 *
 */
static void 
calc_profile (double** Obst_high, 
              double** Obst_dist, 
              double** Raster, 
              double** Offset, 
              const double dx, 
              const double dy, 
              const double xBS, 
              const double yBS,
              const double ZoTransBS, 
              const int xN, 
              const int yN, 
              const double scale, 
              const double radius)
{
	
	int X,Y;
	double x_tmp, y_tmp;
	double temp_dist, max_dist_old, max_dist_los, max_dist_nlos, max_dist, cand_dist, cand_el_dist;
	double offset_dist;
	
	double max_high_old, max_high_los, max_high_nlos, max_high, cand_high, cand_el_high;
	double cand_obstacle, cand_el_obstacle, max_obstacle;
	double EL, EL_temp, EL_max;
	double EL_nlos, EL_nlos_temp, EL_nlos_max;
	
	double interp_high, koef_x, koef_y; 
	
	int LOS_start = 1;
	
	
	//Patrik x_tmp = xBS + dx;
	//Patrik y_tmp = yBS + dy;
	x_tmp = floor(xBS) + 0.5 + dx;
	y_tmp = floor(yBS) + 0.5 + dy;
	//Patrik x_tmp = floor(xBS) + dx;
	//Patrik y_tmp = floor(yBS) + dy;
		
	//Patrik X = Y = 0;
		
	max_high = max_high_old = Raster [(int)xBS][(int)yBS];
	max_dist = max_dist_old = 1;		// 1 bin
	//Patrik cand_high = -999;
	//Patrik cand_dist = 1;
	EL_max = 0;
	EL_nlos_max = 0;
	while ((int)x_tmp >= 0 && (int)x_tmp < xN && (int)y_tmp >= 0 && (int)y_tmp < yN)
	{
		X = (int)x_tmp;
		Y = (int)y_tmp;
			
		//Patrik temp_dist = sqrt(pow(floor(xBS) + 0.5 - x_tmp,2)+pow(floor(yBS) + 0.5 - y_tmp,2));
		temp_dist = sqrt(pow(floor(xBS) - floor(x_tmp),2)+pow(floor(yBS) - floor(y_tmp),2));	//zaokorzena razdalja na bin - pravilno za racun
		offset_dist = sqrt(pow(floor(xBS) + 0.5 - x_tmp,2)+pow(floor(yBS) + 0.5 - y_tmp,2));	//natancna razdalja potrebna za izracun offseta pri dolocanju bina
	
        /*
         * Do not process the whole area, take the radius restriction into account
         */
		if (temp_dist*scale > radius*1000)
        {
			break;
		}
			
// interpolacija DEM-a za natancnejso visino
			if ((x_tmp - (floor(x_tmp) + 0.5)) <= 0){
				if ((X-1) < 0){
					koef_x = 0;	
				}
				else {
					koef_x = (Raster[X][Y] - Raster[X-1][Y]);
				}
			}
			else{
				if ((X+1) >= xN){
					koef_x = 0;	
				}
				else {
					koef_x = (Raster[X+1][Y] - Raster[X][Y]);
				}
			}
			
			if ((y_tmp - (floor(y_tmp) + 0.5)) <= 0){
				if ((Y-1) < 0){
					koef_y = 0;
				}
				else {
					koef_y = (Raster[X][Y] - Raster[X][Y-1]);
				}
			}
			else{
				if ((Y+1) >= yN){
					koef_y = 0;
				}
				else{
					koef_y = (Raster[X][Y+1] - Raster[X][Y]);
				}
			}
		
			
			interp_high = Raster[X][Y] + koef_x * (x_tmp - (floor(x_tmp) + 0.5)) + koef_y * (y_tmp - (floor(y_tmp) + 0.5));
// konec interpolacija DEM-a za natancnejso visino
			
			//Patrik EL = atan(temp_dist * scale / (ZoTransBS - Raster[X][Y]));
			EL = atan(temp_dist * scale / (ZoTransBS - interp_high));
			
			if (EL < 0){
				EL_temp = _PI_ + EL;
			}
			else{
				EL_temp = EL;
			}



// LOS case			
			if (EL_temp > EL_max){
				if (LOS_start == 1){
					max_high_los = max_high;
					max_dist_los = max_dist;
				}
				else{
					max_high_los = max_high_old;
					max_dist_los = max_dist_old;
				}	
	
				
				if (fabs(offset_dist - temp_dist) < Offset[X][Y]){
					Obst_high[X][Y] = max_high_los - (ZoTransBS - max_dist_los * scale / tan(EL));
					Obst_dist[X][Y] = max_dist_los;			//do sedaj prevladuloca ovira
					Offset[X][Y] = fabs(offset_dist - temp_dist);
				}
				EL_max = EL_temp;
					
				//Patrik if (Raster[X][Y] >= max_high){
					//Patrik max_high = cand_high = Raster[X][Y];
					max_high = cand_high = interp_high;
					max_dist = cand_dist = temp_dist;			//nova ovira

				//Patrik }
				
				
				cand_el_high = max_high;			// ponastavi za naslednjo NLOS iteracijo
				cand_el_dist = max_dist;
				EL_nlos_max = 0; // prestudiraj!!
			}
// NLOS case
			else{									// EL_temp <= EL_max
				LOS_start = 0;
				
				
				//Patrik EL_nlos = atan((temp_dist - cand_dist) * scale / (cand_high + 50 - Raster[X][Y]));
				EL_nlos = atan((temp_dist - cand_dist) * scale / (cand_high + 50 - interp_high));
				
				if (EL_nlos < 0){
					EL_nlos_temp = _PI_ + EL_nlos;
				}
				else{
					EL_nlos_temp = EL_nlos;
				}
				
				if (EL_nlos_temp > EL_nlos_max){
					EL_nlos_max = EL_nlos_temp;
				}
				else{					//EL_nlos_temp <= EL_nlos_max
					//Patrik cand_high = Raster[X][Y];
					cand_high = interp_high;
					cand_dist = temp_dist;
					EL_nlos_max = 0; 			// reset kota, ker smo dosegli novo tocko LOS v NLOS obmocju
				}
				
				cand_obstacle = cand_high - ZoTransBS + cand_dist * scale / tan(EL);		// kandidat za novo oviro NLOS
				cand_el_obstacle = cand_el_high - ZoTransBS + cand_el_dist * scale / tan(EL);		// izbran kandidat za novo oviro NLOS
				max_obstacle = max_high - ZoTransBS + max_dist * scale / tan(EL);					// prevladujoca ovira iz LOS
				
				if((cand_obstacle > max_obstacle) && (cand_obstacle > cand_el_obstacle)){
					cand_el_high = cand_high;
					cand_el_dist = cand_dist;		
				}
					
				if (cand_el_obstacle >= max_obstacle){
					max_high_nlos = cand_el_high;
					max_dist_nlos = cand_el_dist;				
				}
				else{
					max_high_nlos = max_high;
					max_dist_nlos = max_dist;					
				}	
					

				if (fabs(offset_dist - temp_dist) < Offset[X][Y]){
					Obst_high[X][Y] = max_high_nlos - (ZoTransBS - max_dist_nlos * scale / tan(EL));
					Obst_dist[X][Y] = max_dist_nlos;
					Offset[X][Y] = fabs(offset_dist - temp_dist);
				}
					
				if (max_high > max_high_nlos){
					max_high_old = max_high;
					max_dist_old = max_dist;
				}
				else{
					max_high_old = max_high_nlos;
					max_dist_old = max_dist_nlos;
				}
			}

		x_tmp = x_tmp + dx;
		y_tmp = y_tmp + dy;	
	}
}



/**
 * Calculates line-of-sight vectors all around the transmitter.-
 *
 */
static int 
DoProfile (double **Obst_high,
           double **Obst_dist,
           double **Offset,
           double ResDist, 
           double **Raster, 
           double xBS, 
           double yBS, 
           double ZoTransBS, 
           int xN, 
           int yN, 
           double scale, 
           double radius)
{
#ifdef _PERFORMANCE_METRICS_
    measure_time ("Line-of-sight");
#endif
	double AZI;
	int ix, iy;
	double dx, dy;
	
    //
    // LOS and obstacle height calculation is executed only once, 
    // because its results are constant throughout the optimization
    //

    // Offset ini
	for (ix = 0; ix < xN; ix++)
    {
		for (iy = 0; iy < yN; iy++)
        {
			Offset[ix][iy]=999;
		}
	}
	
	// Kvadrant I
	for (ix = 0; ix < xN; ix++)
	{
		//Patrik AZI = atan((ix - xBS) / yBS);
		AZI = atan((ix - floor(xBS)) / floor(yBS));
		
		if (cos(AZI) > sin(AZI))
        {
			//Patrik dx = sin(AZI) / cos(AZI);
			//Patrik dy = -cos(AZI) / cos(AZI);
			dx = (ix - floor(xBS)) / floor(yBS);			// tan(AZI)
			dy = -1;
		}
		else
        {
			//Patrik dx = sin(AZI) / sin(AZI);
			//Patrik dy = -cos(AZI) / sin(AZI);
			dx = 1;
			dy = -floor(yBS)/(ix - floor(xBS));				// ctan(AZI)
		}
		calc_profile (Obst_high, Obst_dist, Raster, Offset, dx, dy, xBS, yBS, ZoTransBS, xN, yN, scale, radius);
#ifdef _DEBUG_INFO_
        printf ("DoProfile -> 1st quadrant: %d\n", ix);
#endif
	}

	// Kvadrant III
	for (ix = 0; ix < xN; ix++)
	{
		//Patrik AZI = atan((ix - xBS) / (yN - yBS));
		AZI = atan((ix - floor(xBS)) / (yN - floor(yBS)));
		
		if (cos(AZI) > sin(AZI))
        {
			//Patrik dx = sin(AZI) / cos(AZI);
			//Patrik dy = cos(AZI) / cos(AZI);
			dx = (ix - floor(xBS)) / (yN - floor(yBS));			// tan(AZI)
			dy = 1;
		}
		else
        {
			//Patrik dx = sin(AZI) / sin(AZI);
			//Patrik dy = cos(AZI) / sin(AZI);
			dx = 1;
			dy = (yN - floor(yBS)) / (ix - floor(xBS));				// ctan(AZI)
		}
				
		calc_profile (Obst_high, Obst_dist, Raster, Offset, dx, dy, xBS, yBS, ZoTransBS, xN, yN, scale, radius);
#ifdef _DEBUG_INFO_
        printf ("DoProfile -> 3rd quadrant: %d\n", ix);
#endif
	} 
	
	// Kvadrant II
	for (iy = 0; iy < yN; iy++)
	{
		//Patrik AZI = atan((iy - yBS) / (xN - xBS));
		AZI = atan((iy - floor(yBS)) / (xN - floor(xBS)));
			
		if (cos(AZI) > sin(AZI))
        {
			//Patrik dx = cos(AZI) / cos(AZI);
			//Patrik dy = sin(AZI) / cos(AZI);
			dx = 1;			
			dy = (iy - floor(yBS)) / (xN - floor(xBS));		// tan(AZI)
		}
		else
        {
			//Patrik dx = cos(AZI) / sin(AZI);
			//Patrik dy = sin(AZI) / sin(AZI);
			dx = (xN - floor(xBS)) / (iy - floor(yBS));				// ctan(AZI)
			dy = 1;
		}
		
		calc_profile (Obst_high, Obst_dist, Raster, Offset, dx, dy, xBS, yBS, ZoTransBS, xN, yN, scale, radius);	
#ifdef _DEBUG_INFO_
        printf ("DoProfile -> 2nd quadrant: %d\n", ix);
#endif
	}

	// Kvadrant IV
	for (iy = 0; iy < yN; iy++)
	{
		//Patrik AZI = atan((iy - yBS) / xBS);
		AZI = atan((iy - floor(yBS)) / floor(xBS));
		
		if (cos(AZI) > sin(AZI))
        {
			//Patrik dx = -cos(AZI) / cos(AZI);
			//Patrik dy = sin(AZI) / cos(AZI);
			dx = -1;			
			dy = (iy - floor(yBS)) / floor(xBS);		// tan(AZI)
		}
		else
        {
			//Patrik dx = -cos(AZI) / sin(AZI);
			//Patrik dy = sin(AZI) / sin(AZI);
			dx = -floor(xBS) / (iy - floor(yBS));				// ctan(AZI)
			dy = 1;
		}
		
		calc_profile (Obst_high, Obst_dist, Raster, Offset, dx, dy, xBS, yBS, ZoTransBS, xN, yN, scale, radius);			
#ifdef _DEBUG_INFO_
        printf ("DoProfile -> 4th quadrant: %d\n", ix);
#endif
	}

#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
#endif
	return 0;
}



/**
 * Calculates the coverage prediction for one transmitter, using the E/// model.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters.-
 *
 */
void 
coverage (Parameters    *params,
          Tx_parameters *tx_params)
{
    //
    // calculate the terrain profile from the top of the transmitter,
    // i.e. line-of-sight, only once per transmitter
    // 
    DoProfile (tx_params->m_obst_height,
               tx_params->m_obst_dist,
               tx_params->m_obst_offset,
               1.0,
               tx_params->m_dem,
               tx_params->tx_north_coord_idx,
               tx_params->tx_east_coord_idx,
               tx_params->total_tx_height,
               tx_params->nrows,
               tx_params->ncols,
               params->map_ew_res,
               params->radius);
    //
    // execute the path-loss calculation on CPU or GPU?
    //
    if (params->use_gpu)
    {
        //
        // initialize the OpenCL environment
        //
        init_gpu (params,
                  tx_params);
#ifdef _PERFORMANCE_METRICS_
        measure_time ("E/// on GPU");
#endif
        eric_pathloss_on_gpu (params,
                              tx_params);
    }
    else
    {
#ifdef _PERFORMANCE_METRICS_
        measure_time ("E/// on CPU");
#endif
        eric_pathloss_on_cpu (params,
                              tx_params);
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

}



/**
 * Displays the calculation result in the standard output.
 *
 * params_ptr   a pointer to the structure holding configuration parameters 
 *              which are common to all transmitters;
 */
void *
output_to_stdout (void *params_ptr)
{
    int r, c;

    Parameters *params       = (Parameters *) params_ptr;
    Tx_parameters *tx_params = params->tx_params;
    
    // 
    // prepare the DB server before sending the data
    //
    fprintf (stdout, 
             "CREATE TABLE IF NOT EXISTS pathloss_%s (east float, north float, pl float);\n",
             tx_params->tx_name);
    fprintf (stdout, 
             "TRUNCATE TABLE pathloss_%s;\n",
             tx_params->tx_name);
    fprintf (stdout,
             "\\COPY pathloss_%s (east, north, pl) FROM STDIN WITH DELIMITER '|'\n",
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

            float pl = (float) tx_params->m_loss[r][c];

            if ((!isnan (pl)) && (pl != params->fcell_null_value))
                fprintf (stdout, "%.2f|%.2f|%.5f\n", east_coord,
                                                     north_coord,
                                                     pl);
        }
    }
    //
    // mark end-of transmitter data
    //
    fprintf (stdout, "\\.\n");

    return NULL;
}
