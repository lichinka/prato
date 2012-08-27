#define _ANTENNA_PI_    3.14159265358979

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <grass/gis.h>
#include <grass/glocale.h>

#include "performance/metric.h"



//
// the horizontal and vertical antenna diagram are saved here
// 
static double *horizontal = NULL;
static double *vertical = NULL;

//
// the antenna gain
//
static double gain = -1.0;



/**
 * Reads the antenna gain and directional diagrams from the file
 * given. The data read is saved in the output parameters:
 *
 *      horiz_diag  Horizontal antenna diagram
 *      vert_diag   Vertical antenna diagram
 *
 * This function returns the antenna gain read.-
 */
static double read_antenna_diagram (const char *file_name, 
                                    double *horiz_diag,
                                    double *vert_diag)
{
    int j;
    double gain = -1.0;
    FILE *in;                           // input file descriptor
    char buffer[256], text[32];         // char buffers to read from the file

	//
    // try to open the antenna file for reading
    //
	if ( (in = fopen (file_name, "r")) == NULL )
	    G_fatal_error(_("Unable to open antenna diagram from file <%s>"), file_name);

	//
    // read the gain and find the beginning of horizontal diagram
    //
	double temp_gain;
	while (1)
	{
		if (!fgets (buffer, 250, in)) 
		{	
			G_fatal_error (_("Empty or corrupted antenna diagram file <%s>"), file_name); 
			break;
		}
		sscanf (buffer, "%s %lf", text, &temp_gain);
		if (strcmp (text, "GAIN") == 0)		  
		    gain = temp_gain + 2.15;  // with respect to the isotropic antenna 
		if (strcmp (text, "HORIZONTAL") == 0)
        {
            //
            // we have reached the beggining of the horizontal data
            //
            break;
        }
	}
    //
    // read horizontal data - one angle degree per step
    //
	double angle, loss;
	for (j = 0; j < 360; j++)
	{
		fgets (buffer, 250, in); 
		sscanf (buffer, "%lf %lf", &angle, &loss);
		if (j != (int)angle)
		{
			G_fatal_error (_("Bad antenna diagram format.")); 
			break; 
		}
		horiz_diag[j] = loss;
	}
	//
    // skip one line ("VERTICAL 360")
    //
	fgets (buffer, 250, in); 	

	//
    // read vertical data - one angle degree per step 
    //
	for (j = 0; j < 360; j++)
	{
		fgets (buffer, 250, in); 
		sscanf (buffer, "%lf %lf", &angle, &loss);
		if (j != (int)angle)
		{
			G_fatal_error(_("Bad antenna diagram format.")); 
			break; 
		}
		vert_diag[j] = loss;
	}
    //
    // close the file and return antenna gain read
    //
	fclose (in);
    return gain;
}



/**
 * Calculates the antenna gain on a specific point of the area.
 *
 * d_east       difference between receiver and transmitter eastern 
 *              coordinates;
 * d_north      difference between receiver and transmitter northern
 *              coordinates;
 * total_height transmitter's height above sea level;
 * beam_dir     direction of the antenna beam, in degrees;
 * mech_tilt    mechanical antenna tilt angle, in degress;
 * dem_height   height above sea level from the DEM;
 * rec_height   Rx antenna height above ground level;
 * dist_Tx_Rx   distance between receiver and transmitter;
 * path_loss    path-loss value at the current point.-
 *
 */
static FCELL antenna_influence_on_point (const double d_east,
                                         const double d_north,
                                         const double total_height,
                                         const int beam_dir,
                                         const int mech_tilt,
                                         const FCELL dem_height,
                                         const double rec_height,
                                         const double dist_Tx_Rx,
                                         const FCELL path_loss)
{
    //
    // local variables
    //
    double hor_coor_angle, hor_diag_angle;
    double horizontal_loss;
    double vert_coor_angle, vert_diag_angle;
    double vertical_loss;
    double height_diff_Tx_Rx;
    double temp_angle;

    //
    // determine horizontal angle and loss
    //
    temp_angle = atan (d_east / d_north);
    if (temp_angle < 0)
      temp_angle = - temp_angle;
             
    if (d_north >= 0 && d_east >= 0)
      hor_coor_angle = temp_angle;
    else if (d_north >= 0 && d_east < 0)
      hor_coor_angle = 2*_ANTENNA_PI_ - temp_angle;
    else if (d_north < 0 && d_east < 0)
      hor_coor_angle = _ANTENNA_PI_ + temp_angle;
    else /* (d_north < 0 && d_east >= 0) */
      hor_coor_angle = _ANTENNA_PI_ - temp_angle;

    //
    // convert from radians to degrees
    //
    hor_coor_angle = hor_coor_angle * 180 / _ANTENNA_PI_;  
     
    hor_diag_angle = hor_coor_angle - (double) beam_dir;

    if (hor_diag_angle < 0)
       hor_diag_angle = 360 + hor_diag_angle;

    /* to prevent reading unallocated data (diagram comprises values 0 - 359) */
    temp_angle = ceil (hor_diag_angle);
    if (temp_angle == 360)
        temp_angle = 0;

    /* interpolation */
#ifdef _PERFORMANCE_METRICS_
    memory_access (3, 4);
#endif

    int index = (int) floor (hor_diag_angle);
    horizontal_loss  = horizontal[index];
    horizontal_loss += ((horizontal[(int)temp_angle] - horizontal[index])*(hor_diag_angle - floor(hor_diag_angle)));
    
    /* determine vertical angle and loss */
    height_diff_Tx_Rx = total_height - (double)dem_height - rec_height;

    vert_coor_angle = atan (height_diff_Tx_Rx / (dist_Tx_Rx * 1000));
    vert_coor_angle = vert_coor_angle * 180 / _ANTENNA_PI_;	
  
    if (vert_coor_angle < 0)
        vert_coor_angle = 360 + vert_coor_angle;

    // |-->
    // 3.1.2012 - Vilhar
    // Calculate the impact of mechanical tilt with respect to horizontal angle. 
    // At 0 degrees, the value is the same as the input. At 180 degrees, the 
    // input value is negative. In between, we interpolate. This correction does
    // not contribute essentially, but nevertheless should improve the final
    // result slightly. 

    double mechanicalAntennaTilt_Corrected;

    if (hor_diag_angle >= 0 && hor_diag_angle <= 180)
        mechanicalAntennaTilt_Corrected = (double)mech_tilt * (1 - (hor_diag_angle / 90));
    else if (hor_diag_angle > 180 && hor_diag_angle <= 360)
        mechanicalAntennaTilt_Corrected = (double)mech_tilt * ((hor_diag_angle / 90) - 3);
    else
        G_fatal_error(_("Horizontal angle is not between 0 and 360 degrees.")); 

    // -->|		

    vert_diag_angle = vert_coor_angle - (double)mechanicalAntennaTilt_Corrected;
  
    if (vert_diag_angle < 0)
        vert_diag_angle = 360 + vert_diag_angle;

    /* to prevent reading unallocated data (diagram comprises values 0 - 359) */
    temp_angle = ceil(vert_diag_angle);
    if (temp_angle == 360)
        temp_angle = 0;

    /* interpolation */
#ifdef _PERFORMANCE_METRICS_
    memory_access (3, 4);
#endif
    vertical_loss  = vertical[(int)floor(vert_diag_angle)];
    vertical_loss += ((vertical[(int)temp_angle] - vertical[(int)floor(vert_diag_angle)])*(vert_diag_angle - floor(vert_diag_angle)));
  
    /* finally take into account pathloss for determined diagram angles and antenna gain */
    return (FCELL)((double)path_loss + horizontal_loss + vertical_loss - gain);
}



/**
 * Standard CPU version of the antenna influence algorithm.
 * For a description of the parameters used, see the function
 * 'calculate_antenna_influence'.-
 */
static inline void 
antenna_influence_standard (const double east,
                            const double north,
                            const double total_height,
                            const int beam_dir,
                            const int mech_tilt,
                            const double radius,
                            const double rec_height,
                            const int nrows,
                            const int ncols,
                            const double west_ext,
                            const double north_ext,
                            const double ew_res,
                            const double ns_res,
                            double **dem,
                            double **path_loss)
{
    int r, c;

    //
    // this loop measures around 281.000 MFlops
    // 
    for (r = 0; r < nrows; r++) 
    {
        FCELL f_in, f_out, f_in_dem; 

        // calculate receiver coordinates
        double rec_north = north_ext - (r * ns_res);

        // calculate differences between receiver and transmitter coordinates
        double d_north = rec_north - north;

        //
        // ... and each column 
        // 
        for (c = 0; c < ncols; c++) 
        { 
#ifdef _PERFORMANCE_METRICS_
            memory_access (2, 8);
#endif
            f_in = (FCELL) path_loss[r][c];
            f_in_dem = (FCELL) dem[r][c];

            // calculate receiver coordinates
            double rec_east = west_ext + (c * ew_res);

            // calculate differences between receiver and transmitter coordinates
            double d_east = rec_east - east;
            
            // calculate distance between Tx and Rx
            double dist_Tx_Rx = sqrt (pow (d_east, 2) + 
                                      pow (d_north, 2));
            dist_Tx_Rx = dist_Tx_Rx / 1000;
           
            // If distance between Rx and Tx exceeds given radius, continue with other cells 
            FCELL null_f_out;
            if (dist_Tx_Rx > radius)
            {
                G_set_f_null_value (&null_f_out, 1);   
                f_out = null_f_out;
            }
            else
            {    
                f_out = antenna_influence_on_point (d_east,
                                                    d_north,
                                                    total_height,
                                                    beam_dir,
                                                    mech_tilt,
                                                    f_in_dem,
                                                    rec_height,
                                                    dist_Tx_Rx,
                                                    f_in);
            }
            // 
            // save the result in the output matrix
            //
            path_loss[r][c] = f_out;
        }
    }
}



/**
 * Unrolled-loop CPU version of the antenna influence algorithm.
 * For a description of the parameters used, see the function
 * 'calculate_antenna_influence'.-
 */
static inline void 
antenna_influence_unrolled (const double east,
                            const double north,
                            const double total_height,
                            const int beam_dir,
                            const int mech_tilt,
                            const double radius,
                            const double rec_height,
                            const int nrows,
                            const int ncols,
                            const double west_ext,
                            const double north_ext,
                            const double ew_res,
                            const double ns_res,
                            double **dem,
                            double **path_loss)
{
    int r, c;

    //
    // unroll the loops to improve the number of flops per memory transfer
    //
    // this version measures ~292 MFlops and ~341 MB/s
    //
    for (r = 0; r < (nrows - 1) / 2; r ++)
    {
        // calculate receiver coordinates
        double rec_north [4];

        rec_north[0] = north_ext - ((2*r) * ns_res);
        rec_north[1] = rec_north[0];
        rec_north[2] = north_ext - ((2*r + 1) * ns_res);
        rec_north[3] = rec_north[2];

        // calculate differences between receiver and transmitter coordinates
        double d_north [4], d_north_square [4];

        d_north[0]  = (rec_north[0] - north);
        d_north[1]  = (rec_north[1] - north);
        d_north[2]  = (rec_north[2] - north);
        d_north[3]  = (rec_north[3] - north);
        d_north_square[0] = d_north[0] * d_north[0];
        d_north_square[1] = d_north[1] * d_north[1];
        d_north_square[2] = d_north[2] * d_north[2];
        d_north_square[3] = d_north[3] * d_north[3];

        for (c = 0; c < (ncols - 1) / 2; c ++)
        {
            int i;
            FCELL f_in [4];
            FCELL f_in_dem [4]; 
            FCELL f_out [4];

#ifdef _PERFORMANCE_METRICS_
            memory_access (8, 8);
#endif

            f_in[0] = (FCELL) path_loss[2*r][2*c];
            f_in[1] = (FCELL) path_loss[2*r][2*c + 1];
            f_in[2] = (FCELL) path_loss[2*r + 1][2*c];
            f_in[3] = (FCELL) path_loss[2*r + 1][2*c + 1];

            f_in_dem[0] = (FCELL) dem[2*r][2*c];
            f_in_dem[1] = (FCELL) dem[2*r][2*c + 1];
            f_in_dem[2] = (FCELL) dem[2*r + 1][2*c];
            f_in_dem[3] = (FCELL) dem[2*r + 1][2*c + 1];

            G_set_f_null_value (f_out, 4);

            // calculate receiver coordinates
            double rec_east [4];

            rec_east[0] = west_ext + ((2*c) * ew_res);
            rec_east[1] = west_ext + ((2*c + 1) * ew_res);
            rec_east[2] = rec_east[0];
            rec_east[3] = rec_east[1];

            for (i = 0; i < 4; i ++)
            {
                // calculate distance between receiver and transmitter
                double d_east = rec_east[i] - east;
                double d_east_square = d_east * d_east;
                
                // calculate distance between Tx and Rx in kilometers
                double dist_Tx_Rx = sqrt (d_east_square + d_north_square[i]);
                dist_Tx_Rx /= 1000.0;

                // process this point iif its distance to the Tx falls within
                // the given radius
                if (dist_Tx_Rx <= radius)
                {   
                    f_out[i] = antenna_influence_on_point (d_east,
                                                           d_north[i],
                                                           total_height,
                                                           beam_dir,
                                                           mech_tilt,
                                                           f_in_dem[i],
                                                           rec_height,
                                                           dist_Tx_Rx,
                                                           f_in[i]);
                }
            }
            // 
            // save the result in the output matrix
            //
            path_loss[2*r][2*c]         = f_out[0];
            path_loss[2*r][2*c + 1]     = f_out[1];
            path_loss[2*r + 1][2*c]     = f_out[2];
            path_loss[2*r + 1][2*c + 1] = f_out[3];
        }
    }
    //
    // the last row, for all the columns
    //
    r = nrows - 1;

    // calculate receiver coordinates
    double rec_north = north_ext - (r * ns_res);

    // calculate differences between receiver and transmitter coordinates
    double d_north = rec_north - north;
    double d_north_square = d_north * d_north;

    for (c = 0; c < ncols; c ++)
    {
        FCELL f_in;
        FCELL f_in_dem; 
        FCELL f_out;

#ifdef _PERFORMANCE_METRICS_
        memory_access (2, 8);
#endif

        f_in = (FCELL) path_loss[r][c];

        f_in_dem = (FCELL) dem[r][c];

        G_set_f_null_value (&f_out, 1);

        // calculate receiver coordinates
        double rec_east = west_ext + (c * ew_res);

        // calculate distance between receiver and transmitter coordinates
        double d_east = rec_east - east;
        double d_east_square = d_east * d_east;
        
        // calculate distance between Tx and Rx in kilometers
        double dist_Tx_Rx = sqrt (d_east_square + d_north_square);
        dist_Tx_Rx /= 1000.0;

        // process this point iif its distance to the Tx falls within
        // the given radius
        if (dist_Tx_Rx <= radius)
        {    
            f_out = antenna_influence_on_point (d_east,
                                                d_north,
                                                total_height,
                                                beam_dir,
                                                mech_tilt,
                                                f_in_dem,
                                                rec_height,
                                                dist_Tx_Rx,
                                                f_in);
        }
        // 
        // save the result in the output matrix
        //
        path_loss[r][c] = f_out;
    }
    //
    // the last column, for all the rows
    //
    c = ncols - 1;

    for (r = 0; r < nrows; r ++)
    { 
        // calculate receiver coordinates
        double rec_north = north_ext - (r * ns_res);

        // calculate differences between receiver and transmitter coordinates
        double d_north = rec_north - north;
        double d_north_square = d_north * d_north;

        FCELL f_in;
        FCELL f_in_dem; 
        FCELL f_out;

#ifdef _PERFORMANCE_METRICS_
        memory_access (2, 8);
#endif
        f_in = (FCELL) path_loss[r][c];

        f_in_dem = (FCELL) dem[r][c];

        G_set_f_null_value (&f_out, 1);

        // calculate receiver coordinates
        double rec_east = west_ext + (c * ew_res);

        // calculate distance between receiver and transmitter coordinates
        double d_east = rec_east - east;
        double d_east_square = d_east * d_east;
        
        // calculate distance between Tx and Rx in kilometers
        double dist_Tx_Rx = sqrt (d_east_square + d_north_square);
        dist_Tx_Rx /= 1000.0;
       
        // process this point iif its distance to the Tx falls within
        // the given radius
        if (dist_Tx_Rx <= radius)
        {    
            f_out = antenna_influence_on_point (d_east,
                                                d_north,
                                                total_height,
                                                beam_dir,
                                                mech_tilt,
                                                f_in_dem,
                                                rec_height,
                                                dist_Tx_Rx,
                                                f_in);
        }
        // 
        // save the result in the output matrix
        //
        path_loss[r][c] = f_out;
    }
}



/**
 * Calculates additional gain/pathloss according to the antenna's
 * 3-dimensional diagram.
 *
 * ant_dir      Directory containing antenna diagram file;
 * ant_file     file name of the antenna diagram;
 * east         transmitter's eastern coordinate;
 * north        transmitter's northern coordinate;
 * total_height transmitter's height above sea level;
 * beam_dir     direction of the antenna beam, in degrees;
 * mech_tilt    mechanical antenna tilt angle, in degress;
 * radius       calculation radius around the given transmitter;
 * rec_height   Rx antenna height above ground level;
 * nrows        number of rows within the 2D-matrices;
 * ncols        number of columns within the 2D-matrices;
 * west_ext     western raster extent;
 * north_ext    northern raster extent;
 * ew_res       east/west raster resolution;
 * ns_res       north/south raster resolution;
 * dem          a 2D-matrix containing the digital elevation model 
 *              data;
 * path_loss    a 2D-matrix containing the path-loss values for the
 *              isotrophic antenna
 *              WARNING: the output of this function overwrites this
 *              2D-matrix!
 *
 */
void calculate_antenna_influence (const char *ant_dir,
                                  const char *ant_file,
                                  const double east,
                                  const double north,
                                  const double total_height,
                                  const int beam_dir,
                                  const int mech_tilt,
                                  const double radius,
                                  const double rec_height,
                                  const int nrows,
                                  const int ncols,
                                  const double west_ext,
                                  const double north_ext,
                                  const double ew_res,
                                  const double ns_res,
                                  double **dem,
                                  double **path_loss)
{
    //
    // do we have to load the antenna diagram?
    //
    if ((horizontal == NULL) && (vertical == NULL))
    {
        // allocate memory for the horizontal and vertical diagrams
        horizontal = (double *) calloc (360, sizeof (double));
        vertical = (double *) calloc (360, sizeof (double));

        // get antenna's gain and directional diagrams
        char fileName [1000];
        strcpy (fileName, ant_dir);
        strcat (fileName, "/");
        strcat (fileName, ant_file);
        gain = read_antenna_diagram (fileName,
                                     horizontal,
                                     vertical);
    }
    //antenna_influence_unrolled (east,
    antenna_influence_standard (east,
                                north,
                                total_height,
                                beam_dir,
                                mech_tilt,
                                radius,
                                rec_height,
                                nrows,
                                ncols,
                                west_ext,
                                north_ext,
                                ew_res,
                                ns_res,
                                dem,
                                path_loss);
}

