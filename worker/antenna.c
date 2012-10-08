#include <assert.h>

#include "worker/antenna.h"
#include "worker/common_ocl.h"
#include "performance/metric.h"




/**
 * Reads the antenna gain and directional diagrams from the file
 * given. The data read is saved in the output parameters:
 *
 *      horiz_diag  Horizontal antenna diagram
 *      vert_diag   Vertical antenna diagram
 *
 * This function returns the antenna gain read.-
 */
static double 
read_antenna_diagram (const char *file_name, 
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
    {
        fprintf (stderr, "Unable to open antenna diagram from file <%s>", file_name);
        exit (1);
    }

	//
    // read the gain and find the beginning of horizontal diagram
    //
	double temp_gain;
	while (1)
	{
		if (!fgets (buffer, 250, in)) 
		{	
			fprintf (stderr, "Empty or corrupted antenna diagram file <%s>", file_name); 
			exit (1);
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
	for (j = 0; j < _DIAGRAM_SIZE_; j++)
	{
		fgets (buffer, 250, in); 
		sscanf (buffer, "%lf %lf", &angle, &loss);
		if (j != (int)angle)
		{
			fprintf (stderr, "Bad antenna diagram format."); 
			exit (1);
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
	for (j = 0; j < _DIAGRAM_SIZE_; j++)
	{
		fgets (buffer, 250, in); 
		sscanf (buffer, "%lf %lf", &angle, &loss);
		if (j != (int)angle)
		{
			fprintf (stderr, "Bad antenna diagram format."); 
			exit (1);
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
 * path_loss    path-loss value at the current point;
 * diagram      the antenna diagram and gain.-
 *
 */
static float 
antenna_influence_on_point (const double d_east,
                            const double d_north,
                            const double total_height,
                            const int beam_dir,
                            const int mech_tilt,
                            const float dem_height,
                            const double rec_height,
                            const double dist_Tx_Rx,
                            const float path_loss,
                            const Diagram *diagram)
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
      hor_coor_angle = 2*_PI_ - temp_angle;
    else if (d_north < 0 && d_east < 0)
      hor_coor_angle = _PI_ + temp_angle;
    else /* (d_north < 0 && d_east >= 0) */
      hor_coor_angle = _PI_ - temp_angle;

    //
    // convert from radians to degrees
    //
    hor_coor_angle = hor_coor_angle * 180 / _PI_;  
     
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
    horizontal_loss  = diagram->horizontal[index];
    horizontal_loss += ((diagram->horizontal[(int)temp_angle] - diagram->horizontal[index])*(hor_diag_angle - floor(hor_diag_angle)));
    
    /* determine vertical angle and loss */
    height_diff_Tx_Rx = total_height - (double)dem_height - rec_height;

    vert_coor_angle = atan (height_diff_Tx_Rx / (dist_Tx_Rx * 1000));
    vert_coor_angle = vert_coor_angle * 180 / _PI_;	
  
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
    {
        fprintf (stderr, "Horizontal angle is not between 0 and 360 degrees."); 
	    exit (1);
    }

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
    vertical_loss  = diagram->vertical[(int)floor(vert_diag_angle)];
    vertical_loss += ((diagram->vertical[(int)temp_angle] - diagram->vertical[(int)floor(vert_diag_angle)])*(vert_diag_angle - floor(vert_diag_angle)));
  
    /* finally take into account pathloss for determined diagram angles and antenna gain */
    return (float)((double)path_loss + horizontal_loss + vertical_loss - diagram->gain);
}



/**
 * Standard CPU version of the antenna influence algorithm.
 * For a description of the parameters used, see the function
 * 'calculate_antenna_influence'.
 *
 */
static void 
antenna_influence_cpu (const double east,
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
                       const float  null_value,
                       double **dem,
                       double **path_loss)
{
    int r, c;

    //
    // this loop measures around 281.000 MFlops
    // 
    for (r = 0; r < nrows; r++) 
    {
        float f_in, f_out, f_in_dem; 

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
            f_in = (float) path_loss[r][c];
            f_in_dem = (float) dem[r][c];

            // calculate receiver coordinates
            double rec_east = west_ext + (c * ew_res);

            // calculate differences between receiver and transmitter coordinates
            double d_east = rec_east - east;
            
            // calculate distance between Tx and Rx
            double dist_Tx_Rx = sqrt (pow (d_east, 2) + 
                                      pow (d_north, 2));
            dist_Tx_Rx = dist_Tx_Rx / 1000;
           
            // If distance between Rx and Tx exceeds given radius, continue with other cells 
            if (dist_Tx_Rx > radius)
                f_out = null_value;
            else
                f_out = antenna_influence_on_point (d_east,
                                                    d_north,
                                                    total_height,
                                                    beam_dir,
                                                    mech_tilt,
                                                    f_in_dem,
                                                    rec_height,
                                                    dist_Tx_Rx,
                                                    f_in,
                                                    diagram);
            // 
            // save the result in the output matrix
            //
            path_loss[r][c] = (double) f_out;
        }
    }
}



/**
 * Unrolled-loop CPU version of the antenna influence algorithm.
 * For a description of the parameters used, see the function
 * 'calculate_antenna_influence'.
 *
 * null_value   The value representing NULL on the output map.-
 *
 */
static void 
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
                            const float null_value,
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
            float f_in [4];
            float f_in_dem [4]; 
            float f_out [4];

#ifdef _PERFORMANCE_METRICS_
            memory_access (8, 8);
#endif

            f_in[0] = (float) path_loss[2*r][2*c];
            f_in[1] = (float) path_loss[2*r][2*c + 1];
            f_in[2] = (float) path_loss[2*r + 1][2*c];
            f_in[3] = (float) path_loss[2*r + 1][2*c + 1];

            f_in_dem[0] = (float) dem[2*r][2*c];
            f_in_dem[1] = (float) dem[2*r][2*c + 1];
            f_in_dem[2] = (float) dem[2*r + 1][2*c];
            f_in_dem[3] = (float) dem[2*r + 1][2*c + 1];

            f_out[0] = null_value;
            f_out[1] = null_value;
            f_out[2] = null_value;
            f_out[3] = null_value;

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
                                                           f_in[i],
                                                           diagram);
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
        float f_in;
        float f_in_dem; 
        float f_out;

#ifdef _PERFORMANCE_METRICS_
        memory_access (2, 8);
#endif

        f_in = (float) path_loss[r][c];

        f_in_dem = (float) dem[r][c];

        f_out = null_value;

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
                                                f_in,
                                                diagram);
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

        float f_in;
        float f_in_dem; 
        float f_out;

#ifdef _PERFORMANCE_METRICS_
        memory_access (2, 8);
#endif
        f_in = (float) path_loss[r][c];

        f_in_dem = (float) dem[r][c];

        f_out = null_value;

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
                                                f_in,
                                                diagram);
        }
        // 
        // save the result in the output matrix
        //
        path_loss[r][c] = f_out;
    }
}



/**
 * GPU version of the antenna influence algorithm.
 * For a description of the parameters used, see the function
 * 'calculate_antenna_influence'.
 *
 */
static void 
antenna_influence_gpu (const double east,
                       const double north,
                       const double total_height,
                       const int beam_dir,
                       const int mech_tilt,
                       const double frequency,
                       const double radius,
                       const double rec_height,
                       const int nrows,
                       const int ncols,
                       const double west_ext,
                       const double north_ext,
                       const double ew_res,
                       const double ns_res,
                       const float  null_value,
                       double **dem,
                       double **path_loss)
{
    cl_int            status;
    cl_context        context;
    cl_command_queue *list_queues;
    uint              num_queues = 1;
    cl_event          events [num_queues];
    cl_kernel         kernel;

    // the horizontal and vertical resolution of one raster pixel should match
    assert (ew_res == ns_res);

    // create a new OpenCL environment
    set_opencl_env_multiple_queues (num_queues, 
                                    &list_queues, 
                                    &context);
    // build and activate the kernel
    build_kernel_from_file (&context,
                            "sector_kern",
                            "r.coverage.cl",
                            &kernel);

    // number of processing tiles needed around each transmitter 
    if ((radius*2) < _WORK_ITEMS_PER_DIMENSION_)
    {
        fprintf (stderr, "Cannot execute. Increase radius and try again.");
        exit (1);
    }

    size_t ntile = (radius*2) / _WORK_ITEMS_PER_DIMENSION_ + 1;

    // defines a 2D execution range for the kernel
    size_t global_offsets [] = {0, 0};
    size_t global_sizes [] = {ntile * _WORK_ITEMS_PER_DIMENSION_,
                              ntile * _WORK_ITEMS_PER_DIMENSION_};
    size_t local_sizes [] = {_WORK_ITEMS_PER_DIMENSION_,
                             _WORK_ITEMS_PER_DIMENSION_};
    status = clEnqueueNDRangeKernel (list_queues[0],
                                     kernel,
                                     2,
                                     global_offsets,
                                     global_sizes,
                                     local_sizes,
                                     0,
                                     NULL,
                                     &(events[0]));
    check_error (status, "Set kernel range");

    // kernel parameters 
    set_kernel_arg (&kernel, 
                    0,
                    sizeof (ns_res),
                    &ns_res);
    set_kernel_arg (&kernel, 
                    1,
                    sizeof (north_ext),
                    &north_ext);
    set_kernel_arg (&kernel, 
                    2,
                    sizeof (west_ext),
                    &west_ext);
    set_kernel_arg (&kernel, 
                    3,
                    sizeof (int),
                    &ncols);
    set_kernel_arg (&kernel, 
                    6,
                    sizeof (rec_height),
                    &rec_height);
    set_kernel_arg (&kernel, 
                    7,
                    sizeof (double),
                    &frequency);
    set_kernel_arg (&kernel, 
                    8,
                    sizeof (double),
                    &(diagram->gain));
    set_kernel_arg (&kernel, 
                    9,
                    sizeof (int),
                    &beam_dir);
    set_kernel_arg (&kernel, 
                    10,
                    sizeof (int),
                    &mech_tilt);

    // reserve local memory on the device
    size_t lmem_size = _WORK_ITEMS_PER_DIMENSION_ *
                       _WORK_ITEMS_PER_DIMENSION_ *
                       sizeof (cl_float2);
    set_local_mem (&kernel, 
                   15,
                   lmem_size);

    // create OpenCL buffers
    cl_mem horiz_diag_in_dev;
    cl_mem vert_diag_in_dev;
    cl_mem dem_in_dev;
    cl_mem sect_out_dev;
    size_t diagram_size = 360;
    size_t dem_in_size = nrows * 
                         ncols * 
                         sizeof (dem[0][0]);
    size_t sect_out_size = nrows * 
                           ncols * 
                           sizeof (path_loss[0][0]);

    create_buffer (&context,
                   CL_MEM_READ_ONLY, 
                   diagram_size * sizeof (double),
                   diagram->horizontal,
                   &horiz_diag_in_dev);
    create_buffer (&context,
                   CL_MEM_READ_ONLY, 
                   diagram_size * sizeof (double),
                   diagram->vertical,
                   &vert_diag_in_dev);
    create_buffer (&context,
                   CL_MEM_READ_ONLY, 
                   dem_in_size,
                   dem,
                   &dem_in_dev);
    create_buffer (&context,
                   CL_MEM_WRITE_ONLY, 
                   sect_out_size,
                   path_loss,
                   &sect_out_dev);

    // send input data to the device
    write_buffer (&(list_queues[0]),
                  &horiz_diag_in_dev,
                  CL_TRUE,
                  diagram_size * sizeof (double),
                  diagram->horizontal);
    write_buffer (&(list_queues[0]),
                  &vert_diag_in_dev,
                  CL_TRUE,
                  diagram_size * sizeof (double),
                  diagram->vertical);
    write_buffer (&(list_queues[0]),
                  &dem_in_dev,
                  CL_TRUE,
                  dem_in_size,
                  dem);

    // set the remaining parameters
    set_kernel_arg (&kernel, 
                    11,
                    sizeof (double),
                    &horiz_diag_in_dev);
    set_kernel_arg (&kernel, 
                    12,
                    sizeof (double),
                    &vert_diag_in_dev);
    set_kernel_arg (&kernel, 
                    13,
                    sizeof (double),
                    &dem_in_dev);
    set_kernel_arg (&kernel, 
                    14,
                    sizeof (double),
                    &sect_out_dev);

    printf ("Simulating antenna influence for one sector ...");

    // tile offset within the calculation radius, given in pixel coordinates
    cl_int2 tx_offset;
    tx_offset.s[0] = (int) ((east - west_ext) / ew_res - radius);
    tx_offset.s[1] = (int) ((north_ext - north) / ns_res - radius);
    // transmitter data
    cl_double4 tx_data;
    tx_data.s[0] = (double) east;                 // transmitter coordinate
    tx_data.s[1] = (double) north;                // transmitter coordinate
    tx_data.s[2] = (double) total_height;         // height above sea level
    tx_data.s[3] = (double) 75;                   // FIXME: antenna height

    printf ("Transmitter located at %.0f,%.0f,%.0f,%.0f\n", tx_data.s[0],
                                                            tx_data.s[1],
                                                            tx_data.s[2],
                                                            tx_data.s[3]);
    // set kernel parameters, specific for this transmitter
    set_kernel_arg (&kernel,
                    4,
                    sizeof (double),
                    &tx_data);
    set_kernel_arg (&kernel,
                    5,
                    sizeof (double),
                    &tx_offset);

    // start kernel execution 
    run_kernel (&(list_queues[0]),
                &kernel);

    // sync memory
    read_buffer (&(list_queues[0]),
                 &sect_out_dev,
                 CL_TRUE,
                 sect_out_size,
                 path_loss);

    // deactivate the OpenCL environment
    release_opencl (num_queues, 
                    &list_queues, 
                    &context);
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
 * frequency    transmitter frequency, in Mhz;
 * radius       calculation radius around the given transmitter;
 * rec_height   Rx antenna height above ground level;
 * nrows        number of rows within the 2D-matrices;
 * ncols        number of columns within the 2D-matrices;
 * west_ext     western raster extent;
 * north_ext    northern raster extent;
 * ew_res       east/west raster resolution;
 * ns_res       north/south raster resolution;
 * null_value   the value representing NULL on the output map;
 * dem          a 2D-matrix containing the digital elevation model data;
 * path_loss    a 2D-matrix containing the path-loss values for the isotrophic
 *              antenna;
 *              WARNING: the output of this function overwrites this 2D-matrix!
 *
 */
void calculate_antenna_influence (const char *ant_dir,
                                  const char *ant_file,
                                  const double east,
                                  const double north,
                                  const double total_height,
                                  const int beam_dir,
                                  const int mech_tilt,
                                  const double frequency,
                                  const double radius,
                                  const double rec_height,
                                  const int nrows,
                                  const int ncols,
                                  const double west_ext,
                                  const double north_ext,
                                  const double ew_res,
                                  const double ns_res,
                                  const float null_value,
                                  double **dem,
                                  double **path_loss)
{
    //
    // do we have to load the antenna diagram?
    //
    if (diagram == NULL)
    {
        // allocate memory for the horizontal and vertical diagrams
        diagram = (Diagram *) malloc (sizeof (Diagram));
        diagram->horizontal = (double *) calloc (_DIAGRAM_SIZE_, 
                                                 sizeof (double));
        diagram->vertical = (double *) calloc (_DIAGRAM_SIZE_, 
                                               sizeof (double));

        // get antenna's gain and directional diagrams
        char fileName [1000];
        if (strlen (ant_dir) == 0)
        {
            fprintf (stderr, "ERROR Directory containing antenna files is empty\n");
            exit (1);
        }
        strcpy (fileName, ant_dir);
        strcat (fileName, "/");
        if (strlen (ant_file) == 0)
        {
            fprintf (stderr, "ERROR File name of antenna diagram is empty\n");
            exit (1);
        }
        strcat (fileName, ant_file);

        diagram->gain = read_antenna_diagram (fileName,
                                              diagram->horizontal,
                                              diagram->vertical);
    }

    //antenna_influence_cpu (east,
    antenna_influence_gpu (east,
                           north,
                           total_height,
                           beam_dir,
                           mech_tilt,
                           frequency,
                           radius,
                           rec_height,
                           nrows,
                           ncols,
                           west_ext,
                           north_ext,
                           ew_res,
                           ns_res,
                           null_value,
                           dem,
                           path_loss);
}

