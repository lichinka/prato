#include "worker/antenna.h"




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
 * d_east           difference between receiver and transmitter eastern 
 *                  coordinates;
 * d_north          difference between receiver and transmitter northern
 *                  coordinates;
 * total_height     transmitter's height above sea level;
 * beam_dir         direction of the antenna beam, in degrees;
 * mech_tilt        mechanical antenna tilt angle, in degress;
 * dem_height       height above sea level from the DEM;
 * rec_height       Rx antenna height above ground level;
 * dist_Tx_Rx       distance between receiver and transmitter;
 * diagram          the antenna diagram and gain;
 * main_zone_horiz  indicates the horizontal loss that defines the main
 *                  antenna beam;
 * main_zone_vert   indicates the vertical loss that defines the main
 *                  antenna beam;
 * sec_zone_horiz   indicates the horizontal loss that defines the secondary
 *                  antenna beam;
 * sec_zone_vert    indicates the vertical loss that defines the secondary
 *                  antenna beam;
 * radio_zone       indicates the radio zone to which this point belongs
 *                  (output parameter);
 * antenna_loss     indicates the loss, introduced by the antenna, on this
 *                  point (output parameter);
 *
 */
static void
antenna_influence_on_point (double d_east,
                            double d_north,
                            const double total_height,
                            const int beam_dir,
                            const int mech_tilt,
                            const float dem_height,
                            const double rec_height,
                            double dist_Tx_Rx,
                            const Diagram *diagram,
                            const int main_zone_horiz,
                            const int main_zone_vert,
                            const int sec_zone_horiz,
                            const int sec_zone_vert,
                            char   *radio_zone,
                            double *antenna_loss)
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
    // the arctan cannot be calculated if any of the involved numbers is 0
    //
    if (d_east == 0)
        d_east = 0.01;
    if (d_north == 0)
        d_north = 0.01;
    temp_angle = atan (d_east / d_north);
    if (temp_angle < 0)
      temp_angle = -temp_angle;
             
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
    
    //
    // determine vertical angle and loss
    //
    height_diff_Tx_Rx = total_height - (double)dem_height - rec_height;

    //
    // the arctan cannot be calculated if any of the involved numbers is 0
    //
    if (height_diff_Tx_Rx == 0)
        height_diff_Tx_Rx = 0.001;
    if (dist_Tx_Rx == 0)
        dist_Tx_Rx = 0.001;
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
        vert_diag_angle += 360;
    if (vert_diag_angle > 360)
        vert_diag_angle -= 360;

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
 
    //
    // check if this point is within the main antenna beam, i.e. 
    // horizontally and vertically introduced loss is within threshold
    //
    if ((horizontal_loss <= main_zone_horiz) && (vertical_loss <= main_zone_vert))
        *radio_zone |= _RADIO_ZONE_MAIN_BEAM_ON_;
    //
    // check if this point is within the secondary main antenna beam, i.e. 
    // a larger main antenna beam not including the main one (only the difference)
    //
    else if ((horizontal_loss <= sec_zone_horiz) && (vertical_loss <= sec_zone_vert))
        *radio_zone |= _RADIO_ZONE_SECONDARY_BEAM_ON_;

    //
    // finally save the determined diagram losses and antenna gain 
    //
    *antenna_loss = horizontal_loss + vertical_loss - diagram->gain;
}



/**
 * Standard CPU version of the antenna influence algorithm.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 *
 */
static void 
antenna_influence_cpu (Parameters    *params,
                       Tx_parameters *tx_params)
{
    int r, c;

    //
    // this loop measures around 281.000 MFlops
    // 
    for (r = 0; r < tx_params->nrows; r++) 
    {
        float loss_out, dem_in; 

        // calculate receiver coordinates
        double rec_north = tx_params->map_north - (r * params->map_ns_res);

        // calculate differences between receiver and transmitter coordinates
        double d_north = rec_north - tx_params->tx_north_coord;

        //
        // ... and each column 
        // 
        for (c = 0; c < tx_params->ncols; c++) 
        { 
#ifdef _PERFORMANCE_METRICS_
            memory_access (2, 8);
#endif
            dem_in = (float) tx_params->m_dem[r][c];

            // calculate receiver coordinates
            double rec_east = tx_params->map_west + (c * params->map_ew_res);

            // calculate differences between receiver and transmitter coordinates
            double d_east = rec_east - tx_params->tx_east_coord;
            
            // calculate distance between Tx and Rx (in kilometers)
            double dist_Tx_Rx = sqrt (pow (d_east, 2) + 
                                      pow (d_north, 2));
            dist_Tx_Rx /= 1000;
       
            //
            // ignore pixels exceeding params->radius distance between Rx and Tx or
            // if they are too close (less than 200 mts)
            // 
            if ((dist_Tx_Rx < 0.2) || (dist_Tx_Rx > params->radius))
            {
                loss_out                       = params->null_value;
                tx_params->m_radio_zone[r][c] &= _RADIO_ZONE_MODEL_DISTANCE_OFF_;
            }
            else 
            {
                tx_params->m_radio_zone[r][c] |= _RADIO_ZONE_MODEL_DISTANCE_ON_;
                antenna_influence_on_point (d_east,
                                            d_north,
                                            tx_params->total_tx_height,
                                            tx_params->beam_direction,
                                            tx_params->mechanical_tilt,
                                            dem_in,
                                            params->rx_height_AGL,
                                            dist_Tx_Rx,
                                            tx_params->diagram,
                                            params->main_zone_horiz,
                                            params->main_zone_vert,
                                            params->sec_zone_horiz,
                                            params->sec_zone_vert,
                                            &(tx_params->m_radio_zone[r][c]),
                                            &(tx_params->m_antenna_loss[r][c]));
                loss_out = (float) (tx_params->m_loss[r][c] + 
                                    tx_params->m_antenna_loss[r][c]);
            }
            // 
            // save the result in the output matrix
            //
            tx_params->m_loss[r][c] = (double) loss_out;
        }
    }
}



/**
 * GPU version of the antenna influence algorithm.
 * The losses introduced by the antenna are kept in a separate matrix, ready
 * to be applied over the isotrophic path loss values.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 *
 */
static void 
antenna_influence_gpu (Parameters    *params,
                       Tx_parameters *tx_params)
{
    cl_mem horiz_diag_in_dev, vert_diag_in_dev;

    // activate the kernel
    activate_kernel (tx_params->ocl_obj,
                     "antenna_influence_kern");
    //
    // create the OpenCL buffers, if we are executing this function
    // for the first time
    //
    if (tx_params->m_antenna_loss_dev == NULL)
    {
        //
        // sizes of the OpenCL buffers ...
        //
        size_t diagram_size = _DIAGRAM_SIZE_ *
                               sizeof (tx_params->diagram->horizontal[0]);
        size_t ant_buff_size = tx_params->nrows * 
                               tx_params->ncols * 
                               sizeof (tx_params->m_antenna_loss[0][0]);
        size_t rad_buff_size = tx_params->nrows *
                               tx_params->ncols *
                               sizeof (tx_params->m_radio_zone[0][0]);
        //
        // ... before defining them
        //
        tx_params->m_antenna_loss_dev = (cl_mem *) malloc (sizeof (cl_mem));
        tx_params->m_radio_zone_dev   = (cl_mem *) malloc (sizeof (cl_mem));

        *(tx_params->m_antenna_loss_dev) = create_buffer (tx_params->ocl_obj,
                                                          CL_MEM_READ_WRITE, 
                                                          ant_buff_size);
        *(tx_params->m_radio_zone_dev) = create_buffer (tx_params->ocl_obj,
                                                        CL_MEM_WRITE_ONLY,
                                                        rad_buff_size);
        horiz_diag_in_dev = create_buffer (tx_params->ocl_obj,
                                           CL_MEM_READ_ONLY, 
                                           diagram_size);
        vert_diag_in_dev = create_buffer (tx_params->ocl_obj,
                                          CL_MEM_READ_ONLY, 
                                          diagram_size);
        //
        // send buffers data to the device
        //
        write_buffer (tx_params->ocl_obj,
                      0,
                      &horiz_diag_in_dev,
                      diagram_size,
                      tx_params->diagram->horizontal);
        write_buffer (tx_params->ocl_obj,
                      0,
                      &vert_diag_in_dev,
                      diagram_size,
                      tx_params->diagram->vertical);
        write_buffer (tx_params->ocl_obj,
                      0,
                      tx_params->m_antenna_loss_dev,
                      ant_buff_size,
                      tx_params->m_antenna_loss[0]);
        write_buffer_blocking (tx_params->ocl_obj,
                               0,
                               tx_params->m_radio_zone_dev,
                               rad_buff_size,
                               tx_params->m_radio_zone[0]);
    }
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
    set_kernel_int_arg    (tx_params->ocl_obj,
                           3,
                           &tx_params->ncols);
    set_kernel_double_arg (tx_params->ocl_obj,
                           4,
                           &params->rx_height_AGL);
    set_kernel_double_arg (tx_params->ocl_obj,
                           5,
                           &params->frequency);
    set_kernel_double_arg (tx_params->ocl_obj,
                           6,
                           &(tx_params->diagram->gain));
    set_kernel_int_arg    (tx_params->ocl_obj,
                           7,
                           &tx_params->beam_direction);
    set_kernel_int_arg    (tx_params->ocl_obj,
                           8,
                           &tx_params->mechanical_tilt);
    set_kernel_int_arg    (tx_params->ocl_obj,
                           9,
                           &params->main_zone_horiz);
    set_kernel_int_arg    (tx_params->ocl_obj,
                           10,
                           &params->main_zone_vert);
    set_kernel_int_arg    (tx_params->ocl_obj,
                           11,
                           &params->sec_zone_horiz);
    set_kernel_int_arg    (tx_params->ocl_obj,
                           12,
                           &params->sec_zone_vert);

    // set pointer kernel parameters
    set_kernel_mem_arg (tx_params->ocl_obj,
                        15,
                        &horiz_diag_in_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        16,
                        &vert_diag_in_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        17,
                        tx_params->m_dem_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        18,
                        tx_params->m_radio_zone_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        19,
                        tx_params->m_antenna_loss_dev);

    // reserve local memory on the device
    size_t lmem_size = _WORK_ITEMS_PER_DIMENSION_ *
                       _WORK_ITEMS_PER_DIMENSION_ *
                       (2 * sizeof (double));
    set_local_mem (tx_params->ocl_obj,
                   20,
                   lmem_size);

    //
    // calculation radius and diameter
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
    tx_data.s[3] = (double) -1.0;                       // not used

    // set kernel parameters, specific for this transmitter
    set_kernel_value_arg (tx_params->ocl_obj,
                          13,
                          sizeof (cl_double4),
                          &tx_data);
    set_kernel_value_arg (tx_params->ocl_obj,
                          14,
                          sizeof (cl_int2),
                          &tile_offset);
    //
    // define a 2D execution range for the kernel ...
    //
    size_t global_sizes [2], 
           local_sizes [2];
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
    read_buffer (tx_params->ocl_obj,
                 0,
                 tx_params->m_radio_zone_dev,
                 rad_buff_size,
                 tx_params->m_radio_zone[0]);
    read_buffer (tx_params->ocl_obj,
                 0,
                 tx_params->m_antenna_loss_dev,
                 ant_buff_size,
                 tx_params->m_antenna_loss[0]);
    read_buffer_blocking (tx_params->ocl_obj,
                          0,
                          tx_params->m_loss_dev,
                          ant_buff_size,
                          tx_params->m_loss[0]);
     */
}



/**
 * Applies (sums) the losses introduced by the antenna over the isotrophic 
 * path-loss values.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 *
 */
void
apply_antenna_influence_gpu (Parameters    *params,
                             Tx_parameters *tx_params)
{
    //
    // activate the kernel, 
    // to sum the antenna loss to the isotrophic path-loss
    // 
    activate_kernel (tx_params->ocl_obj,
                     "vector_sum_kern");

    // set pointer kernel parameters
    set_kernel_mem_arg (tx_params->ocl_obj,
                        0,
                        tx_params->m_antenna_loss_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        1,
                        tx_params->m_loss_dev);
    // reserve local memory on the device
    size_t lmem_size = _WORK_ITEMS_PER_DIMENSION_ *
                       _WORK_ITEMS_PER_DIMENSION_ *
                       sizeof (double);
    set_local_mem (tx_params->ocl_obj,
                   2,
                   lmem_size);
    //
    // define a 2D execution range for the kernel ...
    //
    size_t global_sizes [2],
           local_sizes [2];
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
 * Calculates additional gain/pathloss according to the antenna's
 * 3-dimensional diagram.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 *
 * WARNING: the output of this function overwrites 
 *          `tx_params->m_loss`, `tx_params->m_radio_zone` 
 *          and `tx_params->m_antenna_loss`
 *
 */
void
calculate_antenna_influence (Parameters    *params,
                             Tx_parameters *tx_params)
{
    //
    // do we have to load the antenna diagram?
    //
    if (tx_params->diagram == NULL)
    {
        // allocate memory for the horizontal and vertical diagrams
        tx_params->diagram = (Diagram *) malloc (sizeof (Diagram));
        tx_params->diagram->horizontal = (double *) calloc (_DIAGRAM_SIZE_, 
                                                            sizeof (double));
        tx_params->diagram->vertical = (double *) calloc (_DIAGRAM_SIZE_, 
                                                          sizeof (double));

        // get antenna's gain and directional diagrams
        char fileName [1000];
        if (strlen (params->antenna_diagram_dir) == 0)
        {
            fprintf (stderr, "ERROR Directory containing antenna files not given\n");
            exit (1);
        }
        strcpy (fileName, params->antenna_diagram_dir);
        strcat (fileName, "/");
        if (strlen (tx_params->antenna_diagram_file) == 0)
        {
            fprintf (stderr, "ERROR File name of antenna diagram not given\n");
            exit (1);
        }
        strcat (fileName, tx_params->antenna_diagram_file);

        tx_params->diagram->gain = read_antenna_diagram (fileName,
                                                         tx_params->diagram->horizontal,
                                                         tx_params->diagram->vertical);
    }
    // the horizontal and vertical resolution of one raster pixel should match
    assert (params->map_ew_res == params->map_ns_res);

    //
    // calculate the antenna influence
    //
    if (params->use_gpu)
    {
        antenna_influence_gpu (params,
                               tx_params);
        apply_antenna_influence_gpu (params,
                                     tx_params);
     }
    else
        antenna_influence_cpu (params,
                               tx_params);

    /*
     * DEBUG memory
     *
    int r,c;
    for (r=0;r<tx_params->nrows;r++)
    {
        for (c=0;c<tx_params->ncols;c++)
        {
            printf ("%.5f\t", tx_params->m_loss[r][c]);
        }
        printf ("\n");
    }*/
}

