#include <assert.h>
//
// From the performance metrics library
//  
//      https://github.com/lichinka/performance_metrics
//
#include <performance_metric.h>

//
// From the OpenCL common library
//
//      https://github.com/lichinka/ocl_common
//
#include <ocl_common.h>

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
    // finally take into account pathloss for determined diagram angles and antenna gain 
    //
    *antenna_loss = horizontal_loss + vertical_loss - diagram->gain;
}



/**
 * Standard CPU version of the antenna influence algorithm.
 * For a description of the parameters used, see the function
 * 'calculate_antenna_influence'.-
 *
 */
static void 
antenna_influence_cpu (const double tx_east_coord,
                       const double tx_north_coord,
                       const double antenna_height_AGL,
                       const double total_tx_height,
                       const int beam_direction,
                       const int mechanical_tilt,
                       const double frequency,
                       const double radius,  
                       const double rx_height_AGL,
                       const int nrows,       
                       const int ncols,      
                       const double map_west,
                       const double map_north,
                       const double map_ew_res,  
                       const double map_ns_res,  
                       const float  null_value,
                       const int main_zone_horiz,
                       const int main_zone_vert,
                       const int sec_zone_horiz,
                       const int sec_zone_vert,
                       double **m_dem,          
                       double **m_loss,
                       char   **m_radio_zone,
                       double **m_antenna_loss)
{
    int r, c;

    //
    // this loop measures around 281.000 MFlops
    // 
    for (r = 0; r < nrows; r++) 
    {
        float loss_out, dem_in; 

        // calculate receiver coordinates
        double rec_north = map_north - (r * map_ns_res);

        // calculate differences between receiver and transmitter coordinates
        double d_north = rec_north - tx_north_coord;

        //
        // ... and each column 
        // 
        for (c = 0; c < ncols; c++) 
        { 
#ifdef _PERFORMANCE_METRICS_
            memory_access (2, 8);
#endif
            dem_in = (float) m_dem[r][c];

            // calculate receiver coordinates
            double rec_east = map_west + (c * map_ew_res);

            // calculate differences between receiver and transmitter coordinates
            double d_east = rec_east - tx_east_coord;
            
            // calculate distance between Tx and Rx (in kilometers)
            double dist_Tx_Rx = sqrt (pow (d_east, 2) + 
                                      pow (d_north, 2));
            dist_Tx_Rx /= 1000;
       
            //
            // ignore pixels exceeding radius distance between Rx and Tx or
            // if they are too close (less than 200 mts)
            // 
            if ((dist_Tx_Rx < 0.2) || (dist_Tx_Rx > radius))
            {
                loss_out            = null_value;
                m_radio_zone[r][c] &= _RADIO_ZONE_MODEL_DISTANCE_OFF_;
            }
            else 
            {
                m_radio_zone[r][c] |= _RADIO_ZONE_MODEL_DISTANCE_ON_;
                antenna_influence_on_point (d_east,
                                            d_north,
                                            total_tx_height,
                                            beam_direction,
                                            mechanical_tilt,
                                            dem_in,
                                            rx_height_AGL,
                                            dist_Tx_Rx,
                                            diagram,
                                            main_zone_horiz,
                                            main_zone_vert,
                                            sec_zone_horiz,
                                            sec_zone_vert,
                                            &(m_radio_zone[r][c]),
                                            &(m_antenna_loss[r][c]));
                loss_out = (float) (m_loss[r][c] + m_antenna_loss[r][c]);
            }
            // 
            // save the result in the output matrix
            //
            m_loss[r][c] = (double) loss_out;
        }
    }
}





/**
 * GPU version of the antenna influence algorithm.
 * For a description of the parameters used, see the function
 * 'calculate_antenna_influence'.
 *
 */
static void 
antenna_influence_gpu (const double tx_east_coord,
                       const double tx_north_coord,
                       const double antenna_height_AGL,
                       const double total_tx_height,
                       const int beam_direction,
                       const int mechanical_tilt,
                       const double frequency,
                       const double radius,  
                       const double rx_height_AGL,
                       const int nrows,       
                       const int ncols,      
                       const double map_west,
                       const double map_north,
                       const double map_ew_res,  
                       const double map_ns_res,  
                       const float  null_value,
                       GPU_parameters *gpu_params,
                       double **m_dem,          
                       double **m_loss)
{
    // build and activate the kernel
    build_kernel_from_file (gpu_params->ocl_obj,
                            "sector_kern",
                            "r.coverage.cl",
                            "-I. -w");
    //
    // create OpenCL buffers
    //
    size_t diagram_size = 360 *
                          sizeof (diagram->horizontal[0]);
    size_t buffer_size = nrows * 
                         ncols * 
                         sizeof (double);

    cl_mem horiz_diag_in_dev = create_buffer (gpu_params->ocl_obj,
                                              CL_MEM_READ_ONLY, 
                                              diagram_size);
    cl_mem vert_diag_in_dev = create_buffer (gpu_params->ocl_obj,
                                             CL_MEM_READ_ONLY, 
                                             diagram_size);
    //
    // send input data to the device
    //
    cl_event *event;
    write_buffer (gpu_params->ocl_obj,
                  0,
                  &horiz_diag_in_dev,
                  diagram_size,
                  diagram->horizontal);
    event = write_buffer (gpu_params->ocl_obj,
                          0,
                          &vert_diag_in_dev,
                          diagram_size,
                          diagram->vertical);
    clWaitForEvents (1, event);

    // kernel parameters 
    set_kernel_double_arg (gpu_params->ocl_obj,
                           0,
                           &map_ew_res);
    set_kernel_double_arg (gpu_params->ocl_obj,
                           1,
                           &map_north);
    set_kernel_double_arg (gpu_params->ocl_obj,
                           2,
                           &map_west);
    set_kernel_int_arg (gpu_params->ocl_obj,
                        3,
                        &ncols);
    set_kernel_double_arg (gpu_params->ocl_obj,
                           6,
                           &rx_height_AGL);
    set_kernel_double_arg (gpu_params->ocl_obj,
                           7,
                           &frequency);
    set_kernel_double_arg (gpu_params->ocl_obj,
                           8,
                           &(diagram->gain));
    set_kernel_int_arg (gpu_params->ocl_obj,
                        9,
                        &beam_direction);
    set_kernel_int_arg (gpu_params->ocl_obj,
                        10,
                        &mechanical_tilt);

    // reserve local memory on the device
    size_t lmem_size = _WORK_ITEMS_PER_DIMENSION_ *
                       _WORK_ITEMS_PER_DIMENSION_ *
                       sizeof (cl_float2);
    set_local_mem (gpu_params->ocl_obj,
                   15,
                   lmem_size);

    // set the remaining parameters
    set_kernel_mem_arg (gpu_params->ocl_obj,
                        11,
                        &horiz_diag_in_dev);
    set_kernel_mem_arg (gpu_params->ocl_obj,
                        12,
                        &vert_diag_in_dev);
    set_kernel_mem_arg (gpu_params->ocl_obj,
                        13,
                        gpu_params->m_dem_dev);
    set_kernel_mem_arg (gpu_params->ocl_obj,
                        14,
                        gpu_params->m_loss_dev);

    //
    // calculation radius and diameter
    //
    double radius_in_meters = radius * 1000;
    int radius_in_pixels    = (int) (radius_in_meters / map_ew_res);
    int diameter_in_pixels  = 2 * radius_in_pixels;

    //
    // calculation tile offset within the target area, given in pixel coordinates
    //
    cl_int2 tile_offset;
    tile_offset.s[0]  = (int) ((tx_east_coord - map_west) - radius_in_meters);
    tile_offset.s[0] /= map_ew_res;
    tile_offset.s[1]  = (int) ((map_north - tx_north_coord) - radius_in_meters);
    tile_offset.s[1] /= map_ns_res;
    //
    // transmitter data
    //
    cl_double4 tx_data;
    tx_data.s[0] = (double) tx_east_coord;   // transmitter coordinate
    tx_data.s[1] = (double) tx_north_coord;  // transmitter coordinate
    tx_data.s[2] = (double) total_tx_height; // antenna height above sea level
    tx_data.s[3] = (double) -1.0;            // not used

    // set kernel parameters, specific for this transmitter
    set_kernel_value_arg (gpu_params->ocl_obj,
                          4,
                          sizeof (cl_double4),
                          &tx_data);
    set_kernel_value_arg (gpu_params->ocl_obj,
                          5,
                          sizeof (cl_int2),
                          &tile_offset);

    // number of processing tiles needed around each transmitter 

    if (diameter_in_pixels < _WORK_ITEMS_PER_DIMENSION_)
    {
        fprintf (stderr, 
                 "ERROR Cannot execute on GPU. Increase the calculation radius and try again.\n");
        exit (1);
    }
    if ((diameter_in_pixels % _WORK_ITEMS_PER_DIMENSION_) != 0)
    {
        fprintf (stderr, 
                 "ERROR Cannot execute on GPU. Try to set a calculation radius multiple of %d\n",
                 _WORK_ITEMS_PER_DIMENSION_);
        exit (1);
    }
    //
    // define a 2D execution range for the kernel ...
    //
    size_t ntile = diameter_in_pixels / _WORK_ITEMS_PER_DIMENSION_;
    size_t global_sizes [] = {ntile * _WORK_ITEMS_PER_DIMENSION_,
                              ntile * _WORK_ITEMS_PER_DIMENSION_};
    size_t local_sizes [] = {_WORK_ITEMS_PER_DIMENSION_,
                             _WORK_ITEMS_PER_DIMENSION_};

    //
    // ... and execute it
    //
    run_kernel_2D_blocking (gpu_params->ocl_obj,
                            0,
                            NULL,
                            global_sizes,
                            local_sizes);
    //
    // sync memory
    //
    read_buffer_blocking (gpu_params->ocl_obj,
                          0,
                          gpu_params->m_loss_dev,
                          buffer_size,
                          m_loss[0]);
}




/**
 * Calculates additional gain/pathloss according to the antenna's
 * 3-dimensional diagram.
 *
 * use_gpu              A flag indicating whether to use GPU hardware;
 * tx_east_coord        eastern coordinate of the transmitter;
 * tx_north_coord       northern coordinate of the transmitter;
 * antenna_height_AGL   height of the transmitter above ground level;
 * total_tx_height      height of the transmitter above sea level;
 * beam_direction       direction of the antenna beam, ie azimuth, in degrees;
 * mechanical_tilt      mechanical antenna tilt angle, in degress;
 * frequency            transmitter frequency, in Mhz;
 * radius               calculation radius around the given transmitter, in km;
 * rx_height_AGL        receiver's height above ground level;
 * nrows                number of rows within the 2D-matrices;
 * ncols                number of columns within the 2D-matrices;
 * map_west             western coordinate of the given maps, ie matrices;
 * map_north            northern coordinate of the given maps, ie matrices;
 * map_ew_res           east/west map resolution;
 * map_ns_res           north/south map resolution;
 * null_value           value representing NULL data on a map;
 * antenna_diagram_dir  directory containing the files describing the antenna
 *                      diagrams;
 * antenna_diagram_file file containing the antenna diagram description;
 * main_zone_horiz      indicates the horizontal loss that defines the main
 *                      antenna beam;
 * main_zone_vert       indicates the vertical loss that defines the main
 *                      antenna beam;
 * sec_zone_horiz       indicates the horizontal loss that defines the secondary
 *                      antenna beam;
 * sec_zone_vert        indicates the vertical loss that defines the secondary
 *                      antenna beam;
 * gpu_params           parameters for GPU-based execution;
 * m_dem                a 2D-matrix containing the digital elevation model;
 * m_loss               a 2D-matrix containing the path-loss values from the
 *                      isotrophic prediction (output parameter);
 * m_radio_zone         a 2D-matrix containing bit masks that indicate the
 *                      radio zone to which each point belongs 
 *                      (output parameter);
 * m_antenna_loss       a 2D-matrix containing the loss, introduced by the 
 *                      antenna, on every point (output parameter);
 *
 * WARNING: the output of this function overwrites 
 *          `m_loss`, `m_radio_zone` and `m_antenna_loss`
 *
 */
void
calculate_antenna_influence (const char use_gpu,
                             const double tx_east_coord,
                             const double tx_north_coord,
                             const double antenna_height_AGL,
                             const double total_tx_height,
                             const int beam_direction,
                             const int mechanical_tilt,
                             const double frequency,
                             const double radius,  
                             const double rx_height_AGL,
                             const int nrows,       
                             const int ncols,      
                             const double map_west,
                             const double map_north,
                             const double map_ew_res,  
                             const double map_ns_res,  
                             const float  null_value,
                             const char *antenna_diagram_dir,
                             const char *antenna_diagram_file,
                             const int main_zone_horiz,
                             const int main_zone_vert,
                             const int sec_zone_horiz,
                             const int sec_zone_vert,
                             GPU_parameters *gpu_params,
                             double **m_dem,          
                             double **m_loss,
                             char   **m_radio_zone,
                             double **m_antenna_loss)
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
        if (strlen (antenna_diagram_dir) == 0)
        {
            fprintf (stderr, "ERROR Directory containing antenna files not given\n");
            exit (1);
        }
        strcpy (fileName, antenna_diagram_dir);
        strcat (fileName, "/");
        if (strlen (antenna_diagram_file) == 0)
        {
            fprintf (stderr, "ERROR File name of antenna diagram not given\n");
            exit (1);
        }
        strcat (fileName, antenna_diagram_file);

        diagram->gain = read_antenna_diagram (fileName,
                                              diagram->horizontal,
                                              diagram->vertical);
    }
    // the horizontal and vertical resolution of one raster pixel should match
    assert (map_ew_res == map_ns_res);

    //
    // calculate the antenna influence
    //
    if (use_gpu)
        antenna_influence_gpu (tx_east_coord,
                               tx_north_coord,
                               antenna_height_AGL,
                               total_tx_height,
                               beam_direction,
                               mechanical_tilt,
                               frequency,
                               radius,  
                               rx_height_AGL,
                               nrows,       
                               ncols,      
                               map_west,
                               map_north,
                               map_ew_res,  
                               map_ns_res,  
                               null_value,
                               gpu_params,
                               m_dem,          
                               m_loss);
    else
        antenna_influence_cpu (tx_east_coord,
                               tx_north_coord,
                               antenna_height_AGL,
                               total_tx_height,
                               beam_direction,
                               mechanical_tilt,
                               frequency,
                               radius,  
                               rx_height_AGL,
                               nrows,       
                               ncols,      
                               map_west,
                               map_north,
                               map_ew_res,  
                               map_ns_res,  
                               null_value,
                               main_zone_horiz,
                               main_zone_vert,
                               sec_zone_horiz,
                               sec_zone_vert,
                               m_dem,          
                               m_loss,
                               m_radio_zone,
                               m_antenna_loss);
    /*
     * DEBUG memory
     *
    int r,c;
    for (r=0;r<nrows;r++)
    {
        for (c=0;c<ncols;c++)
        {
            printf ("%.5f\t", m_loss[r][c]);
        }
        printf ("\n");
    }*/
}

