#include "worker/gpu.h"




/**
 * Deactivates the OpenCL environment, releasing all associated memory.
 *
 * tx_params    a structure holding transmitter-specific configuration
 *              parameters;
 *
 */
void
release_gpu (Tx_parameters *tx_params)
{ 
    //
    // check the OpenCL environment has been initialized
    //
    if (tx_params->ocl_obj != NULL)
    {
        //
        // release the memory objects on the device
        //
        if (tx_params->m_dem_dev != NULL)
        {
            release_buffer (tx_params->m_dem_dev);
            free (tx_params->m_dem_dev);
            tx_params->m_dem_dev = NULL;
        }
        
        if (tx_params->m_clut_dev != NULL)
        {
            release_buffer (tx_params->m_clut_dev);
            free (tx_params->m_clut_dev);
            tx_params->m_clut_dev = NULL;
        }
        
        if (tx_params->m_field_meas_dev != NULL)
        {
            release_buffer (tx_params->m_field_meas_dev);
            free (tx_params->m_field_meas_dev);
            tx_params->m_field_meas_dev = NULL;
        }
        
        if (tx_params->m_loss_dev != NULL)
        {
            release_buffer (tx_params->m_loss_dev);
            free (tx_params->m_loss_dev);
            tx_params->m_loss_dev = NULL;
        }
        
        if (tx_params->m_antenna_loss_dev != NULL)
        {
            release_buffer (tx_params->m_antenna_loss_dev);
            free (tx_params->m_antenna_loss_dev);
            tx_params->m_antenna_loss_dev = NULL;
        }
        
        if (tx_params->m_radio_zone_dev != NULL)
        {
            release_buffer (tx_params->m_radio_zone_dev);
            free (tx_params->m_radio_zone_dev);
            tx_params->m_radio_zone_dev = NULL;
        }
        
        if (tx_params->m_obst_height_dev != NULL)
        {
            release_buffer (tx_params->m_obst_height_dev);
            free (tx_params->m_obst_height_dev);
            tx_params->m_obst_height_dev = NULL;
        }
        
        if (tx_params->m_obst_dist_dev != NULL)
        {
            release_buffer (tx_params->m_obst_dist_dev);
            free (tx_params->m_obst_dist_dev);
            tx_params->m_obst_dist_dev = NULL;
        }

        if (tx_params->v_partial_sum_dev != NULL)
        {
            release_buffer (tx_params->v_partial_sum_dev);
            free (tx_params->v_partial_sum_dev);
            tx_params->v_partial_sum_dev = NULL;
        }
        
        if (tx_params->v_clutter_loss_dev != NULL)
        {
            release_buffer (tx_params->v_clutter_loss_dev);
            free (tx_params->v_clutter_loss_dev);
            tx_params->v_clutter_loss_dev = NULL;
        }
    }
}



/**
 * Deactivates the OpenCL environment, releasing all associated memory.
 *
 * tx_params    a structure holding transmitter-specific configuration
 *              parameters;
 *
 */
void
close_gpu (Tx_parameters *tx_params)
{
    //
    // check the OpenCL environment has been initialized
    //
    if (tx_params->ocl_obj != NULL)
    {
        release_gpu (tx_params);
        deactivate_opencl (tx_params->ocl_obj);
    }
}



/**
 * Initializes the OpenCL environment that enabled calculation using
 * GPU hardware on the workers, if available.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 * device_hint     sends a hint to the OpenCL backend about which device
 *                  id to select and use; useful for using multiple GPUs.-
 *
 */
void
init_gpu (Parameters    *params,
          Tx_parameters *tx_params,
          const int     device_hint)
{
    //
    // initialize the OpenCL environment if we haven't already;
    //
    if (tx_params->ocl_obj == NULL)
    {
        tx_params->ocl_obj = init_opencl (tx_params->ocl_obj, 
                                          1,
                                          device_hint);
        //
        // build the OpenCL source file (only the first time);
        // in this case, all kernels reside in one source file
        //
        build_kernel_from_file (tx_params->ocl_obj,
                                "r.coverage.cl",
                                "eric_per_tx",
                                "-I. -Werror");
    }
    //
    // create the memory objects on the device, that will be shared 
    // among different function calls; this should minize data 
    // transfers to/from the GPU
    //
    int nelem = tx_params->nrows * tx_params->ncols;

    if (tx_params->m_dem_dev == NULL)
    {
        tx_params->m_dem_dev    = (cl_mem *) malloc (sizeof (cl_mem));
        *(tx_params->m_dem_dev) = create_buffer (tx_params->ocl_obj,
                                                 CL_MEM_READ_ONLY, 
                                                 nelem * sizeof (tx_params->m_dem[0][0]));
    }
    if (tx_params->m_clut_dev == NULL)
    {
        tx_params->m_clut_dev    = (cl_mem *) malloc (sizeof (cl_mem));
        *(tx_params->m_clut_dev) = create_buffer (tx_params->ocl_obj,
                                                  CL_MEM_READ_ONLY, 
                                                  nelem * sizeof (tx_params->m_clut[0][0]));
    }
    if (tx_params->m_loss_dev == NULL)
    {
        tx_params->m_loss_dev    = (cl_mem *) malloc (sizeof (cl_mem));
        *(tx_params->m_loss_dev) = create_buffer (tx_params->ocl_obj,
                                                  CL_MEM_READ_WRITE, 
                                                  nelem * sizeof (tx_params->m_loss[0][0]));
    }
    if (tx_params->m_obst_height_dev == NULL)
    {
        tx_params->m_obst_height_dev    = (cl_mem *) malloc (sizeof (cl_mem));
        *(tx_params->m_obst_height_dev) = create_buffer (tx_params->ocl_obj,
                                                         CL_MEM_READ_ONLY, 
                                                         nelem * sizeof (tx_params->m_obst_height[0][0]));
    }
    if (tx_params->m_obst_dist_dev == NULL)
    {
        tx_params->m_obst_dist_dev    = (cl_mem *) malloc (sizeof (cl_mem));
        *(tx_params->m_obst_dist_dev) = create_buffer (tx_params->ocl_obj,
                                                       CL_MEM_READ_ONLY, 
                                                       nelem * sizeof (tx_params->m_obst_dist[0][0]));
    }
    if (tx_params->v_clutter_loss_dev == NULL)
    {
        tx_params->v_clutter_loss_dev    = (cl_mem *) malloc (sizeof (cl_mem));
        *(tx_params->v_clutter_loss_dev) = create_buffer (tx_params->ocl_obj,
                                                          CL_MEM_READ_ONLY,
                                                          params->clutter_category_count * sizeof (params->clutter_loss));
    }
    //
    // send input data to the device
    // 
    write_buffer (tx_params->ocl_obj,
                  0,
                  tx_params->m_obst_height_dev,
                  nelem * sizeof (tx_params->m_obst_height[0][0]),
                  tx_params->m_obst_height[0]);
    write_buffer (tx_params->ocl_obj,
                  0,
                  tx_params->m_obst_dist_dev,
                  nelem * sizeof (tx_params->m_obst_dist[0][0]),
                  tx_params->m_obst_dist[0]);
    write_buffer (tx_params->ocl_obj,
                  0,
                  tx_params->m_dem_dev,
                  nelem * sizeof (tx_params->m_dem[0][0]),
                  tx_params->m_dem[0]);
    write_buffer (tx_params->ocl_obj,
                  0,
                  tx_params->m_clut_dev,
                  nelem * sizeof (tx_params->m_clut[0][0]),
                  tx_params->m_clut[0]);
    write_buffer_blocking (tx_params->ocl_obj,
                           0,
                           tx_params->m_loss_dev,
                           nelem * sizeof (tx_params->m_loss[0][0]),
                           tx_params->m_loss[0]);
}


/**
 * Defines a 2D execution range for a kernel, consisting of square tiles,
 * which size is based on the execution capabilities of the available GPU 
 * hardware, i.e. the number of concurrent threads it can handle.
 *
 * params   a structure holding configuration parameters which are 
 *          common to all transmitters;
 * global_sizes the two-dimensional global execution range (output parameter);
 * local_sizes  the two-dimensional local execution range (output parameter).-
 *
 */
void 
define_2D_range (const Parameters *params,
                 size_t *global_sizes,
                 size_t *local_sizes)
{
    double radius_in_meters = params->radius * 1000;
    int radius_in_pixels    = (int) (radius_in_meters / params->map_ew_res);
    int diameter_in_pixels  = 2 * radius_in_pixels;
    size_t ntile            = diameter_in_pixels / _WORK_ITEMS_PER_DIMENSION_;

    //
    // number of processing tiles, based on the received radius 
    //
    if (diameter_in_pixels < _WORK_ITEMS_PER_DIMENSION_)
    {
        fprintf (stderr, 
                 "*** ERROR: Cannot execute on GPU. Increase the calculation radius and try again.\n");
        fflush (stderr);
        exit (1);
    }
    if ((diameter_in_pixels % _WORK_ITEMS_PER_DIMENSION_) != 0)
    {
        fprintf (stderr, 
                 "*** ERROR: Cannot execute on GPU. Try setting a calculation radius multiple of %d.\n",
                 _WORK_ITEMS_PER_DIMENSION_);
        fflush (stderr);
        exit (1);
    }

    //
    // define the execution range with global and local dimensions
    //
    global_sizes[0] = ntile * _WORK_ITEMS_PER_DIMENSION_;
    global_sizes[1] = ntile * _WORK_ITEMS_PER_DIMENSION_;
    local_sizes[0]  = _WORK_ITEMS_PER_DIMENSION_;
    local_sizes[1]  = _WORK_ITEMS_PER_DIMENSION_;
}

