#include "worker/gpu.h"


/**
 * Initializes the OpenCL environment that enabled calculation using
 * GPU hardware on the workers, if available.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration
 *              parameters;
 *
 */
void
init_gpu (Parameters *params,
          Tx_parameters *tx_params)
{
    //
    // initialize the OpenCL environment if we haven't already;
    //
    if (tx_params->ocl_obj == NULL)
    {
        //
        // create the memory objects on the device, that will be shared 
        // among different function calls; this should minize data 
        // transfers to/from the GPU
        //
        tx_params->ocl_obj           = (OCL_object *) malloc (sizeof (OCL_object));
        tx_params->m_dem_dev         = (cl_mem *) malloc (sizeof (cl_mem));
        tx_params->m_clut_dev        = (cl_mem *) malloc (sizeof (cl_mem));
        tx_params->m_loss_dev        = (cl_mem *) malloc (sizeof (cl_mem));
        tx_params->m_obst_height_dev = (cl_mem *) malloc (sizeof (cl_mem));
        tx_params->m_obst_dist_dev   = (cl_mem *) malloc (sizeof (cl_mem));
        tx_params->v_clutter_loss_dev= (cl_mem *) malloc (sizeof (cl_mem));

        // initialize the OpenCL platform
        init_opencl (tx_params->ocl_obj, 1);

        //
        // create OpenCL buffers
        //
        int nelem = tx_params->nrows * tx_params->ncols;

        *(tx_params->m_dem_dev)         = create_buffer (tx_params->ocl_obj,
                                                         CL_MEM_READ_ONLY, 
                                                         nelem * sizeof (tx_params->m_dem[0][0]));
        *(tx_params->m_clut_dev)        = create_buffer (tx_params->ocl_obj,
                                                         CL_MEM_READ_ONLY, 
                                                         nelem * sizeof (tx_params->m_clut[0][0]));
        *(tx_params->m_loss_dev)        = create_buffer (tx_params->ocl_obj,
                                                         CL_MEM_READ_WRITE, 
                                                         nelem * sizeof (tx_params->m_loss[0][0]));
        *(tx_params->m_obst_height_dev) = create_buffer (tx_params->ocl_obj,
                                                         CL_MEM_READ_ONLY, 
                                                         nelem * sizeof (tx_params->m_obst_height[0][0]));
        *(tx_params->m_obst_dist_dev)   = create_buffer (tx_params->ocl_obj,
                                                         CL_MEM_READ_ONLY, 
                                                         nelem * sizeof (tx_params->m_obst_dist[0][0]));
        *(tx_params->v_clutter_loss_dev)= create_buffer (tx_params->ocl_obj,
                                                         CL_MEM_READ_ONLY,
                                                         params->clutter_category_count * sizeof (params->clutter_loss));

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
        exit (1);
    }
    if ((diameter_in_pixels % _WORK_ITEMS_PER_DIMENSION_) != 0)
    {
        fprintf (stderr, 
                 "*** ERROR: Cannot execute on GPU. Try setting a calculation radius multiple of %d.\n",
                 _WORK_ITEMS_PER_DIMENSION_);
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

