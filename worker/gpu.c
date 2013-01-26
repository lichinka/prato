#include "worker/gpu.h"


/**
 * Initializes the OpenCL environment that enabled calculation using
 * GPU hardware on the workers, if available.
 *
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 *
 */
void
init_gpu (Tx_parameters *tx_params)
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

