#include "worker/coverage.h"
#include "worker/common_ocl.h"




/**
 * Initializes the OpenCL environment, returning a pointer to the newly
 * created kernel object.-
 *
 *
static OCLKernel* init_ocl_kernel ( )
{
    // create a new kernel object and initialize the OpenCL backend
    OCLKernel *krn = new OCLKernel ("r.hata_gpu.cl");
    krn->init ( );

    // include directory for compiling the OpenCL kernel
    std::string build_options = "-I";
    build_options += _OCL_KERNEL_INCLUDE_DIR_;
    krn->build (build_options.c_str ( ));

    return krn;
}*/



/**
 * Calculates the sector gain for each of the transmitters, based on
 * the provided antenna diagram. The output is saved in separate raster 
 * maps (i.e. one for each tx) to allow later comparisson with the CPU
 * version.
 * Returns the number of elements in each sector matrix. Because the
 * matrices are square, sqrt() gives the size of the matrix in each 
 * dimension.
 *
 * int pixel_res ............. size of one raster pixel (in meters);
 * real raster_north ......... northern limit of the raster area (in meters);
 * real raster_west .......... western limit of the raster area (in meters);
 * int nrows ................. the number of rows in the DEM;
 * int ncols ................. the number of columns in the DEM;
 * cl_real4 *test_tx ......... array containing transmitter coordinates,
 *                             and antenna height, in meters;
 * int ntest_tx .............. number of transmitters in 'test_tx';
 * real rx_height ............ receiver height above ground;
 * real frequency ............ transmitter frequency in Mhz;
 * int radius ................ calculation radius around to the transmitter
 *                             (in raster pixels);
 * int beam_dir .............. antenna beam direction;
 * int mech_tilt ............. whether antenna mechanical tilt is present;
 * char *ant_diag_file ....... file name containing the antenna diagram data;
 * int dem_fd ................ file descriptor to access DEM data;
 * char *outname ............. the base name of all generated rasters;
 * RASTER_MAP_TYPE rtype ..... DEM element data type.-
 *
 *
int sector_calculation_cl (const Parameters     *params,
                           const int pixel_res,
                           const double raster_north,
                           const double raster_west,
                           const int nrows,
                           const int ncols,
                           cl_float4* test_tx, 
                           const unsigned int ntest_tx,
                           const double rx_height, 
                           const double frequency,
                           const int radius,
                           const int beam_dir,
                           const int mech_tilt,
                           const char *ant_diag_file,
                           const int dem_fd,
                           const char *outname,
                           float rtype)
{
    unsigned int      i;
    cl_int            status;
    cl_context        context;
    cl_command_queue *list_queues;
    uint              num_queues = 1;
    cl_event          events [num_queues];
    cl_kernel         kernel;

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
                    sizeof (int),
                    &pixel_res);
    set_kernel_arg (&kernel, 
                    1,
                    sizeof (double),
                    &raster_north);
    set_kernel_arg (&kernel, 
                    2,
                    sizeof (double),
                    &raster_west);
    set_kernel_arg (&kernel, 
                    3,
                    sizeof (int),
                    &ncols);
    set_kernel_arg (&kernel, 
                    6,
                    sizeof (double),
                    &rx_height);
    set_kernel_arg (&kernel, 
                    7,
                    sizeof (double),
                    &frequency);
    set_kernel_arg (&kernel, 
                    8,
                    sizeof (double),
                    &antenna_gain);
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
    const size_t global_range_0 = global_sizes[0];
    const size_t global_range_1 = global_sizes[1];

    // create OpenCL buffers
    cl_mem horiz_diag_in_dev;
    cl_mem vert_diag_in_dev;
    cl_mem dem_in_dev;
    cl_mem sect_out_dev;
    size_t diagram_size = 360;
    size_t dem_in_size = params->nrows * 
                         params->ncols * 
                         sizeof (params->m_dem[0][0]);
    size_t sect_out_size = params->nrows * 
                           params->ncols * 
                           sizeof (params->m_loss[0][0]);
    
    create_buffer (&context,
                   CL_MEM_READ_ONLY, 
                   diagram_size * sizeof (float),
                   &horiz_diag_in_dev);
    create_buffer (&context,
                   CL_MEM_READ_ONLY, 
                   diagram_size * sizeof (float),
                   &vert_diag_in_dev);
    create_buffer (&context,
                   CL_MEM_READ_ONLY, 
                   dem_in_size,
                   params->m_dem,
                   &dem_in_dev);
    create_buffer (&context,
                   CL_MEM_WRITE_ONLY, 
                   sect_out_size,
                   params->m_loss,
                   &sect_out_dev);

    /* send input data to the device
    krn->write_buffer (horiz_diag_in_dev, horiz_diag, 
                       diagram_size * sizeof (float));
    krn->write_buffer (vert_diag_in_dev, vert_diag, 
                       diagram_size * sizeof (float));
    krn->write_buffer (dem_in_dev, dem_in, dem_in_size);

    // set the remaining parameters
    krn->set_arg (11, horiz_diag_in_dev);
    krn->set_arg (12, vert_diag_in_dev);
    krn->set_arg (13, dem_in_dev);
    krn->set_arg (14, sect_out_dev);

    G_message ("Simulating antenna influence per sector for %d transmiter(s) ...", ntest_tx);

    // for each transmitter ...
    for (i = 0; i < ntest_tx; i++)
    {
        // display completion percentage
        G_percent (i, ntest_tx, 1);

        // tile offset within the calculation radius,
        // given in pixel coordinates
        cl_int2 tx_offset;
        tx_offset.s[0] = (int) ((test_tx[i].s[0] - raster_west) / pixel_res - radius);
        tx_offset.s[1] = (int) ((raster_north - test_tx[i].s[1]) / pixel_res - radius);
        // transmitter data
        cl_real4 tx_data;
        tx_data.s[0] = (real) test_tx[i].s[0];          // transmitter coordinate
        tx_data.s[1] = (real) test_tx[i].s[1];          // transmitter coordinate
        tx_data.s[2] = (real) dem_in[(tx_offset.s[1]+radius)*ncols + tx_offset.s[0] + radius]; // height above sea level
        tx_data.s[3] = (real) test_tx[i].s[2];          // antenna height

        G_message ("Transmitter located at %f,%f,%f,%f", tx_data.s[0],
                                                         tx_data.s[1],
                                                         tx_data.s[2],
                                                         tx_data.s[3]);
        // set kernel parameters, specific for this transmitter
        krn->set_arg (4, tx_data);
        krn->set_arg (5, tx_offset);

        // start kernel execution 
        krn->run ( );

        // sync memory
        krn->read_buffer (sect_out_dev, sect_out, sect_out_size);

        // write the calculation output
        sprintf (raster_outname, "%s", outname);
        save_raster_to_file (raster_outname, sect_out,
                             tx_offset.s[0], tx_offset.s[1],
                             global_range_0,
                             global_range_1,
                             ncols, nrows);
    }
    // deactivate the OpenCL environment
    release_opencl (num_queues, 
                    &list_queues, 
                    &context);

    // free allocated resources
    delete krn;
    free (sect_out);
    free (dem_in);

    // number of elements in each sector matrix
    return (int) (sect_out_size / G_raster_size (rtype));
}
*/

