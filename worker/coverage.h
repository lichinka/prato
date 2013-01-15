#ifndef _COVERAGE_CALCULATION_IN_GRASS_H_
#define _COVERAGE_CALCULATION_IN_GRASS_H_

#define _CHAR_BUFFER_SIZE_          1024
#define _COVERAGE_MASTER_RANK_      0
#define _WORKER_IS_IDLE_TAG_        100
#define _WORKER_KEEP_WORKING_TAG_   105
#define _WORKER_SHUTDOWN_TAG_       110

#include <mpi.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <grass/gis.h>
#include <grass/glocale.h>

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

#include "worker/ini.h"



//
// A structure to hold all transmitter-specific configuration parameters
//
struct Tx_parameters
{
    //
    // transmitter-specific parameters
    //
    int     beam_direction;
    int     electrical_tilt;
    int     mechanical_tilt;
    double  antenna_height_AGL;
    double  tx_east_coord;
    double  tx_north_coord;
    double  tx_power;
    double  total_tx_height;
    char    tx_name              [_CHAR_BUFFER_SIZE_];
    char    antenna_diagram_file [_CHAR_BUFFER_SIZE_];
} __attribute__((__packed__));

typedef struct Tx_parameters Tx_parameters;

/**
 * Handles name=value pairs read from the INI file inside 'tx_section'.
 * This is useful for having only one INI file with all the transmitters
 * in it, selecting which one to include with a command line parameter.
 *
 * user_struct  The structure in which the read parameters are kept;
 * section      the current section being read from the INI file;
 * name         the name part of the current line in the INI file;
 * value        the value part of the current line in the INI file;
 * tx_section   the name of the transmitter's section to be parsed.-
 *
 */
static int tx_params_handler (void *user_struct, 
                              const char *section, 
                              const char *name,
                              const char *value,
                              const char *tx_section)
{
    Tx_parameters *pconfig = (Tx_parameters *) user_struct;

    #define MATCH(s,n) strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0

    if (MATCH (tx_section, "cellName"))
        strncpy (pconfig->tx_name, value, _CHAR_BUFFER_SIZE_);
    else if (MATCH (tx_section, "beamDirection"))
    {
        pconfig->beam_direction = atoi (value);
    }
    else if (MATCH (tx_section, "electricalTiltAngle"))
    {
        pconfig->electrical_tilt = atoi (value);
    }
    else if (MATCH (tx_section, "mechanicalTiltAngle"))
    {
        pconfig->mechanical_tilt = atoi (value);
    }
    else if (MATCH (tx_section, "heightAGL"))
    {
        pconfig->antenna_height_AGL = atof (value);
    }
    else if (MATCH (tx_section, "antennaFile"))
        strncpy (pconfig->antenna_diagram_file, value, _CHAR_BUFFER_SIZE_);
    else if (MATCH (tx_section, "positionEast"))
    {
        pconfig->tx_east_coord = atof (value);
    }
    else if (MATCH (tx_section, "positionNorth"))
    {
        pconfig->tx_north_coord = atof (value);
    }
    else if (MATCH (tx_section, "power"))
    {
        pconfig->tx_power = atof (value);
    }
    else 
        return 0;  /* unknown section/name, error */
    return 1;
}

//
// A structure to hold all GPU-specific data
//
struct GPU_parameters
{
    OCL_objects *ocl_obj;
    cl_mem      *m_dem_dev;
    cl_mem      *m_clut_dev;
    cl_mem      *m_loss_dev;
} __attribute__((__packed__));

typedef struct GPU_parameters GPU_parameters;

//
// A structure to hold all common configuration parameters and
// in-run data for the coverage-calculation process
//
struct Parameters
{
    //
    // run-time parameters
    //

    // number of transmitters being processed
    int ntx;

    // parameters of the transmitters being processed
    Tx_parameters *tx_params;

    // NULL value used in maps
    float null_value;

    // 2D matrix containing the digital-elevation-model raster
    double **m_dem;

    // 2D matrix containing the clutter information of the area
    double **m_clut;
    
    // 2D area matrix where the path-loss predictions are saved
    double **m_loss;
    
    // 2D area matrix where the field measurements are saved
    double **m_field_meas;

    //
    // common parameters to all transmitters
    //
    double  radius;
    double  rx_height_AGL;
    double  frequency;
    int     nrows;
    int     ncols;
    double  map_east;
    double  map_west;
    double  map_north;
    double  map_south;
    double  map_ew_res;
    double  map_ns_res;
    char    dem_map             [_CHAR_BUFFER_SIZE_];
    char    clutter_map         [_CHAR_BUFFER_SIZE_];
    char    antenna_diagram_dir [_CHAR_BUFFER_SIZE_];

    // a buffer to save the INI file to memory
    char *  ini_file_content;
    int     ini_file_content_size;

    // a flag to indicate the GPU should be used
    char            use_gpu;
    // GPU-specific data is kept in this structure
    GPU_parameters *gpu_params;
} __attribute__((__packed__));

typedef struct Parameters Parameters;

/**
 * Handles name=value pairs read from the INI file in the [Common] section.
 * Therefore, the 'tx_section' parameter is not taken into account here.
 *
 * user_struct  The structure in which the read parameters are kept;
 * section      the current section being read from the INI file;
 * name         the name part of the current line in the INI file;
 * value        the value part of the current line in the INI file;
 * tx_section   ignored in this context.-
 *
 */
static int common_params_handler (void *user_struct, 
                                  const char *section, 
                                  const char *name,
                                  const char *value,
                                  const char *tx_section)
{
    Parameters *pconfig = (Parameters *) user_struct;

    #define MATCH(s,n) strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0

    if (MATCH ("common", "DEMMapName"))
        strncpy (pconfig->dem_map, value, _CHAR_BUFFER_SIZE_);
    else if (MATCH ("common", "clutterMapName"))
        strncpy (pconfig->clutter_map, value, _CHAR_BUFFER_SIZE_);
    else if (MATCH ("common", "receiverHeightAGL"))
        pconfig->rx_height_AGL = atof (value);
    else if (MATCH ("common", "frequency"))
        pconfig->frequency = atof (value);
    else if (MATCH ("common", "radius"))
        pconfig->radius = atof (value);
    else if (MATCH ("common", "antennaDirectory"))
        strncpy (pconfig->antenna_diagram_dir, value, _CHAR_BUFFER_SIZE_);
    else
        return 0;  /* unknown section/name, error */
    return 1;
}




/**
 * Initializes the coverage calculation by reading the configuration 
 * parameters in the [Common] section of the INI file passed as argument.
 * This function returns a pointer to the newly created parameters structure.
 *
 * ini_file         Pointer to the stream containing the configuration read
 *                  from the INI file;
 * tx_sections_list the names of the sections containing transmitter-specific
 *                  configuration;
 * params           the output parameter, into which everything is saved.-
 *
 */
extern void
init_coverage (FILE       *ini_file,
               char       *tx_sections_list,
               Parameters *params);


/**
 * Calculates the coverage prediction for one transmitter.
 *
 * params           a structure holding configuration parameters which are common
 *                  to all transmitters;
 * tx_params        the output parameter: a structure holding transmitter-specific
 *                  configuration parameters needed for calculation;
 * eric_params      contains the four tunning parameters for the Ericsson 9999 
 *                  model;
 * eric_params_len  the number of parameters within the received vector, four 
 *                  in this case (A0, A1, A2 and A3);
 *
 */
extern void
coverage (const Parameters     *params,
          const Tx_parameters  *tx_params,
          const double         *eric_params, 
          const unsigned int   eric_params_len);


/**
 * Starts a worker process.
 *
 * rank The rank of this worker process.-
 *
 */
extern void 
worker (const int rank,
        MPI_Comm *comm);


/**
 * Displays the calculation result in the standard output.
 *
 * params           a structure holding configuration parameters which are common
 *                  to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 */
extern void 
output_to_stdout (const Parameters *params,
                  const Tx_parameters *tx_params);


#endif
