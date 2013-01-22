#ifndef _COVERAGE_CALCULATION_IN_GRASS_H_
#define _COVERAGE_CALCULATION_IN_GRASS_H_

#define _CHAR_BUFFER_SIZE_          1024
#define _COVERAGE_MASTER_RANK_      0
#define _WORKER_IS_IDLE_TAG_        100
#define _WORKER_KEEP_WORKING_TAG_   105
#define _WORKER_SHUTDOWN_TAG_       110
#define _WORKER_OPTIMIZE_TAG_       115

//
// whether the target point is whithin the main antenna beam or not
//
#define _RADIO_ZONE_MAIN_BEAM_ON_       0x01
#define _RADIO_ZONE_MAIN_BEAM_OFF_      0xfe
//
// whether the target point is within the distance limit of the prediction model
//
#define _RADIO_ZONE_SECONDARY_BEAM_ON_  0x02
#define _RADIO_ZONE_SECONDARY_BEAM_OFF_ 0xfd
//
// whether the target point is within the distance limit of the prediction model
//
#define _RADIO_ZONE_MODEL_DISTANCE_ON_  0x04
#define _RADIO_ZONE_MODEL_DISTANCE_OFF_ 0xfb

//
// dimensions of the search vector - only used for optimization
//
#define _SEARCH_VECTOR_DIMENSIONS_  4



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
    // antenna azimuth
    int     beam_direction;
    int     electrical_tilt;
    int     mechanical_tilt;

    // antenna height above ground level
    double  antenna_height_AGL;

    // antenna height above sea level
    double  total_tx_height;

    // geographical coordinates of this transmitter
    double  tx_east_coord;
    double  tx_north_coord;

    // 2D RADIUS area matrix indeces of the geographical coordinates
    int tx_east_coord_idx;
    int tx_north_coord_idx;

    // transmit pilot power 
    double  tx_power;

    // name of this transmitter
    char    tx_name [_CHAR_BUFFER_SIZE_];

    // full path to the antenna diagram file
    char    antenna_diagram_file [_CHAR_BUFFER_SIZE_];

    // binary file from which the field measurements are read
    char    field_meas_map [_CHAR_BUFFER_SIZE_];

    // number of rows and columns of each of the 2D RADIUS area matrices
    int nrows;
    int ncols;

    // geographical limits of each of the 2D RADIUS area matrices
    double map_north;
    double map_east;
    double map_south;
    double map_west;

    // 2D RADIUS area matrix indeces of the geographical limits
    int map_north_idx;
    int map_east_idx;
    int map_south_idx;
    int map_west_idx;

    // 2D RADIUS area matrix containing the digital-elevation-model raster
    // as received from the master process
    double **m_dem;

    // 2D RADIUS area matrix containing the clutter information 
    // as received from the master process
    double **m_clut;
    
    // 2D RADIUS area matrix where the field measurements are kept
    // as received from the master process
    double **m_field_meas;

    // 2D RADIUS area matrix where the path-loss predictions are saved
    double **m_loss;

    // 2D RADIUS area matrix where the antenna losses are saved
    double **m_antenna_loss;

    // 2D RADIUS area matrix where radio zones are marked; 
    // each element is a bit mask of _RADIO_ZONE_* constants
    char **m_radio_zone;
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
    {
        strncpy (pconfig->tx_name,
                 value,
                 _CHAR_BUFFER_SIZE_);
    }
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
    {
        strncpy (pconfig->antenna_diagram_file, 
                 value, 
                 _CHAR_BUFFER_SIZE_);
    }
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
    else if (MATCH (tx_section, "measurementsMap"))
    {
        strncpy (pconfig->field_meas_map, 
                 value, 
                 _CHAR_BUFFER_SIZE_);
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

    // horizontal and vertical losses that define the main antenna beam
    int     main_zone_horiz;    
    int     main_zone_vert;

    // horizontal and vertical losses that define the secondary antenna beam
    int     sec_zone_horiz;
    int     sec_zone_vert;

    // a buffer to save the INI file to memory
    char *  ini_file_content;
    int     ini_file_content_size;

    // a flag to indicate the framework starts in optimization mode
    char    use_opt;

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
    else if (MATCH ("common", "firstRadioZoneHorizontal"))
        pconfig->main_zone_horiz = atoi (value);
    else if (MATCH ("common", "firstRadioZoneVertical"))
        pconfig->main_zone_vert = atoi (value);
    else if (MATCH ("common", "secondRadioZoneHorizontal"))
        pconfig->sec_zone_horiz = atoi (value);
    else if (MATCH ("common", "secondRadioZoneVertical"))
        pconfig->sec_zone_vert = atoi (value);
    else
        return 0;  /* unknown section/name, error */
    return 1;
}



/**
 *
 * Starts a worker process.
 *
 * rank The rank of this worker process;
 * comm the MPI communicator to use.-
 *
 */
extern void 
worker (const int rank,
        MPI_Comm comm);


/**
 * Calculates the coverage prediction for one transmitter, using the 
 * Ericsson 9999 model.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 *                  configuration parameters needed for calculation;
 * eric_params      contains the four tunning parameters for the Ericsson 9999
 *                  model;
 * eric_params_len  the number of parameters within the received vector, four 
 *                  in this case (A0, A1, A2 and A3);
 *
 */
void 
coverage (const Parameters     *params,
          const Tx_parameters  *tx_params,
          const double         *eric_params, 
          const unsigned int   eric_params_len);


/**
 * Displays the calculation result in the standard output.
 *
 * params           a structure holding configuration parameters which are common
 *                  to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 */
void 
output_to_stdout (const Parameters *params,
                  const Tx_parameters *tx_params);


#endif
