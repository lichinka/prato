#ifndef _COVERAGE_CALCULATION_IN_GRASS_H_
#define _COVERAGE_CALCULATION_IN_GRASS_H_

#include <mpi.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <pthread.h>
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
#include "worker/constants.h"



//
// A structure to hold the antenna diagram and
// runtime data for the antenna-calculation process
//
struct Diagram
{
    //
    // Horizontal antenna diagram
    //
    double *horizontal;
    //
    // Vertical antenna diagram
    //
    double *vertical;
    //
    // Antenna gain
    //
    double  gain;
} __attribute__((__packed__));

typedef struct Diagram Diagram;


//
// A structure to hold transmitter-specific configuration and runtime parameters
//
struct Tx_parameters
{
    // antenna azimuth
    int     beam_direction;
    
    // antenna tilt
    int     electrical_tilt;
    int     mechanical_tilt;

    // antenna diagram
    Diagram *diagram;

    // antenna height above ground level
    double  antenna_height_AGL;

    // antenna height above sea level
    double  total_tx_height;

    // geographical coordinates of this transmitter
    double  tx_east_coord;
    double  tx_north_coord;

    // 2D RADIUS area matrix indices of the geographical coordinates
    int     tx_east_coord_idx;
    int     tx_north_coord_idx;

    // transmit pilot power 
    double  tx_power;

    // name of this transmitter
    char    tx_name [_CHAR_BUFFER_SIZE_];

    // full path to the antenna diagram file
    char    antenna_diagram_file [_CHAR_BUFFER_SIZE_];

    // raster map file from which the field measurements are read
    char    field_meas_map [_CHAR_BUFFER_SIZE_];

    // number of rows and columns of each of the 2D RADIUS area matrices
    int     nrows;
    int     ncols;

    // geographical limits of each of the 2D RADIUS area matrices
    double map_north;
    double map_east;
    double map_south;
    double map_west;

    // the four tunning parameters for the E/// model (A0, A1, A2, A3)
    double eric_params [4];

    // 2D RADIUS area matrix indices of the geographical limits
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

    // the number of valid field measurements,
    // used for calculating the mean square error
    unsigned int field_meas_count;

    // 2D RADIUS area matrix where the path-loss predictions are saved
    double **m_loss;

    // 2D RADIUS area matrix where the antenna losses are saved
    double **m_antenna_loss;

    // 2D RADIUS area matrix where radio zones are marked; 
    // each element is a bit mask of _RADIO_ZONE_* constants
    char **m_radio_zone;

    // 2D RADIUS area matrix where the obstacle heights are saved
    // (used for line-of-sight calculation)
    double **m_obst_height;

    // 2D RADIUS area matrix where the obstacle distances are saved
    // (used for line-of-sight calculation)
    double **m_obst_dist;
       
    // 2D RADIUS area matrix for obstacle distance offsets 
    // (used for line-of-sight calculation)
    double **m_obst_offset;

    //////
    // GPU-specific parameters follow
    //
    
    // a structure holding OpenCL specific data (platform, device, kernel, ...) 
    // as defined in the OCL_Common library
    OCL_object *ocl_obj;

    // GPU mapping of the various 2D RADIUS area matrices defined above
    cl_mem      *m_dem_dev;
    cl_mem      *m_clut_dev;
    cl_mem      *m_field_meas_dev;
    cl_mem      *m_loss_dev;
    cl_mem      *m_antenna_loss_dev;
    cl_mem      *m_radio_zone_dev;

    // GPU matrix holding the heights of the obstacles around the transmitter
    // (only used for the line-of-sight calculation)
    cl_mem      *m_obst_height_dev;

    // GPU matrix holding the distances to the obstacles around the transmitter
    // (only used for the line-of-sight calculation)
    cl_mem      *m_obst_dist_dev;

    // GPU vector where the partial sum by row is saved
    // (only used for the objective-funcion calculation)
    double      *v_partial_sum;
    cl_mem      *v_partial_sum_dev;

    // GPU vector where the clutter categories <-> loss mapping is kept
    cl_mem      *v_clutter_loss_dev;
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
static int
tx_params_handler (void *user_struct, 
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
    else if (MATCH (tx_section, "a0"))
    {
        pconfig->eric_params[0] = atof (value);
    }
    else if (MATCH (tx_section, "a1"))
    {
        pconfig->eric_params[1] = atof (value);
    }
    else if (MATCH (tx_section, "a2"))
    {
        pconfig->eric_params[2] = atof (value);
    }
    else if (MATCH (tx_section, "a3"))
    {
        pconfig->eric_params[3] = atof (value);
    }
    else
        return 0;  /* unknown section/name, error */
    return 1;
}


//
// A structure to hold all common configuration parameters and
// in-run data for the coverage-calculation process
//
struct Parameters
{
    // number of transmitters being processed
    int ntx;

    // parameters of the transmitters being processed
    Tx_parameters *tx_params;

    // NULL value used in DEM and measurement maps
    float fcell_null_value;

    // NULL value used in clutter maps
    int cell_null_value;

    // 2D matrix containing the digital-elevation-model raster
    double **m_dem;

    // 2D matrix containing the clutter information of the area
    double **m_clut;
    
    // 2D area matrix where the path-loss predictions are saved
    double **m_loss;
    
    // 2D area matrix where the field measurements are kept
    // before dispatching them to the workers, e.g. a buffer
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
    
    // the four tunning parameters for the E/// model (A0, A1, A2, A3)
    double eric_params [4];

    // category-based clutter losses are saved here;
    // each category is used as an index within this array to retrieve 
    // the correct loss; the initial values are read from the INI file
    int     clutter_category_count;
    double  clutter_loss        [_CHAR_BUFFER_SIZE_];

    // horizontal and vertical losses that define the main antenna beam
    int     main_zone_horiz;    
    int     main_zone_vert;

    // horizontal and vertical losses that define the secondary antenna beam
    int     sec_zone_horiz;
    int     sec_zone_vert;

    // a buffer to save the INI file to memory
    char *  ini_file_content;
    int     ini_file_content_size;

    // the number of generations to run the optimization algorithm within 
    // the framework, i.e. per-worker optimization mode 
    // (a value of 0 turns this off)
    int     use_opt;

    // the number of generations to run the optimization algorithm, but
    // using the workers just to calculate the objective function value;
    // this happens in every iteration of the optimization algorithm,
    // which runs on the master process 
    // (a value of 0 turns this off)
    int    use_master_opt;

    // a flag to indicate the GPU should be used 
    // on the worker side, if available
    char    use_gpu;
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
static int
common_params_handler (void *user_struct, 
                       const char *section, 
                       const char *name,
                       const char *value,
                       const char *tx_section)
{
    Parameters *pconfig = (Parameters *) user_struct;

    #define MATCH(s,n) strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0

    if (MATCH ("common", "DEMMapName"))
        strncpy (pconfig->dem_map,
                 value,
                 _CHAR_BUFFER_SIZE_);
    else if (MATCH ("common", "clutterMapName"))
        strncpy (pconfig->clutter_map,
                 value,
                 _CHAR_BUFFER_SIZE_);
    else if (MATCH ("common", "receiverHeightAGL"))
        pconfig->rx_height_AGL = atof (value);
    else if (MATCH ("common", "frequency"))
        pconfig->frequency = atof (value);
    else if (MATCH ("common", "radius"))
        pconfig->radius = atof (value);
    else if (MATCH ("common", "antennaDirectory"))
        strncpy (pconfig->antenna_diagram_dir, 
                 value, 
                 _CHAR_BUFFER_SIZE_);
    else if (MATCH ("common", "firstRadioZoneHorizontal"))
        pconfig->main_zone_horiz = atoi (value);
    else if (MATCH ("common", "firstRadioZoneVertical"))
        pconfig->main_zone_vert = atoi (value);
    else if (MATCH ("common", "secondRadioZoneHorizontal"))
        pconfig->sec_zone_horiz = atoi (value);
    else if (MATCH ("common", "secondRadioZoneVertical"))
        pconfig->sec_zone_vert = atoi (value);
    else if (MATCH ("common", "clutterCategoryCount"))
        pconfig->clutter_category_count = atoi (value);
    else if (MATCH ("common", "clutterLoss0"))
        pconfig->clutter_loss[0] = atof (value);
    else if (MATCH ("common", "clutterLoss1"))
        pconfig->clutter_loss[1] = atof (value);
    else if (MATCH ("common", "clutterLoss2"))
        pconfig->clutter_loss[2] = atof (value);
    else if (MATCH ("common", "clutterLoss3"))
        pconfig->clutter_loss[3] = atof (value);
    else if (MATCH ("common", "clutterLoss4"))
        pconfig->clutter_loss[4] = atof (value);
    else if (MATCH ("common", "clutterLoss5"))
        pconfig->clutter_loss[5] = atof (value);
    else if (MATCH ("common", "clutterLoss6"))
        pconfig->clutter_loss[6] = atof (value);
    else if (MATCH ("common", "clutterLoss7"))
        pconfig->clutter_loss[7] = atof (value);
    else if (MATCH ("common", "clutterLoss8"))
        pconfig->clutter_loss[8] = atof (value);
    else if (MATCH ("common", "clutterLoss9"))
        pconfig->clutter_loss[9] = atof (value);
    else if (MATCH ("common", "clutterLoss10"))
        pconfig->clutter_loss[10] = atof (value);
    else if (MATCH ("common", "clutterLoss11"))
        pconfig->clutter_loss[11] = atof (value);
    else if (MATCH ("common", "clutterLoss12"))
        pconfig->clutter_loss[12] = atof (value);
    else if (MATCH ("common", "clutterLoss13"))
        pconfig->clutter_loss[13] = atof (value);
    else if (MATCH ("common", "clutterLoss14"))
        pconfig->clutter_loss[14] = atof (value);
    else if (MATCH ("common", "clutterLoss15"))
        pconfig->clutter_loss[15] = atof (value);
    else
        return 0;       /* unknown section/name, error */
    return 1;
}



/**
 * Allocates a 2D matrix of the specified dimensions as a continous chunk 
 * of memory. Despite this, the matrix can be referenced with two indices
 * as `matrix[i][j]`.
 * This function returns the target pointer to where the memory has been 
 * allocated, i.e. `m_ptr`.
 *
 * nrows    the number of rows of the matrix;
 * ncols    the number of columns in the matrix;
 * m_ptr    target pointer, where the address of the allocated memory
 *          is saved (output parameter).-
 *
 */
static double **
prato_alloc_double_matrix (const int    nrows,
                           const int    ncols,
                           double     **m_ptr)
{
    int r;

    //
    // only allocate new memory if the target pointer is NULL
    //
    if (m_ptr == NULL)
    {
        double *m_ptr_data = (double *) calloc (nrows * ncols, 
                                                sizeof (double));
        m_ptr = (double **) calloc (nrows,
                                    sizeof (double *));
        for (r = 0; r < nrows; r ++)
            m_ptr[r] = &(m_ptr_data[r * ncols]);
    }
    return m_ptr;
}



/**
 * Allocates a 2D matrix of the specified dimensions as a continous chunk 
 * of memory. Despite this, the matrix can be referenced with two indices
 * as `matrix[i][j]`.
 * This function returns the target pointer to where the memory has been 
 * allocated, i.e. `m_ptr`.
 *
 * nrows    the number of rows of the matrix;
 * ncols    the number of columns in the matrix;
 * m_ptr    target pointer, where the address of the allocated memory
 *          is saved (output parameter).-
 *
 */
static char **
prato_alloc_char_matrix (const int    nrows,
                         const int    ncols,
                         char       **m_ptr)
{
    int r;

    //
    // only allocate new memory if the target pointer is NULL
    //
    if (m_ptr == NULL)
    {
        char *m_ptr_data = (char *) calloc (nrows * ncols, 
                                            sizeof (char));
        m_ptr = (char **) calloc (nrows,
                                  sizeof (char *));
        for (r = 0; r < nrows; r ++)
            m_ptr[r] = &(m_ptr_data[r * ncols]);
    }
    return m_ptr;
}



/**
 * Calculates the coverage prediction for one transmitter, using the E/// model.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 * rank             rank of the worker process.-
 *
 */
void 
coverage (Parameters    *params,
          Tx_parameters *tx_params,
          const int rank);


/**
 * Displays the calculation result in the standard output.
 *
 * params_ptr   a pointer to the structure holding configuration parameters 
 *              which are common to all transmitters;
 */
void *
output_to_stdout (void *params_ptr);


#endif
