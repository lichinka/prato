#ifndef _COVERAGE_ANTENNA_H_
#define _COVERAGE_ANTENNA_H_

#define _PI_            3.14159265358979
#define _DIAGRAM_SIZE_  360

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "worker/coverage.h"



//
// A structure to hold all the configuration parameters and
// in-run data for the coverage-calculation process
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

static Diagram *diagram = NULL;

/**
 * Calculates additional gain/pathloss according to the antenna's
 * 3-dimensional diagram.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters.-
 *
 * WARNING: the output of this function overwrites the path-loss matrix!
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
                             GPU_parameters *gpu_params,
                             double **m_dem,          
                             double **m_loss);

#endif
