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
                             double **m_antenna_loss);

#endif
