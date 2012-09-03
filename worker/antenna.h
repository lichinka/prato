#ifndef _COVERAGE_ANTENNA_H_
#define _COVERAGE_ANTENNA_H_

#define _PI_            3.14159265358979
#define _DIAGRAM_SIZE_  360

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>


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
 * ant_dir      Directory containing antenna diagram file;
 * ant_file     file name of the antenna diagram;
 * east         transmitter's eastern coordinate;
 * north        transmitter's northern coordinate;
 * total_height transmitter's height above sea level;
 * beam_dir     direction of the antenna beam, in degrees;
 * mech_tilt    mechanical antenna tilt angle, in degress;
 * radius       calculation radius around the given transmitter;
 * rec_height   Rx antenna height above ground level;
 * nrows        number of rows within the 2D-matrices;
 * ncols        number of columns within the 2D-matrices;
 * west_ext     western raster extent;
 * north_ext    northern raster extent;
 * ew_res       east/west raster resolution;
 * ns_res       north/south raster resolution;
 * null_value   the value representing NULL on the output map;
 * dem          a 2D-matrix containing the digital elevation model data;
 * path_loss    a 2D-matrix containing the path-loss values for the isotrophic
 *              antenna;
 *              WARNING: the output of this function overwrites this 2D-matrix!
 *
 */
void calculate_antenna_influence (const char *ant_dir,
                                  const char *ant_file,
                                  const double east,
                                  const double north,
                                  const double total_height,
                                  const int beam_dir,
                                  const int mech_tilt,
                                  const double radius,
                                  const double rec_height,
                                  const int nrows,
                                  const int ncols,
                                  const double west_ext,
                                  const double north_ext,
                                  const double ew_res,
                                  const double ns_res,
                                  const float null_value,
                                  double **dem,
                                  double **path_loss);

#endif
