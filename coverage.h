#ifndef _COVERAGE_CALCULATION_IN_GRASS_H_
#define _COVERAGE_CALCULATION_IN_GRASS_H_

#define _COVERAGE_MASTER_RANK_  0

#include <mpi.h>
#include <time.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <grass/gis.h>
#include <grass/glocale.h>

#include "ini.h"
#include "performance/metric.h"



//
// A structure to hold all the configuration parameters and
// in-run data for the coverage-calculation process
//
struct Parameters
{
    const char *dem_map;
    const char *clutter_map;
    const char *tx_name;
    int         beam_direction;
    int         electrical_tilt;
    int         mechanical_tilt;
    double      antenna_height_AGL;
    const char *antenna_diagram_dir;
    const char *antenna_diagram_file;
    double      tx_east_coord;
    double      tx_north_coord;
    double      tx_frequency;
    double      tx_power;
    double      radius;
    double      rx_height_AGL;
    
    int     nrows;
    int     ncols;
    double  map_east;
    double  map_west;
    double  map_north;
    double  map_south;
    double  map_ew_res;
    double  map_ns_res;
    double  total_tx_height;
    //
    // NULL value used in maps
    //
    float null_value;
    //
    // 2D matrix containing the digital-elevation-model raster
    //
    double **m_rast;
    //
    // 2D matrix containing the clutter information of the area
    //
    double **m_clut;
    //
    // 2D area matrix where the path-loss predictions are saved
    //
    double **m_loss;
    //
    // 2D area matrix where the field measurements are saved
    //
    double **m_field_meas;
} __attribute__((__packed__));

typedef struct Parameters Parameters;

/**
 * Handles name=value pairs read from the INI file.
 * The [General] section is always parsed, whereas other sections are
 * marked with the 'tx_section' parameter. This is useful for having
 * one INI file with all the transmitters in it, selecting which one
 * to use with a command line parameter.
 *
 * user_struct  The structure in which the read parameters are kept;
 * section      the current section being read from the INI file;
 * name         the name part of the current line in the INI file;
 * value        the value part of the current line in the INI file;
 * tx_section   the name of the transmitter's section to be parsed;
 *              all other transmitter sections are ignored.-
 *
 */
static int params_handler (void *user_struct, 
                           const char *section, 
                           const char *name,
                           const char *value,
                           const char *tx_section)
{
    Parameters *pconfig = (Parameters *) user_struct;

    #define MATCH(s,n) strcasecmp(section, s) == 0 && strcasecmp(name, n) == 0

    if (MATCH ("general", "DEMMapName"))
        pconfig->dem_map = strdup (value);
    else if (MATCH ("general", "clutterMapName"))
        pconfig->clutter_map = strdup (value);
    else if (MATCH ("general", "receiverHeightAGL"))
        pconfig->rx_height_AGL = atof (value);
    else if (MATCH ("general", "frequency"))
        pconfig->tx_frequency = atof (value);
    else if (MATCH ("general", "radius"))
        pconfig->radius = atof (value);
    else if (MATCH ("general", "antennaDirectory"))
        pconfig->antenna_diagram_dir = strdup (value);
    else if (MATCH (tx_section, "cellName"))
        pconfig->tx_name = strdup (value);
    else if (MATCH (tx_section, "beamDirection"))
        pconfig->beam_direction = atoi (value);
    else if (MATCH (tx_section, "electricalTiltAngle"))
        pconfig->electrical_tilt = atoi (value);
    else if (MATCH (tx_section, "mechanicalTiltAngle"))
        pconfig->mechanical_tilt = atoi (value);
    else if (MATCH (tx_section, "heightAGL"))
        pconfig->antenna_height_AGL = atof (value);
    else if (MATCH (tx_section, "antennaFile"))
        pconfig->antenna_diagram_file = strdup (value);
    else if (MATCH (tx_section, "positionEast"))
        pconfig->tx_east_coord = atof (value);
    else if (MATCH (tx_section, "positionNorth"))
        pconfig->tx_north_coord = atof (value);
    else if (MATCH (tx_section, "power"))
        pconfig->tx_power = atof (value);
    else 
        return 0;  /* unknown section/name, error */
    return 1;
}



/**
 * Calculates the area coverage in a serial fashion.
 *
 * eric_params      Contains the four tunning parameters for the 
 *                  Ericsson 9999 model, set by the optimization 
 *                  algorithm;
 * eric_params_len  the number of parameters within the received vector,
 *                  four in this case (A0, A1, A2 and A3);
 * output_raster    the name of the output raster created;
 *                  no output is generated if this parameter is NULL.-
 *
 */
extern void
coverage_serial (const Parameters *params,
                 const double *eric_params, 
                 const unsigned int eric_params_len);

/**
 * Calculate coverage over the received subarea.
 *
 * params   A structure holding all parameters needed for 
 *          calculation;
 * rank     this processes' rank;
 * nrows    number of rows inside the subarea this process has to work on;
 * ncols    number of columns in each of the rows of the subarea.-
 *
 */
extern void process_subarea (const Parameters *params,
                             const int rank,
                             MPI_Comm *comm,
                             const int nrows,
                             const int ncols);

/**
 * Starts a worker process.
 *
 * rank The rank of this worker process.-
 *
 */
extern void 
worker (const int rank,
        MPI_Comm *comm);

#endif
