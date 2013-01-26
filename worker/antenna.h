#ifndef _COVERAGE_ANTENNA_H_
#define _COVERAGE_ANTENNA_H_

#define _DIAGRAM_SIZE_  360

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "worker/coverage.h"



//
// A structure to hold all the antenna diagram and
// in-run data for the antenna-calculation process
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
 *                  parameters;
 *
 * WARNING: the output of this function overwrites 
 *          `tx_params->m_loss`, `tx_params->m_radio_zone` 
 *          and `tx_params->m_antenna_loss`
 *
 */
void
calculate_antenna_influence (Parameters    *params,
                             Tx_parameters *tx_params);

#endif
