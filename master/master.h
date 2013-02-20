#ifndef _MASTER_PROCESS_FOR_COVERAGE_CALCULATION_IN_GRASS_H_
#define _MASTER_PROCESS_FOR_COVERAGE_CALCULATION_IN_GRASS_H_

#include "measurement.h"
#include "worker/worker.h"


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
void
init_coverage (FILE       *ini_file,
               char       *tx_sections_list,
               Parameters *params);


/**
 * Starts coverage calculation over MPI.
 *
 * argc             Number of command line parameters;
 * argv             array containing command line parameters;
 * params           a structure holding all parameters needed for calculation.-
 *
 */
void
coverage_mpi (int argc, 
              char *argv [],
              Parameters *params);

#endif
