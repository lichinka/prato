#ifndef _MASTER_PROCESS_FOR_COVERAGE_CALCULATION_IN_GRASS_H_
#define _MASTER_PROCESS_FOR_COVERAGE_CALCULATION_IN_GRASS_H_



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
 * Starts coverage calculation over MPI.
 *
 * argc             Number of command line parameters;
 * argv             array containing command line parameters;
 * params           a structure holding all parameters needed for calculation.-
 *
 */
extern void
coverage_mpi (int argc, 
              char *argv [],
              Parameters *params);

#endif
