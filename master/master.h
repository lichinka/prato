#ifndef _MASTER_PROCESS_FOR_COVERAGE_CALCULATION_IN_GRASS_H_
#define _MASTER_PROCESS_FOR_COVERAGE_CALCULATION_IN_GRASS_H_

/**
 * Dynamically spawns worker processes to calculate the area coverage using MPI.
 * It returns the number of spawned workers, and a reference to the communicator
 * used to talk to them in the 'worker_comm' parameter.
 *
 * argc             Number of command line parameters;
 * argv             array containing command line parameters;
 * params           a structure holding all parameters needed for calculation;
 * worker_comm      output parameter: the communicator used to talk to the 
 *                  spawned workers.-
 *
 *
int spawn_workers (int argc, 
                   char *argv [],
                   Parameters *params,
                   MPI_Comm *worker_comm);*/


/**
 * Starts coverage calculation over MPI.
 *
 * argc             Number of command line parameters;
 * argv             array containing command line parameters;
 * params           a structure holding all parameters needed for calculation.-
 *
 */
void coverage_mpi (int argc, 
                   char *argv [],
                   Parameters *params);

#endif
