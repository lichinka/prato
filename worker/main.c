#include "worker/coverage.h"




/**
 * Calculates the area coverage using MPI to achieve parallelization.
 * This is the worker part of the algorithm. The worker processes are 
 * dinamically started from the master process.-
 *
 */
int main (int argc, char *argv[])
{
    int rank, size;
 
    MPI_Comm parent;
    MPI_Init (&argc, &argv);
 
    MPI_Comm_get_parent (&parent);
    if (parent == MPI_COMM_NULL)
    {
        fprintf (stderr, "ERROR No master process!\n");
        return -1;
    }
 
    MPI_Comm_remote_size (parent, &size); 
    if (size != 1)
    {
        fprintf (stderr, "ERROR There should be only one master process\n");
        return -1;
    }
  
    MPI_Comm_rank (parent, &rank);
    worker (rank, &parent);

    MPI_Finalize ( );
    return EXIT_SUCCESS;
} 

