#include "worker/coverage.h"



/**
 * Calculates the area coverage using MPI to achieve parallelization.
 * This is the worker part of the algorithm. The worker processes should
 * be concurrently started with `mpirun`.-
 *
 */
int main (int argc, char *argv[])
{
    int rank, size;
 
    MPI_Comm master_comm;
    MPI_Init (&argc, &argv);
 
    master_comm = MPI_COMM_WORLD;
    if (master_comm == MPI_COMM_NULL)
    {
        fprintf (stderr, "ERROR No master process\n");
        return -1;
    }
 
    MPI_Comm_size (master_comm, &size); 
    if (size < 1)
    {
        fprintf (stderr, "ERROR There are not enough worker processes\n");
        return -1;
    }
  
    MPI_Comm_rank (master_comm, &rank);
    worker (rank, master_comm);

    MPI_Finalize ( );
    return EXIT_SUCCESS;
} 

