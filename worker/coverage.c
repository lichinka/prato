#include "coverage.h"
#include "performance/metric.h"




/**
 * Calculate coverage over the received subarea.
 *
 * params   A structure holding all parameters needed for 
 *          calculation;
 * rank     this processes' rank;
 * comm     the MPI communicator to use;
 * nrows    number of rows inside the subarea this process has to work on;
 * ncols    number of columns in each of the rows of the subarea.-
 *
 */
void process_subarea (const Parameters *params,
                      const int rank,
                      MPI_Comm *comm,
                      const int nrows,
                      const int ncols)
{
    int i, j;
    MPI_Datatype subarea;

    //
    // define a datatype for the sub-area each worker receives
    //
    MPI_Type_contiguous (nrows * ncols,
                         MPI_DOUBLE, 
                         &subarea);
    MPI_Type_commit (&subarea);

    //
    // allocate memory for the data received
    //
    double *chunk_data = (double *) calloc (nrows * ncols,
                                            sizeof (double));
    double **chunk = (double **) calloc (nrows,
                                         sizeof (double *));
    for (i = 0; i < nrows; i ++)
        chunk[i] = &(chunk_data[i * ncols]);

    //
    // scatter data around
    //
    MPI_Scatter (&(params->m_rast[0][0]), 
                 1, 
                 subarea, 
                 &(chunk[0][0]), 
                 1, 
                 subarea, 
                 _COVERAGE_MASTER_RANK_, 
                 *comm);

    //
    // print the received data out
    //
    for(i = 0; i < nrows; i++)
        for (j = 0; j < ncols; j++)
            printf("r_%d\t%lf\n", rank, chunk[i][j]);

    MPI_Type_free (&subarea);
    free (chunk);
    free (chunk_data);
}




/**
 * Starts a worker process.
 *
 * rank The rank of this worker process;
 * comm the MPI communicator to use.-
 *
 */
void worker (const int rank,
             MPI_Comm *comm)
{
    int nrows, ncols;

    //
    // synchronize all processes before starting
    //
    MPI_Barrier (*comm);

    //
    // received the row size and the number of rows per worker
    //
    MPI_Bcast (&ncols, 
               1, 
               MPI_INT, 
               _COVERAGE_MASTER_RANK_,
               *comm);

    MPI_Bcast (&nrows, 
               1, 
               MPI_INT, 
               _COVERAGE_MASTER_RANK_,
               *comm);
    printf ("Worker %d received (%d, %d) via broadcast message\n", rank,
                                                                   ncols,
                                                                   nrows);
    /*
    // start processing the subarea
    //
    process_subarea (NULL,
                     rank,
                     comm,
                     nrows, 
                     ncols);*/
}

