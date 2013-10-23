#include "worker/optimize.h"
#include "worker/antenna.h"
#include "worker/eric.h"




/**
 * Objective function evaluation, used by the optimization algorithm.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration 
 *              parameters;
 * radio_zone   radio zone for which the objective function is calculated;
 * sol_vector   solution vector over which the objective function is calculated;
 * comm         the object used to communicate with the workers;
 *
 */
static double
obj_func (Parameters    *params,
          Tx_parameters *tx_params,
          const char     radio_zone,
          double        *sol_vector,
          MPI_Comm      *comm)
{
    double score [params->ntx][2];
    double ret_value [2];
    MPI_Status status;

#ifdef _PERFORMANCE_METRICS_
    measure_time ("Send solution to all workers");
#endif
    //
    // broadcast the new solution to all workers
    //
    MPI_Bcast (sol_vector,
               params->clutter_category_count,
               MPI_DOUBLE,
               _COVERAGE_MASTER_RANK_,
               *comm);
#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
    measure_time ("Gather partial objective-function values");
#endif
    //
    // receive the partial objective-function values from all workers,
    // aggregating the partial values before calculating the total 
    //
    ret_value[0] = 0;
    ret_value[1] = 0;
    int workers_evaluating = 0;
    while (workers_evaluating < params->ntx)
    {
        //
        // receive the partial objective-function value from this worker
        //
        MPI_Recv (&(score[workers_evaluating][0]),
                  2,
                  MPI_DOUBLE,
                  MPI_ANY_SOURCE,
                  MPI_ANY_TAG,
                  *comm,
                  &status);
        if (status.MPI_ERROR)
        {
            int worker_rank = status.MPI_SOURCE;
            fprintf (stderr, 
                     "*** ERROR: Objective-function value incorrectly received from %d. worker\n",
                     worker_rank);
            fflush (stderr);
            exit (1);
        }
        //
        // aggregate the received squared error
        //
        ret_value[0] += score[workers_evaluating][0];
        
        //
        // aggregate the received field-measurement count
        //
        ret_value[1] += score[workers_evaluating][1];

        workers_evaluating ++;
    }
#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
    measure_time ("Build complete objective-function value");
#endif
    //
    // the total mean-squared error
    //
    ret_value[0] /= ret_value[1];

#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
#endif
    return ret_value[0];
}



/**
 * The differential evolution optimization algorithm.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration 
 *              parameters;
 * D            dimension of the search space;
 * NP           population size, try with 20 * `D`;
 * Gmax         number of generations to evolve the solution, try with 1000;
 * CR           cross-over probability factor [0,1], try with 0.9;
 * F            differential weight factor [0, 2], try with 0.9;
 * reset_seed   whether to reset the random seed between runs or not;
 * radio_zone   radio zone under optimization, used for objective calculation;
 * X_low        lower bound of the search vector, should be `D` long;
 * X_up         upper bound of the search vector, should be `D` long;
 * comm         the object used to communicate with the workers;
 *
 */
static void 
de (Parameters     *params,
    Tx_parameters  *tx_params,
    const int       D, 
    const int       NP, 
    const int       Gmax, 
    const double    CR, 
    const double    F,
    const char      reset_seed,
    const char      radio_zone,
    const double   *X_low,
    const double   *X_up,
    MPI_Comm       *comm)
{
    register int i, j, k, r1, r2, r3, jrand, numofFE = 0;
    int index = -1;
    double **popul, **next, **ptr, *iptr, *U, min_value = DBL_MAX, totaltime = 0.0;
    clock_t starttime, endtime;

    //
    // reset the random seed?
    //
    if (reset_seed)
        srand (time (0));
    else
        srand (0);

    // Printing out information about optimization process for the user	

    printf ("*** INFO: DE parameters\t");
    printf ("NP = %d, Gmax = %d, CR = %.2f, F = %.2f\n",
            NP, Gmax, CR, F);

    printf ("*** INFO: Optimization problem dimension is %d.\n", D);

    /* Starting timer    */
    starttime = clock();


    /* Allocating memory for current and next populations, intializing
      current population with uniformly distributed random values and
      calculating value for the objective function	*/

    popul = (double **)malloc(NP*sizeof(double *));
    if (popul == NULL) perror("malloc");

    next = (double **)malloc(NP*sizeof(double *));
    if (next == NULL) perror("malloc");

    for (i=0; i < NP; i++)
    {
        popul[i] = (double *)malloc((D+1)*sizeof(double));
        if (popul[i] == NULL) perror("malloc");

        for (j=0; j < D; j++)
            popul[i][j] = X_low[j] + (X_up[j] - X_low[j])*URAND;

        popul[i][D] = obj_func (params,
                                tx_params,
                                radio_zone,
                                popul[i],
                                comm);
        numofFE++;

        next[i] = (double *)malloc((D+1)*sizeof(double));
        if (next[i] == NULL) perror("malloc");
    }

    /* Allocating memory for a trial vector U	*/

    U = (double *)malloc((D+1)*sizeof(double));
    if (U == NULL) perror("malloc");

    /* The main loop of the algorithm	*/

    for (k=0; k < Gmax; k++)
    {

      for (i=0; i < NP; i++)	/* Going through whole population	*/
      {

         /* Selecting random indeces r1, r2, and r3 to individuals of
            the population such that i != r1 != r2 != r3	*/

         do
         {
            r1 = (int)(NP*URAND);
         } while( r1==i );
         do
         {
            r2 = (int)(NP*URAND);
         } while( r2==i || r2==r1);
         do
         {
            r3 = (int)(NP*URAND);
         } while( r3==i || r3==r1 || r3==r2 );

         jrand = (int)(D*URAND);

         /* Mutation and crossover	*/

         for (j=0; j < D; j++)
         {
            if (URAND < CR || j == jrand)
            {
               U[j] = popul[r3][j] + F*(popul[r1][j] - popul[r2][j]);
            }
            else
               U[j] = popul[i][j];
            //
            // check that all components of the trial vector
            // are within the allowed values
            //
            if (U[j] < X_low[j])
                U[j] = X_low[j];
            if (U[j] > X_up[j])
                U[j] = X_up[j];
         }

         U[D] = obj_func (params,
                          tx_params,
                          radio_zone,
                          U,
                          comm);
         numofFE++;

         printf ("Generation %d/%d\tscore %20.10f\n", k, Gmax, U[D]);

         /* Comparing the trial vector 'U' and the old individual
            'next[i]' and selecting better one to continue in the
            next population.	*/

         if (U[D] <= popul[i][D])
         {
            iptr = U;
            U = next[i];
            next[i] = iptr;
         }
         else
         {
            for (j=0; j <= D; j++)
                next[i][j] = popul[i][j];
         }

      }	/* End of the going through whole population	*/


      /* Pointers of old and new populations are swapped	*/

      ptr = popul;
      popul = next;
      next = ptr;

    }	/* End of the main loop		*/


    /* Stopping timer	*/

    endtime = clock();
    totaltime = (double)(endtime - starttime);

    //
    // output the final population
    //
    for (i=0; i < NP; i++)
    {
        for (j=0; j <= D; j++)
            fprintf (stdout, "%.15f ", popul[i][j]);
        fprintf (stdout, "\n");
    }

    //
    // finding the best individual
    //
    for (i=0; i < NP; i++)
    {
        if (popul[i][D] < min_value)
        {
            min_value = popul[i][D];
            index = i;
        }
    }

    //
    // print out information about optimization process for the user
    //
    printf ("Execution time: %.3f s\n", totaltime / (double)CLOCKS_PER_SEC);
    printf ("Number of objective function evaluations: %d\n", numofFE);

    printf ("[Solution]\n");
    for (i=0; i < D; i++)
      printf("p%d = %.15f\n", i, popul[index][i]);

    printf ("\nObjective function value: ");
    printf ("%.15f\n", popul[index][D]);

    //
    // copy the solution values to the internal structure,
    // used for coverage calculation
    //
    for (i = 0; i < D; i ++)
        params->clutter_loss[i] = popul[index][i];

    /* Freeing dynamically allocated memory	*/ 
    for (i=0; i < NP; i++)
    {
      free(popul[i]);
      free(next[i]);
    }
    free(popul);
    free(next);
    free(U);
}



/**
 * Optimizes the parameters of the prediction model to fit field measurements
 * within different radio zones using a deterministic (analytical) approach
 * for each of the transmitters.
 * It also optimizes the clutter-category losses using differential evolution
 * on the master process, evaluating the objective function on the workers.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration 
 *              parameters;
 * comm         the object used to communicate with the workers;
 *
 */
void 
optimize_on_master (Parameters    *params,
                    Tx_parameters *tx_params,
                    MPI_Comm      *comm)
{
    int i;

    //
    // define lower and upper bounds for each search-vector component,
    // i.e. solutions should be within these limits
    //
    double *search_low = (double *) calloc (params->clutter_category_count,
                                            sizeof (double));
    double *search_up  = (double *) calloc (params->clutter_category_count,
                                            sizeof (double));
    //
    // since we are looking for clutter losses, we define a range 0~40 dB
    //
    for (i = 0; i < params->clutter_category_count; i ++)
    {
        search_low[i] = 0;
        search_up[i]  = 40;
    }
#ifdef _DEBUG_INFO_
    //
    // DEBUG: start only one iteration of the optimization process
    //
    de (params,
        tx_params,
        params->clutter_category_count,
        params->clutter_category_count,
        1,
        0.9,
        0.9,
        0,
        _RADIO_ZONE_MAIN_BEAM_ON_,
        search_low,
        search_up,
        comm);
#else
    //
    // start optimization
    //
    de (params,
        tx_params,
        params->clutter_category_count,
        20 * params->clutter_category_count,
        params->use_master_opt,
        0.9,
        0.9,
        1,
        _RADIO_ZONE_MAIN_BEAM_ON_,
        search_low,
        search_up,
        comm);
#endif
    //
    // tell the workers to shutdown
    //
    double stop_vector [params->clutter_category_count];

    for (i = 0; i < params->clutter_category_count ; i ++)
        stop_vector[i] = _PI_;

    MPI_Bcast (stop_vector,
               params->clutter_category_count,
               MPI_DOUBLE,
               _COVERAGE_MASTER_RANK_,
               *comm);
    //
    // free reserved buffers
    //
    free (search_up);
    free (search_low);
}

