#include "worker/optimize.h"
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
 *
 */
static double
obj_func (Parameters    *params,
          Tx_parameters *tx_params,
          const char     radio_zone,
          const double  *sol_vector)
{
    int r, c, count = 0;
    double ret_value = 0;

    if (params->use_gpu)
    {
#ifdef _PERFORMANCE_METRICS_
        measure_time ("E/// on GPU");
#endif
        eric_pathloss_on_gpu (params,
                              tx_params,
                              sol_vector);

        /*
        // FIXME: move this part to a separate function inside `antenna.c`
        //

        //
        // activate the kernel, 
        // to sum the antenna loss to the isotrophic path-loss
        // 
        activate_kernel (tx_params->ocl_obj,
                         "vector_sum_kern");
        //
        // define a 2D execution range for the kernel ...
        //
        double radius_in_meters = params->radius * 1000;
        int radius_in_pixels    = (int) (radius_in_meters / params->map_ew_res);
        int diameter_in_pixels  = 2 * radius_in_pixels;
        size_t ntile = diameter_in_pixels / _WORK_ITEMS_PER_DIMENSION_;
        size_t global_sizes [] = {ntile * _WORK_ITEMS_PER_DIMENSION_,
                                  ntile * _WORK_ITEMS_PER_DIMENSION_};
        size_t local_sizes [] = {_WORK_ITEMS_PER_DIMENSION_,
                                 _WORK_ITEMS_PER_DIMENSION_};

        // set pointer kernel parameters
        set_kernel_mem_arg (tx_params->ocl_obj,
                            0,
                            tx_params->m_antenna_loss_dev);
        set_kernel_mem_arg (tx_params->ocl_obj,
                            1,
                            tx_params->m_loss_dev);
        // reserve local memory on the device
        lmem_size = _WORK_ITEMS_PER_DIMENSION_ *
                    _WORK_ITEMS_PER_DIMENSION_ *
                    sizeof (double);
        set_local_mem (tx_params->ocl_obj,
                       2,
                       lmem_size);
        //
        // execute the sum kernel
        //
        run_kernel_2D_blocking (tx_params->ocl_obj,
                                0,
                                NULL,
                                global_sizes,
                                local_sizes);
        //
        // sync memory
        //
        read_buffer_blocking (tx_params->ocl_obj,
                              0,
                              tx_params->m_loss_dev,
                              ant_buff_size,
                              tx_params->m_loss[0]);
                              */
    }
    else
    {
#ifdef _PERFORMANCE_METRICS_
        measure_time ("E/// on CPU");
#endif
        //
        // recalculate the isotrophic prediction using the received solution
        //
        struct StructEric IniEric = {tx_params->tx_east_coord_idx,
                                     tx_params->tx_north_coord_idx,
                                     tx_params->antenna_height_AGL, 
                                     params->rx_height_AGL,
                                     tx_params->ncols,
                                     tx_params->nrows,
                                     params->map_ew_res,
                                     params->frequency, 
                                     (double) sol_vector[0],
                                     (double) sol_vector[1],
                                     (double) sol_vector[2], 
                                     (double) sol_vector[3], 
                                     1, 
                                     params->radius};
        EricPathLossSub (tx_params->m_obst_height,
                         tx_params->m_obst_dist,
                         tx_params->m_obst_offset,
                         tx_params->m_dem, 
                         tx_params->m_clut, 
                         tx_params->m_loss, 
                         &IniEric);
    }

#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
    measure_time ("Objective function evaluation on CPU");
#endif
    //
    // for each point in the path loss matrix ...
    //
    for (r = 0; r < tx_params->nrows; r ++)
    {
        for (c = 0; c < tx_params->ncols; c ++)
        {
            //
            // ... apply the previously calculated antenna losses
            //
            tx_params->m_loss[r][c] += tx_params->m_antenna_loss[r][c];

            //
            // look for the target radio zone ...
            //
            char rz = tx_params->m_radio_zone[r][c];
            if ((rz & _RADIO_ZONE_MAIN_BEAM_ON_) > 0)
            {
                //
                // ... and a field measurement there
                //
                float fm = (float) tx_params->m_field_meas[r][c];
                if ((!isnan (fm)) && (fm != params->null_value))
                {
                    count ++;
                    ret_value += ((tx_params->tx_power - tx_params->m_loss[r][c]) - fm) *
                                 ((tx_params->tx_power - tx_params->m_loss[r][c]) - fm);
                }
            }
        }
    }
    //
    // mean square error
    //
    ret_value /= count;
#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
#endif
    return ret_value;
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
    const double   *X_up)
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

    /* Printing out information about optimization process for the user	*/

    printf ("*** INFO: DE parameters\t");
    printf ("NP = %d, Gmax = %d, CR = %.2f, F = %.2f\n",
            NP, Gmax, CR, F);

    printf("*** INFO: Optimization problem dimension is %d.\n", D);

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
                                popul[i]);
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

         /* Selecting random indeces r1, r2, and r3 to individuls of
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
                          U);
         numofFE++;

         //printf ("Generation %d/%d, score %20.10f\n", k, Gmax, U[D]);

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

      printf ("Generation %d/%d\n", k, Gmax);
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
            fprintf (stdout, "%.15e ", popul[i][j]);
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
 * within different radio zones.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration 
 *              parameters;
 *
 */
void 
optimize (Parameters    *params,
          Tx_parameters *tx_params)
{
    //
    // define lower and upper bounds for each search-vector component,
    // i.e. solutions should be within these limits
    //
    double search_low [_SEARCH_VECTOR_DIMENSIONS_] = { 0.0,  0.0, -24.0, 0.0};
    double search_up  [_SEARCH_VECTOR_DIMENSIONS_] = {76.0, 64.0,   0.0, 0.2};

    //
    // calculate the coverage for the first time to initialize all needed structures
    //
    coverage (params,
              tx_params,
              search_low,
              _SEARCH_VECTOR_DIMENSIONS_);
    //
    // start optimization
    //
    de (params,
        tx_params,
        _SEARCH_VECTOR_DIMENSIONS_,
        10 * _SEARCH_VECTOR_DIMENSIONS_,
        10,
        0.9,
        0.9,
        1,
        _RADIO_ZONE_MAIN_BEAM_ON_,
        search_low,
        search_up);
}

