#include "worker/optimize.h"
#include "worker/antenna.h"
#include "worker/eric.h"



/**
 * Objective function evaluation, executed on the GPU, if available.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration 
 *              parameters;
 * radio_zone   radio zone for which the objective function is calculated;
 *
 */
static double
obj_func_on_gpu (Parameters    *params,
                 Tx_parameters *tx_params,
                 const char    radio_zone)
{
    int r, c;
    double ret_value;

    //
    // check if the needed buffers are already on the device
    //
    if (tx_params->m_field_meas_dev == NULL)
    {
        //
        // not yet ... create a device buffer for the field measurements
        //
        size_t fm_buff_size = tx_params->nrows * 
                              tx_params->ncols * 
                              sizeof (tx_params->m_field_meas[0][0]);
        tx_params->m_field_meas_dev = (cl_mem *) malloc (sizeof (cl_mem));
        *(tx_params->m_field_meas_dev) = create_buffer (tx_params->ocl_obj,
                                                        CL_MEM_READ_ONLY, 
                                                        fm_buff_size);
        //
        // create a device buffer for the by-row partial sum 
        //
        size_t ps_buff_size = tx_params->nrows *
                              sizeof (tx_params->m_loss[0][0]);
        tx_params->v_partial_sum = (double *) malloc (ps_buff_size);
        tx_params->v_partial_sum_dev = (cl_mem *) malloc (sizeof (cl_mem));
        *(tx_params->v_partial_sum_dev) = create_buffer (tx_params->ocl_obj,
                                                         CL_MEM_READ_WRITE,
                                                         ps_buff_size);
        //
        // send the buffers data to the device
        //
        write_buffer (tx_params->ocl_obj,
                      0,
                      tx_params->v_partial_sum_dev,
                      ps_buff_size,
                      tx_params->v_partial_sum);
        write_buffer (tx_params->ocl_obj,
                      0,
                      tx_params->m_field_meas_dev,
                      fm_buff_size,
                      tx_params->m_field_meas[0]);
        //
        // count the number of valid field measurements
        // to correctly calculate the mean square error
        //
        tx_params->field_meas_count = 0;
        for (r = 0; r < tx_params->nrows; r ++)
        {
            for (c = 0; c < tx_params->ncols; c ++)
            {
                //
                // use only valid path-loss values
                //
                double pl = tx_params->m_loss[r][c];
                if (!isnan (pl))
                {
                    //
                    // look for the target radio zone ...
                    //
                    char rz = tx_params->m_radio_zone[r][c];
                    if ((rz & radio_zone) > 0)
                    {
                        //
                        // ... and a field measurement there
                        //
                        double fm = tx_params->m_field_meas[r][c];
                        if (!isnan (fm))
                            tx_params->field_meas_count ++;
                    }
                }
            }
        }
    }
    //
    // activate the kernel 
    // 
    activate_kernel (tx_params->ocl_obj,
                     "obj_func_kern");

    // set kernel parameters
    set_kernel_double_arg (tx_params->ocl_obj,
                           0,
                           &tx_params->tx_power);
    set_kernel_value_arg (tx_params->ocl_obj,
                          1,
                          sizeof (radio_zone),
                          &radio_zone);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        2,
                        tx_params->m_field_meas_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        3,
                        tx_params->m_radio_zone_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        4,
                        tx_params->m_loss_dev);
    set_kernel_mem_arg (tx_params->ocl_obj,
                        5,
                        tx_params->v_partial_sum_dev);
    //
    // define a 1D execution range for the kernel ...
    //
    size_t global_size = tx_params->nrows * tx_params->ncols,
           local_size  = tx_params->ncols;

    // reserve local memory on the device
    size_t lmem_size = local_size * 
                       sizeof (double);
    set_local_mem (tx_params->ocl_obj,
                   6,
                   lmem_size);
    //
    // ... and execute it
    //
    run_kernel_1D_blocking (tx_params->ocl_obj,
                            0,
                            NULL,
                            &global_size,
                            &local_size);
    //
    // bring the sum of each row from the GPU;
    // and calculate the average on the CPU
    //
    size_t ps_buff_size = tx_params->nrows * 
                          sizeof (tx_params->v_partial_sum[0]);
    read_buffer_blocking (tx_params->ocl_obj,
                          0,
                          tx_params->v_partial_sum_dev,
                          ps_buff_size,
                          tx_params->v_partial_sum);
    for (r = 0; r < tx_params->nrows; r ++)
        ret_value += tx_params->v_partial_sum[r];
    //
    // the mean square error
    //
    ret_value /= tx_params->field_meas_count;

    return ret_value;
}



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

    //
    // copy the solution values to the internal structure,
    // used for coverage calculation
    //
    for (r = 0; r < params->clutter_category_count; r ++)
        params->clutter_loss[r] = sol_vector[r];

#ifdef _DEBUG_INFO_
    //
    // a fixed solution for testing (score is 101.0050757267 for KPODVI1)
    //
    tx_params->eric_params[0] = 80.36872968;
    tx_params->eric_params[1] = -71.74091439;
    tx_params->eric_params[2] = -23.82273819;
    tx_params->eric_params[3] = 47.83436989;
#endif
    //
    // recalculate the isotrophic prediction using the received solution
    //
    if (params->use_gpu)
    {
#ifdef _PERFORMANCE_METRICS_
        measure_time ("E/// on GPU");
#endif
        eric_pathloss_on_gpu (params,
                              tx_params);
#ifdef _PERFORMANCE_METRICS_
        measure_time (NULL);
        measure_time ("Apply antenna losses on GPU");
#endif
        apply_antenna_influence_gpu (params,
                                     tx_params);
#ifdef _PERFORMANCE_METRICS_
        measure_time (NULL);
        measure_time ("Objective function on GPU");
#endif
        ret_value = obj_func_on_gpu (params,
                                     tx_params,
                                     radio_zone);
    }
    else
    {
#ifdef _PERFORMANCE_METRICS_
        measure_time ("E/// on CPU");
#endif
        eric_pathloss_on_cpu (params,
                              tx_params);
#ifdef _PERFORMANCE_METRICS_
        measure_time (NULL);
        measure_time ("Apply antenna losses and objective function on CPU");
#endif
        //
        // objective function calculation, for 
        // each point in the path loss matrix ...
        //
        for (r = 0; r < tx_params->nrows; r ++)
        {
            for (c = 0; c < tx_params->ncols; c ++)
            {
                //
                // use only valid path-loss values
                //
                double pl = tx_params->m_loss[r][c];
                if (!isnan (pl))
                {
                    //
                    // apply the antenna loss
                    //
                    pl += tx_params->m_antenna_loss[r][c];

                    //
                    // look for the target radio zone ...
                    //
                    char rz = tx_params->m_radio_zone[r][c];
                    if ((rz & radio_zone) > 0)
                    {
                        //
                        // ... and a field measurement there
                        //
                        double fm = tx_params->m_field_meas[r][c];
                        if (!isnan (fm))
                        {
                            count ++;
                            ret_value += (tx_params->tx_power - pl - fm) *
                                         (tx_params->tx_power - pl - fm);
                        }
                    }
                }
            }
        }
        //
        // mean square error
        //
        ret_value /= count;
    }
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
                          U);
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
    // since we are looking for clutter losses, we define a range 0~255 dB
    //
    for (i = 0; i < params->clutter_category_count; i ++)
    {
        search_low[i] = 0;
        search_up[i]  = 255;
    }
    //
    // calculate the coverage for the first time to initialize all needed structures
    //
    coverage (params,
              tx_params);
    //
    // if the first coverage calculation happened on the GPU,
    // we need to refresh the memory buffers on the host
    //
    if (params->use_gpu)
    {
        size_t dbl_buff_size = tx_params->nrows * 
                               tx_params->ncols * 
                               sizeof (tx_params->m_loss[0][0]);
        size_t ch_buff_size = tx_params->nrows * 
                              tx_params->ncols * 
                              sizeof (tx_params->m_radio_zone[0][0]);
        read_buffer (tx_params->ocl_obj,
                     0,
                     tx_params->m_loss_dev,
                     dbl_buff_size,
                     tx_params->m_loss[0]);
        read_buffer (tx_params->ocl_obj,
                     0,
                     tx_params->m_antenna_loss_dev,
                     dbl_buff_size,
                     tx_params->m_antenna_loss[0]);
        read_buffer_blocking (tx_params->ocl_obj,
                              0,
                              tx_params->m_radio_zone_dev,
                              ch_buff_size,
                              tx_params->m_radio_zone[0]);
    }
    //
    // calculate parameter approximation for E/// prediction model
    //
    parameter_fine_tuning (params,
                           tx_params);
    //
    // ... and its objective function value
    //
    double score = obj_func (params,
                             tx_params,
                             _RADIO_ZONE_MAIN_BEAM_ON_,
                             tx_params->eric_params);
    fprintf (stdout, 
             "*** INFO: optimal values for E/// have score %g\n",
             score);

#ifdef _DEBUG_INFO_
    //
    // DEBUG: run only one iteration of the optimization process
    //
    de (params,
        tx_params,
        params->clutter_category_count,
        params->clutter_category_count,
        1,
        0.9,
        0.9,
        1,
        _RADIO_ZONE_MAIN_BEAM_ON_,
        search_low,
        search_up);
#else
    //
    // start optimization
    //
    de (params,
        tx_params,
        params->clutter_category_count,
        20 * params->clutter_category_count,
        200,
        0.9,
        0.9,
        1,
        _RADIO_ZONE_MAIN_BEAM_ON_,
        search_low,
        search_up);
#endif
    //
    // free reserved buffer
    //
    free (search_up);
    free (search_low);
}

