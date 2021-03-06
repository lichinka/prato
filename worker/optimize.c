#include "worker/optimize.h"
#include "worker/antenna.h"
#include "worker/eric.h"




/**
 * Override the clutter losses in order to avoid recalculation of the 
 * isotrophic prediction of the antenna during optimization.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration 
 *              parameters;
 *
 */
static void
override_clutter_losses (Parameters    *params,
                         Tx_parameters *tx_params)
{
    int r, c;

    //
    // substract the already applied clutter loss 
    // in order to have the clean isotrophic prediction 
    //
    for (r = 0; r < tx_params->nrows; r ++)
    {
        for (c = 0; c < tx_params->ncols; c ++)
        {
            //
            // use only valid path-loss values
            //
            if (!isnan (tx_params->m_loss[r][c]))
            {
                //
                // get the clutter loss, based on the category of this point
                //
                int clutter_category = (int) tx_params->m_clut[r][c];

                if (clutter_category != params->cell_null_value)
                    tx_params->m_loss[r][c] -= params->clutter_loss[clutter_category];
            }
        }
    }
}



/**
 * Returns the squared error between the prediction and the field 
 * measurements of this transmitter, i.e. the objective function 
 * evaluation, used by the optimization algorithm.
 * Most of this function is executed on the GPU, if available.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration 
 *              parameters;
 * radio_zone   radio zone for which the objective function is calculated;
 *
 */ 
static double
squared_error_on_gpu (Parameters    *params,
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

    return ret_value;
}



/**
 * Returns the sum of the squared mean and the variance, calculated
 * between the prediction and the field measurements of this transmitter, 
 * i.e., the objective function evaluation used by the optimization algorithm.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration 
 *                  parameters;
 * radio_zone       radio zone for which the objective function is calculated;
 * sol_vector       solution vector over which the objective function is 
 *                  calculated;
 * sol_vector_len   length of the received solution vector;
 * recalculate      whether to recalculate the isotrophic prediction, e.g. for 
 *                  clutter optimization there is no need to recalculate the 
 *                  prediction in every iteration. This is useful to speed-up 
 *                  the process on the CPU;
 * clut_opt         whether the solution being evaluated is clutter losses or
 *                  the E/// parameters.-
 *
 */ 
static double
obj_func (Parameters    *params,
          Tx_parameters *tx_params,
          const char     radio_zone,
          const double  *sol_vector,
          const int      sol_vector_len,
          const char     recalculate,
          const char     clut_opt)
{
    const int Max_meas_count = 102400;
    int r, c;
    double mean = 0.0, variance = 0.0;
    double *temp_variance = (double *) calloc (Max_meas_count, 
                                               sizeof (double));
    double ret_value = 0;

    //
    // copy the solution values to the corresponding internal structure,
    // used for coverage calculation
    //
    if (clut_opt)
    {
        for (r = 0; r < sol_vector_len; r ++)
            params->clutter_loss[r] = sol_vector[r];
    }
    else
    {
        for (r = 0; r < sol_vector_len; r ++)
            tx_params->eric_params[r] = sol_vector[r];
    }

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
        measure_time ("Squared-error function on GPU");
#endif
        ret_value = squared_error_on_gpu (params,
                                          tx_params,
                                          radio_zone);
    }
    else
    {
        if (recalculate)
        {
#ifdef _PERFORMANCE_METRICS_
            measure_time ("E/// on CPU");
#endif
            eric_pathloss_on_cpu (params,
                                  tx_params);
        }
#ifdef _PERFORMANCE_METRICS_
        measure_time (NULL);
        measure_time ("Apply antenna losses and calculate objective function on CPU");
#endif
        //
        // objective function calculation, for 
        // each point in the path loss matrix ...
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
                if (pl > 0.0)
                {
                    //
                    // apply the clutter losses only if the isotrophic 
                    // prediction has not been recalculated
                    //
                    if (!recalculate)
                    {
                        //
                        // get the clutter loss, based on the category of this point
                        //
                        int clutter_category = (int) tx_params->m_clut[r][c];

                        if (clutter_category != params->cell_null_value)
                            pl += params->clutter_loss[clutter_category];
                    }
                    //
                    // apply the antenna loss
                    //
                    pl += tx_params->m_antenna_loss[r][c];

                    //
                    // look for a field measurement there
                    //
                    double fm = tx_params->m_field_meas[r][c];
                    if (!isnan (fm))
                    {
                        //printf ("*** DEBUG: A0...A3 [%.5f, %.5f, %.5f, %.5f], pathloss is [%.5f]\n", tx_params->eric_params[0],
                        //                                                                             tx_params->eric_params[1],
                        //                                                                             tx_params->eric_params[2],
                        //                                                                             tx_params->eric_params[3],
                        //                                                                             pl);
                        if (tx_params->field_meas_count == Max_meas_count)
                        {
                            fprintf (stderr, 
                                     "*** ERROR: cannot save all field measurements. Current array size is [%d].\n",
                                     Max_meas_count);
                            exit (1);
                        }
                        mean += (tx_params->tx_power - pl) - fm;
                        temp_variance[tx_params->field_meas_count] = (tx_params->tx_power - pl) - fm;
                        tx_params->field_meas_count ++;
                    }
                }
             }
        }
        mean /= tx_params->field_meas_count;
        for (r = 0; r < tx_params->field_meas_count; r ++)
            variance += (temp_variance[r] - mean) * 
                        (temp_variance[r] - mean);
        variance /= tx_params->field_meas_count;
        //
        // the unified objective-function value for mean and variance
        //
        ret_value = mean*mean + variance;
        //
        // free allocated memory
        //
        free (temp_variance);
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
 * clut_opt     whether to optimize clutter losses or E/// parameters.-
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
    const char      clut_opt)
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

    printf ("*** INFO: The dimension of the optimization problem is [%d]\n", D);

    /* Starting timer */
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

        //
        // objective-function evaluation
        //
        popul[i][D] = obj_func (params,
                                tx_params,
                                radio_zone,
                                popul[i],
                                D,
                                1,
                                clut_opt);

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

        //
        // objective-function evaluation
        //
        U[D] = obj_func (params,
                         tx_params,
                         radio_zone,
                         U,
                         D,
                         1,
                         clut_opt);
        numofFE++;

        printf ("Generation %d/%d\tscore %20.10f\n", k, Gmax, U[D]);

        /* 
        // DEBUG - print current solution out
        //
        for (j=0; j <= D; j++)
            printf ("%.3f, ", U[j]);
        printf ("\n");*/

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
           for (j=0; j < D; j++)
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
 * Initialize the optimization environment and structures for both 
 * worker and master-based optimization processes.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration 
 *              parameters;
 * search_low   pointer to vector containing the minimum values each 
 *              solution component may take;
 * search_up    pointer to vector containing the maximum values each 
 *              solution component may take;
 * clut_opt     whether to optimize the clutter losses or the E///
 *              parameters.-
 *
 */
static void 
init_optimize (Parameters    *params,
               Tx_parameters *tx_params,
               double       **search_low,
               double       **search_up,
               const char     clut_opt)
{
    int    i;
    size_t search_len;
    double search_min, search_max;

    //
    // define lower and upper bounds for each search-vector component,
    // i.e., solutions should be within these limits
    //
    if (clut_opt)
    {
        search_len = (size_t) params->clutter_category_count;
        search_min = 0.0;
        search_max = 40.0;
    }
    else
    {
        search_len = (size_t) params->clutter_category_count;
        search_min = -100.0;
        search_max = 100.0;
    }
    //
    // apply the above values
    //
    *search_low = (double *) calloc (search_len,
                                     sizeof (double));
    *search_up  = (double *) calloc (search_len,
                                     sizeof (double));
    for (i = 0; i < search_len; i ++)
    { 
        (*search_low)[i] = search_min;
        (*search_up)[i]  = search_max;
    }
    
    //
    // calculate the coverage for the first time to initialize all needed structures
    //
    coverage (params,
              tx_params,
              0);
    //
    // if the first coverage calculation happened on the GPU,
    // we need to refresh the memory buffers on the host
    //
    if (params->use_gpu)
    {
        size_t dbl_buff_size = tx_params->nrows * 
                               tx_params->ncols * 
                               sizeof (tx_params->m_antenna_loss[0][0]);
        size_t ch_buff_size = tx_params->nrows * 
                              tx_params->ncols * 
                              sizeof (tx_params->m_radio_zone[0][0]);
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
}



/**
 * Optimizes the parameters of the prediction model to fit field measurements
 * within different radio zones using a deterministic (analytical) approach.
 * It also optimizes the clutter-category losses using the differential 
 * evolution (DE) metaheuristic algorithm.
 * Optimization happens locally, i.e. within the current worker process.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration 
 *              parameters;
 * clut_opt     whether to optimize the clutter losses or the E/// parameters.-
 *
 */
void 
optimize_on_worker (Parameters    *params,
                    Tx_parameters *tx_params,
                    const char     clut_opt)
{
    //
    // vectors defining the range of each solution component
    //
    double *search_low = NULL;
    double *search_up  = NULL;

    init_optimize (params,
                   tx_params,
                   &search_low,
                   &search_up,
                   clut_opt);
    //
    // calculate parameter approximation for E/// prediction model
    //
    parameter_fine_tuning (params,
                           tx_params);
    //
    // ... and its mean-squared error value
    //
    double score;
    double *sol_vector;
    int    sol_vector_len;

    if (clut_opt)
    {
        sol_vector = params->clutter_loss;
        sol_vector_len = params->clutter_category_count;
    }
    else
    {
        sol_vector = params->eric_params;
        sol_vector_len = 4;
    }
    score = obj_func (params,
                      tx_params,
                      _RADIO_ZONE_MAIN_BEAM_ON_,
                      sol_vector,
                      sol_vector_len,
                      1,
                      clut_opt);
    fprintf (stdout, 
             "*** INFO: optimal values for E/// (%s) have score %g\n",
             tx_params->tx_name,
             score);

    //
    // parameters for the optimization algorithm
    //
    int search_len, population_size;

    if (clut_opt)
        search_len = params->clutter_category_count;
    else
        search_len = 4;

    population_size = 20 * search_len;

#ifdef _DEBUG_INFO_
    //
    // DEBUG: run only one iteration of the optimization process
    //
    de (params,
        tx_params,
        search_len,
        search_len,
        1,
        0.9,
        0.9,
        0,
        _RADIO_ZONE_MAIN_BEAM_ON_,
        search_low,
        search_up,
        clut_opt);
#else
    //
    // start optimization
    //
    de (params,
        tx_params,
        search_len,
        population_size,
        params->use_opt,
        0.9,
        0.9,
        1,
        _RADIO_ZONE_MAIN_BEAM_ON_,
        search_low,
        search_up,
        clut_opt);
#endif
    //
    // free reserved buffer
    //
    free (search_up);
    free (search_low);
}



/**
 * Optimizes the parameters of the prediction model to fit field measurements
 * within different radio zones using a deterministic (analytical) approach.
 * It also calculates the objective-function value using the solutions 
 * received from the master process, where the optimization of the 
 * clutter-category losses (using differential evolution (DE)) runs.
 *
 * params       a structure holding configuration parameters which are 
 *              common to all transmitters;
 * tx_params    a structure holding transmitter-specific configuration 
 *              parameters;
 * comm         the communicator used to receive messages from the master;
 *
 */
void 
optimize_from_master (Parameters    *params,
                      Tx_parameters *tx_params,
                      MPI_Comm      *comm)
{
    int has_finished = 0;

    //
    // vector used to send the squared error and the number of
    // field measurements within the prediction area (used to
    // calculate the mean square error on the master side)
    //
    double score [2];

    //
    // vectors defining the range of each solution component
    //
    double *search_low = NULL;
    double *search_up  = NULL;

    //
    // initialize the optimization environment for this worker
    //
    init_optimize (params,
                   tx_params,
                   &search_low,
                   &search_up,
                   1);
    //
    // calculate parameter approximation for E/// prediction model
    //
    parameter_fine_tuning (params,
                           tx_params);
    //
    // ... and its objective-function value, using the default 
    // clutter-category losses
    //
    score[0] = obj_func (params,
                         tx_params,
                         _RADIO_ZONE_MAIN_BEAM_ON_,
                         params->clutter_loss,
                         params->clutter_category_count,
                         1,
                         1);
    score[1] = (double) tx_params->field_meas_count;
    fprintf (stdout, 
             "*** INFO: optimal values for E/// (%s) have score %g\n",
             tx_params->tx_name,
             score[0]);
    fflush (stdout);
    //
    // a vector to keep the received solution
    //
    double *sol_vector = (double *) calloc (params->clutter_category_count,
                                            sizeof (params->clutter_loss[0]));
    //
    // get ready to evaluate the received solutions
    //
    override_clutter_losses (params,
                             tx_params);
    while (!has_finished)
    {
        //
        // receive the solution vector
        //
        MPI_Bcast (sol_vector,
                   params->clutter_category_count,
                   MPI_DOUBLE,
                   _COVERAGE_MASTER_RANK_,
                   *comm);
        //
        // maybe this worker should stop working?
        //
        if (sol_vector[0] == _PI_)
        {
            int i;
            for (i = 0; i < params->clutter_category_count; i ++)
                if (sol_vector[i] == _PI_)
                    has_finished = 1;
                else
                    has_finished = 0;
        }
        else
        {
            //
            // calculate the objective-function value for this solution ...
            //
            score[0] = obj_func (params,
                                 tx_params,
                                 _RADIO_ZONE_MAIN_BEAM_ON_,
                                 sol_vector,
                                 params->clutter_category_count,
                                 0,
                                 1);
            score[1] = (double) tx_params->field_meas_count;
            //
            // ... and send it back to the master process
            //
            MPI_Send (score,
                      2,
                      MPI_DOUBLE,
                      _COVERAGE_MASTER_RANK_,
                      1,
                      *comm);
        }
    }
    //
    // free reserved vectors
    //
    free (sol_vector);
    free (search_up);
    free (search_low);
}

