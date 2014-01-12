#include <performance_metric.h>

#include "worker/coverage.h"
#include "master/master.h"




/**
 * Reads the received file into memory, returning the number of bytes 
 * saved into the buffer.
 *
 * file_name        The file to be read;
 * content_buffer   the buffer to which the read data is saved.-
 *
 */
static int 
read_file_into_memory (const char *file_name,
                       char *content_buffer)
{
    char *read_ptr = &(content_buffer[0]);
    FILE *fr = fopen (file_name, "r");

    if (fr == NULL)
    {
        fprintf (stderr,
                 "*** ERROR: Cannot open [%s] for reading\n", 
                 file_name);
        fflush (stderr);
        exit (1);
    }

    while (fgets (read_ptr, 1024, fr) != NULL)
       read_ptr += strlen (read_ptr);

    fclose (fr);

    return (strlen (content_buffer) + 1);
}



/**
 * Entry point of the GRASS module.-
 */
int main (int argc, char **argv)
{
    int i;
    struct GModule *module;
    struct Option  *ini_file, *tx_ini_sections, *output,
                   *eric_a0, *eric_a1, *eric_a2, *eric_a3,
                   *worker_opt, *master_opt;
    struct Flag    *use_mpi, *use_gpu;

    //
    // initialize the GIS environment
    //
    G_gisinit (argv[0]);

    //
    // initialize the module itself
    //
    module = G_define_module ( );
    module->keywords = "raster, coverage, E/// 9999, directional antennas";
    module->description = "Parallel coverage-prediction module using the E/// 9999 model and antenna diagrams.-";
  
    //
    // define the module options ...
    //
    ini_file = G_define_option ( );
    ini_file->key = "ini_file";
    ini_file->type = TYPE_STRING;
    ini_file->required = YES;
    ini_file->label = "Full path to the parameters INI file";

    tx_ini_sections = G_define_option ( );
    tx_ini_sections->key = "tx_ini_sections";
    tx_ini_sections->type = TYPE_STRING;
    tx_ini_sections->required = YES;
    tx_ini_sections->label = "A comma-separated list of the names of the INI-file sections containing the transmitters' configuration to be used";

    output = G_define_option ( );
    output->key = "output_raster";
    output->type = TYPE_STRING;
    output->required = NO;
    output->label = "Output raster name. Without it, output is given in standard output.";

    eric_a0 = G_define_option ( );
    eric_a0->key = "A0";
    eric_a0->type = TYPE_DOUBLE;
    eric_a0->required = NO;
    eric_a0->answer = "38.0";
    eric_a0->label = "Value of the A0 parameter for the E/// 9999 model";

    eric_a1 = G_define_option ( );
    eric_a1->key = "A1";
    eric_a1->type = TYPE_DOUBLE;
    eric_a1->required = NO;
    eric_a1->answer = "32.0";
    eric_a1->label = "Value of the A1 parameter for the E/// 9999 model";

    eric_a2 = G_define_option ( );
    eric_a2->key = "A2";
    eric_a2->type = TYPE_DOUBLE;
    eric_a2->required = NO;
    eric_a2->answer = "-12.0";
    eric_a2->label = "Value of the A2 parameter for the E/// 9999 model";

    eric_a3 = G_define_option ( );
    eric_a3->key = "A3";
    eric_a3->type = TYPE_DOUBLE;
    eric_a3->required = NO;
    eric_a3->answer = "0.1";
    eric_a3->label = "Value of the A3 parameter for the E/// 9999 model";

    use_mpi = G_define_flag ( );
    use_mpi->key = 'm';
    use_mpi->description = "Use the MPI implementation";

    use_gpu = G_define_flag ( );
    use_gpu->key = 'g';
    use_gpu->description = "Use the GPU implementation. Implies -m.";

    worker_opt = G_define_option ( );
    worker_opt->key = "t";
    worker_opt->type = TYPE_INTEGER;
    worker_opt->required = NO;
    worker_opt->answer = "0";
    worker_opt->label = "Optimize the clutter-category losses locally, i.e. per worker, for the given number of generations. Implies -m.";

    master_opt = G_define_option ( );
    master_opt->key = "p";
    master_opt->type = TYPE_INTEGER;
    master_opt->required = NO;
    master_opt->answer = "0";
    master_opt->label = "Optimize the clutter-category losses globally, i.e. workers calculate only the objective function, for the given number of generations. Implies -m.";

    //
    // ... and parse them
    //
    if (G_parser (argc, argv) < 0)
	    exit (1);

#ifdef _PERFORMANCE_METRICS_
    measure_time ("Read input data");
#endif
    //
    // allocate the parameters structure
    //
    Parameters *params = (Parameters *) malloc (sizeof (Parameters));

    //
    // parameter values for prediction model
    //
    sscanf (eric_a0->answer, "%lf", &(params->eric_params[0]));
    sscanf (eric_a1->answer, "%lf", &(params->eric_params[1]));
    sscanf (eric_a2->answer, "%lf", &(params->eric_params[2]));
    sscanf (eric_a3->answer, "%lf", &(params->eric_params[3]));

    //
    // number of evaluations in case of optimization
    //
    sscanf (worker_opt->answer, "%d", &(params->use_opt));
    sscanf (master_opt->answer, "%d", &(params->use_master_opt));

    //
    // flag for enabling the GPU implementation on the workers
    //
    params->use_gpu = use_gpu->answer;

    if (params->use_gpu)
    {
        use_mpi->answer = 1;
        fprintf (stdout, 
                 "*** INFO: GPU hardware will be used on the workers, if available\n");
    }
    if (params->use_master_opt)
    {
        params->use_opt = params->use_master_opt;
        use_mpi->answer = 1;
        fprintf (stdout, 
                 "*** INFO: Master-optimization mode enabled\n");
    }
    else if (params->use_opt)
    {
        use_mpi->answer = 1;
        fprintf (stdout, 
                 "*** INFO: Worker-optimization mode enabled\n");
    }

    //
    // read the whole configuration INI file into memory
    //
    char *content = (char *) calloc (1024 * 1024, sizeof (char));
    params->ini_file_content_size = read_file_into_memory (ini_file->answer,
                                                           content);
    params->ini_file_content = (char *) malloc (params->ini_file_content_size);
    memcpy (params->ini_file_content, 
            content, 
            params->ini_file_content_size);

    FILE *ini_file_stream = fmemopen (params->ini_file_content,
                                      params->ini_file_content_size,
                                      "r");
    //
    // initialize the coverage calculation
    //
    init_coverage (ini_file_stream,
                   tx_ini_sections->answer,
                   params);
    //
    // free memory and close input stream
    //
    free (content);
    fclose (ini_file_stream);

#ifdef _PERFORMANCE_METRICS_
    measure_time (NULL);
#endif

    //
    // ... and execute it
    //
    if (use_mpi->answer)
    {
        init_mpi (argc,
                  argv,
                  params);
    }
    else
    {
        if (params->ntx > 1)
            fprintf (stderr, 
                     "*** WARNING: the results of the serial execution mode have not been tested!\n");
        //
        // set the tuning parameters for the prediction model
        //
        sscanf (eric_a0->answer, "%lf", &(params->tx_params->eric_params[0]));
        sscanf (eric_a1->answer, "%lf", &(params->tx_params->eric_params[1]));
        sscanf (eric_a2->answer, "%lf", &(params->tx_params->eric_params[2]));
        sscanf (eric_a3->answer, "%lf", &(params->tx_params->eric_params[3]));

        //
        // calculate the coverage for all transmitters, serially
        //
        for (i = 0; i < params->ntx; i ++)
        {
            Tx_parameters *tx_params = &(params->tx_params[i]);
            coverage (params,
                      tx_params,
                      0);

            //
            // calculation finished, do we have to write the raster output?
            //
            if (output->answer == NULL)
            {
                pthread_t dump_thread;
                //
                // start result dump
                //
                int err = pthread_create (&dump_thread, 
                                          NULL,
                                          &output_to_stdout, 
                                          (void *) params);
                if (err)
                {
                    fprintf (stderr,
                             "*** ERROR: number (%d) while creating result dump thread\n", 
                             err);
                    fflush (stderr);
                    exit (1);
                }
                else
                {
                    //
                    // wait for it to finish
                    //
                    err = pthread_join (dump_thread, NULL);
                }
            }
            else
            {
                int outfd, row, col;
                void *outrast = G_allocate_raster_buf (FCELL_TYPE);    // output buffer
                double path_loss_num;
                FCELL  null_f_out;
                G_set_f_null_value (&null_f_out, 1);   

                // controlling, if we can write the raster 
                if ((outfd = G_open_raster_new (output->answer, FCELL_TYPE)) < 0)
                    G_fatal_error ("Unable to create raster map <%s>", 
                                   output->answer);

                for (row = 0; row < params->nrows; row++)
                {
                    G_percent(row, params->nrows, 2);
                    for (col = 0; col < params->ncols; col++) 
                    {
                        path_loss_num = params->m_loss[row][col];
                        if (path_loss_num == 0)
                            ((FCELL *) outrast)[col] = null_f_out;
                        else
                            ((FCELL *) outrast)[col] = (FCELL)path_loss_num;
                    }
                    // write raster row to output raster map
                    if (G_put_raster_row (outfd, outrast, FCELL_TYPE) < 0)
                        G_fatal_error ("Failed writing raster map <%s>", 
                                       output->answer);
                }
                G_close_cell (outfd);
                G_free (outrast);
            }
        }
    }
    //
    // deallocate memory before exiting
    //
    free_tx_params (params,
                    params->tx_params);
    free (params->ini_file_content);
    free (params->tx_params);
    free (params);

    exit (EXIT_SUCCESS);
}

