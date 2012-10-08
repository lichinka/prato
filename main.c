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
static int read_file_into_memory (const char *file_name,
                                  char *content_buffer)
{
    char *read_ptr = &(content_buffer[0]);
    FILE *fr = fopen (file_name, "r");

    if (fr == NULL)
    {
        fprintf (stderr, "ERROR Cannot open [%s] for reading\n", file_name);
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
int main (int argc, char *argv [])
{
    struct GModule *module;
    struct Option  *ini_file, *tx_ini_sections, *output;
    struct Flag    *use_mpi;

    //
    // initialize the GIS environment
    //
    G_gisinit (argv[0]);

    //
    // initialize the module itself
    //
    module = G_define_module ( );
    module->keywords = _("raster, coverage, E/// 9999, directional antennas");
    module->description = _("Coverage module using the E/// 9999 model and antenna diagrams.-");
  
    //
    // define the module options ...
    //
    ini_file = G_define_option ( );
    ini_file->key = "ini_file";
    ini_file->type = TYPE_STRING;
    ini_file->required = YES;
    ini_file->label = _("Full path to the parameters INI file");	

    tx_ini_sections = G_define_option ( );
    tx_ini_sections->key = "tx_ini_sections";
    tx_ini_sections->type = TYPE_STRING;
    tx_ini_sections->required = YES;
    tx_ini_sections->label = _("A comma-separated list of the names of the INI-file sections containing the transmitters' configuration to be used");

    output = G_define_option ( );
    output->key = "output_raster";
    output->type = TYPE_STRING;
    output->required = NO;
    output->label = _("Output raster name. Without it, output is given in standard output.");

    use_mpi = G_define_flag ( );
    use_mpi->key = 'p';
    use_mpi->description = _("Whether to use the MPI implementation");

    //
    // ... and parse them
    //
    if (G_parser (argc, argv) < 0)
	    exit (EXIT_FAILURE);

    //
    // read the whole configuration INI file into memory
    //
    char content [1024 * 1024];
    Parameters *params = (Parameters *) malloc (sizeof (Parameters));
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
    // initialize the coverage calculation ...
    //
    init_coverage (ini_file_stream,
                   tx_ini_sections->answer,
                   params);
    fclose (ini_file_stream);

    //
    // ... and execute it
    //
    double ericsson_params [4] = {38.0, 32.0, -12.0, 0.1};

#ifdef _PERFORMANCE_METRICS_
    for (i = 0; i < 3; i ++)
#endif
    if (use_mpi->answer)
    {
        coverage_mpi (argc,
                      argv,
                      params);
    }
    else
    {
        if (params->ntx > 1)
            fprintf (stderr, "WARNING Only the first transmitter will be processed\n");
        coverage (params,
                  &(params->tx_params[0]),
                  ericsson_params,
                  4);
        //
        // calculation finished, do we have to write the raster output?
        //
        if (output->answer == NULL)
            output_to_stdout (params);
        else
        {
            int outfd, row, col;
            void *outrast = G_allocate_raster_buf (FCELL_TYPE);    // output buffer
            double path_loss_num;
            FCELL  null_f_out;
            G_set_f_null_value (&null_f_out, 1);   

            // controlling, if we can write the raster 
            if ((outfd = G_open_raster_new (output->answer, FCELL_TYPE)) < 0)
                G_fatal_error(_("Unable to create raster map <%s>"), output->answer);

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
                    G_fatal_error (_("Failed writing raster map <%s>"), output->answer);
            }
            G_close_cell (outfd);
            G_free (outrast);
        }
    }
    //
    // deallocate memory before exiting
    //
    free (params->tx_params);
    free (params->ini_file_content);
    free (&(params->m_loss[0][0]));
    free (params->m_loss);
    free (&(params->m_clut[0][0]));
    free (params->m_clut);
    free (&(params->m_dem[0][0]));
    free (params->m_dem);
    free (params);

    exit (EXIT_SUCCESS);
}

