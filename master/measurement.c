#define _GRASS_DB_DRIVER_       "pg"
#define _GRASS_DB_STRING_       "host=localhost,dbname=umts"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <grass/gis.h>
#include <grass/dbmi.h>
#include <grass/glocale.h>




/**
 * Loads the field-measurements array from a raster map.
 *
 * mapname     the full name of the raster map from which the matrix
 *              should be loaded;
 * field_meas   a 2D-matrix where the loaded data are saved.-
 */
void 
load_field_measurements_from_map (const char *mapname,
                                  double **field_meas)
{
    int row, col, infd, errno;
    char *mapset;
 
    //
    // NULL if the map was not found in any mapset
    //
    mapset = G_find_cell2 (mapname, "");
    if (mapset == NULL)
        G_fatal_error ("Field measurements raster map <%s> not found", 
                       mapname);
    //
    // open raster map
    //
    if ((infd = G_open_cell_old (mapname, mapset)) < 0)
        G_fatal_error ("Unable to open raster map <%s>", 
                       mapname);
    //
    // read map metadata, making sure it matches the current region
    //
    struct Cell_head *metadata = (struct Cell_head *) malloc (sizeof (struct Cell_head));
    struct Cell_head *window   = (struct Cell_head *) malloc (sizeof (struct Cell_head));

	G_get_window (window);
    errno = G_get_cellhd (mapname,
                          mapset,
                          metadata);
    if (errno == 0)
    {
        //
        // compare map metadata with current region
        //
        if ((window->east != metadata->east) ||
            (window->west != metadata->west) ||
            (window->north != metadata->north) ||
            (window->south != metadata->south) ||
            (window->ew_res != metadata->ew_res) ||
            (window->ns_res != metadata->ns_res))
            fprintf (stderr, 
                     "*** WARNING: Metadata of field measurements map does not match with DEM.\n");
        else
        {
            //
            // allocate the reading buffer for DEM 
            //
            void *inrast = G_allocate_raster_buf (FCELL_TYPE);

            //
            // read field measurements into the 2D matrix
            //
            for (row = 0; row < metadata->rows; row ++) 
            {	
                if (G_get_raster_row (infd, inrast, row, FCELL_TYPE) < 0)
                    G_fatal_error ("Unable to read raster map <%s> row %d", 
                                   mapname, 
                                   row);
                for (col = 0; col < metadata->cols; col ++) 
                { 
                    FCELL f_in = ((FCELL *) inrast)[col];
                    field_meas[row][col] = (double) f_in;
                }
            }
            //
            // deallocate reading buffer
            //
            G_free (inrast);
        }
    }
    else
        G_fatal_error ("Unable to open raster map <%s>", 
                       mapname);
    //
    // deallocate memory elements
    //
    free (window);
    free (metadata);
    G_close_cell (infd);
}



/**
 * Loads the field-measurements array from a binary file on disk.
 * The first 96 bits contain the: size of an element (in bytes),
 * number of rows of the matrix, and number of columns of the 
 * matrix, in that order.
 *
 * file_name    the full path of the file from which the matrix
 *              should be loaded;
 * field_meas   a 2D-matrix where the loaded data are saved.-
 */
void load_field_measurements_from_file (const char *file_name,
                                        double **field_meas)
{
    int i;
    FILE *fp = fopen (file_name, "rb");
 
    if (fp == NULL) 
        G_fatal_error ("Failed to open file '%s'", file_name);

    //
    // load matrix metadata
    // 
    int element_size, nrows, ncols;
    fread (&element_size, sizeof (element_size), 1, fp);
    fread (&nrows, sizeof (nrows), 1, fp);
    fread (&ncols, sizeof (ncols), 1, fp);
    if (element_size != (int) sizeof (**field_meas))
        G_fatal_error ("Loaded element size (%d) does not match the allocated one (%ld)", element_size,
                                                                                          sizeof (**field_meas));
    //
    // load the matrix
    //
    for (i = 0; i < nrows; i ++)
        fread (field_meas[i], element_size, ncols, fp);
 
    fclose (fp);
}



/**
 * Dumps the field-measurements array into a binary file on disk.
 * The first 96 bits contain the: size of an element (in bytes),
 * number of rows of the matrix, and number of columns of the 
 * matrix, in that order.
 *
 * file_name    the full path of the file into which the matrix
 *              contents are saved;
 * nrows        number of rows in the 2D matrix;
 * ncols        number of cols in the 2D matrix;
 * field_meas   a 2D-matrix where the field measurements are saved.-
 */
void dump_field_measurements (const char *file_name,
                              const int nrows,
                              const int ncols,
                              double **field_meas)
{
    int i;
    FILE *fp = fopen (file_name, "wb");
 
    if (fp == NULL) 
        G_fatal_error ("Failed to create file '%s'", file_name);

    //
    // first save some data about the matrix, namely
    // the size of an element, the number of rows, and 
    // number of columns, each number being a 32 bit int
    // 
    int element_size = (int) sizeof (**field_meas);
    fwrite (&element_size, sizeof (element_size), 1, fp);
    fwrite (&nrows, sizeof (nrows), 1, fp);
    fwrite (&ncols, sizeof (ncols), 1, fp);
    //
    // now dump the matrix
    //
    for (i = 0; i < nrows; i ++)
        fwrite (field_meas[i], element_size, ncols, fp);
 
    fclose (fp);
}



/**
 * Loads the measurements from the database into a memory array.
 *
 * tx_name      Transmitter's name;
 * tx_east      transmitter's eastern coordinate (m);
 * tx_north     transmitter's northern coordinate (m);
 * radius       calculation radius around the given transmitter (km);
 * raster_west  western-most raster coordinate;
 * raster_north northern-most raster coordinate;
 * ew_res       east/west raster resolution;
 * ns_res       north/south raster resolution;
 * field_meas   a 2D-matrix where the field measurements are saved.-
 *
 */
void load_field_measurements_from_db (const char *tx_name,
                                      const double tx_east,
                                      const double tx_north,
                                      const double radius,
                                      const double raster_west,
                                      const double raster_north,
                                      const double ew_res,
                                      const double ns_res,
                                      double **field_meas)
{
    int more, buf_size = 2048;
    char buf [buf_size];
    dbDriver *driver;
    dbTable  *table;
    dbCursor cursor;
    dbColumn *column;
    dbValue  *value;
    dbString sql;
    //
    // calculate a bounding box around the transmitter, 
    // with the given radius in kilometers
    //
    int low_left_x = (int) ceil (tx_east - radius * 1000);
    int low_left_y = (int) ceil (tx_north - radius * 1000);
    int upper_right_x = (int) ceil (tx_east + radius * 1000);
    int upper_right_y = (int) ceil (tx_north + radius * 1000);

    // open the DB driver 
    driver = db_start_driver_open_database (_GRASS_DB_DRIVER_,
                                            _GRASS_DB_STRING_);
    if (driver == NULL)
        G_fatal_error ("Unable to open database");
    //
    // First we need the record ID of the transmitter used. 
    // This little step improves the speed of the long query by more
    // than 20x, e.g.
    //
    // SELECT id_cell FROM network.cell WHERE name = 'SBANOVA'
    //
    snprintf (buf, buf_size,
              "SELECT id_cell FROM network.cell WHERE name = '%s'",
              tx_name);
    db_init_string (&sql);
    db_set_string  (&sql, buf);
    //printf ("\t*** DEBUG: Executing \n\t%s\n", db_get_string (&sql));
    if (db_open_select_cursor (driver, &sql, &cursor, DB_SEQUENTIAL) != DB_OK) 
        G_fatal_error ("Unable to open select cursor:\n\t%s",
		       db_get_string (&sql));

    table = db_get_cursor_table (&cursor);
    //
    // read the ID of the transmitter record in the database
    //
    int tx_db_id = -1;
    if (db_fetch (&cursor, DB_NEXT, &more) == DB_OK)
    {
        column = db_get_table_column (table, 0);
        int ctype = db_sqltype_to_Ctype (db_get_column_sqltype (column));
        if (ctype != DB_C_TYPE_INT)
            G_fatal_error ("Column type for 'id_cell' is not INT as expected.");
        value = db_get_column_value (column);
        tx_db_id = db_get_value_int (value);
        if (db_close_cursor (&cursor) != DB_OK)
            G_fatal_error ("Cannot close database cursor.");
    }
    else
        G_fatal_error ("Cannot read the transmitter ID from the database.");

    //
    // Select the needed field measurements from the database, e.g.
    //
    // SELECT X(GK.location) - X(GK.location)::int % 100.0, 
    //        Y(GK.location) - Y(GK.location)::int % 100.0, 
    //        AVG(rscp)
    //   FROM (SELECT DISTINCT ON (M.id_coordinate) M.time, M.id_coordinate, rscp
    //           FROM (SELECT time, id_cell, id_coordinate, rscp FROM romes.romes_2012
    //                 UNION
    //                 SELECT time, id_cell, id_coordinate, rscp FROM romes.romes_2011
    //                 UNION
    //                 SELECT time, id_cell, id_coordinate, rscp FROM romes.romes_2010
    //                 UNION
    //                 SELECT time, id_cell, id_coordinate, rscp FROM romes.romes_2009
    //                 UNION
    //                 SELECT time, id_cell, id_coordinate, rscp FROM romes.romes_2008) M
    //           WHERE M.id_cell=6865
    //        ORDER BY M.id_coordinate, M.time DESC) H
    //            JOIN gis.coordinate_gk GK ON GK.id_coordinate = H.id_coordinate 
    //             AND GK.location && setsrid (makebox2d ('POINT (164000 580000)'::geometry,
    //                                                    'POINT (184000 600000)'::geometry), 2710)
    //        GROUP BY X(GK.location) - X(GK.location)::int % 100.0, 
    //                 Y(GK.location) - Y(GK.location)::int % 100.0;
    //
    snprintf (buf, buf_size,
              "SELECT (X(GK.location) - X(GK.location)::int %% %.2f) AS X, \
                      (Y(GK.location) - Y(GK.location)::int %% %.2f) AS Y, \
                      AVG(rscp) AS rscp \
                 FROM (SELECT DISTINCT ON (M.id_coordinate) M.time, M.id_coordinate, rscp \
                         FROM (SELECT time, id_cell, id_coordinate, rscp FROM romes.romes_2012 \
                               UNION \
                               SELECT time, id_cell, id_coordinate, rscp FROM romes.romes_2011 \
                               UNION \
                               SELECT time, id_cell, id_coordinate, rscp FROM romes.romes_2010 \
                               UNION \
                               SELECT time, id_cell, id_coordinate, rscp FROM romes.romes_2009 \
                               UNION \
                               SELECT time, id_cell, id_coordinate, rscp FROM romes.romes_2008) M \
                WHERE M.id_cell = %d \
             ORDER BY M.id_coordinate, M.time DESC) H \
                 JOIN gis.coordinate_gk GK ON GK.location && \
                                              setsrid (makebox2d ('POINT (%d %d)'::geometry, \
                                                                  'POINT (%d %d)'::geometry), 2710) \
                                          AND GK.id_coordinate = H.id_coordinate \
             GROUP BY (X(GK.location) - X(GK.location)::int %% %.2f), \
                      (Y(GK.location) - Y(GK.location)::int %% %.2f)",
                ns_res, ew_res,
                tx_db_id,
                low_left_y, low_left_x,
                upper_right_y, upper_right_x,
                ns_res, ew_res);
    db_set_string  (&sql, buf);
    // printf ("\t*** DEBUG: Executing \n\t%s\n", db_get_string (&sql));
    if (db_open_select_cursor (driver, &sql, &cursor, DB_SEQUENTIAL) != DB_OK) 
        G_fatal_error ("Unable to open select cursor:\n\t'%s'", db_get_string (&sql));

    table = db_get_cursor_table (&cursor);
    //
    // read the field measurements from the database
    //
    while (db_fetch (&cursor, DB_NEXT, &more) == DB_OK && more) 
    {
        // read measurement
        column = db_get_table_column (table, 2);
        int ctype = db_sqltype_to_Ctype (db_get_column_sqltype (column));
        if (ctype != DB_C_TYPE_DOUBLE)
            G_fatal_error ("Measurement (RSCP) column must be of type double");
        value = db_get_column_value (column);
        double m_rscp = db_get_value_double (value);

        //
        // read the geographic coordinates of this measurement
        //
        // X
        double m_x = -1;
        column = db_get_table_column (table, 0);
        ctype = db_sqltype_to_Ctype (db_get_column_sqltype (column));
        if (ctype != DB_C_TYPE_INT && ctype != DB_C_TYPE_DOUBLE)
            G_fatal_error ("Coordinate column must be of type integer or double");
        value = db_get_column_value (column);
        if (ctype == DB_C_TYPE_INT)
            m_x = (double)db_get_value_int (value);
        else
            m_x = db_get_value_double (value);

        // Y
        double m_y = -1;
        column = db_get_table_column (table, 1);
        ctype = db_sqltype_to_Ctype (db_get_column_sqltype (column));
        if (ctype != DB_C_TYPE_INT && ctype != DB_C_TYPE_DOUBLE)
            G_fatal_error ("Coordinate column must be of type integer or double");
        value = db_get_column_value (column);
        if (ctype == DB_C_TYPE_INT)
            m_y = (double)db_get_value_int (value);
        else
            m_y = db_get_value_double (value);

        // convert the geographic coordinates to raster coordinates
        int m_row = (int) ((raster_north - m_x) / ns_res);
        int m_col = (int) ((m_y - raster_west) / ew_res);

        // save the read measurement in the output matrix
        field_meas[m_row][m_col] = m_rscp;
    }
    //
    // close the database connection
    //
    db_close_database_shutdown_driver (driver);
}



/**
 * Tests the correct functioning of dump and load field measurements
 * from a binary file. It also serves as a use case example.
 *
 * file_name    the full path of the file into which the matrix
 *              contents are saved;
 * nrows        number of rows in the 2D matrix;
 * ncols        number of cols in the 2D matrix;
 * field_meas   a 2D-matrix where the field measurements are saved.-
 */
void test_dump_and_load_from_file (const char *file_name,
                                          const int nrows,
                                          const int ncols,
                                          double **field_meas)
{
    int i, j;
    double expected = 0.0;
    dump_field_measurements (file_name,
                             nrows,
                             ncols,
                             field_meas);
    double **m_temp = (double **) G_calloc (nrows, sizeof (double *));
    for (i = 0; i < nrows; i ++)
        m_temp[i] = (double *) G_calloc (ncols, sizeof (double));
    load_field_measurements_from_file (file_name,
                                       m_temp);
    for (i = 0; i < nrows; i ++)
        for (j = 0; j < ncols; j ++)
            expected += (field_meas[i][j] - m_temp[i][j]);
    assert (expected == 0.0);                                  
    G_fatal_error ("Finished with no errors");
}

