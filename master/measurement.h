#ifndef _COVERAGE_MEASUREMENT_H_
#define _COVERAGE_MEASUREMENT_H_

/**
 * Loads the field-measurements array from a raster map.
 *
 * mapname     the full name of the raster map from which the matrix
 *              should be loaded;
 * field_meas   a 2D-matrix where the loaded data are saved.-
 */
void 
load_field_measurements_from_map (const char *mapname,
                                  double **field_meas);

void load_field_measurements_from_db (const char *tx_name,
                                      const double tx_east,
                                      const double tx_north,
                                      const double radius,
                                      const double raster_west,
                                      const double raster_north,
                                      const double ew_res,
                                      const double ns_res,
                                      double **field_meas);

void dump_field_measurements (const char *file_name,
                              const int nrows,
                              const int ncols,
                              double **field_meas);

void load_field_measurements_from_file (const char *file_name,
                                        double **field_meas);

#endif
