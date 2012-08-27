#ifndef _COVERAGE_ANTENNA_H_
#define _COVERAGE_ANTENNA_H_

void calculate_antenna_influence (const char * ant_dir,
                                  const char * ant_file,
                                  const double east,
                                  const double north,
                                  const double total_height,
                                  const int beam_dir,
                                  const int mech_tilt,
                                  const double radius,
                                  const double rec_height,
                                  const int nrows,
                                  const int ncols,
                                  const double west_ext,
                                  const double north_ext,
                                  const double ew_res,
                                  const double ns_res,
                                  double **dem,
                                  double **path_loss);

#endif
