#ifndef _COVERAGE_NORM_H_
#define _COVERAGE_NORM_H_

double calc_norm (const double tx_power,
                  const int nrows,
                  const int ncols,
                  double **field_meas,
                  double **prediction);

#endif
