#include <math.h>
#include <performance_metric.h>



/**
 * Returns the norm difference between field 
 * measurements and the predicted RSCP values.
 *
 * tx_power     Transmitter pilot power, in dBm;
 * nrows        number of rows in the matrix;
 * ncols        number of columns in the matrix;
 * field_meas   a 2D-matrix where the field measurements are saved;
 * prediction   a 2D-matrix containing the predicted RSCP values.-
 *
 */
double calc_norm (const double tx_power,
                  const int nrows,
                  const int ncols,
                  double **field_meas,
                  double **prediction)
{
    double ret_value = -1.0;
    double residual;
    int r, c;

    /**
     * this loops measures around 108.70 MFlops
     *
    for (r = 0; r < nrows; r ++)
    {
        for (c = 0; c < ncols; c ++)
        {
#ifdef _PERFORMANCE_METRICS_
            memory_access (2, 8);
#endif
            double field_measurement = field_meas[r][c];
            double model_prediction  = prediction[r][c];

            // compare predicted RSCP and measured RSCP
            if (!isnan (field_measurement) && field_measurement != 0)
            {
                if (!isnan (model_prediction) && model_prediction != 0)
                {
                    residual  = field_measurement;
                    residual -= tx_power - model_prediction;
                    ret_value += residual * residual;
                }
            }
        }
    }*/

    //
    // unroll the loops to improve the number of flops per memory transfer
    //
    // this unrolled version measures around 142.80 MFlops
    //
    for (r = 0; r < (nrows - 1) / 2; r ++)
    {
        for (c = 0; c < (ncols - 1) / 2; c ++)
        {
#ifdef _PERFORMANCE_METRICS_
            memory_access (8, 8);
#endif

            int i;
            double field_measurement [4];
            double model_prediction [4];

            field_measurement[0] = field_meas[2*r][2*c];
            field_measurement[1] = field_meas[2*r][2*c + 1];
            field_measurement[2] = field_meas[2*r + 1][2*c];
            field_measurement[3] = field_meas[2*r + 1][2*c + 1];

            model_prediction[0]  = prediction[2*r][2*c];
            model_prediction[1]  = prediction[2*r][2*c + 1];
            model_prediction[2]  = prediction[2*r + 1][2*c];
            model_prediction[3]  = prediction[2*r + 1][2*c + 1];

            for (i = 0; i < 4; i ++)
            {
                // compare predicted RSCP and measured RSCP
                if (!isnan (field_measurement[i]) && field_measurement[i] != 0)
                {
                    if (!isnan (model_prediction[i]) && model_prediction[i] != 0)
                    {
                        residual  = field_measurement[i];
                        residual -= tx_power - model_prediction[i];
                        ret_value += residual * residual;
                    }
                }
            }
        }
    }
    //
    // the last row, for all the columns
    //
    r = nrows - 1;
    for (c = 0; c < ncols; c ++)
    {
#ifdef _PERFORMANCE_METRICS_
        memory_access (2, 8);
#endif
        double field_measurement = field_meas[r][c];
        double model_prediction  = prediction[r][c];

        // compare predicted RSCP and measured RSCP
        if (!isnan (field_measurement) && field_measurement != 0)
        {
            if (!isnan (model_prediction) && model_prediction != 0)
            {
                residual  = field_measurement;
                residual -= tx_power - model_prediction;
                ret_value += residual * residual;
            }
        }
    }
    //
    // the last column, for all the rows
    //
    c = ncols - 1;
    for (r = 0; r < nrows; r ++)
    {
#ifdef _PERFORMANCE_METRICS_
        memory_access (2, 8);
#endif
        double field_measurement = field_meas[r][c];
        double model_prediction  = prediction[r][c];

        // compare predicted RSCP and measured RSCP
        if (!isnan (field_measurement) && field_measurement != 0)
        {
            if (!isnan (model_prediction) && model_prediction != 0)
            {
                residual  = field_measurement;
                residual -= tx_power - model_prediction;
                ret_value += residual * residual;
            }
        }
    }

    /*
    // unroll the loops to improve the number of flops per memory transfer
    //
    // this unrolled version measures around 138.90 MFlops
    //
    for (r = 0; r < nrows / 3; r ++)
    {
        for (c = 0; c < (ncols - 1) / 2; c ++)
        {
#ifdef _PERFORMANCE_METRICS_
            memory_access (12, 8);
#endif

            int i;
            double field_measurement [6];
            double model_prediction [6];

            field_measurement[0] = field_meas[3*r][2*c];
            field_measurement[1] = field_meas[3*r][2*c + 1];
            field_measurement[2] = field_meas[3*r + 1][2*c];
            field_measurement[3] = field_meas[3*r + 1][2*c + 1];
            field_measurement[4] = field_meas[3*r + 2][2*c];
            field_measurement[5] = field_meas[3*r + 2][2*c + 1];

            model_prediction[0]  = prediction[3*r][2*c];
            model_prediction[1]  = prediction[3*r][2*c + 1];
            model_prediction[2]  = prediction[3*r + 1][2*c];
            model_prediction[3]  = prediction[3*r + 1][2*c + 1];
            model_prediction[4]  = prediction[3*r + 2][2*c];
            model_prediction[5]  = prediction[3*r + 2][2*c + 1];

            for (i = 0; i < 6; i ++)
            {
                // compare predicted RSCP and measured RSCP
                if (!isnan (field_measurement[i]) && field_measurement[i] != 0)
                {
                    if (!isnan (model_prediction[i]) && model_prediction[i] != 0)
                    {
                        residual  = field_measurement[i];
                        residual -= tx_power - model_prediction[i];
                        ret_value += residual * residual;
                    }
                }
            }
        }
    }
    //
    // the last column, for all the rows
    //
    c = ncols - 1;
    for (r = 0; r < nrows; r ++)
    {
#ifdef _PERFORMANCE_METRICS_
        memory_access (2, 8);
#endif
        double field_measurement = field_meas[r][c];
        double model_prediction  = prediction[r][c];

        // compare predicted RSCP and measured RSCP
        if (!isnan (field_measurement) && field_measurement != 0)
        {
            if (!isnan (model_prediction) && model_prediction != 0)
            {
                residual  = field_measurement;
                residual -= tx_power - model_prediction;
                ret_value += residual * residual;
            }
        }
    }*/
    //
    // return the least-squares value
    //
    return ret_value;
}

