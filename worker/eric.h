#ifndef _COVERAGE_ERICSSON9999_H_
#define _COVERAGE_ERICSSON9999_H_

#include "worker/gpu.h"
#include "worker/coverage.h"



//
// Parameters of the Ericsson 9999 radio propagation model
//
struct StructEric {
	int BSxIndex;       // normalized position of BS -> UTMx/resolution	
	int BSyIndex; 	    // normalized position of BS -> UTMx/resolution
	double BSAntHeight; /* Antenna height of BS [m] */
	double MSAntHeight; /* Antenna height of MS [m] */
	int xN;			    /* dimension of teh input(Raster) and output (PathLoss) */
	int yN;			    /* dimension of teh input(Raster) and output (PathLoss) */
	double scale;		/* Resolution of DEM file */
	double freq;		/* Carrier frequency in MHz */	
	double A0;		    /* Model 9999 parameters */
	double A1;		    /* Model 9999 parameters */
	double A2;		    /* Model 9999 parameters */
	double A3;		    /* Model 9999 parameters */
	double ResDist;		/* Resolution Erricsson model profile calc */
	double radi;		/* Radius of calculation [km] */	
};


/**
 * Calculates the path loss using E/// 9999 model implementation on GPU.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters;
 *
 */
void 
eric_pathloss_on_gpu (Parameters    *params,
                      Tx_parameters *tx_params);

/**
 * Calculates path loss using the E/// 9999 prediction model on CPU.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters.-
 *
 */
void 
eric_pathloss_on_cpu (Parameters    *params,
                      Tx_parameters *tx_params);

/**
 * Fine tunes the A0...A3 E/// parameters to best fit a set of measurements 
 * within a given radio zone. The result is saved in the 
 * `tx_params->eric_params` vector.
 *
 * params           a structure holding configuration parameters which are 
 *                  common to all transmitters;
 * tx_params        a structure holding transmitter-specific configuration
 *                  parameters.-
 *
 */
void
parameter_fine_tuning (Parameters    *params,
                       Tx_parameters *tx_params);

#endif

