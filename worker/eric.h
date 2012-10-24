#ifndef _COVERAGE_ERICSSON9999_H_
#define _COVERAGE_ERICSSON9999_H_

#include <math.h>
#include <stdlib.h>


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


void 
eric_pathloss_on_gpu (const double tx_east_coord,
                      const double tx_north_coord,
                      const int tx_east_idx,
                      const int tx_north_idx,
                      const double antenna_height_AGL,
                      const double total_tx_height,
                      const int beam_direction,
                      const int mechanical_tilt,
                      const double frequency,
                      const double radius,  
                      const double rx_height_AGL,
                      const int nrows,       
                      const int ncols,      
                      const double map_west,
                      const double map_north,
                      const double map_ew_res,  
                      const double map_ns_res,  
                      const float  null_value,
                      double **m_dem,          
                      double **m_clut,
                      double **m_loss);

int 
EricPathLossSub (double **Raster, 
                 double **Clutter, 
                 double **PathLoss, 
                 struct StructEric *IniEric);

#endif

