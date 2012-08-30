#ifndef _COVERAGE_ERICSSON9999_H_
#define _COVERAGE_ERICSSON9999_H_

#include <math.h>
#include <stdlib.h>
#include <grass/gis.h>
#include <grass/glocale.h>


//
// Parameters of the Ericsson 9999 radio propagation model
//
struct StructEric {
	double BSxIndex;    /* normalized position of BS -> UTMx/resolution */	
	double BSyIndex; 	/* normalized position of BS -> UTMx/resolution */
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


int EricPathLossSub (double** Raster, 
                     double** Clutter, 
                     double** PathLoss, 
                     struct StructEric *IniEric);

#endif

