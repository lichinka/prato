
/****************************************************************************
 *
 * MODULE:      r.ericsson
 * AUTHOR(S):   Andrej Vilhar, Tomaz Javornik, Andrej Hrovat Jozef Stefan Institute
 * UPDATE:		Patrik Ritosa, Telekom Slovenije                
 *
 * PURPOSE:     Calculates radio coverage from a single base station 
 *              according to Ericsson model
 *             
 *
 * COPYRIGHT:   (C) 2009 Jozef Stefan Institute
 *
 *
 *****************************************************************************/
/****************************************************************************
 *
 * MODULE:       do_profile function for r.ericsson
 * AUTHOR(S):    Patrik Ritosa, Telekom Slovenije                
 *
 * PURPOSE:      Calculates obstacle high and distance betwen BS and MS
 *             
 *
 * COPYRIGHT:
 *
 *
 *****************************************************************************/

#include "worker/eric.h"
#include "worker/antenna.h"



//
// 2D area matrix where the obstacle heights are saved
//
static double **Obst_high = NULL;

// 
// 2D area matrix where the obstacle distances are saved
//
static double **Obst_dist = NULL;
   
//
// 2D area matrix for some offsets of the line-of-sight calculation
//
static double **Offset = NULL;



/**
 * Calculates the terrain profile for line-of-sight.
 *
 * Obst_high    output 2D matrix of heights;
 * Obst_dist    output 2D matrix of distances;
 * Offset       output 2D matrix of offsets.-
 *
 */
static void calc_profile (double** Obst_high, 
                          double** Obst_dist, 
                          double** Raster, 
                          double** Offset, 
                          const double dx, 
                          const double dy, 
                          const double xBS, 
                          const double yBS,
                          const double ZoTransBS, 
                          const int xN, 
                          const int yN, 
                          const double scale, 
                          const double radius)
{
	
	int X,Y;
	double x_tmp, y_tmp;
	double temp_dist, max_dist_old, max_dist_los, max_dist_nlos, max_dist, cand_dist, cand_el_dist;
	double offset_dist;
	
	double max_high_old, max_high_los, max_high_nlos, max_high, cand_high, cand_el_high;
	double cand_obstacle, cand_el_obstacle, max_obstacle;
	double EL, EL_temp, EL_max;
	double EL_nlos, EL_nlos_temp, EL_nlos_max;
	
	double interp_high, koef_x, koef_y; 
	
	int LOS_start = 1;
	
	
	//Patrik x_tmp = xBS + dx;
	//Patrik y_tmp = yBS + dy;
	x_tmp = floor(xBS) + 0.5 + dx;
	y_tmp = floor(yBS) + 0.5 + dy;
	//Patrik x_tmp = floor(xBS) + dx;
	//Patrik y_tmp = floor(yBS) + dy;
		
	//Patrik X = Y = 0;
		
	max_high = max_high_old = Raster [(int)xBS][(int)yBS];
	max_dist = max_dist_old = 1;		// 1 bin
	//Patrik cand_high = -999;
	//Patrik cand_dist = 1;
	EL_max = 0;
	EL_nlos_max = 0;
	while ((int)x_tmp >= 0 && (int)x_tmp < xN && (int)y_tmp >= 0 && (int)y_tmp < yN)
	{
		X = (int)x_tmp;
		Y = (int)y_tmp;
			
		//Patrik temp_dist = sqrt(pow(floor(xBS) + 0.5 - x_tmp,2)+pow(floor(yBS) + 0.5 - y_tmp,2));
		temp_dist = sqrt(pow(floor(xBS) - floor(x_tmp),2)+pow(floor(yBS) - floor(y_tmp),2));	//zaokorzena razdalja na bin - pravilno za racun
		offset_dist = sqrt(pow(floor(xBS) + 0.5 - x_tmp,2)+pow(floor(yBS) + 0.5 - y_tmp,2));	//natancna razdalja potrebna za izracun offseta pri dolocanju bina
		
		if (temp_dist*scale > radius*1000){
			break;
		}
			
// interpolacija DEM-a za natancnejso visino
			if ((x_tmp - (floor(x_tmp) + 0.5)) <= 0){
				if ((X-1) < 0){
					koef_x = 0;	
				}
				else {
					koef_x = (Raster[X][Y] - Raster[X-1][Y]);
				}
			}
			else{
				if ((X+1) >= xN){
					koef_x = 0;	
				}
				else {
					koef_x = (Raster[X+1][Y] - Raster[X][Y]);
				}
			}
			
			if ((y_tmp - (floor(y_tmp) + 0.5)) <= 0){
				if ((Y-1) < 0){
					koef_y = 0;
				}
				else {
					koef_y = (Raster[X][Y] - Raster[X][Y-1]);
				}
			}
			else{
				if ((Y+1) >= yN){
					koef_y = 0;
				}
				else{
					koef_y = (Raster[X][Y+1] - Raster[X][Y]);
				}
			}
		
			
			interp_high = Raster[X][Y] + koef_x * (x_tmp - (floor(x_tmp) + 0.5)) + koef_y * (y_tmp - (floor(y_tmp) + 0.5));
// konec interpolacija DEM-a za natancnejso visino
			
			//Patrik EL = atan(temp_dist * scale / (ZoTransBS - Raster[X][Y]));
			EL = atan(temp_dist * scale / (ZoTransBS - interp_high));
			
			if (EL < 0){
				EL_temp = _PI_ + EL;
			}
			else{
				EL_temp = EL;
			}



// LOS case			
			if (EL_temp > EL_max){
				if (LOS_start == 1){
					max_high_los = max_high;
					max_dist_los = max_dist;
				}
				else{
					max_high_los = max_high_old;
					max_dist_los = max_dist_old;
				}	
	
				
				if (fabs(offset_dist - temp_dist) < Offset[X][Y]){
					Obst_high[X][Y] = max_high_los - (ZoTransBS - max_dist_los * scale / tan(EL));
					Obst_dist[X][Y] = max_dist_los;			//do sedaj prevladuloca ovira
					Offset[X][Y] = fabs(offset_dist - temp_dist);
				}
				EL_max = EL_temp;
					
				//Patrik if (Raster[X][Y] >= max_high){
					//Patrik max_high = cand_high = Raster[X][Y];
					max_high = cand_high = interp_high;
					max_dist = cand_dist = temp_dist;			//nova ovira

				//Patrik }
				
				
				cand_el_high = max_high;			// ponastavi za naslednjo NLOS iteracijo
				cand_el_dist = max_dist;
				EL_nlos_max = 0; // prestudiraj!!
			}
// NLOS case
			else{									// EL_temp <= EL_max
				LOS_start = 0;
				
				
				//Patrik EL_nlos = atan((temp_dist - cand_dist) * scale / (cand_high + 50 - Raster[X][Y]));
				EL_nlos = atan((temp_dist - cand_dist) * scale / (cand_high + 50 - interp_high));
				
				if (EL_nlos < 0){
					EL_nlos_temp = _PI_ + EL_nlos;
				}
				else{
					EL_nlos_temp = EL_nlos;
				}
				
				if (EL_nlos_temp > EL_nlos_max){
					EL_nlos_max = EL_nlos_temp;
				}
				else{					//EL_nlos_temp <= EL_nlos_max
					//Patrik cand_high = Raster[X][Y];
					cand_high = interp_high;
					cand_dist = temp_dist;
					EL_nlos_max = 0; 			// reset kota, ker smo dosegli novo tocko LOS v NLOS obmocju
				}
				
				cand_obstacle = cand_high - ZoTransBS + cand_dist * scale / tan(EL);		// kandidat za novo oviro NLOS
				cand_el_obstacle = cand_el_high - ZoTransBS + cand_el_dist * scale / tan(EL);		// izbran kandidat za novo oviro NLOS
				max_obstacle = max_high - ZoTransBS + max_dist * scale / tan(EL);					// prevladujoca ovira iz LOS
				
				if((cand_obstacle > max_obstacle) && (cand_obstacle > cand_el_obstacle)){
					cand_el_high = cand_high;
					cand_el_dist = cand_dist;		
				}
					
				if (cand_el_obstacle >= max_obstacle){
					max_high_nlos = cand_el_high;
					max_dist_nlos = cand_el_dist;				
				}
				else{
					max_high_nlos = max_high;
					max_dist_nlos = max_dist;					
				}	
					

				if (fabs(offset_dist - temp_dist) < Offset[X][Y]){
					Obst_high[X][Y] = max_high_nlos - (ZoTransBS - max_dist_nlos * scale / tan(EL));
					Obst_dist[X][Y] = max_dist_nlos;
					Offset[X][Y] = fabs(offset_dist - temp_dist);
				}
					
				if (max_high > max_high_nlos){
					max_high_old = max_high;
					max_dist_old = max_dist;
				}
				else{
					max_high_old = max_high_nlos;
					max_dist_old = max_dist_nlos;
				}
			}

		x_tmp = x_tmp + dx;
		y_tmp = y_tmp + dy;	
	}
}



static int DoProfile (double** Obst_high, 
                      double** Obst_dist, 
                      double** Offset, 
                      double ResDist, 
                      double** Raster, 
                      double xBS, 
                      double yBS, 
                      double ZoTransBS, 
                      int xN, 
                      int yN, 
                      double tiltBS2MS, 
                      double scale, 
                      double radius)

{
	double AZI;
	int ix, iy;
	double dx, dy;
	
    // Offset ini
	for (ix = 0; ix < xN; ix++)
    {
		for (iy = 0; iy < yN; iy++)
        {
			Offset[ix][iy]=999;
		}
	}
	
	
	
	// Kvadrant I
		
	for (ix = 0; ix < xN; ix++)
	{
		//Patrik AZI = atan((ix - xBS) / yBS);
		AZI = atan((ix - floor(xBS)) / floor(yBS));
		
		if (cos(AZI) > sin(AZI))
        {
			//Patrik dx = sin(AZI) / cos(AZI);
			//Patrik dy = -cos(AZI) / cos(AZI);
			dx = (ix - floor(xBS)) / floor(yBS);			// tan(AZI)
			dy = -1;
		}
		else
        {
			//Patrik dx = sin(AZI) / sin(AZI);
			//Patrik dy = -cos(AZI) / sin(AZI);
			dx = 1;
			dy = -floor(yBS)/(ix - floor(xBS));				// ctan(AZI)
		}
		calc_profile (Obst_high, Obst_dist, Raster, Offset, dx, dy, xBS, yBS, ZoTransBS, xN, yN, scale, radius);
	}

	// Kvadrant III

	for (ix = 0; ix < xN; ix++)
	{
		//Patrik AZI = atan((ix - xBS) / (yN - yBS));
		AZI = atan((ix - floor(xBS)) / (yN - floor(yBS)));
		
		if (cos(AZI) > sin(AZI))
        {
			//Patrik dx = sin(AZI) / cos(AZI);
			//Patrik dy = cos(AZI) / cos(AZI);
			dx = (ix - floor(xBS)) / (yN - floor(yBS));			// tan(AZI)
			dy = 1;
		}
		else
        {
			//Patrik dx = sin(AZI) / sin(AZI);
			//Patrik dy = cos(AZI) / sin(AZI);
			dx = 1;
			dy = (yN - floor(yBS)) / (ix - floor(xBS));				// ctan(AZI)
		}
				
		calc_profile (Obst_high, Obst_dist, Raster, Offset, dx, dy, xBS, yBS, ZoTransBS, xN, yN, scale, radius);
	} 
	
	// Kvadrant II

	for (iy = 0; iy < yN; iy++)
	{
		//Patrik AZI = atan((iy - yBS) / (xN - xBS));
		AZI = atan((iy - floor(yBS)) / (xN - floor(xBS)));
			
		if (cos(AZI) > sin(AZI))
        {
			//Patrik dx = cos(AZI) / cos(AZI);
			//Patrik dy = sin(AZI) / cos(AZI);
			dx = 1;			
			dy = (iy - floor(yBS)) / (xN - floor(xBS));		// tan(AZI)
		}
		else
        {
			//Patrik dx = cos(AZI) / sin(AZI);
			//Patrik dy = sin(AZI) / sin(AZI);
			dx = (xN - floor(xBS)) / (iy - floor(yBS));				// ctan(AZI)
			dy = 1;
		}
		
		calc_profile (Obst_high, Obst_dist, Raster, Offset, dx, dy, xBS, yBS, ZoTransBS, xN, yN, scale, radius);	
	}

	// Kvadrant IV

	for (iy = 0; iy < yN; iy++)
	{
		//Patrik AZI = atan((iy - yBS) / xBS);
		AZI = atan((iy - floor(yBS)) / floor(xBS));
		
		if (cos(AZI) > sin(AZI))
        {
			//Patrik dx = -cos(AZI) / cos(AZI);
			//Patrik dy = sin(AZI) / cos(AZI);
			dx = -1;			
			dy = (iy - floor(yBS)) / floor(xBS);		// tan(AZI)
		}
		else
        {
			//Patrik dx = -cos(AZI) / sin(AZI);
			//Patrik dy = sin(AZI) / sin(AZI);
			dx = -floor(xBS) / (iy - floor(yBS));				// ctan(AZI)
			dy = 1;
		}
		
		calc_profile (Obst_high, Obst_dist, Raster, Offset, dx, dy, xBS, yBS, ZoTransBS, xN, yN, scale, radius);			
	}
	
	return 0;
}		// end doProfile






//int EricPathLossSub(double** Raster, double** Clutter, double** PathLoss, double** Obst_high, double** Obst_dist, double** Offset, struct StructEric *IniEric)
int EricPathLossSub (double** Raster, double** Clutter, double** PathLoss, struct StructEric *IniEric)

/*************************************************************************************************
 *
 *		Function EricPathLossSub calculates PathLoss in dB using Ericsson 9999 path loss formula
 *			
 *			**PathLoss:	array of path loss in dB
 *			**Raster:       input DEM file
 *       	**Clutter:      input clutter file
 *			
 *			T.Javornik, Jan. 2010
 * 
 * 		Ericsson 9999 model update:
 * 			- Spherical earth diffraction komponent (JDFR)
 * 			- NLOS komponent (KDFR) 
 * 			- do_profile calculation speed optimisation
 * 
 * 			Patrik Ritosa, May 2012
 *
 *************************************************************************************************/
{
    int i;

	// Ericsson model constants and variables
	int BSxIndex = IniEric->BSxIndex;		    //	normalized position of BS -> UTMx/resolution 
	int BSyIndex = IniEric->BSyIndex;		    //	normalized position of BS -> UTMy/resolution
	double AntHeightBS = IniEric->BSAntHeight;	//	Antenna height of BS [m]
	double AntHeightMS = IniEric->MSAntHeight;	//	Antenna height of MS [m]
	int xN = IniEric->xN;				//	dimension of the input(Raster) and output (PathLoss)
	int yN = IniEric->yN;				//	dimension of the input(Raster) and output (PathLoss)
	double scale = IniEric->scale;			//	Resolution Erricson model
	double A0 = IniEric->A0;				//	Model 9999 parameters
	double A1 = IniEric->A1;				//	Model 9999 parameters
	double A2  = IniEric->A2;			//	Model 9999 parameters
	double A3  = IniEric->A3;			//	Model 9999 parameters
	double freq  = IniEric->freq;			//	carrier frequency
	double ResDist = IniEric->ResDist;		//	distance BS - MS sampling rate [normalized with scale
	double Lambda = 300.0/freq;			//	wave lenght
	double radi = IniEric->radi;			// radius of calculation

	double ZoBS;					// BS and MS height about the sea level
	double ZObs2LOS = 0;
	double DistObs2BS = 0;
	double ZoTransBS,ZoTransMS;
	double log10Zeff;
	double log10DistBS2MSKm;
	double tiltBS2MS;				// (ZoBS-ZoMS)/distBS2MSNorm	
	double PathLossFreq = 0;			// path loss due to carrier frequency
	double PathLossTmp = 0;				// tmp path loss
	int ix; int iy;	
	int DiffX, DiffY;
    double Zeff;			// Difference in X and Y direction
	double PathLossAntHeightBS;
	double DistBS2MSNorm, DistBS2MSKm;		// distance between MS and BS in Km sqrt(x2+y2+z2) * scale / 1000
							// normalized distance between MS and BS in xy plan sqrt(x2+y2)
	double ElevAngCos, Hdot, Ddot, Ddotdot, PathLossDiff, KDFR, Alfa, Fresnel, JDFR;
	
	ZoBS = Raster[(int)BSxIndex][(int)BSyIndex];	// BS height above the sea level calculated from raster DEM file
	ZoTransBS = ZoBS + AntHeightBS;			// BS transmitter height above the sea level

	PathLossFreq = 44.49*log10(freq) - 4.78*pow(log10(freq),2);	// Loss due to carrier frequency
									/*POPRAVJNEO (4.2.2010)*/	
	PathLossAntHeightBS = 3.2*pow(log10(11.75*AntHeightMS),2);


    //
    // LOS and obstacle height calculation is executed only once, 
    // because its results are constant throughout the optimization
    //
    if ((Obst_high == NULL) && (Obst_dist == NULL) && (Offset == NULL))
    {
        // Obstacle height
        Obst_high = (double **) calloc (xN, sizeof(double *));
        for (i = 0; i < xN; i ++)
            Obst_high[i] = (double *) calloc (yN, sizeof (double));

        // Obstacle distance 
        Obst_dist = (double **) calloc (xN, sizeof(double *));
        for (i = 0; i < xN; i ++)
            Obst_dist[i] = (double *) calloc (yN, sizeof (double));
    
        // do_profile offset 
        Offset = (double **) calloc (xN, sizeof(double *));
        for (i = 0; i < xN; i ++)
            Offset[i] = (double *) calloc (yN, sizeof (double));

        // calculate the terrain profile
        DoProfile (Obst_high,Obst_dist,Offset,ResDist,Raster,BSxIndex,BSyIndex,ZoTransBS,xN,yN,tiltBS2MS,scale, radi);
    }

	for (ix = 0; ix < xN; ix++)
	{
		for (iy = 0; iy < yN; iy++)
		{
			// Path Loss due to Hata model
			DiffX = (BSxIndex-ix); DiffY = (BSyIndex-iy);
			// ZoMS = Raster[ix][iy];
			ZoTransMS = Raster[ix][iy]+AntHeightMS;  // ZoMS
			Zeff = ZoTransBS - ZoTransMS;		// ??
			DistBS2MSKm = sqrt (DiffX*DiffX + DiffY*DiffY)*scale/1000; //sqrt(DiffX*DiffX+DiffY*DiffY+Zeff*Zeff)*scale/1000;			
			DistBS2MSNorm = sqrt(DiffX*DiffX+DiffY*DiffY);
//			if(ZoBS <= Raster[ix][iy]) 
//			{
//				Zeff = AntHeightBS;  // ZoMS
//			}
			if (DistBS2MSKm < 0.01)
			{
				DistBS2MSKm = 0.01;
			}

			if (DistBS2MSKm > radi)
		    {    
			      	continue;
		    }

			Zeff = Zeff + (DistBS2MSKm*DistBS2MSKm)/(2 * 4/3 * 6370)*1000; //height correction due to earth sphere
			
			if (- AntHeightBS < Zeff && Zeff < AntHeightBS){
				Zeff = AntHeightBS;		// Preventing Log10(Zeff) to go toward -inf
			}
			
			log10Zeff=log10(abs(Zeff));

			//log10DistBS2MSKm=log10(sqrt(DistBS2MSKm*DistBS2MSKm + Zeff/1000 * Zeff/1000));
			log10DistBS2MSKm=log10(DistBS2MSKm);			

			PathLossTmp = A0 + A1*log10DistBS2MSKm; 

			PathLossTmp = PathLossTmp + A2*log10Zeff + A3*log10DistBS2MSKm*log10Zeff;
			PathLossTmp = PathLossTmp - PathLossAntHeightBS + PathLossFreq;

			// Calc position of the height and position of the highest obstacle

		
			tiltBS2MS = ZoTransBS - ZoTransMS; 	//STARO: tiltBS2MS = Zeff; Zeff je vmes lahko spremenjena /* Sprememba (4.2.2010)*/

			if (DistBS2MSNorm > 0) 
            {
				tiltBS2MS = -tiltBS2MS/DistBS2MSNorm; 
            }
			else 
            {
				tiltBS2MS = 0; 
			}
			
			//DoProfile(&ZObs2LOS,&DistObs2BS,ResDist,Raster,BSxIndex,BSyIndex,ZoTransBS,ix,iy,tiltBS2MS);
			ZObs2LOS = Obst_high[ix][iy];
			DistObs2BS = Obst_dist[ix][iy];
			// Calc path loss due to NLOS conditions
/*Patrik/scale*/	
            ElevAngCos = cos(atan(tiltBS2MS/scale));

			Ddot = DistObs2BS; 
			if (ElevAngCos != 0) 
			{	
				Ddot = DistObs2BS/ElevAngCos;
			}
			Ddotdot = DistBS2MSNorm - Ddot;
			if (ElevAngCos != 0) {
				Ddotdot = (DistBS2MSNorm)/ElevAngCos - Ddot;
			}

//Obstacle height korrection due to earth sphere
			if (Ddot <= Ddotdot){
				ZObs2LOS = ZObs2LOS + (Ddot*scale/1000*Ddot*scale/1000)/(2 * 4/3 * 6370)*1000;
			}
			else{
				ZObs2LOS = ZObs2LOS + (Ddotdot*scale/1000*Ddotdot*scale/1000)/(2 * 4/3 * 6370)*1000;
			}

//Hight correction due to BS2MS line angle
			Hdot = ZObs2LOS*ElevAngCos;

			PathLossDiff = 0;
			KDFR = 0;
			if (Ddot > 0 && Ddotdot > 0) {
				Fresnel=sqrt((Lambda*Ddot*Ddotdot*scale)/(Ddot+Ddotdot)); // First Fresnel elipsoid radius

				PathLossDiff = Hdot/Fresnel;

// NLOS komponent calculation KDFR

				if (PathLossDiff < -0.49 ) {
					KDFR = 0; 
				}
				else if(-0.49 <= PathLossDiff && PathLossDiff < 0.5) {
					KDFR = 6 + 12.2 * PathLossDiff;
				}
				else if(0.5 <= PathLossDiff && PathLossDiff < 2) {
					KDFR = 11.61 * PathLossDiff - 2 * PathLossDiff *PathLossDiff + 6.8;
				}
				else if(2 <= PathLossDiff) {
					KDFR = 16 + 20 * log10(PathLossDiff);
				}

// Alfa correction factor
				
				if (Hdot > Fresnel/2){
					Alfa = 1;
				}
				else if (Fresnel/4 <= Hdot && Hdot <= Fresnel/2){
					Alfa = 4 * Hdot/Fresnel - 1;
				}
				else if (Hdot < Fresnel/4){
					Alfa = 0;
				}
			
			}

//Spherical earth diffraction komponent JDFR
			if (Hdot > 0){
				JDFR = 20 + 0.112 * pow(freq/(16/9),1/3) * (DistBS2MSKm - sqrt(12.73 * 4/3) * ( sqrt(ZoTransBS) + sqrt(ZoTransMS) )  );
			
				if (JDFR < 0){
					JDFR=0;
				}
			}
			else{
				JDFR = 0;
			}


			PathLossTmp = PathLossTmp + sqrt(pow(Alfa*KDFR,2) + pow(JDFR,2));		
			
			// write data to pathloss
			
			PathLoss[ix][iy] = PathLossTmp + Clutter[ix][iy];
		} // end irow
	} // end icol
	return 0;
}

