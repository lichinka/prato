#if defined(cl_amd_fp64) || defined(cl_khr_fp64)

#include "worker/constants.h"



/**
 * Single or double precision arithmetics depend on hardware support
 *
 * Single precision definition
 *
typedef float 	  real;
typedef float2	  real2;
typedef float4    real4;
 *
 */
    
/**
 * Double precision arithmetics
 */
#if defined(cl_amd_fp64)
    #pragma OPENCL EXTENSION cl_amd_fp64 : enable
#elif defined(cl_khr_fp64)
    #pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif

typedef double	real;
typedef double2	real2;
typedef double4 real4;


/**
 * A simple random number number generator as defined in Benedičič et. al. !ref!.-
 */
typedef struct { ulong a, b, c; } random_state;

unsigned long random(random_state *r)
{
    unsigned long old = r->b;
    r->b = r->a * 1103515245 + 12345;
    r->a = (~old ^ (r->b >> 3)) - r->c++;
    return r->b;
}

real random_01(random_state *r)
{
    return (random(r) & 4294967295) / 4294967295.0f;
}

void seed_random(random_state *r, ulong seed)
{
    r->a = seed;
    r->b = 0;
    r->c = 362436;
}



/**
 * Calculates a path-loss value for one thread, using Hata's model
 * for urban areas.-
 *
 * int pixel_res ........ size of one raster pixel (in meters).
 * real raster_north ... northern limit of the raster area (in meters).
 * real raster_west .... western limit of the raster area (in meters).
 * int ncols ............ total number of columns in the area.
 * real4 tx_data ....... coordinates, elevation and anthenna height of 
 *                        the transmitter (all values expressed in meters).
 * int2 tx_offset ....... offset of the current transmitter within the
 *                        area, including the calculation radius (in raster
 *                        pixels).
 * real rx_height ...... height of the receiver above ground (in mts).
 * real frequency ...... transmitter frequency in Mhz.
 * real ahr ............ a constant, part of Hata's formula.
 * __gl float dem_in ... digital elevation model data for the area.
 * __lo real2 pblock ... pointer to local memory, used to speed up the
 *                        whole process. Intermediate numbers, like distance
 *                        and height, are saved here.
 */
real path_loss_hata_urban (const int pixel_res,
                            const real raster_north,
                            const real raster_west,
                            const int ncols, 
                            const real4 tx_data, const int2 tx_offset,
                            const real rx_height, const real frequency, 
                            const real ahr,
                            __global real *dem_in,
                            __local real2 *pblock)
{
    real dist, height_diff;

    // receiver coordinates in pixels 
    int2 px_coord = (int2) ((int)tx_offset.x + get_global_id(0),
                            (int)tx_offset.y + get_global_id(1));
    
    // receiver coordinates in meters
    real2 rec_coord = (real2) ((real) px_coord.x * pixel_res + raster_west,
                                 (real) raster_north - px_coord.y * pixel_res); 
    real2 tx_coord = (real2) ((real) tx_data.x,
                                (real) tx_data.y);
    
    // calculate the distance between transmitter and receiver
    dist = distance (tx_coord, rec_coord);
    
    // correct the distance if it is smaller than 10m
    if (dist < 10.0)
        dist = pixel_res / 2.0;
   
    // element index within the digital elevation model
    uint dem_idx = px_coord.y * ncols + px_coord.x;

    // transform distance to km 
    dist /= 1000.0;
   
    // effective height between tx and rx
    if (tx_data.z > (real)dem_in[dem_idx])
    {
        height_diff = (tx_data.z + tx_data.w) - (real)dem_in[dem_idx] - rx_height;
    }
    else
    {
        height_diff = tx_data.w;
    }
    // cache the distance and height difference
    uint local_idx = get_local_id(1) * get_local_size(0) + 
                     get_local_id(0);
    pblock[local_idx] = (real2) (log10 (dist),
                                 log10 (height_diff));

    // wait for other work items to cache their values
    barrier (CLK_LOCAL_MEM_FENCE);

    // urban path-loss value in dB
    return (real) (69.55 + 26.16 * log10(frequency) -
                   13.82 * pblock[local_idx].y - ahr +
                   (44.9 - 6.55 * pblock[local_idx].y) *
                   pblock[local_idx].x);
}



/**
 * Calculates the coverage in the current area, by searching/analyzing
 * 3D cubes of path loss matrices and interference. The kernel is 
 * limited by the amount of memory on the GPU. It uses a scattered 
 * approach (by Hwu), i.e. output aligned, and local memory to reach
 * as many GFlops as possible.
 *
 * int ntx .................. number of transmitters being processed. This
 *                            is the length of the 'tx_pwr' vector and the
 *                            depth of the 'pl_in' cube.
 * int ncols ................ total number of columns in the area.
 * __gl u char pl_in ........ 3D matrix containing path-loss values in dB
 *                            (i.e. between 0 and 255) for every coordinate
 *                            (2D) and each of the 'ntx' transmitters (3rd
 *                            dimension) in the area.
 * __gl int4 offsets_in ..... vector of tuples containing the offsets for
 *                            each of the path-loss matrices (col-min, 
 *                            row-min, width, height).
 * __gl real qrm_in ........ 2D matrix containing interference values for
 *                            every coordinate in the area.
 * __gl int tx_pwr_in ....... pilot powers for 'ntx' transmitters (in 
 *                            milliwatts).
 * __gl real cov_out ....... coverage matrix where the results are saved.
 * __gl short cov_ctr_out ... coverage count vector, containing the number of
 *                            uncovered raster cells.
 * __gl short cov_coord_out . a vector containing coordinates of uncovered 
 *                            raster cells. Only some, randomly-selected,
 *                            coordinates are saved here.
 * int lmem_ctr_offset ...... the beggining of the uncovered counters within
 *                            local memory.
 * int lmem_coord_offset .... the beggining of the uncovered coordinates 
 *                            within local memory.
 * __lo real pblock ........ local memory used to speed up the whole 
 *                            process. The best received energy at every 
 *                            coordinate is saved here. And flags to count
 *                            the number of uncovered coordinates.
 *
__kernel void coverage_kern (const int ntx,
                             const int ncols,
                             __global unsigned char *pl_in,
                             __global int4 *offsets_in,
                             __global real *qrm_in,
                             __global int *tx_pwr,
                             __global real *cov_out,
                             __global unsigned short *cov_ctr_out,
                             __global unsigned short *cov_coord_out,
                             const int lmem_ctr_offset,
                             const int lmem_coord_offset,
                             __local real *pblock)
{
    int tx;
    int2 ccoord = (int2) (get_global_id(0), get_global_id(1));
    int idx_2d = ccoord.y * ncols + ccoord.x;

    // we will calculate the maximum received energy at each coordinate
    uint idx_local = get_local_id(1) * get_local_size(0) + get_local_id(0);
    // initialize the gamma value for coverage
    pblock[idx_local] = (real) NAN;
    // initialize the uncovered counter
    pblock[lmem_ctr_offset + idx_local] = 0.0f;

    // first transmitter
    // make sure we have path loss data here ...
    if ((ccoord.x >= offsets_in[0].x) &&
       ((ccoord.x - offsets_in[0].x) < offsets_in[0].z))
    {
        if ((ccoord.y >= offsets_in[0].y) &&
           ((ccoord.y - offsets_in[0].y) < offsets_in[0].w))
        {
            int idx_pl = (ccoord.y - offsets_in[0].y) * offsets_in[0].z;
            idx_pl += ccoord.x - offsets_in[0].x;
            
            real pl_value = (real) pl_in[idx_pl];

            // transform the path loss to a linear scale
            if (pl_value > _HATA_MAX_PATH_LOSS_)
                pl_value = _HATA_MAX_PATH_LOSS_;
        
            pl_value = (real) 1.0f - (pl_value /
                       (_HATA_MAX_PATH_LOSS_ - _HATA_MIN_PATH_LOSS_));
          
            // calculate gamma at this coordinate for this transmitter
            pblock[idx_local] = ((tx_pwr[0]/1000.0f) * pl_value) / qrm_in[idx_2d];
        }
    }

    // for every other transmitter
    for (tx = 1; tx < ntx; tx++)
    {
        // make sure we have path loss data here ...
        if ((ccoord.x >= offsets_in[tx].x) &&
           ((ccoord.x - offsets_in[tx].x) < offsets_in[tx].z))
        {
            if ((ccoord.y >= offsets_in[tx].y) &&
               ((ccoord.y - offsets_in[tx].y) < offsets_in[tx].w))
            {
                int idx_pl = offsets_in[tx].z * offsets_in[tx].w * tx;
                idx_pl += (ccoord.y - offsets_in[tx].y) * offsets_in[tx].w;
                idx_pl += ccoord.x - offsets_in[tx].x;
                
                real pl_value = (real) pl_in[idx_pl];

                // transform the path loss to a linear scale
                if (pl_value > _HATA_MAX_PATH_LOSS_)
                    pl_value = _HATA_MAX_PATH_LOSS_;
            
                pl_value = (real) 1.0f - (pl_value /
                           (_HATA_MAX_PATH_LOSS_ - _HATA_MIN_PATH_LOSS_));
              
                // calculate gamma at this coordinate for this transmitter
                pl_value = ((tx_pwr[tx]/1000.0f) * pl_value) / qrm_in[idx_2d];
                
                // keep only the maximum value
                if (pl_value > pblock[idx_local])
                    pblock[idx_local] = pl_value;
            }
        }
    }
    // count uncovered coordinates
    if (pblock[idx_local] < _AG_MINIMUM_GAMMA_COVERAGE_)
    {
        pblock[lmem_ctr_offset + idx_local] = (real) ccoord.x;
        pblock[lmem_coord_offset + idx_local] = (real) ccoord.y;
    }
    // wait for other threads to finish
    barrier (CLK_LOCAL_MEM_FENCE);

    // save all the results at once (coalesced write!)
    cov_out[idx_2d] = pblock[idx_local];

    // save the reduction sum of uncovered pixels
    if (idx_local == 0)
    {
        int cov_ctr = 0;
        int2 cov_coord = (int2) {0, 0};
        int idx_reduced = get_global_id(1) / _HATA_GPU_WITEM_PER_DIM_;
        idx_reduced *= get_global_size(0) / _HATA_GPU_WITEM_PER_DIM_;
        idx_reduced += get_global_id(0) / _HATA_GPU_WITEM_PER_DIM_;
        int local_size = get_local_size(0) * get_local_size(1);

        for (tx = 0; tx < local_size; tx ++)
        {
            if (pblock[lmem_ctr_offset + tx] != 0.0f)
            {
                cov_ctr ++;
                cov_coord = (int2) {(int) pblock[lmem_ctr_offset + tx],
                                    (int) pblock[lmem_coord_offset + tx]};
            }
        }
        cov_ctr_out[idx_reduced] = cov_ctr;
        cov_coord_out[2*idx_reduced] = cov_coord.x;
        cov_coord_out[2*idx_reduced + 1] = cov_coord.y;
    }
}
*/



/**
 * Calculates the interference matrix for an area by accumulating the 
 * path loss for each transmitter using Hata's model for urban areas. 
 * All coordinates are expressed in raster cells.
 *
 * int pixel_res ........ size of one raster pixel (in meters).
 * real raster_north ... northern limit of the raster area (in meters).
 * real raster_west .... western limit of the raster area (in meters).
 * int ncols ............ total number of columns in the area.
 * real4 tx_data ....... coordinates, elevation and anthenna height of 
 *                        the transmitter (coordinates expressed in raster
 *                        cells).
 * int2 tx_offset ....... offset of the current transmitter within the
 *                        area, including the calculation radius (in raster
 *                        cells).
 * real rx_height ...... height of the receiver above ground (in mts).
 * real frequency ...... transmitter frequency in Mhz.
 * real ahr ............ a constant part of Hata's formula.
 * __gl float dem_in ... digital elevation model data for the area.
 * __gl float pl_out ... path loss matrix where the results are saved.
 * __lo real2 pblock ... local memory used to speed up the whole process.
 *                        Intermediate numbers, like distance and height,
 *                        are saved here.
 *
__kernel void hata_urban_interference (const int pixel_res,
                                       const real raster_north,
                                       const real raster_west,
                                       const int ncols, 
                                       const real4 tx_data, const int2 tx_offset,
                                       const real rx_height, const real frequency, 
                                       const real ahr,
                                       __global real *dem_in,
                                       __global real *qrm_out,
                                       __local real2 *pblock)
{
    // calculate the path-loss value on the current location
    real pl_db = path_loss_hata_urban (pixel_res,
                                        raster_north,
                                        raster_west,
                                        ncols,
                                        tx_data, tx_offset,
                                        rx_height, frequency, 
                                        ahr,
                                        dem_in,
                                        pblock);
    // transform the path-loss to a linear scale
    if (pl_db > _HATA_MAX_PATH_LOSS_)
        pl_db = _HATA_MAX_PATH_LOSS_;

    pl_db = (real) 1.0f - (pl_db / 
            (_HATA_MAX_PATH_LOSS_ - _HATA_MIN_PATH_LOSS_));

    // calculate the interference on this point,
    // by accumulating the combination of cell power
    // and path-loss value
    uint element_idx = get_global_id(1) * get_global_size(0) + 
                       get_global_id(0);
    qrm_out[element_idx] += pl_db * _HATA_MAX_CELL_POWER_;
}
*/



/**
 * Calculates the path-loss matrix for one transmitter using 
 * Ericsson 9999 model. All coordinates are expressed in raster cells.
 *
 * int pixel_res ......... size of one raster pixel (in meters).
 * real raster_north .... northern limit of the raster area (in meters).
 * real raster_west ..... western limit of the raster area (in meters).
 * int ncols ............ total number of columns in the area.
 * real4 tx_data ........ coordinates, elevation and anthenna height of 
 *                        the transmitter (coordinates expressed in raster
 *                        cells).
 * int2 tx_offset ....... offset of the current transmitter within the
 *                        area, including the calculation radius (in raster
 *                        cells).
 * real rx_height ...... height of the receiver above ground (in mts).
 * real frequency ...... transmitter frequency in Mhz.
 * real ahr ............ a constant part of Hata's formula.
 * __gl float dem_in ... digital elevation model data for the area.
 * __gl float pl_out ... path loss matrix where the results are saved.
 * __lo real2 pblock ... local memory used to speed up the whole process.
 *                        Intermediate numbers, like distance and height,
 *                        are saved here.
 */
__kernel void eric_per_tx (const real pixel_res,
                           const real raster_north,
                           const real raster_west,
                           const real4 tx_data, 
                           const int2 tile_offset,
                           const real rx_height,
                           const real frequency,
                           const real4 ericsson_params,
                           __global real *obst_height_in,
                           __global real *obst_dist_in,
                           __global real *dem_in,
                           __global real *clut_cat_in,
                           __global real *clut_loss_in,
                           __global real *pl_out,
                           __local  real *pblock)
{
    uint element_idx = get_global_id(1) * get_global_size(0) + 
                       get_global_id(0);
    uint local_idx   = get_local_id(1) * get_local_size(0) + 
                       get_local_id(0);

    real Lambda = 300.0 / frequency;			//	wave lenght
    real tiltBS2MS;
    real ElevAngCos, Hdot, Ddot, Ddotdot, PathLossDiff, KDFR, Alfa, Fresnel, JDFR;
    
    // receiver coordinates in pixels 
    int2 rx_coord_in_pixels  = (int2) ((int)tile_offset.x + get_global_id (0),
                                       (int)tile_offset.y + get_global_id (1));
    
    // receiver coordinates in meters
    real2 rx_coord_in_meters = (real2) ((real) rx_coord_in_pixels.x * pixel_res + raster_west,
                               (real) raster_north - rx_coord_in_pixels.y * pixel_res); 
    real2 tx_coord_in_meters = (real2) ((real) tx_data.x,
                                        (real) tx_data.y);
    //
    // access the needed data from global memory
    //
    real ZObs2LOS            = obst_height_in[element_idx];
    real DistObs2BS          = obst_dist_in[element_idx];
    real rx_corrected_height = dem_in[element_idx];

    //
    // cache the partial path-loss value
    //
    int clutter_category = (int) clut_cat_in[element_idx];
    pblock[local_idx]    = (real) clut_loss_in[clutter_category];

    // wait for other work items to cache their values
    barrier (CLK_LOCAL_MEM_FENCE);

    rx_corrected_height += rx_height;

    //
    // calculate the distance between the Tx and Rx in kilometers
    //
    real dist_tx_rx_km;
    dist_tx_rx_km  = distance (tx_coord_in_meters, 
                               rx_coord_in_meters);
    dist_tx_rx_km /= 1000;

    //
    // calculate the distance between the Tx and Rx in pixels
    //
    real dist_tx_rx_px;
    dist_tx_rx_px = (tx_data.x - rx_coord_in_meters.x) *
                    (tx_data.x - rx_coord_in_meters.x);
    dist_tx_rx_px += (tx_data.y - rx_coord_in_meters.y) *
                     (tx_data.y - rx_coord_in_meters.y);
    dist_tx_rx_px = sqrt (dist_tx_rx_px) / pixel_res;

    //
    // prevent prediction of receiver points that are too near the antenna
    //
    if (dist_tx_rx_km < 0.01)
    {
        dist_tx_rx_km = 0.01;
        dist_tx_rx_px = round ((dist_tx_rx_km * 1000) / pixel_res);
    }

    //
    // HOA, HEBK calculation
    //
    real HOA, HEBK;
    
    //
    // effective antenna height
    //
    HEBK  = tx_data.z - rx_corrected_height;

    //
    // height correction due to earth curvature
    //
    HEBK += (dist_tx_rx_km * dist_tx_rx_km) / ((6370 * 8000) / 3);

    //
    // prevent the influence of this term exceeding the limit 
    // of the prediction model; `tx_data.w` is the height of the 
    // antenna above ground level
    //
    if ((HEBK < tx_data.w) && (HEBK > -tx_data.w))
        HEBK = tx_data.w;
    else
        HEBK = fabs (HEBK);

    HOA  = ericsson_params.x;
    HOA += ericsson_params.y * log10 (dist_tx_rx_km);
    HOA += ericsson_params.z * log10 (HEBK);
    HOA += ericsson_params.w * log10 (dist_tx_rx_km) * log10 (HEBK);
    HOA -=   3.2 * log10 (11.75 * rx_height) * log10 (11.75 * rx_height);
    HOA += 44.49 * log10 (frequency);
    HOA -=  4.78 * log10 (frequency) * log10 (frequency);

    // partial path-loss value
    pblock[local_idx] += HOA;

    //
    // calculate the height and position of the highest obstacle
    //
    tiltBS2MS = tx_data.z - rx_corrected_height;

    if (dist_tx_rx_px > 0) 
        tiltBS2MS = -tiltBS2MS / dist_tx_rx_px;
    else 
        tiltBS2MS = 0; 
    

    // Calculate path loss due to NLOS conditions
    ElevAngCos = cos ( atan (tiltBS2MS / pixel_res) );

    Ddot = DistObs2BS; 
    if (ElevAngCos != 0) 
        Ddot = DistObs2BS/ElevAngCos;

    Ddotdot = dist_tx_rx_px - Ddot;
    if (ElevAngCos != 0)
        Ddotdot = dist_tx_rx_px / ElevAngCos - Ddot;

    // Obstacle height korrection due to earth sphere
    if (Ddot <= Ddotdot)
        ZObs2LOS = ZObs2LOS + (Ddot*pixel_res/1000*Ddot*pixel_res/1000)/(2 * 4/3 * 6370)*1000;
    else
        ZObs2LOS = ZObs2LOS + (Ddotdot*pixel_res/1000*Ddotdot*pixel_res/1000)/(2 * 4/3 * 6370)*1000;

    // Height correction due to BS2MS line angle
    Hdot = ZObs2LOS*ElevAngCos;

    PathLossDiff = 0;
    KDFR = 0;
    if (Ddot > 0 && Ddotdot > 0)
    { 
        // First Fresnel elipsoid radius
        Fresnel = sqrt ( (Lambda*Ddot*Ddotdot*pixel_res) / (Ddot+Ddotdot) );

        PathLossDiff = Hdot/Fresnel;

        // NLOS komponent calculation KDFR
        if (PathLossDiff < -0.49 )
            KDFR = 0; 
        else if (-0.49 <= PathLossDiff && PathLossDiff < 0.5)
            KDFR = 6 + 12.2 * PathLossDiff;
        else if (0.5 <= PathLossDiff && PathLossDiff < 2)
            KDFR = 11.61 * PathLossDiff - 2 * PathLossDiff *PathLossDiff + 6.8;
        else if (2 <= PathLossDiff)
            KDFR = 16 + 20 * log10(PathLossDiff);
       
        // Alfa correction factor
        if (Hdot > Fresnel/2)
            Alfa = 1;
        else if (Fresnel/4 <= Hdot && Hdot <= Fresnel/2)
            Alfa = 4 * Hdot/Fresnel - 1;
        else if (Hdot < Fresnel/4)
            Alfa = 0;
    }

    // Spherical earth diffraction komponent JDFR
    if (Hdot > 0)
    {
        JDFR = 20 + 0.112 * pow(frequency/(16/9),1/3) * (dist_tx_rx_km - sqrt(12.73 * 4/3) * ( sqrt(tx_data.z) + sqrt(rx_corrected_height) )  );
    
        if (JDFR < 0)
            JDFR=0;
    }
    else
        JDFR = 0;

    pblock[local_idx] += sqrt(pow(Alfa*KDFR,2) + pow(JDFR,2));

    // wait for other work items to finish their calculations
    barrier (CLK_LOCAL_MEM_FENCE);

    //
    // output all at once
    //
    pl_out[element_idx] = pblock[local_idx];
}




/**
 * Calculates the losses introduced by the antenna, taking into account its
 * horizontal and vertical diagrams.
 * All coordinates are expressed in raster pixels.
 *
 * pixel_res ................. size of one raster pixel (in meters);
 * raster_north .............. northern limit of the raster area (in meters);
 * raster_west ............... western limit of the raster area (in meters);
 * ncols ..................... total number of columns in the area;
 * rx_height ................. height of the receiver above ground (in mts);
 * frequency ................. transmitter frequency in Mhz;
 * gain ...................... gain (in dB), relative to the isotropic antenna;
 * beam_dir .................. direction of the antenna beam;
 * mech_tilt ................. whether antenna mechanical tilt is present;
 * main_zone_horiz ........... indicates the horizontal loss that defines the
 *                             main antenna beam;
 * main_zone_vert ............ indicates the vertical loss that defines the 
 *                             main antenna beam;
 * sec_zone_horiz ............ indicates the horizontal loss that defines the
 *                             secondary antenna beam;
 * sec_zone_vert ............. indicates the vertical loss that defines the
 *                             secondary antenna beam;
 * tx_data ................... coordinates, elevation and anthenna height of 
 *                             the transmitter (coordinates expressed in meters);
 * tx_offset ................. offset of the current transmitter within the
 *                             area, including the calculation radius (in raster
 *                             pixels);
 * __gl float horiz_diag_in .. horizontal antenna diagram;
 * __gl float vert_diag_in ... vertical antenna diagram;
 * __gl float dem_in ......... digital elevation model data for the area;
 * __gl float radio_zone_out . indicates the radio zone to which each point belongs;
 * __gl float ant_loss_out ... indicates the loss, introduced by the antenna, on every point;
 * __lo real2 pblock ......... pointer to local memory, used to speed up the
 *                             whole process. Intermediate numbers, like 
 *                             distance and height, are saved here.
 */
__kernel void antenna_influence_kern (const real pixel_res,
                                      const real raster_north,
                                      const real raster_west,
                                      const int ncols, 
                                      const real rx_height,
                                      const real frequency,
                                      const real gain, 
                                      const int beam_dir,
                                      const int mech_tilt,
                                      const int main_zone_horiz,
                                      const int main_zone_vert,
                                      const int sec_zone_horiz,
                                      const int sec_zone_vert,
                                      const real4 tx_data, 
                                      const int2 tile_offset,
                                      __global real *horiz_diag_in,
                                      __global real *vert_diag_in,
                                      __global real *dem_in,
                                      __global char *radio_zone_out,
                                      __global real *ant_loss_out,
                                      __local real2 *pblock)
{
    real height_diff;
    real angle;
    real horiz_angle, horiz_loss;
    real vert_angle, vert_loss;
    real dist_Tx_Rx_in_km;
    real2 dist;
    uint element_idx = get_global_id(1) * get_global_size(0) + 
                       get_global_id(0);
 
    // receiver coordinates in pixels 
    int2 rx_coord_in_pixels  = (int2) ((int)tile_offset.x + get_global_id(0),
                                       (int)tile_offset.y + get_global_id(1));
    
    // receiver coordinates in meters
    real2 rx_coord_in_meters = (real2) ((real) rx_coord_in_pixels.x * pixel_res + raster_west,
                               (real) raster_north - rx_coord_in_pixels.y * pixel_res); 
    real2 tx_coord_in_meters = (real2) ((real) tx_data.x,
                                        (real) tx_data.y);
    
    // calculate the distance using the raster coordinates
    dist_Tx_Rx_in_km  = distance (tx_coord_in_meters, 
                                  rx_coord_in_meters);
    dist_Tx_Rx_in_km /= 1000;
    dist.x = rx_coord_in_meters.x - tx_coord_in_meters.x;
    dist.y = rx_coord_in_meters.y - tx_coord_in_meters.y;
    
    // calculate the horizontal angle between transmitter and receiver;
    // the arctan cannot be calculated if any of the involved numbers is 0
    if (dist.x == 0)
        dist.x = 0.001;
    if (dist.y == 0)
        dist.y = 0.001;
    angle = atan (dist.x / dist.y);
    if (angle < 0)
        angle = -angle;

    // take care of the orientation
    if ((dist.x >= 0) && (dist.y >= 0))
        horiz_angle = angle;
    else if ((dist.x < 0) && (dist.y >= 0))
        horiz_angle = 2*_PI_ - angle;
    else if ((dist.x < 0) && (dist.y < 0))
        horiz_angle = _PI_ + angle;
    else // (dist.y < 0 && dist.x >= 0)
        horiz_angle = _PI_ - angle;

    // convert from radians to degrees and substract the beam direction
    horiz_angle  = (horiz_angle * 180 / _PI_) - beam_dir;
    horiz_angle += (horiz_angle < 0) ? 360 : 0;

    // avoid to reading unallocated data (antenna diagram comprises 0 to 359 deg.)
    angle = ceil (horiz_angle);
    angle = (angle == 360) ? 0 : angle;
    
    // interpolation
    int index = (int) floor (horiz_angle);
    horiz_loss  = horiz_diag_in[index];
    horiz_loss += ((horiz_diag_in[(int) angle] - horiz_diag_in[index]) *
                   (horiz_angle - floor (horiz_angle)));

    //
    // continue with the vertical angle and loss calculation
    //
    // element index within the digital elevation model
    uint dem_idx = rx_coord_in_pixels.y * ncols + rx_coord_in_pixels.x;

    // effective height between Tx and Rx
    height_diff = tx_data.z - dem_in[dem_idx] - rx_height;

    vert_angle  = atan (height_diff / (dist_Tx_Rx_in_km * 1000));
    vert_angle  = vert_angle * 180 / _PI_;

    // impact of the mechanical tilt with respect to the horizontal angle
    real mech_tilt_impact;

    if (horiz_angle >= 0 && horiz_angle <= 180)
    {
        mech_tilt_impact = (real) mech_tilt * (1 - (horiz_angle / 90));
    }
    else if (horiz_angle > 180 && horiz_angle <= 360)
    {
        mech_tilt_impact = (real) mech_tilt * ((horiz_angle / 90) - 3);
    }
    // adjust the vertical angle with the mechanical tilt impact
    vert_angle = vert_angle - mech_tilt_impact;
    vert_angle += (vert_angle < 0) ? 360 : 0;

    // prevent reading unallocated data (antenna diagram comprises 0 to 359 deg.)
    angle = ceil (vert_angle);
    angle = (angle == 360) ? 0 : angle;

    // interpolation
    vert_loss = vert_diag_in[(int) floor (vert_angle)] + 
                ((vert_diag_in[(int) angle] - 
                  vert_diag_in[(int) floor (vert_angle)]) * 
                  (vert_angle - floor (vert_angle)));

    // cache the horizontal and vertical losses
    uint local_idx = get_local_id(1) * get_local_size(0) + 
                     get_local_id(0);
    pblock[local_idx] = (real2) (horiz_loss,
                                 vert_loss);

    // wait for other work items to cache their values
    barrier (CLK_LOCAL_MEM_FENCE);

    // mark the radio zone based on the distances defined by E/// propagation;
    // we should mark the pixels which are too close for calculation 
    // (less than 200 mts), and those within distances for which the 
    // propagation model is defined
    if (dist_Tx_Rx_in_km < 0.01)
        radio_zone_out[element_idx] &= _RADIO_ZONE_MODEL_DISTANCE_OFF_;
    else 
        radio_zone_out[element_idx] |= _RADIO_ZONE_MODEL_DISTANCE_ON_;

    // check if this point is within the main antenna beam, i.e. 
    // horizontally and vertically introduced loss is within threshold
    if ((pblock[local_idx].x <= (real) main_zone_horiz) && 
        (pblock[local_idx].y <= (real) main_zone_vert))
        radio_zone_out[element_idx] |= _RADIO_ZONE_MAIN_BEAM_ON_;
    // check if this point is within the secondary main antenna beam, i.e. 
    // a larger main antenna beam not including the main one (only the difference)
    else if ((pblock[local_idx].x <= (real) sec_zone_horiz) && 
             (pblock[local_idx].y <= (real) sec_zone_vert))
        radio_zone_out[element_idx] |= _RADIO_ZONE_SECONDARY_BEAM_ON_;

    // combine pathloss with determined diagram angles and antenna gain
    ant_loss_out[element_idx] = (real) (horiz_loss + vert_loss - gain);
}



/**
 * Sums the contents of the first vector to the second one.
 *
 * __gl vect_in ........ the first vector;
 * __gl vect_out ....... the second and output vector, to which the first
 *                       one will be sumed;
 * __lo pblock ......... pointer to local memory, used to speed up the
 *                       whole process. Intermediate results are saved
 *                       here.
 */
__kernel void 
vector_sum_kern (__global real *vect_in,
                 __global real *vect_out,
                 __local  real *pblock)
{
    uint element_idx = get_global_id(1) * get_global_size(0) + 
                       get_global_id(0);
    uint local_idx = get_local_id(1) * get_local_size(0) + 
                     get_local_id(0);

    // cache the sum to achieve coallesced memory access
    pblock[local_idx]  = vect_in[element_idx];
    pblock[local_idx] += vect_out[element_idx];

    // wait for other work items to cache their values
    barrier (CLK_LOCAL_MEM_FENCE);

    // output all at once
    vect_out[element_idx] = pblock[local_idx];
}



/**
 * Calculate the objective function, i.e. mean square error of the path-loss
 * prediction against field measurements per radio zone.
 *
 * tx_power ............... transmitting power at the antenna, in dBm;
 * __gl m_field_meas_in ... the field-measurements matrix;
 * __gl m_radio_zone_in ... the radio-zone matrix;
 * __gl m_loss_in ......... the path-loss matrix;
 * __gl v_partial_sum_out . the vector where the sum of the elements in each
 *                          row is saved;
 * __lo pblock ............ pointer to local memory, used to speed up the
 *                          whole process. Intermediate results are saved
 *                          here.
 */
__kernel void 
obj_func_kern (const real     tx_power,
               const char     radio_zone,
               __global real *m_field_meas_in,
               __global char *m_radio_zone_in,
               __global real *m_loss_in,
               __global real *v_partial_sum_out,
               __local  real *pblock)
{
    uint element_idx = get_global_id (0);
    uint local_idx = get_local_id (0);

    // get the data from global memory
    char rz = m_radio_zone_in[element_idx];
    pblock[local_idx]  = tx_power;
    pblock[local_idx] -= m_loss_in[element_idx];
    pblock[local_idx] -= m_field_meas_in[element_idx];

    // make sure the calculated value is valid
    if (isnan (pblock[local_idx]))
    {
        // value is not a number
        pblock[local_idx] = 0;
    }
    else
    {
        // make sure we are within the target radio zone
        if ((rz & radio_zone) > 0)
            // within the radio zone
            pblock[local_idx] *= pblock[local_idx];
        else
            // out of the radio zone
            pblock[local_idx] = 0;
    }
    // wait for other work items to calculate their values
    barrier (CLK_LOCAL_MEM_FENCE);

    // calculate the sum of this row with only one thread
    if (local_idx == 0)
    {
        uint i, local_size = get_local_size (0);
        
        // accumulate the sum of the rest of the elements
        for (i = local_idx; i < local_size; i ++)
            pblock[local_idx] += pblock[i];

        // save the calculated value
        uint out_element_idx = element_idx / local_size;
        v_partial_sum_out[out_element_idx] = pblock[local_idx];
    }
}

#endif

