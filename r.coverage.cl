#include "prot.h"

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
 * Double precision definition:
 */
#ifdef cl_khr_fp64
	#pragma OPENCL EXTENSION cl_khr_fp64 : enable
#endif
#ifdef cl_amd_fp64
	#pragma OPENCL EXTENSION cl_amd_fp64 : enable
#endif
typedef double	real;
typedef double2	real2;
typedef double4 real4;

#pragma OPENCL EXTENSION cl_amd_printf : enable

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
                            __global float *dem_in,
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
 * Implements an autonomous agent optimizing the service area.
 * Agent behavior is explained in the README file.
 *
 * int ntx .................. number of transmitters being processed. This
 *                            is the length of the 'tx_pwr' vector and the
 *                            depth of the 'pl_in' cube.
 * int ncols ................ total number of columns in the area.
 * int4 mdim ................ area widthin which the agents may move.
 * int uncovered count ...... number of uncovered raster cells.
 * int uncov_coord_length ... number of elements in the vector 'uncov_coord_in'.
 * long random_seed ......... used to initialize the random number generator.
 * __gl u char pl_in ........ 3D matrix containing path-loss values in dB
 *                            (i.e. between 0 and 255) for every coordinate
 *                            (2D) and each of the 'ntx' transmitters (3rd
 *                            dimension) in the area.
 * __gl int4 offsets_in ..... vector of tuples containing the offsets for
 *                            each of the path-loss matrices (col-min, 
 *                            row-min, width, height).
 * __gl real qrm_in ........ 2D matrix containing interference values for
 *                            every coordinate in the area.
 * __gl real cov_in ........ 2D matrix containing the service area coverage.
 * __gl short uncov_coord_in. a vector containing coordinates of uncovered 
 *                            raster cells.
 * __gl int tx_pwr .......... pilot powers for 'ntx' transmitters (in 
 *                            milliwatts).
 * __lo real pblock ........ local memory used to speed up the whole 
 *                            process. The new pilot powers with the 
 *                            transmitter IDs are saved here.
 */
__kernel void agent_kern (const int ntx,
                          const int ncols,
                          const int4 mdim,
                          const int uncovered_count,
                          const int uncov_coord_length,
                          const unsigned long random_seed,
                          __global unsigned char *pl_in,
                          __global int4 *offsets_in,
                          __global real *qrm_in,
                          __global real *cov_in,
                          __global unsigned short *uncov_coord_in,
                          __global int *tx_pwr,
                          __local unsigned short *pblock)
{
    int tx;
    int idx_thread = get_local_id(0);
    int2 ccoord = (int2) (-1, -1);

    // initialize random number generator
    random_state randstate;
    seed_random (&randstate, random_seed);

    // is this a special agent?
    if ((idx_thread % 2) == 0)
    {
        // randomly select a non-covered pixel
        // from the uncovered coordinates vector
        if (uncovered_count > 0)
        {
            int idx_coord;
            do
            {
                idx_coord = (int) (random(&randstate) % uncov_coord_length);
                idx_coord --;
                ccoord.x = (int) uncov_coord_in[2*idx_coord];
                ccoord.y = (int) uncov_coord_in[2*idx_coord+1];
            } while (ccoord.x == 0);
        }
    }
    // already have a coordinate?
    if (ccoord.x < 0)
    {
        // select a random coordinate from the whole area
        ccoord.x = (int) (random(&randstate) % (mdim.z - mdim.x));
        ccoord.x += mdim.x;
        ccoord.y = (int) (random(&randstate) % (mdim.w - mdim.y));
        ccoord.y += mdim.y;
    }
    int tx_id;
    real tx_new_pwr;
    int idx_2d = ccoord.y * ncols + ccoord.x;

    // is this coordinate under network coverage?
    if (cov_in[idx_2d] < _AG_MINIMUM_GAMMA_COVERAGE_)
    {
        // coordinate NOT covered, analyze signals
        real2 max_signal = (real2) (-1.0f, _AG_EMPTY_VALUE_);

        max_signal.x = (real) (random(&randstate) % ntx);

        /*
        for (tx = 0; tx < ntx; tx++)
        {
            // is this transmitter still adjustable?
            if (tx_pwr[tx] < _AG_MAX_PILOT_POWER_MW_)
            {
                int idx_pl = (ccoord.y - offsets_in[tx].y) * offsets_in[tx].z;
                idx_pl += ccoord.x - offsets_in[tx].x;
                
                real pl_value = (real) pl_in[idx_pl];

                // transform the path loss to a linear scale
                if (pl_value > _HATA_MAX_PATH_LOSS_)
                    pl_value = _HATA_MAX_PATH_LOSS_;
            
                pl_value = (real) 1.0f - (pl_value /
                           (_HATA_MAX_PATH_LOSS_ - _HATA_MIN_PATH_LOSS_));
              
                // calculate gamma at this coordinate for this transmitter
                pl_value = ((tx_pwr[tx]/1000.0f) * pl_value) / qrm_in[idx_2d];

                // keep the maximum signal
                if (max_signal.y < pl_value)
                {
                    max_signal.x = (real) tx;
                    max_signal.y = pl_value;
                }
            }
        }*/
        // raise pilot power of the transmitter with max rcv signal
        tx_id = (int) max_signal.x;
        tx_new_pwr = _AG_INCREASE_PILOT_DB_ / 10.0f;
        tx_new_pwr = pow(10, tx_new_pwr);
    }
    else
    {
        // coordinate covered, analyze signals
        real2 min_signal = (real2) (-1.0f, _AG_EMPTY_VALUE_ * -1);
      
        min_signal.x = (real) (random(&randstate) % ntx);

        /*for (tx = 0; tx < ntx; tx++)
        {
            // is this transmitter still adjustable?
            if (tx_pwr[tx] > 0)
            {
                int idx_pl = (ccoord.y - offsets_in[tx].y) * offsets_in[tx].z;
                idx_pl += ccoord.x - offsets_in[tx].x;
                
                real pl_value = (real) pl_in[idx_pl];

                // transform the path loss to a linear scale
                if (pl_value > _HATA_MAX_PATH_LOSS_)
                    pl_value = _HATA_MAX_PATH_LOSS_;
            
                pl_value = (real) 1.0f - (pl_value /
                           (_HATA_MAX_PATH_LOSS_ - _HATA_MIN_PATH_LOSS_));
         
                // calculate gamma at this coordinate for this transmitter
                pl_value = ((tx_pwr[tx]/1000.0f) * pl_value) / qrm_in[idx_2d];

                // keep the minimum signal, if it is still adjustable
                if (min_signal.y > pl_value)
                {
                    min_signal.x = (real) tx;
                    min_signal.y = pl_value;
                }
            }
        }*/
        // lower pilot power of the transmitter with min rcv signal
        tx_id = (int) min_signal.x;
        tx_new_pwr = _AG_DECREASE_PILOT_DB_ / 10.0f;
        tx_new_pwr = pow(10, tx_new_pwr);
    }
    // save it to local memory
    pblock[2*idx_thread] = (unsigned short) tx_id;
    pblock[2*idx_thread+1] = (unsigned short) (tx_new_pwr*tx_pwr[tx_id]);

    // wait for all threads to finish
    barrier (CLK_LOCAL_MEM_FENCE);

    // save the results
    if (idx_thread == 0)
    {
        // for each agent ...
        for (tx = 0; tx < get_local_size(0); tx++)
        {
            tx_id = pblock[2*tx];
            int ty, new_pwr = pblock[2*tx+1];
            // for each other agent ...
            for (ty = 0; ty < get_local_size(0); ty++)
            {
                if ((tx != ty) && (pblock[2*ty] == tx_id))
                {
                    if (pblock[2*ty+1] > new_pwr)
                    {
                        new_pwr = -1;
                        break;
                    }
                }
            }
            // apply this change?
            if (new_pwr != -1)
            {
                // FIXME multiply this value by -1 and correct the coverage kernel!
                tx_pwr[tx_id] = new_pwr; // * -1;
            }
        }
    }
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
 */
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
 */
__kernel void hata_urban_interference (const int pixel_res,
									   const real raster_north,
									   const real raster_west,
                                       const int ncols, 
                                       const real4 tx_data, const int2 tx_offset,
                                       const real rx_height, const real frequency, 
                                       const real ahr,
                                       __global float *dem_in,
                                       __global float *qrm_out,
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



/**
 * Calculates the path-loss matrix for one transmitter using Hata's model
 * for urban areas. All coordinates are expressed in raster cells.
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
__kernel void hata_urban_pl_per_tx (const int pixel_res,
									const real raster_north,
									const real raster_west,
                                    const int ncols, 
                                    const real4 tx_data, const int2 tx_offset,
                                    const real rx_height, const real frequency, 
                                    const real ahr,
                                    __global float *dem_in,
                                    __global float *pl_out,
                                    __local real2 *pblock)
{
    uint element_idx = get_global_id(1) * get_global_size(0) + 
                       get_global_id(0);
    pl_out[element_idx] = (float) path_loss_hata_urban (pixel_res,
														raster_north,
														raster_west,
														ncols,
														tx_data, tx_offset,
														rx_height, frequency, 
														ahr,
														dem_in,
														pblock);
}



/**
 * Calculates the sectorization path-loss, as per antenna diagram.
 * All coordinates are expressed in raster pixels.
 *
 * int pixel_res ............. size of one raster pixel (in meters);
 * real raster_north ......... northern limit of the raster area (in meters);
 * real raster_west .......... western limit of the raster area (in meters);
 * int ncols ................. total number of columns in the area;
 * real4 tx_data ............. coordinates, elevation and anthenna height of 
 *                             the transmitter (coordinates expressed in meters);
 * int2 tx_offset ............ offset of the current transmitter within the
 *                             area, including the calculation radius (in raster
 *                             pixels);
 * real rx_height ............ height of the receiver above ground (in mts);
 * real frequency ............ transmitter frequency in Mhz;
 * real gain ................. gain (in dB), relative to the isotropic antenna;
 * int  beam_dir ............. direction of the antenna beam;
 * int  mech_tilt ............ whether antenna mechanical tilt is present;
 * __gl float horiz_diag_in .. horizontal antenna diagram;
 * __gl float vert_diag_in ... vertical antenna diagram;
 * __gl float dem_in ......... digital elevation model data for the area;
 * __gl float sect_out ....... sectorized data where the results are saved;
 * __lo real2 pblock ......... pointer to local memory, used to speed up the
 *                             whole process. Intermediate numbers, like 
 *                             distance and height, are saved here.
 */
__kernel void sector_kern (const int pixel_res,
                           const real raster_north,
                           const real raster_west,
                           const int ncols, 
                           const real4 tx_data, const int2 tx_offset,
                           const real rx_height,
                           const real frequency,
                           const real gain, 
                           const int beam_dir,
                           const int mech_tilt,
                           __global float *horiz_diag_in,
                           __global float *vert_diag_in,
                           __global float *dem_in,
                           __global float *sect_out,
                           __local real2 *pblock)
{
    real height_diff, angle;
    real horiz_angle, horiz_loss;
    real vert_angle, vert_loss;
    real linear_dist;
    real2 dist;
    uint element_idx = get_global_id(1) * get_global_size(0) + 
                       get_global_id(0);
 
    // receiver coordinates in pixels 
    int2 px_coord = (int2) ((int)tx_offset.x + get_global_id(0),
                            (int)tx_offset.y + get_global_id(1));
    
    // receiver coordinates in meters
    real2 rec_coord = (real2) ((real) px_coord.x * pixel_res + raster_west,
    						   (real) raster_north - px_coord.y * pixel_res); 
    real2 tx_coord = (real2) ((real) tx_data.x,
    						  (real) tx_data.y);
    
    // calculate the distance using the raster coordinates
    linear_dist = distance (tx_coord, rec_coord);
    dist.x = rec_coord.x - tx_data.x;
    dist.y = rec_coord.y - tx_data.y;
    
    // calculate the horizontal angle between transmitter and receiver
    // angle = fabs (atan2 (dist.x, dist.y));
    angle = fabs (atan (dist.x / dist.y));

	// take care of the orientation
	if ((dist.x >= 0) && (dist.y >= 0))
	{
		horiz_angle = angle;
	}
	else if ((dist.x < 0) && (dist.y >= 0))
	{
		horiz_angle = 2*PI - angle;
	}
	else if ((dist.x < 0) && (dist.y < 0))
	{
		horiz_angle = PI + angle;
	}
	else
	{
		horiz_angle = PI - angle;
	}
    // convert from radians to degrees and substract the beam direction
    horiz_angle  = (horiz_angle * 180 / PI) - beam_dir;
    horiz_angle += (horiz_angle < 0) ? 360 : 0;

    // to prevent reading unallocated data (antenna diagram comprises 0 to 359 deg.)
    angle = ceil (horiz_angle);
    angle = (angle == 360) ? 0 : angle;
    
    // interpolation
    horiz_loss = horiz_diag_in[(int) floor (horiz_angle)] +
    		     ((horiz_diag_in[(int) angle] - 
    		       horiz_diag_in[(int) floor (horiz_angle)]) *
    		       (horiz_angle - floor (horiz_angle)));

    //
    // continue with the vertical angle and loss calculation
    //
    // element index within the digital elevation model
    uint dem_idx = px_coord.y * ncols + px_coord.x;

    // effective height between Tx and Rx
    height_diff = (tx_data.z + tx_data.w) - dem_in[dem_idx] - rx_height;

    vert_angle  = atan (height_diff / linear_dist);
    vert_angle  = vert_angle * 180 / PI;

    /* calculate the impact of mechanical tilt with respect
    // to the horizontal angle
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
    */
    vert_angle += (vert_angle < 0) ? 360 : 0;
    // to prevent reading unallocated data (antenna diagram comprises 0 to 359 deg.)
    angle = ceil (vert_angle);
    angle = (angle == 360) ? 0 : angle;

    // interpolation
    vert_loss = vert_diag_in[(int) floor (vert_angle)] + 
                ((vert_diag_in[(int) angle] - 
                  vert_diag_in[(int) floor (vert_angle)]) * 
                  (vert_angle - floor (vert_angle)));

    /* cache the horizontal and vertical losses
    uint local_idx = get_local_id(1) * get_local_size(0) + 
                     get_local_id(0);
    pblock[local_idx] = (real2) (log10 (horiz_loss),
                                  log10 (vert_loss));

    // wait for other work items to cache their values
    barrier (CLK_LOCAL_MEM_FENCE);

    // resulting sectorization
    sect_out[element_idx] = (float)(pl_dbm + pblock.x + pblock.y - gain);
    */

    real ahr = (1.1 * log10 (frequency) - 0.7) * rx_height -
                (1.56 * log10 (frequency) - 0.8);

    real pl_dbm = (real) path_loss_hata_urban (pixel_res,
                                              raster_north,
                                              raster_west,
                                              ncols,
                                              tx_data, tx_offset,
                                              rx_height, frequency, 
                                              ahr,
                                              dem_in,
                                              pblock);
    // combine pathloss with determined diagram angles and antenna gain
    sect_out[element_idx] = (float)(pl_dbm + horiz_loss + vert_loss - gain);
}

