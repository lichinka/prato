#ifndef _COVERAGE_CONSTANTS_FOR_GRASS_H_
#define _COVERAGE_CONSTANTS_FOR_GRASS_H_

//
// value for Pi
//
#define _PI_    3.14159265358979

//
// character buffers are 1k by default
//
#define _CHAR_BUFFER_SIZE_          1024

//
// MPI constants for master-worker communication
//
#define _COVERAGE_MASTER_RANK_      0
#define _WORKER_IS_IDLE_TAG_        100
#define _WORKER_KEEP_WORKING_TAG_   105
#define _WORKER_SHUTDOWN_TAG_       110
#define _WORKER_OPTIMIZE_TAG_       115

//
// whether the target point is whithin the main antenna beam or not
//
#define _RADIO_ZONE_MAIN_BEAM_ON_       0x01
#define _RADIO_ZONE_MAIN_BEAM_OFF_      0xfe
//
// whether the target point is within the distance limit of the prediction model
//
#define _RADIO_ZONE_SECONDARY_BEAM_ON_  0x02
#define _RADIO_ZONE_SECONDARY_BEAM_OFF_ 0xfd
//
// whether the target point is within the distance limit of the prediction model
//
#define _RADIO_ZONE_MODEL_DISTANCE_ON_  0x04
#define _RADIO_ZONE_MODEL_DISTANCE_OFF_ 0xfb

//
// dimensions of the search vector - only used for optimization
//
#define _SEARCH_VECTOR_DIMENSIONS_  4

#endif
