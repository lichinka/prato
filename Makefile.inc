#
# Common definitions used in all Makefiles
#
#GRASS_TOPDIR        = /usr/local/src/grass-6.4.2/dist.x86_64-unknown-linux-gnu
GRASS_TOPDIR        = /home/luka/grass-6.4.2
#CUDA_DIR            = $(HOME)/AMD_SDK
CUDA_DIR            = $(ATISTREAMSDKROOT)
PAPI_DIR            = $(HOME)
OCL_COMMON_LIB_DIR  = ../ocl_common
PERFORMANCE_LIB_DIR = ../performance

CC          = /home/luka/bin/gcc
#PRJCFLAGS   = -D_DEBUG_INFO_ -D_PERFORMANCE_METRICS_ -fno-strict-aliasing -mtune=native -g -O0 -Wall -Wno-unused-function -pipe
#PRJCFLAGS   = -D_PERFORMANCE_METRICS_ -fno-strict-aliasing -mtune=native -g -O0 -Wall -Wno-unused-function -pipe
PRJCFLAGS   = -fno-strict-aliasing -mtune=native -g -O0 -Wall -Wno-unused-function -pipe
PRJINCS     = -I$(PAPI_DIR)/include -I$(CUDA_DIR)/include -I$(GRASS_TOPDIR)/include $(shell mpicc --showme:compile) -I.
PRJLIBS     = -L$(PAPI_DIR)/lib -lpapi -L$(CUDA_DIR)/lib/x86_64 $(shell mpicc --showme:link)
RM          = rm
ECHO        = echo

SHELL       = /bin/bash


.SILENT :
