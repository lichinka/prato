#
# Common definitions used in all Makefiles
#
#GRASS_TOPDIR        = /usr/local/grass-6.4.2
GRASS_TOPDIR        = /usr/local/src/grass-6.4.2/dist.x86_64-unknown-linux-gnu
#CUDA_DIR            = /usr/local/cuda
CUDA_DIR             = $(HOME)/AMD_SDK
OCL_COMMON_LIB_DIR  = ../ocl_common
PERFORMANCE_LIB_DIR = ../performance

CC          = gcc
PRJCFLAGS   = -D_PERFORMANCE_METRICS_ -fno-strict-aliasing -mtune=native -g -O0 -Wall -Wno-unused-function -pipe
#PRJCFLAGS   = -fno-strict-aliasing -mtune=native -g -O0 -Wall -Wno-unused-function -pipe
PRJINCS     = -I$(HOME)/papi/include -I$(CUDA_DIR)/include -I$(GRASS_TOPDIR)/include $(shell mpicc --showme:compile) -I.
PRJLIBS     = $(shell mpicc --showme:link)
RM          = rm
ECHO        = echo

SHELL       = /bin/bash


.SILENT :
