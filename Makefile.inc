#
# Common definitions used in all Makefiles
#
GRASS_TOPDIR        = /usr/local/grass-6.4.2
CUDA_DIR            = /usr/local/cuda
OCL_COMMON_LIB_DIR  = ../ocl_common
PERFORMANCE_LIB_DIR = ../performance

CC          = gcc
PRJCFLAGS   = -D_PERFORMANCE_METRICS_ -fno-strict-aliasing -mtune=native -g -O0 -Wall -Wno-unused-function -pipe
PRJINCS     = -I$(CUDA_DIR)/include -I$(GRASS_TOPDIR)/include $(shell mpicc --showme:compile) -I.
PRJLIBS     = -L$(GRASS_TOPDIR)/lib -lgrass_gis -lgrass_datetime -lgrass_dbmibase -lgrass_dbmiclient -lpapi -lm -Wl,-rpath=$(GRASS_TOPDIR)/lib $(shell mpicc --showme:link)
RM          = rm
ECHO        = echo

SHELL       = /bin/bash


.SILENT :
