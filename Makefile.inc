# common definitions used in all Makefiles
GRASS_TOPDIR = /usr/lib/grass64
CUDA_DIR = /usr/local/cuda

CC          = gcc
PRJCFLAGS   = -D_PERFORMANCE_METRICS_ -fno-strict-aliasing -mtune=native -g -O0 -Wall -Wno-unused-function -pipe
PRJINCS     = -I. -I$(CUDA_DIR)/include -I$(GRASS_TOPDIR)/include $(shell mpicc --showme:compile)
PRJLIBS     = -L$(GRASS_TOPDIR)/lib -lgrass_gis -lgrass_datetime -lgrass_dbmibase -lgrass_dbmiclient -lpapi -lm -Wl,-rpath=$(GRASS_TOPDIR)/lib $(shell mpicc --showme:link)
RM          = rm
ECHO        = echo

SHELL       = /bin/bash


.SILENT :
