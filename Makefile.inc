# common definitions used in all Makefiles
GRASS_TOPDIR = /usr/local/grass-6.4.2

CC          = gcc
#PRJCFLAGS   = -D_PERFORMANCE_METRICS_ -fno-strict-aliasing -mtune=native -O -Wall -pipe
PRJCFLAGS   = -fno-strict-aliasing -mtune=native -O -Wall -pipe
PRJINCS     = -I. -I${HOME}/include -I$(GRASS_TOPDIR)/include $(shell mpicc --showme:compile)
PRJLIBS     = -L${HOME}/lib -L$(GRASS_TOPDIR)/lib -lgrass_gis -lgrass_datetime -lgrass_dbmibase -lgrass_dbmiclient -lpapi -lm -Wl,-rpath=$(GRASS_TOPDIR)/lib $(shell mpicc --showme:link)
PRJSRCS     = $(wildcard *.c)
PRJOBJS     = $(PRJSRCS:.c=.o)
RM          = rm
ECHO        = echo

SHELL       = /bin/bash


.SILENT :
