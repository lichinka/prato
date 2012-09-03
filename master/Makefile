include ../Makefile.inc

CFLAGS  = $(PRJCFLAGS)
INCS    = $(PRJINCS) -I..
LIBS    = $(PRJLIBS)
SRCS    = $(wildcard *.c)
OBJS    = $(SRCS:.c=.o)

all: $(OBJS)

.c.o:
	$(ECHO) $(CC) $(CFLAGS) $(INCS) -c $< -o $@
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

clean:
	$(ECHO) Cleaning up in master ...
	$(ECHO) $(RM) -f $(OBJS)
	$(RM) -f $(OBJS)

