include ../Makefile.inc

CFLAGS = $(PRJCFLAGS)
INCS   = $(PRJINCS) -I..
OBJS   = $(PRJOBJS)

all: $(OBJS)

.c.o:
	$(ECHO) $(CC) $(CFLAGS) $(INCS) -c $< -o $@
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

clean:
	$(ECHO) Cleaning up in performance ...
	$(ECHO) $(RM) -f $(OBJS)
	$(RM) -f $(OBJS)
