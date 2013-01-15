include Makefile.inc

DIRS 	= master worker
BIN		= r.coverage
CFLAGS 	= $(PRJCFLAGS) 
INCS 	= $(PRJINCS) -I$(OCL_COMMON_LIB_DIR) -I$(PERFORMANCE_LIB_DIR)
LIBS 	= $(PRJLIBS) -L. -loclcommon -lperformance -lworker -lOpenCL
SRCS    = $(wildcard *.c) $(wildcard master/*.c)
OBJS    = $(SRCS:.c=.o)
OBJLIBS = libworker.so

all: $(BIN)

$(BIN): $(OBJS) $(OBJLIBS)
	$(ECHO) Building $(BIN) ...
	$(ECHO) $(CC) $(CFLAGS) $(INCS) -o $@ $(OBJS) $(LIBS)
	$(CC) $(CFLAGS) $(INCS) -o $@ $(OBJS) $(LIBS)

.c.o:
	$(ECHO) $(CC) $(CFLAGS) $(INCS) -c $< -o $@
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

libworker.so: force_look
	$(ECHO) Looking into worker : $(MAKE) $(MFLAGS)
	cd worker; $(MAKE) $(MFLAGS)

clean:
	$(ECHO) Cleaning up ...
	$(ECHO) $(RM) -f $(BIN) $(OBJS) $(OBJLIBS)
	-$(RM) -f $(BIN) $(OBJS) $(OBJLIBS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look: 
	true
