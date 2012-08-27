include Makefile.inc

DIRS 	= performance
BIN		= r.coverage
CFLAGS 	= $(PRJCFLAGS) 
INCS 	= $(PRJINCS)
LIBS 	= $(PRJLIBS)
SRCS 	= $(PRJSRCS)
OBJS 	= $(PRJOBJS)

all: $(BIN)

$(BIN): $(DIRS) $(OBJS) 
	$(ECHO) Building $(BIN) ...
	$(ECHO) $(OBJS)
	$(ECHO) $(CC) $(CFLAGS) $(INCS) -o $@ $(OBJS) $(DIRS:=/*.o) $(LIBS)
	$(CC) $(CFLAGS) $(INCS) -o $@ $(OBJS) $(DIRS:=/*.o) $(LIBS)

.c.o:
	$(ECHO) $(CC) $(CFLAGS) $(INCS) -c $< -o $@
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

$(DIRS): force_look
	$(ECHO) Looking into $@ : $(MAKE) $(MFLAGS)
	cd $@; $(MAKE) $(MFLAGS)

clean:
	$(ECHO) Cleaning up ...
	-$(RM) -f $(BIN) $(OBJS)
	-for d in $(DIRS); do (cd $$d; $(MAKE) clean ); done

force_look: 
	true
