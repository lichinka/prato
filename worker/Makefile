include ../Makefile.inc

CFLAGS  = $(PRJCFLAGS) -fPIC
INCS    = $(PRJINCS) -I../$(OCL_COMMON_LIB_DIR) -I../$(PERFORMANCE_LIB_DIR) -I..
LIBS    = $(PRJLIBS) -L.. -loclcommon -lperformance -lOpenCL
SRCS    = $(wildcard *.c)
OBJS    = $(SRCS:.c=.o)
BIN		= worker
DIR		= worker
OBJLIBS = lib$(DIR).so

all: $(OBJLIBS) $(BIN)

lib$(DIR).so : $(OBJS)
	$(ECHO) $(CC) $(CFLAGS) $(INCS) -shared -Wl,-soname,$@ -o ../$@ $(OBJS) $(LIBS)
	$(CC) $(CFLAGS) $(INCS) -shared -Wl,-soname,$@ -o ../$@ $(OBJS) $(LIBS)

$(BIN): main.o
	$(ECHO) Building $(BIN) ...
	$(ECHO) $(CC) $(CFLAGS) $(INCS) -o $@ $< $(LIBS) -L.. -lworker
	$(CC) $(CFLAGS) $(INCS) -o $@ $< $(LIBS) -L.. -lworker

.c.o:
	$(ECHO) $(CC) $(CFLAGS) $(INCS) -c $< -o $@
	$(CC) $(CFLAGS) $(INCS) -c $< -o $@

clean:
	$(ECHO) Cleaning up in $(DIR) ...
	$(ECHO) $(RM) -f $(OBJS) main.o $(BIN)
	$(RM) -f $(OBJS) main.o $(BIN)

