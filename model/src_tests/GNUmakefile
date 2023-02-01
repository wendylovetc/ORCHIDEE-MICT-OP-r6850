#
# Make sure to compile Orchidee in the first place
# 
# Notes:
#    pFunit is compiled in serial. Tests are not parallel.

BASE_FFLAGS=$(shell grep BASE_FFLAGS ../arch.fcm | cut -c 20-)
DEB_FFLAGS=$(shell grep DEBUG_FFLAGS ../arch.fcm | cut -c 20-)
ALL_FFLAGS=$(BASE_FFLAGS) $(DEB_FFLAGS)

FCOMP=$(shell grep COMPILER ../arch.fcm | cut -c 10-)
FLINK=$(shell grep LINK ../arch.fcm | cut -c 10-)
ORCHIDEE=$(shell pwd)/..
ORCH_INC=$(ORCHIDEE)/lib
ORCH_LIB=$(ORCHIDEE)/lib 
MODIPSL=$(ORCHIDEE)../.. 
# Get from config.fcm (arch dependent)
ALL_INCDIR=$(shell grep INCDIR ../config.fcm | cut -c 10-)
ALL_LIBDIR=$(shell grep LIBDIR ../config.fcm | cut -c 10-)

# pFunit makefiles
include $(PFUNIT)/include/base.mk

# Get source files
SRCS = $(wildcard *.pf)
OBJS = $(SRCS:.pf=.o)


#tests: $(OBJS)

#testSuites.inc: $(SRCS) 

# Parse files
%.F90: %.pf
	$(PFUNIT)/bin/pFUnitParser.py $< $@

# Compile tests
%.o: %.F90
	$(FCOMP) $(ALL_FFLAGS) -c $< -I$(PFUNIT)/mod -I$(PFUNIT)/include -I$(ORCH_INC) $(ALL_INCDIR) 

	
tests.x: $(OBJS) 
	$(FLINK) -o $@ -I$(PFUNIT)/mod -I$(PFUNIT)/mod -I$(PFUNIT)/include -I$(ORCH_INC) $(ALL_INCDIR)  \
	$(PFUNIT)/include/driver.F90 \
	./$^ $(LIBS) $(ALL_FFLAGS) -L$(PFUNIT)/lib -lpfunit -L$(ORCH_LIB) -lorchidee $(ALL_LIBDIR) -lioipsl


clean:
	rm -f *.F90 tests.x *.o *.mod
