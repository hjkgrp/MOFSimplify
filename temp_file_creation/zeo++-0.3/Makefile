# This Makefile compiles Zeo++ into an executable named network
# To see instructions of implementation please see README

#Use GNU C++ Compilier
CC = g++

#Flags to compile with
CFLAGS = -g #-fPIC
CFLAGS_DYLIB = -fPIC

UNAME_S := $(shell uname -s)
ifeq ($(UNAME_S),Linux)
    LDFLAGS_DYLIB = -shared
else ifeq ($(UNAME_S),Darwin)
    LDFLAGS_DYLIB = -dynamiclib
endif

#Libraries to use while compiling
LIB = -lvoro++

# The relative of path of the main library source files
VOROINCLDIR = -Ivoro++/src 
VOROLINKDIR = -Lvoro++/src

#Object files to be created for network
MAIN_OBJ = main.o
FB_OBJ = framework_builder.o
MTA_OBJ = molecule_to_abstract.o
OBJS = network.o \
  networkinfo.o \
  networkio.o \
  instructions.o \
  networkanalysis.o \
  area_and_volume.o \
  segment.o \
  string_additions.o \
  symbcalc.o \
  ray.o \
  graphstorage.o \
  networkstorage.o \
  voronoicell.o \
  channel.o \
  mindist.o \
  feature.o \
  v_network.o \
  geometry.o \
  holograms.o \
  grid.o \
  psd.o \
  sphere_approx.o \
  networkaccessibility.o \
  poreinfo.o \
  cluster.o \
  net.o \
  rmsd.o \
  symmetry.o \
  cycle.o \
  arguments.o \
  material.o \
  OMS.o
DY_OBJS=${OBJS}
STATIC_OBJS=${OBJS}

TOT_OBJS=$(MAIN_OBJ) $(FB_OBJ) $(MTA_OBJ) $(OBJS)
SRC=$(patsubst %.o,%.cc,$(TOT_OBJS))

INCS = network.h \
  heap.h \
  general.h \
  networkstorage.h \
  networkinfo.h \
  mindist.h \
  networkio.h \
  geometry.h \
  graphstorage.h \
  graphstorage.h \
  instructions.h \
  voronoicell.h \
  channel.h \
  segment.h \
  networkanalysis.h \
  area_and_volume.h \
  holograms.h \
  string_additions.h \
  symbcalc.h \
  ray.h \
  feature.h \
  v_network.h \
  grid.h \
  psd.h \
  sphere_approx.h \
  networkaccessibility.h \
  poreinfo.h \
  cluster.h \
  net.h \
  rmsd.h \
  symmetry.h \
  cycle.h \
  arguments.h \
  zeojobs.h \
  material.h \
  OMS.hh

# List of executables
EXECUTABLES = network framework_builder molecule_to_abstract

DY_LIB = libzeo++.so

STATIC_LIB = libzeo++.a

# Makefile rules
all: network framework_builder molecule_to_abstract

include Makefile.dep

depend:
	$(CC) -MG $(CFLAGS) -MM $(SRC) | sed 's/voro++\.hh //g' > Makefile.dep

network: ${OBJS} ${MAIN_OBJ}
	@echo
	@echo Linking $@
	$(CC) $(VOROINCLDIR) $(VOROLINKDIR) -o $@ $(CFLAGS) ${OBJS} ${MAIN_OBJ} $(LIB)

framework_builder: ${OBJS} ${FB_OBJ}
	@echo
	@echo Linking $@
	$(CC) $(VOROINCLDIR) $(VOROLINKDIR) -o $@ $(CFLAGS) ${OBJS} ${FB_OBJ} $(LIB)

molecule_to_abstract: ${OBJS} ${MTA_OBJ}
	@echo
	@echo Linking $@
	$(CC) $(VOROINCLDIR) $(VOROLINKDIR) -o $@ $(CFLAGS) ${OBJS} ${MTA_OBJ} $(LIB)

statlib: ${STATIC_LIB}
${STATIC_LIB}: ${STATIC_OBJS}
	@echo
	@echo Linking $@
	ar rs $@ $^ #$(VOROINCLDIR) $(VOROLINKDIR) $(LIB)

dylib: CFLAGS += -fPIC 
dylib: ${DY_LIB}

${DY_LIB}: ${DY_OBJS}
	@echo
	@echo Linking $@
	$(CC) $(VOROINCLDIR) $(VOROLINKDIR) -o $@ $(LDFLAGS_DYLIB) ${OBJS} $(LIB)

%.o: %.cc
	$(CC) $(VOROINCLDIR) $(CFLAGS) -c $<

clean: 
	rm -f $(TOT_OBJS) *.err $(EXECUTABLES)

.PHONY: all clean depend

