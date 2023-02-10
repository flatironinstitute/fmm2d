PROJECT = int2-test-mbhfmm2d

HOST = gcc
#HOST = gcc-dbg
#HOST = gcc-openmp
#HOST = intel

# FC - fortran compiler
# FFLAGS - fortran compiler flags

ifeq ($(HOST),gcc)
    FC=gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -std=legacy 
endif

ifeq ($(HOST),gcc-dbg)
    FC=gfortran 
    FFLAGS=-g -pg -std=legacy 
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy
endif

ifeq ($(HOST),intel)
    FC = ifort
    FFLAGS=-fPIC -O2 -qopenmp
endif

# Test objects
#

HFMM = ../../src/helmholtz
LFMM = ../../src/laplace
MBHFMM = ../../src/modified-biharmonic
COM = ../../src/common


.PHONY: all clean list


OBJECTS =  test_mbhfmm2d.o \
  $(COM)/prini.o \
  $(COM)/hkrand.o \
  $(COM)/dlaran.o \
  $(COM)/pts_tree2d.o \
  $(COM)/tree_routs2d.o \
  $(COM)/cumsum.o \
  $(HFMM)/hfmm2d.o \
  $(HFMM)/hfmm2dwrap.o \
  $(HFMM)/hfmm2dwrap_vec.o \
  $(HFMM)/hndiv2d.o \
  $(COM)/fmmcommon2d.o \
  $(COM)/cdjseval2d.o \
  $(COM)/dfft_threadsafe.o \
  $(COM)/next235.o \
  $(HFMM)/helmrouts2d.o \
  $(HFMM)/helmkernels2d.o \
  $(COM)/hank103.o \
  $(HFMM)/h2dcommon.o \
  $(HFMM)/wideband2d.o \
  $(HFMM)/h2dterms.o \
  $(LFMM)/l2dterms.o \
  $(MBHFMM)/mbhkernels2d.o \
  $(MBHFMM)/mbhrouts2d.o \
  $(MBHFMM)/mbhgreen2d.o \
  $(MBHFMM)/mbhfmm2d.o 


#

%.o : %.f %.h
	$(FC) $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(PROJECT)
	$(FC) $(FFLAGS) -o $(PROJECT) $(OBJECTS)
	./$(PROJECT)

clean:
	rm -f $(OBJECTS)
	rm -f $(PROJECT)

list: $(SOURCES)
	$(warning Requires:  $^)





  
