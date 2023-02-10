PROJECT = int2-hfmm2d

HOST = gcc
HOST = gcc-openmp
HOST = gcc-openmp-x86
#HOST = intel

# FC - fortran compiler
# FFLAGS - fortran compiler flags

ifeq ($(HOST),gcc)
    FC=gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -std=legacy 
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=native -fopenmp -std=legacy
endif

ifeq ($(HOST),gcc-openmp-x86)
    FC = gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=x86-64 -fopenmp -std=legacy
endif

ifeq ($(HOST),intel)
    FC = ifort
    FFLAGS=-fPIC -O2 -qopenmp
endif

# Test objects
#

HFMM = ../../src/helmholtz
COM = ../../src/common


.PHONY: all clean list


OBJECTS =  test_hfmm2d_mps.o \
  $(COM)/prini.o \
  $(COM)/hkrand.o \
  $(COM)/dlaran.o \
  $(COM)/pts_tree2d.o \
  $(COM)/tree_routs2d.o \
  $(COM)/cumsum.o \
  $(HFMM)/hfmm2d.o \
  $(HFMM)/hfmm2d_mps.o \
  $(HFMM)/hfmm2dwrap.o \
  $(HFMM)/hfmm2dwrap_vec.o \
  $(COM)/fmmcommon2d.o \
  $(COM)/cdjseval2d.o \
  $(COM)/dfft_threadsafe.o \
  $(COM)/next235.o \
  $(HFMM)/helmrouts2d.o \
  $(HFMM)/helmkernels2d.o \
  $(COM)/hank103.o \
  $(HFMM)/h2dcommon.o \
  $(HFMM)/hndiv2d.o \
  $(HFMM)/wideband2d.o \
  $(HFMM)/h2dterms.o 


#

%.o : %.f %.h
	$(FC) $(FFLAGS) $< -o $@
%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@

all: $(OBJECTS)
	rm -f $(PROJECT)
	$(FC) $(FFLAGS) -o $(PROJECT) $(OBJECTS)
	./$(PROJECT)

clean:
	rm -f $(OBJECTS)
	rm -f $(PROJECT)

list: $(SOURCES)
	$(warning Requires:  $^)





  
