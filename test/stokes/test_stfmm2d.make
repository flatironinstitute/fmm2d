PROJECT = int2-stfmm2d

HOST = gcc
HOST = gcc-openmp

arch = native
# arch = x86-64

# FC - fortran compiler
# FFLAGS - fortran compiler flags

ifeq ($(HOST),gcc)
    FC=gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=${arch} -std=legacy 
endif

ifeq ($(HOST),gcc-openmp)
    FC = gfortran 
    FFLAGS=-fPIC -O3 -funroll-loops -march=${arch} -fopenmp -std=legacy
endif

# Test objects
#

BHFMM = ../../src/biharmonic
STFMM = ../../src/stokes
COM = ../../src/common


.PHONY: all clean list


OBJECTS =  test_stfmm2d.o \
  $(COM)/prini.o \
  $(COM)/hkrand.o \
  $(COM)/dlaran.o \
  $(COM)/pts_tree2d.o \
  $(COM)/tree_routs2d.o \
  $(COM)/cumsum.o \
  $(BHFMM)/bhfmm2d.o \
  $(COM)/fmmcommon2d.o \
  $(BHFMM)/bhrouts2d.o \
  $(STFMM)/stfmm2d.o \
  $(STFMM)/stokkernels2d.o \
  $(BHFMM)/bhkernels2d.o \
  $(BHFMM)/bhndiv2d.o \
  $(BHFMM)/bh2dterms.o 


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





  
