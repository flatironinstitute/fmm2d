PROJECT = int2-bhrouts2d

HOST = gcc
#HOST = gcc-openmp
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

ifeq ($(HOST),intel)
    FC = ifort
    FFLAGS=-fPIC -O2 -qopenmp
endif

# Test objects
#

BHFMM = ../../src/biharmonic
COM = ../../src/common


.PHONY: all clean list


OBJECTS =  test_bhrouts2d.o \
  $(COM)/prini.o \
  $(COM)/hkrand.o \
  $(COM)/dlaran.o \
  $(COM)/fmmcommon2d.o \
  $(BHFMM)/bhrouts2d.o \
  $(BHFMM)/bhkernels2d.o \
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





  
