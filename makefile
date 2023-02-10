# Makefile for FMM2d
# # This is the only makefile; there are no makefiles in subdirectories.
# Users should not need to edit this makefile (doing so would make it
# hard to stay up to date with repo version). Rather in order to
# change OS/environment-specific compilers and flags, create 
# the file make.inc, which overrides the defaults below (which are 
# for ubunutu linux/gcc system). 


# compiler, and linking from C, fortran
CC=gcc
CXX=g++
FC=gfortran


# set compiler flags for c and fortran
FFLAGS= -fPIC -O3 -march=native -funroll-loops -std=legacy -w 
CFLAGS= -std=c99 
CFLAGS+= $(FFLAGS) 
CXXFLAGS= -std=c++11 -DSCTL_PROFILE=-1
CXXFLAGS+=$(FFLAGS)

# set linking libraries
CLIBS = -lgfortran -lm -ldl 
LIBS = -lm

# extra flags for multithreaded: C/Fortran, MATLAB
OMPFLAGS =-fopenmp
OMPLIBS =-lgomp 

# Python Exetucable
PYTHON=python


# flags for MATLAB MEX compilation..
MFLAGS=-compatibleArrayDims -DMWF77_UNDERSCORE1 
MWFLAGS=-c99complex 
MOMPFLAGS = -D_OPENMP

# location of MATLAB's mex compiler
MEX=mex

# For experts, location of Mwrap executable
MWRAP=../../mwrap/mwrap
MEXLIBS=-lm -lstdc++ -ldl -lgfortran

FMM_INSTALL_DIR=$(PREFIX)
ifeq ($(PREFIX),)
	FMM_INSTALL_DIR = ${HOME}/lib
endif

DYLIBS = $(LIBS)

LIBNAME=libfmm2d
DYNAMICLIB = $(LIBNAME).so
STATICLIB = $(LIBNAME).a
LIMPLIB = $(DYNAMICLIB)

LLINKLIB = -lfmm2d



# For your OS, override the above by placing make variables in make.inc
-include make.inc

# additional compile flags for FAST_KER
ifeq ($(FAST_KER),ON)
  LIBS += -lstdc++
  DYLIBS += -lstdc++
  CLIBS += -lstdc++
  FFLAGS += -lstdc++
  CFLAGS += -lstdc++
  OMP = ON
endif


# multi-threaded libs & flags needed
ifneq ($(OMP),OFF)

CFLAGS += $(OMPFLAGS)
FFLAGS += $(OMPFLAGS)
MFLAGS += $(MOMPFLAGS)
LIBS += $(OMPLIBS)
DYLIBS += $(OMPLIBS)
MEXLIBS += $(OMPLIBS)

endif


# vectorized kernel directory
SRCDIR = ./vec-kernels/src
INCDIR = ./vec-kernels/include
LIBDIR = lib-static

# objects to compile
#
# Common objects
COM = src/common
COMOBJS = $(COM)/cdjseval2d.o $(COM)/dfft_threadsafe.o \
	$(COM)/fmmcommon2d.o $(COM)/next235.o $(COM)/prini.o \
	$(COM)/tree_routs2d.o $(COM)/pts_tree2d.o \
	$(COM)/cumsum.o $(COM)/hank103.o

# Helmholtz objects
HELM = src/helmholtz
HOBJS = $(HELM)/h2dcommon.o $(HELM)/h2dterms.o \
	$(HELM)/helmrouts2d.o $(HELM)/hfmm2d.o $(HELM)/hfmm2dwrap.o \
	$(HELM)/wideband2d.o $(HELM)/hfmm2dwrap_vec.o \
	$(HELM)/hndiv2d.o $(HELM)/hfmm2d_ndiv.o $(HELM)/hfmm2d_mps.o \
	$(HELM)/hfmm2d_mps_ndiv.o

# laplace objects
LAP = src/laplace
LOBJS = $(LAP)/l2dterms.o \
	$(LAP)/laprouts2d.o $(LAP)/lfmm2d.o $(LAP)/lfmm2dwrap.o \
	$(LAP)/lfmm2dwrap_vec.o \
	$(LAP)/cfmm2d.o $(LAP)/cfmm2dwrap.o \
	$(LAP)/cfmm2dwrap_vec.o \
	$(LAP)/rfmm2d.o $(LAP)/rfmm2dwrap.o \
	$(LAP)/rfmm2dwrap_vec.o $(LAP)/lndiv2d.o \
	$(LAP)/rfmm2d_ndiv.o $(LAP)/lfmm2d_ndiv.o \
	$(LAP)/cfmm2d_ndiv.o 

BH = src/biharmonic
BHOBJS = $(BH)/bh2dterms.o \
	$(BH)/bhrouts2d.o $(BH)/bhfmm2d.o $(BH)/bhndiv2d.o \
	$(BH)/bhfmm2dwrap.o

ST = src/stokes
STOBJS = $(ST)/stfmm2d.o 


ifneq ($(FAST_KER),ON)
LOBJS += $(LAP)/lapkernels2d.o
LOBJS += $(LAP)/rlapkernels2d.o
LOBJS += $(LAP)/cauchykernels2d.o
HOBJS += $(HELM)/helmkernels2d.o
BHOBJS += $(BH)/bhkernels2d.o
STOBJS += $(ST)/stokkernels2d.o
endif

ifeq ($(FAST_KER),ON)
LOBJS += $(LAP)/rlapkernels2d.o
LOBJS += $(LAP)/cauchykernels2d.o
LOBJS += $(LAP)/lapkernels2d.o
HOBJS += $(HELM)/helmkernels2d.o
BHOBJS += $(BH)/bhkernels2d.o
STOBJS += $(ST)/stokkernels2d.o
COMOBJS+= $(SRCDIR)/libkernels.o
endif

# Test objects
TOBJS = $(COM)/hkrand.o $(COM)/dlaran.o

OBJS = $(COMOBJS) $(HOBJS) $(LOBJS) $(BHOBJS) $(STOBJS) 

.PHONY: usage install lib test all python python-dist 

default: usage

cxxkernel: $(CXXOBJ)

$(SRCDIR)/libkernels.o: $(SRCDIR)/libkernels.cpp
		$(CXX) $(CXXFLAGS) -I$(INCDIR) -c $^ -o $@

usage:
	@echo "Makefile for FMM2d. Specify what to make:"
	@echo "  make install - compile and install the main library"
	@echo "  make install PREFIX=(INSTALL_DIR) - compile and install the main library at custom location given by PREFIX"
	@echo "  make lib - compile the main library (in lib/ and lib-static/)"
	@echo "  make test - compile and run validation tests (will take a couple of mins)"
	@echo "  make matlab - compile matlab interfaces"
	@echo "  make python - compile and test python interfaces"
	@echo "  make python-dist    compile python interfaces for distribution"
	@echo "  make objclean - removal all object files, preserving lib & MEX"
	@echo "  make clean - also remove lib, MEX, py, and demo executables"
	@echo "For faster (multicore) making, append the flag -j"
	@echo "  'make [task] OMP=OFF' for single-threaded"


# implicit rules for objects (note -o ensures writes to correct dir)
%.o: %.cpp %.h
	$(CXX) -c $(CXXFLAGS) $< -o $@
%.o: %.c %.h
	$(CC) -c $(CFLAGS) $< -o $@
%.o: %.f %.h
	$(FC) -c $(FFLAGS) $< -o $@
%.o: %.f90 
	$(FC) -c $(FFLAGS) $< -o $@

# build the library...
lib: $(STATICLIB) $(DYNAMICLIB)
ifneq ($(OMP),OFF)
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, multithread versions"
else
	@echo "$(STATICLIB) and $(DYNAMICLIB) built, single-threaded versions"
endif


install: $(STATICLIB) $(DYNAMICLIB)
	echo $(FMM_INSTALL_DIR)
	mkdir -p $(FMM_INSTALL_DIR)
	cp -f lib/$(DYNAMICLIB) $(FMM_INSTALL_DIR)/
	cp -f lib-static/$(STATICLIB) $(FMM_INSTALL_DIR)/
	[ ! -f lib/$(LIMPLIB) ] || cp lib/$(LIMPLIB) $(FMM_INSTALL_DIR)/
	@echo "Make sure to include " $(FMM_INSTALL_DIR) " in the appropriate path variable"
	@echo "    LD_LIBRARY_PATH on Linux"
	@echo "    PATH on windows"
	@echo "    DYLD_LIBRARY_PATH on Mac OSX (not needed if default installation directory is used"
	@echo " "
	@echo "In order to link against the dynamic library, use -L"$(FMM_INSTALL_DIR) " -lfmm2d"


$(STATICLIB): $(OBJS) 
	ar rcs $(STATICLIB) $(OBJS)
	mv $(STATICLIB) lib-static/
$(DYNAMICLIB): $(OBJS) 
	$(FC) -shared -fPIC $(OBJS) -o $(DYNAMICLIB) $(DYLIBS) 
	mv $(DYNAMICLIB) lib/
	[ ! -f $(LIMPLIB) ] || mv $(LIMPLIB) lib/


# testing routines
#
test: $(STATICLIB) $(TOBJS) test/hfmm2d test/hfmm2d_vec test/lfmm2d test/lfmm2d_vec \
		test/cfmm2d test/cfmm2d_vec test/rfmm2d test/rfmm2d_vec test/bhfmm2d test/hfmm2d_mps \
		test/hfmm2d_hf 
	(cd test/helmholtz; ./run_helmtest.sh)
	(cd test/laplace; ./run_laptest.sh)
	(cd test/biharmonic; ./run_bhtest.sh)
	cat print_testreshelm.txt
	cat print_testreslap.txt
	cat print_testresbh.txt
	rm print_testreshelm.txt
	rm print_testreslap.txt
	rm print_testresbh.txt

test/hfmm2d:
	$(FC) $(FFLAGS) test/helmholtz/test_hfmm2d.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/helmholtz/int2-test-hfmm2d $(LIBS)

test/hfmm2d_hf:
	$(FC) $(FFLAGS) test/helmholtz/test_hfmm2d_hf.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/helmholtz/int2-test-hfmm2d-hf $(LIBS)

test/hfmm2d_vec:
	$(FC) $(FFLAGS) test/helmholtz/test_hfmm2d_vec.f $(TOBJS) $(COMOBJS) $(HOBJS) -o test/helmholtz/int2-test-hfmm2d-vec  $(LIBS)

test/hfmm2d_mps:
	$(FC) $(FFLAGS) test/helmholtz/test_hfmm2d_mps.f90 $(TOBJS) $(COMOBJS) $(HOBJS) -o test/helmholtz/int2-test-hfmm2d-mps  $(LIBS)

test/lfmm2d:
	$(FC) $(FFLAGS) test/laplace/test_lfmm2d.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/laplace/int2-test-lfmm2d $(LIBS)

test/lfmm2d_vec:
	$(FC) $(FFLAGS) test/laplace/test_lfmm2d_vec.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/laplace/int2-test-lfmm2d-vec $(LIBS) 

test/cfmm2d:
	$(FC) $(FFLAGS) test/laplace/test_cfmm2d.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/laplace/int2-test-cfmm2d $(LIBS)

test/cfmm2d_vec:
	$(FC) $(FFLAGS) test/laplace/test_cfmm2d_vec.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/laplace/int2-test-cfmm2d-vec $(LIBS) 

test/rfmm2d:
	$(FC) $(FFLAGS) test/laplace/test_rfmm2d.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/laplace/int2-test-rfmm2d $(LIBS)

test/rfmm2d_vec:
	$(FC) $(FFLAGS) test/laplace/test_rfmm2d_vec.f $(TOBJS) $(COMOBJS) $(LOBJS) -o test/laplace/int2-test-rfmm2d-vec $(LIBS) 

test/bhfmm2d:
	$(FC) $(FFLAGS) test/biharmonic/test_bhfmm2d.f $(TOBJS) $(COMOBJS) $(BHOBJS) -o test/biharmonic/int2-test-bhfmm2d $(LIBS)


#python
python: $(STATICLIB)
	cd python && \
	FMM_FLIBS='$(LIBS) $(OMPFLAGS)' FMM_FFLAGS='$(FFLAGS)' $(PYTHON) -m pip install -e .

python-dist: $(STATICLIB)
	cd python && \
	FMM_FLIBS='$(LIBS) $(OMPFLAGS)' FMM_FFLAGS='$(FFLAGS)' $(PYTHON) setup.py bdist_wheel

# matlab..
MWRAPFILE = fmm2d
GATEWAY = $(MWRAPFILE)

matlab:	$(STATICLIB) matlab/$(GATEWAY).c 
	$(MEX) matlab/$(GATEWAY).c lib-static/$(STATICLIB) $(MFLAGS) \
	-output matlab/fmm2d $(MEXLIBS) 


mex:  $(STATICLIB)
	cd matlab; $(MWRAP) $(MWFLAGS) -list -mex $(GATEWAY) -mb $(MWRAPFILE).mw;\
	$(MWRAP) $(MWFLAGS) -mex $(GATEWAY) -c $(GATEWAY).c $(MWRAPFILE).mw;\
	$(MEX) $(GATEWAY).c ../lib-static/$(STATICLIB) $(MFLAGS) -output $(MWRAPFILE) \
	$(MEXLIBS); \

clean: objclean
	rm -f lib-static/*.a lib/*.so lib/*.dll lib/*.lib
	rm -f test/laplace/int2-*
	rm -f test/helmholtz/int2-*

objclean: 
	rm -f $(OBJS) $(COBJS) $(TOBJS)
	rm -f test/laplace/*.o test/helmholtz/*.o
