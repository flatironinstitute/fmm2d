#!/bin/bash

set -e -x

# Replace native compilation flags with more generic ones.
cp make.inc.manylinux make.inc

# Clean up the build and make the library.
make clean
make lib

# Test to make sure everything is ok.
make test

# Needed for pip install to work
export FMM2D_DIR=$(pwd)
# Needed for auditwheel to find the dynamic libraries
export LD_LIBRARY_PATH=${FMM2D_DIR}/lib:${LD_LIBRARY_PATH}
export FMM_FLIBS="-lm -lgomp -fopenmp"
export FMM_FFLAGS="-fPIC -O3 -funroll-loops -march=x86-64 -mtune=generic -msse4 -fcx-limited-range -std=legacy -w -fopenmp"

# Explicitly list Python versions to build
versions=("cp36-cp36m"
          "cp37-cp37m"
          "cp38-cp38"
          "cp39-cp39"
          "cp310-cp310"
          "cp311-cp311"
          "pp37-pypy37_pp73"
          "pp38-pypy38_pp73"
          "pp39-pypy39_pp73")

pys=()
for version in "${versions[@]}"; do
    pys+=("/opt/python/${version}/bin")
done

# build wheel
for pybin in "${pys[@]}"; do
    "${pybin}/pip" install --upgrade pip
    "${pybin}/pip" install auditwheel wheel twine numpy
    "${pybin}/pip" wheel ./python -w python/wheelhouse
done

# fix wheel
for whl in python/wheelhouse/fmm2dpy-*.whl; do
    auditwheel repair "$whl" -w python/wheelhouse/
done

# test wheel
for pybin in "${pys[@]}"; do
    "${pybin}/pip" install fmm2dpy -f ./python/wheelhouse/
    "${pybin}/python" ./python/test/test_bhfmm.py
    "${pybin}/python" ./python/test/test_cfmm.py
    "${pybin}/python" ./python/test/test_hfmm.py
    "${pybin}/python" ./python/test/test_lfmm.py
    "${pybin}/python" ./python/test/test_rfmm.py
done
