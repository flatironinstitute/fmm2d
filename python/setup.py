import setuptools
import string
import os
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from sys import platform

pkg_name = "fmm2dpy"

list_files=[]
list_files.append('../src/helmholtz/hfmm2dwrap.f')

FLIBS = os.getenv('FMM_FLIBS')
FLIBS = FLIBS.rstrip().split(' ')
FLIBS = list(filter(None,FLIBS))
FLIBS.append('../lib-static/libfmm2d.a')

helm = []
helm.append('hfmm2d_t_c_p')


ext_helm = Extension(
    name='fmm2d',
    sources=list_files,
    f2py_options=['only:']+helm+[':'],
#    extra_f77_compile_args=FFLAGS,
    extra_f90_compile_args=["-std=legacy"],
    extra_f77_compile_args=["-std=legacy"],
    extra_link_args=FLIBS
)

## TODO: fill in the info below
setup(
    name=pkg_name,
    version="0.1.0",
    author="Manas Rachh",
    author_email="mrachh@flatironinstitute.org",
    description="This pacakge contains basic routines for fmm in 2d",
    url="",
    packages=setuptools.find_packages(),
    install_requires=[
        "numpy",
        "pytest"
    ],
    ext_modules=[ext_helm],
    classifiers=(
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    )    
)
