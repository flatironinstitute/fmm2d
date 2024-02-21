import setuptools
import string
import os
from numpy.distutils.core import setup
from numpy.distutils.core import Extension
from sys import platform

pkg_name = "fmm2dpy"

## TODO: this should be automatically populated using "read directory, or whatever"
## TODO: fix problem with relative location for executable

list_helm=['hfmm2dwrap.f','hfmm2dwrap_vec.f','helmkernels2d.f']
list_lap=['rfmm2dwrap.f','rfmm2dwrap_vec.f','rlapkernels2d.f','lfmm2dwrap.f','lfmm2dwrap_vec.f','lapkernels2d.f','cfmm2dwrap.f','cfmm2dwrap_vec.f','cauchykernels2d.f']
list_bh=['bhfmm2dwrap.f','bhkernels2d.f']
list_stk=['stfmm2d.f','stokkernels2d.f']
list_common=[]

FLIBS = os.getenv('FMM_FLIBS')
FLIBS = FLIBS.rstrip().split(' ')
FLIBS = list(filter(None, FLIBS))
FLIBS.append('../lib-static/libfmm2d.a')
FFLAGS = os.getenv('FMM_FFLAGS')
FFLAGS = FFLAGS.rstrip().split(' ')
FFLAGS = list(filter(None, FFLAGS))



c_opts = ['_c','_d','_cd']
c_opts2 = ['c','d','cd']
st_opts = ['_s','_t','_st']
p_optsh = ['_p','_g']
p_optsh2 = ['p','g']

p_optsl = ['_p','_g','_h']
p_optsl2 = ['p','g','h']

list_int_helm = []
list_int_helm_vec = []
list_int_helm_dir = []

list_int_lap = []
list_int_lap_vec = []
list_int_lap_dir = []

for st in st_opts:
    for cd in c_opts:
        for pg in p_optsh:
            list_int_helm.append('hfmm2d'+st+cd+pg)
            list_int_helm_vec.append('hfmm2d'+st+cd+pg+'_vec')
        for pg in p_optsl:
            list_int_lap.append('rfmm2d'+st+cd+pg)
            list_int_lap.append('lfmm2d'+st+cd+pg)
            list_int_lap.append('cfmm2d'+st+cd+pg)
            list_int_lap_vec.append('rfmm2d'+st+cd+pg+'_vec')
            list_int_lap_vec.append('lfmm2d'+st+cd+pg+'_vec')
            list_int_lap_vec.append('cfmm2d'+st+cd+pg+'_vec')

list_int_bh = []
list_int_bh_dir = []
for cd in c_opts2:
    for pg in p_optsh2:
        list_int_helm_dir.append('h2d_direct'+cd+pg)
        list_int_bh_dir.append('bh2d_direct'+cd+pg)
    for pg in p_optsl2:
        list_int_lap_dir.append('r2d_direct'+cd+pg)
        list_int_lap_dir.append('l2d_direct'+cd+pg)
        list_int_lap_dir.append('c2d_direct'+cd+pg)
list_int_bh.append('bhfmm2dwrap_guru')

ext_helm = Extension(
    name='fmm2dpy.hfmm2d_fortran',
    sources=['../src/helmholtz/'+item for item in list_helm]+['../src/common/'+item for item in list_common],
    f2py_options=['only:']+list_int_helm+list_int_helm_vec+list_int_helm_dir+[':'],
    extra_f90_compile_args=FFLAGS,
    extra_f77_compile_args=FFLAGS,
    extra_link_args=FLIBS
)

ext_lap = Extension(
    name='fmm2dpy.lfmm2d_fortran',
    sources=['../src/laplace/'+item for item in list_lap]+['../src/common/'+item for item in list_common],
    f2py_options=['only:']+list_int_lap+list_int_lap_vec+list_int_lap_dir+[':'],
    extra_f90_compile_args=FFLAGS,
    extra_f77_compile_args=FFLAGS,
    extra_link_args=FLIBS
)

ext_bh = Extension(
    name='fmm2dpy.bhfmm2d_fortran',
    sources=['../src/biharmonic/'+item for item in list_bh]+['../src/common/'+item for item in list_common],
    f2py_options=['only:']+list_int_bh+list_int_bh_dir+[':'],
    extra_f90_compile_args=FFLAGS,
    extra_f77_compile_args=FFLAGS,
    extra_link_args=FLIBS
)

ext_st = Extension(
    name='fmm2dpy.stfmm2d_fortran',
    sources=['../src/stokes/stfmm2d.f']+['../src/stokes/stokkernels2d.f'],
    f2py_options=['only:']+['stfmm2d']+['st2ddirectstokg']+['st2ddirectstokstrsg']+[':'],
    extra_link_args=FLIBS
)


## TODO: fill in the info below
setup(
    name=pkg_name,
    python_requires='>=3.0.0',
    version="1.1.0",
    author="Travis Askham, Zydrunas Gimbutas, Leslie Greengard, Libin Lu, Michael O'Neil, Manas Rachh, and Vladimir Rokhlin",
    author_email="mrachh@flatironinstitute.org",
    description="This pacakge contains basic routines for Laplace and Helmholtz fast multipole methods in two dimensions",
    long_description=open('../README.md').read(),
    long_description_content_type='text/markdown',
    url="https://github.com/flatironinstitute/fmm2d",
    packages=['fmm2dpy'],
    install_requires=[
        "numpy",
        "pytest"
    ],
    ext_modules=[ext_helm,ext_lap,ext_bh,ext_st],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ]
)
