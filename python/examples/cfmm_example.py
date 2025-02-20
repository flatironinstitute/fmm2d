#!/usr/bin/env python

import fmm2dpy as fmm
import numpy as np


#
#  This is a sample code to demonstrate how to use
#  the fmm libraries
#

# sample with one density, sources to sources,
# charge interactions, and potential only
#
n = 2000000
nd = 1
sources = np.random.uniform(0,1,(2,n))
eps = 10**(-5)

charges = np.random.uniform(0,1,n) + 1j*np.random.uniform(0,1,n)
out = fmm.lfmm2d(eps=eps,sources=sources,charges=charges,pg=1)


# sample with a vector of densities, sources to 
# sources and targets, dipole interactions, 
# potential and gradietns

nd = 3
nt = 1870
targ = np.random.uniform(0,1,(2,nt))
dipstr = np.random.uniform(0,1,(nd,n)) +  1j*np.random.uniform(0,1,n)
out2 = fmm.cfmm2d(eps=eps,sources=sources,dipstr=dipstr,\
    targets=targ,nd=nd,pg=2,pgt=2)
