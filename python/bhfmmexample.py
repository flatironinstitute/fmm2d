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
n = 100000
nd = 1
sources = np.random.uniform(0,1,(2,n))
eps = 10**(-5)

charges = np.random.uniform(0,1,(2,n)) + 1j*np.random.uniform(0,1,(2,n))
dipoles = np.random.uniform(0,1,(3,n)) + 1j*np.random.uniform(0,1,(3,n))
out = fmm.bhfmm2d(eps=eps,sources=sources,charges=charges,dipoles=dipoles,pg=1)


# sample with a vector of densities, sources to 
# sources and targets, dipole interactions, 
# potential and gradietns

nd = 3
nt = 1870
targ = np.random.uniform(0,1,(2,nt))
dipoles = np.random.uniform(0,1,(nd,3,n)) +  1j*np.random.uniform(0,1,(nd,3,n))
out2 = fmm.bhfmm2d(eps=eps,sources=sources,dipoles=dipoles,\
    targets=targ,nd=nd,pg=2,pgt=2)
