import fmm2d as fmm
import numpy as np


n = 2000
sources = np.random.uniform(0,1,(2,n))

m = 1000
targets = np.random.uniform(0,1,(2,m))

zk = 1.1+ 1j*0
charges = np.random.uniform(0,1,n) + 1j*np.random.uniform(0,1,n)

eps = 10**(-5)

pottarg,ier = fmm.hfmm2d_t_c_p(eps,zk,sources,charges,targets)

print(pottarg[1:3])
