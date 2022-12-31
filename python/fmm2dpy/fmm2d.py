from . import hfmm2d_fortran as hfmm
from . import lfmm2d_fortran as lfmm
from . import bhfmm2d_fortran as bhfmm
import numpy as np
import numpy.linalg as la


class Output():
    pot = None
    grad = None
    hess = None
    pottarg = None
    gradtarg = None
    hesstarg = None
    ier = 0

def hfmm2d(*,eps,zk,sources,charges=None,dipstr=None,dipvec=None,
          targets=None,pg=0,pgt=0,nd=1):
    r"""
      This subroutine computes the N-body Helmholtz interactions
      in three dimensions where the interaction kernel is given by e^{ikr}/r 
      and its gradients. 

      .. math::

          u(x) = \sum_{j=1}^{N} c_{j} H_{0}^{(1)}(k \|x-x_{j}\|) - d_{j} v_{j}\cdot \\nabla \left( H_{0}^{(1)}(k \|x-x_{j}\|) \\right)  \, ,

      where $c_{j}$ are the charge densities, $d_{j}$ are the dipole densities, 
      $v_{j}$ are the dipole orientation vectors, and 
      $x_{j}$ are the source locations.

      When $x=x_{m}$, the term corresponding to $x_{m}$ is dropped from the
      sum

      Args:
        eps (float): precision requested
        zk (complex): Helmholtz parameter 
        sources (float(2,n)): source locations ($x_{j}$)
        charges (complex(nd,n) or complex(n)): charge densities ($c_{j}$)
        dipstr (complex(nd,n) or complex(n)): dipole densities ($d_{j}$)
        dipvec (float(nd,2,n) or float(2,n)): dipole orientation vectors ($v_{j}$)
        targets (float(2,nt)): target locations (x)
        pg (integer): source eval flag. Potential at sources evaluated if pg = 1. Potenial and gradient at sources evaluated if pg=2. Potential, gradient and hessian evaluated at sources if pg=3
        pgt (integer): target eval flag. Potential at targets evaluated if pgt = 1. Potenial and gradient at targets evaluated if pgt=2. Potential, gradient and hessian evaluated at targets if pgt=3
        nd (integer): number of densities

      Returns:
        Returns an object of type Output (out) with the following variables

        out.pot: potential at source locations if requested
        out.grad: gradient at source locations if requested
        out.hess: hessian at source locations if requested
        out.pottarg: potential at target locations if requested
        out.gradtarg: gradient at target locations if requested
        out.hesstarg: hessian at target locations if requested


      Example:
        see hfmmexample.py
    r"""
    
    out = Output()
    assert sources.shape[0] == 2, "The first dimension of sources must be 2"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1

    ifcharge = 0
    ifdipole = 0
    iftarg = 0
    if(pg == 0 and pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
        ifcharge = 1
    if(dipvec is not None or dipstr is not None):
        if nd == 1 and ns>1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            assert dipvec.shape[0] == 2 and dipvec.shape[1] == ns, "dipole vectors must be of shape [2,number of sources]"
        if nd == 1 and ns==1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            assert dipvec.shape[0] == 2, "dipole vectors must be of shape [2,number of sources]"
        if nd>1:
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 2 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,2,ns] where nd is number of densities, and ns is number of sources"
            assert dipstr.shape[0] == nd and dipstr.shape[1]==ns, "Dipole strengths must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1
    if(targets is not None):
        assert targets.shape[0] == 2, "The first dimension of targets must be 2"
        iftarg = 1
    if(iftarg == 0 or pgt != 1 or pgt !=2):
        if(pg == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.ier = hfmm.hfmm2d_s_c_p_vec(eps,zk,sources,charges,nd)
            if(nd == 1):
                out.pot,out.ier = hfmm.hfmm2d_s_c_p(eps,zk,sources,charges)
        if(pg == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.ier = hfmm.hfmm2d_s_c_g_vec(eps,zk,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = hfmm.hfmm2d_s_c_g(eps,zk,sources,charges)

        if(pg == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.ier = hfmm.hfmm2d_s_d_p_vec(eps,zk,sources,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.ier = hfmm.hfmm2d_s_d_p(eps,zk,sources,dipstr,dipvec)
        if(pg == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.ier = hfmm.hfmm2d_s_d_g_vec(eps,zk,sources,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = hfmm.hfmm2d_s_d_g(eps,zk,sources,dipstr,dipvec)

        if(pg == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.ier = hfmm.hfmm2d_s_cd_p_vec(eps,zk,sources,charges,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.ier = hfmm.hfmm2d_s_cd_p(eps,zk,sources,charges,dipstr,dipvec)
        if(pg == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.ier = hfmm.hfmm2d_s_cd_g_vec(eps,zk,sources,charges,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = hfmm.hfmm2d_s_cd_g(eps,zk,sources,charges,dipstr,dipvec)
    
    if(pg !=1 and pg !=2 and targets is not None):
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.ier = hfmm.hfmm2d_t_c_p_vec(eps,zk,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = hfmm.hfmm2d_t_c_p(eps,zk,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = hfmm.hfmm2d_t_c_g_vec(eps,zk,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = hfmm.hfmm2d_t_c_g(eps,zk,sources,charges,targets)
        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.ier = hfmm.hfmm2d_t_d_p_vec(eps,zk,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = hfmm.hfmm2d_t_d_p(eps,zk,sources,dipstr,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = hfmm.hfmm2d_t_d_g_vec(eps,zk,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = hfmm.hfmm2d_t_d_g(eps,zk,sources,dipstr,dipvec,targets)
        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.ier = hfmm.hfmm2d_t_cd_p_vec(eps,zk,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = hfmm.hfmm2d_t_cd_p(eps,zk,sources,charges,dipstr,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = hfmm.hfmm2d_t_cd_g_vec(eps,zk,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = hfmm.hfmm2d_t_cd_g(eps,zk,sources,charges,dipstr,dipvec,targets)
    

    if((pg == 1 or pg == 2) and targets is not None):
        assert pg == pgt, "if both potential or potential at gradient are requested at sources and targets, then the same pg must be equal to pgt"
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = hfmm.hfmm2d_st_c_p_vec(eps,zk,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = hfmm.hfmm2d_st_c_p(eps,zk,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = hfmm.hfmm2d_st_c_g_vec(eps,zk,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = hfmm.hfmm2d_st_c_g(eps,zk,sources,charges,targets)
        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = hfmm.hfmm2d_st_d_p_vec(eps,zk,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = hfmm.hfmm2d_st_d_p(eps,zk,sources,dipstr,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = hfmm.hfmm2d_st_d_g_vec(eps,zk,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = hfmm.hfmm2d_st_d_g(eps,zk,sources,dipstr,dipvec,targets)
        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = hfmm.hfmm2d_st_cd_p_vec(eps,zk,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = hfmm.hfmm2d_st_cd_p(eps,zk,sources,charges,dipstr,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = hfmm.hfmm2d_st_cd_g_vec(eps,zk,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = hfmm.hfmm2d_st_cd_g(eps,zk,sources,charges,dipstr,dipvec,targets)

    return out


def rfmm2d(*,eps,sources,charges=None,dipstr=None,dipvec=None,
          targets=None,pg=0,pgt=0,nd=1):
    r"""
      This subroutine computes the N-body Laplace interactions
      in two dimensions where the interaction kernel is given by log(r) 
      and its gradients. 


      .. math:: 

          u(x) = \sum_{j=1}^{N} c_{j} * log\(\|x-x_{j}\|\) + d_{j}v_{j} \cdot \\nabla( log(\|x-x_{j}\|) )  \, ,

      where $c_{j}$ are the charge densities, $d_{j}$ are the dipole strengths,
      $v_{j}$ are the dipole orientation vectors, and 
      $x_{j}$ are the source locations.

      When $x=x_{m}$, the term corresponding to $x_{m}$ is dropped from the
      sum


      Args:
        eps: float
             precision requested

        sources: float(2,n)   
               source locations (x_{j})
        charges: float(nd,n) or float(n)
               charge densities (c_{j})
        dipstr: float(nd,n) or float(n)
               dipole densities (d_{j})
        dipvec: float(nd,2,n) or float(2,n)
               dipole orientation vectors (v_{j})
        targets: float(2,nt)
                target locations (x)
        pg:  integer
               source eval flag
               potential at sources evaluated if pg = 1
               potenial and gradient at sources evaluated if pg=2
               potential, gradient and hessian at sources evaluated if pg=3

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
               potential, gradient and hessian at targets evaluated if pgt=3
        
        nd:   integer
               number of densities

      Returns:
        out.pot: potential at source locations if requested
        out.grad: gradient at source locations if requested
        out.hess: hessian at source locations if requested
        out.pottarg: potential at target locations if requested
        out.gradtarg: gradient at target locations if requested
        out.hesstarg: hessian at target locations if requested
      
      Example:
        see rfmmexample.py
    r"""

    out = Output()
    assert sources.shape[0] == 2, "The first dimension of sources must be 2"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    ifcharge = 0
    ifdipole = 0
    iftarg = 0
    if(pg == 0 and pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipvec is not None or dipstr is not None):
        if nd == 1 and ns>1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            assert dipvec.shape[0] == 2 and dipvec.shape[1] == ns, "dipole vectors must be of shape [2,number of sources]"
        if nd == 1 and ns==1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            assert dipvec.shape[0] == 2, "dipole vectors must be of shape [2,number of sources]"
        if nd>1:
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 2 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,2,ns] where nd is number of densities, and ns is number of sources"
            assert dipstr.shape[0] == nd and dipstr.shape[1]==ns, "Dipole strengths must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1
    if(targets is not None):
        assert targets.shape[0] == 2, "The first dimension of targets must be 2"
        iftarg = 1
#
# sources -> sources routines
#
    if(iftarg == 0 or pgt != 1 or pgt !=2 or pgt !=3):
        if(pg == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.ier = lfmm.rfmm2d_s_c_p_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.ier = lfmm.rfmm2d_s_c_p(eps,sources,charges)
        if(pg == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.ier = lfmm.rfmm2d_s_c_g_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = lfmm.rfmm2d_s_c_g(eps,sources,charges)
        if(pg == 3 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.rfmm2d_s_c_h_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.rfmm2d_s_c_h(eps,sources,charges)


        if(pg == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.ier = lfmm.rfmm2d_s_d_p_vec(eps,sources,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.ier = lfmm.rfmm2d_s_d_p(eps,sources,dipstr,dipvec)
        if(pg == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.ier = lfmm.rfmm2d_s_d_g_vec(eps,sources,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = lfmm.rfmm2d_s_d_g(eps,sources,dipstr,dipvec)
        if(pg == 3 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.rfmm2d_s_d_h_vec(eps,sources,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.rfmm2d_s_d_h(eps,sources,dipstr,dipvec)


        if(pg == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.ier = lfmm.rfmm2d_s_cd_p_vec(eps,sources,charges,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.ier = lfmm.rfmm2d_s_cd_p(eps,sources,charges,dipstr,dipvec)
        if(pg == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.ier = lfmm.rfmm2d_s_cd_g_vec(eps,sources,charges,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = lfmm.rfmm2d_s_cd_g(eps,sources,charges,dipstr,dipvec)
        if(pg == 3 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.rfmm2d_s_cd_h_vec(eps,sources,charges,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.rfmm2d_s_cd_h(eps,sources,charges,dipstr,dipvec)

#
# sources -> targets routines
#


    if(pg !=1 and pg !=2 and pg !=3 and targets is not None):
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.ier = lfmm.rfmm2d_t_c_p_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = lfmm.rfmm2d_t_c_p(eps,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.rfmm2d_t_c_g_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.rfmm2d_t_c_g(eps,sources,charges,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.rfmm2d_t_c_h_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.rfmm2d_t_c_h(eps,sources,charges,targets)


        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.ier = lfmm.rfmm2d_t_d_p_vec(eps,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = lfmm.rfmm2d_t_d_p(eps,sources,dipstr,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.rfmm2d_t_d_g_vec(eps,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.rfmm2d_t_d_g(eps,sources,dipstr,dipvec,targets)
        if(pgt == 3 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.rfmm2d_t_d_h_vec(eps,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.rfmm2d_t_d_h(eps,sources,dipstr,dipvec,targets)


        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.ier = lfmm.rfmm2d_t_cd_p_vec(eps,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = lfmm.rfmm2d_t_cd_p(eps,sources,charges,dipstr,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.rfmm2d_t_cd_g_vec(eps,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.rfmm2d_t_cd_g(eps,sources,charges,dipstr,dipvec,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.rfmm2d_t_cd_h_vec(eps,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.rfmm2d_t_cd_h(eps,sources,charges,dipstr,dipvec,targets)
    
#
# sources to sources + targets
#
    if((pg == 1 or pg == 2 or pg == 3) and targets is not None):
        assert pg == pgt, "if output is requested at both sources and targets, then the same pg must be equal to pgt"
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = lfmm.rfmm2d_st_c_p_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = lfmm.rfmm2d_st_c_p(eps,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.rfmm2d_st_c_g_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.rfmm2d_st_c_g(eps,sources,charges,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.rfmm2d_st_c_h_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.rfmm2d_st_c_h(eps,sources,charges,targets)


        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = lfmm.rfmm2d_st_d_p_vec(eps,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = lfmm.rfmm2d_st_d_p(eps,sources,dipstr,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.rfmm2d_st_d_g_vec(eps,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.rfmm2d_st_d_g(eps,sources,dipstr,dipvec,targets)
        if(pgt == 3 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.rfmm2d_st_d_h_vec(eps,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.rfmm2d_st_d_h(eps,sources,dipstr,dipvec,targets)


        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = lfmm.rfmm2d_st_cd_p_vec(eps,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = lfmm.rfmm2d_st_cd_p(eps,sources,charges,dipstr,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.rfmm2d_st_cd_g_vec(eps,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.rfmm2d_st_cd_g(eps,sources,charges,dipstr,dipvec,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.rfmm2d_st_cd_h_vec(eps,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.rfmm2d_st_cd_h(eps,sources,charges,dipstr,dipvec,targets)

    return out


def lfmm2d(*,eps,sources,charges=None,dipstr=None,dipvec=None,
          targets=None,pg=0,pgt=0,nd=1):
    r"""
      This subroutine computes the N-body Laplace interactions
      in two dimensions where the interaction kernel is given by log(r) 
      and its gradients. 


      .. math:: 

          u(x) = \sum_{j=1}^{N} c_{j} * log\(\|x-x_{j}\|\) + d_{j}v_{j} \cdot \\nabla( log(\|x-x_{j}\|) )  \, ,

      where $c_{j}$ are the charge densities, $d_{j}$ are the dipole strengths,
      $v_{j}$ are the dipole orientation vectors, and 
      $x_{j}$ are the source locations.

      When $x=x_{m}$, the term corresponding to $x_{m}$ is dropped from the
      sum


      Args:
        eps: float
             precision requested

        sources: float(2,n)   
               source locations (x_{j})
        charges: complex(nd,n) or complex(n)
               charge densities (c_{j})
        dipstr: complex(nd,n) or complex(n)
               dipole densities (d_{j})
        dipvec: float(nd,2,n) or float(2,n)
               dipole orientation vectors (v_{j})
        targets: float(2,nt)
                target locations (x)
        pg:  integer
               source eval flag
               potential at sources evaluated if pg = 1
               potenial and gradient at sources evaluated if pg=2
               potential, gradient and hessian at sources evaluated if pg=3

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
               potential, gradient and hessian at targets evaluated if pgt=3
        
        nd:   integer
               number of densities

      Returns:
        out.pot: potential at source locations if requested
        out.grad: gradient at source locations if requested
        out.hess: hessian at source locations if requested
        out.pottarg: potential at target locations if requested
        out.gradtarg: gradient at target locations if requested
        out.hesstarg: hessian at target locations if requested
      
      Example:
        see lfmmexample.py
    r"""

    out = Output()
    assert sources.shape[0] == 2, "The first dimension of sources must be 2"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    ifcharge = 0
    ifdipole = 0
    iftarg = 0
    if(pg == 0 and pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipvec is not None or dipstr is not None):
        if nd == 1 and ns>1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            assert dipvec.shape[0] == 2 and dipvec.shape[1] == ns, "dipole vectors must be of shape [2,number of sources]"
        if nd == 1 and ns==1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            assert dipvec.shape[0] == 2, "dipole vectors must be of shape [2,number of sources]"
        if nd>1:
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 2 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,2,ns] where nd is number of densities, and ns is number of sources"
            assert dipstr.shape[0] == nd and dipstr.shape[1]==ns, "Dipole strengths must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1
    if(targets is not None):
        assert targets.shape[0] == 2, "The first dimension of targets must be 2"
        iftarg = 1
#
# sources -> sources routines
#
    if(iftarg == 0 or pgt != 1 or pgt !=2 or pgt !=3):
        if(pg == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.ier = lfmm.lfmm2d_s_c_p_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.ier = lfmm.lfmm2d_s_c_p(eps,sources,charges)
        if(pg == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.ier = lfmm.lfmm2d_s_c_g_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = lfmm.lfmm2d_s_c_g(eps,sources,charges)
        if(pg == 3 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.lfmm2d_s_c_h_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.lfmm2d_s_c_h(eps,sources,charges)


        if(pg == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.ier = lfmm.lfmm2d_s_d_p_vec(eps,sources,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.ier = lfmm.lfmm2d_s_d_p(eps,sources,dipstr,dipvec)
        if(pg == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.ier = lfmm.lfmm2d_s_d_g_vec(eps,sources,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = lfmm.lfmm2d_s_d_g(eps,sources,dipstr,dipvec)
        if(pg == 3 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.lfmm2d_s_d_h_vec(eps,sources,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.lfmm2d_s_d_h(eps,sources,dipstr,dipvec)


        if(pg == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.ier = lfmm.lfmm2d_s_cd_p_vec(eps,sources,charges,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.ier = lfmm.lfmm2d_s_cd_p(eps,sources,charges,dipstr,dipvec)
        if(pg == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.ier = lfmm.lfmm2d_s_cd_g_vec(eps,sources,charges,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = lfmm.lfmm2d_s_cd_g(eps,sources,charges,dipstr,dipvec)
        if(pg == 3 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.lfmm2d_s_cd_h_vec(eps,sources,charges,dipstr,dipvec,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.lfmm2d_s_cd_h(eps,sources,charges,dipstr,dipvec)

#
# sources -> targets routines
#


    if(pg !=1 and pg !=2 and pg !=3 and targets is not None):
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.ier = lfmm.lfmm2d_t_c_p_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = lfmm.lfmm2d_t_c_p(eps,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.lfmm2d_t_c_g_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.lfmm2d_t_c_g(eps,sources,charges,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm2d_t_c_h_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm2d_t_c_h(eps,sources,charges,targets)


        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.ier = lfmm.lfmm2d_t_d_p_vec(eps,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = lfmm.lfmm2d_t_d_p(eps,sources,dipstr,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.lfmm2d_t_d_g_vec(eps,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.lfmm2d_t_d_g(eps,sources,dipstr,dipvec,targets)
        if(pgt == 3 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm2d_t_d_h_vec(eps,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm2d_t_d_h(eps,sources,dipstr,dipvec,targets)


        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.ier = lfmm.lfmm2d_t_cd_p_vec(eps,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = lfmm.lfmm2d_t_cd_p(eps,sources,charges,dipstr,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.lfmm2d_t_cd_g_vec(eps,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.lfmm2d_t_cd_g(eps,sources,charges,dipstr,dipvec,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm2d_t_cd_h_vec(eps,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm2d_t_cd_h(eps,sources,charges,dipstr,dipvec,targets)
    
#
# sources to sources + targets
#
    if((pg == 1 or pg == 2 or pg == 3) and targets is not None):
        assert pg == pgt, "if output is requested at both sources and targets, then the same pg must be equal to pgt"
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = lfmm.lfmm2d_st_c_p_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = lfmm.lfmm2d_st_c_p(eps,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.lfmm2d_st_c_g_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.lfmm2d_st_c_g(eps,sources,charges,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm2d_st_c_h_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm2d_st_c_h(eps,sources,charges,targets)


        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = lfmm.lfmm2d_st_d_p_vec(eps,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = lfmm.lfmm2d_st_d_p(eps,sources,dipstr,dipvec,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.lfmm2d_st_d_g_vec(eps,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.lfmm2d_st_d_g(eps,sources,dipstr,dipvec,targets)
        if(pgt == 3 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm2d_st_d_h_vec(eps,sources,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm2d_st_d_h(eps,sources,dipstr,dipvec,targets)


        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = lfmm.lfmm2d_st_cd_p_vec(eps,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = lfmm.lfmm2d_st_cd_p(eps,sources,charges,dipstr,dipvec,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.lfmm2d_st_cd_g_vec(eps,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.lfmm2d_st_cd_g(eps,sources,charges,dipstr,dipvec,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm2d_st_cd_h_vec(eps,sources,charges,dipstr,dipvec,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.lfmm2d_st_cd_h(eps,sources,charges,dipstr,dipvec,targets)

    return out


def cfmm2d(*,eps,sources,charges=None,dipstr=None,
          targets=None,pg=0,pgt=0,nd=1):
    r"""
      This subroutine computes the N-body Laplace interactions with Cauchy kernel
      in two dimensions where the interaction kernel is given by log(r)
      and its gradients. 


      .. math:: 

          u(x) = \sum_{j=1}^{N} c_{j} * log\(\|x-x_{j}\|\) + d_{j}/(x-x_{j})  \, ,

      where $c_{j}$ are the charge densities, $d_{j}$ are the dipole strengths,
      and $x_{j}$ are the source locations.

      When $x=x_{m}$, the term corresponding to $x_{m}$ is dropped from the
      sum


      Args:
        eps: float
             precision requested

        sources: float(2,n)   
               source locations (x_{j})
        charges: complex(nd,n) or complex(n)
               charge densities (c_{j})
        dipstr: complex(nd,n) or complex(n)
               dipole densities (d_{j})
        targets: float(2,nt)
                target locations (x)
        pg:  integer
               source eval flag
               potential at sources evaluated if pg = 1
               potenial and gradient at sources evaluated if pg=2
               potential, gradient and hessian at sources evaluated if pg=3

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
               potential, gradient and hessian at targets evaluated if pgt=3
        
        nd:   integer
               number of densities

      Returns:
        out.pot: potential at source locations if requested
        out.grad: gradient at source locations if requested
        out.hess: hessian at source locations if requested
        out.pottarg: potential at target locations if requested
        out.gradtarg: gradient at target locations if requested
        out.hesstarg: hessian at target locations if requested
      
      Example:
        see cfmmexample.py
    r"""

    out = Output()
    assert sources.shape[0] == 2, "The first dimension of sources must be 2"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    ifcharge = 0
    ifdipole = 0
    iftarg = 0
    if(pg == 0 and pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
            charges = charges.reshape(1,ns)
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipstr is not None):
        if nd == 1 and ns>1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
        if nd == 1 and ns==1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
        if nd>1:
            assert dipstr.shape[0] == nd and dipstr.shape[1]==ns, "Dipole strengths must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1
    if(targets is not None):
        assert targets.shape[0] == 2, "The first dimension of targets must be 2"
        iftarg = 1
#
# sources -> sources routines
#
    if(iftarg == 0 or pgt != 1 or pgt !=2 or pgt !=3):
        if(pg == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.ier = lfmm.cfmm2d_s_c_p_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.ier = lfmm.cfmm2d_s_c_p(eps,sources,charges)
        if(pg == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.ier = lfmm.cfmm2d_s_c_g_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = lfmm.cfmm2d_s_c_g(eps,sources,charges)
        if(pg == 3 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.cfmm2d_s_c_h_vec(eps,sources,charges,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.cfmm2d_s_c_h(eps,sources,charges)


        if(pg == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.ier = lfmm.cfmm2d_s_d_p_vec(eps,sources,dipstr,nd)
            if(nd == 1):
                out.pot,out.ier = lfmm.cfmm2d_s_d_p(eps,sources,dipstr)
        if(pg == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.ier = lfmm.cfmm2d_s_d_g_vec(eps,sources,dipstr,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = lfmm.cfmm2d_s_d_g(eps,sources,dipstr)
        if(pg == 3 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.cfmm2d_s_d_h_vec(eps,sources,dipstr,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.cfmm2d_s_d_h(eps,sources,dipstr)


        if(pg == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.ier = lfmm.cfmm2d_s_cd_p_vec(eps,sources,charges,dipstr,nd)
            if(nd == 1):
                out.pot,out.ier = lfmm.cfmm2d_s_cd_p(eps,sources,charges,dipstr)
        if(pg == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.ier = lfmm.cfmm2d_s_cd_g_vec(eps,sources,charges,dipstr,nd)
            if(nd == 1):
                out.pot,out.grad,out.ier = lfmm.cfmm2d_s_cd_g(eps,sources,charges,dipstr)
        if(pg == 3 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.cfmm2d_s_cd_h_vec(eps,sources,charges,dipstr,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.ier = lfmm.cfmm2d_s_cd_h(eps,sources,charges,dipstr)

#
# sources -> targets routines
#


    if(pg !=1 and pg !=2 and pg !=3 and targets is not None):
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.ier = lfmm.cfmm2d_t_c_p_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = lfmm.cfmm2d_t_c_p(eps,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.cfmm2d_t_c_g_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.cfmm2d_t_c_g(eps,sources,charges,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.cfmm2d_t_c_h_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.cfmm2d_t_c_h(eps,sources,charges,targets)


        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.ier = lfmm.cfmm2d_t_d_p_vec(eps,sources,dipstr,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = lfmm.cfmm2d_t_d_p(eps,sources,dipstr,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.cfmm2d_t_d_g_vec(eps,sources,dipstr,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.cfmm2d_t_d_g(eps,sources,dipstr,targets)
        if(pgt == 3 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.cfmm2d_t_d_h_vec(eps,sources,dipstr,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.cfmm2d_t_d_h(eps,sources,dipstr,targets)


        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.ier = lfmm.cfmm2d_t_cd_p_vec(eps,sources,charges,dipstr,targets,nd)
            if(nd == 1):
                out.pottarg,out.ier = lfmm.cfmm2d_t_cd_p(eps,sources,charges,dipstr,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.cfmm2d_t_cd_g_vec(eps,sources,charges,dipstr,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.ier = lfmm.cfmm2d_t_cd_g(eps,sources,charges,dipstr,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.cfmm2d_t_cd_h_vec(eps,sources,charges,dipstr,targets,nd)
            if(nd == 1):
                out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.cfmm2d_t_cd_h(eps,sources,charges,dipstr,targets)
    
#
# sources to sources + targets
#
    if((pg == 1 or pg == 2 or pg == 3) and targets is not None):
        assert pg == pgt, "if output is requested at both sources and targets, then the same pg must be equal to pgt"
        if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = lfmm.cfmm2d_st_c_p_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = lfmm.cfmm2d_st_c_p(eps,sources,charges,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.cfmm2d_st_c_g_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.cfmm2d_st_c_g(eps,sources,charges,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 0):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.cfmm2d_st_c_h_vec(eps,sources,charges,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.cfmm2d_st_c_h(eps,sources,charges,targets)


        if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = lfmm.cfmm2d_st_d_p_vec(eps,sources,dipstr,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = lfmm.cfmm2d_st_d_p(eps,sources,dipstr,targets)
        if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.cfmm2d_st_d_g_vec(eps,sources,dipstr,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.cfmm2d_st_d_g(eps,sources,dipstr,targets)
        if(pgt == 3 and ifcharge == 0 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.cfmm2d_st_d_h_vec(eps,sources,dipstr,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.cfmm2d_st_d_h(eps,sources,dipstr,targets)


        if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.pottarg,out.ier = lfmm.cfmm2d_st_cd_p_vec(eps,sources,charges,dipstr,targets,nd)
            if(nd == 1):
                out.pot,out.pottarg,out.ier = lfmm.cfmm2d_st_cd_p(eps,sources,charges,dipstr,targets)
        if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.cfmm2d_st_cd_g_vec(eps,sources,charges,dipstr,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.pottarg,out.gradtarg,out.ier = lfmm.cfmm2d_st_cd_g(eps,sources,charges,dipstr,targets)
        if(pgt == 3 and ifcharge == 1 and ifdipole == 1):
            if(nd > 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.cfmm2d_st_cd_h_vec(eps,sources,charges,dipstr,targets,nd)
            if(nd == 1):
                out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = lfmm.cfmm2d_st_cd_h(eps,sources,charges,dipstr,targets)

    return out


def bhfmm2d(*,eps,sources,charges=None,dipoles=None,
          targets=None,pg=0,pgt=0,nd=1):
    r"""
      This subroutine computes the N-body biharmonic interactions 
      in two dimensions where the interaction kernel is related to the
      biharmonic greens function r^2 log (r) and its derivatives


      .. math:: 

          u(x) = \sum_{j=1}^{N} c_{j} * log\(\|x-x_{j}\|\) + 
          \overline{c}_{j} (x-x_{j})/(\overline{x-x_{j}) + d_{j,1}/(x-x_{j}) + 
          d_{j,2}/(\overline{x-x_{j}}) - 
          \overline{d_{j,1}} (x-x_{j})/(\overline{x-x_{j}})^2\, ,

      where $c_{j}$ are the charge densities, $d_{j,1}$, $d_{j,2}$ are the dipole strengths,
      and $x_{j}$ are the source locations.

      When $x=x_{m}$, the term corresponding to $x_{m}$ is dropped from the
      sum


      Args:
        eps: float
             precision requested

        sources: float(2,n)   
               source locations (x_{j})
        charges: complex(nd,n) or complex(n)
               charge densities (c_{j})
        dipoles: complex(nd,2,n) or complex(2,n)
               dipole densities (d_{j,1}, d_{j,2})
        targets: float(2,nt)
                target locations (x)
        pg:  integer
               source eval flag
               potential at sources evaluated if pg = 1
               potenial and gradient at sources evaluated if pg=2

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
        
        nd:   integer
               number of densities

      Returns:
        out.pot: potential at source locations if requested
        out.grad: gradient at source locations if requested
        out.pottarg: potential at target locations if requested
        out.gradtarg: gradient at target locations if requested
      
      Example:
        see bhfmmexample.py
    r"""

    out = Output()
    assert sources.shape[0] == 2, "The first dimension of sources must be 2"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    ifcharge = 0
    ifdipole = 0
    iftarg = 0
    if(pg == 0 and pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
            charges = charges.reshape(1,ns)
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
        charges = charges.reshape([nd,ns])
    else:
       charges = np.zeros([nd,ns],dtype='complex') 

    if(dipoles is not None):
        if nd == 1 and ns>1:
            assert dipoles.shape[0] == 2 and dipoles.shape[1] == ns, "Dipole must of shape [2, number of sources]"
        if nd == 1 and ns==1:
            assert dipoles.shape[0] == 2, "Dipole must of shape [2, number of sources]"
        if nd>1:
            assert dipoles.shape[0] == nd and dipoles.shape[1] == 2 and dipoles.shape[2]==ns, "Dipole must of shape [nd,2, number of sources]"
        dipoles = dipoles.reshape([nd,2,ns])
        ifdipole = 1
    else:
        dipoles = np.zeros([nd,2,ns],dtype='complex')
    if(targets is not None):
        assert targets.shape[0] == 2, "The first dimension of targets must be 2"
        iftarg = 1
        nt = np.shape(targets)[1]
    else:
        targets = np.zeros([2,0],dtype='double')
        nt = 0 
    iper = 0
    out.pot,out.grad,out.hess,out.pottarg,out.gradtarg,out.hesstarg,out.ier = bhfmm.bhfmm2dwrap_guru(eps,sources,ifcharge,charges,ifdipole,dipoles,iper,pg,targets,pgt)
    out.hess = None
    out.hesstarg = None
    if(nd == 1):
        if(pgt>0):
            out.pottarg = out.pottarg.reshape(nt,)
        if(pgt==2):
            out.gradtarg = out.gradtarg.reshape(2,nt)
        if(pg>0):
            out.pot = out.pot.reshape(ns,)
        if(pg==2):
            out.grad = out.grad.reshape(2,ns)

    if(pg<2):
        out.grad = None
    if(pgt<2):
        out.gradtarg = None

    return out

def h2ddir(*,zk,sources,targets,charges=None,dipstr=None,dipvec=None,
          pgt=0,nd=1,thresh=1e-16):
    r"""
      This subroutine computes the N-body Helmholtz interactions
      in three dimensions where the interaction kernel is given by $H_{0}^{(1)(kr)$ 
      and its gradients. 


      .. math::


          u(x) = \sum_{j=1}^{N} c_{j} H_{0}^{(1)}(k \|x-x_{j}\|) - d_{j} v_{j}\cdot \\nabla \left( H_{0}^{(1)}(k \|x-x_{j}\|) \\right)  \, ,

      where $c_{j}$ are the charge densities, $d_{j}$ are dipole densities 
      $v_{j}$ are the dipole orientation vectors, and 
      $x_{j}$ are the source locations.

      When |x-x_{m}| \leq thresh, the term corresponding to $x_{m}$ is dropped from the
      sum


      Args:
        eps: float   
               precision requested
        zk: complex
               Helmholtz parameter - k
        sources: float(2,n)   
               source locations (x_{j})
        charges: complex(nd,n) or complex(n)
                charge densities (c_{j})
        dipstr: complex(nd,n) or complex(n)
                dipole densities (d_{j})
        dipole orientation vectors: float(nd,2,n) or complex(2,n)
                dipole orientation vectors (v_{j})
        targets: float(2,nt)
                target locations (x)

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
               potential, gradient and hessians at targets evaluated if pgt=3
        
        nd:   integer
               number of densities
        thresh: contribution of source x_i, at location x ignored if |x-x_i|<=thresh

      Returns:
        out.pottarg  - potential at target locations if requested
        out.gradtarg - gradient at target locations if requested
        out.hesstarg - hessian at target locations if requested
              
      Example:
        see hfmmexample.py
    r"""

    out = Output()
    assert sources.shape[0] == 2, "The first dimension of sources must be 2"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    ifcharge = 0
    ifdipole = 0
    if(pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
            charges = charges.reshape(1,ns)
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
        ifcharge = 1
    if(dipvec is not None or dipstr is not None):
        if nd == 1 and ns>1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            assert dipvec.shape[0] == 2 and dipvec.shape[1] == ns, "dipole vectors must be of shape [2,number of sources]"
            dipvec = dipvec.reshape(1,2,ns)
            dipstr = dipstr.reshape(1,ns)
        if nd == 1 and ns==1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            assert dipvec.shape[0] == 2, "dipole vectors must be of shape [2,number of sources]"
            dipvec = dipvec.reshape(1,2,ns)
            dipstr = dipstr.reshape(1,ns)
        if nd>1:
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 2 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,2,ns] where nd is number of densities, and ns is number of sources"
            assert dipstr.shape[0] == nd and dipstr.shape[1]==ns, "Dipole strengths must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1

    assert targets.shape[0] == 2, "The first dimension of targets must be 2"
    nt = targets.shape[1]
    if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
        out.pottarg = hfmm.h2d_directcp(zk,sources,charges,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
        out.pottarg,out.gradtarg = hfmm.h2d_directcg(zk,sources,charges,targets,thresh)
    if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
        out.pottarg = hfmm.h2d_directdp(zk,sources,dipstr,dipvec,targets,thresh)
    if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
        out.pottarg,out.gradtarg = hfmm.h2d_directdg(zk,sources,dipstr,dipvec,targets,thresh)
    if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
        out.pottarg = hfmm.h2d_directcdp(zk,sources,charges,dipstr,dipvec,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
        out.pottarg,out.gradtarg = hfmm.h2d_directcdg(zk,sources,charges,dipstr,dipvec,targets,thresh)

    if(nd == 1):
        if(ifcharge==1):
            charges = charges.reshape(ns,)
        if(ifdipole==1):
            dipvec = dipvec.reshape(2,ns)
            dipstr = dipstr.reshape(ns,)
        if(pgt>0):
            out.pottarg = out.pottarg.reshape(nt,)
        if(pgt==2):
            out.gradtarg = out.gradtarg.reshape(2,nt)


    return out


def r2ddir(*,sources,targets,charges=None,dipstr=None,dipvec=None,
          pgt=0,nd=1,thresh=1e-16):
    r"""
      This subroutine computes the N-body Laplace interactions
      in two dimensions where the interaction kernel is given by $log(r)$ 
      and its gradients. 


      .. math::

          u(x) = \sum_{j=1}^{N} c_{j} * log\(\|x-x_{j}\|\) + d_{j}v_{j} \cdot \\nabla( log(\|x-x_{j}\|) )  \, ,

      where $c_{j}$ are the charge densities, $d_{j}$ are the dipole strengths 
      $v_{j}$ are the dipole orientation vectors, and 
      $x_{j}$ are the source locations.

      When |x-x_{m}|leq thresh, the term corresponding to $x_{m}$ is dropped from the
      sum


      Args:
        sources: float(2,n)   
               source locations (x_{j})
        charges: float(nd,n) or float(n)
                charge densities (c_{j})
        dipstr: float(nd,n) or float(n)
                dipole densities  (d_{j})
        dipvec: float(nd,2,n) or float(2,n)
                dipole orientation vectors (v_{j})
        targets: float(2,nt)
                target locations (x)

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
               potenial, gradient, and hessians at targets evaluated if pgt=3
        
        nd:   integer
               number of densities
        thresh: contribution of source x_i, at location x ignored if |x-x_i|<=thresh

      Returns:
        out.pottarg  - potential at target locations if requested
        out.gradtarg - gradient at target locations if requested
        out.hesstarg - hessian at target locations if requested
              
      Example:
        see rfmmexample.py

    r"""

    out = Output()
    assert sources.shape[0] == 2, "The first dimension of sources must be 2"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    ifcharge = 0
    ifdipole = 0
    if(pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
            charges = charges.reshape(1,ns)
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipvec is not None or dipstr is not None):
        if nd == 1 and ns>1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            assert dipvec.shape[0] == 2 and dipvec.shape[1] == ns, "dipole vectors must be of shape [2,number of sources]"
            dipvec = dipvec.reshape(1,2,ns)
            dipstr = dipstr.reshape(1,ns)
        if nd == 1 and ns==1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            assert dipvec.shape[0] == 2, "dipole vectors must be of shape [2,number of sources]"
            dipvec = dipvec.reshape(1,2,ns)
            dipstr = dipstr.reshape(1,ns)
        if nd>1:
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 2 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,2,ns] where nd is number of densities, and ns is number of sources"
            assert dipstr.shape[0] == nd and dipstr.shape[1]==ns, "Dipole strengths must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1

    assert targets.shape[0] == 2, "The first dimension of targets must be 2"
    nt = targets.shape[1]
    if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
        out.pottarg = lfmm.r2d_directcp(sources,charges,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
        out.pottarg,out.gradtarg = lfmm.r2d_directcg(sources,charges,targets,thresh)
    if(pgt == 3 and ifcharge == 1 and ifdipole == 0):
        out.pottarg,out.gradtarg,out.hesstarg = lfmm.r2d_directch(sources,charges,targets,thresh)
    if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
        out.pottarg = lfmm.r2d_directdp(sources,dipstr,dipvec,targets,thresh)
    if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
        out.pottarg,out.gradtarg = lfmm.r2d_directdg(sources,dipstr,dipvec,targets,thresh)
    if(pgt == 3 and ifcharge == 0 and ifdipole == 1):
        out.pottarg,out.gradtarg,out.hesstarg = lfmm.r2d_directdh(sources,dipstr,dipvec,targets,thresh)
    if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
        out.pottarg = lfmm.r2d_directcdp(sources,charges,dipstr,dipvec,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
        out.pottarg,out.gradtarg = lfmm.r2d_directcdg(sources,charges,dipstr,dipvec,targets,thresh)
    if(pgt == 3 and ifcharge == 1 and ifdipole == 1):
        out.pottarg,out.gradtarg,out.hesstarg = lfmm.r2d_directcdh(sources,charges,dipstr,dipvec,targets,thresh)

    if(nd == 1):
        if(ifcharge == 1):
            charges = charges.reshape(ns,)
        if(ifdipole ==1): 
            dipvec = dipvec.reshape(2,ns)
            dipstr = dipstr.reshape(ns,)
        if(pgt>0):
            out.pottarg = out.pottarg.reshape(nt,)
        if(pgt==2):
            out.gradtarg = out.gradtarg.reshape(2,nt)
        if(pgt==3):
            out.hesstarg = out.hesstarg.reshape(3,nt)

    return out


def l2ddir(*,sources,targets,charges=None,dipstr=None,dipvec=None,
          pgt=0,nd=1,thresh=1e-16):
    r"""
      This subroutine computes the N-body Laplace interactions
      in two dimensions where the interaction kernel is given by $log(r)$ 
      and its gradients. 


      .. math::

          u(x) = \sum_{j=1}^{N} c_{j} * log\(\|x-x_{j}\|\) + d_{j}v_{j} \cdot \\nabla( log(\|x-x_{j}\|) )  \, ,

      where $c_{j}$ are the charge densities, $d_{j}$ are the dipole strengths 
      $v_{j}$ are the dipole orientation vectors, and 
      $x_{j}$ are the source locations.

      When |x-x_{m}|leq thresh, the term corresponding to $x_{m}$ is dropped from the
      sum


      Args:
        sources: float(2,n)   
               source locations (x_{j})
        charges: complex(nd,n) or complex(n)
                charge densities (c_{j})
        dipstr: complex(nd,n) or complex(n)
                dipole densities  (d_{j})
        dipvec: float(nd,2,n) or float(2,n)
                dipole orientation vectors (v_{j})
        targets: float(2,nt)
                target locations (x)

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
               potenial, gradient, and hessians at targets evaluated if pgt=3
        
        nd:   integer
               number of densities
        thresh: contribution of source x_i, at location x ignored if |x-x_i|<=thresh

      Returns:
        out.pottarg  - potential at target locations if requested
        out.gradtarg - gradient at target locations if requested
        out.hesstarg - hessian at target locations if requested
              
      Example:
        see lfmmexample.py

    r"""

    out = Output()
    assert sources.shape[0] == 2, "The first dimension of sources must be 2"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    ifcharge = 0
    ifdipole = 0
    if(pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
            charges = charges.reshape(1,ns)
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipvec is not None or dipstr is not None):
        if nd == 1 and ns>1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            assert dipvec.shape[0] == 2 and dipvec.shape[1] == ns, "dipole vectors must be of shape [2,number of sources]"
            dipvec = dipvec.reshape(1,2,ns)
            dipstr = dipstr.reshape(1,ns)
        if nd == 1 and ns==1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            assert dipvec.shape[0] == 2, "dipole vectors must be of shape [2,number of sources]"
            dipvec = dipvec.reshape(1,2,ns)
            dipstr = dipstr.reshape(1,ns)
        if nd>1:
            assert dipvec.shape[0] == nd and dipvec.shape[1] == 2 and dipvec.shape[2] == ns, "Dipole vectors must be of shape [nd,2,ns] where nd is number of densities, and ns is number of sources"
            assert dipstr.shape[0] == nd and dipstr.shape[1]==ns, "Dipole strengths must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1

    assert targets.shape[0] == 2, "The first dimension of targets must be 2"
    nt = targets.shape[1]
    if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
        out.pottarg = lfmm.l2d_directcp(sources,charges,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
        out.pottarg,out.gradtarg = lfmm.l2d_directcg(sources,charges,targets,thresh)
    if(pgt == 3 and ifcharge == 1 and ifdipole == 0):
        out.pottarg,out.gradtarg,out.hesstarg = lfmm.l2d_directch(sources,charges,targets,thresh)
    if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
        out.pottarg = lfmm.l2d_directdp(sources,dipstr,dipvec,targets,thresh)
    if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
        out.pottarg,out.gradtarg = lfmm.l2d_directdg(sources,dipstr,dipvec,targets,thresh)
    if(pgt == 3 and ifcharge == 0 and ifdipole == 1):
        out.pottarg,out.gradtarg,out.hesstarg = lfmm.l2d_directdh(sources,dipstr,dipvec,targets,thresh)
    if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
        out.pottarg = lfmm.l2d_directcdp(sources,charges,dipstr,dipvec,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
        out.pottarg,out.gradtarg = lfmm.l2d_directcdg(sources,charges,dipstr,dipvec,targets,thresh)
    if(pgt == 3 and ifcharge == 1 and ifdipole == 1):
        out.pottarg,out.gradtarg,out.hesstarg = lfmm.l2d_directcdh(sources,charges,dipstr,dipvec,targets,thresh)

    if(nd == 1):
        if(ifcharge == 1):
            charges = charges.reshape(ns,)
        if(ifdipole ==1): 
            dipvec = dipvec.reshape(2,ns)
            dipstr = dipstr.reshape(ns,)
        if(pgt>0):
            out.pottarg = out.pottarg.reshape(nt,)
        if(pgt==2):
            out.gradtarg = out.gradtarg.reshape(2,nt)
        if(pgt==3):
            out.hesstarg = out.hesstarg.reshape(3,nt)

    return out


def c2ddir(*,sources,targets,charges=None,dipstr=None,
          pgt=0,nd=1,thresh=1e-16):
    r"""
      This subroutine computes the N-body Laplace interactions
      in two dimensions where the interaction kernel is given by $log(r)$ 
      and its gradients. 


      .. math::

          u(x) = \sum_{j=1}^{N} c_{j} * log\(\|x-x_{j}\|\) + d_{j}/(x-x_{j})  \, ,

      where $c_{j}$ are the charge densities, $d_{j}$ are the dipole strengths 
      and $x_{j}$ are the source locations.

      When |x-x_{m}|leq thresh, the term corresponding to $x_{m}$ is dropped from the
      sum


      Args:
        sources: float(2,n)   
               source locations (x_{j})
        charges: complex(nd,n) or complex(n)
                charge densities (c_{j})
        dipstr: complex(nd,n) or complex(n)
                dipole densities  (d_{j})
        targets: float(2,nt)
                target locations (x)

        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
               potenial, gradient, and hessians at targets evaluated if pgt=3
        
        nd:   integer
               number of densities
        thresh: contribution of source x_i, at location x ignored if |x-x_i|<=thresh

      Returns:
        out.pottarg  - potential at target locations if requested
        out.gradtarg - gradient at target locations if requested
        out.hesstarg - hessian at target locations if requested
              
      Example:
        see cfmmexample.py

    r"""

    out = Output()
    assert sources.shape[0] == 2, "The first dimension of sources must be 2"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    ifcharge = 0
    ifdipole = 0
    if(pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
            charges = charges.reshape(1,ns)
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        ifcharge = 1
    if(dipstr is not None):
        if nd == 1 and ns>1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            dipstr = dipstr.reshape(1,ns)
        if nd == 1 and ns==1:
            assert dipstr.shape[0] == ns, "Dipole strengths must be same length as second dimension of sources"
            dipstr = dipstr.reshape(1,ns)
        if nd>1:
            assert dipstr.shape[0] == nd and dipstr.shape[1]==ns, "Dipole strengths must be of shape [nd,ns] where nd is number of densities, and ns is number of sources"
        ifdipole = 1

    assert targets.shape[0] == 2, "The first dimension of targets must be 2"
    nt = targets.shape[1]
    if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
        out.pottarg = lfmm.c2d_directcp(sources,charges,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
        out.pottarg,out.gradtarg = lfmm.c2d_directcg(sources,charges,targets,thresh)
    if(pgt == 3 and ifcharge == 1 and ifdipole == 0):
        out.pottarg,out.gradtarg,out.hesstarg = lfmm.c2d_directch(sources,charges,targets,thresh)
    if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
        out.pottarg = lfmm.c2d_directdp(sources,dipstr,targets,thresh)
    if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
        out.pottarg,out.gradtarg = lfmm.c2d_directdg(sources,dipstr,targets,thresh)
    if(pgt == 3 and ifcharge == 0 and ifdipole == 1):
        out.pottarg,out.gradtarg,out.hesstarg = lfmm.c2d_directdh(sources,dipstr,targets,thresh)
    if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
        out.pottarg = lfmm.c2d_directcdp(sources,charges,dipstr,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
        out.pottarg,out.gradtarg = lfmm.c2d_directcdg(sources,charges,dipstr,targets,thresh)
    if(pgt == 3 and ifcharge == 1 and ifdipole == 1):
        out.pottarg,out.gradtarg,out.hesstarg = lfmm.c2d_directcdh(sources,charges,dipstr,targets,thresh)

    if(nd == 1):
        if(ifcharge == 1):
            charges = charges.reshape(ns,)
        if(ifdipole ==1): 
            dipstr = dipstr.reshape(ns,)
        if(pgt>0):
            out.pottarg = out.pottarg.reshape(nt,)
        if(pgt==2):
            out.gradtarg = out.gradtarg.reshape(nt,)
        if(pgt==3):
            out.hesstarg = out.hesstarg.reshape(nt,)

    return out


def bh2ddir(*,sources,targets,charges=None,dipoles=None,
          pgt=0,nd=1,thresh=1e-16):

    r"""
      This subroutine computes the N-body biharmonic interactions 
      in two dimensions where the interaction kernel is related to the
      biharmonic greens function r^2 log (r) and its derivatives


      .. math:: 

          u(x) = \sum_{j=1}^{N} c_{j} * log\(\|x-x_{j}\|\) + 
          \overline{c}_{j} (x-x_{j})/(\overline{x-x_{j}) + d_{j,1}/(x-x_{j}) - 
          d_{j,2}/(\overline{x-x_{j}}) - 
          \overline{d_{j,1}} (x-x_{j})/(\overline{x-x_{j}})^2\, ,

      where $c_{j}$ are the charge densities, $d_{j,1}$, $d_{j,2}$ are the dipole strengths,
      and $x_{j}$ are the source locations.

      When $x=x_{m}$, the term corresponding to $x_{m}$ is dropped from the
      sum


      Args:
        eps: float
             precision requested

        sources: float(2,n)   
               source locations (x_{j})
        charges: complex(nd,n) or complex(n)
               charge densities (c_{j})
        dipoles: complex(nd,2,n) or complex(2,n)
               dipole densities (d_{j,1}, d_{j,2})
        targets: float(2,nt)
                target locations (x)
        pgt:  integer
               target eval flag
               potential at targets evaluated if pgt = 1
               potenial and gradient at targets evaluated if pgt=2
        
        nd:   integer
               number of densities
        thresh: contribution of source x_i, at location x ignored if |x-x_i|<=thresh

      Returns:
        out.pottarg: potential at target locations if requested
        out.gradtarg: gradient at target locations if requested
      
      Example:
        see bhfmmexample.py
    r"""


    out = Output()
    assert sources.shape[0] == 2, "The first dimension of sources must be 2"
    if(np.size(np.shape(sources))==2):
        ns = sources.shape[1]
    if(np.size(np.shape(sources))==1):
        ns = 1
    ifcharge = 0
    ifdipole = 0
    if(pgt == 0):
        print("Nothing to compute, set either pg or pgt to non-zero")
        return out
    ifcharge = 0
    ifdipole = 0
    iftarg = 0
    if charges is not None:
        if nd == 1:
            assert charges.shape[0] == ns, "Charges must be same length as second dimension of sources"
            charges = charges.reshape(1,ns)
        if nd>1:
            assert charges.shape[0] == nd and charges.shape[1]==ns, "Charges must be of shape [nd,ns] where nd is number of densities, and ns is number of sources" 
        charges = charges.reshape([nd,ns])
        ifcharge = 1
    else:
       charges = np.zeros([nd,ns],dtype='complex') 

    if(dipoles is not None):
        if nd == 1 and ns>1:
            assert dipoles.shape[0] == 2 and dipoles.shape[1] == ns, "Dipole must of shape [2, number of sources]"
        if nd == 1 and ns==1:
            assert dipoles.shape[0] == 2, "Dipole must of shape [2, number of sources]"
        if nd>1:
            assert dipoles.shape[0] == nd and dipoles.shape[1] == 2 and dipoles.shape[2]==ns, "Dipole must of shape [nd,2, number of sources]"
        dipoles = dipoles.reshape([nd,2,ns])
        ifdipole = 1
    else:
        dipoles = np.zeros([nd,2,ns],dtype='complex')

    assert targets.shape[0] == 2, "The first dimension of targets must be 2"
    nt = targets.shape[1]
#
    if(pgt == 1 and ifcharge == 1 and ifdipole == 0):
        out.pottarg = bhfmm.bh2d_directcp(sources,charges,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 0):
        out.pottarg,out.gradtarg = bhfmm.bh2d_directcg(sources,charges,targets,thresh)
    if(pgt == 1 and ifcharge == 0 and ifdipole == 1):
        out.pottarg = bhfmm.bh2d_directdp(sources,dipoles,targets,thresh)
    if(pgt == 2 and ifcharge == 0 and ifdipole == 1):
        out.pottarg,out.gradtarg = bhfmm.bh2d_directdg(sources,dipoles,targets,thresh)
    if(pgt == 1 and ifcharge == 1 and ifdipole == 1):
        out.pottarg = bhfmm.bh2d_directcdp(sources,charges,dipoles,targets,thresh)
    if(pgt == 2 and ifcharge == 1 and ifdipole == 1):
        out.pottarg,out.gradtarg = bhfmm.bh2d_directcdg(sources,charges,dipoles,targets,thresh)

    if(nd == 1):
        if(pgt>0):
            out.pottarg = out.pottarg.reshape(nt,)
        if(pgt==2):
            out.gradtarg = out.gradtarg.reshape(2,nt,)

    return out



def comperr(*,ntest,out,outex,pg=0,pgt=0,nd=1,cauchy=0):
    r = 0
    err = 0
    if(nd == 1):
        if(pg > 0):
            if(cauchy==0):
                r = r+la.norm(outex.pot[0:ntest])**2
                err = err+la.norm(outex.pot[0:ntest]-out.pot[0:ntest])**2
            else:
                r = r+la.norm(outex.pot.real[0:ntest])**2
                err = err+la.norm(outex.pot.real[0:ntest]-out.pot.real[0:ntest])**2
        if(pg >= 2):
            if(cauchy==0):
                g = out.grad[:,0:ntest].reshape(2*ntest,)
                gex = outex.grad[:,0:ntest].reshape(2*ntest,)
            else:
                g = out.grad[0:ntest].reshape(ntest,)
                gex = outex.grad[0:ntest].reshape(ntest,)
            r = r +la.norm(gex)**2
            err = err+la.norm(gex-g)**2
        if( pg >= 3):
            if(cauchy==0):
                h = out.hess[:,0:ntest].reshape(3*ntest,)
                hhex = outex.hess[:,0:ntest].reshape(3*ntest,)
            else:
                h = out.hess[0:ntest].reshape(ntest,)
                hhex = outex.hess[0:ntest].reshape(ntest,)
            r = r + la.norm(hhex)**2
            err = err + la.norm(hhex-h)**2
        if(pgt > 0):
            if(cauchy==0):
                r = r+la.norm(outex.pottarg[0:ntest])**2
                err = err+la.norm(outex.pottarg[0:ntest]-out.pottarg[0:ntest])**2
            else:
                r = r+la.norm(outex.pottarg.real[0:ntest])**2
                err = err+la.norm(outex.pottarg.real[0:ntest]-out.pottarg.real[0:ntest])**2
        if(pgt >= 2):
            if(cauchy==0):
                g = out.gradtarg[:,0:ntest].reshape(2*ntest,)
                gex = outex.gradtarg[:,0:ntest].reshape(2*ntest,)
            else:
                g = out.gradtarg[0:ntest].reshape(ntest,)
                gex = outex.gradtarg[0:ntest].reshape(ntest,)
            r = r +la.norm(gex)**2
            err = err+la.norm(gex-g)**2
        if( pgt >= 3):
            if(cauchy==0):
                h = out.hesstarg[:,0:ntest].reshape(3*ntest,)
                hhex = outex.hesstarg[:,0:ntest].reshape(3*ntest,)
            else:
                h = out.hesstarg[0:ntest].reshape(ntest,)
                hhex = outex.hesstarg[0:ntest].reshape(ntest,)
            r = r + la.norm(hhex)**2
            err = err + la.norm(hhex-h)**2
    if(nd > 1):
        if(pg > 0):
            if(cauchy==0):
                p = out.pot[:,0:ntest].reshape(nd*ntest,)
                pex = outex.pot[:,0:ntest].reshape(nd*ntest,)
            else:
                p = out.pot.real[:,0:ntest].reshape(nd*ntest,)
                pex = outex.pot.real[:,0:ntest].reshape(nd*ntest,)
            r = r+la.norm(pex)**2
            err = err+la.norm(p-pex)**2
        if(pg >= 2):
            if(cauchy==0):
                g = out.grad[:,:,0:ntest].reshape(2*nd*ntest,)
                gex = outex.grad[:,:,0:ntest].reshape(2*nd*ntest,)
            else:
                g = out.grad[:,0:ntest].reshape(nd*ntest,)
                gex = outex.grad[:,0:ntest].reshape(nd*ntest,)
            r = r +la.norm(gex)**2
            err = err+la.norm(gex-g)**2
        if( pg >= 3):
            if(cauchy==0):
                h = out.hess[:,:,0:ntest].reshape(3*nd*ntest,)
                hhex = outex.hess[:,:,0:ntest].reshape(3*nd*ntest,)
            else:
                h = out.hess[:,0:ntest].reshape(nd*ntest,)
                hhex = outex.hess[:,0:ntest].reshape(nd*ntest,)
            r = r + la.norm(hhex)**2
            err = err + la.norm(hhex-h)**2
        if(pgt > 0):
            if(cauchy==0):
                p = out.pottarg[:,0:ntest].reshape(nd*ntest,)
                pex = outex.pottarg[:,0:ntest].reshape(nd*ntest,)
            else:
                p = out.pottarg.real[:,0:ntest].reshape(nd*ntest,)
                pex = outex.pottarg.real[:,0:ntest].reshape(nd*ntest,)
            r = r+la.norm(pex)**2
            err = err+la.norm(p-pex)**2
        if(pgt >= 2):
            if(cauchy==0):
                g = out.gradtarg[:,:,0:ntest].reshape(2*nd*ntest,)
                gex = outex.gradtarg[:,:,0:ntest].reshape(2*nd*ntest,)
            else:
                g = out.gradtarg[:,0:ntest].reshape(nd*ntest,)
                gex = outex.gradtarg[:,0:ntest].reshape(nd*ntest,)
            r = r +la.norm(gex)**2
            err = err+la.norm(gex-g)**2
        if( pgt >= 3):
            if(cauchy==0):
                h = out.hesstarg[:,:,0:ntest].reshape(3*nd*ntest,)
                hhex = outex.hesstarg[:,:,0:ntest].reshape(3*nd*ntest,)
            else:
                h = out.hesstarg[:,0:ntest].reshape(nd*ntest,)
                hhex = outex.hesstarg[:,0:ntest].reshape(nd*ntest,)
            r = r + la.norm(hhex)**2
            err = err + la.norm(hhex-h)**2
    err = np.sqrt(err/r)
    return err
