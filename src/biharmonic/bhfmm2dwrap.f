

      subroutine bhfmm2dwrap_guru(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dip,iper,ifpgh,pot,grad,
     2            hess,nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
cf2py intent(in) nd,eps,ns,sources,ifcharge,charge,ifdipole,dip
cf2py intent(in) iper,ifpgh,nt,targ,ifpghtarg
cf2py intent(out) pot,grad,hess,pottarg,gradtarg,hesstarg,ier
c
c----------------------------------------------
c
c   This subroutine evaluates biharmonic sums related to stokes
c   and elasticity due to point ``charges'' and ``dipoles''
c
c
c   Notes:
c     * charges are complex, and dipoles are described by vectors
c       in C^2
c     * In the legacy interface, dipoles were input as two
c       complex arrays dip1 and dip2, as opposed to a single
c       complex array (nd,2,n) where n is the number of sources.
c       Similarly, the gradients in the legacy interface
c       were two complex arrays, grada, and gradaa; while
c       in the new interface, the gradient is a (nd,2,n) 
c       complex array where (nd,1,n) is the derivative with
c       respect to z and the (nd,2,n) component is the derivative
c       with respect to z_bar
c     * In this subroutine and the rest the terms potential
c       and velocity are interchangably used
c     
c
c    vel(z) = \sum charge1_j 2*\log|z-z_j| + 
c                (charge2_j) (z-z_j)/(z-z_j)_bar+
c                dip1/(z-z_j) + dip3/(z-z_j)_bar+
c                dip2 (z-z_j)/(z-z_j)^2_bar
c
c    The velocity for stokes is related to the goursat functions
c    as vel = \phi(z) + z (d/dz (\phi))_bar+\psi(z)
c 
c   The goursat functions are given by
c
c   \phi(z) = \sum charge_j log(z-z_j)+
c                     dip1/(z-z_j)
c
c   \psi(z) = \sum dippar2_bar/(z-z_j)+
c               z_j_bar dippar1/(z-z_j)^2 +
c
c    In all of the sums above, the terms for which |z-z_{j}| <= 2^-51|b|
c    where |b| is the size of the root box will be ignored.
c
c   Gradient here has
c    
c
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   ifcharge      : flag for including charge interactions
c                   charge interactions included if ifcharge =1
c                   not included otherwise
c   charge(nd,2,ns) : charge strengths
c                   charge(j,1,i) is the charge1 strength of the
c                   jth charge density at source location i, and
c                   charge(j,2,i) is the charge2 stregth of the 
c                   jth charge density at source location i
c   ifdipole      : flag for including dipole interactions
c                   dipole interactions included if ifcharge =1
c                   not included otherwise
c   dip(nd,3,ns)  : dipole strengths
c                   dip(j,1,i) is dip1 for the jth dipole density
c                   at source location i, dip(j,2,i) is dip2
c                   for the jth dipole density at source location i, and
c                   dip(j,3,i) is dip3 for the jth dipole density
c                   at source locaation i
c   iper          : flag for periodic implmentations. Currently unused
c   ifpgh         : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed (hessians not supported currently)
c   nt            : number of targets
c   targ(2,nt)    : target locations
c   ifpghtarg     : flag for computing pottarg/gradtarg/hesstarg
c                   ifpghtarg = 1, only potential is computed at targets
c                   ifpghtarg = 2, potential and gradient are 
c                   computed at targets
c                   ifpghtarg = 3, potential, gradient, and hessian are 
c                   computed at targets (hessians not supported
c                     currently)
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : velocity at the source locations
c   grad(nd,3,*)    : gradients (d/dz)_projz, (d/dz)_projzbar, d/dzbar at the source locations
c   hess(nd,3,*)    : hessian  at the source locations (currently unused)
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,3,*): gradient (d/dz)_projz, (d/dz)_projzbar, d/dzbar at the target locations
c   hesstarg(nd,3,*): hessian at the target locations (currently unused)
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 omp_get_wtime
      real *8 eps
      integer ns,nt
      integer ifcharge,ifdipole
      integer iper,ier,ifpgh,ifpghtarg
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,2,ns),dip(nd,3,ns)

      complex *16 pot(nd,ns),grad(nd,3,ns), hess(nd,3,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,3,nt)
      complex *16 hesstarg(nd,3,nt)

      call bhfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dip,iper,ifpgh,pot,grad,
     2            hess,nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
