c
c     This file contains the Stokes FMM wrappers
c
c**********************************************************************
c      
c     We take the following conventions for the Stokes kernels
c
c     For a source y and target x, let r_i = x_i-y_i
c     and let r = sqrt(r_1^2 + r_2^2)
c
c     The Stokeslet, G_{ij}, and its associated pressure tensor, P_j,
c     (without the 1/2pi scaling) are
c
c     G_{ij}(x,y) = (r_i r_j)/(2r^2) - delta_{ij}log(r)/(2)
c     P_j(x,y) = r_j/r^2
c
c     The (Type I) stresslet, T_{ijk}, and its associated pressure
c     tensor, PI_{jk}, (without the 1/2pi scaling) are
c     
c     T_{ijk}(x,y) = -2 r_i r_j r_k/ r^4
c     PI_{jk} = - delta_{jk}/r^2 + 2 r_j r_k/r^4      
c
c
c
      

      subroutine stfmm2d(nd, eps, 
     $                 nsource, source,
     $                 ifstoklet, stoklet, ifstrslet, strslet, strsvec,
     $                 ifppreg, pot, pre, grad, ntarg, targ, 
     $                 ifppregtarg, pottarg, pretarg, gradtarg,ier)

cf2py  intent(in) nd,eps
cf2py  intent(in) nsource,source
cf2py  intent(in) ifstoklet,stoklet
cf2py  intent(in) ifstrslet,strslet,strsvec
cf2py  intent(in) ifppreg,ifppregtarg
cf2py  intent(in) ntarg,targ
cf2py  intent(out) pot,pre,grad
cf2py  intent(out) pottarg,pretarg,gradtarg
cf2py  intent(out) ier
c
c     Stokes FMM in R^{3}: evaluate all pairwise particle
c     interactions (ignoring self-interactions) and
c     interactions with targs.
c      
c     This routine computes sums of the form
c
c       u(x) = sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k
c
c     where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
c     stresslet charge, and nu^{(m)} is the stresslet orientation
c     (note that each of these is a 3 vector per source point y^{(m)}).
c     For x a source point, the self-interaction in the sum is omitted. 
c
c     Optionally, the associated pressure p(x) and gradient grad u(x)
c     are returned
c
c       p(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
c          + sum_m T_{ijk}(x,y^{(m)}) PI_{jk} mu^{(m)}_j nu^{(m)}_k
c
c       grad u(x) = grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
c                + sum_m T_{ijk}(x,y^{(m)}) mu^{(m)}_j nu^{(m)}_k]
c
c-----------------------------------------------------------------------
c     INPUT PARAMETERS:
c     
c   nd:    in: integer
c              number of densities
c   
c   eps:   in: double precision
c              requested precision
c
c   nsource in: integer  
c               number of sources
c
c   source  in: double precision (2,nsource)
c               source(k,j) is the kth component of the jth
c               source locations
c
c   ifstoklet  in: integer  
c               Stokeslet charge computation flag
c               ifstoklet = 1   =>  include Stokeslet contribution
c                                   otherwise do not
c 
c   stoklet in: double precision (nd,2,nsource) 
c               Stokeslet charge strengths (sigma vectors above)
c
c   ifstrslet in: integer
c               stresslet computation flag
c               ifstrslet = 1   =>  include standard stresslet
c                                   (type I)
c
c            NOT YET IMPLEMENTED
c      
c               ifstrslet = 2   =>  include symmetric stresslet
c                                   (type II)
c               ifstrslet = 3   =>  include rotlet
c               ifstrslet = 4   =>  include Stokes doublet
c                      otherwise do not include
c
c   strslet  in: double precision (nd,2,nsource) 
c               stresslet strengths (mu vectors above)
c
c   strsvec  in: double precision (nd,2,nsource)   
c               stresslet orientations (nu vectors above)
c
c     ifppreg    in: integer      
c               flag for evaluating potential, gradient, and pressure
c               at the sources
c               ifppreg = 1, only potential
c               ifppreg = 2, potential and pressure
c               ifppreg = 3, potential, pressure, and gradient 
c      
c   ntarg   in: integer  
c              number of targs 
c
c   targ    in: double precision (2,ntarg)
c             targ(k,j) is the kth component of the jth
c             targ location
c      
c   ifppregtarg in: integer
c                flag for evaluating potential, gradient, and pressure
c                at the targets
c                ifppregtarg = 1, only potential
c                ifppregtarg = 2, potential and pressure
c                ifppregtarg = 3, potential, pressure, and gradient
c
c-----------------------------------------------------------------------
c
c   OUTPUT parameters:
c
c   pot   out: double precision(nd,2,nsource) 
c           velocity at the source locations
c      
c   pre   out: double precision(nd,nsource)
c           pressure at the source locations
c      
c         GRADIENT NOT IMPLEMENTED
c   grad   out: double precision(nd,2,2,nsource) 
c              gradient of velocity at the source locations
c              grad(l,i,j,k) is the ith component of the
c              gradient of the jth component of the velocity
c              for the lth density at the kth source location
c     
c   pottarg   out: double precision(nd,2,ntarg) 
c               velocity at the targets
c      
c   pretarg   out: double precision(nd,ntarg)
c               pressure at the targets
c      
c   gradtarg   out: double precision(nd,2,2,ntarg) 
c               gradient of velocity at the targets
c               gradtarg(l,i,j,k) is the ith component of the
c               gradient of the jth component of the velocity
c               for the lth density at the kth target
c     ier     out: integer
c               error flag
c
c     TODO: implement other stresslet options and gradient
c------------------------------------------------------------------
      implicit none
      integer nd, ifstoklet, ifstrslet, ntarg
      double precision eps
      integer nsource, ifppreg, ifppregtarg
      double precision source(2, nsource), targ(2, ntarg)
      double precision stoklet(nd, 2, nsource), strslet(nd, 2, nsource)
      double precision strsvec(nd, 2, nsource)
      double precision pot(nd, 2, nsource), pre(nd, nsource)
      double precision grad(nd, 2, 2, nsource)
      double precision pottarg(nd, 2, ntarg), pretarg(nd,ntarg),
     1     gradtarg(nd, 2, 2, ntarg)     
c
c  CONTINUE FROM HERE
c
c


c     local
      complex *16, allocatable :: charge(:,:,:),dip(:,:,:)
      complex *16, allocatable :: potl(:,:),gradl(:,:,:),zsum(:)
      complex *16, allocatable :: pottargl(:,:),gradtargl(:,:,:)
      
      complex *16 hesstmp(10),ima,zn,zd1,zd2

      integer ndl, ifchargel, ifdipolel, ifpghl, ifpghtargl

      integer i, j, ii, ifppreg1, l, npt, ier,iper
      data ima/(0.0d0,1.0d0)/
      

      ifdipolel = 0
      ifchargel = 0

      if (ifstoklet .eq. 1) ifchargel = 1
      if (ifstrslet .eq. 1) ifdipolel = 1
      


c     allocate necessary arrays
      
      allocate(charge(nd,2,nsource),dip(nd,3,nsource),
     1     potl(nd,nsource),pottargl(nd,ntarg),
     2     gradl(nd,3,nsource),gradtargl(nd,3,ntarg),stat=ier)
      allocate(zsum(nd))
      if(ier .ne. 0) then
         print *, "In stfmm2d: cannot allocate Laplace call storage"
         print *, "nd =",nd
         print *, "nsource =",nsource
         print *, "ntarg =",ntarg       
         ier = 4
         return
      endif

c     set-up appropriate vector charge and dipole arrays
      do j=1,nd
        zsum(j) = 0
      enddo
ccc$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,zd1,zd2) REDUCTION(+:zsum)
      do i = 1,nsource
         do j = 1,nd
           charge(j,1,i) = 0
           charge(j,2,i) = 0
           dip(j,1,i) = 0
           dip(j,2,i) = 0
           dip(j,3,i) = 0
         enddo

         if(ifstoklet.eq.1) then
           do j=1,nd
             charge(j,1,i) = (-ima*stoklet(j,1,i) + stoklet(j,2,i))/4
             charge(j,2,i) = dconjg(charge(j,1,i))
             zsum(j) = zsum(j) - charge(j,1,i)
           enddo
         endif

         if(ifstrslet.eq.1) then
           do j=1,nd
             zd1 = -(-strslet(j,2,i) + ima*strslet(j,1,i))/2
             zd2 =  (-strsvec(j,2,i) + ima*strsvec(j,1,i))
             dip(j,1,i) = -ima*zd1*zd2 
             dip(j,2,i) = -dconjg(dip(j,1,i))
             dip(j,3,i) =  ima*(zd1*dconjg(zd2) + dconjg(zd1)*zd2) 
           enddo
         endif
      enddo
ccc$OMP END PARALLEL DO      
      ifpghl = 0
      ifpghtargl = 0
      if(ifppreg.eq.1) ifpghl = 1
      if(ifppreg.ge.2) ifpghl = 2

      if(ifppregtarg.eq.1) ifpghtargl = 1
      if(ifppregtarg.ge.2) ifpghtargl = 2




c     call biharmonic FMM

      iper = 0
      ier = 0
      
      call bhfmm2d(nd,eps,nsource,source,ifchargel,charge,
     1  ifdipolel,dip,iper,ifpghl,potl,gradl,hesstmp,ntarg,
     2  targ,ifpghtargl,pottargl,gradtargl,hesstmp,ier)

      if(ifppreg.ge.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)      
        do i=1,nsource
          do j=1,nd
            pot(j,1,i) =  imag(potl(j,i)+zsum(j)+charge(j,1,i)) 
            pot(j,2,i) = -real(potl(j,i)+zsum(j)+charge(j,1,i)) 
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifppreg.ge.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
        do i=1,nsource
          do j=1,nd
            pre(j,i) = -4*imag(gradl(j,1,i))
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifppreg.ge.3) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
        do i=1,nsource
          do j=1,nd
             grad(j,1,1,i) =  imag(gradl(j,3,i))
             grad(j,2,2,i) = -imag(gradl(j,3,i))
             grad(j,2,1,i) =  real(2*gradl(j,1,i)-gradl(j,3,i))
             grad(j,1,2,i) = -real(2*gradl(j,1,i)+gradl(j,3,i))
          enddo
        enddo
C$OMP END PARALLEL DO
      endif




      if(ifppregtarg.ge.1) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)      
        do i=1,ntarg
          do j=1,nd
            pottarg(j,1,i) =  imag(pottargl(j,i)+zsum(j)) 
            pottarg(j,2,i) = -real(pottargl(j,i)+zsum(j)) 
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifppregtarg.ge.2) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
        do i=1,ntarg
          do j=1,nd
            pretarg(j,i) = -4*imag(gradtargl(j,1,i))
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

      if(ifppregtarg.ge.3) then
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
        do i=1,ntarg
          do j=1,nd
             gradtarg(j,1,1,i) =  imag(gradtargl(j,3,i))
             gradtarg(j,2,2,i) = -imag(gradtargl(j,3,i))
             gradtarg(j,2,1,i) =  real(2*gradtargl(j,1,i)-
     1           gradtargl(j,3,i))
             gradtarg(j,1,2,i) = -real(2*gradtargl(j,1,i)+
     1           gradtargl(j,3,i))
          enddo
        enddo
C$OMP END PARALLEL DO
      endif

      
      return
      end

