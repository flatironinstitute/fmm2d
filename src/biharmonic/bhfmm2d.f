cc Copyright (C) 2018-2019: Leslie Greengard, Zydrunas Gimbutas, 
cc and Manas Rachh
c     c Contact: greengard@cims.nyu.edu
cc
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$

      subroutine bhfmm2d(nd,eps,ns,sources,ifcharge,charge,
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
c       complex array (nd,3,n) where n is the number of sources.
c       Similarly, the gradients in the legacy interface
c       were two complex arrays, grada, and gradaa; while
c       in the new interface, the gradient is a (nd,3,n) 
c       complex array where (nd,1,n) is the derivative with
c       respect to z projection onto z component, 
c       the (nd,2,n) component is the derivative with respect to z
c       projection onto zbar component
c       and the (nd,3,n) component is the derivative
c       with respect to z_bar
c     * In this subroutine and the rest the terms potential
c       and velocity are interchangably used
c     
c
c    vel(z) = \sum 2*charge1_j *\log|z-z_j| + 
c                (charge2_j) (z-z_j)/(z-z_j)_bar+
c                dip1_j/(z-z_j) + dip2_j/(z-z_j)_bar-
c                dip3_j (z-z_j)/(z-z_j)^2_bar
c
c    For the stokes case, if charge2_j = charge1_j_bar, and dip3 =
c    dip1_bar, then the velocity is related to the goursat functions
c    as vel = \phi(z) + z (d/dz (\phi))_bar+\psi(z), where
c 
c    The goursat functions are given by
c
c    \phi(z) = \sum charge_j log(z-z_j)+
c                     dip1/(z-z_j)
c
c    \psi(z) = \sum dip2_bar/(z-z_j)+
c               z_j_bar dip1/(z-z_j)^2 +
c
c    In all of the sums above, the terms for which |z-z_{j}| <= 2^-51|b|
c    where |b| is the size of the root box will be ignored.
c
c    grad(:,1,:) = d/dz vel(z)_proj{z} = 
c                     charge1_j/(z-z_j) - dip1_j/(z-z_j)^2
c    grad(:,2,:) = d/dz vel(z)_proj{zbar} = 
c                     charge2_j/(z-z_j)_bar - dip3_j/(z-z_j)^2_bar
c    grad(:,3,:) = d/dzbar vel(z) = 
c                    -charge2_j (z-z_j)/(z-z_j)^2_bar - dip2_j
c                             (z-z_j)^2_bar +
c                     2dip3_j(z-z_j)/(z-z_j)^3_bar
c
c   INPUT PARAMETERS:
c   nd              : number of expansions
c   eps             : FMM precision requested
c   ns              : number of sources
c   sources(2,ns)   : source locations
c   ifcharge        : flag for including charge interactions
c                     charge interactions included if ifcharge =1
c                     not included otherwise
c   charge(nd,2,ns) : charge strengths
c                     charge(j,1,i) is the charge1 strength of the
c                     jth charge density at source location i,
c                     and charge2(j,1,i) is the charge2 strength
c   ifdipole        : flag for including dipole interactions
c                     dipole interactions included if ifcharge =1
c                     not included otherwise
c   dip(nd,3,ns)    : dipole strengths
c                     dip(j,1,i) is dip1 for the jth dipole density
c                     at source location i, dip(j,2,i) is dip2
c                     for the jth dipole density at source location i,
c                     and dip(j,3,i) is dip3
c   iper            : flag for periodic implmentations. Currently unused
c   ifpgh           : flag for computing pot/grad/hess
c                     ifpgh = 1, only potential is computed
c                     ifpgh = 2, potential and gradient are computed
c                     ifpgh = 3, potential, gradient, and hessian 
c                     are computed (hessians not supported currently)
c   nt              : number of targets
c   targ(2,nt)      : target locations
c   ifpghtarg       : flag for computing pottarg/gradtarg/hesstarg
c                     ifpghtarg = 1, only potential is computed at targets
c                     ifpghtarg = 2, potential and gradient are 
c                     computed at targets
c                     ifpghtarg = 3, potential, gradient, and hessian are 
c                     computed at targets (hessians not supported
c                     currently)
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : velocity at the source locations
c   grad(nd,3,*)    : gradients (d/dz)_proj{z}, 
c                     (d/dz)_proj{zbar}, d/dzbar at the source locations
c   hess(nd,3,*)    : hessian  at the source locations (currently unused)
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,3,*): gradient (d/dz)_proj{z}, 
c                     (d/dz)_proj{zbar}, d/dzbar at the target locations
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
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,2,ns),dip(nd,3,ns)

      complex *16 pot(nd,ns),grad(nd,3,ns), hess(nd,3,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,3,nt)
      complex *16 hesstarg(nd,3,nt)

c
cc      Tree variables
c
      integer, allocatable :: itree(:)
      integer iptr(8)
      integer iper,nlmin,nlmax,ifunif
      real *8, allocatable :: tcenters(:,:),boxsize(:)
      integer nexpc,ntj
      real *8 expc(2)
      real *8 scj
      complex *16 jexps(100)
      integer idivflag,nlevels,nboxes,ndiv
      integer ltree

c
cc     sorted arrays
c
      integer, allocatable :: isrc(:),isrcse(:,:)
      integer, allocatable :: itarg(:),itargse(:,:),iexpcse(:,:)
      real *8, allocatable :: sourcesort(:,:)
      real *8, allocatable :: targsort(:,:)
      complex *16, allocatable :: chargesort(:,:,:),dipsort(:,:,:)
      complex *16, allocatable :: potsort(:,:),gradsort(:,:,:),
     1                             hesssort(:,:,:)
      complex *16, allocatable :: pottargsort(:,:),gradtargsort(:,:,:),
     1                            hesstargsort(:,:,:)

c
cc     additional fmm variables

      integer lmptot
      real *8, allocatable :: rscales(:)
      integer, allocatable :: nterms(:),iaddr(:,:)
      real *8, allocatable :: rmlexp(:)
      complex *16, allocatable :: mptemp(:)

c
cc      temporary variables
c
      integer i,ilev,lmptmp,nmax,idim
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg,ifprint,ier
      real *8 time1,time2,pi,done

      done = 1
      pi = atan(done)*4.0d0


      nexpc = 0

c
c    Need to fix ndiv in biharmonic FMM
c   
c

      nlevels = 0
      nboxes = 0

      call bhndiv2d(eps,ns,nt,ifcharge,ifdipole,ifpgh,
     1  ifpghtarg,ndiv,idivflag)

      ltree = 0
      nlmin = 0
      nlmax = 51
      ifunif = 0
      iper = 0

      ifprint = 0
c
cc      call the tree memory management
c       code to determine number of boxes,
c       number of levels and length of tree
c

      
      call pts_tree_mem(sources,ns,targ,nt,idivflag,ndiv,nlmin,nlmax,
     1  ifunif,iper,nlevels,nboxes,ltree)


      allocate(itree(ltree))
      allocate(boxsize(0:nlevels))
      allocate(tcenters(2,nboxes))

c
c       call the tree code
c

      call pts_tree_build(sources,ns,targ,nt,idivflag,ndiv,nlmin,nlmax,
     1  ifunif,iper,nlevels,nboxes,ltree,itree,iptr,tcenters,boxsize)

      allocate(isrc(ns),isrcse(2,nboxes))
      allocate(itarg(nt),itargse(2,nboxes),iexpcse(2,nboxes))
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nboxes
        iexpcse(1,i) = 1
        iexpcse(2,i) = 0
      enddo
C$OMP END PARALLEL DO      

      call pts_tree_sort(ns,sources,itree,ltree,nboxes,nlevels,iptr,
     1   tcenters,isrc,isrcse)

      call pts_tree_sort(nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1     tcenters,itarg,itargse)

      
      allocate(sourcesort(2,ns))
      allocate(targsort(2,nt))


      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        allocate(chargesort(nd,2,ns),dipsort(nd,3,1))
      endif
      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        allocate(chargesort(nd,2,1),dipsort(nd,3,ns))
      endif
      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        allocate(chargesort(nd,2,ns),dipsort(nd,3,ns))
      endif

      if(ifpgh.eq.1) then
        allocate(potsort(nd,ns),gradsort(nd,3,1),hesssort(nd,3,1))
      else if(ifpgh.eq.2) then
        allocate(potsort(nd,ns),gradsort(nd,3,ns),hesssort(nd,3,1))
      else if(ifpgh.eq.3) then
        allocate(potsort(nd,ns),gradsort(nd,3,ns),hesssort(nd,3,ns))
      else
        allocate(potsort(nd,1),gradsort(nd,3,1),hesssort(nd,3,1))
      endif

      
      if(ifpghtarg.eq.1) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,3,1),
     1     hesstargsort(nd,3,1))
      else if(ifpghtarg.eq.2) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,3,nt),
     1      hesstargsort(nd,3,1))
      else if(ifpghtarg.eq.3) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,3,nt),
     1     hesstargsort(nd,3,nt))
      else
        allocate(pottargsort(nd,1),gradtargsort(nd,3,1),
     1     hesstargsort(nd,3,1))
      endif
      
c
cc      initialize potentials,hessians,gradients
c


      if(ifpgh.eq.1) then
        do i=1,ns
          do idim=1,nd
            potsort(idim,i) = 0
          enddo
        enddo
      endif

      if(ifpgh.eq.2) then
        do i=1,ns
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
            gradsort(idim,3,i) = 0
          enddo
        enddo
      endif

      if(ifpgh.eq.3) then
        do i=1,ns
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
            gradsort(idim,3,i) = 0
            hesssort(idim,1,i) = 0
            hesssort(idim,2,i) = 0
            hesssort(idim,3,i) = 0
          enddo
        enddo
      endif


      if(ifpghtarg.eq.1) then
        do i=1,nt
          do idim=1,nd
            pottargsort(idim,i) = 0
          enddo
        enddo
      endif

      if(ifpghtarg.eq.2) then
        do i=1,nt
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
            gradtargsort(idim,3,i) = 0
          enddo
        enddo
      endif

      if(ifpghtarg.eq.3) then
        do i=1,nt
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
            gradtargsort(idim,3,i) = 0
            hesstargsort(idim,1,i) = 0
            hesstargsort(idim,2,i) = 0
            hesstargsort(idim,3,i) = 0
          enddo
        enddo
      endif


c
cc      compute scaling factor for multipole/local expansions
c       and lengths of multipole and local expansions
c
      allocate(rscales(0:nlevels),nterms(0:nlevels))

      nmax = 0
      do i=0,nlevels
        rscales(i) = boxsize(i)
        call bh2dterms(eps,nterms(i),ier)
        nterms(i) = nterms(i) 
        if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo

      if(ifprint.eq.1) call prinf('nmax=*',nmax,1)
      if(ifprint.eq.1) call prinf('nterms=*',nterms,nlevels+1)

c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c

      allocate(iaddr(2,nboxes))

      lmptmp = (nmax+1)*nd*5
      allocate(mptemp(lmptmp))

c     reorder sources
c
      call dreorderf(2,ns,sources,sourcesort,isrc)
      if(ifcharge.eq.1) 
     1     call dreorderf(4*nd,ns,charge,chargesort,isrc)
      if(ifdipole.eq.1) then
         call dreorderf(6*nd,ns,dip,dipsort,isrc)
      endif

c
cc     reorder targets
c
      call dreorderf(2,nt,targ,targsort,itarg)



c
c
c     allocate memory need by multipole, local expansions at all
c     levels
c     irmlexp is pointer for workspace need by various fmm routines,
c
      call bh2dmpalloc(nd,itree,iaddr,nlevels,lmptot,
     1     nterms)
      if(ifprint .eq. 1) call prinf(' lmptot is *',lmptot,1)

      allocate(rmlexp(lmptot),stat=ier)


c
cc     call the main fmm routine
c

c     Memory allocation is complete. 
c     Call main fmm routine
c
      call cpu_time(time1)
C$      time1=omp_get_wtime()
      call bhfmm2dmain(nd,eps,
     $   ns,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,dipsort,
     $   nt,targsort,nexpc,expc,
     $   iaddr,rmlexp,mptemp,lmptmp,
     $   itree,ltree,iptr,ndiv,nlevels,
     $   nboxes,iper,boxsize,rscales,tcenters,itree(iptr(1)),
     $   isrcse,itargse,iexpcse,nterms,ntj,
     $   ifpgh,potsort,gradsort,hesssort,
     $   ifpghtarg,pottargsort,gradtargsort,
     $   hesstargsort,jexps,scj)
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)


c
cc      resort the output arrays in input order
c

      if(ifpgh.eq.1) then
        call dreorderi(2*nd,ns,potsort,pot,isrc)
      endif

      if(ifpgh.eq.2) then
        call dreorderi(2*nd,ns,potsort,pot,isrc)
        call dreorderi(6*nd,ns,gradsort,grad,isrc)
      endif

      if(ifpgh.eq.3) then
        call dreorderi(2*nd,ns,potsort,pot,isrc)
        call dreorderi(6*nd,ns,gradsort,grad,isrc)
        call dreorderi(6*nd,ns,hesssort,hess,isrc)
      endif

      if(ifpghtarg.eq.1) then
        call dreorderi(2*nd,nt,pottargsort,pottarg,itarg)
      endif

      if(ifpghtarg.eq.2) then
        call dreorderi(2*nd,nt,pottargsort,pottarg,itarg)
        call dreorderi(6*nd,nt,gradtargsort,gradtarg,itarg)
      endif

      if(ifpghtarg.eq.3) then
        call dreorderi(2*nd,nt,pottargsort,pottarg,itarg)
        call dreorderi(6*nd,nt,gradtargsort,gradtarg,itarg)
        call dreorderi(6*nd,nt,hesstargsort,hesstarg,itarg)
      endif


      return
      end
c
c
c
c
c
      subroutine bhfmm2dmain(nd,eps,
     $     nsource,sourcesort,
     $     ifcharge,chargesort,
     $     ifdipole,dipsort,
     $     ntarget,targetsort,nexpc,expcsort,
     $     iaddr,rmlexp,mptemp,lmptmp,
     $     itree,ltree,iptr,ndiv,nlevels, 
     $     nboxes,iper,boxsize,rscales,centers,laddr,
     $     isrcse,itargse,iexpcse,nterms,ntj,
     $     ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg,
     $     jsort,scjsort)
c
c
c              
c   Generalized biharmonic FMM in R^2: evaluate all pairwise particle
c   interactions (ignoring self-interaction). 
c   
c   bhfmm2d: chargesort and dipsort are complex valued,
c   z are complex numbers.
c
c   See documentation above for definitions of velocity 
c 
c
c   All the source/target/expansion center related quantities
c   are assumed to be tree-sorted
c
c-----------------------------------------------------------------------
c  INPUT PARAMETERS:
c    - nd: integer
c        number of charge densities
c    - eps:  real *8
c        FMM precision requested
c    - nsource: integer
c        number of sources
c    - sourcesort: real *8 (2,nsource)
c        tree sorted source locations
c    - ifcharge: integer
c        charge computation flag
c        ifcharge = 1 => include charge contribution
c                           otherwise do not
c    - chargesort: complex *16 (nd,2,nsource): 
c        charge strengths
c    - ifdipole:  integer
c        dipole computation flag
c        ifdipole = 1 => include dipole contribution
c                            otherwise do not
c    - dipsort: complex *16 (nd,3,nsource): 
c        dipole strengths
c    - ntarget: integer
c        number of targets
c    - targetsort: real *8 (2,ntarget)
c        target locations
c    - nexpc: integer
c        number of expansion centers
c    - expcsort: real *8 (2,nexpc)
c        expansion center locations
c    - iaddr: integer (2,nboxes)
c        pointer in rmlexp where multipole and local expansions for 
c        each box is stored
c        * iaddr(1,ibox) is the starting index in rmlexp for the 
c          multipole expansion of ibox, and
c        * iaddr(2,ibox) is the starting index in rmlexp for the
c          local expansion of ibox
c    - mptemp: complex *16(lmptmp): 
c        temporary multipole/local expansion
c        (may not be needed in new setting)
c    - lmptmp: integer
c        length of temporary expansion
c    - itree: integer (ltree)
c        This array contains all the information about the tree
c        See to pts_tree2d.f for documentation
c    - ltree: integer
c        length of itree
c    - iptr: integer(8)
c         iptr is a collection of pointers 
c         which points to where different elements 
c         of the tree are stored in the itree array
c    - ndiv: integer
c        Max number of points per box
c    - nlevels: integer
c        number of levels in the tree
c    - nboxes: integer
c        number of boxes in the tree
c    - boxsize: real*8 (0:nlevels)
c        boxsize(i) is the size of the box from end to end
c        at level i
c    - iper: integer
c        flag for periodic implementation
c    - centers: real *8(2,nboxes)
c        array containing the centers of all the boxes
c    - isrcse: integer(2,nboxes)
c        starting and ending location of sources in ibox
c        in sorted list of sources
c    - itargse: integer(2,nboxes)
c        starting and ending location of targets in ibox
c        in sorted list of sources
c    - iexpcse: integer(2,nboxes)
c        starting and ending location of expansion centers
c        in ibox in sorted list of sources
c    - nterms: integer(0:nlevels) 
c        length of multipole and local expansions at various levels
c    - ntj: integer
c        order of the output expansions (unused)
c    - ifpgh: integer
c        flag for evaluating potential/gradients/hessians at sources.
c        * ifpgh = 1, only potentials will be evaluated
c        * ifpgh = 2, potentials/gradients will be evaluated
c        * ifpgh = 3, potentials/gradients/hessians will be evaluated
c    - ifpghtarg: integer
c        flag for evaluating potential/gradients/hessians at targets.
c        * ifpghtarg = 1, only potentials will be evaluated
c        * ifpghtarg = 2, potentials/gradients will be evaluated
c        * ifpghtarg = 3, potentials/gradients/hessians will be evaluated
c
c  OUTPUT PARAMETERS:
c    - pot: complex *16 (nd,nsource)
c        potential at the source locations
c    - grad: complex *16 (nd,3,nsource)
c        gradient at the source locations (d/dz_proj(z),
c        d/dz_proj(zbar),d/dzbar)
c    - hess: complex *16 (nd,3,nsource)
c        currently unused
c    - pottarg: complex *16 (nd,ntarget)
c        potential at the target locations
c    - gradtarg: complex *16 (nd,3,ntarget)
c        gradient at the target locations (d/dz_proj(z),
c        d/dz_proj(zbar),d/dzbar)
c    - hesstarg: complex *16 (nd,3,ntarget)
c        currently unused
c    - jexps: complex *16 (nd,5,0:ntj,nexpc)
c        currently unused
c    - scjsort: real *8 (nexpc)
c        currently unused
c------------------------------------------------------------------

      implicit none

      integer nd

      integer iper

      integer nsource,ntarget,nexpc
      integer ndiv,nlevels,ntj

      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 eps

      real *8 sourcesort(2,nsource)

      complex *16 chargesort(nd,2,*)
      complex *16 dipsort(nd,3,*)

      real *8 targetsort(2,ntarget)
      complex *16 jsort(nd,5,0:ntj,*)

      real *8 expcsort(2,*)

      complex *16 pot(nd,*)
      complex *16 grad(nd,3,*)
      complex *16 hess(nd,3,*)

      complex *16 pottarg(nd,*)
      complex *16 gradtarg(nd,3,*)
      complex *16 hesstarg(nd,3,*)

      integer iaddr(2,nboxes),lmptmp
      real *8 rmlexp(*)
      complex *16 mptemp(lmptmp)
       
      real *8 timeinfo(10)
      real *8 timelev(0:200)
      real *8 centers(2,*)

      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer iptr(8),ltree
      integer itree(ltree)
      integer nboxes
      integer isrcse(2,nboxes),itargse(2,nboxes)
      integer iexpcse(2,nboxes)
      real *8 rscales(0:nlevels),boxsize(0:nlevels)

      real *8 scjsort(*)

      real *8 thresh

      integer nterms_eval(4,0:200)

c     temp variables
      integer i,j,k,l,idim
      integer ibox,jbox,ilev,npts
      integer nchild,nlist1,nlist2,nlist3,nlist4
      integer mnlist1,mnlist2,mnlist3,mnlist4
      integer, allocatable :: list1(:,:),list2(:,:),list3(:,:)
      integer, allocatable :: list4(:,:)
      integer, allocatable :: nlist1s(:),nlist2s(:),nlist3s(:)
      integer, allocatable :: nlist4s(:)

      integer istart,iend,istarts,iends
      integer isstart,isend,jsstart,jsend
      integer jstart,jend
      integer istarte,iende,istartt,iendt

      integer ifprint

      integer ifhesstarg,nn
      real *8 d,time1,time2,omp_get_wtime
      real *8 tt1,tt2
      complex *16 pottmp,gradtmp,hesstmp
      
      real *8, allocatable :: carray(:,:)
      integer ldc

      double precision pi
      
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c      
        ifprint=0

        pi = 4*atan(1.0d0)
c

        do i=0,nlevels
          timelev(i) = 0
        enddo

        ldc = 100
        allocate(carray(0:ldc,0:ldc))

        call init_carray(carray,ldc)

c
c        compute list info
c
        call computemnlists(nlevels,nboxes,itree,ltree,iptr,centers,
     1    boxsize,iper,mnlist1,mnlist2,mnlist3,mnlist4)
        allocate(nlist1s(nboxes),list1(mnlist1,nboxes))
        allocate(nlist2s(nboxes),list2(mnlist2,nboxes))
        allocate(nlist3s(nboxes),list3(mnlist3,nboxes))
        allocate(nlist4s(nboxes),list4(mnlist4,nboxes))
        
        call computelists(nlevels,nboxes,itree,ltree,iptr,centers,
     1    boxsize,iper,mnlist1,nlist1s,list1,mnlist2,nlist2s,list2,
     2    mnlist3,nlist3s,list3,mnlist4,nlist4s,list4)
c
c
c     ... set the expansion coefficients to zero
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(idim,i,j)
      do i=1,nexpc
         do j = 0,ntj
           do idim=1,nd
             jsort(idim,1,j,i)=0
             jsort(idim,2,j,i)=0
             jsort(idim,3,j,i)=0
             jsort(idim,4,j,i)=0
             jsort(idim,5,j,i)=0
           enddo
         enddo
      enddo
C$OMP END PARALLEL DO
C
c       
        do i=1,10
          timeinfo(i)=0
        enddo
c
c       ... set all multipole and local expansions to zero
c
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            call bh2dmpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
            call bh2dmpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
         enddo
C$OMP END PARALLEL DO         
       enddo

c     Set scjsort
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,i)
         do ibox = laddr(1,ilev), laddr(2,ilev)
            nchild = itree(iptr(4)+ibox-1)
            if(nchild.eq.0) then
                istart = iexpcse(1,ibox)
                iend = iexpcse(2,ibox)
                do i=istart,iend
                   scjsort(i) = rscales(ilev)
                enddo
            endif
         enddo
C$OMP END PARALLEL DO         
      enddo
       
c
c
      if(ifprint .ge. 1) 
     $   call prinf('=== STEP 1 (form mp) ====*',i,0)
        call cpu_time(time1)
C$        time1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions

      do ilev = 2,nlevels
C
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox=laddr(1,ilev),laddr(2,ilev)
             nchild = itree(iptr(4)+ibox-1)
             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend-istart+1
c              Check if current box is a leaf box            
             if(nchild.eq.0.and.npts.gt.0) then
                 call bh2dformmpc(nd,rscales(ilev),
     1             sourcesort(1,istart),npts,chargesort(1,1,istart),
     2             centers(1,ibox),nterms(ilev),
     3             rmlexp(iaddr(1,ibox)))
             endif
          enddo
C$OMP END PARALLEL DO 
        endif

        if(ifdipole.eq.1.and.ifcharge.eq.0) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox=laddr(1,ilev),laddr(2,ilev)
             nchild = itree(iptr(4)+ibox-1)
             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend-istart+1
c              Check if current box is a leaf box            
             if(nchild.eq.0.and.npts.gt.0) then
                call bh2dformmpd(nd,rscales(ilev),
     1          sourcesort(1,istart),npts,dipsort(1,1,istart),
     2          centers(1,ibox),
     3          nterms(ilev),rmlexp(iaddr(1,ibox))) 
             endif
          enddo
C$OMP END PARALLEL DO 
        endif

        if(ifdipole.eq.1.and.ifcharge.eq.1) then
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
          do ibox=laddr(1,ilev),laddr(2,ilev)
             nchild = itree(iptr(4)+ibox-1)
             istart = isrcse(1,ibox)
             iend = isrcse(2,ibox)
             npts = iend-istart+1
c             Check if current box is a leaf box            
             if(nchild.eq.0.and.npts.gt.0) then
                call bh2dformmpcd(nd,rscales(ilev),
     1             sourcesort(1,istart),npts,chargesort(1,1,istart),
     2             dipsort(1,1,istart),
     3             centers(1,ibox),
     4             nterms(ilev),rmlexp(iaddr(1,ibox))) 
             endif
          enddo
C$OMP END PARALLEL DO 
        endif
      enddo


      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(1)=time2-time1

      if(ifprint.ge.1)
     $      call prinf('=== STEP 2 (form lo) ====*',i,0)
      call cpu_time(time1)
C$        time1=omp_get_wtime()
      do ilev = 2,nlevels
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
           do ibox = laddr(1,ilev),laddr(2,ilev)
              npts = 0
              if(ifpghtarg.gt.0) then
                 istart = itargse(1,ibox)
                 iend = itargse(2,ibox)
                 npts = npts + iend-istart+1
              endif

              istart = iexpcse(1,ibox)
              iend = iexpcse(2,ibox)
              npts = npts + iend-istart+1

              if(ifpgh.gt.0) then
                 istart = isrcse(1,ibox)
                 iend = isrcse(2,ibox)
                 npts = npts + iend-istart+1
              endif
              
              if (npts .gt. 0) then
                 do i=1,nlist4s(ibox)
                    jbox = list4(i,ibox)
                    istart = isrcse(1,jbox)
                    iend = isrcse(2,jbox)
                    npts = iend-istart+1
                    
                    call bh2dformtac(nd,rscales(ilev),
     1                   sourcesort(1,istart),npts,
     2                   chargesort(1,1,istart),centers(1,ibox),
     3                   nterms(ilev),rmlexp(iaddr(2,ibox)))
                 enddo
              endif
           enddo
C$OMP END PARALLEL DO        
        endif
        if(ifcharge.eq.0.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
           do ibox = laddr(1,ilev),laddr(2,ilev)
              npts = 0
              if(ifpghtarg.gt.0) then
                 istart = itargse(1,ibox)
                 iend = itargse(2,ibox)
                 npts = npts + iend-istart+1
              endif

              istart = iexpcse(1,ibox)
              iend = iexpcse(2,ibox)
              npts = npts + iend-istart+1

              if(ifpgh.gt.0) then
                 istart = isrcse(1,ibox)
                 iend = isrcse(2,ibox)
                 npts = npts + iend-istart+1
              endif
              
              if (npts .gt. 0) then
                 do i=1,nlist4s(ibox)
                    jbox = list4(i,ibox)
                    istart = isrcse(1,jbox)
                    iend = isrcse(2,jbox)
                    npts = iend-istart+1

                    call bh2dformtad(nd,rscales(ilev),
     1                   sourcesort(1,istart),npts,
     2                   dipsort(1,1,istart),
     3                   centers(1,ibox),nterms(ilev),
     4                   rmlexp(iaddr(2,ibox)))
                 enddo
              endif
          enddo
C$OMP END PARALLEL DO        
        endif
        if(ifcharge.eq.1.and.ifdipole.eq.1) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,nlist4,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
           do ibox = laddr(1,ilev),laddr(2,ilev)
              npts = 0
              if(ifpghtarg.gt.0) then
                 istart = itargse(1,ibox)
                 iend = itargse(2,ibox)
                 npts = npts + iend-istart+1
              endif

              istart = iexpcse(1,ibox)
              iend = iexpcse(2,ibox)
              npts = npts + iend-istart+1

              if(ifpgh.gt.0) then
                 istart = isrcse(1,ibox)
                 iend = isrcse(2,ibox)
                 npts = npts + iend-istart+1
              endif
              
              if (npts .gt. 0) then              
                 do i=1,nlist4s(ibox)
                    jbox = list4(i,ibox)
                    istart = isrcse(1,jbox)
                    iend = isrcse(2,jbox)
                    npts = iend-istart+1
                    
                    call bh2dformtacd(nd,rscales(ilev),
     1                   sourcesort(1,istart),npts,
     2                   chargesort(1,1,istart),dipsort(1,1,istart),
     3                   centers(1,ibox),
     3                   nterms(ilev),rmlexp(iaddr(2,ibox)))
                 enddo
              endif
          enddo
C$OMP END PARALLEL DO        
        endif
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1

cc      call prin2('carray=*',carray,(ldc+1)*(ldc+1))

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 3 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
      do ilev=nlevels-1,1,-1

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          nchild = itree(iptr(4)+ibox-1)
          do i=1,nchild
            jbox = itree(iptr(5)+4*(ibox-1)+i-1)
            istart = isrcse(1,jbox)
            iend = isrcse(2,jbox)
            npts = iend-istart+1
            if(npts.gt.0) then
              call bh2dmpmp(nd,rscales(ilev+1),
     1             centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2             nterms(ilev+1),rscales(ilev),centers(1,ibox),
     3             rmlexp(iaddr(1,ibox)),nterms(ilev),carray,ldc)
            endif
          enddo
        enddo
C$OMP END PARALLEL DO    
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(3)=time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 4 (mp to loc) ===*',i,0)
c      ... step 3, convert multipole expansions into local
c       expansions

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      do ilev = 2,nlevels

       tt1 = second()
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istart,iend,npts,i,nlist2)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          npts = 0
          if(ifpghtarg.gt.0) then
            istart = itargse(1,ibox)
            iend = itargse(2,ibox)
            npts = npts + iend-istart+1
          endif

          istart = iexpcse(1,ibox)
          iend = iexpcse(2,ibox)
          npts = npts + iend-istart+1

          if(ifpgh.gt.0) then
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = npts + iend-istart+1
          endif

          if(npts.gt.0) then
            do i=1,nlist2s(ibox)
              jbox = list2(i,ibox) 
              call bh2dmploc(nd,rscales(ilev),
     $          centers(1,jbox),rmlexp(iaddr(1,jbox)),nterms(ilev),
     2          rscales(ilev),centers(1,ibox),rmlexp(iaddr(2,ibox)),
     3          nterms(ilev),carray,ldc)
            enddo
          endif
        enddo
C$OMP END PARALLEL DO        
       tt2 = second()
       timelev(ilev) = tt2-tt1
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(4) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 5 (split loc) ===*',i,0)

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      do ilev = 1,nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          nchild = itree(iptr(4)+ibox-1)
          istart = iexpcse(1,ibox)
          iend = iexpcse(2,ibox) 
          npts = iend - istart + 1


          if(ifpghtarg.gt.0) then
            istart = itargse(1,ibox) 
            iend = itargse(2,ibox) 
            npts = npts + iend-istart+1
          endif

          if(ifpgh.gt.0) then
            istart = isrcse(1,ibox) 
            iend = isrcse(2,ibox) 
            npts = npts + iend-istart+1
          endif

          if(npts.gt.0) then
            do i=1,nchild
              jbox = itree(iptr(5)+4*(ibox-1)+i-1)
              call bh2dlocloc(nd,rscales(ilev),centers(1,ibox),
     1          rmlexp(iaddr(2,ibox)),nterms(ilev),rscales(ilev+1),
     2          centers(1,jbox),rmlexp(iaddr(2,jbox)),nterms(ilev+1),
     3          carray,ldc)
            enddo
          endif
        enddo
C$OMP END PARALLEL DO        
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(5) = time2-time1

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      if(ifprint.ge.1)
     $    call prinf('=== Step 6 (mp eval) ===*',i,0)

cc      call prinf('ifpgh=*',ifpgh,1)
cc      call prinf('ifpghtarg=*',ifpghtarg,1)
cc      call prinf('laddr=*',laddr,2*(nlevels+1))
      do ilev=1,nlevels-1
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,nlist3,istart,iend,npts,j,i)
C$OMP$PRIVATE(jbox)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox=laddr(1,ilev),laddr(2,ilev)
c              evalute multipole expansion at all targets
          istart = itargse(1,ibox)
          iend = itargse(2,ibox) 
          npts = iend-istart+1

          if(ifpghtarg.eq.1) then
            do i=1,nlist3s(ibox)
              jbox = list3(i,ibox) 
                  
              call bh2dmpevalp(nd,rscales(ilev+1),
     1         centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2         nterms(ilev+1),targetsort(1,istart),npts,
     3         pottarg(1,istart))
            enddo
          endif
          if(ifpghtarg.eq.2) then
            do i=1,nlist3s(ibox)
              jbox = list3(i,ibox)
              call bh2dmpevalg(nd,rscales(ilev+1),
     1          centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2          nterms(ilev+1),targetsort(1,istart),npts,
     3          pottarg(1,istart),gradtarg(1,1,istart))
            enddo
          endif


c              evalute multipole expansion at all sources
          istart = isrcse(1,ibox)
          iend = isrcse(2,ibox) 
          npts = iend-istart+1
            

          if(ifpgh.eq.1) then
            do i=1,nlist3s(ibox)
              jbox = list3(i,ibox) 
              call bh2dmpevalp(nd,rscales(ilev+1),
     1           centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2           nterms(ilev+1),sourcesort(1,istart),npts,
     3           pot(1,istart))
            enddo
          endif
          if(ifpgh.eq.2) then
            do i=1,nlist3s(ibox)
              jbox = list3(i,ibox) 
              call bh2dmpevalg(nd,rscales(ilev+1),
     1           centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2           nterms(ilev+1),sourcesort(1,istart),npts,
     3           pot(1,istart),grad(1,1,istart))
            enddo
          endif

        enddo
C$OMP END PARALLEL DO     
      enddo

 1000 continue    


      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(6) = time2-time1



      if(ifprint.ge.1)
     $    call prinf('=== step 7 (eval lo) ===*',i,0)

c     ... step 7, evaluate all local expansions
      call cpu_time(time1)
C$    time1=omp_get_wtime()
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,i,npts,nchild)
C$OMP$SCHEDULE(DYNAMIC)
        do ibox = laddr(1,ilev),laddr(2,ilev)
          nchild = itree(iptr(4)+ibox-1)
          if(nchild.eq.0) then
c
cc               evaluate local expansion
c                at targets
            istart = itargse(1,ibox) 
            iend = itargse(2,ibox)
            npts = iend-istart + 1
            if(ifpghtarg.eq.1) then
              call bh2dtaevalp(nd,rscales(ilev),
     1              centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2              nterms(ilev),targetsort(1,istart),npts,
     3              pottarg(1,istart))
            endif
            if(ifpghtarg.eq.2) then
              call bh2dtaevalg(nd,rscales(ilev),
     1          centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2          nterms(ilev),targetsort(1,istart),npts,
     3          pottarg(1,istart),gradtarg(1,1,istart))
            endif

c
cc                evaluate local expansion at sources

            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1
            if(ifpgh.eq.1) then
              call bh2dtaevalp(nd,rscales(ilev),
     1           centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2           nterms(ilev),sourcesort(1,istart),npts,
     3           pot(1,istart))
            endif
            if(ifpgh.eq.2) then
              call bh2dtaevalg(nd,rscales(ilev),
     1           centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2           nterms(ilev),sourcesort(1,istart),npts,
     3           pot(1,istart),grad(1,1,istart))
            endif
          endif
        enddo
C$OMP END PARALLEL DO        
      enddo

      call cpu_time(time2)
C$    time2 = omp_get_wtime()      
      timeinfo(7) = time2 - time1

      if(ifprint .ge. 1)
     $     call prinf('=== STEP 8 (direct) =====*',i,0)

c
cc     set threshold for ignoring interactions with 
c      |r| < thresh
c
      thresh = boxsize(0)*2.0d0**(-51)

c
cc
      call cpu_time(time1)
C$    time1=omp_get_wtime()  
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istartt,iendt,i,jstart,jend,istarte,iende)
C$OMP$PRIVATE(nlist1,istarts,iends)
C$OMP$SCHEDULE(DYNAMIC)  
         do ibox = laddr(1,ilev),laddr(2,ilev)

            istartt = itargse(1,ibox)
            iendt = itargse(2,ibox)


            istarte = iexpcse(1,ibox)
            iende = iexpcse(2,ibox)

            istarts = isrcse(1,ibox)
            iends = isrcse(2,ibox)

            do i =1,nlist1s(ibox)
               jbox = list1(i,ibox) 

               jstart = isrcse(1,jbox)
               jend = isrcse(2,jbox)

                
               call bhfmm2dpart_direct(nd,jstart,jend,istartt,
     1         iendt,sourcesort,ifcharge,chargesort,ifdipole,
     2         dipsort,targetsort,ifpghtarg,pottarg,
     3         gradtarg,hesstarg,thresh)
         
               call bhfmm2dpart_direct(nd,jstart,jend,istarts,iends,
     1         sourcesort,ifcharge,chargesort,ifdipole,
     2         dipsort,sourcesort,ifpgh,pot,grad,hess,
     3         thresh)
            enddo   
         enddo
C$OMP END PARALLEL DO         
      enddo

      call cpu_time(time2)
C$    time2=omp_get_wtime()  
      timeinfo(8) = time2-time1
      if(ifprint.ge.1) call prin2('timeinfo=*',timeinfo,8)
      d = 0
      do i = 1,8
         d = d + timeinfo(i)
      enddo

      if(ifprint.ge.1) call prin2('sum(timeinfo)=*',d,1)
      if(ifprint.ge.1) call prin2('timlev=*',timelev,nlevels+1)

      return
      end
c
c
c
c
c
c------------------------------------------------------------------     
      subroutine bhfmm2dpart_direct(nd,istart,iend,jstart,jend,
     $     source,ifcharge,charge,ifdipole,dip,
     $     targ,ifpgh,pot,grad,hess,thresh)
c--------------------------------------------------------------------
c     This subroutine adds the contribuition due to sources
c     istart to iend in the source array at the expansion centers
c     jstart to jend in the target array to the computed velocities
c     and gradients. Note that contributions for sources
c     within thresh of the targets are not added to the potential
c     
c
c     INPUT arguments
c-------------------------------------------------------------------
c     nd           in: integer
c                  number of charge densities
c
c     istart       in:Integer
c                  Starting index in source array whose expansions
c                  we wish to add
c
c     iend         in:Integer
c                  Last index in source array whose expansions
c                  we wish to add
c
c     jstart       in: Integer
c                  First index in target array at which we
c                  wish to update the potential and gradients
c 
c     jend         in:Integer
c                  Last index in target array at which we wish
c                  to update the potential and gradients
c
c     source       in: real *8(2,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: complex *16 (nd,2,ns)
c                  Charge at the source locations
c
c     ifdipole     in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dip         in: complex *16(nd,3,ns)
c                 dipole strengths at the source locations
c
c     targ        in: real *8(2,nt)
c                 target locations
c
c     ifpgh        in: Integer
c                  Flag for computing the potential/gradient/hessian.
c                  ifpgh = 1, only potential is computed
c                  ifpgh = 2, potential/gradient are computed
c                  ifpgh = 3, potential/gradient/hessian are computed
c
c     thresh       in: real *8
c                  threshold for computing interactions
c                  if |r| < threshold, then interactions are
c                  not included
c
c
c------------------------------------------------------------
c     OUTPUT
c
c   Updated velocity and gradients at the targets
c   pot : potential at the targets
c   grad: gradient at the targets (d/dz_projz, d/dz_projzbar,d/dzbar)
c   hess: Hessian at the targets
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j,i
        integer ifcharge,ifdipole,nt

        integer nd



        real *8 source(2,*)
        complex *16 charge(nd,2,*),dip(nd,3,*)

        integer ifpgh
        real *8 targ(2,*),thresh
        
c
        complex *16 pot(nd,*)
        complex *16 grad(nd,3,*)
        complex *16 hess(nd,3,*)

c
        ns = iend - istart + 1
        nt = jend - jstart + 1
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          if(ifpgh.eq.1) then
             call bh2d_directcp(nd,source(1,istart),ns,
     1        charge(1,1,istart),targ(1,jstart),nt,pot(1,jstart),thresh)
          endif

          if(ifpgh.eq.2) then
            call bh2d_directcg(nd,source(1,istart),ns,
     1         charge(1,1,istart),targ(1,jstart),nt,pot(1,jstart),
     1         grad(1,1,jstart),thresh)
          endif
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             call bh2d_directdp(nd,source(1,istart),ns,
     1         dip(1,1,istart),
     2         targ(1,jstart),nt,pot(1,jstart),thresh)
          endif

          if(ifpgh.eq.2) then
               call bh2d_directdg(nd,source(1,istart),ns,
     1            dip(1,1,istart),
     2            targ(1,jstart),nt,pot(1,jstart),
     2            grad(1,1,jstart),thresh)
          endif
        endif

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
               call bh2d_directcdp(nd,source(1,istart),ns,
     1            charge(1,1,istart),dip(1,1,istart),
     2            targ(1,jstart),nt,pot(1,jstart),thresh)
          endif

          if(ifpgh.eq.2) then
               call bh2d_directcdg(nd,source(1,istart),ns,
     1            charge(1,1,istart),dip(1,1,istart),
     2            targ(1,jstart),nt,pot(1,jstart),
     2            grad(1,1,jstart),
     2            thresh)
          endif
        endif


c
        return
        end
c------------------------------------------------------------------    
      subroutine bh2dmpalloc(nd,laddr,iaddr,nlevels,lmptot,
     1                          nterms)
c     This subroutine determines the size of the array
c     to be allocated for the multipole expansions
c     iaddr(1,i) points to the starting location of the multipole
c     expansion of box i and iaddr(2,i) points to the local
c     expansion of box i
c  
c     Input arguments
c     nd          in: integer
c                 number of expansions
c
c     laddr       in: Integer(2,0:nlevels)
c                 indexing array providing access to boxes at each
c                 level
c
c     nlevels     in: Integer
c                 Total numner of levels
c     
c     nterms      in: Integer(0:nlevels)
c                 Number of terms requried in expansions at each
c                 level
c
c------------------------------------------------------------------
c     Output arguments
c     iaddr       out: Integer(2,nboxes)
c                 Points the multipole and local expansions in box i
c 
c     lmptot      out: Integer
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,nterms(0:nlevels),nd,nsig,nt1,nt2,next235
      integer iaddr(2,*), lmptot, laddr(2,0:nlevels)
      integer ibox,i,iptr,istart,nn,itmp
      real *8 ddn
c
      istart = 1
      do i = 0,nlevels

         nn = (nterms(i)+1)*2*nd*5
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole expansion         
c
           itmp = ibox - laddr(1,i)
           iaddr(1,ibox) = istart + itmp*nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
       enddo
c
c            Allocate memory for the local expansion
c
       do i=0,nlevels
         nn = (nterms(i)+1)*2*nd*5
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
             itmp = ibox - laddr(1,i)
             iaddr(2,ibox) = istart + itmp*nn 
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
      lmptot = istart

      return
      end
c----------------------------------------------------------------     
