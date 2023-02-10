cc Copyright (C) 2018-2019: Leslie Greengard, Zydrunas Gimbutas, 
cc and Manas Rachh
cc Contact: greengard@cims.nyu.edu
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
      subroutine hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   ifcharge      : flag for including charge interactions
c                   charge interactions included if ifcharge =1
c                   not included otherwise
c   charge(nd,ns)    : charge strengths
c   ifdipole      : flag for including dipole interactions
c                   dipole interactions included if ifcharge =1
c                   not included otherwise
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orienstation vectors
c   iper          : flag for periodic implmentations. Currently unused
c   ifpgh         : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed
c   nt            : number of targets
c   targ(2,nt)    : target locations
c   ifpghtarg     : flag for computing pottarg/gradtarg/hesstarg
c                   ifpghtarg = 1, only potential is computed at targets
c                   ifpghtarg = 2, potential and gradient are 
c                   computed at targets
c                   ifpghtarg = 3, potential, gradient, and hessian are 
c                   computed at targets
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : potential at the source locations
c   grad(nd,2,*)    : gradients at the source locations
c   hess(nd,3,*)    : hessian at the source locations
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,2,*): gradient at the target locations
c   hesstarg(nd,3,*): hessian at the target locations
c   ier             : error code
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      integer iper
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,*)
      complex *16 charge(nd,*),dipstr(nd,*)

      complex *16 pot(nd,*),grad(nd,2,*),hess(nd,3,*)
      complex *16 pottarg(nd,*),gradtarg(nd,2,*),hesstarg(nd,3,*)

c
cc      Tree variables
c
      integer, allocatable :: itree(:)
      integer iptr(8)
      integer nlmin,nlmax,ifunif
      real *8, allocatable :: tcenters(:,:),boxsize(:)
      integer nexpc,ntj
      real *8 expc(2)
      real *8 scj
      complex *16 jexps(100)
      integer idivflag,nlevels,nboxes,ndiv
      integer ltree

      real *8, allocatable :: radsrc(:)
      real *8 radexp

c
cc     sorted arrays
c
      integer, allocatable :: isrc(:),isrcse(:,:)
      integer, allocatable :: itarg(:),itargse(:,:),iexpcse(:,:)
      real *8, allocatable :: sourcesort(:,:)
      real *8, allocatable :: targsort(:,:)
      complex *16, allocatable :: chargesort(:,:),dipstrsort(:,:)
      real *8, allocatable :: dipvecsort(:,:,:)
      complex *16, allocatable :: potsort(:,:),gradsort(:,:,:),
     1                             hesssort(:,:,:)
      complex *16, allocatable :: pottargsort(:,:),gradtargsort(:,:,:),
     1                              hesstargsort(:,:,:)

c
cc     additional fmm variables

      integer lmptot
      real *8, allocatable :: rscales(:)
      integer, allocatable :: nterms(:),iaddr(:,:)
      real *8, allocatable :: rmlexp(:)
      complex *16, allocatable :: mptemp(:)

      real *8 timeinfo(8)
      integer ifnear

c
cc      temporary variables
c
      integer i,ilev,lmptmp,nmax,idim
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg,ifprint,iert
      real *8 time1,time2,pi,done
      real *8 omp_get_wtime

      done = 1
      pi = atan(done)*4.0d0


      nexpc = 0

      nlevels = 0
      nboxes = 0



      call hndiv2d(eps,ns,nt,ifcharge,ifdipole,ifpgh,
     1  ifpghtarg,ndiv,idivflag)

      ltree = 0
      nlmin = 0
      nlmax = 51
      ifunif = 0
      iper = 0

      ifprint = 1

c
c  turn on computation of list 1
c
      ifnear = 1

c
c  initialize timeinfo
c
      do i=1,8
        timeinfo(i) = 0
      enddo

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
     1   tcenters,itarg,itargse)
      allocate(sourcesort(2,ns))
      allocate(targsort(2,nt))


      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        allocate(chargesort(nd,ns),dipstrsort(nd,1),dipvecsort(nd,2,1))
      endif
      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        allocate(chargesort(nd,1),dipstrsort(nd,ns),dipvecsort(nd,2,ns))
      endif
      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        allocate(chargesort(nd,ns),dipstrsort(nd,ns),
     1     dipvecsort(nd,2,ns))
      endif

      if(ifpgh.eq.1) then
        allocate(potsort(nd,ns),gradsort(nd,2,1),hesssort(nd,3,1))
      else if(ifpgh.eq.2) then
        allocate(potsort(nd,ns),gradsort(nd,2,ns),hesssort(nd,3,1))
      else if(ifpgh.eq.3) then
        allocate(potsort(nd,ns),gradsort(nd,2,ns),hesssort(nd,3,ns))
      else
        allocate(potsort(nd,1),gradsort(nd,2,1),hesssort(nd,3,1))
      endif

      
      if(ifpghtarg.eq.1) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,2,1),
     1     hesstargsort(nd,3,1))
      else if(ifpghtarg.eq.2) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,2,nt),
     1      hesstargsort(nd,3,1))
      else if(ifpghtarg.eq.3) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,2,nt),
     1     hesstargsort(nd,3,nt))
      else
        allocate(pottargsort(nd,1),gradtargsort(nd,2,1),
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
          enddo
        enddo
      endif

      if(ifpgh.eq.3) then
        do i=1,ns
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,1,i) = 0
            gradsort(idim,2,i) = 0
            hesssort(idim,1,i) = 0
            hesssort(idim,2,i) = 0
            hesssort(idim,3,i) = 0
          enddo
        enddo
      endif


      if(ifpghtarg.eq.1) then
        do i=1,nt
          do idim=1,nd
            pottarg(idim,i) = 0
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
          enddo
        enddo
      endif

      if(ifpghtarg.eq.3) then
        do i=1,nt
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,1,i) = 0
            gradtargsort(idim,2,i) = 0
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
      ier = 0
      do i=0,nlevels
        rscales(i) = min(abs(zk*boxsize(i)/(2.0d0*pi)),1.0d0)
        call h2dterms(boxsize(i),zk,eps,nterms(i),ier)
        nterms(i) = nterms(i) 
        if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo

      if(ifprint.eq.1) call prinf('nmax=*',nmax,1)
      if(ifprint.eq.1) call prinf('nterms=*',nterms,nlevels+1)

c       
c     Multipole and local expansions and diag forms
c     will be held in workspace
c     in locations pointed to by array iaddr(4,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c

      allocate(iaddr(4,nboxes))

      lmptmp = (2*nmax+1)*nd
      allocate(mptemp(lmptmp))

c     reorder sources
c
      call dreorderf(2,ns,sources,sourcesort,isrc)
      if(ifcharge.eq.1) 
     1    call dreorderf(2*nd,ns,charge,chargesort,isrc)
      if(ifdipole.eq.1) then
         call dreorderf(2*nd,ns,dipstr,dipstrsort,isrc)
         call dreorderf(2*nd,ns,dipvec,dipvecsort,isrc)
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
      call h2dmpalloc(nd,itree(iptr(1)),iaddr,nlevels,lmptot,
     1    nterms)
      if(ifprint .eq. 1) call prinf(' lmptot is *',lmptot,1)

      ier = 0
      allocate(rmlexp(lmptot),stat=iert)
      if(iert.ne.0) then
         print *, "Cannot allocate mpole expansion workspace"
         print *, "lmptot=", lmptot
         ier = 4
         return
      endif



c
cc     call the main fmm routine
c

c     Memory allocation is complete. 
c     Call main fmm routine
c
      call cpu_time(time1)
C$      time1=omp_get_wtime()
      call hfmm2dmain(nd,eps,
     $   zk,ns,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,dipstrsort,dipvecsort,
     $   nt,targsort,nexpc,expc,
     $   iaddr,rmlexp,mptemp,lmptmp,
     $   itree,ltree,iptr,ndiv,nlevels,
     $   nboxes,iper,boxsize,rscales,tcenters,itree(iptr(1)),
     $   isrcse,itargse,iexpcse,nterms,ntj,
     $   ifpgh,potsort,gradsort,hesssort,
     $   ifpghtarg,pottargsort,gradtargsort,
     $   hesstargsort,jexps,scj,ifnear,timeinfo,ier)
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)
      if(ier.ne.0) return


c
cc      resort the output arrays in input order
c

      if(ifpgh.eq.1) then
        call dreorderi(2*nd,ns,potsort,pot,isrc)
      endif

      if(ifpgh.eq.2) then
        call dreorderi(2*nd,ns,potsort,pot,isrc)
        call dreorderi(4*nd,ns,gradsort,grad,isrc)
      endif

      if(ifpgh.eq.3) then
        call dreorderi(2*nd,ns,potsort,pot,isrc)
        call dreorderi(4*nd,ns,gradsort,grad,isrc)
        call dreorderi(6*nd,ns,hesssort,hess,isrc)
      endif

cc      call prini(6,13)
cc      call prin2('eps = *', eps, 1)
cc      call prin2('after hfmm2dmain, pottargsort = *', pottargsort, 30)
cc      stop
      
      if(ifpghtarg.eq.1) then
        call dreorderi(2*nd,nt,pottargsort,pottarg,itarg)
      endif

      if(ifpghtarg.eq.2) then
        call dreorderi(2*nd,nt,pottargsort,pottarg,itarg)
        call dreorderi(4*nd,nt,gradtargsort,gradtarg,itarg)
      endif

      if(ifpghtarg.eq.3) then
        call dreorderi(2*nd,nt,pottargsort,pottarg,itarg)
        call dreorderi(4*nd,nt,gradtargsort,gradtarg,itarg)
        call dreorderi(6*nd,nt,hesstargsort,hesstarg,itarg)
      endif


      return
      end

      subroutine hfmm2dmain(nd,eps,
     $     zk,nsource,sourcesort,
     $     ifcharge,chargesort,
     $     ifdipole,dipstrsort,dipvecsort,
     $     ntarget,targetsort,nexpc,expcsort,
     $     iaddr,rmlexp,mptemp,lmptmp,
     $     itree,ltree,iptr,ndiv,nlevels, 
     $     nboxes,iper,boxsize,rscales,centers,laddr,
     $     isrcse,itargse,iexpcse,nterms,ntj,
     $     ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg,
     $     jsort,scjsort,ifnear,timeinfo,ier)
c   Helmholtz FMM in R^2: evaluate all pairwise particle
c   interactions (ignoring self-interaction) 
c   and interactions with targets.
c
c   We use H_0(kr)*(i/4) for the Green's function.
c   Self-interactions are not included
c
c   h2d: charge and dipstr are complex valued, x in \R^2
c
c   \phi(x_i) = (i/4)\sum_{j\ne i} charge_j H^{(1)}_0(k |x_i - x_j|)
c   + dipstr_j (dipvec_j \dot (x_i - x_j)) H^{(1)}_1(k |x_i - x_j|*
c                                          k/|x_i-x_j|
c
c   All the source/target/expansion center related quantities
c   are assumed to be tree-sorted
c
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
c
c   nd:   number of charge densities
c
c   eps:  FMM precision requested
c
c   zk: complex *16, Helmholtz parameter
c
c   nsource:     integer:  number of sources
c   sourcesort: real *8 (2,ns):  source locations
c
c   ifcharge:  charge computation flag
c              ifcharge = 1   =>  include charge contribution
c                                     otherwise do not
c   chargesort: complex *16 (nsource): charge strengths
c
c   ifdipole:  dipole computation flag
c              ifdipole = 1   =>  include dipole contribution
c                                     otherwise do not
c   dipstrsort: complex *16 (nsource): dipole strengths
c   dipvecsort: real *8 (2,nsource): dip orientation
c   ntarget: integer:  number of targets
c   targetsort: real *8 (2,ntarget):  target locations
c   nexpc: number of expansion centers
c   expcsort: real *8 (2,nexpc): expansion center locations
c   iaddr: (4,nboxes): pointer in rmlexp where multipole
c                      and local expansions and diag forms for each
c                      box are stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      multipole expansion of ibox
c                      iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the local expansion of ibox
c                      iaddr(3,ibox) is the
c                      starting index in rmlexp
c                      for the outgoing diag form of ibox
c                      iaddr(4,ibox) is the
c                      starting index in rmlexp
c                      for the incoming diag form of ibox
c  mptemp: (lmptmp): temporary multipole/local expansion
c                        (may not be needed in new setting)
c  lmptmp: length of temporary expansion
c   
c
c
c   itree    in: integer (ltree)
c             This array contains all the information
c             about the tree
c             Refer to pts_tree2d.f 
c
c   ltree    in: integer
c            length of tree
c
c    iptr in: integer(8)
c             iptr is a collection of pointers 
c             which points to where different elements 
c             of the tree are stored in the itree array
c
c     ndiv    in: integer
c             Max number of chunks per box
c
c     nlevels in: integer
c             number of levels in the tree
c
c     
c     nboxes  in: integer
c             number of boxes in the tree
c
c     boxsize in: real*8 (0:nlevels)
c             boxsize(i) is the size of the box from end to end
c             at level i
c   
c     iper    in: integer
c             flag for periodic implementation. Currently unused
c
c     centers in: real *8(2,nboxes)
c                 array containing the centers of all the boxes
c
c     isrcse in: integer(2,nboxes)
c               starting and ending location of sources in ibox
c                in sorted list of sources
c
c     itargse in: integer(2,nboxes)
c               starting and ending location of targets in ibox
c                in sorted list of sources
c
c     iexpcse in: integer(2,nboxes)
c               starting and ending location of expansion centers
c               in ibox in sorted list of sources
c
c     nterms: (0:nlevels) length of multipole and local expansions
c              at various levels
c     ntj     in: integer
c             order of the output expansions
c
c     ifpgh  in: integer
c             flag for evaluating potential/gradients/hessians 
c             at sources.
c             ifpgh = 1, only potentials will be evaluated
c             ifpgh = 2, potentials/gradients will be evaluated
c             ifpgh = 3, potentials/gradients/hessians will be evaluated
c
c     ifpghtarg  in: integer
c             flag for evaluating potential/gradients/hessians 
c             at targets.
c             ifpghtarg = 1, only potentials will be evaluated
c             ifpghtarg = 2, potentials/gradients will be evaluated
c             ifpghtarg = 3, potentials/gradients/hessians will be evaluated
c
c
c   OUTPUT
c
c   Expansions at the targets
c   jexps : coeffs for local expansion
c   scj: scaling parameter for the expansions
c
c   pot: potential at the source locations
c   grad: gradient at the source locations
c   hess: gradient at the source locations
c  
c   pottarg: potential at the target locations
c   gradtarg: gradient at the target locations
c   hesstarg: gradient at the target locations
c------------------------------------------------------------------

      implicit none

      integer nd

      complex *16 zk
      real *8 zi

      integer nsource,ntarget,nexpc,ier
      integer ndiv,nlevels,ntj,ix,iy

      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      real *8 eps
      integer iper,ifnear

      real *8 sourcesort(2,nsource)

      complex *16 chargesort(nd,*)
      complex *16 dipstrsort(nd,*)
      real *8 dipvecsort(nd,2,*)

      real *8 targetsort(2,ntarget)
      complex *16 jsort(nd,-ntj:ntj,*)

      real *8 expcsort(2,*)

      complex *16 pot(nd,*)
      complex *16 grad(nd,2,*)
      complex *16 hess(nd,3,*)

      complex *16 pottarg(nd,*)
      complex *16 gradtarg(nd,2,*)
      complex *16 hesstarg(nd,3,*)

      integer iaddr(4,nboxes),lmptmp
      real *8 rmlexp(*)
      complex *16 mptemp(lmptmp)
       
      real *8 timeinfo(8)
      real *8 timelev(0:200)
      real *8 centers(2,*)

      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer iptr(8),ltree
      integer isrcse(2,nboxes),itargse(2,nboxes),iexpcse(2,nboxes)
      integer itree(*)
      integer nboxes
      real *8 rscales(0:nlevels),boxsize(0:nlevels)

      real *8 scjsort(*)

      real *8 zkiupbound,thresh,dn,dx,dy
      real *8 c1(2),c2(2)

      integer nterms_eval(4,0:200)

c     temp variables
      integer i,j,k,l,idim
      integer ibox,jbox,ilev,npts,next235,ilevhf
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
      complex *16 pottmp,gradtmp(2),hesstmp(3)

      integer :: ni, nsig
      double complex, allocatable :: wsave(:)
      double complex, allocatable :: transvecall(:,:,:)
      double complex, allocatable :: transvecmpmp(:,:)
      double complex, allocatable :: sig(:,:)
      double complex, allocatable :: sig2(:,:,:)
      
      double precision dlam, pi, boxlam
      
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c      
        ifprint=1

        pi = 4*atan(1.0d0)
c
cc
c           upper limit for zk along imaginary axis
        zkiupbound = 40.0d0
        zi = imag(zk)

        do i=0,nlevels
          timelev(i) = 0
        enddo

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
         do j = -ntj,ntj
           do idim=1,nd
             jsort(idim,j,i)=0
           enddo
         enddo
      enddo
C$OMP END PARALLEL DO
C
c       
        do i=1,8
          timeinfo(i)=0
        enddo
c
c       ... set all multipole and local expansions to zero
c
      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nn,dn)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            call h2dmpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
            call h2dmpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
            dn = 2*(nterms(ilev)+nterms(ilev)) + 1
            nn = next235(dn)
            call h2dsigzero(nd,rmlexp(iaddr(3,ibox)),nn)
            call h2dsigzero(nd,rmlexp(iaddr(4,ibox)),nn)
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
c
      if(ifprint .ge. 1) 
     $   call prinf('=== STEP 1 (form mp) ====*',i,0)
        call cpu_time(time1)
C$        time1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions

      do ilev = 2,nlevels
       if(zi*boxsize(ilev).lt.zkiupbound) then
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
                  call h2dformmpc(nd,zk,rscales(ilev),
     1               sourcesort(1,istart),npts,chargesort(1,istart),
     2               centers(1,ibox),nterms(ilev),
     3               rmlexp(iaddr(1,ibox)))
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
                  call h2dformmpd(nd,zk,rscales(ilev),
     1            sourcesort(1,istart),npts,dipstrsort(1,istart),
     2            dipvecsort(1,1,istart),centers(1,ibox),
     3            nterms(ilev),rmlexp(iaddr(1,ibox))) 
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
c              Check if current box is a leaf box            
               if(nchild.eq.0.and.npts.gt.0) then
                  call h2dformmpcd(nd,zk,rscales(ilev),
     1               sourcesort(1,istart),npts,chargesort(1,istart),
     2               dipstrsort(1,istart),
     3               dipvecsort(1,1,istart),centers(1,ibox),
     4               nterms(ilev),rmlexp(iaddr(1,ibox))) 
               endif
            enddo
C$OMP END PARALLEL DO 
         endif
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
       if(zi*boxsize(ilev).lt.zkiupbound) then
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

               if(npts.gt.0) then
                  
                  do i=1,nlist4s(ibox)
                     jbox = list4(i,ibox)
                     istart = isrcse(1,jbox)
                     iend = isrcse(2,jbox)
                     npts = iend-istart+1
                     
                     call h2dformtac(nd,zk,rscales(ilev),
     1                    sourcesort(1,istart),npts,
     2                    chargesort(1,istart),centers(1,ibox),
     3                    nterms(ilev),rmlexp(iaddr(2,ibox)))
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

               if(npts.gt.0) then
                  
                  do i=1,nlist4s(ibox)
                     jbox = list4(i,ibox)
                     istart = isrcse(1,jbox)
                     iend = isrcse(2,jbox)
                     npts = iend-istart+1

                     call h2dformtad(nd,zk,rscales(ilev),
     1                    sourcesort(1,istart),npts,
     2                    dipstrsort(1,istart),dipvecsort(1,1,istart),
     3                    centers(1,ibox),nterms(ilev),
     4                    rmlexp(iaddr(2,ibox)))
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

               if(npts.gt.0) then
                  
                  do i=1,nlist4s(ibox)
                     jbox = list4(i,ibox)
                     istart = isrcse(1,jbox)
                     iend = isrcse(2,jbox)
                     npts = iend-istart+1

                     call h2dformtacd(nd,zk,rscales(ilev),
     1                    sourcesort(1,istart),npts,
     2                    chargesort(1,istart),dipstrsort(1,istart),
     3                    dipvecsort(1,1,istart),centers(1,ibox),
     3                    nterms(ilev),rmlexp(iaddr(2,ibox)))
                  enddo
               endif
            enddo
C$OMP END PARALLEL DO        
         endif
       endif
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1
c
c     Identify level ilevhf where HF regine begins
c
      ilevhf = 0
      do ilev=nlevels-1,1,-1
         if(zi*boxsize(ilev).lt.zkiupbound) then
            dlam = zk
            dlam = 1/(dlam/(2*pi))                 
            boxlam = boxsize(ilev)/dlam
            if(boxlam.gt.16.0d0) then
               ilevhf = ilev
               goto 111
            endif
          endif
      enddo
111   continue
      if(ifprint.ge.1) call prinf(' ilevhf is *',ilevhf,1)

      if(ifprint .ge. 1)
     $      call prinf('=== STEP 3a (merge mp low freq) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
      do ilev=nlevels-1,ilevhf+1,-1
c
        if(zi*boxsize(ilev).lt.zkiupbound) then
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
                   call h2dmpmp(nd,zk,rscales(ilev+1),
     1               centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2               nterms(ilev+1),rscales(ilev),centers(1,ibox),
     3               rmlexp(iaddr(1,ibox)),nterms(ilev))
                 endif
              enddo
           enddo
C$OMP END PARALLEL DO    
        endif
      enddo
      if(ifprint .ge. 1)
     $      call prinf('=== STEP 3b (merge mp high freq) ====*',i,0)
c
c    convert to diag form at level ilev
c

      do ilev=ilevhf,1,-1

        if(zi*boxsize(ilev).lt.zkiupbound) then
           dn = 2*(nterms(ilev)+nterms(ilev+1)) + 1
           nsig = next235(dn)
           allocate(wsave(4*nsig+100))
           allocate(transvecmpmp(nsig,4))
           call zffti(nsig, wsave)
           if(ifprint.ge.1) print *, "Doing mpmp using hf"
           c2(1) = 0.0d0
           c2(2) = 0.0d0
           do jbox=1,4
              k=2
              if (jbox.le.2) k=1
              c1(1) = 0.25d0*boxsize(ilev)*(-1)**jbox
              c1(2) = 0.25d0*boxsize(ilev)*(-1)**k
              call h2d_mkmpshift(zk,c1,nterms(ilev+1),
     1           c2,nterms(ilev),nsig,wsave,transvecmpmp(1,jbox))
           enddo

           allocate(sig(nd,nsig))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,nchild,i)
C$OMP$PRIVATE(jbox,istart,iend,npts,sig)
           do ibox = laddr(1,ilev),laddr(2,ilev)
              nchild = itree(iptr(4)+ibox-1)
              call h2dsigzero(nd,sig,nsig)
              do i=1,nchild
                 jbox = itree(iptr(5)+4*(ibox-1)+i-1)
                 istart = isrcse(1,jbox)
                 iend = isrcse(2,jbox)
                 npts = iend-istart+1
                 if(npts.gt.0) then
                   call h2dmpmphf(nd,zk,rscales(ilev+1),
     1               centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2               nterms(ilev+1),rscales(ilev),centers(1,ibox),
     3               sig,nterms(ilev),nsig,wsave,
     4                transvecmpmp(1,i))
                 endif
              enddo
              call h2d_sig2exp(nd,nsig,sig,wsave,nterms(ilev),
     1          rmlexp(iaddr(1,ibox)))
           enddo
C$OMP END PARALLEL DO           
           deallocate(wsave)
           deallocate(transvecmpmp,sig)
        endif
      enddo
c
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

        ni = nterms(ilev)
        dn = 2*(ni + ni)+1
        nsig = next235(dn)
        allocate(wsave(4*nsig+100))        
        allocate(transvecall(nsig,-3:3,-3:3))
        dlam = zk
        dlam = 1/(dlam/(2*pi))                 
        boxlam = boxsize(ilev)/dlam
        call zffti(nsig,wsave)
        if(boxlam.gt.16.and.ifprint.ge.1) print *, "in high freq"
c
c   precompute mp to sig for all boxes 
c       
        if (boxlam .gt. 16.0d0) then
C$OMP PARALLEL DO DEFAULT(SHARED)          
           do ibox = laddr(1,ilev),laddr(2,ilev)
             call h2d_mptosig(nd,nterms(ilev),nsig,
     1            rmlexp(iaddr(1,ibox)),rmlexp(iaddr(3,ibox)),
     2            wsave)
           enddo
C$OMP END PARALLEL DO  

c
c  evaluate diagonal shifts for translation operators 
c
           c1(1) = 0.0d0
           c1(2) = 0.0d0
           do ix = -3,3
             do iy = -3,3
               c2(1) = ix*boxsize(ilev)
               c2(2) = iy*boxsize(ilev)
               call h2d_mkm2ltrans(zk,c1,nterms(ilev),
     1           c2,nterms(ilev),nsig,wsave,transvecall(1,ix,iy))
             enddo
           enddo
        endif

        call cpu_time(tt1)
C$    tt1=omp_get_wtime()
c
       if(zi*boxsize(ilev).lt.zkiupbound) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,istart,iend,npts,i,nlist2,ix,dx,iy,dy)
C$OMP$SCHEDULE(DYNAMIC)
c
         do ibox = laddr(1,ilev),laddr(2,ilev)
            npts = 0
c
c     if in high freq regime, transform multipole expansions
c     to diagonal form and store in rmlexp(3,*).
c
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

                  if (boxlam .gt. 16.0d0) then
                    dx = centers(1,ibox)-centers(1,jbox)
                    ix = nint(dx/boxsize(ilev))
                    dy = centers(2,ibox)-centers(2,jbox)
                    iy = nint(dy/boxsize(ilev))
                    call h2d_diagtrans(nd,nsig,rmlexp(iaddr(3,jbox)),
     1                   transvecall(1,ix,iy),rmlexp(iaddr(4,ibox)))
                  else
                    call h2dmploc(nd,zk,rscales(ilev),
     $                  centers(1,jbox),
     1                  rmlexp(iaddr(1,jbox)),nterms(ilev),
     2                  rscales(ilev),centers(1,ibox),
     3                  rmlexp(iaddr(2,ibox)),nterms(ilev))
                  endif
               enddo
            endif
         enddo
C$OMP END PARALLEL DO        
         if (boxlam .gt. 16.0d0) then
C$OMP PARALLEL DO DEFAULT(SHARED)         
            do ibox = laddr(1,ilev),laddr(2,ilev)
              call h2d_sig2exp(nd,nsig,
     3                rmlexp(iaddr(4,ibox)),wsave,nterms(ilev),
     3                rmlexp(iaddr(2,ibox)))
             enddo
C$OMP END PARALLEL DO             
          endif
        endif

        call cpu_time(tt2)
C$    tt2=omp_get_wtime()
        timelev(ilev) = tt2-tt1

        deallocate(wsave)
        deallocate(transvecall)
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(4) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== Step 5 (split loc) ===*',i,0)

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      do ilev = 1,nlevels-1
        dn = 2*(nterms(ilev)+nterms(ilev+1)) + 1
        nsig = next235(dn)
        allocate(wsave(4*nsig+100))
        allocate(sig(nd,nsig))
        allocate(transvecmpmp(nsig,4))
        call zffti(nsig, wsave)
        c1(1) = 0.0d0
        c1(2) = 0.0d0
        do jbox=1,4
           k=2
           if (jbox.le.2) k=1
            c2(1) = 0.25d0*boxsize(ilev)*(-1)**jbox
            c2(2) = 0.25d0*boxsize(ilev)*(-1)**k
            call h2d_mkmpshift(zk,c1,nterms(ilev+1),
     1           c2,nterms(ilev),nsig,wsave,transvecmpmp(1,jbox))
        enddo
        dlam = zk
        dlam = 1/(dlam/(2*pi))                 
        boxlam = boxsize(ilev)/dlam
C$OMP PARALLEL DO DEFAULT(SHARED)        
        do ibox = laddr(1,ilev),laddr(2,ilev)
            call h2d_mptosig(nd,nterms(ilev),nsig,
     1           rmlexp(iaddr(2,ibox)),rmlexp(iaddr(4,ibox)),wsave)
        enddo
C$OMP END PARALLEL DO        
c
       if(zi*boxsize(ilev).lt.zkiupbound) then
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

                  if (boxlam .gt. 16.0d0) then
                   call h2dloclochf(nd,zk,rscales(ilev),
     1                  centers(1,ibox),rmlexp(iaddr(4,ibox)),
     2                  nterms(ilev),nsig,rscales(ilev+1),
     3                  centers(1,jbox),rmlexp(iaddr(2,jbox)),
     4                  nterms(ilev+1),transvecmpmp(1,i),wsave)
                  else
                    call h2dlocloc(nd,zk,rscales(ilev),
     1                  centers(1,ibox),
     1                  rmlexp(iaddr(2,ibox)),nterms(ilev),
     2                  rscales(ilev+1),centers(1,jbox),
     3                  rmlexp(iaddr(2,jbox)),nterms(ilev+1))
                  endif
                  
               enddo
            endif
         enddo
C$OMP END PARALLEL DO        
       endif
       deallocate(wsave)
       deallocate(transvecmpmp)
       deallocate(sig)
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(5) = time2-time1

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      if(ifprint.ge.1)
     $    call prinf('=== Step 6 (mp eval) ===*',i,0)

      do ilev=1,nlevels-1
       if(zi*boxsize(ilev+1).lt.zkiupbound) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,j,i)
C$OMP$PRIVATE(jbox,dlam,boxlam)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            istart = iexpcse(1,ibox)
            iend = iexpcse(2,ibox)
            do j=istart,iend
               do i=1,nlist3s(ibox)
                  jbox = list3(i,ibox)
c                 shift multipole expansion directly to box
c                 for all expansion centers
                  dlam = zk
                  dlam = 1/(dlam/(2*pi))                 
                  boxlam = boxsize(ilev)/dlam
                    call h2dmploc(nd,zk,rscales(ilev+1),
     $                  centers(1,jbox),
     1                  rmlexp(iaddr(1,jbox)),nterms(ilev+1),scjsort(j),
     2                  expcsort(1,j),jsort(1,-ntj,j),ntj)
                  
               enddo
            enddo

c              evalute multipole expansion at all targets
            istart = itargse(1,ibox)
            iend = itargse(2,ibox)
            npts = iend-istart+1

            if(ifpghtarg.eq.1) then
               do i=1,nlist3s(ibox)
                  jbox = list3(i,ibox)
                  
                  call h2dmpevalp(nd,zk,rscales(ilev+1),
     1            centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2            nterms(ilev+1),targetsort(1,istart),npts,
     3            pottarg(1,istart))
               enddo
            endif
            if(ifpghtarg.eq.2) then
               do i=1,nlist3s(ibox)
                  jbox = list3(i,ibox)
                  call h2dmpevalg(nd,zk,rscales(ilev+1),
     1            centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2            nterms(ilev+1),targetsort(1,istart),npts,
     3            pottarg(1,istart),
     4            gradtarg(1,1,istart))
               enddo
            endif
            if(ifpghtarg.eq.3) then
               do i=1,nlist3s(ibox)
                  jbox = list3(i,ibox)

                  call h2dmpevalh(nd,zk,rscales(ilev+1),
     1            centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2            nterms(ilev+1),targetsort(1,istart),npts,
     3            pottarg(1,istart),
     3            gradtarg(1,1,istart),hesstarg(1,1,istart))
               enddo
            endif


c              evalute multipole expansion at all sources
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1
            

            if(ifpgh.eq.1) then
               do i=1,nlist3s(ibox)
                  jbox = list3(i,ibox)
                  call h2dmpevalp(nd,zk,rscales(ilev+1),
     1                centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2                nterms(ilev+1),sourcesort(1,istart),npts,
     3                pot(1,istart))
               enddo
            endif
            if(ifpgh.eq.2) then
               do i=1,nlist3s(ibox)
                  jbox = list3(i,ibox)
                  call h2dmpevalg(nd,zk,rscales(ilev+1),
     1                centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2                nterms(ilev+1),sourcesort(1,istart),npts,
     3                pot(1,istart),grad(1,1,istart))
               enddo
            endif
            if(ifpgh.eq.3) then
               do i=1,nlist3s(ibox)
                  jbox = list3(i,ibox)
                  call h2dmpevalh(nd,zk,rscales(ilev+1),
     1                centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2                nterms(ilev+1),sourcesort(1,istart),npts,
     3                pot(1,istart),grad(1,1,istart),hess(1,1,istart))
               enddo
            endif

         enddo
C$OMP END PARALLEL DO     
       endif
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
       if(zi*boxsize(ilev).lt.zkiupbound) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,i,npts,nchild)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild = itree(iptr(4)+ibox-1)
            if(nchild.eq.0) then
               istart = iexpcse(1,ibox)
               iend = iexpcse(2,ibox)
               do i=istart,iend
                  call h2dlocloc(nd,zk,rscales(ilev),
     $            centers(1,ibox),
     1            rmlexp(iaddr(2,ibox)),nterms(ilev),scjsort(i),
     2            expcsort(1,i),jsort(1,-ntj,i),ntj)
               enddo
c
cc               evaluate local expansion
c                at targets
               istart = itargse(1,ibox)
               iend = itargse(2,ibox)
               npts = iend-istart + 1
               if(ifpghtarg.eq.1) then
                  call h2dtaevalp(nd,zk,rscales(ilev),
     1                centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2                nterms(ilev),targetsort(1,istart),npts,
     3                pottarg(1,istart))
               endif
               if(ifpghtarg.eq.2) then
                  call h2dtaevalg(nd,zk,rscales(ilev),
     1                centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2                nterms(ilev),targetsort(1,istart),npts,
     3                pottarg(1,istart),gradtarg(1,1,istart))
               endif
               if(ifpghtarg.eq.3) then
                  call h2dtaevalh(nd,zk,rscales(ilev),
     1                centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2                nterms(ilev),targetsort(1,istart),npts,
     3                pottarg(1,istart),gradtarg(1,1,istart),
     4                hesstarg(1,1,istart))
               endif

c
cc                evaluate local expansion at sources

               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox)
               npts = iend-istart+1
               if(ifpgh.eq.1) then
                 call h2dtaevalp(nd,zk,rscales(ilev),
     1              centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2              nterms(ilev),sourcesort(1,istart),npts,
     3              pot(1,istart))
               endif
               if(ifpgh.eq.2) then
                 call h2dtaevalg(nd,zk,rscales(ilev),
     1              centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2              nterms(ilev),sourcesort(1,istart),npts,
     3              pot(1,istart),grad(1,1,istart))
               endif
               if(ifpgh.eq.3) then
                 call h2dtaevalh(nd,zk,rscales(ilev),
     1              centers(1,ibox),rmlexp(iaddr(2,ibox)),
     2              nterms(ilev),sourcesort(1,istart),npts,
     3              pot(1,istart),grad(1,1,istart),hess(1,1,istart))
               endif
            endif
         enddo
C$OMP END PARALLEL DO        
       endif
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
      thresh = abs(zk)*boxsize(0)*2.0d0**(-51)


cc      call prin2('thresh=*',thresh,1)
c
cc
      call cpu_time(time1)
C$    time1=omp_get_wtime() 

      if(ifnear.eq.0) goto 1233
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

               call hfmm2dexpc_direct(nd,jstart,jend,istarte,
     1         iende,zk,rscales,nlevels, 
     2         sourcesort,ifcharge,chargesort,ifdipole,dipstrsort,
     3         dipvecsort,expcsort,jsort,scjsort,ntj)

                
               call hfmm2dpart_direct(nd,jstart,jend,istartt,
     1         iendt,zk,sourcesort,ifcharge,chargesort,ifdipole,
     2         dipstrsort,dipvecsort,targetsort,ifpghtarg,pottarg,
     3         gradtarg,hesstarg,thresh)
         
               call hfmm2dpart_direct(nd,jstart,jend,istarts,iends,
     1         zk,sourcesort,ifcharge,chargesort,ifdipole,
     2         dipstrsort,dipvecsort,sourcesort,ifpgh,pot,grad,hess,
     3         thresh)
            enddo   
         enddo
C$OMP END PARALLEL DO         
      enddo
 1233 continue
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
      subroutine hfmm2dexpc_direct(nd,istart,iend,jstart,jend,
     $     zk,rscales,nlevels,source,ifcharge,charge,ifdipole,dipstr,
     $     dipvec,targ,jexps,scj,ntj)
c--------------------------------------------------------------------
c     This subroutine adds the local expansions due to sources
c     istart to iend in the source array at the expansion centers
c     jstart to jend in the target array to the existing local
c     expansions
c
c     INPUT arguments
c-------------------------------------------------------------------
c     nd           in: integer
c                   number of expansions
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
c                  wish to compute the expansions
c 
c     jend         in:Integer
c                  Last index in target array at which we wish
c                  to compute the expansions
c 
c     zk           in: complex *16
c                  Helmholtz parameter
c
c     rscales       in: real*8(0:nlevels)
c                  Scale of expansions formed at all levels
c
c     nlevels      in:Integer
c                  Number of levels in the tree structure
c
c     source       in: real *8(2,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: complex *16
c                  Charge at the source locations
c
c     ifdipole     in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dipstr        in: complex *16(ns)
c                   dip strengths at the source locations
c
c     dipvec       in: complex *16(ns)
c                  dip orientation vectors at the source locations
c
c     targ        in: real *8(2,nexpc)
c                 Expansion center locations
c
c     scj         in: real *8(nexpc)
c                 scaling parameter for expansions
c
c     ntj         in: Integer
c                 Number of terms in expansion
c------------------------------------------------------------
c     OUTPUT
c
c   Updated expansions at the targets
c   jexps : coeffs for local expansions
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j
        integer ifcharge,ifdipole,ier,nd
        complex *16 zk
        real *8 source(2,*)
        real *8 rscales(0:nlevels)
        complex *16 charge(nd,*),dipstr(nd,*)
        real *8 dipvec(nd,2,*)
        real *8 targ(2,*)
        real *8 scj(*)

        integer nlevels,ntj
c
        complex *16 jexps(nd,-ntj:ntj,*)
        
c
        ns = iend - istart + 1
        do j=jstart,jend
           if(ifcharge.eq.1.and.ifdipole.eq.0) then
              call h2dformtac(nd,zk,scj(j),
     1        source(1,istart),charge(1,istart),ns,targ(1,j),
     2        ntj,jexps(1,-ntj,j))
           endif

           if(ifdipole.eq.1.and.ifcharge.eq.0) then
               call h2dformtad(nd,zk,scj(j),
     1         source(1,istart),dipstr(1,istart),dipvec(1,1,istart),
     2         ns,targ(1,j),ntj,jexps(1,-ntj,j))
           endif        
           if(ifdipole.eq.1.and.ifcharge.eq.1) then
               call h2dformtacd(nd,zk,scj(j),
     1         source(1,istart),charge(1,istart),dipstr(1,istart),
     2         dipvec(1,1,istart),
     2         ns,targ(1,j),ntj,jexps(1,-ntj,j))
           endif        
        enddo
c
        return
        end
c------------------------------------------------------------------     
      subroutine hfmm2dpart_direct(nd,istart,iend,jstart,jend,
     $     zk,source,ifcharge,charge,ifdipole,dipstr,dipvec,
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
c     zk           in: complex *16
c                  Complex helmholtz parameter
c
c     source       in: real *8(2,ns)
c                  Source locations
c
c     ifcharge     in: Integer
c                  flag for including expansions due to charges
c                  The expansion due to charges will be included
c                  if ifcharge == 1
c
c     charge       in: complex *16
c                  Charge at the source locations
c
c     ifdipole     in: Integer
c                 flag for including expansions due to dipoles
c                 The expansion due to dipoles will be included
c                 if ifdipole == 1
c
c     dipstr        in: complex *16(ns)
c                 dipole strengths at the source locations
c
c     dipvec        in: complex *16(ns)
c                 dipole orientation vectors at the source locations
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
c                  if |zk*r| < threshold, then interactions are
c                  not included
c
c
c------------------------------------------------------------
c     OUTPUT
c
c   Updated velocity and gradients at the targets
c   pot : potential at the targets
c   grad: gradient at the targets
c   hess: Hessian at the targets
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j,i
        integer ifcharge,ifdipole

        integer nd


        complex *16 zk

        real *8 source(2,*)
        complex *16 charge(nd,*),dipstr(nd,*)
        real *8 dipvec(nd,2,*)

        integer ifpgh
        real *8 targ(2,*),thresh
        
c
        complex *16 pot(nd,*)
        complex *16 grad(nd,2,*)
        complex *16 hess(nd,3,*)

c
        ns = iend - istart + 1
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          if(ifpgh.eq.1) then
             do j=jstart,jend
               call h2d_directcp(nd,zk,source(1,istart),ns,
     1            charge(1,istart),targ(1,j),1,pot(1,j),thresh)
             enddo
          endif

          if(ifpgh.eq.2) then
             do j=jstart,jend
               call h2d_directcg(nd,zk,source(1,istart),ns,
     1            charge(1,istart),targ(1,j),1,pot(1,j),grad(1,1,j),
     2            thresh)
             enddo
          endif
          if(ifpgh.eq.3) then
             do j=jstart,jend
               call h2d_directch(nd,zk,source(1,istart),ns,
     1            charge(1,istart),targ(1,j),1,pot(1,j),grad(1,1,j),
     2            hess(1,1,j),thresh)
             enddo
          endif
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             do j=jstart,jend
               call h2d_directdp(nd,zk,source(1,istart),ns,
     1            dipstr(1,istart),dipvec(1,1,istart),
     2            targ(1,j),1,pot(1,j),thresh)
             enddo
          endif

          if(ifpgh.eq.2) then
             do j=jstart,jend
               call h2d_directdg(nd,zk,source(1,istart),ns,
     1            dipstr(1,istart),dipvec(1,1,istart),
     2            targ(1,j),1,pot(1,j),grad(1,1,j),
     2            thresh)
             enddo
          endif
          if(ifpgh.eq.3) then
             do j=jstart,jend
               call h2d_directdh(nd,zk,source(1,istart),ns,
     1            dipstr(1,istart),dipvec(1,1,istart),targ(1,j),
     2            1,pot(1,j),grad(1,1,j),
     2            hess(1,1,j),thresh)
             enddo
          endif
        endif

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          if(ifpgh.eq.1) then
             do j=jstart,jend
               call h2d_directcdp(nd,zk,source(1,istart),ns,
     1            charge(1,istart),dipstr(1,istart),dipvec(1,1,istart),
     2            targ(1,j),1,pot(1,j),thresh)
             enddo
          endif

          if(ifpgh.eq.2) then
             do j=jstart,jend
               call h2d_directcdg(nd,zk,source(1,istart),ns,
     1            charge(1,istart),dipstr(1,istart),dipvec(1,1,istart),
     2            targ(1,j),1,pot(1,j),grad(1,1,j),
     2            thresh)
             enddo
          endif
          if(ifpgh.eq.3) then
             do j=jstart,jend
               call h2d_directcdh(nd,zk,source(1,istart),ns,
     1            charge(1,istart),dipstr(1,istart),dipvec(1,1,istart),
     2            targ(1,j),1,pot(1,j),grad(1,1,j),
     2            hess(1,1,j),thresh)
             enddo
          endif
        endif


c
        return
        end
c------------------------------------------------------------------    
      subroutine h2dmpalloc(nd,laddr,iaddr,nlevels,lmptot,
     1                          nterms)
c     This subroutine determines the size of the array
c     to be allocated for the multipole expansions
c     iaddr(1,i) points to the starting location of the multipole
c     expansion of box i 
c     iaddr(2,i) points to the local
c     expansion of box i
c     iaddr(3,i) points to the outgoing diag form of box i
c     iaddr(4,i) points to the incoming diag form of box i
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
c     iaddr       out: Integer(4,nboxes)
c                 Points the multipole and local expansions in box i
c 
c     lmptot      out: Integer
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,nterms(0:nlevels),nd,nsig,nt1,nt2,next235
      integer iaddr(4,*), lmptot, laddr(2,0:nlevels)
      integer ibox,i,iptr,istart,nn,itmp
      real *8 ddn,dn
c
      istart = 1
      do i = 0,nlevels
         nn = (2*nterms(i)+1)*2*nd
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
      do i=0,nlevels
         nn = (2*nterms(i)+1)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the local expansion
c
            itmp = ibox - laddr(1,i)
            iaddr(2,ibox) = istart + itmp*nn 
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
      do i=0,nlevels
         dn = 2*(nterms(i)+nterms(i)) + 1
         nn = 2*nd*next235(dn)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the local expansion
c
            itmp = ibox - laddr(1,i)
            iaddr(3,ibox) = istart + itmp*nn 
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
      do i=0,nlevels
         dn = 2*(nterms(i)+nterms(i)) + 1
         nn = 2*nd*next235(dn)
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the local expansion
c
            itmp = ibox - laddr(1,i)
            iaddr(4,ibox) = istart + itmp*nn 
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn
      enddo
c
      lmptot = istart

      return
      end
