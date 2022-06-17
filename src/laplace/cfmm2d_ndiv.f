cc Copyright (C) 2018-2019: Le=",nlevelsslie Greengard, Zydrunas Gimbutas, 
cc and Manas Rachh
c     c Contact: greengard@cims.nyu.edu
cc
cc    convert to Cauchy FMM - Travis Askham 2021/07/07
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

      subroutine cfmm2d_ndiv(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ndiv,idivflag,ifnear,timeinfo,ier)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
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
c   ndiv          : subdivision criterion, number of points of box
c   idivflag      : subdivision criterion, whether to use sources, targets
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : potential at the source locations
c   grad(nd,*)    : gradients (d/dz) at the source locations
c   hess(nd,*)    : hessian (d^2/dz^2) at the source locations
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,*): gradient (d/dz) at the target locations
c   hesstarg(nd,*): hessian (d^2/dz^2) at the target locations
c   timeinfo        : time distribution
c   ier             : error code
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
      complex *16 charge(nd,*),dipstr(nd,*)

      complex *16 pot(nd,*),grad(nd,*),hess(nd,*)
      complex *16 pottarg(nd,*),gradtarg(nd,*),hesstarg(nd,*)

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

      integer ifnear
      real *8 timeinfo(8)

c
cc     sorted arrays
c
      integer, allocatable :: isrc(:),isrcse(:,:)
      integer, allocatable :: itarg(:),itargse(:,:),iexpcse(:,:)
      real *8, allocatable :: sourcesort(:,:)
      real *8, allocatable :: targsort(:,:)
      complex *16, allocatable :: chargesort(:,:),dipstrsort(:,:)
      complex *16, allocatable :: potsort(:,:),gradsort(:,:),
     1                             hesssort(:,:)
      complex *16, allocatable :: pottargsort(:,:),gradtargsort(:,:),
     1                              hesstargsort(:,:)

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
c    Need to fix ndiv in Laplace FMM
c   
c

      nlevels = 0
      nboxes = 0


      ltree = 0
      nlmin = 0
      nlmax = 51
      ifunif = 0
      iper = 0

      ifprint = 0

c
c  initialize timeinfo
c
      do i=1,8
        timeinfo(i) = 0
      enddo

c
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
        allocate(chargesort(nd,ns),dipstrsort(nd,1))
      endif
      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        allocate(chargesort(nd,1),dipstrsort(nd,ns))
      endif
      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        allocate(chargesort(nd,ns),dipstrsort(nd,ns))
      endif

      if(ifpgh.eq.1) then
        allocate(potsort(nd,ns),gradsort(nd,1),hesssort(nd,1))
      else if(ifpgh.eq.2) then
        allocate(potsort(nd,ns),gradsort(nd,ns),hesssort(nd,1))
      else if(ifpgh.eq.3) then
        allocate(potsort(nd,ns),gradsort(nd,ns),hesssort(nd,ns))
      else
        allocate(potsort(nd,1),gradsort(nd,1),hesssort(nd,1))
      endif

      
      if(ifpghtarg.eq.1) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,1),
     1     hesstargsort(nd,1))
      else if(ifpghtarg.eq.2) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,nt),
     1      hesstargsort(nd,1))
      else if(ifpghtarg.eq.3) then
        allocate(pottargsort(nd,nt),gradtargsort(nd,nt),
     1     hesstargsort(nd,nt))
      else
        allocate(pottargsort(nd,1),gradtargsort(nd,1),
     1     hesstargsort(nd,1))
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
            gradsort(idim,i) = 0
          enddo
        enddo
      endif

      if(ifpgh.eq.3) then
        do i=1,ns
          do idim=1,nd
            potsort(idim,i) = 0
            gradsort(idim,i) = 0
            hesssort(idim,i) = 0
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
            gradtargsort(idim,i) = 0
            gradtargsort(idim,i) = 0
          enddo
        enddo
      endif

      if(ifpghtarg.eq.3) then
        do i=1,nt
          do idim=1,nd
            pottargsort(idim,i) = 0
            gradtargsort(idim,i) = 0
            hesstargsort(idim,i) = 0
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
        rscales(i) = boxsize(i)

        call l2dterms(eps,nterms(i),ier)
        nterms(i) = nterms(i) 
        if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo

      if(ifprint.eq.1) call prinf('nmax=*',nmax,1)
      if(ifprint.eq.1) call prinf('nterms=*',nterms,nlevels+1)

c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(3,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c

      allocate(iaddr(2,nboxes))

      lmptmp = (nmax+1)*nd
      allocate(mptemp(lmptmp))

c     reorder sources
c
      call dreorderf(2,ns,sources,sourcesort,isrc)
      if(ifcharge.eq.1) 
     1     call dreorderf(2*nd,ns,charge,chargesort,isrc)
      if(ifdipole.eq.1) then
         call dreorderf(2*nd,ns,dipstr,dipstrsort,isrc)
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
      call l2dmpalloc(nd,itree,iaddr,nlevels,lmptot,
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
      call cfmm2dmain(nd,eps,
     $   ns,sourcesort,
     $   ifcharge,chargesort,
     $   ifdipole,dipstrsort,
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


c
cc      resort the output arrays in input order
c

      if(ifpgh.eq.1) then
        call dreorderi(2*nd,ns,potsort,pot,isrc)
      endif

      if(ifpgh.eq.2) then
        call dreorderi(2*nd,ns,potsort,pot,isrc)
        call dreorderi(2*nd,ns,gradsort,grad,isrc)
      endif

      if(ifpgh.eq.3) then
        call dreorderi(2*nd,ns,potsort,pot,isrc)
        call dreorderi(2*nd,ns,gradsort,grad,isrc)
        call dreorderi(2*nd,ns,hesssort,hess,isrc)
      endif

cc      call prini(6,13)
cc      call prin2('eps = *', eps, 1)
cc      call prin2('after lfmm2dmain, pottargsort = *', pottargsort, 30)
cc      stop
      
      if(ifpghtarg.eq.1) then
        call dreorderi(2*nd,nt,pottargsort,pottarg,itarg)
      endif

      if(ifpghtarg.eq.2) then
        call dreorderi(2*nd,nt,pottargsort,pottarg,itarg)
        call dreorderi(2*nd,nt,gradtargsort,gradtarg,itarg)
      endif

      if(ifpghtarg.eq.3) then
        call dreorderi(2*nd,nt,pottargsort,pottarg,itarg)
        call dreorderi(2*nd,nt,gradtargsort,gradtarg,itarg)
        call dreorderi(2*nd,nt,hesstargsort,hesstarg,itarg)
      endif


      return
      end
c
c
c
