cc Copyright (C) 2018-2021: Leslie Greengard, Zydrunas Gimbutas, 
cc  Manas Rachh, Travis Askham
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
      subroutine mbhfmm2d(nd,eps,beta,ns,sources,ifcharge,charge,
     1     ifdipole,dipstr,dipvec,ifquadpole,quadstr,quadvec,ifoctpole,
     2     octstr,octvec,iper,ifpgh,pot,grad,hess,nt,targ,ifpghtarg,
     3     pottarg,gradtarg,hesstarg,ier)
c----------------------------------------------
c     FMM for the modified biharmonic Green's
c     function
c      
c         G(x,y) = (S_beta(r) - S_0(r))/beta^2
c
c     where S_beta(r) = K_0(beta r)/(2*pi) is the fundamental solution
c     of the Yukawa equation and S_0(r) = log(r)/(2*pi) is the
c     fundamental solution of the Laplace equation, with
c     r = sqrt( (x(1)-y(1))**2 + (x(2)-y(2))**2 )
c
c     The expansions represent the sum
c
c     u(x) = sum_j charge(j)*G(x,y(j))
c     + dipstr(j)*(G_{y1}(x,y(j))*dipvec(1,j) + G_{y2}*dipvec(2,j))
c     + quadstr(j)*(G_{y1y1}(x,y(j))*quadvec(1,j)
c             + G_{y1y2}*quadvec(2,j) + G_{y2y2}*quadvec(3,j))
c     + octstr(j)*(G_{y1y1y1}(x,y(j))*octvec(1,j)
c             + G_{y1y1y2}*octvec(2,j) + G_{y1y2y2}*octvec(3,j)
c             + G_{y2y2y2}*octvec(4,j))
c
c     INPUT PARAMETERS:
c      
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
c                   dipole interactions included if ifdipole =1
c                   not included otherwise
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   ifquadpole      : flag for including quadrupole interactions
c                   quadrupole interactions included if ifquadpole =1
c                   not included otherwise
c   quadstr(nd,ns)    : quadrupole strengths
c   quadvec(nd,3,ns)  : quadrupole orientation vectors
c   ifoctpole      : flag for including octopole interactions
c                   octopole interactions included if ifoctpole =1
c                   not included otherwise
c   octstr(nd,ns)    : octopole strengths
c   octvec(nd,4,ns)  : octopole orientation vectors
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
      integer nd, ifcharge, ifdipole, ifquadpole, ifoctpole
      integer iper
      real *8 eps
      real *8 beta
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,*)
      real *8 charge(nd,*),dipstr(nd,*)
      real *8 quadvec(nd,3,*),quadstr(nd,*)
      real *8 octvec(nd,4,*),octstr(nd,*)      

      real *8 pot(nd,*),grad(nd,2,*),hess(nd,3,*)
      real *8 pottarg(nd,*),gradtarg(nd,2,*),hesstarg(nd,3,*)

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
      integer idivflag,nlevels,nboxes,ndiv
      integer ltree

      real *8, allocatable :: radsrc(:)
      real *8 radexp

c
cc     sorted arrays
c
      integer, allocatable :: isrc(:),isrcse(:,:)
      integer, allocatable :: itarg(:),itargse(:,:)
      real *8, allocatable :: sourcesort(:,:)
      real *8, allocatable :: targsort(:,:)
      
      complex *16, allocatable :: mbhmps(:,:,:), ymps(:,:,:)
      complex *16, allocatable :: mbhmpssort(:,:,:), ympssort(:,:,:)   

      real *8, allocatable :: potsort(:,:),gradsort(:,:,:),
     1     hesssort(:,:,:)
      real *8, allocatable :: pottargsort(:,:),gradtargsort(:,:,:),
     1     hesstargsort(:,:,:)

c
cc     additional fmm variables

      integer lmptot, ldc
      real *8, allocatable :: rscales(:)
      integer, allocatable :: nterms(:),iaddr(:,:)
      real *8, allocatable :: rmlexp(:),carray(:,:)

c
cc      temporary variables
c
      integer i,ilev,lmptmp,nmax,idim,j,l
      integer ifpgh,ifpghtarg,ifprint,iert
      integer ntermsmps,ndtmp,nterms1
      real *8 time1,time2,pi,done
      real *8 omp_get_wtime
      complex *16 zk, eye
      data eye /(0.0d0,1.0d0)/

      done = 1
      pi = atan(done)*4.0d0

      zk = beta*eye

      nexpc = 0

      nlevels = 0
      nboxes = 0
      idivflag =0
      ndiv = 20
      ltree = 0
      nlmin = 0
      nlmax = 51
      ifunif = 0
      iper = 0

      ifprint = 1

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
      allocate(itarg(nt),itargse(2,nboxes))

      call pts_tree_sort(ns,sources,itree,ltree,nboxes,nlevels,iptr,
     1   tcenters,isrc,isrcse)

      call pts_tree_sort(nt,targ,itree,ltree,nboxes,nlevels,iptr,
     1   tcenters,itarg,itargse)
      allocate(sourcesort(2,ns))
      allocate(targsort(2,nt))

      ntermsmps=0
      if (ifdipole .eq. 1) ntermsmps = 1
      if (ifquadpole .eq. 1) ntermsmps = 2
      if (ifoctpole .eq. 1) ntermsmps = 3

      call prinf('ntermsmp *',ntermsmps,1)

      allocate(mbhmps(nd,0:ntermsmps,ns),
     1     ymps(nd,0:ntermsmps,ns))

      allocate(mbhmpssort(nd,0:ntermsmps,ns),
     1     ympssort(nd,0:ntermsmps,ns))

      do i = 1,ns
         do j = 0,ntermsmps
            do l = 1,nd
               mbhmps(l,j,i)=0
               ymps(l,j,i)=0
            enddo
         enddo
      enddo
      
      call mbh2dconvtomp_vec(nd,beta,ns,ifcharge,charge,
     1     ifdipole,dipstr,dipvec,ifquadpole,quadstr,quadvec,
     2     ifoctpole,octstr,octvec,ntermsmps,mbhmps,ymps)

      

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
        call l2dterms(eps,nterms1,ier)
        nterms(i) = nterms1
        if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo

      if(ifprint.eq.1) call prinf('nmax=*',nmax,1)
      if(ifprint.eq.1) call prinf('nterms=*',nterms,nlevels+1)

      ldc = nmax+5
      allocate(carray(0:ldc,0:ldc))
      call mbh2d_init_carray(carray,ldc)

      
c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(4,nboxes).
c
c       ... allocate iaddr and temporary arrays
c

      allocate(iaddr(4,nboxes))


c     reorder sources
c
      call dreorderf(2,ns,sources,sourcesort,isrc)
      ndtmp = nd*2*(ntermsmps+1)
      call dreorderf(ndtmp,ns,mbhmps,mbhmpssort,isrc)
      call dreorderf(ndtmp,ns,ymps,ympssort,isrc)      

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
      call mbh2dmpalloc(nd,itree(iptr(1)),iaddr,nlevels,lmptot,
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
      call mbhfmm2dmain(nd,eps,
     $   beta,ns,sourcesort,
     $   ntermsmps,mbhmpssort,ympssort,
     $   nt,targsort,
     $   iaddr,rmlexp,carray,ldc,
     $   itree,ltree,iptr,ndiv,nlevels,
     $   nboxes,iper,boxsize,rscales,tcenters,itree(iptr(1)),
     $   isrcse,itargse,nterms,
     $   ifpgh,potsort,gradsort,hesssort,
     $   ifpghtarg,pottargsort,gradtargsort,
     $   hesstargsort,ier)
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)
      if(ier.ne.0) return


c
cc      resort the output arrays in input order
c

      if(ifpgh.eq.1) then
        call dreorderi(nd,ns,potsort,pot,isrc)
      endif

      if(ifpgh.eq.2) then
        call dreorderi(nd,ns,potsort,pot,isrc)
        call dreorderi(2*nd,ns,gradsort,grad,isrc)
      endif

      if(ifpgh.eq.3) then
        call dreorderi(nd,ns,potsort,pot,isrc)
        call dreorderi(2*nd,ns,gradsort,grad,isrc)
        call dreorderi(3*nd,ns,hesssort,hess,isrc)
      endif

cc      call prini(6,13)
cc      call prin2('eps = *', eps, 1)
c      call prin2('after mbhfmm2dmain, potsort = *', potsort,
c     1     2*nt)
cc      stop
      
      if(ifpghtarg.eq.1) then
        call dreorderi(nd,nt,pottargsort,pottarg,itarg)
      endif

      if(ifpghtarg.eq.2) then
        call dreorderi(nd,nt,pottargsort,pottarg,itarg)
        call dreorderi(2*nd,nt,gradtargsort,gradtarg,itarg)
      endif

      if(ifpghtarg.eq.3) then
        call dreorderi(nd,nt,pottargsort,pottarg,itarg)
        call dreorderi(2*nd,nt,gradtargsort,gradtarg,itarg)
        call dreorderi(3*nd,nt,hesstargsort,hesstarg,itarg)
      endif


      return
      end

      subroutine mbhfmm2dmain(nd,eps,
     $     beta,nsource,source,
     $     ntermsmps,mbhmps,ymps,
     $     ntarget,target,
     $     iaddr,rmlexp,carray,ldc,
     $     itree,ltree,iptr,ndiv,nlevels, 
     $     nboxes,iper,boxsize,rscales,centers,laddr,
     $     isrcse,itargse,nterms,
     $     ifpgh,pot,grad,hess,
     $     ifpghtarg,pottarg,gradtarg,hesstarg,
     $     ier)
c     2D FMM for the modified biharmonic Green's
c     function
c      
c         G(x,y) = (S_beta(r) - S_0(r))/beta^2
c
c     where S_beta(r) = K_0(beta r)/(2*pi) is the fundamental solution
c     of the Yukawa equation and S_0(r) = log(r)/(2*pi) is the
c     fundamental solution of the Laplace equation, with
c     r = sqrt( (x(1)-y(1))**2 + (x(2)-y(2))**2 )
c
c     This FMM takes charges to be multipolar expansions
c     of a given length.
c
c   All the source/target related quantities
c   are assumed to be tree-sorted
c
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
c
c   nd:   number of charge densities
c
c   eps:  FMM precision requested
c
c   beta: real *8, modified biharmonic parameter
c
c     nsource:     integer:  number of sources
c     source: real *8 (2,ns):  source locations
c     ntermsmps:   integer: order of multipolar sources
c     mbhmps:   complex *16(nd,0:ntermsmps,nsource): difference function
c                   type expansion at each source
c     ymps:     complex *16(nd,0:ntermsmps,nsource): Yukawa (modified
c     bessel K) type expansion at each source
c     ntarget: integer:  number of targets
c     target: real *8 (2,ntarget):  target locations
c     iaddr: (4,nboxes): pointer in rmlexp where multipole
c                      and local expansions for each
c                      box is stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      difference typ multipole expansion of ibox
c                      iaddr(2,ibox) is the
c                      starting index in rmlexp for the 
c                      Yukawa (Bessel K) multipole expansion of ibox
c                      iaddr(3,ibox) is the
c                      starting index in rmlexp for the 
c                      difference type local expansion of ibox
c                      and iaddr(4,ibox) is the
c                      starting index in rmlexp
c                      for the Laplace (power series)
c                      local expansion of ibox
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

      real *8 zi, beta

      integer nsource,ntarget,ier
      integer ndiv,nlevels,ntermsmps

      integer ifpgh,ifpghtarg
      real *8 eps
      integer iper

      real *8 source(2,nsource)

      complex *16 :: mbhmps(nd,0:ntermsmps,nsource)
      complex *16 :: ymps(nd,0:ntermsmps,nsource)     
      
      real *8 target(2,ntarget)

      real *8 pot(nd,*)
      real *8 grad(nd,2,*)
      real *8 hess(nd,3,*)

      real *8 pottarg(nd,*)
      real *8 gradtarg(nd,2,*)
      real *8 hesstarg(nd,3,*)

      integer iaddr(4,nboxes),lmptmp
      real *8 rmlexp(*)
       
      real *8 timeinfo(10)
      real *8 timelev(0:200)
      real *8 centers(2,*)

      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer iptr(8),ltree
      integer isrcse(2,nboxes),itargse(2,nboxes)
      integer itree(*)
      integer nboxes, ldc
      real *8 rscales(0:nlevels),boxsize(0:nlevels)
      real *8 carray(0:ldc,0:ldc)

      real *8 zkiupbound,thresh

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
      real *8 pottmp,gradtmp(2),hesstmp(3)
      
      double precision dlam, pi, boxlam
      
c     ifprint is an internal information printing flag. 
c     Suppressed if ifprint=0.
c     Prints timing breakdown and other things if ifprint=1.
c     Prints timing breakdown, list information, and other things if ifprint=2.
c      
      ifprint=1

      pi = 4*atan(1.0d0)
c     
c     c
c     upper limit for zk along imaginary axis
      zkiupbound = 40.0d0
      zi = beta

      do i=0,nlevels
         timelev(i) = 0
      enddo

c
c        compute list info
c     
      call computemnlists(nlevels,nboxes,itree,ltree,iptr,centers,
     1     boxsize,iper,mnlist1,mnlist2,mnlist3,mnlist4)
      allocate(nlist1s(nboxes),list1(mnlist1,nboxes))
      allocate(nlist2s(nboxes),list2(mnlist2,nboxes))
      allocate(nlist3s(nboxes),list3(mnlist3,nboxes))
      allocate(nlist4s(nboxes),list4(mnlist4,nboxes))
      
      call computelists(nlevels,nboxes,itree,ltree,iptr,centers,
     1     boxsize,iper,mnlist1,nlist1s,list1,mnlist2,nlist2s,list2,
     2     mnlist3,nlist3s,list3,mnlist4,nlist4s,list4)


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
            call mbh2dmpzero(nd,rmlexp(iaddr(1,ibox)),nterms(ilev))
            call mbh2dmpzero(nd,rmlexp(iaddr(2,ibox)),nterms(ilev))
            call mbh2dmpzero(nd,rmlexp(iaddr(3,ibox)),nterms(ilev))
            call mbh2dmpzero(nd,rmlexp(iaddr(4,ibox)),nterms(ilev))
         enddo
C$OMP END PARALLEL DO         
      enddo
             
c
      if(ifprint .ge. 1) 
     $     call prinf('=== STEP 1 (form mp) ====*',i,0)
      call cpu_time(time1)
C$        time1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions

      do ilev = 2,nlevels
         if(zi*boxsize(ilev).lt.zkiupbound) then
C     

C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox=laddr(1,ilev),laddr(2,ilev)
               nchild = itree(iptr(4)+ibox-1)
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox)
               npts = iend-istart+1
c     Check if current box is a leaf box            
               if(nchild.eq.0.and.npts.gt.0) then
                  call mbh2dformmpmp_vec(nd,beta,rscales(ilev),
     1                 source(1,istart),npts,mbhmps(1,0,istart),
     1                 ymps(1,0,istart),ntermsmps,centers(1,ibox),
     1                 nterms(ilev),rmlexp(iaddr(1,ibox)),
     1                 rmlexp(iaddr(2,ibox)))
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
       if(zi*boxsize(ilev).lt.zkiupbound) then
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
                   
                   call mbh2dformtamp_vec(nd,beta,rscales(ilev),
     1                  source(1,istart),npts,mbhmps(1,0,istart),
     2                  ymps(1,0,istart),ntermsmps,centers(1,ibox),
     3                  nterms(ilev),rmlexp(iaddr(3,ibox)),
     4                  rmlexp(iaddr(4,ibox)))

                enddo
             endif
          enddo
C$OMP END PARALLEL DO        

       endif
      enddo
      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(2)=time2-time1

      if(ifprint .ge. 1)
     $      call prinf('=== STEP 3 (merge mp) ====*',i,0)
      call cpu_time(time1)
C$    time1=omp_get_wtime()
c
      do ilev=nlevels-1,1,-1

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
                     call mbh2dmpmp_vec(nd,beta,rscales(ilev+1),
     1                    centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2                    rmlexp(iaddr(2,jbox)),nterms(ilev+1),
     3                    rscales(ilev),centers(1,ibox),
     4                    rmlexp(iaddr(1,ibox)),rmlexp(iaddr(2,ibox)),
     5                    nterms(ilev))
                  endif
               enddo
            enddo
C$OMP END PARALLEL DO    
         endif
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

       if(zi*boxsize(ilev).lt.zkiupbound) then
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

            if(ifpgh.gt.0) then
              istart = isrcse(1,ibox)
              iend = isrcse(2,ibox)
              npts = npts + iend-istart+1
            endif

            if(npts.gt.0) then
               do i=1,nlist2s(ibox)
                  jbox = list2(i,ibox) 

                  call mbh2dmploc_vec(nd,beta,rscales(ilev),
     1                 centers(1,jbox),rmlexp(iaddr(1,jbox)),
     2                 rmlexp(iaddr(2,jbox)),nterms(ilev),
     3                 rscales(ilev),centers(1,ibox),
     4                 rmlexp(iaddr(3,ibox)),
     5                 rmlexp(iaddr(4,ibox)),nterms(ilev))
                  
               enddo
            endif
         enddo
C$OMP END PARALLEL DO        
       endif

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
         if(zi*boxsize(ilev).lt.zkiupbound) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,nchild,istart,iend,npts)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
               nchild = itree(iptr(4)+ibox-1)
               
               npts = 0

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

                     call mbh2dlocloc_vec(nd,beta,rscales(ilev),
     1                    centers(1,ibox),rmlexp(iaddr(3,ibox)),
     2                    rmlexp(iaddr(4,ibox)),nterms(ilev),
     3                    rscales(ilev+1),centers(1,jbox),
     4                    rmlexp(iaddr(3,jbox)),rmlexp(iaddr(4,jbox)),
     5                    nterms(ilev+1),carray,ldc)
                     
                  enddo
               endif
            enddo
C$OMP END PARALLEL DO        
         endif
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
       if(zi*boxsize(ilev+1).lt.zkiupbound) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,j,i)
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
                  
                  call mbh2dmpevalp_vec(nd,beta,rscales(ilev+1),
     1                 centers(1,jbox),rmlexp(iaddr(1,jbox)),
     1                 rmlexp(iaddr(2,jbox)),
     2                 nterms(ilev+1),target(1,istart),npts,
     3                 pottarg(1,istart))
               enddo
            endif
            if(ifpghtarg.eq.2) then
               do i=1,nlist3s(ibox)
                  jbox = list3(i,ibox)
                  call mbh2dmpevalg_vec(nd,beta,rscales(ilev+1),
     1                 centers(1,jbox),rmlexp(iaddr(1,jbox)),
     1                 rmlexp(iaddr(2,jbox)),
     2                 nterms(ilev+1),target(1,istart),npts,
     3                 pottarg(1,istart),
     4                 gradtarg(1,1,istart))
               enddo
            endif
            if(ifpghtarg.eq.3) then
               do i=1,nlist3s(ibox)
                  jbox = list3(i,ibox)

                  call mbh2dmpevalh_vec(nd,beta,rscales(ilev+1),
     1                 centers(1,jbox),rmlexp(iaddr(1,jbox)),
     1                 rmlexp(iaddr(2,jbox)),
     2                 nterms(ilev+1),target(1,istart),npts,
     3                 pottarg(1,istart),
     3                 gradtarg(1,1,istart),hesstarg(1,1,istart))
               enddo
            endif


c              evalute multipole expansion at all sources
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1
            

            if(ifpgh.eq.1) then
               do i=1,nlist3s(ibox)
                  jbox = list3(i,ibox)
                  call mbh2dmpevalp_vec(nd,beta,rscales(ilev+1),
     1                 centers(1,jbox),rmlexp(iaddr(1,jbox)),
     1                 rmlexp(iaddr(2,jbox)),
     2                nterms(ilev+1),source(1,istart),npts,
     3                pot(1,istart))
               enddo
            endif
            if(ifpgh.eq.2) then
               do i=1,nlist3s(ibox)
                  jbox = list3(i,ibox)
                  call mbh2dmpevalg_vec(nd,beta,rscales(ilev+1),
     1                 centers(1,jbox),rmlexp(iaddr(1,jbox)),
     1                 rmlexp(iaddr(2,jbox)),
     2                 nterms(ilev+1),source(1,istart),npts,
     3                 pot(1,istart),grad(1,1,istart))
               enddo
            endif
            if(ifpgh.eq.3) then
               do i=1,nlist3s(ibox)
                  jbox = list3(i,ibox)
                  call mbh2dmpevalh_vec(nd,beta,rscales(ilev+1),
     1                 centers(1,jbox),rmlexp(iaddr(1,jbox)),
     1                 rmlexp(iaddr(2,jbox)),
     2                 nterms(ilev+1),source(1,istart),npts,
     3                 pot(1,istart),grad(1,1,istart),hess(1,1,istart))
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
C$OMP$PRIVATE(ibox,istart,iend,i,npts)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox = laddr(1,ilev),laddr(2,ilev)
            nchild = itree(iptr(4)+ibox-1)
            if(nchild.eq.0) then
c     c               evaluate local expansion
c                at targets
               istart = itargse(1,ibox)
               iend = itargse(2,ibox)
               npts = iend-istart + 1
               if(ifpghtarg.eq.1) then
                  call mbh2dtaevalp_vec(nd,beta,rscales(ilev),
     1                 centers(1,ibox),rmlexp(iaddr(3,ibox)),
     1                 rmlexp(iaddr(4,ibox)),
     2                 nterms(ilev),target(1,istart),npts,
     3                 pottarg(1,istart))
               endif
               if(ifpghtarg.eq.2) then
                  call mbh2dtaevalg_vec(nd,beta,rscales(ilev),
     1                 centers(1,ibox),rmlexp(iaddr(3,ibox)),
     2                 rmlexp(iaddr(4,ibox)),
     2                nterms(ilev),target(1,istart),npts,
     3                pottarg(1,istart),gradtarg(1,1,istart))
               endif
               if(ifpghtarg.eq.3) then
                  call mbh2dtaevalh_vec(nd,beta,rscales(ilev),
     1                centers(1,ibox),rmlexp(iaddr(3,ibox)),
     2                 rmlexp(iaddr(4,ibox)),
     2                 nterms(ilev),target(1,istart),npts,
     3                pottarg(1,istart),gradtarg(1,1,istart),
     4                hesstarg(1,1,istart))
               endif

c
cc                evaluate local expansion at sources

               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox)
               npts = iend-istart+1
               if(ifpgh.eq.1) then
                 call mbh2dtaevalp_vec(nd,beta,rscales(ilev),
     1                 centers(1,ibox),rmlexp(iaddr(3,ibox)),
     1                 rmlexp(iaddr(4,ibox)),
     2                 nterms(ilev),source(1,istart),npts,
     3                 pot(1,istart))
              endif
              if(ifpgh.eq.2) then
                 call mbh2dtaevalg_vec(nd,beta,rscales(ilev),
     1                centers(1,ibox),rmlexp(iaddr(3,ibox)),
     1                rmlexp(iaddr(4,ibox)),
     2                nterms(ilev),source(1,istart),npts,
     3                pot(1,istart),grad(1,1,istart))
              endif
              if(ifpgh.eq.3) then
                 call mbh2dtaevalh_vec(nd,beta,rscales(ilev),
     1                centers(1,ibox),rmlexp(iaddr(3,ibox)),
     1                rmlexp(iaddr(4,ibox)),
     2                nterms(ilev),source(1,istart),npts,
     3                pot(1,istart),grad(1,1,istart),hess(1,1,istart))
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
      thresh = boxsize(0)*1.0d-16
      call prin2('thresh=*',thresh,1)
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

            istarts = isrcse(1,ibox)
            iends = isrcse(2,ibox)

            do i =1,nlist1s(ibox)
               jbox = list1(i,ibox) 

               jstart = isrcse(1,jbox)
               jend = isrcse(2,jbox)

               call mbhfmm2dmps_direct_vec(nd,jstart,jend,istartt,
     1              iendt,beta,source,mbhmps,ymps,ntermsmps,
     2              target,ifpghtarg,pottarg,
     3              gradtarg,hesstarg,thresh)
         
               call mbhfmm2dmps_direct_vec(nd,jstart,jend,istarts,
     1              iends,beta,source,mbhmps,ymps,ntermsmps,
     2              source,ifpgh,pot,
     3              grad,hess,thresh)
         
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
c------------------------------------------------------------------     
      subroutine mbhfmm2dmps_direct_vec(nd,istart,iend,jstart,jend,
     $     beta,source,mbhmps,ymps,ntermsmps,
     $     targ,ifpgh,pot,grad,hess,thresh)
c--------------------------------------------------------------------
c     This subroutine adds the contribuition due to multipolar sources
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
c     beta         in: real *8
c                  modified biharmonic parameter
c
c     source       in: real *8(2,ns)
c                  Source locations
c     
c     mbhmps       in: complex *16 (nd,0:ntermsmps,ns)
c                  difference kernel type expansion at each source 
c
c     ymps         in: complex *16 (nd,0:ntermsmps,ns)
c                  Yukawa (Bessel K) type expansion at each source 
c
c     ntermsmps    in: integer, order of the multipolar charges at
c                  each source
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
c   grad: gradient at the targets
c   hess: Hessian at the targets
c-------------------------------------------------------               
        implicit none
c
        integer istart,iend,jstart,jend,ns,j,i
        integer ifcharge,ifdipole

        integer nd,ntermsmps

        real *8 beta

        real *8 source(2,*)
        complex *16 mbhmps(nd,0:ntermsmps,*),ymps(nd,0:ntermsmps,*)

        integer ifpgh
        real *8 targ(2,*),thresh
        
c
        real *8  pot(nd,*)
        real *8  grad(nd,2,*)
        real *8  hess(nd,3,*)

        integer nt

c
        ns = iend - istart + 1
        nt = jend - jstart + 1
        
        if(ifpgh.eq.1) then
           call mbh2d_directmpsp_vec(nd,beta,source(1,istart),ns,
     1          mbhmps(1,0,istart),ymps(1,0,istart),ntermsmps,
     1          targ(1,jstart),nt,pot(1,jstart),thresh)
        endif

        if(ifpgh.eq.2) then
           call mbh2d_directmpsg_vec(nd,beta,source(1,istart),ns,
     1          mbhmps(1,0,istart),ymps(1,0,istart),ntermsmps,
     1          targ(1,jstart),nt, pot(1,jstart),grad(1,1,jstart),
     1          thresh)
        endif
        if(ifpgh.eq.3) then
           call mbh2d_directmpsh_vec(nd,beta,source(1,istart),ns,
     1          mbhmps(1,0,istart),ymps(1,0,istart),ntermsmps,
     1          targ(1,jstart),nt,pot(1,jstart),grad(1,1,jstart),
     1          hess(1,1,jstart),thresh)
        endif

c
        return
        end
c------------------------------------------------------------------    
      subroutine mbh2dmpalloc(nd,laddr,iaddr,nlevels,lmptot,
     1     nterms)
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
c     iaddr       out: Integer(4,nboxes)
c     Points to the multipole and local expansions in box i
c     iaddr(1,ibox) - difference type multipole exp
c     iaddr(2,ibox) - Yukawa (Bessel K) type multipole exp
c     iaddr(3,ibox) - difference type local exp
c     iaddr(4,ibox) - Laplace (power series) type multipole exp
c 
c     lmptot      out: Integer
c                 Total length of expansions array required
c------------------------------------------------------------------

      implicit none
      integer nlevels,nterms(0:nlevels),nd,nsig,nt1,nt2,next235
      integer iaddr(4,*), lmptot, laddr(2,0:nlevels)
      integer ibox,i,iptr,istart,nn,itmp
      real *8 ddn
c
      istart = 1
      do i = 0,nlevels

         nn = (nterms(i)+1)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
c     Allocate memory for the multipole expansion         
c
           itmp = ibox - laddr(1,i)
           iaddr(1,ibox) = istart + itmp*nn*2
           iaddr(2,ibox) = istart + itmp*nn*2 + nn
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn*2
       enddo
c
c            Allocate memory for the local expansion
c
       do i=0,nlevels
         nn = (nterms(i)+1)*2*nd
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,itmp)
         do ibox = laddr(1,i),laddr(2,i)
             itmp = ibox - laddr(1,i)
             iaddr(3,ibox) = istart + itmp*nn*2
             iaddr(4,ibox) = istart + itmp*nn*2 + nn              
         enddo
C$OMP END PARALLEL DO         
         istart = istart + (laddr(2,i)-laddr(1,i)+1)*nn*2
      enddo
      lmptot = istart

      return
      end
c----------------------------------------------------------------     
c
C***********************************************************************
      subroutine mbh2dmpzero(nd,mpole,nterms)
      implicit none
C***********************************************************************
c
c     This subroutine sets a multipole expansion of the given
c     order to zero.
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :   number of expansions
c     nterms :   order of multipole expansions
C---------------------------------------------------------------------
c     OUTPUT:
c
c     mpole  :   coeffs for the expansion set to zero.
C---------------------------------------------------------------------
      integer n,nterms,nd,idim
      complex *16 mpole(nd,nterms+1)
c
      do n=1,(nterms+1)
        do idim=1,nd
          mpole(idim,n)=0.0d0
        enddo
      enddo
      return
      end
c
c
