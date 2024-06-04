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
      subroutine hfmm2d_mps(nd,eps,zk,nmpole,cmpole,rmpole,
     1  mterms,ntot,mpole,impole,local,ier)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd               : number of expansions
c   eps              : FMM precision requested
c   zk               : Helmholtz parameter
c   nmpole           : number of multipole expansions 
c   cmpole(2,nmpole) : center of mpole expansions 
c   rmpole(nmpole)   : scaling factor for multipole expansion 
c   mterms(nmpole)   : number of terms in multipole and local expansions 
c   ntot             : total length of multipole expansion array,
c                        should be \sum_{j} (2*mterms(j) + 1) in size*nd
c   mpole(ntot)      : coefficients of multipole expansion
c   impole(nmpole)   : indexing array for mpole, the ith expansion is at
c                       location mpole(1,impole(i))
c
c   OUTPUT PARAMETERS
c   local(nd,ntot)  : local expansion at each center, due to all 
c                     incoming multipole expansions (self is ignored).
c                     the orders are same as for the incoming multipole
c   ier             : error code
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      complex *16 zk
      integer nmpole,ier,mterms(nmpole),impole(nmpole),ntot
      real *8 cmpole(2,nmpole),rmpole(nmpole)
      double complex :: mpole(ntot),local(ntot)

c
cc      Tree variables
c
      integer, allocatable :: itree(:)
      real *8 targ(2)
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
      integer, allocatable :: isrc(:),isrcse(:,:)
      integer iper

c
cc     sorted arrays
c
      integer :: lmpole,mt,ilen
      integer, allocatable :: mtermssort(:),impolesort(:)
      real *8, allocatable :: cmpolesort(:,:),rmpolesort(:)
      double complex, allocatable :: mpolesort(:),localsort(:)

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
      integer i,ilev,lmptmp,nmax,idim,j,ijk,ntarg,l,iloc
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg,ifprint,iert
      real *8 time1,time2,pi,done
      real *8 omp_get_wtime

      done = 1
      pi = atan(done)*4.0d0

      
      nexpc = 0
      ntarg = 0
      nlevels = 0
      nboxes = 0

      idivflag = 0
      ndiv = 1



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

      
      call pts_tree_mem(cmpole,nmpole,targ,ntarg,idivflag,ndiv,
     1  nlmin,nlmax,ifunif,iper,nlevels,nboxes,ltree)


      allocate(itree(ltree))
      allocate(boxsize(0:nlevels))
      allocate(tcenters(2,nboxes))

c
c       call the tree code
c

      call pts_tree_build(cmpole,nmpole,targ,ntarg,idivflag,
     1  ndiv,nlmin,nlmax,ifunif,iper,nlevels,nboxes,ltree,itree,
     2  iptr,tcenters,boxsize)

      allocate(isrc(nmpole),isrcse(2,nboxes))

      call pts_tree_sort(nmpole,cmpole,itree,ltree,nboxes,nlevels,iptr,
     1   tcenters,isrc,isrcse)

      allocate(cmpolesort(2,nmpole),rmpolesort(nmpole))
      allocate(mpolesort(ntot),impolesort(nmpole))
      allocate(mtermssort(nmpole))
      allocate(localsort(ntot))
      
      call dreorderf(2,nmpole,cmpole,cmpolesort,isrc)
      call dreorderf(1,nmpole,rmpole,rmpolesort,isrc)
      call ireorderf(1,nmpole,mterms,mtermssort,isrc)

      impolesort(1) = 1
c
c  This loop needs to be openmped
c
      do i=1,nmpole
        mt = mtermssort(i)
        ilen = 2*mt+1
        ijk = 1
        do j=1,ilen
          do l=1,nd
            mpolesort(impolesort(i)+ijk-1) =
     1         mpole(impole(isrc(i))+ijk-1)
             ijk = ijk + 1
          enddo
        enddo
        if(i.lt.nmpole) impolesort(i+1) = impolesort(i) + nd*ilen
      enddo
c
cc      compute scaling factor for multipole/local expansions
c       and lengths of multipole and local expansions
c
      allocate(rscales(0:nlevels),nterms(0:nlevels))

      nmax = 0
      ier = 0
      do i=0,nlevels
        rscales(i) = min(abs(zk*boxsize(i)/(2.0d0*pi)),1.0d0)
cc        rscales(i) = 1
        call h2dterms(boxsize(i),zk,eps,nterms(i),ier)
        nterms(i) = nterms(i) 
        if(nterms(i).gt.nmax) nmax = nterms(i)
      enddo

      if(ifprint.eq.1) call prinf('nmax=*',nmax,1)
      if(ifprint.eq.1) call prinf('nterms=*',nterms,nlevels+1)
      if(ifprint.ge.1) call prin2('rscales=*',rscales,nlevels+1)

c       
c     Multipole and local expansions will be held in workspace
c     in locations pointed to by array iaddr(2,nboxes).
c
c     iiaddr is pointer to iaddr array, itself contained in workspace.
c     imptemp is pointer for single expansion (dimensioned by nmax)
c
c       ... allocate iaddr and temporary arrays
c

      allocate(iaddr(4,nboxes))

      lmptmp = (2*nmax+1)*nd
      allocate(mptemp(lmptmp))



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
      call cpu_time(time1)
C$      time1=omp_get_wtime()
      call hfmm2dmain_mps(nd,eps,
     $   zk,nmpole,cmpolesort,rmpolesort,mtermssort,ntot,
     $   mpolesort,impolesort,localsort,iaddr,rmlexp,lmptot,
     $   mptemp,lmptmp,itree,ltree,iptr,ndiv,nlevels,nboxes,
     $   iper,boxsize,rscales,tcenters,itree(iptr(1)),
     $   isrcse,nterms,timeinfo,ier)
      call cpu_time(time2)
C$        time2=omp_get_wtime()
      if( ifprint .eq. 1 ) call prin2('time in fmm main=*',
     1   time2-time1,1)
      if(ier.ne.0) return
c
c  now unsort the local expansion 
c
      do i=1,nmpole
        mt = mtermssort(i)
        ilen = 2*mt+1
        ijk = 1
        do j=1,ilen
          do l=1,nd
            local(impole(isrc(i))+ijk-1) = 
     1        localsort(impolesort(i)+ijk-1)
            ijk = ijk + 1
          enddo
        enddo
      enddo


      return
      end

      subroutine hfmm2dmain_mps(nd,eps,
     $     zk,nmpole,cmpolesort,rmpolesort,mtermssort,
     $     ntot,mpolesort,impolesort,localsort,iaddr,rmlexp,
     $     lmptot,mptemp,lmptmp,itree,ltree,iptr,ndiv,
     $     nlevels,nboxes,iper,boxsize,rscales,centers,
     $     laddr,isrcse,nterms,timeinfo,ier)
c   Helmholtz multipole scattering FMM in R^2: 
c   evaluate all pairwise multipole
c   interactions (ignoring self-interaction)
c   and compute local expansions
c
c-----------------------------------------------------------------------
c   INPUT PARAMETERS:
c
c   nd:   number of charge densities
c
c   eps:  FMM precision requested
c
c   zk: complex *16, Helmholtz parameter
c   nmpole: number of multipole expansions
c   cmpolesort: real *8 (2,nmpole)
c       location of centers
c   rmpolesort: real *8(nmpole)
c      scaling factors for multipole and local expansions
c   mtermssort: integer(nmpole)
c      number of terms in multipole expansions
c   ntot: integer
c      total length of multipole expansion array
c   mpolesort: complex *16 (ntot)
c      sorted multipole expanson array
c   impolesort: integer (nmpole)
c      how to index into multipole expansionarray
c   iaddr: (4,nboxes): pointer in rmlexp where multipole
c                      and local expansions for each
c                      box is stored
c                      iaddr(1,ibox) is the
c                      starting index in rmlexp for the 
c                      multipole expansion of ibox
c                      and iaddr(2,ibox) is the
c                      starting index in rmlexp
c                      for the local expansion of ibox
c  rmlexp: real *8 (lmptot)
c      multipole expansion for tree hierarchy
c  lmptot: integer
c      length of multipole expansion array for quad tree
c  mptemp: (lmptmp): temporary multipole/local expansion
c                        (may not be needed in new setting)
c  lmptmp: length of temporary expansion
c
c   itree    in: integer (ltree)
c             This array contains all the information
c             about the tree
c             Refer to pts_tree2d.f 
c   ltree    in: integer
c            length of tree
c    iptr in: integer(8)
c             iptr is a collection of pointers 
c             which points to where different elements 
c             of the tree are stored in the itree array
c     ndiv    in: integer
c             Max number of chunks per box
c     nlevels in: integer
c             number of levels in the tree
c     nboxes  in: integer
c             number of boxes in the tree
c     iper    in: integer
c             flag for periodic implementation. Currently unused
c
c     boxsize in: real*8 (0:nlevels)
c             boxsize(i) is the size of the box from end to end
c             at level i
c     rscales: real *8(0:nlevels)
c         scaling parameter for multipole and local expansions
c         at various levels
c     centers in: real *8(2,nboxes)
c                 array containing the centers of all the boxes
c     laddr: integer (2,0:nlevels)
c        laddr(1,i), laddr(2,i) are box indices
c        which start and end at level i
c     
c     isrcse in: integer(2,nboxes)
c               starting and ending location of sources in ibox
c                in sorted list of sources
c     nterms: (0:nlevels) length of multipole and local expansions
c              at various levels
c   OUTPUT
c      localsort: complex *16 (ntot)
c        local expansions at source locations
c      timeinfo: real *8 (8)
c        time taken in various steps of fmm
c      ier: integer
c        error code, ier=0 corresponds to successful execution
c      
c------------------------------------------------------------------
      implicit none

      integer nd

      complex *16 zk
      real *8 zi

      integer nmpole,ier
      integer ndiv,nlevels,ix,iy

      real *8 eps
      integer iper,ifnear

      real *8 cmpolesort(2,nmpole),rmpolesort(nmpole)
      integer mtermssort(nmpole),ntot
      complex *16 mpolesort(ntot),localsort(ntot)
      integer impolesort(nmpole)

      integer iaddr(4,nboxes),lmptmp,lmptot
      real *8 rmlexp(*)
      complex *16 mptemp(lmptmp)
       
      real *8 timeinfo(8)
      real *8 timelev(0:200)
      real *8 centers(2,*)

      integer laddr(2,0:nlevels)
      integer nterms(0:nlevels)
      integer iptr(8),ltree
      integer isrcse(2,nboxes)
      integer itree(*)
      integer nboxes
      real *8 rscales(0:nlevels),boxsize(0:nlevels)


      real *8 zkiupbound,thresh,dn,dx,dy
      real *8 c1(2),c2(2)

      integer nterms_eval(4,0:200)

c     temp variables
      integer i,j,k,l,idim,iloc
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
      integer istarte,iende,istartt,iendt,mt

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
c        call prinf('mnlist4=*',mnlist4,1)
c        call prinf('nlist4s=*',nlist4s,nboxes)
c        call prinf('mnlist3=*',mnlist3,1)
c        call prinf('nlist3s=*',nlist3s,nboxes)


C
c       
      do i=1,8
        timeinfo(i)=0
      enddo

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ntot
        localsort(i) = 0
      enddo
C$OMP END PARALLEL DO

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

c
c
c
      if(ifprint .ge. 1) 
     $   call prinf('=== STEP 1 (shift mp) ====*',i,0)
        call cpu_time(time1)
C$        time1=omp_get_wtime()
c
c       ... step 1, locate all charges, assign them to boxes, and
c       form multipole expansions

      do ilev = 2,nlevels
         if(zi*boxsize(ilev).lt.zkiupbound) then
C
C$OMP PARALLEL DO DEFAULT (SHARED)
C$OMP$PRIVATE(ibox,nchild,istart,iend,npts,i)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox=laddr(1,ilev),laddr(2,ilev)
               nchild = itree(iptr(4)+ibox-1)
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox)
               npts = iend-istart+1
c              Check if current box is a leaf box            
               if(nchild.eq.0.and.npts.gt.0) then
                  do i=istart,iend
                    call h2dmpmp(nd,zk,rmpolesort(i),
     1                 cmpolesort(1,i),mpolesort(impolesort(i)),
     2                 mtermssort(i),rscales(ilev),centers(1,ibox),
     2                 rmlexp(iaddr(1,ibox)),nterms(ilev))
                  enddo
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
C$OMP$PRIVATE(ibox,jbox,istart,iend,npts,i,j,jstart,jend)
C$OMP$SCHEDULE(DYNAMIC)
            do ibox = laddr(1,ilev),laddr(2,ilev)
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox)
               npts = iend-istart+1
               if(npts.ge.0) then
                  do i=1,nlist4s(ibox)
                     jbox = list4(i,ibox)
                     jstart = isrcse(1,jbox)
                     jend = isrcse(2,jbox)
                     do j=jstart,jend
                       call h2dmploc(nd,zk,rmpolesort(j),
     $                    cmpolesort(1,j),mpolesort(impolesort(j)),
     $                    mtermssort(j),rscales(ilev),centers(1,ibox),
     $                    rmlexp(iaddr(2,ibox)),nterms(ilev))
                     enddo
                  enddo
               endif
            enddo
C$OMP END PARALLEL DO        
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
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = npts + iend-istart+1

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
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1
            nchild = itree(iptr(4)+ibox-1)

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




c
c
c

      call cpu_time(time1)
C$    time1=omp_get_wtime()
      if(ifprint.ge.1)
     $    call prinf('=== Step 6 (mp eval) ===*',i,0)

cc      call prinf('ifpgh=*',ifpgh,1)
cc      call prinf('ifpghtarg=*',ifpghtarg,1)
cc      call prinf('laddr=*',laddr,2*(nlevels+1))
      mt = 2*mtermssort(9)+1
      do ilev=1,nlevels-1
       if(zi*boxsize(ilev+1).lt.zkiupbound) then
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,istart,iend,npts,j,i,mt)
C$OMP$PRIVATE(jbox)
C$OMP$SCHEDULE(DYNAMIC)
         do ibox=laddr(1,ilev),laddr(2,ilev)
            istart = isrcse(1,ibox)
            iend = isrcse(2,ibox)
            npts = iend-istart+1
            do i=1,nlist3s(ibox)
               jbox = list3(i,ibox)

               do j=istart,iend
c                 shift multipole expansion directly to box
c                 for all expansion centers
                  call h2dmploc(nd,zk,rscales(ilev+1),
     $             centers(1,jbox),rmlexp(iaddr(1,jbox)),nterms(ilev+1),
     $             rmpolesort(j),cmpolesort(1,j),
     $             localsort(impolesort(j)),mtermssort(j))
               
               enddo
            enddo
         enddo
C$OMP END PARALLEL DO     
       endif
      enddo

 1000 continue    


      call cpu_time(time2)
C$    time2=omp_get_wtime()
      timeinfo(6) = time2-time1

      if(ifprint.ge.1)
     $    call prinf('=== step 7 (LOC to CEN) ===*',i,0)

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
               istart = isrcse(1,ibox)
               iend = isrcse(2,ibox)
               do i=istart,iend
                  call h2dlocloc(nd,zk,rscales(ilev),
     $            centers(1,ibox),
     1            rmlexp(iaddr(2,ibox)),nterms(ilev),
     2            rmpolesort(i),cmpolesort(1,i),
     $            localsort(impolesort(i)),mtermssort(i))
               enddo
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

      do ilev = 0,nlevels
C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,jbox,i,jstart,jend,iloc)
C$OMP$PRIVATE(nlist1,istarts,iends,d,j)
C$OMP$SCHEDULE(DYNAMIC)  
         do ibox = laddr(1,ilev),laddr(2,ilev)
            istarts = isrcse(1,ibox)
            iends = isrcse(2,ibox)
            do iloc = istarts,iends

              do i =1,nlist1s(ibox)
                 jbox = list1(i,ibox) 

                 jstart = isrcse(1,jbox)
                 jend = isrcse(2,jbox)
                 do j=jstart,jend
                   d = (cmpolesort(1,j)-cmpolesort(1,iloc))**2 + 
     1               (cmpolesort(2,j)-cmpolesort(2,iloc))**2
                   d = sqrt(d)
                   if(d.gt.thresh) then
                     call h2dmploc(nd,zk,rmpolesort(j),
     $                 cmpolesort(1,j),
     $                 mpolesort(impolesort(j)),mtermssort(j),
     $                 rmpolesort(iloc),cmpolesort(1,iloc),
     $                 localsort(impolesort(iloc)),mtermssort(iloc))
                   endif
                 enddo
              enddo
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
