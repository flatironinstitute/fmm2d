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
      subroutine hfmm2d_mps_direct(nd,eps,zk,nmpole,cmpole,rmpole,
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
c   timeinfo(8)     : time taken in various steps
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
      real *8 timeinfo(8)

      call hfmm2d_mps_ndiv(nd,eps,zk,nmpole,cmpole,rmpole,
     1  mterms,ntot,mpole,impole,local,nmpole,timeinfo,ier)



      return
      end
c
c
c
c
c

      subroutine hfmm2d_mps_ndiv(nd,eps,zk,nmpole,cmpole,rmpole,
     1  mterms,ntot,mpole,impole,local,ndiv,timeinfo,ier)
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
c   ndiv             : subdivision criterion, set it to nmpole to do
c                      things directly
c
c   OUTPUT PARAMETERS
c   local(nd,ntot)  : local expansion at each center, due to all 
c                     incoming multipole expansions (self is ignored).
c                     the orders are same as for the incoming multipole
c   timeinfo(8)     : time taken in various steps
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
