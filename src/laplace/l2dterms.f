cc Copyright (C) 2009-2011: Leslie Greengard and Zydrunas Gimbutas
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
c    $Date: 2010-12-28 14:28:30 -0500 (Tue, 28 Dec 2010) $
c    $Revision: 1571 $
c
c
c-----------------------------------------------------------------------------
c
c      l2dterms - determine number of terms in mpole expansions 
c
c      l2dterms_list2 - build the number of terms table for all boxes 
c           in list 2
c
c      l2dterms_list2w - build the number of terms table for all boxes 
c           in list 2, worst case error in multipole to local translation 
c
c      l2dterms_list2ew - build the number of terms table for all boxes
c           in extended list 2, worst case error in multipole to local
c           translation
c
c-----------------------------------------------------------------------------
c
c
c
      subroutine l2dterms(eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions.
c
c     The method is based on examining the decay of \rho^n / r^{n+1}
c     for \rho a worst case source and r a worst case target.
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      ier = 0
c
      ntmax = 1000
c       
      z1 = 1.5d0
      do i = 0,ntmax
         hfun(i) = 1.0d0/(z1**(i+1))
      enddo
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
      z2 = dsqrt(2d0)/2.d0
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          return
        endif
c
      enddo
      return
      end
c
c
c
c
c
      subroutine l2dterms_far(eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk=0. 
c
c     The method is based on examining the decay of h_n * j_n.
c
c     This routine assumes slightly larger separation of boxes: the
c     first unit box is located at the origin and the second box is
c     located at (3,0).
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      ier = 0
c
      ntmax = 1000
c
      z1 = 2.5d0
      do i = 0,ntmax
         hfun(i) = 1.0d0/(z1**(i+1))
      enddo
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
      z2 = dsqrt(2d0)/2.d0
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
ccc      call prin2(' jfun is *',jfun,2*ntmax+2)
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          return
        endif
c
      enddo
      return
      end
c
c
c
c
c
      subroutine l2dterms_list2(eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in list 2
c
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      integer nterms_table(2:3,0:3)
      integer itable(-3:3,-3:3)
c
      ier = 0
c
      do 1800 ii=2,3
      do 1800 jj=0,3
c
        dx=ii
        dy=jj
c       
        if( dx .gt. 0 ) dx=dx-.5
        if( dy .gt. 0 ) dy=dy-.5
c
        rr=sqrt(dx*dx+dy*dy)
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(3.0d0)/2*5,1)
c
      ntmax = 1000
c
      z1 = rr
      do i = 0,ntmax
         hfun(i) = 1.0d0/( z1**(i+1))
      enddo
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
      z2 = dsqrt(2d0)/2.d0
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          goto 1600
        endif
      enddo
 1600   continue
c
        nterms_table(ii,jj)=nterms
c
 1800   continue
c
ccc        call prinf('nterms=*',nterms_table,2*4*4)
c
c       build the rank table for all boxes in list 2
c
        do i=-3,3
        do j=-3,3
        itable(i,j)=0
        enddo
        enddo
c
        do 2200 i=-3,3
        do 2200 j=-3,3
c
        if( abs(i) .gt. 1 ) then
        itable(i,j)=nterms_table(abs(i),abs(j))
        else if( abs(j) .gt. 1) then
        itable(i,j)=nterms_table(abs(j),abs(i))
        endif
c
 2200   continue
c
      return
      end
c
c
c
c
c
      subroutine l2dterms_list2w(eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in list 2
c     Estimate worst case multipole to local translation operator errors.
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      integer nterms_table(2:3,0:3)
      integer itable(-3:3,-3:3)
c
      ier = 0
c
      do 1800 ii=2,3
      do 1800 jj=0,3
c
        dx=ii
        dy=jj
c       
c        if( dx .gt. 0 ) dx=dx-.5
c        if( dy .gt. 0 ) dy=dy-.5
c
        rr=sqrt(dx*dx+dy*dy)
        rr=rr-sqrt(2.0d0)/2
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(3.0d0)/2*5,1)
c
      ntmax = 1000
c
      z1 = rr
      do i = 0,ntmax
         hfun(i) = 1.0d0/( z1**(i+1))
      enddo
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
      z2 = dsqrt(2d0)/2.d0
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          goto 1600
        endif
      enddo
 1600   continue
c
        nterms_table(ii,jj)=nterms
c
 1800   continue
c
ccc        call prinf('nterms=*',nterms_table,2*4*4)
c
c       build the rank table for all boxes in list 2
c
        do i=-3,3
        do j=-3,3
        itable(i,j)=0
        enddo
        enddo
c
        do 2200 i=-3,3
        do 2200 j=-3,3
c
        if( abs(i) .gt. 1 ) then
        itable(i,j)=nterms_table(abs(i),abs(j))
        else if( abs(j) .gt. 1) then
        itable(i,j)=nterms_table(abs(j),abs(i))
        endif
c
 2200   continue
c
      return
      end
c
c
c
c
c
      subroutine l2dterms_list2e(eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in extended list 2
c
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      integer nterms_table(2:7,0:7)
      integer itable(-7:7,-7:7)
c
      ier = 0
c
      do 1800 ii=2,7
      do 1800 jj=0,7
c
        dx=ii
        dy=jj
c       
        if( dx .gt. 0 ) dx=dx-.5
        if( dy .gt. 0 ) dy=dy-.5
c
        rr=sqrt(dx*dx+dy*dy)
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(2.0d0)/2*5,1)
c
      ntmax = 1000
c
      z1 = rr
      do i = 0,ntmax
         hfun(i) = 1.0d0/( z1**(i+1))
      enddo
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
      z2 = dsqrt(2d0)/2.d0
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          goto 1600
        endif
      enddo
 1600   continue
c
        nterms_table(ii,jj)=nterms
c
 1800   continue
c
ccc        call prinf('nterms=*',nterms_table,2*4*4)
c
c       build the rank table for all boxes in extended list 2
c
        do i=-7,7
        do j=-7,7
        itable(i,j)=0
        enddo
        enddo
c
        do 2200 i=-7,7
        do 2200 j=-7,7
c
        if( abs(i) .gt. 2 ) then
        itable(i,j)=nterms_table(abs(i),abs(j))
        else if( abs(j) .gt. 2) then
        itable(i,j)=nterms_table(abs(j),abs(i))
        endif
c
 2200   continue
c
      return
      end
c
c
c
c
c
      subroutine l2dterms_list2ew(eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in extended list 2
c     Estimate worst case multipole to local translation operator errors.
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      integer nterms_table(2:7,0:7)
      integer itable(-7:7,-7:7)
c
      ier = 0
c
      do 1800 ii=2,7
      do 1800 jj=0,7
c
        dx=ii
        dy=jj
c       
c        if( dx .gt. 0 ) dx=dx-.5
c        if( dy .gt. 0 ) dy=dy-.5
c
        rr=sqrt(dx*dx+dy*dy)
        rr=rr-sqrt(2.0d0)/2
ccc        call prin2('rr=*',rr,1)
ccc        call prin2('rr=*',sqrt(2.0d0)/2*5,1)
c
      ntmax = 1000
c
      z1 = rr
      do i = 0,ntmax
         hfun(i) = 1.0d0/( z1**(i+1))
      enddo
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
      z2 = dsqrt(2d0)/2.d0
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          goto 1600
        endif
      enddo
 1600   continue
c
        nterms_table(ii,jj)=nterms
c
 1800   continue
c
ccc        call prinf('nterms=*',nterms_table,2*4*4)
c
c       build the rank table for all boxes in extended list 2
c
        do i=-7,7
        do j=-7,7
        itable(i,j)=0
        enddo
        enddo
c
        do 2200 i=-7,7
        do 2200 j=-7,7
c
        if( abs(i) .gt. 2 ) then
        itable(i,j)=nterms_table(abs(i),abs(j))
        else if( abs(j) .gt. 2) then
        itable(i,j)=nterms_table(abs(j),abs(i))
        endif
c
 2200   continue
c
      return
      end
c
c
c
c
c
      subroutine l2dterms_eval(itype, eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "size" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:2000), ht0,
     1             ht1, ht2, ztmp,
     1             hfun(0:2000)
c
      ier = 0
      ntmax = 1000
c
      z1 = 1.5d0
      do i = 0,ntmax
         hfun(i) = 1.0d0/(z1**(i+1))
      enddo
c
ccc      call prin2(' hfun is *',hfun,2*ntmax+2)
c
        z2 = dsqrt(2d0)/2.d0
c
c       corners included
        if( itype .eq. 1 ) z2 = dsqrt(2d0)/2.d0
c       edges included, no corners
        if( itype .eq. 2 ) z2 = dsqrt(1d0)/2.d0
c       center only
        if( itype .eq. 3 ) z2 = 1.0d0/2.d0
c       center only, small interior sphere
        if( itype .eq. 4 ) z2 = 0.8d0/2.d0
c
c
      do i = 0,ntmax
         jfun(i) = z2**i
      enddo
ccc      call prin2(' jfun is *',jfun,2*ntmax+2)
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        if(xtemp1 .lt. eps)then
          nterms = j
          return
        endif
c
      enddo
      return
      end

