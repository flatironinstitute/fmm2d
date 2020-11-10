cc Copyright (C) 2009: Leslie Greengard and Zydrunas Gimbutas
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
c    $Date: 2011-02-22 17:34:23 -0500 (Tue, 22 Feb 2011) $
c    $Revision: 1670 $
c
c
c-----------------------------------------------------------------------------
c
c      h2dterms - determine number of terms in mpole expansions for box
c           of size "bsize" with Helmholtz parameter zk.
c
c      h2dterms_list2 - build the number of terms table for all boxes 
c           in list 2
c
c      h2dterms_list2e - build the number of terms table for all boxes 
c           in extended list 2
c
c-----------------------------------------------------------------------------
c
c
c
      subroutine h2dterms(bsize, zk, eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c     Determine number of terms in mpole expansions for box of size
c     "bsize" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of H_n * J_n.
c
c     Maximum number of terms is 10000, which 
c     works for boxes up to 1000 wavelengths in size     
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:20000), ht0,
     1             ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:20000), fhder(0:1)
c
      ier = 0
      pi2 = 8*datan(1.0d0)
c
      z1 = (zk*bsize)*1.5d0
c
c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c       
ccc        if( abs(z1) .lt. 1d-6 ) z1=1d-6
c
      ntmax = 10000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*bsize) .lt. pi2) rscale = cdabs(zk*bsize)
      call h2dall(ntmax,z1,rscale,hfun,ifder,fhder)
c
      z2 = (zk*bsize) * dsqrt(2d0)/2.d0
c
      call jbessel2d(ntmax, z2, rscale, jfun, ifder, fjder)
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      xtemp2 = cdabs(jfun(1)*hfun(1))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        xtemp2 = cdabs(jfun(j-1)*hfun(j-1))
        xtemp = (xtemp1+xtemp2)*cdabs(hfun(0))
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          return
        endif
c
      enddo
c
c       ... computational box is too big, set nterms to 10001
c
        ier = 13
        nterms=10001
c
      return
      end
c
c
c-------------------------------------------------------------
      subroutine h2dterms2(bsize, zk, eps, nterms, np, ier)
      implicit real *8 (a-h,o-z)
c
c     Determine number of terms in mpole expansions for box of size
c     "bsize" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of H_n+np * J_n.
c
c     Maximum number of terms is 10000, which 
c     works for boxes up to 1000 wavelengths in size     
c
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, z4,jfun(0:20000), jfun2(0:20000),
     1             ht0,ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:20000), fhder(0:1)
c
      ier = 0
c
      z1 = (zk*bsize)*1.5d0
c
c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c       
ccc        if( abs(z1) .lt. 1d-6 ) z1=1d-6
c
      ntmax = 10000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*bsize) .lt. 1.0d0) rscale = cdabs(zk*bsize)
      call h2dall(ntmax,z1,rscale,hfun,ifder,fhder)
c
      z2 = (zk*bsize) * dsqrt(2d0)/2.d0
c
      call jbessel2d(ntmax, z2, rscale, jfun, ifder, fjder)

      z4 = 0.001*zk/2
cc       z4 = 5.0d0/zk
      call prin2('z4=*',z4,2)
      call jbessel2d(ntmax,z4,rscale,jfun2,ifder,fjder)

cc      call prin2('jfun2=*',jfun2,20)
c
      xtemp1 = cdabs(jfun(0)*hfun(np)*jfun2(np))
      xtemp2 = cdabs(jfun(1)*hfun(np+1)*jfun2(np))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax-np
        xtemp1 = cdabs(jfun(j)*hfun(j+np)*jfun2(np))
        xtemp2 = cdabs(jfun(j-1)*hfun(j+np-1)*jfun2(np))
        xtemp = xtemp1+xtemp2
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          return
        endif
c
      enddo
c
c       ... computational box is too big, set nterms to 10001
c
        ier = 13
        nterms=10001
c
      return
      end
c------------------------------------------------------------------------      
c
c
c
      subroutine h2dterms_far(bsize, zk, eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "bsize" with Helmholtz parameter zk. 
c
c     The method is based on examining the decay of H_n * J_n.
c
c     This routine aasumes slightly larger separation of boxes: the
c     first unit box is located at the origin and the second box is
c     located at (3,0,0).
c
c     Maximum number of terms is 10000, which 
c     works for boxes up to 1000 wavelengths in size     
c
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:20000), ht0,
     1             ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:20000), fhder(0:1)
c
      ier = 0
c
      z1 = (zk*bsize)*2.5d0
c
c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c       
ccc        if( abs(z1) .lt. 1d-6 ) z1=1d-6
c
      ntmax = 10000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*bsize) .lt. 1.0d0) rscale = cdabs(zk*bsize)
      call h2dall(ntmax,z1,rscale,hfun,ifder,fhder)
c
      z2 = (zk*bsize) * dsqrt(2d0)/2.d0
c
      call jbessel2d(ntmax, z2, rscale, jfun, ifder, fjder)
ccc      call prin2(' jfun is *',jfun,2*ntmax+2)
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      xtemp2 = cdabs(jfun(1)*hfun(1))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        xtemp2 = cdabs(jfun(j-1)*hfun(j-1))
        xtemp = xtemp1+xtemp2
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          return
        endif
      enddo
c
c       ... computational box is too big, set nterms to 10001
c
        ier = 13
        nterms=10001
c
      return
      end
c
c
c
c
c
      subroutine h2dterms_list2(bsize, zk, eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "bsize" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in list 2
c
c     Maximum number of terms is 1000, which 
c     works for boxes up to 160 wavelengths in size     
c
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:20000), ht0,
     1             ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:20000), fhder(0:1)
c
      integer nterms_table(2:3,0:3)
      integer itable(-3:3,-3:3)
c
      ier = 0
c
c
c       first box in x direction
c
ccc      z1 = (zk*size)*1.5d0
c
c       second box in x direction
c
ccc      z1 = (zk*size)*2.5d0
c
c       on the diagonal
c
c      z1 = (zk*size)*dsqrt(2d0)/2.d0*3
c      z1 = (zk*size)*dsqrt(2d0)/2.d0*5
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
ccc        call prin2('rr=*',sqrt(2.0d0)/2*5,1)
c
        z1 = (zk*bsize)*rr

c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c  
ccc      if( abs(z1) .lt. 1d-6 ) z1=1d-6
c
      ntmax = 10000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*bsize) .lt. 1.0d0) rscale = cdabs(zk*bsize)
      call h2dall(ntmax,z1,rscale,hfun,ifder,fhder)
c
      z2 = (zk*bsize) * dsqrt(2d0)/2.d0
c
      call jbessel2d(ntmax, z2, rscale, jfun, ifder, fjder)
ccc      call prin2(' jfun is *',jfun,2*ntmax+2)
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      xtemp2 = cdabs(jfun(1)*hfun(1))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        xtemp2 = cdabs(jfun(j-1)*hfun(j-1))
        xtemp = xtemp1+xtemp2
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          goto 1600
        endif
      enddo
c
c       ... computational box is too big, set nterms to 10001
c
        ier = 13
        nterms=10001
c
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
      subroutine h2dterms_list2e(bsize, zk, eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "bsize" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Build nterms table for all boxes in extended list 2
c
c     Maximum number of terms is 10000, which 
c     works for boxes up to 1000 wavelengths in size     
c
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:20000), ht0,
     1             ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:20000), fhder(0:1)
c
      integer nterms_table(2:7,0:7)
      integer itable(-7:7,-7:7)
c
      ier = 0
c
c
c       first box in x direction
c
ccc      z1 = (zk*size)*1.5d0
c
c       second box in x direction
c
ccc      z1 = (zk*size)*2.5d0
c
c       on the diagonal
c
c      z1 = (zk*size)*dsqrt(2d0)/2.d0*3
c      z1 = (zk*size)*dsqrt(2d0)/2.d0*5
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
        z1 = (zk*bsize)*rr

c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c  
ccc      if( abs(z1) .lt. 1d-6 ) z1=1d-6
c
      ntmax = 10000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*bsize) .lt. 1.0d0) rscale = cdabs(zk*bsize)
      call h2dall(ntmax,z1,rscale,hfun,ifder,fhder)
c
      z2 = (zk*bsize) * dsqrt(2d0)/2.d0
c
      call jbessel2d(ntmax, z2, rscale, jfun, ifder, fjder)
ccc      call prin2(' jfun is *',jfun,2*ntmax+2)
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      xtemp2 = cdabs(jfun(1)*hfun(1))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        xtemp2 = cdabs(jfun(j-1)*hfun(j-1))
        xtemp = xtemp1+xtemp2
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          goto 1600
        endif
      enddo
c
c       ... computational box is too big, set nterms to 10001
c
        ier = 13
        nterms=10001
c
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
      subroutine h2dterms_eval(itype, bsize, zk, eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for box of size
c     "bsize" with Helmholtz parameter zk.
c
c     The method is based on examining the decay of h_n * j_n.
c
c     Maximum number of terms is 10000, which 
c     works for boxes up to 1000 wavelengths in size     
c
c
c-----------------------------------------------------------------------------
c
      complex *16  zk, z1, z2, z3, jfun(0:20000), ht0,
     1             ht1, ht2, fjder(0:1), ztmp,
     1             hfun(0:20000), fhder(0:1)
c
      ier = 0
c
      z1 = (zk*bsize)*1.5d0
c
c       the code will run out memory if frequency is too small 
c       set frequency to something more reasonable, nterms is 
c       approximately the same for all small frequencies
c       
ccc      if( abs(z1) .lt. 1d-6 ) z1=1d-6
c
      ntmax = 1000
      ifder = 0
      rscale = 1.0d0
      if (cdabs(zk*bsize) .lt. 1.0d0) rscale = cdabs(zk*bsize)
      call h2dall(ntmax,z1,rscale,hfun,ifder,fhder)
c
        z2 = (zk*bsize) * dsqrt(2d0)/2.d0
c
c       corners included
        if( itype .eq. 1 ) z2 = (zk*bsize) * dsqrt(2d0)/2.d0
c       edges included, no corners
        if( itype .eq. 2 ) z2 = (zk*bsize) * dsqrt(1d0)/2.d0
c       center only
        if( itype .eq. 3 ) z2 = (zk*bsize) * 1.0d0/2.d0
c       center only, small interior sphere
        if( itype .eq. 4 ) z2 = (zk*bsize) * 0.8d0/2.d0
c
      call jbessel2d(ntmax, z2, rscale, jfun, ifder, fjder)
ccc      call prin2(' jfun is *',jfun,2*ntmax+2)
c
      xtemp1 = cdabs(jfun(0)*hfun(0))
      xtemp2 = cdabs(jfun(1)*hfun(1))
      xtemp0 = xtemp1+xtemp2
      nterms = 1
      do j = 2, ntmax
        xtemp1 = cdabs(jfun(j)*hfun(j))
        xtemp2 = cdabs(jfun(j-1)*hfun(j-1))
        xtemp = xtemp1+xtemp2
        if(xtemp .lt. eps*xtemp0)then
          nterms = j + 1
          return
        endif
c
      enddo
c
c       ... computational box is too big, set nterms to 10001
c
        ier = 13
        nterms=10001
c
      return
      end
c
c
c
