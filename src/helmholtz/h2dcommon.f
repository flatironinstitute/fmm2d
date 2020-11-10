cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
c
c  Modified: (C) 2019, Leslie Greengard
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
c      This file contains the basic subroutines for 
c      forming and evaluating multipole (partial wave) expansions
c      in two dimensions.
c
c      Remarks on scaling conventions.
c
c      1)  Hankel and Bessel functions are consistently scaled as
c       	hvec(n)= H_n(z)*rscale^(n)
c       	jvec(n)= J_n(z)/rscale^(n)
c
c          rscale should be of the order of |z| if |z| < 1. Otherwise,
c          rscale should be set to 1.
c
c      2) The potential scaling is that required of the delta function
c      response:  pot = H_0(k*r)*(eye/4)
c
c-----------------------------------------------------------------------
c
c
c
c
c
c**********************************************************************
      subroutine h2cart2polar(zat,r,theta)
      implicit none
c**********************************************************************
c
c     Convert from Cartesian to polar coordinates.
c-----------------------------------------------------------------------
c     INPUT:
c
c     zat   :  Cartesian vector
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     r     :  |zat|
c     theta :  angle of (zat(1),zat(2)) subtended with 
c              respect to x-axis
c-----------------------------------------------------------------------
      real *8 zat(2),r,theta
c 
      r= sqrt(zat(1)**2+zat(2)**2)
      if( abs(zat(1)) .eq. 0 .and. abs(zat(2)) .eq. 0 ) then
      theta = 0
      else
      theta = datan2(zat(2),zat(1))
      endif
      return
      end
c
c
c
c**********************************************************************
      subroutine h2dall(nterms,z,rscale,hvec,ifder,hder)
      implicit none
c**********************************************************************
c
c     This subroutine computes scaled versions of the Hankel 
c     functions H_n of orders 0 to nterms.
c
c       	hvec(n)= H_n(z)*rscale^(n)
c
c     The parameter SCALE is useful when |z| < 1, in which case
c     it damps out the rapid growth of H_n as n increases. In such 
c     cases, we recommend setting 
c                                 
c               rscale approx |z|
c
c     or something close. If |z| > 1, set scale = 1.
c
c     If the flag IFDER is set to one, it also computes the 
c     derivatives of h_n.
c
c		hder(n)= H_n'(z)*rscale^(n)
c
c     NOTE: If |z| < 1.0d-200, the subroutine returns zero.
c     
c-----------------------------------------------------------------------
c     INPUT:
c
c     nterms :  highest order of the Hankel functions to be computed.
c     z      :  argument of the Hankel functions.
c     rscale :  scaling parameter discussed above
c     ifder  :  1 => compute hder array, otherwise do not
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     hvec   :  the vector of Hankel functions 
c     hder   :  the derivatives of the Hankel functions 
c-----------------------------------------------------------------------
      integer i,ifexpon,ifder,nterms
      real *8 thresh,done,rscale,scal2,dtmp 
      complex *16 hvec(0:*),hder(0:*)
      complex *16 zk2,z,zinv,ztmp,fhextra,h0,h1
c
      data thresh/1.0d-200/,done/1.0d0/
c
c     If |z| < thresh, return zeros.
c
      if (abs(z).lt.thresh) then
         do i=0,nterms
            hvec(i)=0
            hder(i)=0
         enddo
         return
      endif
c
c     Otherwise, get H_0 and H_1 via hank103 and the rest via
c     recursion.
c
      ifexpon=1
      call hank103(z,h0,h1,ifexpon)
      hvec(0)=h0
      hvec(1)=h1*rscale
c
c
c     From Abramowitz and Stegun (9.1.27)
c
c       H_{n-1}(z) + H_{n+1}(z) = 2*n/z * H_n(z)
c
c     With scaling:
c
c       H_{n-1}(z) *rscale + H_{n+1}(z)/rscale = 2*n/z * H_n(z)
c       H_{n+1}(z) = rscale*(2*n/z*H_n(z) - H_{n-1}(z)*rscale)
c
      scal2=rscale*rscale
      zinv=rscale/z
      do i=1,nterms-1
         dtmp=2*i
         ztmp=zinv*dtmp
         hvec(i+1)=ztmp*hvec(i)-scal2*hvec(i-1)
      enddo
c
c
c     From Abramowitz and Stegun (9.1.27)
c
c     H_{n}'(z)= H_{n-1}(z) - (n)/z * H_n(z)
c
c     With scaling:
c
c     hder(n)=scale* hvec(n-1) - n/z * hvec(n)
c
      if (ifder.eq.1) then
c
         hder(0)=-hvec(1)/rscale
         zinv=1.0d0/z
         do i=1,nterms
            dtmp=(i)
            ztmp=zinv*dtmp
            hder(i)=rscale*hvec(i-1)-ztmp*hvec(i)
         enddo
c
      endif
c
      return
      end
c
c
c 
c
c
C***********************************************************************
      subroutine h2dmpzero(nd,mpole,nterms)
      implicit none
C***********************************************************************
c
c     This subroutine sets a multipole expansion to zero.
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :   number of expansions
c     nterms :   order of multipole expansion
C---------------------------------------------------------------------
c     OUTPUT:
c
c     mpole  :   coeffs for the expansion set to zero.
C---------------------------------------------------------------------
      integer n,nterms,nd,idim
      complex *16 mpole(nd,-nterms:nterms)
c
      do n=-nterms,nterms
        do idim=1,nd
          mpole(idim,n)=0.0d0
        enddo
      enddo
      return
      end
c
c
