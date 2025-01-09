c-----------------------------------------------------------------------------
c
c      bh2dterms - determine number of terms in mpole expansions 
c
c      bh2dterms_far - determine number of terms in mpole expansions 
c                      for separation by 2
c
c      bh2dterms_list2 - build the number of terms table for all boxes 
c           in list 2
c
c-----------------------------------------------------------------------------
c
c
c
      subroutine bh2dterms(eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for biharmonic equation
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
      subroutine bh2dterms_far(eps, nterms, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for biharmonic
c     equation
c
c     The method is based on examining the decay of \rho^n / r^{n+1}
c     for \rho a worst case source and r a worst case target.
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
      subroutine bh2dterms_list2(eps, itable, ier)
      implicit real *8 (a-h,o-z)
c
c
c     Determine number of terms in mpole expansions for 
c     biharmonic equation
c
c     The method is based on examining the decay of \rho^n / r^{n+1}
c     for \rho a worst case source and r a worst case target.
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
