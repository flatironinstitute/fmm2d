cc Copyright (C) 2017: Travis Askham
cc email: askhamwhat@gmail.com      
cc 
cc This software is being released under a modified FreeBSD license
cc (see licenses folder in home directory). 
      
      subroutine modbhgreen_all(beta,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess,ifder3,der3,ifder4,der4,ifder5,der5)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Evaluates the Green function 
c
c     G(x,y) = (S_beta(r) - S_0(r))/beta^2 ,
c
c     where S_beta is the fundamental solution of the Yukawa equation
c     and S_0 is the fundamental solution of the Laplace equation and 
c     r = sqrt( (x_1-y_1)^2 + (x_2-y_2)^2),
c     and its derivatives stably (with the derivatives on the
c     "x" or "target" variable).
c
c     INPUTS
c
c     beta - real *8, modified helmholtz parameter
c     zx(2) - real *8, coordinates of target
c     zy(2) - real *8, coordinates of source
c     ifpot - if equal to 1, compute potential
c     ifgrad - if equal to 1, compute gradient
c     ifhess - if equal to 1, compute hessian
c     ifder3 - if equal to 1, compute 3rd order derivatives
c     ifder4 - if equal to 1, compute 4th order derivatives
c     ifder5 - if equal to 1, compute 5th order derivatives
c
c     OUTPUTS
c
c     note: for derivatives, the ordering is such that,
c     for a given order, the first entry has all of the 
c     derivatives on the x1 variables, the second entry
c     has one derivative on the x2 variable and the rest 
c     on x1, and so on. e.g. 
c
c     hess(1) = d_x1 d_x1 G(x,y)
c     hess(2) = d_x1 d_x2 G(x,y)
c     hess(3) = d_x2 d_x2 G(x,y)
c
c     pot - real *8, potential, if requested
c     grad(2) - real *8, gradient, if requested
c     hess(3) - real *8, hessian, if requested
c     der3(4) - real *8, 3rd order derivatives, if requested.
c     der4(5) - real *8, 4th order derivatives, if requested.
c     der5(6) - real *8, 5th order derivatives, if requested.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, zx(2), zy(2), pot, grad(2), hess(3)
      real *8 der3(4), der4(5), der5(6)
      integer ifpot, ifgrad, ifhess, ifder3, ifder4, ifder5
c     local variables
      real *8 dx, dy, r, r2, r3, r4, r5, r6
      real *8 dx2, dx3, dx4, dx5, dy2, dy3, dy4, dy5
      real *8 g0, g1, g2, g3, g4, g5, betah, dtemp
      real *8 g0temp, g1temp, g2temp, g3temp
      integer if1, if2, if3, if0, ifders
      real *8 diffs(0:5), rscale, pih
      integer nmax, ifder, j
      complex *16 z, eye, hvec(0:10), hder(0:10)
      real *8 kvec(0:10), ders(0:5)
      data eye /(0.0d0,1.0d0)/

      dx = zx(1) - zy(1)
      dy = zx(2) - zy(2)

      dx2 = dx*dx
      dx3 = dx2*dx
      dx4 = dx3*dx
      dx5 = dx4*dx
      dy2 = dy*dy
      dy3 = dy2*dy
      dy4 = dy3*dy
      dy5 = dy4*dy

      r2 = dx2 + dy2
      r = dsqrt(r2)
      r3 = r2*r
      r4 = r3*r
      r5 = r4*r

      nmax = 5

c     get values of difference kernel and higher modes
      
      rscale = 1.0d0
      ifders = 0
c      call diffslogbk(r,beta,rscale,diffs,nmax)
      call diffslogbk_fast(r,beta,rscale,diffs,ifders,ders,kvec,nmax)

c     get values of modified Bessel functions

      z = eye*r*beta
      ifder = 0
      rscale = 1.0d0
c      call h2dall(nmax,z,rscale,hvec,ifder,hder)

      pih = 2.0d0*atan(1.0d0)
      
c      kvec(0) = -dimag(hvec(0))*pih
c      kvec(1) = -dreal(hvec(1))*pih
c      kvec(2) = dimag(hvec(2))*pih
c      kvec(3) = dreal(hvec(3))*pih

c     combine to get values of derivatives of G 
c     with respect to r

      dtemp = 1.0d0/(4.0d0*pih*beta**2)
      betah = 0.5d0*beta

      g0 = diffs(0)*dtemp
      dtemp = -dtemp*beta
      g1 = diffs(1)*dtemp
      dtemp = -dtemp*betah
      g2 = (diffs(2)+kvec(0))*dtemp
      dtemp = -dtemp*betah
      g3 = (diffs(3)+3.0d0*kvec(1))*dtemp
      dtemp = -dtemp*betah
      g4 = (diffs(4)+3.0d0*kvec(0)+4.0d0*kvec(2))*dtemp
      dtemp = -dtemp*betah
      g5 = (diffs(5)+10.0d0*kvec(1)+5.0d0*kvec(3))*dtemp
      
      if0 = 1
      if1 = 1
      if2 = 1
      if3 = 1
c      call difflogbk(r,beta,if0,g0temp,if1,g1temp,if2,g2temp,if3,g3temp)

c      write(*,*) g0, g0temp/(4.0d0*pih*beta**2)
c      write(*,*) g1, g1temp/(4.0d0*pih*beta**2)
c      write(*,*) g2, g2temp/(4.0d0*pih*beta**2)
c      write(*,*) g3, g3temp/(4.0d0*pih*beta**2)
c      write(*,*) g4
c      write(*,*) g5

c     evaluate potential and derivatives

      if (ifpot .eq. 1) then
         pot = g0
      endif
      
      if (ifgrad .eq. 1) then
         grad(1) = dx*g1/r
         grad(2) = dy*g1/r
      endif

      if (ifhess .eq. 1) then
         hess(1) = dx2*g2/r2+g1*(1.0d0/r-dx2/r3)
         hess(2) = dx*dy*(g2/r2-g1/r3)
         hess(3) = dy2*g2/r2+g1*(1.0d0/r-dy2/r3)
      endif

      if (ifder3 .eq. 1) then
         der3(1) = (dx3*g3+3.0d0*dy2*dx*(g2/r-g1/r2))/r3
         der3(2) = dx2*dy*(g3/r3-3.0d0*(g2/r4-g1/r5)) + 
     1        dy*(g2/r2-g1/r3)
         der3(3) = dx*dy2*(g3/r3-3.0d0*(g2/r4-g1/r5)) + 
     1        dx*(g2/r2-g1/r3)
         der3(4) = (dy3*g3+3.0d0*dx2*dy*(g2/r-g1/r2))/r3
      endif

      if (ifder4 .eq. 1) then
         der4(1) = (dx4*(g4-6.0d0*g3/r+15.0d0*(g2/r2-g1/r3)))/r4 +
     1        (dx2*6.0d0*(g3-3.0d0*(g2/r-g1/r2)))/r3 +
     2        (3.0d0*(g2-g1/r))/r2
         der4(2) = (dx3*dy*(g4-6.0d0*g3/r+15.0d0*(g2/r2-g1/r3)))/r4 + 
     1        (3.0d0*dx*dy*(g3-3.0d0*(g2/r-g1/r2)))/r3
         der4(3) = dx2*dy2*(g4-6.0d0*g3/r+15.0d0*g2/r2-15.0d0*g1/r3)/r4
     1        + g3/r - 2.0d0*g2/r2 + 2.0d0*g1/r3
         der4(4) = dx*dy3*(g4-6.0d0*g3/r+15.0d0*(g2/r2-g1/r3))/r4 + 
     1        3.0d0*dx*dy*(g3-3.0d0*(g2/r-g1/r2))/r3
         der4(5) = dy4*(g4-6.0d0*g3/r+15.0d0*(g2/r2-g1/r3))/r4 +
     1        dy2*6.0d0*(g3-3.0d0*(g2/r-g1/r2))/r3 +
     2        3.0d0*(g2-g1/r)/r2        
      endif

      if (ifder5 .eq. 1) then
         der5(1) = (dx5*g5+10.0d0*dy2*dx3*g4/r +
     1        (15.0d0*dy4*dx-30.0d0*dy2*dx3)*g3/r2 +
     2        (60.0d0*dy2*dx3-45.0d0*dy4*dx)*g2/r3 +
     3        (45.0d0*dy4*dx-60.0d0*dy2*dx3)*g1/r4)/r5
         der5(2) = (dy*dx4*g5+(6.0d0*dy3*dx2-4.0d0*dy*dx4)*g4/r +
     1        (3.0d0*dy5+12.0d0*dy*dx4-30.0d0*dy3*dx2)*g3/r2 +
     2        (72.0d0*dy3*dx2-9.0d0*dy5-24.0d0*dy*dx4)*g2/r3 +
     3        (9.0d0*dy5-72.0d0*dy3*dx2+24.0d0*dy*dx4)*g1/r4)/r5
         der5(3) = (dy2*dx3*g5+(3.0d0*dy4*dx-6.0d0*dy2*dx3+dx5)*g4/r +
     1        (27.0d0*dy2*dx3-15.0d0*dy4*dx-3.0d0*dx5)*g3/r2 +
     2        (36.0d0*dy4*dx-63.0d0*dy2*dx3+6.0d0*dx5)*g2/r3 +
     3        (63.0d0*dy2*dx3-36.0d0*dy4*dx-6.0d0*dx5)*g1/r4)/r5
         der5(4) = (dx2*dy3*g5+(3.0d0*dx4*dy-6.0d0*dx2*dy3+dy5)*g4/r +
     1        (27.0d0*dx2*dy3-15.0d0*dx4*dy-3.0d0*dy5)*g3/r2 +
     2        (36.0d0*dx4*dy-63.0d0*dx2*dy3+6.0d0*dy5)*g2/r3 +
     3        (63.0d0*dx2*dy3-36.0d0*dx4*dy-6.0d0*dy5)*g1/r4)/r5
         der5(5) = (dx*dy4*g5+(6.0d0*dx3*dy2-4.0d0*dx*dy4)*g4/r +
     1        (3.0d0*dx5+12.0d0*dx*dy4-30.0d0*dx3*dy2)*g3/r2 +
     2        (72.0d0*dx3*dy2-9.0d0*dx5-24.0d0*dx*dy4)*g2/r3 +
     3        (9.0d0*dx5-72.0d0*dx3*dy2+24.0d0*dx*dy4)*g1/r4)/r5
         der5(6) = (dy5*g5+10.0d0*dx2*dy3*g4/r +
     1        (15.0d0*dx4*dy-30.0d0*dx2*dy3)*g3/r2 +
     2        (60.0d0*dx2*dy3-45.0d0*dx4*dy)*g2/r3 +
     3        (45.0d0*dx4*dy-60.0d0*dx2*dy3)*g1/r4)/r5
      endif


      return
      end

      subroutine modbhgreen(beta,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Evaluates the Green's function (S_beta - S_0)/beta^2 
c     and its derivatives stably (with the derivatives on the
c     "x" or "target" variable).
c
c     INPUTS
c
c     beta - real *8, modified helmholtz parameter
c     zx(2) - real *8, coordinates of target
c     zy(2) - real *8, coordinates of source
c     ifpot - if equal to 1, compute potential
c     ifgrad - if equal to 1, compute gradient
c     ifhess - if equal to 1, compute hessian
c
c     OUTPUTS
c
c     pot - real *8, potential, if requested
c     grad(2) - real *8, gradient, if requested
c     hess(3) - real *8, hessian, if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, zx(2), zy(2), pot, grad(2), hess(3)
      integer ifpot, ifgrad, ifhess
c     local variables
      real *8 der3(4), der4(5), der5(6)
      integer ifder3, ifder4, ifder5

      ifder3 = 0
      ifder4 = 0
      ifder5 = 0
      
      call modbhgreen_all(beta,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess,ifder3,der3,ifder4,der4,ifder5,der5)

      return
      end

      subroutine modbhgreend1(beta,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess,dir1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Evaluates a directional derivative of the Green's function 
c     G = (S_beta - S_0)/beta^2 at the source ("y" variable)
c     and its derivatives stably (with the derivatives on the
c     "x" or "target" variable).
c
c     The potential is (nabla_y * dir1) G
c     The gradient is grad_x (nabla_y * dir1) G
c     The hessian is hess_x (nabla_y * dir1) G
c
c     INPUTS
c
c     beta - real *8, modified helmholtz parameter
c     zx(2) - real *8, coordinates of target
c     zy(2) - real *8, coordinates of source
c     ifpot - if equal to 1, compute potential
c     ifgrad - if equal to 1, compute gradient
c     ifhess - if equal to 1, compute hessian
c
c     OUTPUTS
c
c     pot - real *8, potential, if requested
c     grad(2) - real *8, gradient, if requested
c     hess(3) - real *8, hessian, if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, zx(2), zy(2), pot, grad(2), hess(3), dir1(2)
      integer ifpot, ifgrad, ifhess
c     local variables
      real *8 potloc, gradloc(2), hessloc(3), der3(4), der4(5)
      real *8 der5(6)
      integer ifder3, ifder4, ifder5, ifpotloc, ifgradloc, ifhessloc

      ifpotloc = 0
      ifgradloc = 0
      ifhessloc = 0
      ifder3 = 0
      ifder4 = 0
      ifder5 = 0

      if (ifpot .eq. 1) ifgradloc = 1
      if (ifgrad .eq. 1) ifhessloc = 1
      if (ifhess .eq. 1) ifder3 = 1

      call modbhgreen_all(beta,zx,zy,ifpotloc,potloc,ifgradloc,gradloc,
     1     ifhessloc,hessloc,ifder3,der3,ifder4,der4,ifder5,der5)

      if (ifpot .eq. 1) then
         pot = -dir1(1)*gradloc(1)-dir1(2)*gradloc(2)
      endif

      if (ifgrad .eq. 1) then
         grad(1) = -dir1(1)*hessloc(1) - dir1(2)*hessloc(2)
         grad(2) = -dir1(1)*hessloc(2) - dir1(2)*hessloc(3)
      endif

      if (ifhess .eq. 1) then
         hess(1) = -dir1(1)*der3(1) - dir1(2)*der3(2)
         hess(2) = -dir1(1)*der3(2) - dir1(2)*der3(3)
         hess(3) = -dir1(1)*der3(3) - dir1(2)*der3(4)
      endif
      
      return
      end

      subroutine modbhgreend2(beta,zx,zy,ifpot,pot,ifgrad,grad,
     1     ifhess,hess,dir1,dir2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     Evaluates the Green's function (S_beta - S_0)/beta^2 
c     and its derivatives stably (with the derivatives on the
c     "x" or "target" variable).
c
c     INPUTS
c
c     beta - real *8, modified helmholtz parameter
c     zx(2) - real *8, coordinates of target
c     zy(2) - real *8, coordinates of source
c     ifpot - if equal to 1, compute potential
c     ifgrad - if equal to 1, compute gradient
c     ifhess - if equal to 1, compute hessian
c
c     OUTPUTS
c
c     pot - real *8, potential, if requested
c     grad(2) - real *8, gradient, if requested
c     hess(3) - real *8, hessian, if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, zx(2), zy(2), pot, grad(2), hess(3), dir1(2)
      real *8 dir2(2)
      integer ifpot, ifgrad, ifhess
c     local variables
      real *8 r, r2, dx, dy, g0, g1, g2, g3, temp1, dx2, dy2
      real *8 r3, pi, r4, r5
      real *8 g, gx, gy, gxx, gxy, gyy, gxxx, gxxy, gxyy, gyyy
      real *8 gxxxx, gxxxy, gxxyy, gxyyy, gyyyy
      real *8 quadvec(3)
      integer if0, if1, if2, if3
      integer ifpot1,ifgrad1,ifhess1,ifder3,ifder4,ifder5
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)

      pi = 4.0d0*datan(1.0d0)

c     convert two directional derivatives to xx, xy, yy

      quadvec(1) = dir1(1)*dir2(1)
      quadvec(2) = dir1(1)*dir2(2) + dir1(2)*dir2(1)
      quadvec(3) = dir1(2)*dir2(2)

      ifpot1 = 0
      ifgrad1 = 0
      ifhess1 = 1
      ifder3 = 1
      ifder4 = 1
      ifder5 = 0

      call modbhgreen_all(beta,zx,zy,ifpot1,potloc,ifgrad1,gradloc,
     1     ifhess1,hessloc,ifder3,der3,ifder4,der4,ifder5,der5)

      if (ifpot .eq. 1) then
         pot = quadvec(1)*hessloc(1)+quadvec(2)*hessloc(2)
     1        +quadvec(3)*hessloc(3)
      endif

      if (ifgrad .eq. 1) then
         grad(1) = quadvec(1)*der3(1)+quadvec(2)*der3(2)
     1        +quadvec(3)*der3(3)
         grad(2) = quadvec(1)*der3(2)+quadvec(2)*der3(3)
     1        +quadvec(3)*der3(4)
      endif

      if (ifhess .eq. 1) then
         hess(1) = quadvec(1)*der4(1)+quadvec(2)*der4(2)
     1        +quadvec(3)*der4(3)
         hess(2) = quadvec(1)*der4(2)+quadvec(2)*der4(3)
     1        +quadvec(3)*der4(4)
         hess(3) = quadvec(1)*der4(3)+quadvec(2)*der4(4)
     1        +quadvec(3)*der4(5)
      endif

      return
      end

      subroutine difflogbk(x,beta,if0,g0,if1,g1,if2,g2,if3,g3)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine evaluates the difference kernel
c
c     K_0(beta x) - (-log(x))                     (1)
c
c     where K_0 is the modified Bessel function. 
c
c     INPUT
c
c     if0, if1, if2, if3 - integer, flags. If ifj = 1, then the jth 
c                          x derivative of (1) is returned in the 
c                          variable gj
c     
c     x - real. Argument of the function as in (1)
c
c     beta - real. Yukawa parameter as in (1)
c
c     OUTPUT
c
c     g0, g1, g2, g3 - real. Values of the derivatives of the
c                      difference kernel (1), if requested
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 x, beta, g0, g1, g2, g3
      integer if0, if1, if2, if3
c     local variables
      real *8 yh, yh2, yh3, yh2p, lnyh, y
      real *8 betah, lnbetah
      real *8 i0, i1, i2, i3
      real *8 k0ps, k1ps, k2ps, k3ps
      real *8 k0, k1, k2, k3
      real *8 psi0, psi1, psi2, psi3
      real *8 dtemp, dfac2
      complex *16 z, zk(4)
      integer p
      real *8 gamma
      data gamma /0.5772156649015328606d+00/
      integer pmax
      data pmax / 12 /

      y = x*beta
      yh = 0.5d0*y
      yh2 = yh*yh
      yh3 = yh2*yh
      lnyh = dlog(yh)

      dfac2 = 1
      i0 = 1
      i1 = 1
      i2 = 1/2.0d0
      i3 = 1/6.0d0
      psi0 = -gamma
      psi1 = psi0+1/2.0d0
      psi2 = psi1+1/4.0d0
      psi3 = psi2+1/6.0d0
      k0ps = i0*psi0
      k1ps = i1*psi1
      k2ps = i2*psi2
      k3ps = i3*psi3

      yh2p = yh2
      do p = 1,pmax
         dtemp = yh2p/dfac2
         psi0 = psi0+1.0d0/p
         i0 = i0+dtemp
         k0ps = k0ps+dtemp*psi0

         dtemp = dtemp/(p+1)
         psi1 = psi0+0.5d0/(p+1)
         i1 = i1+dtemp
         k1ps = k1ps+dtemp*psi1

         dtemp = dtemp/(p+2)
         psi2 = psi1+0.5d0/(p+2)
         i2 = i2+dtemp
         k2ps = k2ps+dtemp*psi2

         dtemp = dtemp/(p+3)
         psi3 = psi2+0.5d0/(p+3)
         i3 = i3+dtemp
         k3ps = k3ps+dtemp*psi3
         
         yh2p = yh2p*yh2
         dfac2 = dfac2*(p+1)*(p+1)

      enddo

      i1 = yh*i1
      i2 = yh2*i2
      i3 = yh3*i3

      k0 = k0ps-lnyh*i0
      k1 = -yh*k1ps+lnyh*i1+1/y
      k2 = yh2*k2ps-lnyh*i2+0.5d0/yh2-0.5d0
      k3 = -yh3*k3ps+lnyh*i3+1/yh3-0.5d0/yh+0.25d0*yh
         
      betah = beta*0.5d0
      lnbetah = dlog(betah)

c     g0 = k0+log(x)
      g0 = k0ps-lnyh*(i0-1)-lnbetah

c     g1 = -beta*k1+1/x
      g1 = beta*(yh*k1ps-lnyh*i1)

c     g2 = beta^2/2*(k0+k2)-1/x^2
      g2 = beta**2*0.5d0*(k0ps+yh2*k2ps-lnyh*(i0+i2)-0.5d0)

c     g3 = -beta^3/4*(3*k1+k3)+2/x^3
      g3 = beta**3*(0.75d0*yh*k1ps+0.25d0*yh3*k3ps
     1             -lnyh*(0.75d0*i1+0.25d0*i3)
     2             -0.25d0/yh-0.25d0*0.25d0*yh)

      return
      end


      subroutine diffslogbk(x,beta,rscale,diffs,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine evaluates the difference kernel
c
c     K_0(beta x) - (-log(x))                     (1)
c
c     where K_0 is the modified Bessel function. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 x, beta, diffs(0:n), rscale
      integer n
c     local variables
      real *8 yh, yh2, yhp, yh2p, lnyh, y, rscale2, rscalep
      real *8 betah, lnbetah
      real *8 i0, i1, i2, i3
      real *8 k0ps, k1ps, k2ps, k3ps
      real *8 k0, k1, k2, k3
      real *8 psitemp1, psitemp2
      real *8 dtemp, dfac2
      real *8 psi(60+1000), kps(0:1000), ips(0:1000)
      real *8 fac(0:1000)
      complex *16 z, zk(4), eye
      data eye / (0.0d0,1.0d0) /
      complex *16 hvec(0:1000), hder(0:1000)
      real *8 kvec(0:1000)
      integer p, j, isign, ifder
      real *8 gamma, pih
      data gamma /0.5772156649015328606d+00/
      integer pmax
      data pmax / 12 /

      y = x*beta
      yh = 0.5d0*y
      lnyh = dlog(yh)

      if ( yh .le. 1.0d0 ) then

         yh2 = yh*yh

         psi(1) = -gamma
         do j = 2,pmax+n+1
            psi(j) = psi(j-1)+1.0d0/(j-1)
         enddo
         
         do j = 0,n
            kps(j) = 0.0d0
            ips(j) = 0.0d0
            diffs(j) = 0.0d0
         enddo

         ips(0) = -1.0d0

         yh2p = 1.0d0
         do p = 0,pmax
            dtemp = yh2p
            psitemp1 = psi(p+1)
            do j = 0,n
               psitemp2 = psi(p+1+j)
               kps(j) = kps(j)+(psitemp1+psitemp2)*dtemp
               ips(j) = ips(j)+dtemp
               dtemp = dtemp/(j+p+1.0d0)
            enddo
            yh2p = yh2p*yh2/((p+1.0d0)**2)
         enddo

         yhp = yh
         isign = -1
         rscalep = rscale
         do j = 1,n
            diffs(j) = isign*(0.5d0*yhp*kps(j)-lnyh*ips(j)*yhp)*rscalep
            yhp = yh*yhp
            rscalep = rscalep*rscale
            isign = -isign
         enddo

         diffs(0) = 0.5d0*kps(0)-lnyh*ips(0)-log(0.5d0*beta)

         fac(0) = 1.0d0
         fac(1) = 1.0d0
         
         do j = 2,n
            fac(j) = j*fac(j-1)
         enddo

         rscale2 = rscale*rscale
         do j = 1,n
            dtemp = -0.5d0*(yh/rscale)**(2-j)*rscale2
            do p = 1,j-1
               diffs(j) = diffs(j)+fac(j-p-1)*dtemp/fac(p)
               dtemp = -dtemp*yh2
            enddo
         enddo

      else

         z = eye*y
         ifder = 0
         call h2dall(n+5,z,rscale,hvec,ifder,hder)

         pih = 2.0d0*atan(1.0d0)
         
         do j = 0,n+1,4
            kvec(j) = -dimag(hvec(j))*pih
            kvec(j+1) = -dreal(hvec(j+1))*pih
            kvec(j+2) = dimag(hvec(j+2))*pih
            kvec(j+3) = dreal(hvec(j+3))*pih
         enddo

         diffs(0) = kvec(0) + log(x)

         dfac2 = 1.0d0/yh

         dtemp = 1.0d0/y

         rscalep = rscale
         do j = 1,n
            diffs(j) = (kvec(j) - dtemp)*rscalep
            dtemp = dtemp*dfac2*j
            rscalep = rscalep*rscale
         enddo
      
      endif

      return
      end


      subroutine diffslogbk_fast(x,beta,rscale,diffs,ifders,ders,
     1     kvec,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine evaluates the difference kernel
c
c     K_0(beta x) - (-log(x))                     (1)
c
c     where K_0 is the modified Bessel function. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 x, beta, diffs(0:n), rscale, ders(0:n), kvec(0:n+4)
      integer n, ifders
c     local variables
      real *8 yh, yh2, yhp, yh2p, lnyh, y, rscale2
      real *8 betah, lnbetah
      real *8 i0, i1, i2, i3
      real *8 k0ps, k1ps, k2ps, k3ps
      real *8 k0, k1, k2, k3
      real *8 psitemp1, psitemp2
      real *8 dtemp, dfac2
      real *8 psi(60+2), kps0, kps1, ips0, ips1
      real *8 diffsnp1
      complex *16 hvec(0:1000), hder(0:1000)
      complex *16 eye
      integer ifder
      data eye / (0.0d0,1.0d0) /
      real *8 fac(0:1000)
      complex *16 z, zk(4)
      integer p, j, isign
      real *8 gamma, pih
      data gamma /0.5772156649015328606d+00/

      integer pmax
      parameter (pmax = 12)
      real *8 coeffs0(0:pmax), coeffs1(0:pmax)

      data coeffs0 / -1.1544313298030657d0,    
     1     0.84556867019693427d0,     
     2     1.8455686701969343d0,     
     3     2.5122353368636010d0,     
     4     3.0122353368636010d0,     
     5     3.4122353368636009d0,     
     6     3.7455686701969344d0,     
     7     4.0312829559112204d0,     
     8     4.2812829559112204d0,     
     9     4.5035051781334428d0,    
     $     4.7035051781334429d0,     
     1     4.8853233599516246d0,     
     2     5.0519900266182916d0 /

      data coeffs1 /  -0.15443132980306573d0,     
     1     1.3455686701969343d0,     
     2     2.1789020035302675d0,     
     3     2.7622353368636010d0,     
     4     3.2122353368636007d0,     
     5     3.5789020035302679d0,     
     6     3.8884258130540772d0,     
     7     4.1562829559112204d0,     
     8     4.3923940670223320d0,     
     9     4.6035051781334424d0,     
     $     4.7944142690425338d0,     
     1     4.9686566932849576d0,     
     2     5.1289131035413682d0 /


      y = x*beta
      yh = 0.5d0*y
      yh2 = yh*yh
      lnyh = dlog(yh)
      betah = beta/2.0d0

      z = eye*y
      ifder = 0


      if (yh .gt. 1.0d0 .or. n .gt. 1 .or. ifders .eq. 1) then
         call h2dall(n+4,z,rscale,hvec,ifder,hder)
         pih = 2.0d0*atan(1.0d0)
         do j = 0,n+1,4
            kvec(j) = -dimag(hvec(j))*pih
            kvec(j+1) = -dreal(hvec(j+1))*pih
            kvec(j+2) = dimag(hvec(j+2))*pih
            kvec(j+3) = dreal(hvec(j+3))*pih
         enddo
      endif

      if (yh .le. 1.0d0) then

         kps0 = 0.0d0
         ips0 = -1.0d0
c         ips0 = 0.0d0
         kps1 = 0.0d0
         ips1 = 0.0d0

         yh2p = 1.0d0
         do p = 0,pmax
            kps0 = kps0+coeffs0(p)*yh2p
            ips0 = ips0+yh2p
            yh2p = yh2p/(p+1.0d0)
            kps1 = kps1+coeffs1(p)*yh2p
            ips1 = ips1+yh2p
            yh2p = yh2p*yh2/(p+1.0d0)
         enddo

          diffs(0) = 0.5d0*kps0-lnyh*ips0 - dlog(betah)
c         diffs(0) = 0.5d0*kps0-lnyh*ips0
c         diffs(0) = 0.5d0*kps0-lnyh*ips0 + dlog(x)
         diffs(1) = (lnyh*ips1*yh-0.5d0*yh*kps1)*rscale

      else

         diffs(0) = kvec(0) + dlog(x)

         dtemp = 1.0d0/y
         diffs(1) = kvec(1) - dtemp*rscale

      endif

      rscale2 = rscale*rscale

      dtemp = diffs(1)
      do j = 2,n
         if (dabs(dtemp) .lt. 1.0d-200) dtemp = 0.0d0
         dtemp = 2.0d0*(j-1.0d0)*dtemp*(rscale/y) 
     1        + kvec(j-2)*rscale2
         diffs(j) = dtemp
      enddo

      if (ifders .eq. 1) then

         diffsnp1 = 2.0d0*n*diffs(n)*(rscale/y) + kvec(n-1)*rscale2
         betah = 0.5d0*beta

         ders(0) = -beta*diffs(1)/rscale

         do j = 1,n-1
            ders(j) = -betah*(diffs(j+1)/rscale+kvec(j-1)*rscale)
         enddo
         
         ders(n) = -betah*(diffsnp1/rscale+kvec(n-1)*rscale)

      endif

      return
      end


      subroutine diffszkik(x,beta,rscale,diffs,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine evaluates the difference kernel
c
c     I_0(beta x) - 1                  (1)
c
c     (I_n(beta x) - beta^n x^n/ (2^n n!))/rscale^n
c
c     where I_0 is the modified Bessel function. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 x, beta, diffs(0:n), rscale
      integer n
c     local variables
      real *8 yh, yh2, yhp, yh2p, y
      real *8 dtemp
      real *8, allocatable :: ivals(:)
      complex *16, allocatable :: fjs(:)
      real *8 fjder(10)
      integer, allocatable :: iscale(:)
      complex *16 z, eye
      data eye / (0.0d0,1.0d0) /
      integer p, j, ifder, ntop, lwfjs, ier
      integer pmax
      data pmax / 12 /

      y = x*beta
      yh = 0.5d0*y

      if ( yh .le. 1.0d0 ) then

         yh2 = yh*yh

         do j = 0,n
            diffs(j) = 0.0d0
         enddo

         yh2p = yh2
         do p = 1,pmax
            dtemp = yh2p
            do j = 0,n
               diffs(j) = diffs(j) + dtemp
               dtemp = dtemp/(j+p+1.0d0)
            enddo
            yh2p = yh2p*yh2/((p+1.0d0)**2)
         enddo

         yh2p = 1.0d0
         do j = 0,n
            diffs(j) = yh2p*diffs(j)
            yh2p = yh2p*yh/rscale
         enddo

      else

         lwfjs = n + 5 + 4*n + 100

         allocate(fjs(0:lwfjs),iscale(lwfjs),ivals(0:n+5))

         z = eye*y
         ifder = 0
c         call jfuns2d(ier,n,z,rscale,fjs,ifder,fjder,
c     1	      lwfjs,iscale,ntop)
         call jbessel2d(n,z,rscale,fjs,ifder,fjder)

         do j = 0,n+1,4
            ivals(j) = dreal(fjs(j))
            ivals(j+1) = dimag(fjs(j+1))
            ivals(j+2) = -dreal(fjs(j+2))
            ivals(j+3) = -dimag(fjs(j+3))
         enddo

         diffs(0) = ivals(0) - 1.0d0

         dtemp = yh/rscale
         do j = 1,n
            diffs(j) = ivals(j) - dtemp
            dtemp = dtemp*yh/(rscale*(j+1))
         enddo
      
      endif

      return
      end


      subroutine diffszkik_fast(x,beta,rscale,diffs,ifders,ders,ivec,n)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This subroutine evaluates the difference kernel
c
c     I_0(beta x) - 1                  (1)
c
c     (I_n(beta x) - beta^n x^n/ (2^n n!))/rscale^n
c
c     where I_0 is the modified Bessel function. 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 x, beta, diffs(0:n), rscale, ders(0:n), ivec(0:n)
      integer n, ifders
c     local variables
      real *8 yh, yh2, yhp, yh2p, y, betah
      real *8 dtemp
      real *8, allocatable :: ivals(:)
      complex *16, allocatable :: fjs(:)
      real *8 fjder(10)
      integer, allocatable :: iscale(:)
      complex *16 z, eye
      data eye / (0.0d0,1.0d0) /
      integer p, j, ifder, ntop, lwfjs, ier
      integer pmax
      data pmax / 12 /

      y = x*beta
      yh = 0.5d0*y

      lwfjs = n + 5 + 4*n + 500
      
      allocate(fjs(0:lwfjs),iscale(lwfjs),ivals(0:n+5))
      
      z = eye*y
      ifder = 0
      if (abs(z) .gt. 300) then
         do j = 0,n+5
            fjs(j) = 0.0d0
         enddo
      else
c         call jfuns2d(ier,n+5,z,rscale,fjs,ifder,fjder,
c     1        lwfjs,iscale,ntop)
         call jbessel2d(n+5,z,rscale,fjs,ifder,fjder)
      endif

c      call prin2('fjs *',fjs,2*n+10)

      do j = 0,n+1,4
         ivals(j) = dreal(fjs(j))
         ivals(j+1) = dimag(fjs(j+1))
         ivals(j+2) = -dreal(fjs(j+2))
         ivals(j+3) = -dimag(fjs(j+3))
      enddo

      do j = 0,n
         ivec(j) = ivals(j)
      enddo


      if ( yh .le. 1.0d0 ) then

         yh2 = yh*yh

         do j = 0,n
            diffs(j) = 0.0d0
         enddo

         yh2p = yh2
         do p = 1,pmax
            dtemp = yh2p
            do j = 0,n
               diffs(j) = diffs(j) + dtemp
               if (dtemp .lt. 1.0d-200) dtemp = 0.0d0
               dtemp = dtemp/(j+p+1.0d0)
            enddo
            if (yh2p .lt. 1.0d-200) yh2p = 0.0d0
            yh2p = yh2p*yh2/((p+1.0d0)**2)
         enddo

         yh2p = 1.0d0
         do j = 0,n
            diffs(j) = yh2p*diffs(j)
            if (yh2p .lt. 1.0d-200) yh2p = 0.0d0
            yh2p = yh2p*yh/rscale
         enddo

         if (ifders .eq. 1) then

            betah = 0.5d0*beta
            ders(0) = beta*ivals(1)*rscale

            do j = 1,n
               ders(j) = betah*(diffs(j-1)/rscale+rscale*ivals(j+1))
            enddo
         endif

      else


         diffs(0) = ivals(0) - 1.0d0

         dtemp = yh/rscale
         do j = 1,n
            diffs(j) = ivals(j) - dtemp
            dtemp = dtemp*yh/(rscale*(j+1))
            if (dtemp .gt. 1.0d250) dtemp = 0.0d0
         enddo

         if (ifders .eq. 1) then
            
            betah = 0.5d0*beta
            ders(0) = beta*ivals(1)*rscale

            do j = 1,n
               ders(j) = betah*(diffs(j-1)/rscale+rscale*ivals(j+1))
            enddo
         endif
      
      endif


      return
      end

      
