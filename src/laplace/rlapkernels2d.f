c
c     Laplace interactions
c      
c     real-valued charge, dipstr, pot, grad, hess
c     real-valued dipvec
c      
c**********************************************************************
      subroutine r2d_directcp_vec(nd,sources,ns,charge,
     $           targ,pot,thresh)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c
c     pot(ii)  = \sum_j log(|targ-source(*,j)|)*charge(ii,j)
c
c     The potential is not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ          :   location of the target
c     thresh        :   threshold for computing potential
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)   (real *8)      : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(2,ns),targ(2),xdiff,ydiff,rr,r
      real *8 thresh,rtmp,thresh2
      real *8 pot(nd)
      real *8 charge(nd,ns)

      thresh2 = thresh*thresh
c
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.le.thresh2) goto 1000

         rtmp = log(rr)/2 
         do ii=1,nd
            pot(ii) = pot(ii) + rtmp*charge(ii,i)
         enddo
 1000    continue         
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine r2d_directcg_vec(nd,sources,ns,charge,targ,pot,
     1             grad,thresh)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c
c     pot(ii)  = \sum_j log(|targ-source(*,j)|)*charge(ii,j)
c
c     grad(1:2,ii)  = \nabla pot(ii)
c
c     The potential and gradient are not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ          :   location of the target
c     thresh        :   threshold for computing potential and gradient
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (real *8)      : potential is incremented
c     grad(nd,2)  (real *8)      : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(2,ns),targ(2)
      real *8 xdiff,ydiff,rr,r,thresh,rtmp,thresh2
      real *8 pot(nd),grad(nd,2)
      real *8 dx,dy
      real *8 charge(nd,ns)

      thresh2 = thresh*thresh
c
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         rtmp = log(rr)/2
         dx = xdiff/rr
         dy = ydiff/rr
c
         do ii = 1,nd
            pot(ii) = pot(ii) + rtmp*charge(ii,i)
            grad(ii,1) = grad(ii,1) + dx*charge(ii,i)
            grad(ii,2) = grad(ii,2) + dy*charge(ii,i)
         enddo
 1000 continue         
      enddo
      
      return
      end
c
c
c
c
c
c
c
c
c**********************************************************************
      subroutine r2d_directch_vec(nd,sources,ns,charge,targ,
     1           pot,grad,hess,thresh)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c
c     pot(ii)  = \sum_j log(|targ-source(*,j)|)*charge(ii,j)
c
c     grad  = d/dz(pot)
c     hess = d^2/dz^2(pot)
c
c     The potential and gradient are not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ          :   location of the target
c     thresh        :   threshold for computing potential,gradient,and
c                       hessian
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (real *8)      : potential is incremented
c     grad(nd,2)  (real *8)      : gradient is incremented
c     hess(nd,3)  (real *8)      : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd
      real *8 sources(2,ns),targ(2)
      real *8 xdiff,ydiff,rr,r,thresh,thresh2
      real *8 pot(nd),grad(nd,2),hess(nd,3)
      real *8 rtmp, xdiff2, ydiff2
      real *8 rr2, dx, dy, dxx, dxy, dyy
      real *8 charge(nd,ns)

      thresh2 = thresh*thresh

c
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         xdiff2 = xdiff*xdiff
         ydiff2 = ydiff*ydiff
         rr=xdiff2+ydiff2
         rr2 = rr*rr
         if(rr.lt.thresh2) goto 1000
         rtmp = log(rr)/2
         dx = xdiff/rr
         dy = ydiff/rr
         dxx = (rr - 2*xdiff2)/rr2
         dxy = -2*xdiff*ydiff/rr2
         dyy = (rr - 2*ydiff2)/rr2
         
c
         do ii = 1,nd
            pot(ii) = pot(ii) + rtmp*charge(ii,i)
            grad(ii,1) = grad(ii,1) + dx*charge(ii,i)
            grad(ii,2) = grad(ii,2) + dy*charge(ii,i)
            hess(ii,1) = hess(ii,1) + dxx*charge(ii,i)
            hess(ii,2) = hess(ii,2) + dxy*charge(ii,i)
            hess(ii,3) = hess(ii,3) + dyy*charge(ii,i)  
         enddo
 1000 continue         
      enddo
      

      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine r2d_directdp_vec(nd,sources,ns,dipstr,dipvec,
     $           targ,pot,thresh)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c     pot(ii)  = \sum_j dipstr(ii,j)* dipstr(ii,:,j) \cdot
c                               \nabla_src \log |targ(:)-sources(:,j)|
c
c     The potential is not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     dipstr(nd,ns) :   dipole strengths
c     dipvec(nd,2,ns) :   dipole orientations
c     targ          :   location of the target
c     thresh        :   threshold for computing potential
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)   (real *8)      : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(2,ns),targ(2),xdiff,ydiff,rr,r
      real *8 thresh,thresh2,p1,p2
      real *8 pot(nd)
      real *8 dipstr(nd,ns)
      real *8 dipvec(nd,2,ns)
c
      thresh2 = thresh*thresh
      
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr = xdiff*xdiff + ydiff*ydiff
         if(rr.le.thresh2) goto 1000

         p1 = -xdiff/rr
         p2 = -ydiff/rr         

         do ii=1,nd
            pot(ii) = pot(ii) + dipstr(ii,i)*(dipvec(ii,1,i)*p1
     1           + dipvec(ii,2,i)*p2)
         enddo
 1000    continue         
      enddo
      return
      end
c
c
c
c
c
c
c**********************************************************************
      subroutine r2d_directdg_vec(nd,sources,ns,dipstr,dipvec,
     1     targ,pot,grad,thresh)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c     pot(ii)  = \sum_j dipstr(ii,j)* dipstr(ii,:,j) \cdot
c                               \nabla_src \log |targ(:)-sources(:,j)|
c
c     grad  = \nabla(pot)
c
c     The potential and gradient are not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     dipstr(nd,ns) :   charge strengths
c     dipvec(nd,2,ns) :   dipole orientations
c     targ          :   location of the target
c     thresh        :   threshold for computing potential and gradient
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (real *8)    : potential is incremented
c     grad(nd,2)  (real *8)      : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(2,ns),targ(2)
      real *8 xdiff,ydiff,rr,r,thresh,rtmp
      real *8 pot(nd),grad(nd,2)
      real *8 dipstr(nd,ns),d1,d2
      real *8 dx1,dx2,dy1,dy2
      real *8 dipvec(nd,2,ns)
      real *8 xdiff2,ydiff2,p1,p2,thresh2,rr2
c
      thresh2 = thresh*thresh
      
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         xdiff2 = xdiff*xdiff
         ydiff2 = ydiff*ydiff
         rr = xdiff2 + ydiff2
         if(rr.le.thresh2) goto 1000 

         rr2 = rr*rr
         
         p1 = -xdiff/rr
         p2 = -ydiff/rr

         dx1 = -(rr-2*xdiff2)/rr2
         dx2 = 2*xdiff*ydiff/rr2
         dy1 = dx2
         dy2 = -(rr-2*ydiff2)/rr2

c         
         do ii = 1,nd
            d1 = dipstr(ii,i)*dipvec(ii,1,i)
            d2 = dipstr(ii,i)*dipvec(ii,2,i)            
            pot(ii) = pot(ii) + d1*p1 + d2*p2
            grad(ii,1) = grad(ii,1) + d1*dx1 + d2*dx2
            grad(ii,2) = grad(ii,2) + d1*dy1 + d2*dy2 
         enddo
 1000 continue         
      enddo
      
      return
      end
c
c
c
c
c
c
c
c
c**********************************************************************
      subroutine r2d_directdh_vec(nd,sources,ns,dipstr,dipvec,targ,
     1           pot,grad,hess,thresh)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c     pot(ii)  = \sum_j dipstr(ii,j)* dipstr(ii,:,j) \cdot
c                               \nabla_src \log |targ(:)-sources(:,j)|
c
c     grad  = \nabla (pot)
c     hess = \nabla \nabla (pot) (dxx,dxy,dyy)
c
c     The potential and gradient are not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     dipstr(nd,ns) :   charge strengths
c     dipvec(nd,2,ns) :   dipole orientations      
c     targ          :   location of the target
c     thresh        :   threshold for computing potential,gradient,and
c                       hessian
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (real *8)      : potential is incremented
c     grad(nd,2)  (real *8)      : gradient is incremented
c     hess(nd,3)  (real *8)      : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd
      real *8 sources(2,ns),targ(2)
      real *8 xdiff,ydiff,rr,r,thresh
      real *8 pot(nd),grad(nd,2),hess(nd,3)
      real *8 rtmp
      real *8 dipstr(nd,ns),d1,d2
      real *8 dipvec(nd,2,ns),dx1,dx2,dy1,dy2
      real *8 xdiff2,ydiff2,p1,p2,dxx1,dxx2,dxy1,dxy2,dyy1,dyy2
      real *8 rr2,rr3,thresh2

c
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         xdiff2 = xdiff*xdiff
         ydiff2 = ydiff*ydiff
         rr = xdiff2 + ydiff2
         if(rr.le.thresh2) goto 1000 

         rr2 = rr*rr
         rr3 = rr2*rr
         
         p1 = -xdiff/rr
         p2 = -ydiff/rr

         dx1 = -(rr-2*xdiff2)/rr2
         dx2 = 2*xdiff*ydiff/rr2
         dy1 = dx2
         dy2 = -(rr-2*ydiff2)/rr2

         dxx1 = -2*xdiff*(xdiff2-3*ydiff2)/rr3
         dxx2 = 2*ydiff*(ydiff2-3*xdiff2)/rr3
         dxy1 = dxx2
         dxy2 = -dxx1
         dyy1 = dxy2
         dyy2 = -dxx2

c         
         do ii = 1,nd
            d1 = dipstr(ii,i)*dipvec(ii,1,i)
            d2 = dipstr(ii,i)*dipvec(ii,2,i)            
            pot(ii) = pot(ii) + d1*p1 + d2*p2
            grad(ii,1) = grad(ii,1) + d1*dx1 + d2*dx2
            grad(ii,2) = grad(ii,2) + d1*dy1 + d2*dy2
            hess(ii,1) = hess(ii,1) + d1*dxx1 + d2*dxx2
            hess(ii,2) = hess(ii,2) + d1*dxy1 + d2*dxy2
            hess(ii,3) = hess(ii,3) + d1*dyy1 + d2*dyy2 
         enddo
 1000 continue         
      enddo
      

      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine r2d_directcdp_vec(nd,sources,ns,charge,
     $           dipstr,dipvec,targ,pot,thresh)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c
c     pot(ii)  = \sum_j log(|z_{j}|)*charge(ii,j)+
c                   dipstr(ii,j)* dipstr(ii,:,j) \cdot
c                      \nabla_src \log |targ(:)-sources(:,j)|
c
c     The potential is not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     dipstr(nd,ns) :   dipole strengths
c     dipvec(nd,2,ns) :   dipole orientations
c     targ          :   location of the target
c     thresh        :   threshold for computing potential
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)   (real *8)      : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(2,ns),targ(2),xdiff,ydiff,rr,r
      real *8 thresh,rtmp,thresh2
      real *8 pot(nd)
      real *8 dipvec(nd,2,ns)
      real *8 charge(nd,ns),dipstr(nd,ns)
      real *8 xdiff2,ydiff2,p1,p2

      thresh2 = thresh*thresh
c
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.le.thresh2) goto 1000

         rtmp = log(rr)/2.0d0

         p1 = -xdiff/rr
         p2 = -ydiff/rr         

         do ii=1,nd
            pot(ii) = pot(ii) + rtmp*charge(ii,i)
            pot(ii) = pot(ii) + dipstr(ii,i)*(dipvec(ii,1,i)*p1
     1           + dipvec(ii,2,i)*p2)
         enddo
 1000    continue         
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine r2d_directcdg_vec(nd,sources,ns,charge,dipstr,dipvec,
     1     targ,pot,grad,thresh)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c     pot(ii)  = \sum_j log(|z_{j}|)*charge(ii,j)+
c                   dipstr(ii,j)* dipstr(ii,:,j) \cdot
c                      \nabla_src \log |targ(:)-sources(:,j)|
c
c     grad  = \nabla(pot)
c
c     The potential and gradient are not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     dipstr(nd,ns) :   dipole strengths
c     dipvec(nd,2,ns) :   dipole orientations
c     targ          :   location of the target
c     thresh        :   threshold for computing potential and gradient
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (real *8)      : potential is incremented
c     grad(nd,2)  (real *8)      : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd
      real *8 sources(2,ns),targ(2)
      real *8 xdiff,ydiff,rr,r,thresh,rtmp,thresh2
      real *8 pot(nd),grad(nd,2)
      real *8 charge(nd,ns),dipstr(nd,ns)
      real *8 d1,d2
      real *8 dipvec(nd,2,ns)
      real *8 xdiff2,ydiff2,p1,p2,dxx1,dxx2,dxy1,dxy2,dyy1,dyy2
      real *8 rr2,rr3,dx,dy,dx1,dx2,dy1,dy2

      thresh2 = thresh*thresh
c
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         xdiff2 = xdiff*xdiff
         ydiff2 = ydiff*ydiff
         rr = xdiff2 + ydiff2
         if(rr.lt.thresh2) goto 1000
         rtmp = log(rr)/2

         rr2 = rr*rr
         
         dx = xdiff/rr
         dy = ydiff/rr

         p1 = -dx
         p2 = -dy

         dx1 = -(rr-2*xdiff2)/rr2
         dx2 = 2*xdiff*ydiff/rr2
         dy1 = dx2
         dy2 = -(rr-2*ydiff2)/rr2
         
c
         do ii = 1,nd
            pot(ii) = pot(ii) + rtmp*charge(ii,i)
            grad(ii,1) = grad(ii,1) + dx*charge(ii,i)
            grad(ii,2) = grad(ii,2) + dy*charge(ii,i) 
            
            d1 = dipstr(ii,i)*dipvec(ii,1,i)
            d2 = dipstr(ii,i)*dipvec(ii,2,i)

            pot(ii) = pot(ii) + d1*p1 + d2*p2
            grad(ii,1) = grad(ii,1) + d1*dx1 + d2*dx2
            grad(ii,2) = grad(ii,2) + d1*dy1 + d2*dy2 
         enddo
 1000 continue         
      enddo
      
      return
      end
c
c
c
c
c
c
c
c
c**********************************************************************
      subroutine r2d_directcdh_vec(nd,sources,ns,charge,dipstr,dipvec,
     1     targ,pot,grad,hess,thresh)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target point TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c
c     pot(ii)  = \sum_j log(|z_{j}|)*charge(ii,j)+
c                   dipstr(ii,j)* dipstr(ii,:,j) \cdot
c                      \nabla_src \log |targ(:)-sources(:,j)|
c
c     grad  = \nabla(pot)
c     hess = \nabla\nabla (pot)    (dxx,dxy,dyy)
c
c     The potential and gradient are not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     dipstr(nd,ns) :   dipole strengths
c     targ          :   location of the target
c     thresh        :   threshold for computing potential,gradient,and
c                       hessian
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (real *8)      : potential is incremented
c     grad(nd,2)  (real *8)      : gradient is incremented
c     hess(nd,3)  (real *8)      : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd
      real *8 sources(2,ns),targ(2)
      real *8 xdiff,ydiff,rr,r,thresh,thresh2
      real *8 pot(nd),grad(nd,2),hess(nd,3)
      real *8 rtmp
      real *8 charge(nd,ns),dipstr(nd,ns),d1,d2
      real *8 dipvec(nd,2,ns),dx1,dx2,dy1,dy2
      real *8 xdiff2,ydiff2,p1,p2,dxx1,dxx2,dxy1,dxy2,dyy1,dyy2
      real *8 rr2,rr3,dxx,dxy,dyy,dx,dy

      thresh2 = thresh*thresh

c
      do i = 1,ns
         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         xdiff2 = xdiff*xdiff
         ydiff2 = ydiff*ydiff
         rr = xdiff2 + ydiff2
         if(rr.lt.thresh2) goto 1000
         rtmp = log(rr)/2

         rr2 = rr*rr
         rr3 = rr2*rr
         
         dx = xdiff/rr
         dy = ydiff/rr
         dxx = (rr - 2*xdiff2)/rr2
         dxy = -2*xdiff*ydiff/rr2
         dyy = (rr - 2*ydiff2)/rr2

         p1 = -dx
         p2 = -dy

         dx1 = -dxx
         dx2 = -dxy
         dy1 = dx2
         dy2 = -dyy

         dxx1 = -2*xdiff*(xdiff2-3*ydiff2)/rr3
         dxx2 = 2*ydiff*(ydiff2-3*xdiff2)/rr3
         dxy1 = dxx2
         dxy2 = -dxx1
         dyy1 = dxy2
         dyy2 = -dxx2

         do ii = 1,nd
            
            pot(ii) = pot(ii) + rtmp*charge(ii,i)
            grad(ii,1) = grad(ii,1) + dx*charge(ii,i)
            grad(ii,2) = grad(ii,2) + dy*charge(ii,i)            
            hess(ii,1) = hess(ii,1) + dxx*charge(ii,i)
            hess(ii,2) = hess(ii,2) + dxy*charge(ii,i)
            hess(ii,3) = hess(ii,3) + dyy*charge(ii,i)

            d1 = dipstr(ii,i)*dipvec(ii,1,i)
            d2 = dipstr(ii,i)*dipvec(ii,2,i)

            pot(ii) = pot(ii) + d1*p1 + d2*p2
            grad(ii,1) = grad(ii,1) + d1*dx1 + d2*dx2
            grad(ii,2) = grad(ii,2) + d1*dy1 + d2*dy2
            hess(ii,1) = hess(ii,1) + d1*dxx1 + d2*dxx2
            hess(ii,2) = hess(ii,2) + d1*dxy1 + d2*dxy2
            hess(ii,3) = hess(ii,3) + d1*dyy1 + d2*dyy2 
         enddo
 1000 continue         
      enddo
      

      return
      end
c
c
c
