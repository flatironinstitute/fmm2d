c
c**********************************************************************
      subroutine c2d_directcp(nd,sources,ns,charge,
     $           targ,nt,pot,thresh)
cf2py  intent(in) nd
cf2py  intent(in) sources,ns,charge,targ,nt,thresh
cf2py  intent(out) pot
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target points TARGET, due to a vector of 
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
c     targ(2,nt)    :   location of the targets
c     nt            :   number of targets
c     thresh        :   threshold for computing potential
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)   (complex *16)      : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt),xdiff,ydiff,rr,r
      real *8 thresh,rtmp,thresh2
      complex *16 pot(nd,nt),z,ima4inv,ztmp
      complex *16 charge(nd,ns)

      thresh2 = thresh*thresh
c
      do j = 1,nt
         do i = 1,ns
            xdiff=targ(1,j)-sources(1,i)
            ydiff=targ(2,j)-sources(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.le.thresh2) goto 1000

            rtmp = log(rr)/2 
            do ii=1,nd
               pot(ii,j) = pot(ii,j) + rtmp*charge(ii,i)
            enddo
 1000       continue         
         enddo
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine c2d_directcg(nd,sources,ns,charge,targ,nt,pot,
     1             grad,thresh)
cf2py  intent(in) nd
cf2py  intent(in) sources,ns,charge,targ,nt,thresh
cf2py  intent(out) pot,grad
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target points TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c
c     pot(ii)  = \sum_j log(|targ-source(*,j)|)*charge(ii,j)
c
c     grad  = d/dz(pot)
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
c     targ(2,nt)    :   location of the targets
c     nt            :   number of targets
c     thresh        :   threshold for computing potential and gradient
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)     (complex *16)      : potential is incremented
c     grad(nd,nt)  (complex *16)      : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt)
      real *8 xdiff,ydiff,rr,r,thresh,rtmp,thresh2
      complex *16 pot(nd,nt),grad(nd,nt)
      complex *16 z,ztmp,zinv
      complex *16 charge(nd,ns)

      thresh2 = thresh*thresh
c
      do j = 1,nt
         do i = 1,ns
            xdiff=targ(1,j)-sources(1,i)
            ydiff=targ(2,j)-sources(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
            rtmp = log(rr)/2
            z=dcmplx(xdiff,ydiff)
            zinv = 1/z
c
            do ii = 1,nd
               pot(ii,j) = pot(ii,j) + rtmp*charge(ii,i)
               grad(ii,j) = grad(ii,j) + zinv*charge(ii,i)
            enddo
 1000       continue         
         enddo
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
      subroutine c2d_directch(nd,sources,ns,charge,targ,nt,
     1           pot,grad,hess,thresh)
cf2py  intent(in) nd
cf2py  intent(in) sources,ns,charge,targ,nt,thresh
cf2py  intent(out) pot,grad,hess
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target points TARGET, due to a vector of 
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
c     targ(2,nt)    :   location of the targets
c     nt            :   number of targets
c     thresh        :   threshold for computing potential,gradient,and
c                       hessian
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)     (complex *16)      : potential is incremented
c     grad(nd,nt)  (complex *16)      : gradient is incremented
c     hess(nd,nt)  (complex *16)      : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt)
      real *8 xdiff,ydiff,rr,r,thresh,thresh2
      complex *16 pot(nd,nt),grad(nd,nt),hess(nd,nt)
      real *8 rtmp
      complex *16 z,zinv,ztmp
      complex *16 charge(nd,ns)

      thresh2 = thresh*thresh

c
      do j = 1,nt
         do i = 1,ns
            xdiff=targ(1,j)-sources(1,i)
            ydiff=targ(2,j)-sources(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
            rtmp = log(rr)/2
            z=dcmplx(xdiff,ydiff)
            zinv = 1/z
            ztmp = -zinv*zinv
c
            do ii = 1,nd
               pot(ii,j) = pot(ii,j) + rtmp*charge(ii,i)
               grad(ii,j) = grad(ii,j) + zinv*charge(ii,i)
               hess(ii,j) = hess(ii,j) + ztmp*charge(ii,i)
            enddo
 1000       continue         
         enddo
      enddo
      

      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine c2d_directdp(nd,sources,ns,dipstr,
     $           targ,nt,pot,thresh)
cf2py  intent(in) nd
cf2py  intent(in) sources,ns,dipstr,targ,nt,thresh
cf2py  intent(out) pot
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target points TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c      z_j = cmplx(targ - source(:,j))
c
c     pot(ii)  = \sum_j dipstr(ii,j)/z_j
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
c     targ(2,nt)    :   location of the targets
c     nt            :   number of targets
c     thresh        :   threshold for computing potential
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)   (complex *16)      : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt),xdiff,ydiff,rr,r
      real *8 thresh,rtmp
      complex *16 pot(nd,nt),z,zinv
      complex *16 dipstr(nd,ns)
c
      do j = 1,nt
         do i = 1,ns
            xdiff=targ(1,j)-sources(1,i)
            ydiff=targ(2,j)-sources(2,i)
            z = dcmplx(xdiff,ydiff)
            if(abs(z).le.thresh) goto 1000

            zinv = 1/z 
            do ii=1,nd
               pot(ii,j) = pot(ii,j) + zinv*dipstr(ii,i)
            enddo
 1000       continue         
         enddo
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine c2d_directdg(nd,sources,ns,dipstr,targ,nt,pot,
     1             grad,thresh)
cf2py  intent(in) nd
cf2py  intent(in) sources,ns,dipstr,targ,nt,thresh
cf2py  intent(out) pot,grad
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target points TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c      z_j = cmplx(targ - source(:,j))
c
c     pot(ii)  = \sum_j dipstr(ii,j)/z_j
c
c     grad  = d/dz(pot)
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
c     targ(2,nt)    :   location of the targets
c     nt            :   number of targets
c     thresh        :   threshold for computing potential and gradient
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)     (complex *16)    : potential is incremented
c     grad(nd,nt)  (complex *16)      : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt)
      real *8 xdiff,ydiff,rr,r,thresh,rtmp
      complex *16 pot(nd,nt),grad(nd,nt)
      complex *16 z,zinv,ztmp2
      complex *16 dipstr(nd,ns)
c
      do j = 1,nt
         do i = 1,ns
            xdiff=targ(1,j)-sources(1,i)
            ydiff=targ(2,j)-sources(2,i)
            z=dcmplx(xdiff,ydiff) 
            if(abs(z).lt.thresh) goto 1000
            zinv = 1/z
            ztmp2 = -zinv*zinv
c
            do ii = 1,nd
               pot(ii,j) = pot(ii,j) + zinv*dipstr(ii,i)
               grad(ii,j) = grad(ii,j) + ztmp2*dipstr(ii,i)
            enddo
 1000       continue         
         enddo
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
      subroutine c2d_directdh(nd,sources,ns,dipstr,targ,nt,
     1           pot,grad,hess,thresh)
cf2py  intent(in) nd
cf2py  intent(in) sources,ns,dipstr,targ,nt,thresh
cf2py  intent(out) pot,grad,hess
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target points TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c      z_j = cmplx(targ - source(:,j))
c
c     pot(ii)  = \sum_j dipstr(ii,j)/z_j
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
c     dipstr(nd,ns) :   charge strengths
c     targ(2,nt)    :   location of the targets
c     nt            :   number of targets
c     thresh        :   threshold for computing potential,gradient,and
c                       hessian
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)     (complex *16)      : potential is incremented
c     grad(nd,nt)  (complex *16)      : gradient is incremented
c     hess(nd,nt)  (complex *16)      : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt)
      real *8 xdiff,ydiff,rr,r,thresh
      complex *16 pot(nd,nt),grad(nd,nt),hess(nd,nt)
      real *8 rtmp
      complex *16 z,zinv,ztmp2,ztmp3
      complex *16 dipstr(nd,ns)

c
      do j = 1,nt
         do i = 1,ns
            xdiff=targ(1,j)-sources(1,i)
            ydiff=targ(2,j)-sources(2,i)
            z = dcmplx(xdiff,ydiff)
            if(abs(z).lt.thresh) goto 1000
            zinv = 1/z
            ztmp2 = -zinv*zinv
            ztmp3 = -2*ztmp2*zinv
c
            do ii = 1,nd
               pot(ii,j) = pot(ii,j) + zinv*dipstr(ii,i)
               grad(ii,j) = grad(ii,j) + ztmp2*dipstr(ii,i)
               hess(ii,j) = hess(ii,j) + ztmp3*dipstr(ii,i)
            enddo
 1000       continue         
         enddo
      enddo
      

      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine c2d_directcdp(nd,sources,ns,charge,
     $           dipstr,targ,nt,pot,thresh)
cf2py  intent(in) nd
cf2py  intent(in) sources,ns,charge,dipstr,targ,nt,thresh
cf2py  intent(out) pot
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target points TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c
c     z_j = cmplx(targ - source(:,j))
c     pot(ii)  = \sum_j log(|z_{j}|)*charge(ii,j)+
c           dipstr(ii,j)/z_{j}
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
c     targ(2,nt)    :   location of the targets
c     nt            :   number of targets
c     thresh        :   threshold for computing potential
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)   (complex *16)      : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt),xdiff,ydiff,rr,r
      real *8 thresh,rtmp,thresh2
      complex *16 pot(nd,nt),z,ztmp,zinv
      complex *16 charge(nd,ns),dipstr(nd,ns)

      thresh2 = thresh*thresh
c
      do j = 1,nt
         do i = 1,ns
            xdiff=targ(1,j)-sources(1,i)
            ydiff=targ(2,j)-sources(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.le.thresh2) goto 1000

            rtmp = log(rr)/2.0d0
            z = dcmplx(xdiff,ydiff)
            zinv = 1/z
            do ii=1,nd
               pot(ii,j) = pot(ii,j) + rtmp*charge(ii,i)
     1                               + zinv*dipstr(ii,i)
            enddo
 1000       continue         
         enddo
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine c2d_directcdg(nd,sources,ns,charge,dipstr,targ,nt,pot,
     1             grad,thresh)
cf2py  intent(in) nd
cf2py  intent(in) sources,ns,charge,dipstr,targ,nt,thresh
cf2py  intent(out) pot,grad
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target points TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c     pot(ii)  = \sum_j log(|z_{j}|)*charge(ii,j)+
c           dipstr(ii,j)/z_{j}
c
c     grad  = d/dz(pot)
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
c     dipstr(nd,ns) :   charge strengths
c     targ(2,nt)    :   location of the targets
c     nt            :   number of targets
c     thresh        :   threshold for computing potential and gradient
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)     (complex *16)      : potential is incremented
c     grad(nd,nt)  (complex *16)      : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt)
      real *8 xdiff,ydiff,rr,r,thresh,rtmp,thresh2
      complex *16 pot(nd,nt),grad(nd,nt)
      complex *16 z,ztmp,zinv,ztmp2
      complex *16 charge(nd,ns),dipstr(nd,ns)

      thresh2 = thresh*thresh
c
      do j = 1,nt
         do i = 1,ns
            xdiff=targ(1,j)-sources(1,i)
            ydiff=targ(2,j)-sources(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
            rtmp = log(rr)/2
            z=dcmplx(xdiff,ydiff)
            zinv = 1/z
            ztmp2 = -zinv*zinv
c
            do ii = 1,nd
               pot(ii,j) = pot(ii,j) + rtmp*charge(ii,i)
     1                               + zinv*dipstr(ii,i)
               grad(ii,j) = grad(ii,j) + zinv*charge(ii,i)
     1                                 + ztmp2*dipstr(ii,i)
            enddo
 1000       continue         
         enddo
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
      subroutine c2d_directcdh(nd,sources,ns,charge,dipstr,targ,nt,
     1           pot,grad,hess,thresh)
cf2py  intent(in) nd
cf2py  intent(in) sources,ns,charge,dipstr,targ,nt,thresh
cf2py  intent(out) pot,grad,hess
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target points TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     We use the unscaled version of log
c     response: i.e., log|z|
c     
c
c     pot(ii)  = \sum_j log(|z_{j}|)*charge(ii,j)+
c           dipstr(ii,j)/z_{j}
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
c     dipstr(nd,ns) :   dipole strengths
c     targ(2,nt)    :   location of the targets
c     nt            :   number of targets
c     thresh        :   threshold for computing potential,gradient,and
c                       hessian
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)     (complex *16)      : potential is incremented
c     grad(nd,nt)  (complex *16)      : gradient is incremented
c     hess(nd,nt)  (complex *16)      : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt)
      real *8 xdiff,ydiff,rr,r,thresh,thresh2
      complex *16 pot(nd,nt),grad(nd,nt),hess(nd,nt)
      real *8 rtmp
      complex *16 z,zinv,ztmp2,ztmp3
      complex *16 charge(nd,ns),dipstr(nd,ns)

      thresh2 = thresh*thresh

c
      do j = 1,nt
         do i = 1,ns
            xdiff=targ(1,j)-sources(1,i)
            ydiff=targ(2,j)-sources(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
            rtmp = log(rr)/2
            z=dcmplx(xdiff,ydiff)
            zinv = 1/z
            ztmp2 = -zinv*zinv
            ztmp3 = -2*ztmp2*zinv
c
            do ii = 1,nd
               pot(ii,j) = pot(ii,j) + rtmp*charge(ii,i)
     1                               + zinv*dipstr(ii,i)
               grad(ii,j) = grad(ii,j) + zinv*charge(ii,i)
     1                                 + ztmp2*dipstr(ii,i)
               hess(ii,j) = hess(ii,j) + ztmp2*charge(ii,i)
     1                                 + ztmp3*dipstr(ii,i)
            enddo
 1000       continue         
         enddo
      enddo
      

      return
      end
c
c
c
