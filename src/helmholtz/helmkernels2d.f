c
c**********************************************************************
      subroutine h2d_directcp(nd,wavek,sources,ns,charge,
     $           targ,nt,pot,thresh)
cf2py  intent(in) wavek,nd
cf2py  intent(in) ns,sources,charge,targ,nt,thresh
cf2py  intent(out) pot
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target points TARGETs, due to a vector of 
c     charges at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c
c     pot(ii)  = \sum_j H_0(k*|targ-source(*,j)|)*(eye/4)*charge(ii,j)
c
c     The potential is not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     |zk|*R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     wavek         :   Helmholtz parameter
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ(2,nt)    :   location of the targets
c     thresh        :   threshold for computing potential
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)   (complex *16)      : potential is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt),xdiff,ydiff,rr,r
      real *8 thresh
      complex *16 wavek,pot(nd,nt),z,ima4inv,ima,h0,h1
      complex *16 charge(nd,ns)
      data ima/(0.0d0,1.0d0)/
c
      ima4inv=ima/4
c
      do j=1,nt
        do i = 1,ns
           xdiff=targ(1,j)-sources(1,i)
           ydiff=targ(2,j)-sources(2,i)
           rr=xdiff*xdiff+ydiff*ydiff
           r=sqrt(rr)
           z=wavek*r
           if(abs(z).le.thresh) goto 1000
           ifexpon = 1
           call hank103(z, h0, h1, ifexpon)
           do ii=1,nd
              pot(ii,j) = pot(ii,j) + h0*charge(ii,i)*ima4inv
           enddo
 1000      continue         
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
      subroutine h2d_directcg(nd,wavek,sources,ns,charge,targ,nt,pot,
     1             grad,thresh)
cf2py  intent(in) wavek,nd
cf2py  intent(in) ns,sources,charge,targ,nt,thresh
cf2py  intent(out) pot,grad
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT and gradients GRAD
c     at the target points TARGETs, due to a vector of 
c     charges at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c
c     pot(ii)  = \sum_j H_0(k*|targ-source(*,j)|)*(eye/4)*charge(ii,j)
c
c     grad  = gradient(pot)
c
c     The potential and gradient are not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     |zk|*R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     wavek         :   Helmholtz parameter
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ(2,nt)    :   location of the targets
c     thresh        :   threshold for computing potential and gradient
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)     (complex *16)      : potential is incremented
c     grad(nd,2,nt)  (complex *16)      : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd,nt,j
      real *8 sources(2,ns),targ(2,nt)
      real *8 xdiff,ydiff,rr,r,thresh
      complex *16 wavek,pot(nd,nt),grad(nd,2,nt)
      complex *16 h0,h1,cd,z,ima,ima4inv,cdd
      complex *16 charge(nd,ns)
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      ima4inv=ima/4
c
      do j=1,nt
        do i = 1,ns
          xdiff=targ(1,j)-sources(1,i)
          ydiff=targ(2,j)-sources(2,i)
          rr=xdiff*xdiff+ydiff*ydiff
          r=sqrt(rr)
          z=wavek*r
          if(abs(z).lt.thresh) goto 1000
c
          ifexpon = 1
          call hank103(z, h0, h1, ifexpon)
          cdd = -h1*(wavek*ima4inv/r)
c
          do ii = 1,nd
            pot(ii,j) = pot(ii,j) + h0*charge(ii,i)*ima4inv
            cd = cdd*charge(ii,i)
            grad(ii,1,j) = grad(ii,1,j) + cd*xdiff
            grad(ii,2,j) = grad(ii,2,j) + cd*ydiff
         enddo
 1000 continue         
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
      subroutine h2d_directch(nd,wavek,sources,ns,charge,targ,nt,
     1           pot,grad,hess,thresh)
cf2py  intent(in) wavek,nd
cf2py  intent(in) ns,sources,charge,nt,targ,thresh
cf2py  intent(out) pot,grad,hess
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT, gradients GRAD
c     and Hessians at the target points TARGET, due to a vector of 
c     charges at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c
c     pot(ii)  = \sum_j H_0(k*|targ-source(*,j)|)*(eye/4)*charge(ii,j)
c
c     grad  = gradient(pot)
c     hess = Hessian(pot)
c
c     The potential,gradient and hessian are not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     |zk|*R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     wavek         :   Helmholtz parameter
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     targ(2,nt)    :   location of the target
c     nt            :   number of targets
c     thresh        :   threshold for computing potential,gradient,and
c                       hessian
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)     (complex *16)      : potential is incremented
c     grad(nd,2,nt)  (complex *16)      : gradient is incremented
c     hess(nd,3,nt)  (complex *16)      : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns,ifexpon,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt)
      real *8 xdiff,ydiff,rr,r,thresh
      complex *16 wavek,pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)
      complex *16 h0,h1,h2z,cd,z,ima,ima4inv,cdd,cdd2,hf1,hf2,hf3
      complex *16 charge(nd,ns)
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      ima4inv=ima/4
c
      do j = 1,nt
        do i = 1,ns
          xdiff=targ(1,j)-sources(1,i)
          ydiff=targ(2,j)-sources(2,i)
          rr=xdiff*xdiff+ydiff*ydiff
          r=sqrt(rr)
          z=wavek*r
          if(abs(z).lt.thresh) goto 1000
c 
          ifexpon = 1
          call hank103(z, h0, h1, ifexpon)
          cdd = -h1*(wavek*ima4inv/r)
          cdd2 = (wavek*ima4inv/r)/rr
          h2z=(-z*h0+2*h1)
          hf1 = (h2z*xdiff*xdiff-rr*h1)
          hf2 = (h2z*xdiff*ydiff      )
          hf3 = (h2z*ydiff*ydiff-rr*h1)
c
          do ii = 1,nd
            pot(ii,j) = pot(ii,j) + h0*charge(ii,i)*ima4inv
            cd = cdd*charge(ii,i)
            grad(ii,1,j) = grad(ii,1,j) + cd*xdiff
            grad(ii,2,j) = grad(ii,2,j) + cd*ydiff
c
            cd = cdd2*charge(ii,i)
            hess(ii,1,j) = hess(ii,1,j) + cd*hf1
            hess(ii,2,j) = hess(ii,2,j) + cd*hf2
            hess(ii,3,j) = hess(ii,3,j) + cd*hf3
          enddo
 1000 continue         
        enddo
      enddo

      return
      end
c
c
c
c
c**********************************************************************
      subroutine h2d_directdp(nd,wavek,sources,ns,dipstr,dipvec,
     $           targ,nt,pot,thresh)
cf2py  intent(in) wavek,nd
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ,thresh
cf2py  intent(out) pot
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials POT
c     at the target point TARGET, due to a collection of 
c     dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c     pot(ii)  = \sum_j [\grad(H_0(k*|targ-source(*,j)|)) \cdot 
c                       dipvec(ii,*,j)]
c                      *(eye/4)*dipstr(ii,j)
c
c     The potential is not computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     |zk|*R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     wavek         :   Helmholtz parameter
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     dipstr(nd,ns) :   dipole strengths
c     dipvec(nd,2,ns) :   dipole vectors
c     targ(2,nt)    :   location of the target
c     nt            :   number of targets
c     thresh        :   threshold for computing potential
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)   (complex *16)      : calculated potential
c---------------------------------------------------------------------
      integer i,ns, ifexpon,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt),dipvec(nd,2,ns)
      real *8 xdiff,ydiff,rr,r,ctheta,stheta,dotprod
      real *8 thresh
      complex *16 wavek,pot(nd,nt)
      complex *16 h0,h1,h2,h3,cd,z,ima,ima4inv,cdd
      complex *16 dipstr(nd,ns)
      data ima/(0.0d0,1.0d0)/
c
      ima4inv=ima/4
c
      do j = 1,nt
        do i = 1,ns
          xdiff=targ(1,j)-sources(1,i)
          ydiff=targ(2,j)-sources(2,i)
          rr=xdiff*xdiff+ydiff*ydiff
          r=sqrt(rr)
          z=wavek*r
          if(abs(z).lt.thresh) goto 1000
c
          ifexpon = 1
          call hank103(z, h0, h1, ifexpon)
          cdd = h1/r*wavek*ima4inv
c
          ctheta=xdiff/r
          stheta=ydiff/r
          h2 = 2*h1/z-h0
          h3 = 4*h2/z-h1
          do ii = 1,nd
            cd = cdd*dipstr(ii,i)
            dotprod = xdiff*dipvec(ii,1,i)+ydiff*dipvec(ii,2,i)
            pot(ii,j) = pot(ii,j)+cd*dotprod
          enddo
 1000 continue         
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
c**********************************************************************
      subroutine h2d_directdg(nd,wavek,sources,ns,dipstr,dipvec,
     $           targ,nt,pot,grad,thresh)
cf2py  intent(in) wavek,nd,nt
cf2py  intent(in) ns,sources,dipstr,dipvec,targ,thresh
cf2py  intent(out) pot,grad
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potential POT and gradient GRAD
c     at the target point TARGET, due to a collection of 
c     dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c     
c     pot(ii)  = \sum_j [\grad(H_0(k*|targ-source(*,j)|)) 
c                      \cdot dipvec(ii,*,j)]
c                      *(eye/4)*dipstr(ii,j)
c
c     The potential and gradient are not computed 
c     if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     |zk|*R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     wavek         :   Helmholtz parameter
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     dipstr(nd,ns) :   dipole strengths
c     dipvec(nd,2,ns):   dipole vectors
c     targ(2,nt)    :   location of the target
c     nt            :   number of targets
c     thresh        :   threshold for computing potential and gradient
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)     (complex *16)   : potential is incremented
c     grad(nd,2)  (complex *16)   : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns, ifexpon,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt),dipvec(nd,2,ns)
      real *8 xdiff,ydiff,rr,r,ctheta,stheta,dotprod
      real *8 thresh
      complex *16 wavek,pot(nd,nt),grad(nd,2,nt)
      complex *16 h0,h1,h2,h3,cd,z,ima,ima4inv,hxx,hxy,hyy,cdd,cdd2
      complex *16 hf1,hf2,hf3
      complex *16 dipstr(nd,ns)
      data ima/(0.0d0,1.0d0)/
c
      ima4inv=ima/4

c
      do j = 1,nt
        do i = 1,ns
          xdiff=targ(1,j)-sources(1,i)
          ydiff=targ(2,j)-sources(2,i)
          rr=xdiff*xdiff+ydiff*ydiff
          r=sqrt(rr)
          z=wavek*r
          if(abs(z).lt.thresh) goto 1000
c
          ifexpon = 1
          call hank103(z, h0, h1, ifexpon)
c
          ctheta=xdiff/r
          stheta=ydiff/r
          h2 = 2*h1/z-h0
          h3 = 4*h2/z-h1
          cdd = h1/r*wavek*ima4inv
          cdd2 = -wavek**2*ima4inv
          hf1 = (h2*(ctheta*ctheta-0.5d0)-h0/2)
          hf2 = (h2*ctheta*stheta             )
          hf3 = (h2*(stheta*stheta-0.5d0)-h0/2)
          do ii = 1,nd
            cd = cdd*dipstr(ii,i)
            dotprod = xdiff*dipvec(ii,1,i)+ydiff*dipvec(ii,2,i)
            pot(ii,j) = pot(ii,j)+cd*dotprod
            cd = cdd2*dipstr(ii,i)
            hxx = cd*hf1
            hxy = cd*hf2
            hyy = cd*hf3
            grad(ii,1,j) = grad(ii,1,j) + hxx*dipvec(ii,1,i) + 
     1          hxy*dipvec(ii,2,i)
            grad(ii,2,j) = grad(ii,2,j) + hxy*dipvec(ii,1,i) + 
     1          hyy*dipvec(ii,2,i)
          enddo
 1000 continue         
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
c**********************************************************************
      subroutine h2d_directdh(nd,wavek,sources,ns,dipstr,dipvec,
     $           targ,nt,pot,grad,hess,thresh)
cf2py  intent(in) wavek,nd,nt
cf2py  intent(in) ns,sources,dipstr,dipvec,targ,thresh
cf2py  intent(out) pot,grad,hess
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potential POT, gradient GRAD and
c     Hessian HESS at the target point TARGET, due to a collection of 
c     dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c     
c     pot(ii)  = \sum_j [\grad(H_0(k*|targ-source(*,j)|)) 
c                     \cdot dipvec(ii,*,j)]
c                      *(eye/4)*dipstr(ii,j)
c
c     The potential,gradient, and hessian are not 
c     computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     |zk|*R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     wavek         :   Helmholtz parameter
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     dipstr(nd,ns) :   dipole strengths
c     dipvec(nd,2,ns) :   dipole vectors
c     targ          :   location of the target
c     thresh        :   threshold for computing potential, gradient,
c                        and hessian
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)     (complex *16)   : potential is incremented
c     grad(nd,2)  (complex *16)   : gradient is incremented
c     hess(nd,3)  (complex *16)   : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns, ifexpon,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt),dipvec(nd,2,ns)
      real *8 xdiff,ydiff,rr,r,ctheta,stheta,dotprod
      real *8 thresh
      complex *16 wavek,pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)
      complex *16 h0,h1,h2,h3,cd,z,ima,ima4inv,hxx,hxy,hyy
      complex *16 hxxx,hxxy,hxyy,hyyy,cdd,cdd2,hf1,hf2,hf3
      complex *16 hxxx1,hxxy1,hxyy1,hyyy1
      complex *16 dipstr(nd,ns)
      data ima/(0.0d0,1.0d0)/
c
      ima4inv=ima/4
c
      do j = 1,nt
        do i = 1,ns
          xdiff=targ(1,j)-sources(1,i)
          ydiff=targ(2,j)-sources(2,i)
          rr=xdiff*xdiff+ydiff*ydiff
          r=sqrt(rr)
          z=wavek*r
          if(abs(z).lt.thresh) goto 1000
c
          ifexpon = 1
          call hank103(z, h0, h1, ifexpon)
c
          ctheta=xdiff/r
          stheta=ydiff/r
          h2 = 2*h1/z-h0
          h3 = 4*h2/z-h1
          cdd = h1/r*wavek*ima4inv
          cdd2 = -wavek**2*ima4inv
          hf1 = (h2*(ctheta*ctheta-0.5d0)-h0/2)
          hf2 = (h2*ctheta*stheta             )
          hf3 = (h2*(stheta*stheta-0.5d0)-h0/2)
          hxxx1 = (-h1/2*(-1.5d0)
     $         -h3/2*(+ctheta*ctheta-0.5d0-stheta*stheta)
     $         )*ctheta
          hxxy1 = (-h1/2*(-0.5d0)
     $         -h3/2*(1.5d0*ctheta*ctheta-0.5d0*stheta*stheta)
     $         )*stheta
          hxyy1 = (-h1/2*(-0.5d0)
     $         -h3/2*(1.5d0*stheta*stheta-0.5d0*ctheta*ctheta)
     $         )*ctheta
          hyyy1 = (-h1/2*(-1.5d0)
     $         -h3/2*(+stheta*stheta-0.5d0-ctheta*ctheta)
     $         )*stheta

          do ii = 1,nd
            dotprod = xdiff*dipvec(ii,1,i)+ydiff*dipvec(ii,2,i)
            cd = cdd*dipstr(ii,i)
            pot(ii,j) = pot(ii,j)+cd*dotprod
            cd = cdd2*dipstr(ii,i)
            hxx = cd*hf1
            hxy = cd*hf2
            hyy = cd*hf3
            grad(ii,1,j) = grad(ii,1,j) + hxx*dipvec(ii,1,i) + 
     $         hxy*dipvec(ii,2,i)
            grad(ii,2,j) = grad(ii,2,j) + hxy*dipvec(ii,1,i) + 
     $         hyy*dipvec(ii,2,i)
c
            cd = -wavek**3*dipstr(ii,i)*ima4inv
            hxxx = cd*hxxx1
            hxxy = cd*hxxy1
            hxyy = cd*hxyy1
            hyyy = cd*hyyy1
            hess(ii,1,j) = hess(ii,1,j) + hxxx*dipvec(ii,1,i) + 
     $          hxxy*dipvec(ii,2,i)
            hess(ii,2,j) = hess(ii,2,j) + hxxy*dipvec(ii,1,i) + 
     $           hxyy*dipvec(ii,2,i)
            hess(ii,3,j) = hess(ii,3,j) + hxyy*dipvec(ii,1,i) + 
     $           hyyy*dipvec(ii,2,i)
          enddo
 1000 continue         
        enddo
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine h2d_directcdp(nd,wavek,sources,ns,
     $     charge,dipstr,dipvec,targ,nt,pot,thresh)
cf2py  intent(in) wavek,nd,nt
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,targ,thresh
cf2py  intent(out) pot
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potential POT
c     at the target point TARGET, due to a collection of 
c     charges and dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c     pot  = H_0(k*r)*(eye/4)
c
c     pot(ii)  = \sum_j H_0(k*|targ-source(*,j)|)*(eye/4)*charge(ii,j)
c     
c              + \sum_j [\grad(H_0(k*|targ-source(*,j)|)) 
c                    \cdot dipvec(ii,*,j)]
c                      *(eye/4)*dipstr(ii,j)
c
c     The potential is not 
c     computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     |zk|*R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     wavek         :   Helmholtz parameter
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     dipstr        :   dipole strengths
c     dipvec(nd,2,ns):   dipole vectors
c     targ          :   location of the target
c     thresh        :   threshold for computing potential
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd)   (complex *16)      : potential is incremented
c---------------------------------------------------------------------
      integer i,ns, ifexpon,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt),dipvec(nd,2,ns)
      real *8 xdiff,ydiff,rr,r,ctheta,stheta,dotprod
      real *8 thresh
      complex *16 wavek,pot(nd,nt)
      complex *16 h0,h1,h2,h3,cd,z,ima,ima4inv,cdd
      complex *16 dipstr(nd,ns),charge(nd,ns)
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      ima4inv=ima/4
c
      do j = 1,nt
        do i = 1,ns
c
          xdiff=targ(1,j)-sources(1,i)
          ydiff=targ(2,j)-sources(2,i)
          rr=xdiff*xdiff+ydiff*ydiff
          r=sqrt(rr)
          z=wavek*r
          if(abs(z).lt.thresh) goto 1000

          ifexpon = 1
          call hank103(z, h0, h1, ifexpon)
c
          ctheta=xdiff/r
          stheta=ydiff/r
          h2 = 2*h1/z-h0
          h3 = 4*h2/z-h1
          cdd = h1/r*wavek*ima4inv
          do ii = 1,nd
            cd = cdd*dipstr(ii,i)
            pot(ii,j) = pot(ii,j) + h0*charge(ii,i)*ima4inv
            dotprod = xdiff*dipvec(ii,1,i)+ydiff*dipvec(ii,2,i)
            pot(ii,j) = pot(ii,j) + cd*dotprod
          enddo
 1000 continue         
        enddo
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine h2d_directcdg(nd,wavek,sources,ns,
     $     charge,dipstr,dipvec,targ,nt,pot,grad,thresh)
cf2py  intent(in) wavek,nd,nt
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,targ,thresh
cf2py  intent(out) pot,grad
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potential POT
c     at the target point TARGET, due to a collection of 
c     charges and dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c     pot  = H_0(k*r)*(eye/4)
c
c     pot(ii)  = \sum_j H_0(k*|targ-source(*,j)|)*(eye/4)*charge(ii,j)
c     
c              + \sum_j [\grad(H_0(k*|targ-source(*,j)|)) 
c                    \cdot dipvec(ii,*,j)]
c                      *(eye/4)*dipstr(ii,j)
c
c
c     The potential and gradient are not 
c     computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     |zk|*R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     wavek         :   Helmholtz parameter
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     dipstr(nd,ns) :   dipole strengths
c     dipvec(nd,2,ns):   dipole vectors
c     targ(2,nt)    :   location of the target
c     nt            :   number of targets
c     thresh        :   threshold for computing potential and gradient
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)    (complex *16)      : potential is incremented
c     grad(nd,2,nt) (complex *16)      : gradient is incremented
c---------------------------------------------------------------------
      integer i,ns, ifexpon,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt),dipvec(nd,2,ns)
      real *8 xdiff,ydiff,rr,r,ctheta,stheta,dotprod
      real *8 thresh
      complex *16 wavek,pot(nd,nt),grad(nd,2,nt)
      complex *16 h0,h1,h2,h3,cd,z,ima,ima4inv,cdd,cdd2
      complex *16 hxx,hxy,hyy,hf1,hf2,hf3
      complex *16 dipstr(nd,ns),charge(nd,ns)
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      ima4inv=ima/4
c
      do j = 1,nt
        do i = 1,ns
c
          xdiff=targ(1,j)-sources(1,i)
          ydiff=targ(2,j)-sources(2,i)
          rr=xdiff*xdiff+ydiff*ydiff
          r=sqrt(rr)
          z=wavek*r
          if(abs(z).lt.thresh) goto 1000
          ifexpon = 1
          call hank103(z, h0, h1, ifexpon)
          cdd = -h1*(wavek*ima4inv/r)
          do ii = 1,nd
            pot(ii,j) = pot(ii,j) + h0*charge(ii,i)*ima4inv
            cd = cdd*charge(ii,i)
            grad(ii,1,j) = grad(ii,1,j) + cd*xdiff
            grad(ii,2,j) = grad(ii,2,j) + cd*ydiff
         enddo
c
         ctheta=xdiff/r
         stheta=ydiff/r
         h2 = 2.0d0*h1/z-h0
         h3 = 4.0d0*h2/z-h1
         cdd = h1/r*wavek*ima4inv
         cdd2 = -wavek**2*ima4inv
         hf1 = (h2*(ctheta*ctheta-0.5d0)-h0/2)
         hf2 = (h2*ctheta*stheta             )
         hf3 = (h2*(stheta*stheta-0.5d0)-h0/2)
         do ii = 1,nd
            cd = cdd*dipstr(ii,i)
            dotprod = xdiff*dipvec(ii,1,i)+ydiff*dipvec(ii,2,i)
            pot(ii,j) = pot(ii,j)+cd*dotprod
            cd = cdd2*dipstr(ii,i)
            hxx = cd*hf1
            hxy = cd*hf2
            hyy = cd*hf3
            grad(ii,1,j) = grad(ii,1,j) + hxx*dipvec(ii,1,i) + 
     $          hxy*dipvec(ii,2,i)
            grad(ii,2,j) = grad(ii,2,j) + hxy*dipvec(ii,1,i) + 
     $          hyy*dipvec(ii,2,i)
          enddo
 1000 continue         
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
c**********************************************************************
      subroutine h2d_directcdh(nd,wavek,sources,ns,
     $     charge,dipstr,dipvec,targ,nt,pot,grad,hess,thresh)
cf2py  intent(in) wavek,nd,nt
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,targ,thresh
cf2py  intent(out) pot,grad,hess
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potential POT
c     at the target point TARGET, due to a collection of 
c     charges and dipoles at SOURCE(2,ns). 
c     The scaling is that required of the delta function
c     response: i.e.,
c     
c     pot  = H_0(k*r)*(eye/4)
c
c     pot(ii)  = \sum_j H_0(k*|targ-source(*,j)|)*(eye/4)*charge(ii,j)
c     
c              + \sum_j [\grad(H_0(k*|targ-source(*,j)|)) 
c                 \cdot dipvec(ii,*,j)]
c                      *(eye/4)*dipstr(ii,j)
c
c     The potential, gradient and hessian are not 
c     computed if |r| < thresh
c     (Recommended value for threshold in an FMM is 
c     |zk|*R*eps, where R is the size of the computation
c     domain and eps is machine precision
c     
c---------------------------------------------------------------------
c     INPUT:
c
c     wavek         :   Helmholtz parameter
c     sources(2,ns) :   location of the sources
c     ns            :   number of sources
c     charge(nd,ns) :   charge strengths
c     dipstr(nd,ns) :   dipole strengths
c     dipvec(nd,2,ns) :   dipole vectors
c     targ(2,nt)    :   location of the target
c     nt            :   number of targets
c     thresh        :   threshold for computing potential, gradient,
c                       and hessian
c---------------------------------------------------------------------
c     OUTPUT:
c
c     pot(nd,nt)     (complex *16)      : potential is incremented
c     grad(nd,2,nt)  (complex *16)      : gradient is incremented
c     hess(nd,3,nt)  (complex *16)      : Hessian is incremented
c---------------------------------------------------------------------
      integer i,ns, ifexpon,ii,nd,j,nt
      real *8 sources(2,ns),targ(2,nt),dipvec(nd,2,ns)
      real *8 xdiff,ydiff,rr,r,ctheta,stheta,dotprod
      real *8 thresh
      complex *16 wavek,pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)
      complex *16 h0,h1,h2,h3,cd,z,ima,ima4inv,cdd,cdd2
      complex *16 hxx,hxy,hyy,h2z
      complex *16 hxxx,hxxy,hxyy,hyyy
      complex *16 hxxx1,hxxy1,hxyy1,hyyy1,hf1,hf2,hf3
      complex *16 dipstr(nd,ns),charge(nd,ns)
      data ima/(0.0d0,1.0d0)/
c
c ... Calculate offsets and distance
c
      ima4inv=ima/4
c
      do j = 1,nt
        do i = 1,ns
c
          xdiff=targ(1,j)-sources(1,i)
          ydiff=targ(2,j)-sources(2,i)
          rr=xdiff*xdiff+ydiff*ydiff
          r=sqrt(rr)
          z=wavek*r
          if(abs(z).lt.thresh) goto 1000
          ifexpon = 1
          call hank103(z, h0, h1, ifexpon)
          do ii = 1,nd
            pot(ii,j) = pot(ii,j) + h0*charge(ii,i)*ima4inv
            cd = -h1*(wavek*charge(ii,i)*ima4inv/r)
            grad(ii,1,j) = grad(ii,1,j) + cd*xdiff
            grad(ii,2,j) = grad(ii,2,j) + cd*ydiff
            cd = (wavek*charge(ii,i)*ima4inv/r)/rr
            h2z=(-z*h0+2*h1)
            hess(ii,1,j) = hess(ii,1,j) + cd*(h2z*xdiff*xdiff-rr*h1)
            hess(ii,2,j) = hess(ii,2,j) + cd*(h2z*xdiff*ydiff      )
            hess(ii,3,j) = hess(ii,3,j) + cd*(h2z*ydiff*ydiff-rr*h1)
          enddo

          ctheta=xdiff/r
          stheta=ydiff/r
          h2 = 2*h1/z-h0
          h3 = 4*h2/z-h1
          cdd = h1/r*wavek*ima4inv
          cdd2 = -wavek**2*ima4inv
          hf1 = (h2*(ctheta*ctheta-0.5d0)-h0/2)
          hf2 = (h2*ctheta*stheta             )
          hf3 = (h2*(stheta*stheta-0.5d0)-h0/2)
          hxxx1 = (-h1/2*(-1.5d0)
     $         -h3/2*(+ctheta*ctheta-0.5d0-stheta*stheta)
     $         )*ctheta
          hxxy1 = (-h1/2*(-0.5d0)
     $         -h3/2*(1.5d0*ctheta*ctheta-0.5d0*stheta*stheta)
     $         )*stheta
          hxyy1 = (-h1/2*(-0.5d0)
     $         -h3/2*(1.5d0*stheta*stheta-0.5d0*ctheta*ctheta)
     $         )*ctheta
          hyyy1 = (-h1/2*(-1.5d0)
     $         -h3/2*(+stheta*stheta-0.5d0-ctheta*ctheta)
     $         )*stheta
c
          do ii = 1,nd
            cd = cdd*dipstr(ii,i)
            dotprod = xdiff*dipvec(ii,1,i)+ydiff*dipvec(ii,2,i)
            pot(ii,j) = pot(ii,j)+cd*dotprod
            cd = cdd2*dipstr(ii,i)
            hxx = cd*hf1
            hxy = cd*hf2
            hyy = cd*hf3
            grad(ii,1,j) = grad(ii,1,j) + hxx*dipvec(ii,1,i) + 
     $         hxy*dipvec(ii,2,i)
            grad(ii,2,j) = grad(ii,2,j) + hxy*dipvec(ii,1,i) + 
     $         hyy*dipvec(ii,2,i)
            cd = -wavek**3*dipstr(ii,i)*ima4inv
            hxxx = cd*hxxx1
            hxxy = cd*hxxy1
            hxyy = cd*hxyy1
            hyyy = cd*hyyy1
            hess(ii,1,j) = hess(ii,1,j) + hxxx*dipvec(ii,1,i) +
     $          hxxy*dipvec(ii,2,i)
            hess(ii,2,j) = hess(ii,2,j) + hxxy*dipvec(ii,1,i) + 
     $          hxyy*dipvec(ii,2,i)
            hess(ii,3,j) = hess(ii,3,j) + hxyy*dipvec(ii,1,i) + 
     $           hyyy*dipvec(ii,2,i)
          enddo
 1000 continue         
        enddo
      enddo
      return
      end
c
c
