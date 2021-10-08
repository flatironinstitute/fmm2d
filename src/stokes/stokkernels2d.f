c
c     This file contains the Stokes direct kernel evaluators in two dimensions
c
c**********************************************************************
c      
c     We take the following conventions for the Stokes kernels
c
c     For a source y and target x, let r_i = x_i-y_i
c     and let r = sqrt(r_1^2 + r_2^2)
c
c     The Stokeslet, G_{ij}, and its associated pressure tensor, P_j,
c     (without the 1/2pi scaling) are
c
c     G_{ij}(x,y) = (r_i r_j)/(2r^2) - delta_{ij}log(r)/(2)
c     P_j(x,y) = r_j/r^2
c
c     The (Type I) stresslet, T_{ijk}, and its associated pressure
c     tensor, PI_{jk}, (without the 1/2pi scaling) are
c     
c     T_{ijk}(x,y) = -2 r_i r_j r_k/ r^4
c     PI_{jk} = - delta_{jk}/r^2 + 2 r_j r_k/r^4      
c

      subroutine st2ddirectstokg(nd,sources,stoklet,ns,targ,nt,
     1     pot,pre,grad,thresh)
c
c     This subroutine evaluates the potential and gradient due
c     to a collection of Stokeslet and stresslet sources and adds
c     to existing quantities (see definitions at top of file).
c
c       pot(x) = pot(x) + sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
c
c       pre(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
c
c       grad(x) = Grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
c      
c     where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
c     stresslet charge, and nu^{(m)} is the stresslet orientation
c     (note that each of these is a 2 vector per source point y^{(m)}).
c     For x a source point, the self-interaction in the sum is omitted. 
c
c
c-----------------------------------------------------------------------
      
c     INPUT:
c     
c     nd in: integer
c        number of densities
c     
c     nsource in: integer  
c        number of sources
c
c     source  in: double precision (2,nsource)
c        source(k,j) is the kth component of the jth
c        source location
c
c     stoklet in: double precision (nd,2,nsource) 
c        Stokeslet charge strengths (sigma vectors above)
c
c     ntarg   in: integer  
c        number of targs 
c
c     targ    in: double precision (2,ntarg)
c        targ(k,j) is the kth component of the jth
c        targ location
c     
c     thresh in: double precision
c        threshold for updating potential,
c        potential at target won't be updated if
c        |t - s| <= thresh, where t is the target
c        location and, and s is the source location 
c
c-----------------------------------------------------------------------
c
c   OUTPUT:
c
c     pot out: double precision(nd,2,ntarg) 
c        velocity at the targets
c      
c     pre out: double precision(nd,ntarg)
c        pressure at the targets
c      
c     grad out: double precision(nd,2,2,ntarg) 
c        gradient of velocity at the targets
c        gradtarg(l,i,j,k) is the ith component of the
c        gradient of the jth component of the velocity
c        for the lth density at the kth target
c
c------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,stoklet,ns
cf2py intent(in) targ,nt,thresh
cf2py intent(out) pot,pre,grad

      integer nd, ns, nt, istress
      real *8 sources(2,ns),targ(2,nt)
      real *8 stoklet(nd,2,ns)
      real *8 pot(nd,2,nt),pre(nd,nt),grad(nd,2,2,nt)
      real *8 thresh
      
c     local

      real *8 zdiff(2), d1, d2, d3, tempx1, tempx2, tempx3
      real *8 pl, pl2,pv, dmu(2), dnu(2), temp, r, r2, r3, r4, r5
      real *8 dmunu,rtmp      
      real *8 threshsq

      integer i, j, idim, l

      threshsq = thresh**2


c     stokeslet contribution
      
      do i = 1,nt
         do j = 1,ns
            zdiff(1) = targ(1,i)-sources(1,j)
            zdiff(2) = targ(2,i)-sources(2,j)

            r2 = zdiff(1)**2 + zdiff(2)**2 
            if (r2 .lt. threshsq) goto 10

            rtmp = zdiff(1)**2 - zdiff(2)**2
            r4 = r2*r2
            
            do idim = 1,nd

               pot(idim,1,i) = pot(idim,1,i)-
     1            stoklet(idim,1,j)*log(r2)/(4)
               pot(idim,2,i) = pot(idim,2,i)-
     1             stoklet(idim,2,j)*log(r2)/(4)
               
               pl = (zdiff(1)*stoklet(idim,1,j) +
     1              zdiff(2)*stoklet(idim,2,j))

               pl2 = zdiff(1)*stoklet(idim,2,j) + 
     1               zdiff(2)*stoklet(idim,1,j)
               
               pot(idim,1,i) = pot(idim,1,i) + zdiff(1)*pl/r2/2
               pot(idim,2,i) = pot(idim,2,i) + zdiff(2)*pl/r2/2

               grad(idim,1,1,i) = grad(idim,1,1,i) - rtmp*pl/r4/2 
               grad(idim,2,2,i) = grad(idim,2,2,i) + rtmp*pl/r4/2
               grad(idim,2,1,i) = grad(idim,2,1,i) + 
     1           (pl2*rtmp - 
     2             stoklet(idim,1,j)*4*zdiff(2)*zdiff(1)**2)/2/r4

               grad(idim,1,2,i) = grad(idim,1,2,i) - 
     1           (pl2*rtmp + 
     2             stoklet(idim,2,j)*4*zdiff(1)*zdiff(2)**2)/2/r4

               pre(idim,i) = pre(idim,i) + pl/r2

               
            enddo
 10         continue
         enddo
      enddo

      
      return
      end
c
c
c
c
c

      subroutine st2ddirectstokstrsg(nd,sources,ifstoklet,
     1     stoklet,
     1     istress, strslet, strsvec, ns,targ,nt,pot,pre,
     2     grad,thresh)
c
c     This subroutine evaluates the potential and gradient due
c     to a collection of Stokeslet and stresslet sources and adds
c     to existing quantities (see definitions at top of file).
c
c       pot(x) = pot(x) + sum_m G_{ij}(x,y^{(m)}) sigma^{(m)}_j
c
c       pre(x) = sum_m P_j(x,y^m) sigma^{(m)}_j
c
c       grad(x) = Grad[sum_m G_{ij}(x,y^m) sigma^{(m)}_j
c      
c     where sigma^{(m)} is the Stokeslet charge, mu^{(m)} is the
c     stresslet charge, and nu^{(m)} is the stresslet orientation
c     (note that each of these is a 2 vector per source point y^{(m)}).
c     For x a source point, the self-interaction in the sum is omitted. 
c
c
c-----------------------------------------------------------------------
      
c     INPUT:
c     
c     nd in: integer
c        number of densities
c     
c     nsource in: integer  
c        number of sources
c
c     source  in: double precision (2,nsource)
c        source(k,j) is the kth component of the jth
c        source location
c
c     stoklet in: double precision (nd,2,nsource) 
c        Stokeslet charge strengths (sigma vectors above)
c
c     istress in: integer
c        stresslet computation flag
c           istress = 1   =>  include standard stresslet
c                               (type I)
c     
c           NOT YET IMPLEMENTED
c      
c           ifstress = 2   =>  include symmetric stresslet
c                                   (type II)
c           ifstress = 3   =>  include rotlet
c           ifstress = 4   =>  include Stokes doublet
c                      otherwise do not include
c
c     strslet  in: double precision (nd,2,nsource) 
c        stresslet strengths (mu vectors above)
c
c     strsvec  in: double precision (nd,2,nsource)   
c        stresslet orientations (nu vectors above)
c
c     ntarg   in: integer  
c        number of targs 
c
c     targ    in: double precision (2,ntarg)
c        targ(k,j) is the kth component of the jth
c        targ location
c     
c     thresh in: double precision
c        threshold for updating potential,
c        potential at target won't be updated if
c        |t - s| <= thresh, where t is the target
c        location and, and s is the source location 
c
c-----------------------------------------------------------------------
c
c   OUTPUT:
c
c     pot out: double precision(nd,2,ntarg) 
c        velocity at the targets
c      
c     pre out: double precision(nd,ntarg)
c        pressure at the targets
c      
c     grad out: double precision(nd,2,2,ntarg) 
c        gradient of velocity at the targets
c        gradtarg(l,i,j,k) is the ith component of the
c        gradient of the jth component of the velocity
c        for the lth density at the kth target
c
c------------------------------------------------------------------
      implicit none
cf2py intent(in) nd,sources,stoklet,ns
cf2py intent(in) targ,nt,thresh
cf2py intent(out) pot,pre,grad

      integer nd, ns, nt, istress,ifstoklet
      real *8 sources(2,ns),targ(2,nt)
      real *8 stoklet(nd,2,ns),strslet(nd,2,ns),strsvec(nd,2,ns)
      real *8 pot(nd,2,nt),pre(nd,nt),grad(nd,2,2,nt)
      real *8 thresh
      
c     local

      real *8 zdiff(2), d1, d2, d3, tempx1, tempx2, tempx3
      real *8 pl, pl2,pv, dmu(2), dnu(2), temp, r, r2, r3, r4, r6
      real *8 dmunu,rtmp      
      real *8 threshsq

      integer i, j, idim, l
      if(ifstoklet.eq.1) call st2ddirectstokg(nd,sources,stoklet,
     1     ns,targ,nt,pot,pre,grad,thresh)

      threshsq = thresh**2
      if(istress.ne.1) goto 100
      

c     stresslet contribution
      
      do i = 1,nt
         do j = 1,ns
            zdiff(1) = targ(1,i)-sources(1,j)
            zdiff(2) = targ(2,i)-sources(2,j)

            r2 = zdiff(1)**2 + zdiff(2)**2 
            if (r2 .lt. threshsq) goto 10

            rtmp = zdiff(1)**2 - zdiff(2)**2
            r4 = r2*r2
            r6 = r4*r2
            
            do idim = 1,nd

               
               pl = (zdiff(1)*strslet(idim,1,j) +
     1              zdiff(2)*strslet(idim,2,j))

               pl2 = zdiff(1)*strsvec(idim,1,j) + 
     1               zdiff(2)*strsvec(idim,2,j)
               
               pv = strslet(idim,1,j)*strsvec(idim,1,j) + 
     1             strslet(idim,2,j)*strsvec(idim,2,j)

               
               pot(idim,1,i) = pot(idim,1,i) - 2*zdiff(1)*pl*pl2/r4
               pot(idim,2,i) = pot(idim,2,i) - 2*zdiff(2)*pl*pl2/r4

               grad(idim,1,1,i) = grad(idim,1,1,i) - 2*pl*pl2/r4 - 
     1            2*zdiff(1)*strslet(idim,1,j)*pl2/r4 -
     2            2*zdiff(1)*strsvec(idim,1,j)*pl/r4 +
     3            8*zdiff(1)*zdiff(1)*pl*pl2/r6
               grad(idim,2,1,i) = grad(idim,2,1,i) - 
     1            2*zdiff(1)*strslet(idim,2,j)*pl2/r4 -
     2            2*zdiff(1)*strsvec(idim,2,j)*pl/r4 +
     3            8*zdiff(1)*zdiff(2)*pl*pl2/r6
               grad(idim,1,2,i) = grad(idim,1,2,i) - 
     1            2*zdiff(2)*strslet(idim,1,j)*pl2/r4 -
     2            2*zdiff(2)*strsvec(idim,1,j)*pl/r4 +
     3            8*zdiff(1)*zdiff(2)*pl*pl2/r6
               grad(idim,2,2,i) = grad(idim,2,2,i) - 2*pl*pl2/r4 - 
     1            2*zdiff(2)*strslet(idim,2,j)*pl2/r4 -
     2            2*zdiff(2)*strsvec(idim,2,j)*pl/r4 +
     3            8*zdiff(2)*zdiff(2)*pl*pl2/r6

               pre(idim,i) = pre(idim,i) - 4*pl*pl2/r4 + 2*pv/r2 

               
            enddo
 10         continue
         enddo
      enddo
 100  continue
      
      return
      end

