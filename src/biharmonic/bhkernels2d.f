c
c      In this file, the subroutines compute the potential, and
c      gradients due to a collection of biharmonic charges and dipoles  
c
c      The complex velocity at zt = xt + i yt 
c      due to a charge a complex charge c1,c2 at zs = xs+i ys is given by
c
c
c      vel = 2*c1 log|zt-zs| + c2 (zt-zs)/(zt_bar - zs_bar)
c            
c      The complex velocity due to a dipole with complex parameters
c      d1,d2,d3 is given by
c
c
c      vel = d1/(zt-zs) + d2 (zt-zs)/(zt_bar-zs_bar)^2 + 
c             d3/(zt_bar - zs_bar)
c
c      Analytic component of the gradient grada_z = d/dz(vel(z))_proj(z)
c      Analytic component of the gradient grada_zbar = d/dz(vel(z))_proj(zbar)
c
c      Anti analytic component of the gradient (gradaa) = d/dzbar (vel(z))) 
c
c
c

      subroutine bh2d_directcp(nd,sources,ns,charges,targ,nt,vel,
     1         thresh)
c
cf2py intent(in) nd,sources,ns,charges,targ,thresh,nt
cf2py intent(out) vel
c********************************************************************
c      This subroutine INCREMENTS the complex velocity VEL 
c      at the targets at targ(2) due to a collection of ns
c      charges at sources(2,ns)
c----------------------------------------------------------------------
c      INPUT parameters
c      nd           : number of densities
c      sources(2,ns): location of the sources
c      ns           : number of sources
c      charges(2,ns): charge strength
c      targ         : target location
c      nt           : number of targets
c--------------------------------------------------------------------
c      OUTPUT
c      vel          : Complex velocity at target
c
c--------------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(2,ns), targ(2,nt)
      complex *16 vel(nd,nt),zs,zt,zdis
      complex *16 charges(nd,2,ns)
      complex *16 zdis1,zdis2
      complex *16 eye

      eye = dcmplx(0,1.0d0)
      do j=1,nt
        zt = dcmplx(targ(1,j),targ(2,j))
        do i=1,ns
           zs = dcmplx(sources(1,i),sources(2,i))
           zdis = zt-zs
           if(abs(zdis).le.thresh) goto 1111
           zdis1 = 1.0d0/zdis
           do idim=1,nd
             vel(idim,j)=vel(idim,j)+
     1              2*charges(idim,1,i)*log(cdabs(zdis))+
     1              charges(idim,2,i)*dconjg(zdis1)*zdis
           enddo
 1111 continue
         enddo
      enddo
   
      return
      end
c
c
c
c
c

      subroutine bh2d_directcg(nd,sources,ns,charges,
     1         targ,nt,vel,grad,thresh)
cf2py intent(in) nd,sources,ns,charges,targ,thresh,nt
cf2py intent(out) vel,grad
c********************************************************************
c      This subroutine INCREMENTS the complex velocity VEL and its
c      gradients GRADA_Z, GRADA_ZBAR, GRADAA, at the
c      targets at targ(2) due to a collection of ns
c      charges at sources(2,ns)
c      
c----------------------------------------------------------------------
c      INPUT parameters
c      nd           : number of densities
c      sources(2,ns): location of the sources
c      ns           : number of sources
c      charges(2,ns): charge strength
c      targ         : target location
c      nt           : number of targets
c--------------------------------------------------------------------
c      OUTPUT
c      vel          : Complex velocity at target
c      grad         : Complex gradient 
c
c--------------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(2,ns), targ(2,nt)
      complex *16 vel(nd,nt),zs,zt,zdis
      complex *16 charges(nd,2,ns)
      complex *16 zdis1,zdis2
      complex *16 grad(nd,3,nt),eye

      eye = dcmplx(0,1.0d0)
      do j=1,nt
        zt = dcmplx(targ(1,j),targ(2,j))
        do i=1,ns
           zs = dcmplx(sources(1,i),sources(2,i))
           zdis = zt-zs
           if(abs(zdis).le.thresh) goto 1111
           zdis1 = 1.0d0/zdis
           zdis2=zdis1**2
           do idim=1,nd
             vel(idim,j)=vel(idim,j)+
     1               2*charges(idim,1,i)*log(cdabs(zdis))+
     1               charges(idim,2,i)*dconjg(zdis1)*zdis
     
             grad(idim,1,j)=grad(idim,1,j) + charges(idim,1,i)*zdis1
             
             grad(idim,2,j)=grad(idim,2,j) +
     1           charges(idim,2,i)*dconjg(zdis1)

             grad(idim,3,j)=grad(idim,3,j) + 
     1           charges(idim,1,i)*dconjg(zdis1)
             grad(idim,3,j)=grad(idim,3,j) - 
     1          charges(idim,2,i)*dconjg(zdis2)*zdis
          enddo
 1111 continue
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

      subroutine bh2d_directdp(nd,sources,ns,dippar,
     1         targ,nt,vel,thresh)
cf2py intent(in) nd,sources,ns,dippar,targ,thresh,nt
cf2py intent(out) vel
c********************************************************************
c      This subroutine INCREMENTS the complex velocity VEL 
c      at the
c      targets at targ(2) due to a collection of ns
c      dipoles at sources(2,ns)
c----------------------------------------------------------------------
c      INPUT parameters
c      nd           : number of densities
c      sources(2,ns): location of the sources
c      ns           : number of sources
c      dippar(3,ns) : dipole parameters 
c      targ         : target location
c      nt           : number of targets
c--------------------------------------------------------------------
c      OUTPUT
c      vel          : Complex velocity at target
c
c--------------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns,nt
      real *8 sources(2,ns), targ(2,nt)
      complex *16 vel(nd,nt),zs,zt,zdis
      complex *16 dippar(nd,3,ns)
      complex *16 zdis1,zdis2
      complex *16 eye

      eye = dcmplx(0,1.0d0)
      do j=1,nt
        zt = dcmplx(targ(1,j),targ(2,j))
        do i=1,ns
           zs = dcmplx(sources(1,i),sources(2,i))
           zdis = zt-zs
           if(abs(zdis).le.thresh) goto 1111
           zdis1 = 1.0d0/zdis
           zdis2=zdis1**2
           do idim=1,nd
             vel(idim,j)=vel(idim,j)+dippar(idim,1,i)*zdis1 + 
     1         dippar(idim,3,i)*dconjg(zdis1)
             vel(idim,j)=vel(idim,j)+dippar(idim,2,i)*dconjg(zdis2)*zdis
          enddo
 1111 continue
        enddo
      enddo
   
      return
      end
c
c
c
c
c

      subroutine bh2d_directdg(nd,sources,ns,dippar,
     1         targ,nt,vel,grad,thresh)
cf2py intent(in) nd,sources,ns,dippar,targ,thresh,nt
cf2py intent(out) vel,grad
c********************************************************************
c      This subroutine INCREMENTS the complex velocity VEL and its
c      gradients GRADA_Z, GRADA_ZBAR, GRADAA, at the
c      targets at targ(2) due to a collection of ns
c      dipoles at sources(2,ns)
c----------------------------------------------------------------------
c      INPUT parameters
c      nd           : number of densities
c      sources(2,ns): location of the sources
c      ns           : number of sources
c      dippar(3,ns) : dipole parameters 
c      targ         : target location
c      nt           : number of targets
c--------------------------------------------------------------------
c      OUTPUT
c      vel          : Complex velocity at target
c      grad         : Complex gradient 
c
c--------------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns,nt
      real *8 sources(2,ns), targ(2,nt)
      complex *16 vel(nd,nt),zs,zt,zdis
      complex *16 dippar(nd,3,ns)
      complex *16 zdis1,zdis2
      complex *16 grad(nd,3,nt),eye


      eye = dcmplx(0.0d0,1.0d0)
      do j=1,nt
        zt = dcmplx(targ(1,j),targ(2,j))
        do i=1,ns
           zs = dcmplx(sources(1,i),sources(2,i))
           zdis = zt-zs
           if(abs(zdis).le.thresh) goto 1111
           zdis1 = 1.0d0/zdis
           zdis2=zdis1**2
           do idim=1,nd
             vel(idim,j)=vel(idim,j)+dippar(idim,1,i)*zdis1 + 
     1         dippar(idim,3,i)*dconjg(zdis1)
             vel(idim,j)=vel(idim,j)+dippar(idim,2,i)*dconjg(zdis2)*zdis
             grad(idim,1,j)=grad(idim,1,j)-dippar(idim,1,i)*(zdis2)
              
             grad(idim,2,j)=grad(idim,2,j)+
     1           dippar(idim,2,i)*dconjg(zdis2)

             grad(idim,3,j)=grad(idim,3,j)-
     1           dippar(idim,3,i)*dconjg(zdis2)
             grad(idim,3,j)=grad(idim,3,j)-
     1           2*dippar(idim,2,i)*dconjg(zdis2*zdis1)*zdis
          enddo
 1111 continue
        enddo
      enddo
   
      return
      end

c
c
c

      subroutine bh2d_directcdp(nd,sources,ns,charges,dippar,
     1         targ,nt,vel,thresh)
cf2py intent(in) nd,sources,ns,charges,dippar,targ,thresh,nt
cf2py intent(out) vel
c********************************************************************
c      This subroutine INCREMENTS the complex velocity VEL 
c      at the
c      targets at targ(2) due to a collection of ns
c      charges and dipoles at sources(2,ns)
c----------------------------------------------------------------------
c      INPUT parameters
c      nd           : number of densities
c      sources(2,ns): location of the sources
c      ns           : number of sources
c      charges      : charge strength
c      dippar       : dipole parameters 
c      targ         : target location
c      nt           : number of targets
c--------------------------------------------------------------------
c      OUTPUT
c      vel          : Complex velocity at target
c
c--------------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(2,ns), targ(2,nt)
      complex *16 vel(nd,nt),zs,zt,zdis
      complex *16 charges(nd,2,ns),dippar(nd,3,ns)
      complex *16 zdis1,zdis2
      complex *16 eye

      eye = dcmplx(0,1.0d0)
      do j=1,nt
        zt = dcmplx(targ(1,j),targ(2,j))
        do i=1,ns
           zs = dcmplx(sources(1,i),sources(2,i))
           zdis = zt-zs
           if(abs(zdis).le.thresh) goto 1111
           zdis1 = 1.0d0/zdis
           zdis2=zdis1**2
           do idim=1,nd
             vel(idim,j)=vel(idim,j)+
     1            2*charges(idim,1,i)*log(cdabs(zdis))+
     1            charges(idim,2,i)*dconjg(zdis1)*zdis
     
             vel(idim,j)=vel(idim,j)+dippar(idim,1,i)*zdis1 + 
     1         dippar(idim,3,i)*dconjg(zdis1)
             vel(idim,j)=vel(idim,j)+dippar(idim,2,i)*dconjg(zdis2)*zdis
          enddo
 1111 continue
        enddo
      enddo
   
      return
      end
c
c
c
c
c

      subroutine bh2d_directcdg(nd,sources,ns,charges,dippar,
     1         targ,nt,vel,grad,thresh)
cf2py intent(in) nd,sources,ns,charges,dippar,targ,thresh,nt
cf2py intent(out) vel,grad
c********************************************************************
c      This subroutine INCREMENTS the complex velocity VEL and its
c      gradients GRADA, GRADAA, at the
c      targets at targ(2) due to a collection of ns
c      charges and dipoles at sources(2,ns)
c----------------------------------------------------------------------
c      INPUT parameters
c      nd           : number of densities
c      sources(2,ns): location of the sources
c      ns           : number of sources
c      charges      : charge strength
c      dippar       : dipole parameters 
c      targ         : target location
c--------------------------------------------------------------------
c      OUTPUT
c      vel          : Complex velocity at target
c      grad         : Complex gradient 
c
c--------------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns,nt
      real *8 sources(2,ns), targ(2,nt)
      complex *16 vel(nd,nt),zs,zt,zdis
      complex *16 charges(nd,2,ns),dippar(nd,3,ns)
      complex *16 zdis1,zdis2
      complex *16 grad(nd,3,nt),eye

      eye = dcmplx(0,1.0d0)
      do j=1,nt
        zt = dcmplx(targ(1,j),targ(2,j))
        do i=1,ns
           zs = dcmplx(sources(1,i),sources(2,i))
           zdis = zt-zs
           if(abs(zdis).le.thresh) goto 1111
           zdis1 = 1.0d0/zdis
           zdis2=zdis1**2
           do idim=1,nd
             vel(idim,j)=vel(idim,j)+
     1           2*charges(idim,1,i)*log(cdabs(zdis))+
     1           charges(idim,2,i)*dconjg(zdis1)*zdis
     
             vel(idim,j)=vel(idim,j)+dippar(idim,1,i)*zdis1 + 
     1         dippar(idim,3,i)*dconjg(zdis1)
             vel(idim,j)=vel(idim,j)+dippar(idim,2,i)*dconjg(zdis2)*zdis

             grad(idim,1,j)=grad(idim,1,j)+charges(idim,1,i)*zdis1
             grad(idim,1,j)=grad(idim,1,j)-dippar(idim,1,i)*(zdis2)

             grad(idim,2,j)=grad(idim,2,j) +
     1           charges(idim,2,i)*dconjg(zdis1)
             grad(idim,2,j)=grad(idim,2,j)+
     1           dippar(idim,2,i)*dconjg(zdis2)

             grad(idim,3,j)=grad(idim,3,j)+
     1          charges(idim,1,i)*dconjg(zdis1)
             grad(idim,3,j)=grad(idim,3,j)-
     1          charges(idim,2,i)*dconjg(zdis2)*zdis
             grad(idim,3,j)=grad(idim,3,j)-
     1          dippar(idim,3,i)*dconjg(zdis2)
             grad(idim,3,j)=grad(idim,3,j)-
     1         2*dippar(idim,2,i)*dconjg(zdis2*zdis1)*zdis
          enddo
 1111 continue
        enddo
      enddo
   
      return
      end

