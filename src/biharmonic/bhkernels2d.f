c
c
c
c

      subroutine bh2d_directcp_vec(nd,sources,ns,charges,targ,vel,
     1         thresh)
c********************************************************************
c      This subroutine INCREMENTS the complex velocity VEL 
c      at the targets at targ(2) due to a collection of ns
c      charges at sources(2,ns)
c      
c      In this subroutine the following rescaled versions of charges and
c      dipoles are used to compute the velocity field due to them
c
c      The complex velocity at zt = xt + i yt 
c      due to a charge a complex charge c1 at zs = xs+i ys is given by
c
c
c      vel = 2*c1 log|zt-zs| + c1_bar (zt-zs)/(zt_bar - zs_bar)
c            
c      The complex velocity due to a dipole with complex parameters
c      c2,c3 is given by
c
c
c      vel = c2/(zt-zs) - c2_bar (zt-zs)/(zt_bar-zs_bar)^2 + 
c             c3/(zt_bar - zs_bar)
c
c      The corresponding goursat functions (\phi \psi) for derivation 
c      of the complex velocity=\phi + z d/dz (\phi)_bar + \psi_bar
c
c      \phi(z) = c1 log (z-zs) + d1/(z-zs)
c      \psi(z) = d2_bar/(z-zs) + zs_bar d1/(z-zs)^2 + c1_bar log(z-zs)-
c               zs_bar c1 /(z-zs)
c
c      Analytic component of the gradient (grada)= d/dz (\phi(z))
c
c      Anti analytic component of the gradient (gradaa) = 
c      z (d^2/ dz^2 (\phi))_bar + (d/dz (\psi))_bar
c----------------------------------------------------------------------
c      INPUT parameters
c      nd           : number of densities
c      sources(2,ns): location of the sources
c      ns           : number of sources
c      charges      : charge strength
c      targ         : target location
c--------------------------------------------------------------------
c      OUTPUT
c      vel          : Complex velocity at target
c
c--------------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(2,ns), targ(2)
      complex *16 vel(nd),zs,zt,zdis
      complex *16 charges(nd,ns)
      complex *16 zdis1,zdis2
      complex *16 eye

      eye = dcmplx(0,1.0d0)
      zt = dcmplx(targ(1),targ(2))
      do i=1,ns
         zs = dcmplx(sources(1,i),sources(2,i))
         zdis = zt-zs
         if(abs(zdis).le.thresh) goto 1111
         zdis1 = 1.0d0/zdis
         do idim=1,nd
           vel(idim)=vel(idim)+2*charges(idim,i)*log(cdabs(zdis))+
     1                     dconjg(charges(idim,i)*zdis1)*zdis
        enddo
 1111 continue
      enddo
   
      return
      end
c
c
c
c
c

      subroutine bh2d_directcg_vec(nd,sources,ns,charges,
     1         targ,vel,grad,thresh)
c********************************************************************
c      This subroutine INCREMENTS the complex velocity VEL and its
c      gradients GRADA, GRADAA, at the
c      targets at targ(2) due to a collection of ns
c      charges at sources(2,ns)
c      
c      In this subroutine the following rescaled versions of charges and
c      dipoles are used to compute the velocity field due to them
c
c      The complex velocity at zt = xt + i yt 
c      due to a charge a complex charge c1 at zs = xs+i ys is given by
c
c
c      vel = 2*c1 log|zt-zs| + c1_bar (zt-zs)/(zt_bar - zs_bar)
c            
c      The complex velocity due to a dipole with complex parameters
c      c2,c3 is given by
c
c
c      vel = c2/(zt-zs) - c2_bar (zt-zs)/(zt_bar-zs_bar)^2 + 
c             c3/(zt_bar - zs_bar)
c
c      The corresponding goursat functions (\phi \psi) for derivation 
c      of the complex velocity=\phi + z d/dz (\phi)_bar + \psi_bar
c
c      \phi(z) = c1 log (z-zs) + d1/(z-zs)
c      \psi(z) = d2_bar/(z-zs) + zs_bar d1/(z-zs)^2 + c1_bar log(z-zs)-
c               zs_bar c1 /(z-zs)
c
c      Analytic component of the gradient (grada)= d/dz (\phi(z))
c
c      Anti analytic component of the gradient (gradaa) = 
c      z (d^2/ dz^2 (\phi))_bar + (d/dz (\psi))_bar
c----------------------------------------------------------------------
c      INPUT parameters
c      nd           : number of densities
c      sources(2,ns): location of the sources
c      ns           : number of sources
c      charges      : charge strength
c      targ         : target location
c--------------------------------------------------------------------
c      OUTPUT
c      vel          : Complex velocity at target
c      grad         : Complex gradient (d/dz, d/dzbar)
c
c--------------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(2,ns), targ(2)
      complex *16 vel(nd),zs,zt,zdis
      complex *16 charges(nd,ns)
      complex *16 zdis1,zdis2
      complex *16 grad(nd,2),eye

      eye = dcmplx(0,1.0d0)
      zt = dcmplx(targ(1),targ(2))
      do i=1,ns
         zs = dcmplx(sources(1,i),sources(2,i))
         zdis = zt-zs
         if(abs(zdis).le.thresh) goto 1111
         zdis1 = 1.0d0/zdis
         zdis2=zdis1**2
         do idim=1,nd
           vel(idim)=vel(idim)+2*charges(idim,i)*log(cdabs(zdis))+
     1                     dconjg(charges(idim,i)*zdis1)*zdis
     
           grad(idim,1)=grad(idim,1)+charges(idim,i)*zdis1
           grad(idim,2)=grad(idim,2)+charges(idim,i)*dconjg(zdis1)
           grad(idim,2)=grad(idim,2)-dconjg(charges(idim,i)*zdis2)*zdis
        enddo
 1111 continue
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

      subroutine bh2d_directdp_vec(nd,sources,ns,dippar,
     1         targ,vel,thresh)
c********************************************************************
c      This subroutine INCREMENTS the complex velocity VEL 
c      at the
c      targets at targ(2) due to a collection of ns
c      dipoles at sources(2,ns)
c      
c      In this subroutine the following rescaled versions of charges and
c      dipoles are used to compute the velocity field due to them
c
c      The complex velocity at zt = xt + i yt 
c      due to a charge a complex charge c1 at zs = xs+i ys is given by
c
c
c      vel = 2*c1 log|zt-zs| + c1_bar (zt-zs)/(zt_bar - zs_bar)
c            
c      The complex velocity due to a dipole with complex parameters
c      c2,c3 is given by
c
c
c      vel = c2/(zt-zs) - c2_bar (zt-zs)/(zt_bar-zs_bar)^2 + 
c             c3/(zt_bar - zs_bar)
c
c      The corresponding goursat functions (\phi \psi) for derivation 
c      of the complex velocity=\phi + z d/dz (\phi)_bar + \psi_bar
c
c      \phi(z) = c1 log (z-zs) + d1/(z-zs)
c      \psi(z) = d2_bar/(z-zs) + zs_bar d1/(z-zs)^2 + c1_bar log(z-zs)-
c               zs_bar c1 /(z-zs)
c
c      Analytic component of the gradient (grada)= d/dz (\phi(z))
c
c      Anti analytic component of the gradient (gradaa) = 
c      z (d^2/ dz^2 (\phi))_bar + (d/dz (\psi))_bar
c----------------------------------------------------------------------
c      INPUT parameters
c      nd           : number of densities
c      sources(2,ns): location of the sources
c      ns           : number of sources
c      dippar       : dipole parameters (corresponding to c2,c3 in
c                     above expression)
c      targ         : target location
c--------------------------------------------------------------------
c      OUTPUT
c      vel          : Complex velocity at target
c
c--------------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(2,ns), targ(2)
      complex *16 vel(nd),zs,zt,zdis
      complex *16 dippar(nd,2,ns)
      complex *16 zdis1,zdis2
      complex *16 eye

      eye = dcmplx(0,1.0d0)
      zt = dcmplx(targ(1),targ(2))
      do i=1,ns
         zs = dcmplx(sources(1,i),sources(2,i))
         zdis = zt-zs
         if(abs(zdis).le.thresh) goto 1111
         zdis1 = 1.0d0/zdis
         zdis2=zdis1**2
         do idim=1,nd
           vel(idim)=vel(idim)+dippar(idim,1,i)*zdis1 + 
     1         dippar(idim,2,i)*dconjg(zdis1)
           vel(idim)=vel(idim)-dconjg(dippar(idim,1,i)*zdis2)*zdis
        enddo
 1111 continue
      enddo
   
      return
      end
c
c
c
c
c

      subroutine bh2d_directdg_vec(nd,sources,ns,charges,dippar,
     1         targ,vel,grad,thresh)
c********************************************************************
c      This subroutine INCREMENTS the complex velocity VEL and its
c      gradients GRADA, GRADAA, at the
c      targets at targ(2) due to a collection of ns
c      dipoles at sources(2,ns)
c      
c      In this subroutine the following rescaled versions of charges and
c      dipoles are used to compute the velocity field due to them
c
c      The complex velocity at zt = xt + i yt 
c      due to a charge a complex charge c1 at zs = xs+i ys is given by
c
c
c      vel = 2*c1 log|zt-zs| + c1_bar (zt-zs)/(zt_bar - zs_bar)
c            
c      The complex velocity due to a dipole with complex parameters
c      c2,c3 is given by
c
c
c      vel = c2/(zt-zs) - c2_bar (zt-zs)/(zt_bar-zs_bar)^2 + 
c             c3/(zt_bar - zs_bar)
c
c      The corresponding goursat functions (\phi \psi) for derivation 
c      of the complex velocity=\phi + z d/dz (\phi)_bar + \psi_bar
c
c      \phi(z) = c1 log (z-zs) + d1/(z-zs)
c      \psi(z) = d2_bar/(z-zs) + zs_bar d1/(z-zs)^2 + c1_bar log(z-zs)-
c               zs_bar c1 /(z-zs)
c
c      Analytic component of the gradient (grada)= d/dz (\phi(z))
c
c      Anti analytic component of the gradient (gradaa) = 
c      z (d^2/ dz^2 (\phi))_bar + (d/dz (\psi))_bar
c----------------------------------------------------------------------
c      INPUT parameters
c      nd           : number of densities
c      sources(2,ns): location of the sources
c      ns           : number of sources
c      dippar       : dipole parameters (corresponding to c2,c3 in
c                     above expression)
c      targ         : target location
c--------------------------------------------------------------------
c      OUTPUT
c      vel          : Complex velocity at target
c      grad         : Complex gradient (d/dz, d/dzbar)
c
c--------------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(2,ns), targ(2)
      complex *16 vel(nd),zs,zt,zdis
      complex *16 dippar(nd,2,ns)
      complex *16 zdis1,zdis2
      complex *16 grad(nd,2),eye

      eye = dcmplx(0,1.0d0)
      zt = dcmplx(targ(1),targ(2))
      do i=1,ns
         zs = dcmplx(sources(1,i),sources(2,i))
         zdis = zt-zs
         if(abs(zdis).le.thresh) goto 1111
         zdis1 = 1.0d0/zdis
         zdis2=zdis1**2
         do idim=1,nd
           vel(idim)=vel(idim)+dippar(idim,1,i)*zdis1 + 
     1         dippar(idim,2,i)*dconjg(zdis1)
           vel(idim)=vel(idim)-dconjg(dippar(idim,1,i)*zdis2)*zdis
           grad(idim,1)=grad(idim,1)-dippar(idim,1,i)*(zdis2)
           grad(idim,2)=grad(idim,2)-dippar(idim,2,i)*dconjg(zdis2)
           grad(idim,2)=grad(idim,2)+
     1         2*dconjg(dippar(idim,1,i)*zdis2*zdis1)*zdis
        enddo
 1111 continue
      enddo
   
      return
      end

c
c
c

      subroutine bh2d_directcdp_vec(nd,sources,ns,charges,dippar,
     1         targ,vel,thresh)
c********************************************************************
c      This subroutine INCREMENTS the complex velocity VEL 
c      at the
c      targets at targ(2) due to a collection of ns
c      charges and dipoles at sources(2,ns)
c      
c      In this subroutine the following rescaled versions of charges and
c      dipoles are used to compute the velocity field due to them
c
c      The complex velocity at zt = xt + i yt 
c      due to a charge a complex charge c1 at zs = xs+i ys is given by
c
c
c      vel = 2*c1 log|zt-zs| + c1_bar (zt-zs)/(zt_bar - zs_bar)
c            
c      The complex velocity due to a dipole with complex parameters
c      c2,c3 is given by
c
c
c      vel = c2/(zt-zs) - c2_bar (zt-zs)/(zt_bar-zs_bar)^2 + 
c             c3/(zt_bar - zs_bar)
c
c      The corresponding goursat functions (\phi \psi) for derivation 
c      of the complex velocity=\phi + z d/dz (\phi)_bar + \psi_bar
c
c      \phi(z) = c1 log (z-zs) + d1/(z-zs)
c      \psi(z) = d2_bar/(z-zs) + zs_bar d1/(z-zs)^2 + c1_bar log(z-zs)-
c               zs_bar c1 /(z-zs)
c
c      Analytic component of the gradient (grada)= d/dz (\phi(z))
c
c      Anti analytic component of the gradient (gradaa) = 
c      z (d^2/ dz^2 (\phi))_bar + (d/dz (\psi))_bar
c----------------------------------------------------------------------
c      INPUT parameters
c      nd           : number of densities
c      sources(2,ns): location of the sources
c      ns           : number of sources
c      charges      : charge strength
c      dippar       : dipole parameters (corresponding to c2,c3 in
c                     above expression)
c      targ         : target location
c--------------------------------------------------------------------
c      OUTPUT
c      vel          : Complex velocity at target
c
c--------------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(2,ns), targ(2)
      complex *16 vel(nd),zs,zt,zdis
      complex *16 charges(nd,ns),dippar(nd,2,ns)
      complex *16 zdis1,zdis2
      complex *16 eye

      eye = dcmplx(0,1.0d0)
      zt = dcmplx(targ(1),targ(2))
      do i=1,ns
         zs = dcmplx(sources(1,i),sources(2,i))
         zdis = zt-zs
         if(abs(zdis).le.thresh) goto 1111
         zdis1 = 1.0d0/zdis
         zdis2=zdis1**2
         do idim=1,nd
           vel(idim)=vel(idim)+2*charges(idim,i)*log(cdabs(zdis))+
     1                     dconjg(charges(idim,i)*zdis1)*zdis
     
           vel(idim)=vel(idim)+dippar(idim,1,i)*zdis1 + 
     1         dippar(idim,2,i)*dconjg(zdis1)
           vel(idim)=vel(idim)-dconjg(dippar(idim,1,i)*zdis2)*zdis
        enddo
 1111 continue
      enddo
   
      return
      end
c
c
c
c
c

      subroutine bh2d_directcdg_vec(nd,sources,ns,charges,dippar,
     1         targ,vel,grad,thresh)
c********************************************************************
c      This subroutine INCREMENTS the complex velocity VEL and its
c      gradients GRADA, GRADAA, at the
c      targets at targ(2) due to a collection of ns
c      charges and dipoles at sources(2,ns)
c      
c      In this subroutine the following rescaled versions of charges and
c      dipoles are used to compute the velocity field due to them
c
c      The complex velocity at zt = xt + i yt 
c      due to a charge a complex charge c1 at zs = xs+i ys is given by
c
c
c      vel = 2*c1 log|zt-zs| + c1_bar (zt-zs)/(zt_bar - zs_bar)
c            
c      The complex velocity due to a dipole with complex parameters
c      c2,c3 is given by
c
c
c      vel = c2/(zt-zs) - c2_bar (zt-zs)/(zt_bar-zs_bar)^2 + 
c             c3/(zt_bar - zs_bar)
c
c      The corresponding goursat functions (\phi \psi) for derivation 
c      of the complex velocity=\phi + z d/dz (\phi)_bar + \psi_bar
c
c      \phi(z) = c1 log (z-zs) + d1/(z-zs)
c      \psi(z) = d2_bar/(z-zs) + zs_bar d1/(z-zs)^2 + c1_bar log(z-zs)-
c               zs_bar c1 /(z-zs)
c
c      Analytic component of the gradient (grada)= d/dz (\phi(z))
c
c      Anti analytic component of the gradient (gradaa) = 
c      z (d^2/ dz^2 (\phi))_bar + (d/dz (\psi))_bar
c----------------------------------------------------------------------
c      INPUT parameters
c      nd           : number of densities
c      sources(2,ns): location of the sources
c      ns           : number of sources
c      charges      : charge strength
c      dippar       : dipole parameters (corresponding to c2,c3 in
c                     above expression)
c      targ         : target location
c--------------------------------------------------------------------
c      OUTPUT
c      vel          : Complex velocity at target
c      grad         : Complex gradient (d/dz, d/dzbar)
c
c--------------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns
      real *8 sources(2,ns), targ(2)
      complex *16 vel(nd),zs,zt,zdis
      complex *16 charges(nd,ns),dippar(nd,2,ns)
      complex *16 zdis1,zdis2
      complex *16 grad(nd,2),eye

      eye = dcmplx(0,1.0d0)
      zt = dcmplx(targ(1),targ(2))
      do i=1,ns
         zs = dcmplx(sources(1,i),sources(2,i))
         zdis = zt-zs
         if(abs(zdis).le.thresh) goto 1111
         zdis1 = 1.0d0/zdis
         zdis2=zdis1**2
         do idim=1,nd
           vel(idim)=vel(idim)+2*charges(idim,i)*log(cdabs(zdis))+
     1                     dconjg(charges(idim,i)*zdis1)*zdis
     
           vel(idim)=vel(idim)+dippar(idim,1,i)*zdis1 + 
     1         dippar(idim,2,i)*dconjg(zdis1)
           vel(idim)=vel(idim)-dconjg(dippar(idim,1,i)*zdis2)*zdis
           grad(idim,1)=grad(idim,1)+charges(idim,i)*zdis1
           grad(idim,1)=grad(idim,1)-dippar(idim,1,i)*(zdis2)
           grad(idim,2)=grad(idim,2)+charges(idim,i)*dconjg(zdis1)
           grad(idim,2)=grad(idim,2)-dconjg(charges(idim,i)*zdis2)*zdis
           grad(idim,2)=grad(idim,2)-dippar(idim,2,i)*dconjg(zdis2)
           grad(idim,2)=grad(idim,2)+
     1         2*dconjg(dippar(idim,1,i)*zdis2*zdis1)*zdis
        enddo
 1111 continue
      enddo
   
      return
      end

