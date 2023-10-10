c    This file contains basic subroutines for forming and
c    evaluating multipole expansions for Stokes flow
c    in 2D
c
c-----------------------------------------------------------------
c     Expansion formation subroutines
c
c     bh2dformmpc: creates multipole expansion (outgoing) due to 
c                     a collection of charge sources
c     bh2dformmpd: creates multipole expansion (outgoing) due to 
c                     a collection of dipole sources
c     bh2dformmpcd: creates multipole expansion (outgoing) due to 
c                     a collection of charge and dipole sources
c     bh2dformta: creates local expansion (incoming) due to 
c                     a collection of charge sources
c     bh2dformtad: creates multipole expansion (incoming) due to 
c                     a collection of dipole sources
c     bh2dformtacd: creates multipole expansion (incoming) due to 
c                     a collection of charge and dipole sources
c
c------------------------------------------------------------------     
c     Multipole and local translation operators
c
c     bh2dmpmp: Converts multipole expansion to a multipole expansion
c     bh2dmploc: Converts multipole expansion to local expansion
c     bh2dlocloc: converts local expansion to local expansion
c
c---------------------------------------------------------------------
c     Evaluation subroutines
c 
c     bh2dmpevalp: Computes the complex velocity due to a 
c                      multipole expansion at a collection of targets
c     bh2dmpevalg: Computes the complex velocity and its gradients 
c                      due to a multipole expansion at a collection 
c                      of targets
c     bh2dtaevalp: Computes the complex velocity due to a 
c                      local expansion at a collection of targets
c     bh2dtaevalg: Computes the complex velocity and its gradients 
c                      due to a local expansion at a collection 
c                      of targets
c----------------------------------------------------------------------
c     Utility subroutines
c     
c     bh2dzero: utility to intialize multipole/local coefficients to
c                   zero
c    
c*********************************************************************
c
c

      subroutine bh2dmpevalp(nd,rscale,center,mpole,nterms,
     1                     ztarg,ntarg,vel)
c********************************************************************
c     This subroutine evaluates the multipole expansion about
c     the "center" due to the multipole expansion
c     at the target ztarg
c
c     The multipole expansion has 5 set of coefficients
c     vel = \sum mpole(1,k)/z^k + \sum mpole(2,k)/z_bar^k +
c            k                k
c      
c         z \sum mpole(3,k)/z_bar^k + Re(mp4(0) log(z) + \sum mpole(4,k)/z^k)
c              k                                       k
c
c          + i Re(mp5(0) log(z) + \sum mpole(5,k)/z^k)
c                                  k
c      z = (ztarg - center)/rscale
c----------------------------------------------------------------------
c      INPUT parameters
c      rscale        : scaling parameter
c      center        : expansion center
c      mpole         : coeffs for multipole expansions
c      mpole(1,:)    : coeffs for multipole expansion - laurent
c      mpole(2,:)    : coeffs for multipole expansion - antilaurent 
c      mpole(3,:)    : coeffs for multipole expansion - z antilaurent
c      mpole(4,:)    : coeffs for multipole expansion - logsource1
c      mpole(5,:)    : coeffs for multipole expansion - logsource2
c      nterms        : Number of terms in multipole expansion
c      ztarg         : target location
c      ntarg         : number of targets
c---------------------------------------------------------------------
c     OUTPUT
c     vel - Complex velocity
c*********************************************************************

      implicit real *8 (a-h,o-z)
      integer nterms,ifgrada,ifgradaa,nmax
      real *8 rscale,center(2),ztarg(2,ntarg)
      complex *16 zc,zt,zdis,eye
      complex *16 mpole(nd,5,0:nterms)
      complex *16 vel(nd,ntarg)
      complex *16 ztemp,ztemp1
      complex *16 zdisinv,zpow(1:nterms+5)


      nmax = nterms+3
      rinv=1.0d0/rscale
      zc = dcmplx(center(1),center(2))
      eye = dcmplx(0.0d0,1.0d0)
      do k=1,ntarg
        zt = dcmplx(ztarg(1,k),ztarg(2,k))
        zdis = zt - zc
        zdisinv = 1.0d0/zdis
        ztemp = rscale*zdisinv
        zpow(1)=ztemp
        do i=2,nmax
          zpow(i)=zpow(i-1)*ztemp
        enddo
        do idim=1,nd
          vel(idim,k) = vel(idim,k)+
     1        (mpole(idim,4,0)+eye*mpole(idim,5,0))*log(cdabs(zdis))
        enddo
   
        do i=1,nterms
          do idim=1,nd
            vel(idim,k) = vel(idim,k) + 
     1        mpole(idim,1,i)*zpow(i) + mpole(idim,2,i)*dconjg(zpow(i))
            vel(idim,k) = vel(idim,k) + 
     1        mpole(idim,3,i)*dconjg(zpow(i))*zdis
            vel(idim,k) = vel(idim,k) + 
     1        dreal(mpole(idim,4,i)*zpow(i))+
     2         eye*dreal(mpole(idim,5,i)*zpow(i))
          enddo
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

      subroutine bh2dmpevalg(nd,rscale,center,mpole,nterms,
     1                     ztarg,ntarg,vel,grad)
c********************************************************************
c     This subroutine evaluates the multipole expansion about
c     the "center" due to the multipole expansion
c     at the target ztarg
c
c     The multipole expansion has 5 set of coefficients
c     vel = \sum mpole(1,k)/z^k + \sum mpole(2,k)/z_bar^k +
c            k                k
c      
c         z \sum mpole(3,k)/z_bar^k + Re(mp4(0) log(z) + \sum mpole(4,k)/z^k)
c              k                                       k
c
c          + i Re(mp5(0) log(z) + \sum mpole(5,k)/z^k)
c                                  k
c      z = (ztarg - center)/rscale
c
c      note that we use log(|z|) = 1/2(log(z) + log(zbar))
c
c      grada_z = gradient (analytic comoponent of vel, projected onto z
c              component), i.e.
c              terms involving mpole1, mpole4, mpole5 and ignoring z \sum_k
c              mpole(3,k)/z_bar^k
c      grada_zbar = gradient (analytic component of grad vel, projected onto
c              zbar component), i.e.
c              terms involving mpole3 and ignoring mpole1 terms
c      gradaa = gradient(all but analytic component of vel)
c             = gradient(anti analytic + z times anti analytic 
c                         component of vel)
c
c----------------------------------------------------------------------
c      INPUT parameters
c      rscale        : scaling parameter
c      center        : expansion center
c      mpole         : coeffs for multipole expansions
c      mpole(1,:)    : coeffs for multipole expansion - laurent
c      mpole(2,:)    : coeffs for multipole expansion - antilaurent 
c      mpole(3,:)    : coeffs for multipole expansion - z antilaurent
c      mpole(4,:)    : coeffs for multipole expansion - logsource1
c      mpole(5,:)    : coeffs for multipole expansion - logsource2
c      nterms        : Number of terms in multipole expansion
c      ztarg         : target location
c      ntarg         : number of targets
c---------------------------------------------------------------------
c     OUTPUT
c     vel - Complex velocity
c     grad(1:3) - gradients (grada_z, grada_zbar, gradaa) 
c*********************************************************************

      implicit real *8 (a-h,o-z)
      integer nterms,ifgrada,ifgradaa,nmax
      real *8 rscale,center(2),ztarg(2,ntarg)
      complex *16 zc,zt,zdis,eye
      complex *16 mpole(nd,5,0:nterms)
      complex *16 vel(nd,ntarg),grad(nd,3,ntarg)
      complex *16 ztemp,ztemp1
      complex *16 zdisinv,zpow(1:nterms+5)


      nmax = nterms+3
      rinv=1.0d0/rscale
      zc = dcmplx(center(1),center(2))
      eye = dcmplx(0.0d0,1.0d0)
      do k=1,ntarg
        zt = dcmplx(ztarg(1,k),ztarg(2,k))
        zdis = zt - zc
        zdisinv = 1.0d0/zdis
        ztemp = rscale*zdisinv
        zpow(1)=ztemp
        do i=2,nmax
          zpow(i)=zpow(i-1)*ztemp
        enddo
        do idim=1,nd
          vel(idim,k) = vel(idim,k)+
     1        (mpole(idim,4,0)+eye*mpole(idim,5,0))*log(cdabs(zdis))
          grad(idim,1,k) = grad(idim,1,k) + 
     1        0.5d0*(mpole(idim,4,0)+eye*mpole(idim,5,0))*zpow(1)*rinv
          grad(idim,3,k)= grad(idim,3,k) + 
     1       0.5d0*(mpole(idim,4,0)+eye*mpole(idim,5,0))*
     2       dconjg(zpow(1))*rinv
        enddo
   
        do i=1,nterms
          do idim=1,nd
            vel(idim,k) = vel(idim,k) + 
     1        mpole(idim,1,i)*zpow(i) + mpole(idim,2,i)*dconjg(zpow(i))
            vel(idim,k) = vel(idim,k) + 
     1        mpole(idim,3,i)*dconjg(zpow(i))*zdis
            vel(idim,k) = vel(idim,k) + 
     1        dreal(mpole(idim,4,i)*zpow(i))+
     2         eye*dreal(mpole(idim,5,i)*zpow(i))

            grad(idim,1,k) = grad(idim,1,k) - 
     1         mpole(idim,1,i)*zpow(i+1)*i*rinv
            grad(idim,1,k) = grad(idim,1,k) - 
     1         0.5*mpole(idim,4,i)*i*zpow(i+1)*rinv
            grad(idim,1,k) = grad(idim,1,k) - 
     1         eye*0.5*mpole(idim,5,i)*i*zpow(i+1)*rinv
            
            grad(idim,2,k) = grad(idim,2,k) + 
     1           mpole(idim,3,i)*dconjg(zpow(i))

            grad(idim,3,k) = grad(idim,3,k) - 
     1         mpole(idim,2,i)*dconjg(zpow(i+1))*i*rinv
            grad(idim,3,k) = grad(idim,3,k) - 
     1         zdis*mpole(idim,3,i)*dconjg(zpow(i+1))*i*rinv
            grad(idim,3,k) = grad(idim,3,k) - 
     1         dconjg(0.5*mpole(idim,4,i)*i*zpow(i+1))*rinv
            grad(idim,3,k) = grad(idim,3,k) + 
     1         dconjg(eye*0.5*mpole(idim,5,i)*i*zpow(i+1))*rinv
          enddo
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
c********************************************************************
      subroutine bh2dtaevalp(nd,rscale,center,mpole,nterms,
     1                      ztarg,ntarg,vel)
c********************************************************************
c     This subroutine evaluates the multipole expansion about
c     the "center" due to the local expansion
c     at the target ztarg
c
c     The multipole expansion has 5 set of coefficients
c     vel = \sum mpole(1,k) z^k + \sum mpole(2,k) z_bar^k +
c            k                k
c      
c         z \sum mpole(3,k) z_bar^k + Re(\sum mpole(4,k) z^k)
c              k                       k
c
c          + i Re( \sum mpole(5,k) z^k)
c                    k
c      z = (ztarg - center)/rscale
c----------------------------------------------------------------------
c      INPUT parameters
c      rscale        : scaling parameter
c      center        : expansion center
c      mpole         : coeffs for local expansion
c      mpole(1,:)    : coeffs for multipole expansion - taylor
c      mpole(2,:)    : coeffs for multipole expansion - antitaylor 
c      mpole(3,:)    : coeffs for multipole expansion - z antitaylor
c      mpole(4,:)    : coeffs for multipole expansion - logsource1
c      mpole(5,:)    : coeffs for multipole expansion - logsource2
c      nterms        : Number of terms in multipole expansion
c      ztarg         : target location
c      ntarg         : number of targets
c---------------------------------------------------------------------
c     OUTPUT
c     vel - Complex velocity
c*********************************************************************

      implicit real *8 (a-h,o-z)
      integer nterms,ifgrada,ifgradaa
      real *8 rscale,center(2),ztarg(2,ntarg)
      complex *16 zc,zt,zdis,eye
      complex *16 mpole(nd,5,0:nterms)
      complex *16 vel(nd,ntarg)
      complex *16 ztemp,ztemp1,zpow(0:nterms)

      rinv=1.0d0/rscale
      zc = dcmplx(center(1),center(2))
      eye = dcmplx(0.0d0,1.0d0)

      do k=1,ntarg
        zt = dcmplx(ztarg(1,k),ztarg(2,k))
        zdis = zt - zc
        ztemp1 = zdis*rinv
        zpow(0)=1.0d0
        do i=1,nterms
          zpow(i)=zpow(i-1)*ztemp1
        enddo
      
        do i=0,nterms
          do idim=1,nd
            vel(idim,k) = vel(idim,k) + mpole(idim,1,i)*zpow(i) + 
     1          mpole(idim,2,i)*dconjg(zpow(i))
            vel(idim,k) = vel(idim,k) + 
     1          mpole(idim,3,i)*dconjg(zpow(i))*zdis
            vel(idim,k) = vel(idim,k) + 
     1          dreal(mpole(idim,4,i)*zpow(i)) + 
     2          eye*dreal(mpole(idim,5,i)*zpow(i))
          enddo
        enddo
      enddo

   
      return
      end

c
c
c
c
c********************************************************************
      subroutine bh2dtaevalg(nd,rscale,center,mpole,nterms,
     1                      ztarg,ntarg,vel,grad)
c********************************************************************
c     This subroutine evaluates the multipole expansion about
c     the "center" due to the local expansion
c     at the target ztarg
c
c     The multipole expansion has 5 set of coefficients
c     vel = \sum mpole(1,k) z^k + \sum mpole(2,k) z_bar^k +
c            k                k
c      
c         z \sum mpole(3,k) z_bar^k + Re(\sum mpole(4,k) z^k)
c              k                       k
c
c          + i Re( \sum mpole(5,k) z^k)
c                    k
c      z = (ztarg - center)/rscale
c
c      note that we use log(|z|) = 1/2(log(z) + log(zbar))
c
c      grada_z = gradient (analytic comoponent of vel, projected onto z
c              component), i.e.
c              terms involving mpole1, mpole4, mpole5 and ignoring z \sum_k
c              mpole(3,k)/z_bar^k
c      grada_zbar = gradient (analytic component of grad vel, projected onto
c              zbar component), i.e.
c              terms involving mpole3 and ignoring mpole1 terms
c      gradaa = gradient(all but analytic component of vel)
c             = gradient(anti analytic + z times anti analytic 
c                         component of vel)
c
c----------------------------------------------------------------------
c      INPUT parameters
c      rscale        : scaling parameter
c      center        : expansion center
c      mpole         : coeffs for local expansion
c      mpole(1,:)    : coeffs for multipole expansion - taylor
c      mpole(2,:)    : coeffs for multipole expansion - antitaylor 
c      mpole(3,:)    : coeffs for multipole expansion - z antitaylor
c      mpole(4,:)    : coeffs for multipole expansion - logsource1
c      mpole(5,:)    : coeffs for multipole expansion - logsource2
c      nterms        : Number of terms in multipole expansion
c      ztarg         : target location
c      ntarg         : number of targets
c---------------------------------------------------------------------
c     OUTPUT
c     vel - Complex velocity
c     grad - gradient (grada_z, grada_zbar, gradaa) 
c*********************************************************************

      implicit real *8 (a-h,o-z)
      integer nterms,ifgrada,ifgradaa
      real *8 rscale,center(2),ztarg(2,ntarg)
      complex *16 zc,zt,zdis,eye
      complex *16 mpole(nd,5,0:nterms)
      complex *16 vel(nd,ntarg),grad(nd,3,ntarg)
      complex *16 ztemp,ztemp1,zpow(0:nterms)

      rinv=1.0d0/rscale
      zc = dcmplx(center(1),center(2))
      eye = dcmplx(0.0d0,1.0d0)

      do k=1,ntarg
        zt = dcmplx(ztarg(1,k),ztarg(2,k))
        zdis = zt - zc
        ztemp1 = zdis*rinv
        zpow(0)=1.0d0
        do i=1,nterms
          zpow(i)=zpow(i-1)*ztemp1
        enddo
      
        do i=0,nterms
          do idim=1,nd
            vel(idim,k) = vel(idim,k) + mpole(idim,1,i)*zpow(i) + 
     1          mpole(idim,2,i)*dconjg(zpow(i))
            vel(idim,k) = vel(idim,k) + 
     1          mpole(idim,3,i)*dconjg(zpow(i))*zdis
            vel(idim,k) = vel(idim,k) + 
     1          dreal(mpole(idim,4,i)*zpow(i)) + 
     2          eye*dreal(mpole(idim,5,i)*zpow(i))

          enddo
        enddo

        do idim=1,nd
          grad(idim,2,k) = grad(idim,2,k) + mpole(idim,3,0)
        enddo
        do i=1,nterms
          do idim=1,nd
            grad(idim,1,k) = grad(idim,1,k) + 
     1         mpole(idim,1,i)*zpow(i-1)*i*rinv
            grad(idim,1,k) = grad(idim,1,k) + 
     1         0.5*mpole(idim,4,i)*i*zpow(i-1)*rinv
            grad(idim,1,k) = grad(idim,1,k) + 
     1         eye*0.5*mpole(idim,5,i)*i*zpow(i-1)*rinv

            grad(idim,2,k) = grad(idim,2,k) + 
     1         mpole(idim,3,i)*dconjg(zpow(i))

            grad(idim,3,k) = grad(idim,3,k) + 
     1          mpole(idim,2,i)*dconjg(zpow(i-1))*i*rinv
            grad(idim,3,k) = grad(idim,3,k) + 
     1          mpole(idim,3,i)*dconjg(zpow(i-1))*i*zdis*rinv
            grad(idim,3,k) = grad(idim,3,k) + 
     1         dconjg(0.5*mpole(idim,4,i)*i*zpow(i-1))*rinv
            grad(idim,3,k) = grad(idim,3,k) - 
     1         dconjg(eye*0.5*mpole(idim,5,i)*i*zpow(i-1))*rinv
          enddo
        enddo
      enddo

   
      return
      end

c*******************************************************************
c     EXPANSION FORMATION
c*******************************************************************
      subroutine bh2dformmpd(nd,rscale,sources,ns,dip,
     1           center,nterms,mpole)
c*****************************************************************
c     This subroutine computes the multipole expansion about
c     the "center" due to ns sources located at sources(2,*)
c
c      The complex velocity at zt = xt + i yt 
c      due to a charge a complex charge (c1,c2) at zs = xs+i ys is given by
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
c     The multipole expansion has 5 set of coefficients
c     vel = \sum mpole(1,k)/z^k + \sum mpole(2,k)/z_bar^k +
c            k                      k
c      
c         z \sum mpole(3,k)/z_bar^k + Re(mpole(4,0) log(z) + \sum mpole(4,k)/z^k)
c              k                                               k
c
c          + i Re(mpole(5,0) log(z) + \sum mpole(5,k)/z^k)
c                                       k
c------------------------------------------------------------------
c     INPUT parameters
c      nd            : number of densities
c      rscale        : scaling parameter
c      sources(2,ns) : coordinates of sources
c      dip(3,ns)    : dipole parameters (corresponding to c2,c3 above)
c      ns           : number of sources
c      center(2)    : expansion center
c      nterms       : Order of multipole expansion
c------------------------------------------------------------------
c      OUTPUT
c      mpole         : coeffs for multipole expansion
c      mpole(1,:)    : coeffs for multipole expansion - laurent
c      mpole(2,:)    : coeffs for multipole expansion - antilaurent 
c      mpole(3,:)    : coeffs for multipole expansion - z antilaurent
c      mpole(4,:)    : coeffs for multipole expansion - logsource1
c      mpole(5,:)    : coeffs for multipole expansion - logsource2
c-----------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns,nterms
      complex *16 mpole(nd,5,0:nterms)
      complex *16 dip(nd,3,ns)
      real *8 sources(2,ns),center(2),rscale
      complex *16 zc,z1,z2,zs,zdis,ztemp,zdisc,ztempc
      complex *16 ztempc2
      complex *16 zt1,zt2,zt3

      ier = 0
      zc = dcmplx(center(1),center(2))
      rinv=1.0d0/rscale
    
      do i=1,ns
         zs=dcmplx(sources(1,i),sources(2,i))
         zdis = (zs-zc)/rscale
         ztemp= zs-zc
         zdisc=dconjg(zdis)

         if(abs(zdis).le.1.0d-16) then
            do idim=1,nd
                mpole(idim,1,1) = mpole(idim,1,1) + dip(idim,1,i)*rinv
                mpole(idim,2,1) = mpole(idim,2,1) + dip(idim,3,i)*rinv
                mpole(idim,3,2) = mpole(idim,3,2) + 
     1              dip(idim,2,i)*rinv**2
            enddo
         endif

         if(abs(zdis).gt.1.0d-16) then
            ztempc=1.0d0/dconjg(ztemp)
            ztempc2=ztempc*ztempc
            do j=1,nterms
               do idim=1,nd
                  zt2 = dip(idim,2,i)*ztempc2
                  zt3 = dip(idim,3,i)*ztempc
c           Expansion corresponding to d1/z-zs
                  mpole(idim,1,j)=mpole(idim,1,j)+
     1               dip(idim,1,i)*zdis/ztemp
c           Expansion corresponding to d2(z-zs)/(z_bar - zs_bar)^2
                  mpole(idim,2,j)=mpole(idim,2,j)-zt2*(j-1)*zdisc*ztemp
                  mpole(idim,3,j)=mpole(idim,3,j)+zt2*(j-1)*zdisc
c           Expansion corresponding to d3/(z_bar - zs_bar)
                  mpole(idim,2,j)=mpole(idim,2,j)+zt3*zdisc
               enddo 
               zdis=zdis*ztemp*rinv
               zdisc=zdisc/ztempc*rinv
            enddo
         endif
      enddo

      return
      end
c
c
c
      subroutine bh2dformmpc(nd,rscale,sources,ns,c1,
     1           center,nterms,mpole)
c*****************************************************************
c     This subroutine computes the multipole expansion about
c     the "center" due to ns sources located at sources(2,*)
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
c     The multipole expansion has 5 set of coefficients
c     vel = \sum mpole(1,k)/z^k + \sum mpole(2,k)/z_bar^k +
c            k                      k
c      
c         z \sum mpole(3,k)/z_bar^k + Re(mpole(4,0) log(z) + \sum mpole(4,k)/z^k)
c              k                                               k
c
c          + i Re(mpole(5,0) log(z) + \sum mpole(5,k)/z^k)
c                                       k
c------------------------------------------------------------------
c     INPUT parameters
c      nd           : number of densities
c      rscale       : scaling parameter
c      sources(2,ns): coordinates of sources
c      c1(2,ns)     : charge strength
c      ns           : number of sources
c      center(2)    : expansion center
c      nterms       : Order of multipole expansion
c------------------------------------------------------------------
c      OUTPUT
c      mpole         : coeffs for multipole expansion
c      mpole(1,:)    : coeffs for multipole expansion - laurent
c      mpole(2,:)    : coeffs for multipole expansion - antilaurent 
c      mpole(3,:)    : coeffs for multipole expansion - z antilaurent
c      mpole(4,:)    : coeffs for multipole expansion - logsource1
c      mpole(5,:)    : coeffs for multipole expansion - logsource2
c-----------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns,nterms
      complex *16 mpole(nd,5,0:nterms)
      complex *16 c1(nd,2,ns)
      real *8 sources(2,ns),center(2),rscale
      complex *16 zc,z1,z2,zs,zdis,ztemp,zdisc,ztempc
      complex *16 ztempc2
      complex *16 zt1,zt2,zt3

      ier = 0
      zc = dcmplx(center(1),center(2))
      rinv=1.0d0/rscale
    
      do i=1,ns
         zs=dcmplx(sources(1,i),sources(2,i))
         zdis = (zs-zc)/rscale
         ztemp= zs-zc
         zdisc=dconjg(zdis)

         if(abs(zdis).le.1.0d-16) then
            do idim=1,nd
                mpole(idim,4,0) = mpole(idim,4,0) + real(2*c1(idim,1,i))
                mpole(idim,5,0) = mpole(idim,5,0) + imag(2*c1(idim,1,i))
c       Expansion corresponding to c2 z-zs/(z_bar - zs_bar)
                mpole(idim,3,1) = mpole(idim,3,1) + 
     1             c1(idim,2,i)*rinv
            enddo
         endif

         if(abs(zdis).gt.1.0d-16) then
            ztempc=1.0d0/dconjg(ztemp)
            ztempc2=ztempc*ztempc
            do idim=1,nd
               mpole(idim,4,0)=mpole(idim,4,0)+dreal(2*c1(idim,1,i))
               mpole(idim,5,0)=mpole(idim,5,0)+dimag(2*c1(idim,1,i))
            enddo
            do j=1,nterms
               do idim=1,nd
                  zt1 = c1(idim,2,i)*ztempc
c          Expansion corresponding to 2 c1 log |z-zs|
                 mpole(idim,4,j)=mpole(idim,4,j)-
     1               dreal(2*c1(idim,1,i))*zdis/j
                 mpole(idim,5,j)=mpole(idim,5,j)-
     1               dimag(2*c1(idim,1,i))*zdis/j
c          Expansion corresponding to c2 z-zs/(z_bar - zs_bar)
                 mpole(idim,2,j)=mpole(idim,2,j)-zt1*zdisc*ztemp
                 mpole(idim,3,j)=mpole(idim,3,j)+zt1*zdisc
               enddo

               zdis=zdis*ztemp*rinv
               zdisc=zdisc/ztempc*rinv
            enddo
         endif
      enddo

      return
      end
c
c
c
c
c
c********************************************************************
      subroutine bh2dformmpcd(nd,rscale,sources,ns,c1,dip,
     1           center,nterms,mpole)
c*****************************************************************
c     This subroutine computes the multipole expansion about
c     the "center" due to ns sources located at sources(2,*)
c
c      The complex velocity at zt = xt + i yt 
c      due to a charge a complex charge c1,c2 at zs = xs+i ys is given by
c
c
c      vel = 2*c1 log|zt-zs| + c2 (zt-zs)/(zt_bar - zs_bar)
c            
c      The complex velocity due to a dipole with complex parameters
c      c2,c3 is given by
c
c
c      vel = d1/(zt-zs) + d2 (zt-zs)/(zt_bar-zs_bar)^2 + 
c             d3/(zt_bar - zs_bar)
c
c     The multipole expansion has 5 set of coefficients
c     vel = \sum mpole(1,k)/z^k + \sum mpole(2,k)/z_bar^k +
c            k                      k
c      
c         z \sum mpole(3,k)/z_bar^k + Re(mpole(4,0) log(z) + \sum mpole(4,k)/z^k)
c              k                                               k
c
c          + i Re(mpole(5,0) log(z) + \sum mpole(5,k)/z^k)
c                                       k
c------------------------------------------------------------------
c     INPUT parameters
c      nd           : number of densities
c      rscale       : scaling parameter
c      sources(2,ns): coordinates of sources
c      c1(2,ns)     : charge strength
c      dip(3,ns)    : dipole parameters (corresponding to c2,c3 above)
c      ns           : number of sources
c      center(2)    : expansion center
c      nterms       : Order of multipole expansion
c------------------------------------------------------------------
c      OUTPUT
c      mpole         : coeffs for multipole expansion
c      mpole(1,:)    : coeffs for multipole expansion - laurent
c      mpole(2,:)    : coeffs for multipole expansion - antilaurent 
c      mpole(3,:)    : coeffs for multipole expansion - z antilaurent
c      mpole(4,:)    : coeffs for multipole expansion - logsource1
c      mpole(5,:)    : coeffs for multipole expansion - logsource2
c-----------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns,nterms
      complex *16 mpole(nd,5,0:nterms)
      complex *16 c1(nd,2,ns),dip(nd,3,ns)
      real *8 sources(2,ns),center(2),rscale
      complex *16 zc,z1,z2,zs,zdis,ztemp,zdisc,ztempc
      complex *16 ztempc2
      complex *16 zt1,zt2,zt3

      ier = 0
      zc = dcmplx(center(1),center(2))
      rinv=1.0d0/rscale
    
      do i=1,ns
         zs=dcmplx(sources(1,i),sources(2,i))
         zdis = (zs-zc)/rscale
         ztemp= zs-zc
         zdisc=dconjg(zdis)

         if(abs(zdis).le.1.0d-16) then
            do idim=1,nd
                mpole(idim,4,0) = mpole(idim,4,0) + real(2*c1(idim,1,i))
                mpole(idim,5,0) = mpole(idim,5,0) + imag(2*c1(idim,1,i))
c       Expansion corresponding to c2 z-zs/(z_bar - zs_bar)
                mpole(idim,3,1) = mpole(idim,3,1) + 
     1             c1(idim,2,i)*rinv
                mpole(idim,1,1) = mpole(idim,1,1) + dip(idim,1,i)*rinv
                mpole(idim,2,1) = mpole(idim,2,1) + dip(idim,3,i)*rinv
                mpole(idim,3,2) = mpole(idim,3,2) + 
     1              dip(idim,2,i)*rinv**2
            enddo
         endif

         if(abs(zdis).gt.1.0d-16) then
            ztempc=1.0d0/dconjg(ztemp)
            ztempc2=ztempc*ztempc
            do idim=1,nd
               mpole(idim,4,0)=mpole(idim,4,0)+dreal(2*c1(idim,1,i))
               mpole(idim,5,0)=mpole(idim,5,0)+dimag(2*c1(idim,1,i))
            enddo
            do j=1,nterms
               do idim=1,nd
                  zt1 = c1(idim,2,i)*ztempc
c          Expansion corresponding to 2 c1 log |z-zs|
                 mpole(idim,4,j)=mpole(idim,4,j)-
     1               dreal(2*c1(idim,1,i))*zdis/j
                 mpole(idim,5,j)=mpole(idim,5,j)-
     1               dimag(2*c1(idim,1,i))*zdis/j
c          Expansion corresponding to c2 z-zs/(z_bar - zs_bar)
                 mpole(idim,2,j)=mpole(idim,2,j)-zt1*zdisc*ztemp
                 mpole(idim,3,j)=mpole(idim,3,j)+zt1*zdisc
               enddo

               do idim=1,nd
                  zt2 = dip(idim,2,i)*ztempc2
                  zt3 = dip(idim,3,i)*ztempc
c           Expansion corresponding to d1/z-zs
                  mpole(idim,1,j)=mpole(idim,1,j)+
     1               dip(idim,1,i)*zdis/ztemp
c           Expansion corresponding to d2(z-zs)/(z_bar - zs_bar)^2
                  mpole(idim,2,j)=mpole(idim,2,j)-zt2*(j-1)*zdisc*ztemp
                  mpole(idim,3,j)=mpole(idim,3,j)+zt2*(j-1)*zdisc
c           Expansion corresponding to d3/(z_bar - zs_bar)
                  mpole(idim,2,j)=mpole(idim,2,j)+zt3*zdisc
               enddo 
               zdis=zdis*ztemp*rinv
               zdisc=zdisc/ztempc*rinv
            enddo
         endif
      enddo

      return
      end
c
c
c
c
c
c
c********************************************************************
      subroutine bh2dformtac(nd,rscale,sources,ns,c1,
     1       center,nterms,mpole)
c*****************************************************************
c     This subroutine computes the local expansion about
c     the "center" due to ns sources located at sources(2,*)
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
c     The multipole expansion has 5 set of coefficients
c     vel = \sum mp1(k) z^k + \sum mp2(k) z_bar^k +
c            k                k
c      
c         z \sum mp3(k) z_bar^k + Re(mp4(0) log(z) + \sum mp4(k) z^k)
c              k                                       k
c
c          + i Re(mp5(0) log(z) + \sum mp5(k) z^k)
c                                  k
c------------------------------------------------------------------
c     INPUT parameters
c
c      rscale        : scaling parameter
c      sources(2,ns) : coordinates of sources
c      c1(2,ns)      : charge strength
c      ns            : number of sources
c      center(2)     : expansion center
c      nterms        : Order of multipole expansion
c------------------------------------------------------------------
c      OUTPUT
c      mpole         : coeffs for local expansion
c      mpole(1,:)    : coeffs for multipole expansion - laurent
c      mpole(2,:)    : coeffs for multipole expansion - antilaurent 
c      mpole(3,:)    : coeffs for multipole expansion - z antilaurent
c      mpole(4,:)    : coeffs for multipole expansion - logsource1
c      mpole(5,:)    : coeffs for multipole expansion - logsource2
c-----------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns,nterms
      complex *16 mpole(nd,5,0:nterms)
      complex *16 c1(nd,2,ns)
      real *8 sources(2,ns),center(2),rscale
      complex *16 zc,z1,z2,zs,zdis,ztemp,zdisc,ztempc

      zc = dcmplx(center(1),center(2))

      do i=1,ns
         zs=dcmplx(sources(1,i),sources(2,i))
         zdis = 1.0d0
         ztemp= 1.0d0/(zs-zc)
         zdisc=dconjg(zdis)
         ztempc=dconjg(ztemp)
         do j=0,nterms
            do idim=1,nd
              if(j.eq.0) then
                 mpole(idim,4,j)=mpole(idim,4,j) + 
     1               dreal(2*c1(idim,1,i))*log(cdabs(1.0d0/ztemp))
                 mpole(idim,5,j)=mpole(idim,5,j) + 
     1               dimag(2*c1(idim,1,i))*log(cdabs(1.0d0/ztemp))
              else
                mpole(idim,4,j)=mpole(idim,4,j) - 
     1              dreal(2*c1(idim,1,i))*zdis/j
                mpole(idim,5,j)=mpole(idim,5,j) - 
     1              dimag(2*c1(idim,1,i))*zdis/j
              endif
              mpole(idim,2,j)=mpole(idim,2,j) + 
     1            c1(idim,2,i)*zdisc*ztempc/ztemp
              mpole(idim,3,j)=mpole(idim,3,j) - 
     1           c1(idim,2,i)*zdisc*ztempc
 
            enddo
            zdis=zdis*ztemp*rscale
            zdisc=zdisc*ztempc*rscale
         enddo
      enddo

      return
      end

c
c
c
c
c********************************************************************
      subroutine bh2dformtad(nd,rscale,sources,ns,dip,
     1       center,nterms,mpole)
c*****************************************************************
c     This subroutine computes the local expansion about
c     the "center" due to ns sources located at sources(2,*)
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
c     The multipole expansion has 5 set of coefficients
c     vel = \sum mp1(k) z^k + \sum mp2(k) z_bar^k +
c            k                k
c      
c         z \sum mp3(k) z_bar^k + Re(mp4(0) log(z) + \sum mp4(k) z^k)
c              k                                       k
c
c          + i Re(mp5(0) log(z) + \sum mp5(k) z^k)
c                                  k
c------------------------------------------------------------------
c     INPUT parameters
c
c      rscale        : scaling parameter
c      sources(2,ns) : coordinates of sources
c      dip(3,n s)    : dipole parameters (corresponding to c2,c3 above)
c      ns            : number of sources
c      center(2)     : expansion center
c      nterms        : Order of multipole expansion
c------------------------------------------------------------------
c      OUTPUT
c      mpole         : coeffs for local expansion
c      mpole(1,:)    : coeffs for multipole expansion - laurent
c      mpole(2,:)    : coeffs for multipole expansion - antilaurent 
c      mpole(3,:)    : coeffs for multipole expansion - z antilaurent
c      mpole(4,:)    : coeffs for multipole expansion - logsource1
c      mpole(5,:)    : coeffs for multipole expansion - logsource2
c-----------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns,nterms
      complex *16 mpole(nd,5,0:nterms)
      complex *16 dip(nd,3,ns)
      real *8 sources(2,ns),center(2),rscale
      complex *16 zc,z1,z2,zs,zdis,ztemp,zdisc,ztempc

      zc = dcmplx(center(1),center(2))

      do i=1,ns
         zs=dcmplx(sources(1,i),sources(2,i))
         zdis = 1.0d0
         ztemp= 1.0d0/(zs-zc)
         zdisc=dconjg(zdis)
         ztempc=dconjg(ztemp)
         do j=0,nterms
            do idim=1,nd
              mpole(idim,1,j)=mpole(idim,1,j)-dip(idim,1,i)*zdis*ztemp
              mpole(idim,2,j)=mpole(idim,2,j)-dip(idim,3,i)*zdisc*ztempc
              mpole(idim,2,j)=mpole(idim,2,j)-
     1             dip(idim,2,i)*(j+1)*zdisc*ztempc*ztempc/ztemp
              mpole(idim,3,j)=mpole(idim,3,j) + 
     1            dip(idim,2,i)*(j+1)*zdisc*ztempc*ztempc
            
            enddo
            zdis=zdis*ztemp*rscale
            zdisc=zdisc*ztempc*rscale
         enddo
      enddo

      return
      end

c
c
c
c
c********************************************************************
      subroutine bh2dformtacd(nd,rscale,sources,ns,c1,dip,
     1       center,nterms,mpole)
c*****************************************************************
c     This subroutine computes the local expansion about
c     the "center" due to ns sources located at sources(2,*)
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
c     The multipole expansion has 5 set of coefficients
c     vel = \sum mp1(k) z^k + \sum mp2(k) z_bar^k +
c            k                k
c      
c         z \sum mp3(k) z_bar^k + Re(mp4(0) log(z) + \sum mp4(k) z^k)
c              k                                       k
c
c          + i Re(mp5(0) log(z) + \sum mp5(k) z^k)
c                                  k
c------------------------------------------------------------------
c     INPUT parameters
c
c      rscale       : scaling parameter
c      sources(2,ns): coordinates of sources
c      c1(2,ns)     : charge strength
c      d1(3,ns)     : dipole strenghts
c      ns           : number of sources
c      center(2)    : expansion center
c      nterms       : Order of multipole expansion
c------------------------------------------------------------------
c      OUTPUT
c      mpole         : coeffs for local expansion
c      mpole(1,:)    : coeffs for multipole expansion - laurent
c      mpole(2,:)    : coeffs for multipole expansion - antilaurent 
c      mpole(3,:)    : coeffs for multipole expansion - z antilaurent
c      mpole(4,:)    : coeffs for multipole expansion - logsource1
c      mpole(5,:)    : coeffs for multipole expansion - logsource2
c-----------------------------------------------------------------

      implicit real *8 (a-h,o-z)
      integer ns,nterms
      complex *16 mpole(nd,5,0:nterms)
      complex *16 c1(nd,2,ns),dip(nd,3,ns)
      real *8 sources(2,ns),center(2),rscale
      complex *16 zc,z1,z2,zs,zdis,ztemp,zdisc,ztempc

      zc = dcmplx(center(1),center(2))

      do i=1,ns
         zs=dcmplx(sources(1,i),sources(2,i))
         zdis = 1.0d0
         ztemp= 1.0d0/(zs-zc)
         zdisc=dconjg(zdis)
         ztempc=dconjg(ztemp)
         do j=0,nterms
            do idim=1,nd
              if(j.eq.0) then
                 mpole(idim,4,j)=mpole(idim,4,j) + 
     1               dreal(2*c1(idim,1,i))*log(cdabs(1.0d0/ztemp))
                 mpole(idim,5,j)=mpole(idim,5,j) + 
     1               dimag(2*c1(idim,1,i))*log(cdabs(1.0d0/ztemp))
              else
                mpole(idim,4,j)=mpole(idim,4,j) - 
     1              dreal(2*c1(idim,1,i))*zdis/j
                mpole(idim,5,j)=mpole(idim,5,j) - 
     1              dimag(2*c1(idim,1,i))*zdis/j
              endif
              mpole(idim,2,j)=mpole(idim,2,j) + 
     1            c1(idim,2,i)*zdisc*ztempc/ztemp
              mpole(idim,3,j)=mpole(idim,3,j) - 
     1           c1(idim,2,i)*zdisc*ztempc
 
              mpole(idim,1,j)=mpole(idim,1,j)-dip(idim,1,i)*zdis*ztemp
              mpole(idim,2,j)=mpole(idim,2,j)-dip(idim,3,i)*zdisc*ztempc
              mpole(idim,2,j)=mpole(idim,2,j)-
     1             dip(idim,2,i)*(j+1)*zdisc*ztempc*ztempc/ztemp
              mpole(idim,3,j)=mpole(idim,3,j) + 
     1            dip(idim,2,i)*(j+1)*zdisc*ztempc*ztempc
            
            enddo
            zdis=zdis*ztemp*rscale
            zdisc=zdisc*ztempc*rscale
         enddo
      enddo

      return
      end

c********************************************************************
c     TRANSLATION OPERATORS
c*******************************************************************

      subroutine bh2dlocloc(nd,rscale1,c1,hexp,nterms1,
     1       rscale2,c2,jexp,nterms2,carray,ldc)
c******************************************************************
c     Converts local expansion to local expansion
c     Given original expansions and precomputed binomial coefficients
c
c     Given local expansion in the following form
c     vel = \sum h1(k) z1^k + \sum h2(k) z1_bar^k +
c            k                k
c      
c     z1 \sum h3(k) z1_bar^k + Re(h4(0) log(z1) + \sum h4(k) z1^k)
c              k                                       k
c
c          + i Re(h5(0) log(z1) + \sum h5(k) z1^k)
c     z1 = z-c1
c
c     NOTE: The subroutine will work if nterms<1000, as 1000, is
c     a hardwired number for computing zpow in the computation

c-----------------------------------------------------------------
c     INPUT
c     rscale1 : scaling parameter for original expansion 
c     c1      : center of original local expansion
c     hexp    : coefficients of original taylor expansion
c     hexp(1,:) : coefficients of original taylor expansion
c     hexp(2,:) : coefficients of original antitaylor expansion
c     hexp(3,:) : coefficients of original z1 (antitaylor expansion)
c     hexp(4,:) : coefficients of original real part logsource1
c     hexp(5,:) : coefficients of original imag part logsource2
c     nterms1 : number of terms in original expansion
c     rscale2 : scaling parameter for shifted expansion
c     c2      : center of new local expansion
c     nterms2 : number of terms of the shifted expansion
c     carray  : Precomputed binomial coefficients
c     ldc     : size of carray
c-----------------------------------------------------------------
c     OUTPUT
c     jexp    : coefficients of shifted taylor expansion
c     jexp(1,:) : coefficients of shifted taylor expansion
c     jexp(2,:) : coefficients of shifted antitaylor expansion
c     jexp(3,:) : coefficients of shifted z1 (antitaylor expansion)
c     jexp(4,:) : coefficients of shifted real part logsource1
c     jexp(5,:) : coefficients of shifted imag part logsource2
c---------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms1,nterms2,ldc
      complex *16 hexp(nd,5,0:nterms1),jexp(nd,5,0:nterms2)
      complex *16 h1a(nd,0:nterms1),h2a(nd,0:nterms1),h3a(nd,0:nterms1)
      complex *16 h4a(nd,0:nterms1),h5a(nd,0:nterms1)
      complex *16 j1(nd,0:nterms2),j2(nd,0:nterms2),j3(nd,0:nterms2)
      complex *16 j4(nd,0:nterms2),j5(nd,0:nterms2)
      complex *16 zc1,zc2,zpow1(0:1000),zpow2(0:1000),zdis,zdisinv
      real *8 carray(0:ldc,0:ldc),rscale1,rscale2,c1(2),c2(2)

      zc1 = dcmplx(c1(1),c1(2))
      zc2 = dcmplx(c2(1),c2(2))
      zdis = zc2-zc1
      zdisinv = 1.0d0/zdis
      zpow1(0)=1.0d0
      zpow2(0)=1.0d0
      rinv1 = 1.0d0/rscale1

      do i=0,nterms2
        do idim=1,nd
          j1(idim,i) = 0
          j2(idim,i) = 0
          j3(idim,i) = 0
          j4(idim,i) = 0
          j5(idim,i) = 0
        enddo
      enddo
      do i=1,max(nterms1,nterms2)
           zpow1(i)=zpow1(i-1)*zdis*rinv1
           zpow2(i)=zpow2(i-1)*rscale2*zdisinv
      enddo
      
      do i=0,nterms1
        do idim=1,nd
          h1a(idim,i)=hexp(idim,1,i)*zpow1(i)
          h2a(idim,i)=(hexp(idim,2,i)+
     1        hexp(idim,3,i)*zdis)*dconjg(zpow1(i))
          h3a(idim,i)=hexp(idim,3,i)*dconjg(zpow1(i))
          h4a(idim,i)=hexp(idim,4,i)*zpow1(i)
          h5a(idim,i)=hexp(idim,5,i)*zpow1(i)
        enddo
      enddo

      do i=0,nterms2
         do j=i,nterms1
           do idim=1,nd
             j1(idim,i)=j1(idim,i)+h1a(idim,j)*carray(j,i)
             j2(idim,i)=j2(idim,i)+h2a(idim,j)*carray(j,i)
             j3(idim,i)=j3(idim,i)+h3a(idim,j)*carray(j,i)
             j4(idim,i)=j4(idim,i)+h4a(idim,j)*carray(j,i)
             j5(idim,i)=j5(idim,i)+h5a(idim,j)*carray(j,i)
           enddo
        enddo
        do idim=1,nd
          j1(idim,i)=j1(idim,i)*zpow2(i)
          j2(idim,i)=j2(idim,i)*dconjg(zpow2(i))
          j3(idim,i)=j3(idim,i)*dconjg(zpow2(i))
          j4(idim,i)=j4(idim,i)*zpow2(i)
          j5(idim,i)=j5(idim,i)*zpow2(i)
          jexp(idim,1,i) = jexp(idim,1,i) + j1(idim,i)
          jexp(idim,2,i) = jexp(idim,2,i) + j2(idim,i)
          jexp(idim,3,i) = jexp(idim,3,i) + j3(idim,i)
          jexp(idim,4,i) = jexp(idim,4,i) + j4(idim,i)
          jexp(idim,5,i) = jexp(idim,5,i) + j5(idim,i)
        enddo 
      enddo

      return
      end
c******************************************************************
      subroutine bh2dmpmp(nd,rscale1,c1,hexp,nterms1,
     1       rscale2,c2,jexp,nterms2,carray,ldc)
c******************************************************************
c     Converts multipole expansion to multipole expansion
c     Given original expansions and precomputed binomial coefficients
c
c     Given multipole expansion in the following form
c     vel = \sum h1(k) /z1^k + \sum h2(k) /z1_bar^k +
c            k                k
c      
c     z1 \sum h3(k) /z1_bar^k + Re(h4(0) log(z1) + \sum h4(k) /z1^k)
c              k                                       k
c
c          + i Re(h5(0) log(z1) + \sum h5(k) /z1^k)
c     z1 = z-c1
c
c     NOTE: The subroutine will work if nterms<1000, as 1000, is
c     a hardwired number for computing zpow in the computation

c-----------------------------------------------------------------
c     INPUT
c     nd      : number of densities
c     rscale1 : scaling parameter of original multipole expansion
c     c1      : center of original mutlipole expansion
c     hexp    : coefficients of original multipole expansion
c     hexp(1,:) : coefficients of original laurent expansion
c     hexp(2,:) : coefficients of original antilaurent expansion
c     hexp(3,:) : coefficients of original z1 *antilaurent expansion
c     hexp(4,:) : coefficients of original real part logsource1
c     hexp(5,:) : coefficients of original imag part logsource2
c     nterms1 : number of terms in original expansion
c     rscale2 : scaling parameter of shifted multipole expansion
c     c2      : center of new multipole expansion
c     nterms2 : number of terms of the shifted expansion
c     carray  : precomputed binomial coefficieints array
c-----------------------------------------------------------------
c     OUTPUT
c     jexp    : coefficients of shifted multipole expansion
c     jexp(1,:) : coefficients of shifted laurent expansion
c     jexp(2,:) : coefficients of shifted antilaurent expansion
c     jexp(3,:) : coefficients of shifted z1 *antilaurent expansion
c     jexp(4,:) : coefficients of shifted real part logsource1
c     jexp(5,:) : coefficients of shifted imag part logsource2
c---------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms1,nterms2,ldc
      complex *16 hexp(nd,5,0:nterms1),jexp(nd,5,0:nterms2) 
      complex *16 h1a(nd,0:nterms1),h2a(nd,0:nterms1),h3a(nd,0:nterms1)
      complex *16 h4a(nd,0:nterms1),h5a(nd,0:nterms1)
      complex *16 j1(nd,0:nterms2),j2(nd,0:nterms2),j3(nd,0:nterms2)
      complex *16 j4(nd,0:nterms2),j5(nd,0:nterms2)
      complex *16 zc1,zc2,zpow1(0:1000),zdis,zdisinv
      complex *16 zpow2(0:1000)
      real *8 carray(0:ldc,0:ldc)
      real *8 rscale1,rscale2,c1(2),c2(2),rinv2

      nmax=max(nterms1,nterms2)
      
      zc1 = dcmplx(c1(1),c1(2))
      zc2 = dcmplx(c2(1),c2(2))
      zdis = zc1-zc2
      zdisinv = 1.0d0/zdis
      rinv2 = 1.0d0/rscale2
      zpow1(0)=1.0d0
      zpow2(0)=1.0d0
      do i=1,nmax
           zpow1(i)=zpow1(i-1)*zdisinv*rscale1
           zpow2(i)=zpow2(i-1)*zdis*rinv2
      enddo
      
      do i=0,nterms2
        do idim=1,nd
         j1(idim,i)=0.0d0
         j2(idim,i)=0.0d0
         j3(idim,i)=0.0d0
         j4(idim,i)=0.0d0
         j5(idim,i)=0.0d0
       enddo
      enddo

c     Handling the log term in the expansion
      do idim=1,nd
        jexp(idim,4,0)=jexp(idim,4,0) + hexp(idim,4,0)
        jexp(idim,5,0)=jexp(idim,5,0) + hexp(idim,5,0)
      enddo
      do i=1,nterms2
        do idim=1,nd
          j4(idim,i)=j4(idim,i)-hexp(idim,4,0)/i
          j5(idim,i)=j5(idim,i)-hexp(idim,5,0)/i
        enddo
      enddo

      do i=1,min(nterms2,nterms1)
        do idim=1,nd
          h1a(idim,i)=hexp(idim,1,i)*zpow1(i)
          h2a(idim,i)=(hexp(idim,2,i)-
     1        hexp(idim,3,i)*zdis)*dconjg(zpow1(i))
          h3a(idim,i)=hexp(idim,3,i)*dconjg(zpow1(i))
          h4a(idim,i)=hexp(idim,4,i)*zpow1(i)
          h5a(idim,i)=hexp(idim,5,i)*zpow1(i)
        enddo
      enddo

      do i=1,nterms2
         do j=1,min(i,nterms1)
           do idim=1,nd
             j1(idim,i)=j1(idim,i)+h1a(idim,j)*carray(i-1,j-1)
             j2(idim,i)=j2(idim,i)+h2a(idim,j)*carray(i-1,j-1)
             j3(idim,i)=j3(idim,i)+h3a(idim,j)*carray(i-1,j-1)
             j4(idim,i)=j4(idim,i)+h4a(idim,j)*carray(i-1,j-1)
             j5(idim,i)=j5(idim,i)+h5a(idim,j)*carray(i-1,j-1)
           enddo
        enddo
        do idim=1,nd
          j1(idim,i)=j1(idim,i)*zpow2(i)
          j2(idim,i)=j2(idim,i)*dconjg(zpow2(i))
          j3(idim,i)=j3(idim,i)*dconjg(zpow2(i))
          j4(idim,i)=j4(idim,i)*zpow2(i)
          j5(idim,i)=j5(idim,i)*zpow2(i)
          jexp(idim,1,i) = jexp(idim,1,i) + j1(idim,i)
          jexp(idim,2,i) = jexp(idim,2,i) + j2(idim,i)
          jexp(idim,3,i) = jexp(idim,3,i) + j3(idim,i)
          jexp(idim,4,i) = jexp(idim,4,i) + j4(idim,i)
          jexp(idim,5,i) = jexp(idim,5,i) + j5(idim,i)
        enddo
      enddo

      return
      end
c*******************************************************************
      subroutine bh2dmploc(nd,rscale1,c1,hexp,nterms1,
     1       rscale2,c2,jexp,nterms2,carray,ldc)
c******************************************************************
c     Converts multipole expansion to local expansion
c     Given original expansions and precomputed binomial coeffients
c
c     Given multipole expansion in the following form
c     vel = \sum h1(k) /z1^k + \sum h2(k) /z1_bar^k +
c            k                k
c      
c     z1 \sum h3(k) /z1_bar^k + Re(h4(0) log(z1) + \sum h4(k) /z1^k)
c              k                                       k
c
c          + i Re(h5(0) log(z1) + \sum h5(k) /z1^k)
c     z1 = z-c1
c
c     NOTE: The subroutine will work if nterms<1000, as 1000, is
c     a hardwired number for computing zpow in the computation

c-----------------------------------------------------------------
c     INPUT
c     nd      : number of densities
c     rscale1 : scaling parameter of original multipole expansion 
c     c1      : center of original mutlipole expansion
c     hexp    : coefficients of original multipole expansion
c     hexp(1,:) : coefficients of original laurent expansion
c     hexp(2,:) : coefficients of original antilaurent expansion
c     hexp(3,:) : coefficients of original z1 *antilaurent expansion
c     hexp(4,:) : coefficients of original real part logsource1
c     hexp(5,:) : coefficients of original imag part logsource2
c     nterms1 : number of terms of original expansion
c     rscale2 : scaling parameter of shifted local expansion
c     c2      : center of new local expansion
c     nterms2 : number of terms of the shifted expansion
c-----------------------------------------------------------------
c     OUTPUT
c     jexp    : coefficients of shifted taylor expansion
c     jexp(1,:) : coefficients of shifted taylor expansion
c     jexp(2,:) : coefficients of shifted antaylor expansion
c     jexp(3,:) : coefficients of shifted z1 *antitaylor expansion
c     jexp(4,:) : coefficients of shifted real part logsource1
c     jexp(5,:) : coefficients of shifted imag part logsource2
c---------------------------------------------------------------
      implicit real *8 (a-h,o-z)
      integer nterms1,nterms2,ldc
      complex *16 hexp(nd,5,0:nterms1),jexp(nd,5,0:nterms2)
      complex *16 h1a(nd,0:nterms1),h2a(nd,0:nterms1),h3a(nd,0:nterms1)
      complex *16 h4a(nd,0:nterms1),h5a(nd,0:nterms1)
      complex *16 j1(nd,0:nterms2),j2(nd,0:nterms2),j3(nd,0:nterms2)
      complex *16 j4(nd,0:nterms2),j5(nd,0:nterms2)
      complex *16 zc1,zc2,zdis,zdis1
      complex *16 zpow1(0:1000)
      complex *16 zzz,zzz1
      complex *16 zpow2(0:1000)
      real *8 carray(0:ldc,0:ldc)
      real *8 rscale1,rscale2,c1(2),c2(2)

      nmax=nterms1+nterms2
      
      zc1 = dcmplx(c1(1),c1(2))
      zc2 = dcmplx(c2(1),c2(2))
      zdis = 1/(zc1-zc2)
      zdis1 = (zc1-zc2)
      zpow1(0)=1.0d0
      zpow2(0)=1.0d0
      do i=1,nmax
           zpow1(i)=-zpow1(i-1)*zdis*rscale1
           zpow2(i)=zpow2(i-1)*zdis*rscale2
      enddo
      
      do i=0,nterms2
        do idim=1,nd
         j1(idim,i)=0.0d0
         j2(idim,i)=0.0d0
         j3(idim,i)=0.0d0
         j4(idim,i)=0.0d0
         j5(idim,i)=0.0d0
        enddo
      enddo

c     Handling the log term in the expansion
      do idim=1,nd
        j4(idim,0)=hexp(idim,4,0)*log(cdabs(zdis1))
        j5(idim,0)=hexp(idim,5,0)*log(cdabs(zdis1))
      enddo
      do i=1,nterms2
        do idim=1,nd
          j4(idim,i)=j4(idim,i)-hexp(idim,4,0)/i
          j5(idim,i)=j5(idim,i)-hexp(idim,5,0)/i
        enddo
      enddo
c
      do i=1,nterms1
        do idim=1,nd
          h1a(idim,i)=hexp(idim,1,i)*zpow1(i)
          h2a(idim,i)=(hexp(idim,2,i)-
     1        hexp(idim,3,i)*zdis1)*dconjg(zpow1(i))
          h3a(idim,i)=hexp(idim,3,i)*dconjg(zpow1(i))
          h4a(idim,i)=hexp(idim,4,i)*zpow1(i)
          h5a(idim,i)=hexp(idim,5,i)*zpow1(i)
        enddo
      enddo
      do i=0,nterms2
         do j=1,nterms1
           do idim=1,nd
             j1(idim,i)=j1(idim,i)+h1a(idim,j)*carray(i+j-1,i)
             j2(idim,i)=j2(idim,i)+h2a(idim,j)*carray(i+j-1,i)
             j3(idim,i)=j3(idim,i)+h3a(idim,j)*carray(i+j-1,i)
             j4(idim,i)=j4(idim,i)+h4a(idim,j)*carray(i+j-1,i)
             j5(idim,i)=j5(idim,i)+h5a(idim,j)*carray(i+j-1,i)
           enddo
        enddo
        do idim=1,nd
          j1(idim,i)=j1(idim,i)*zpow2(i)
          j2(idim,i)=j2(idim,i)*dconjg(zpow2(i))
          j3(idim,i)=j3(idim,i)*dconjg(zpow2(i))
          j4(idim,i)=j4(idim,i)*zpow2(i)
          j5(idim,i)=j5(idim,i)*zpow2(i)
          jexp(idim,1,i) = jexp(idim,1,i) + j1(idim,i)
          jexp(idim,2,i) = jexp(idim,2,i) + j2(idim,i)
          jexp(idim,3,i) = jexp(idim,3,i) + j3(idim,i)
          jexp(idim,4,i) = jexp(idim,4,i) + j4(idim,i)
          jexp(idim,5,i) = jexp(idim,5,i) + j5(idim,i)
        enddo
      enddo

      return
      end
c*******************************************************************

      subroutine bh2dmpzero(nd,mpole,nterms)
      implicit real *8 (a-h,o-z)
      complex *16 mpole(nd,5,0:nterms)
      
      do i=0,nterms
        do j=1,5
          do idim=1,nd
            mpole(idim,j,i) = 0
          enddo
        enddo
      enddo

      return
      end
