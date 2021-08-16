cc Copyright (C) 2017: Travis Askham, Leslie Greengard,
cc Zydrunas Gimbutas
cc email: askhamwhat@gmail.com      
cc 
cc This software is being released under a modified FreeBSD license
cc (see licenses folder in home directory). 

c      This file contains the basic subroutines for 
c      forming and evaluating multipole (partial wave) expansions
c      in two dimensions.
c
c-----------------------------------------------------------------------
c         **************  IMPORTANT NOTE **************  
c
c      THE CONVENTION IN THIS SUBROUTINE LIBRARY IS THAT ALL COMPUTATIONS
C      USING EXPANSIONS ARE INCREMENTAL. 
c      That is, the multipole and local expansions 
c      on output are incremented from their input values, not overwritten.
c      The same is true for translation routines.
c
c-----------------------------------------------------------------------
c
c      MBH2DCONVTOMP_VEC: converts charge, dipole, quadrupole, octopole
c      			  sources to multipolar
c
c      MBH2DFORMMPMP_VEC:  creates multipole expansions (outgoing) due to 
c                       a collection of multipolar sources.
c      MBH2DMPEVALP_VEC:  computes potentials from multipole expansions
c                       at a collection of targets
c      MBH2DMPEVALG_VEC:  computes potentials/gradients from multipole
c                          expansions at a collection of targets
c      MBH2DMPEVALH_VEC:  computes potentials/gradients/Hessians 
c                         from multipole expansions at a collection of
c                         targets
c      MBH2DFORMTAMP_VEC:  creates local expansions due to 
c                       a collection of multipolar sources
c      MBH2DTAEVALP_VEC:  computes potentials from local expansions
c                       at a collection of targets
c      MBH2DTAEVALG_VEC:  computes potentials/gradients from local expansions
c                       at a collection of targets
c      MBH2DTAEVALH_VEC:  computes potentials/gradients/Hessians 
c                       from local expansions at a collection of targets
c 
c      MBH2DMPMP_VEC:     Translates center of multipole expansions 
c      MBH2DLOCLOC_VEC:   Translates center of local expansions 
c      MBH2DMPLOC_VEC:    Converts multipole expansions to local expansions
c
c     MBH2DMPZERO_VEC:   utility to initialize multipole coefficients
c                           to zero.
C
c-----------------------------------------------------------------------


      subroutine mbh2dconvtomp_vec(nd,beta,ns,ifcharge,charge,
     1     ifdipole,dipstr,dipvec,ifquad,quadstr,quadvec,ifoct,octstr,
     2     octvec,nterms,mbhmpole,ympole)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     Converts up to octopole strength sources to a short, equivalent
c     multipolar expansion
c
c     input:
c
c     nd - integer, number of expansions at each source location
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     ifcharge - integer, flag, equals 1 if charges are present
c     ifdipole - integer, flag, equals 1 if dipoles are present
c     ifquad - integer, flag, equals 1 if quadrupoles are present
c     ifoct - integer, flag, equals 1 if octopoles are present            
c     charge - real *8 (ns) charge strengths
c     nterms - integer, 0:nterms is the second dimension of
c                the output mbhmpole,ympole arrays
c                requirements:
c            for charge nterms >= 0
c            for dipole nterms >= 1
c            for quadrupoles nterms >= 2
c            for octopoles nterms >= 3      
c     
c     output:
c
c     NOTE: these arrays are *incremented* on output
c
c     mbhmpole - complex *16 (nd,0:nterms,ns), difference kernel part
c                 of equivalent multipolar expansions
c     ympole - complex *16 (nd,0:nterms,ns), Yukawa part of equivalent
c                 multipolar expansions
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ns, nterms, ifcharge, ifdipole, ifquad, ifoct, nd
      real *8 beta, charge(nd,*), dipstr(nd,*)
      real *8 dipvec(nd,2,*)
      real *8 quadstr(nd,*), quadvec(nd,3,*), octstr(nd,*)
      real *8 octvec(nd,4,*)
      complex *16 mbhmpole(nd,0:nterms,*), ympole(nd,0:nterms,*)
c     local variables
      integer i, j
      complex *16 ima
      real *8 pi
c     
      data ima/(0.0d0,1.0d0)/
c     

      pi = 4.0d0*datan(1.0d0)

      do i = 1,ns

         if (ifcharge .eq. 1) then
            do j = 1,nd
               mbhmpole(j,0,i) = mbhmpole(j,0,i) +
     1              charge(j,i)/(2.0d0*pi*beta**2)
            enddo
         endif

         if (ifdipole .eq. 1) then
            do j = 1,nd
               mbhmpole(j,1,i) = mbhmpole(j,1,i) +
     1              dipstr(j,i)*(dipvec(j,1,i) +
     1              ima*dipvec(j,2,i))/(2.0d0*pi*beta)
            enddo
         endif

         if (ifquad .eq. 1) then
            do j = 1,nd
               mbhmpole(j,2,i) = mbhmpole(j,2,i) +
     1              quadstr(j,i)*(quadvec(j,1,i) +
     1              quadvec(j,2,i)*ima -
     2              quadvec(j,3,i))/(4.0d0*pi)
               ympole(j,0,i) = ympole(j,0,i) +
     1              quadstr(j,i)*(quadvec(j,1,i) +
     1              quadvec(j,3,i))/(4.0d0*pi)
            enddo
         endif

         if (ifoct .eq. 1) then
            do j = 1,nd
               mbhmpole(j,3,i) = mbhmpole(j,3,i) +
     1              beta*octstr(j,i)*(octvec(j,1,i) +
     1              ima*octvec(j,2,i) - octvec(j,3,i) -
     2              ima*octvec(j,4,i))/(8.0d0*pi)
               ympole(j,1,i) = ympole(j,1,i) +
     1              beta*octstr(j,i)*(3.0d0*octvec(j,1,i) +
     1              ima*octvec(j,2,i) + octvec(j,3,i) +
     2              3.0d0*ima*octvec(j,4,i))/(8.0d0*pi)
            enddo
         endif
         
      enddo
      

      return
      end


      

      subroutine mbh2dmpevalp_vec(nd,beta,rscale,center,
     1     mbhmpole,ympole,nterms,ztarg,ntarg,
     2     pot)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor of the expansions
c     center(2)          : expansion center
c     mbhmpole(0:nterms) : coefficients for difference-type functions
c     ympole(nterms)     : coefficients for modified Bessel functions
c     nterms             : number of terms in expansions
c     ztarg(2,ntarg)     : coordinates of targets
c     ntarg              : number of targets
c
c     OUTPUT:
c      
c     pot(ntarg)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhmpole(nd,0:nterms), ympole(nd,0:nterms)
      real *8 rscale, center(2), ztarg(2,ntarg), beta
      real *8 pot(nd,ntarg)
      integer nterms, ntarg, nd
c     local
      real *8 zdiff(2), r, theta
      real *8 ders(0:1), pih
      real *8, allocatable :: diffs(:), kvec(:)
      complex *16 z, eye, ztemp1, ztemp2
      complex *16, allocatable :: mptemp1(:),mptemp2(:)      
      real *8 dc, ds, dc2, ds2, dcx, dsx, dcy, dsy
      integer ifder, ifders, l, j, i
      data eye /(0.0d0,1.0d0)/

      allocate(mptemp1(0:nterms+6),mptemp2(0:nterms+6),
     1     diffs(0:nterms+6),kvec(0:nterms+6))

      do i = 1,ntarg
         zdiff(1) = ztarg(1,i)-center(1)
	 zdiff(2) = ztarg(2,i)-center(2)
      	 call h2cart2polar(zdiff,r,theta)
         
c     get values of difference functions
         ifders = 0
         call diffslogbk_fast(r,beta,rscale,diffs,ifders,ders,kvec,
     1        nterms+2)       
         
         mptemp1(0)=diffs(0)
         mptemp2(0)=kvec(0)
         ztemp2=exp(eye*theta)
         ztemp1=ztemp2
	 do j=1,nterms+2
   	    mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1           diffs(j)*dimag(ztemp1))
	    mptemp2(j)=dcmplx(kvec(j)*dreal(ztemp1),
     1           kvec(j)*dimag(ztemp1))
            ztemp1 = ztemp1*ztemp2
	 enddo
         
c     evaluate
         do j = 1,nd
            pot(j,i) = pot(j,i) +
     1           dreal(mptemp1(0))*dreal(mbhmpole(j,0)) +
     1           dreal(mptemp2(0))*dreal(ympole(j,0))
         enddo
            
         do l = 1,nterms
            do j=1,nd
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp1(l))*dreal(mbhmpole(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhmpole(j,l))
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp2(l))*dreal(ympole(j,l)) +
     1              dimag(mptemp2(l))*dimag(ympole(j,l))
            enddo
         enddo
      enddo

      return 
      end

      subroutine mbh2dmpevalg_vec(nd,beta,rscale,center,
     1     mbhmpole,ympole,nterms,ztarg,ntarg,
     2     pot,grad)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor of the expansions
c     center(2)          : expansion center
c     mbhmpole(0:nterms) : coefficients for difference-type functions
c     ympole(nterms)     : coefficients for modified Bessel functions
c     nterms             : number of terms in expansions
c     ztarg(2,ntarg)     : coordinates of targets
c     ntarg              : number of targets
c
c     OUTPUT:
c      
c     pot(ntarg)      : value of potential at targets
c     grad(2,ntarg)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhmpole(nd,0:nterms), ympole(nd,0:nterms)
      real *8 rscale, center(2), ztarg(2,ntarg), beta
      real *8 pot(nd,ntarg)
      real *8 grad(nd,2,ntarg)
      integer nterms, ntarg, nd
c     local
      real *8 zdiff(2), r, theta
      real *8 ders(0:1), pih
      real *8, allocatable :: diffs(:), kvec(:)
      complex *16 z, eye, ztemp1, ztemp2
      complex *16, allocatable :: mptemp1(:),mptemp2(:)
      complex *16, allocatable :: ympolex(:,:), ympoley(:,:)
      complex *16, allocatable :: mbhmpolex(:,:), mbhmpoley(:,:)      
      real *8 dc, ds, dc2, ds2, dcx, dsx, dcy, dsy
      integer ifder, ifders, l, j, i
      data eye /(0.0d0,1.0d0)/
c     
      allocate(ympolex(nd,0:nterms+1),ympoley(nd,0:nterms+1),
     1     mbhmpolex(nd,0:nterms+1),mbhmpoley(nd,0:nterms+1))
      do i=0,nterms+1
         do j = 1,nd
            ympolex(j,i)=0
            ympoley(j,i)=0
            mbhmpolex(j,i)=0
            mbhmpoley(j,i)=0
         enddo
      enddo
      
      do j=1,nd
         ympolex(j,1) = -beta/rscale*ympole(j,0)
         ympoley(j,1) = -beta/rscale*ympole(j,0)*eye

         mbhmpolex(j,1) = -beta/rscale*mbhmpole(j,0)
         mbhmpoley(j,1) = -beta/rscale*mbhmpole(j,0)*eye
      enddo
      
      do i=1,nterms
         do j = 1,nd
            ympolex(j,i-1) = ympolex(j,i-1) 
     1           -beta/2.0d0*ympole(j,i)*rscale
            ympolex(j,i+1) = ympolex(j,i+1) 
     1           -beta/2.0d0*ympole(j,i)/rscale
            ympoley(j,i-1) = ympoley(j,i-1) 
     1           +beta/2.0d0*ympole(j,i)*rscale*eye
            ympoley(j,i+1) = ympoley(j,i+1)
     1           -beta/2.0d0*ympole(j,i)/rscale*eye

            ympolex(j,i-1) = ympolex(j,i-1) 
     1           -beta/2.0d0*mbhmpole(j,i)*rscale
            
            mbhmpolex(j,i+1) = mbhmpolex(j,i+1) 
     1           -beta/2.0d0*mbhmpole(j,i)/rscale
            ympoley(j,i-1) = ympoley(j,i-1) 
     1           +beta/2.0d0*mbhmpole(j,i)*rscale*eye
            mbhmpoley(j,i+1) = mbhmpoley(j,i+1) 
     1           -beta/2.0d0*mbhmpole(j,i)/rscale*eye
         enddo
      enddo

      allocate(mptemp1(0:nterms+6),mptemp2(0:nterms+6),
     1     diffs(0:nterms+6),kvec(0:nterms+6))

      do i = 1,ntarg
         zdiff(1) = ztarg(1,i)-center(1)
	 zdiff(2) = ztarg(2,i)-center(2)
      	 call h2cart2polar(zdiff,r,theta)
         
c     get values of difference functions
         ifders = 0
         call diffslogbk_fast(r,beta,rscale,diffs,ifders,ders,kvec,
     1        nterms+2)       
         
         mptemp1(0)=diffs(0)
         mptemp2(0)=kvec(0)
         ztemp2=exp(eye*theta)
         ztemp1=ztemp2
	 do j=1,nterms+2
   	    mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1           diffs(j)*dimag(ztemp1))
	    mptemp2(j)=dcmplx(kvec(j)*dreal(ztemp1),
     1           kvec(j)*dimag(ztemp1))
            ztemp1 = ztemp1*ztemp2
	 enddo
         
c     evaluate
         do j = 1,nd
            pot(j,i) = pot(j,i) +
     1           dreal(mptemp1(0))*dreal(mbhmpole(j,0)) +
     1           dreal(mptemp2(0))*dreal(ympole(j,0))
         enddo
            
         do l = 1,nterms
            do j=1,nd
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp1(l))*dreal(mbhmpole(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhmpole(j,l))
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp2(l))*dreal(ympole(j,l)) +
     1              dimag(mptemp2(l))*dimag(ympole(j,l))
            enddo
         enddo
         do j = 1,nd
            grad(j,1,i) = grad(j,1,i) + 
     1           dreal(mptemp1(0))*dreal(mbhmpolex(j,0)) +
     1           dreal(mptemp2(0))*dreal(ympolex(j,0))
            grad(j,2,i) = grad(j,2,i) + 
     1           dreal(mptemp1(0))*dreal(mbhmpoley(j,0)) +
     1           dreal(mptemp2(0))*dreal(ympoley(j,0))
         enddo
         do l = 1,nterms+1
            do j = 1,nd
               grad(j,1,i) = grad(j,1,i) +
     1              dreal(mptemp1(l))*dreal(mbhmpolex(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhmpolex(j,l))
               grad(j,1,i) = grad(j,1,i) +
     1              dreal(mptemp2(l))*dreal(ympolex(j,l)) +
     1              dimag(mptemp2(l))*dimag(ympolex(j,l))
               grad(j,2,i) = grad(j,2,i) +
     1              dreal(mptemp1(l))*dreal(mbhmpoley(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhmpoley(j,l))
               grad(j,2,i) = grad(j,2,i) +
     1              dreal(mptemp2(l))*dreal(ympoley(j,l)) +
     1              dimag(mptemp2(l))*dimag(ympoley(j,l))
            enddo
         enddo
      enddo

      return 
      end

      subroutine mbh2dmpevalh_vec(nd,beta,rscale,center,
     1     mbhmpole,ympole,nterms,ztarg,ntarg,
     2     pot,grad,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor of the expansions
c     center(2)          : expansion center
c     mbhmpole(0:nterms) : coefficients for difference-type functions
c     ympole(nterms)     : coefficients for modified Bessel functions
c     nterms             : number of terms in expansions
c     ztarg(2,ntarg)     : coordinates of targets
c     ntarg              : number of targets
c
c     OUTPUT:
c      
c     pot(ntarg)      : value of potential at targets
c     grad(2,ntarg)   : value of gradient at targets
c     hess(3,ntarg)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhmpole(nd,0:nterms), ympole(nd,0:nterms)
      real *8 rscale, center(2), ztarg(2,ntarg), beta
      real *8 pot(nd,ntarg)
      real *8 grad(nd,2,ntarg)
      real *8 hess(nd,3,ntarg)
      integer nterms, ntarg, nd
c     local
      real *8 zdiff(2), r, theta
      real *8 ders(0:1), pih
      real *8, allocatable :: diffs(:), kvec(:)
      complex *16 z, eye, ztemp1, ztemp2
      complex *16, allocatable :: mptemp1(:),mptemp2(:)
      complex *16, allocatable :: ympolex(:,:), ympoley(:,:)
      complex *16, allocatable :: mbhmpolex(:,:), mbhmpoley(:,:)
      complex *16, allocatable :: ympolexx(:,:), ympolexy(:,:),
     1     ympoleyy(:,:)
      complex *16, allocatable :: mbhmpolexx(:,:),mbhmpolexy(:,:),
     1     mbhmpoleyy(:,:)      
      real *8 dc, ds, dc2, ds2, dcx, dsx, dcy, dsy
      integer ifder, ifders, l, j, i
      data eye /(0.0d0,1.0d0)/
c     
      allocate(ympolex(nd,0:nterms+1),ympoley(nd,0:nterms+1),
     1     mbhmpolex(nd,0:nterms+1),mbhmpoley(nd,0:nterms+1))
      do i=0,nterms+1
         do j = 1,nd
            ympolex(j,i)=0
            ympoley(j,i)=0
            mbhmpolex(j,i)=0
            mbhmpoley(j,i)=0
         enddo
      enddo
      
      do j=1,nd
         ympolex(j,1) = -beta/rscale*ympole(j,0)
         ympoley(j,1) = -beta/rscale*ympole(j,0)*eye

         mbhmpolex(j,1) = -beta/rscale*mbhmpole(j,0)
         mbhmpoley(j,1) = -beta/rscale*mbhmpole(j,0)*eye
      enddo
      
      do i=1,nterms
         do j = 1,nd
            ympolex(j,i-1) = ympolex(j,i-1) 
     1           -beta/2.0d0*ympole(j,i)*rscale
            ympolex(j,i+1) = ympolex(j,i+1) 
     1           -beta/2.0d0*ympole(j,i)/rscale
            ympoley(j,i-1) = ympoley(j,i-1) 
     1           +beta/2.0d0*ympole(j,i)*rscale*eye
            ympoley(j,i+1) = ympoley(j,i+1)
     1           -beta/2.0d0*ympole(j,i)/rscale*eye

            ympolex(j,i-1) = ympolex(j,i-1) 
     1           -beta/2.0d0*mbhmpole(j,i)*rscale
            
            mbhmpolex(j,i+1) = mbhmpolex(j,i+1) 
     1           -beta/2.0d0*mbhmpole(j,i)/rscale
            ympoley(j,i-1) = ympoley(j,i-1) 
     1           +beta/2.0d0*mbhmpole(j,i)*rscale*eye
            mbhmpoley(j,i+1) = mbhmpoley(j,i+1) 
     1           -beta/2.0d0*mbhmpole(j,i)/rscale*eye
         enddo
      enddo
c
      allocate(ympolexx(nd,0:nterms+2),ympolexy(nd,0:nterms+2),
     1     ympoleyy(nd,0:nterms+2),
     1     mbhmpolexx(nd,0:nterms+2),mbhmpolexy(nd,0:nterms+2),
     1     mbhmpoleyy(nd,0:nterms+2))

      do i=0,nterms+2
         do j=1,nd
            ympolexx(j,i)=0
            ympolexy(j,i)=0
            ympoleyy(j,i)=0
            mbhmpolexx(j,i)=0
            mbhmpolexy(j,i)=0
            mbhmpoleyy(j,i)=0
         enddo
      enddo

      do j = 1,nd
         ympolexx(j,1) = -beta/1.0d0/rscale*dreal(ympolex(j,0))
         ympolexy(j,1) = -beta/1.0d0/rscale*dreal(ympolex(j,0))*eye
         ympoleyy(j,1) = -beta/1.0d0/rscale*dreal(ympoley(j,0))*eye

         mbhmpolexx(j,1) = -beta/1.0d0/rscale*dreal(mbhmpolex(j,0))
         mbhmpolexy(j,1) = 
     1        -beta/1.0d0/rscale*dreal(mbhmpolex(j,0))*eye
         mbhmpoleyy(j,1) = 
     1        -beta/1.0d0/rscale*dreal(mbhmpoley(j,0))*eye
      enddo
      
      do i=1,nterms+1
         do j=1,nd
            ympolexx(j,i-1) = ympolexx(j,i-1)
     1  	 -beta/2.0d0*ympolex(j,i)*rscale
            ympolexx(j,i+1) = ympolexx(j,i+1)
     1           -beta/2.0d0*ympolex(j,i)/rscale
            ympolexy(j,i-1) = ympolexy(j,i-1) 
     1           +beta/2.0d0*ympolex(j,i)*rscale*eye
            ympolexy(j,i+1) = ympolexy(j,i+1) 
     1           -beta/2.0d0*ympolex(j,i)/rscale*eye
            ympoleyy(j,i-1) = ympoleyy(j,i-1) 
     1           +beta/2.0d0*ympoley(j,i)*rscale*eye
            ympoleyy(j,i+1) = ympoleyy(j,i+1) 
     1           -beta/2.0d0*ympoley(j,i)/rscale*eye

            ympolexx(j,i-1) = ympolexx(j,i-1) 
     1           -beta/2.0d0*mbhmpolex(j,i)*rscale
            mbhmpolexx(j,i+1) = mbhmpolexx(j,i+1) 
     1           -beta/2.0d0*mbhmpolex(j,i)/rscale
            ympolexy(j,i-1) = ympolexy(j,i-1) 
     1           +beta/2.0d0*mbhmpolex(j,i)*rscale*eye
            mbhmpolexy(j,i+1) = mbhmpolexy(j,i+1) 
     1           -beta/2.0d0*mbhmpolex(j,i)/rscale*eye
            ympoleyy(j,i-1) = ympoleyy(j,i-1) 
     1           +beta/2.0d0*mbhmpoley(j,i)*rscale*eye
            mbhmpoleyy(j,i+1) = mbhmpoleyy(j,i+1) 
     1           -beta/2.0d0*mbhmpoley(j,i)/rscale*eye
         enddo
      enddo

      allocate(mptemp1(0:nterms+6),mptemp2(0:nterms+6),
     1     diffs(0:nterms+6),kvec(0:nterms+6))

      do i = 1,ntarg
         zdiff(1) = ztarg(1,i)-center(1)
	 zdiff(2) = ztarg(2,i)-center(2)
      	 call h2cart2polar(zdiff,r,theta)
         
c     get values of difference functions
         ifders = 0
         call diffslogbk_fast(r,beta,rscale,diffs,ifders,ders,kvec,
     1        nterms+2)       
         
         mptemp1(0)=diffs(0)
         mptemp2(0)=kvec(0)
         ztemp2=exp(eye*theta)
         ztemp1=ztemp2
	 do j=1,nterms+2
   	    mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1           diffs(j)*dimag(ztemp1))
	    mptemp2(j)=dcmplx(kvec(j)*dreal(ztemp1),
     1           kvec(j)*dimag(ztemp1))
            ztemp1 = ztemp1*ztemp2
	 enddo
         
c     evaluate
         do j = 1,nd
            pot(j,i) = pot(j,i) +
     1           dreal(mptemp1(0))*dreal(mbhmpole(j,0)) +
     1           dreal(mptemp2(0))*dreal(ympole(j,0))
         enddo
            
         do l = 1,nterms
            do j=1,nd
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp1(l))*dreal(mbhmpole(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhmpole(j,l))
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp2(l))*dreal(ympole(j,l)) +
     1              dimag(mptemp2(l))*dimag(ympole(j,l))
            enddo
         enddo
         do j = 1,nd
            grad(j,1,i) = grad(j,1,i) + 
     1           dreal(mptemp1(0))*dreal(mbhmpolex(j,0)) +
     1           dreal(mptemp2(0))*dreal(ympolex(j,0))
            grad(j,2,i) = grad(j,2,i) + 
     1           dreal(mptemp1(0))*dreal(mbhmpoley(j,0)) +
     1           dreal(mptemp2(0))*dreal(ympoley(j,0))
         enddo
         do l = 1,nterms+1
            do j = 1,nd
               grad(j,1,i) = grad(j,1,i) +
     1              dreal(mptemp1(l))*dreal(mbhmpolex(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhmpolex(j,l))
               grad(j,1,i) = grad(j,1,i) +
     1              dreal(mptemp2(l))*dreal(ympolex(j,l)) +
     1              dimag(mptemp2(l))*dimag(ympolex(j,l))
               grad(j,2,i) = grad(j,2,i) +
     1              dreal(mptemp1(l))*dreal(mbhmpoley(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhmpoley(j,l))
               grad(j,2,i) = grad(j,2,i) +
     1              dreal(mptemp2(l))*dreal(ympoley(j,l)) +
     1              dimag(mptemp2(l))*dimag(ympoley(j,l))
            enddo
         enddo
         do j=1,nd
            hess(j,1,i) = hess(j,1,i) +
     1           dreal(mptemp1(0))*dreal(mbhmpolexx(j,0)) +
     1           dreal(mptemp2(0))*dreal(ympolexx(j,0))
            hess(j,2,i) = hess(j,2,i) +
     1           dreal(mptemp1(0))*dreal(mbhmpolexy(j,0)) +
     1           dreal(mptemp2(0))*dreal(ympolexy(j,0))
            hess(j,3,i) = hess(j,3,i) +
     1           dreal(mptemp1(0))*dreal(mbhmpoleyy(j,0)) +
     1           dreal(mptemp2(0))*dreal(ympoleyy(j,0))
         enddo
         do l = 1,nterms+1
            do j = 1,nd
               hess(j,1,i) = hess(j,1,i) +
     1              dreal(mptemp1(l))*dreal(mbhmpolexx(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhmpolexx(j,l))
               hess(j,1,i) = hess(j,1,i) + 
     1              dreal(mptemp2(l))*dreal(ympolexx(j,l)) +
     1              dimag(mptemp2(l))*dimag(ympolexx(j,l))
               hess(j,2,i) = hess(j,2,i) + 
     1              dreal(mptemp1(l))*dreal(mbhmpolexy(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhmpolexy(j,l))
               hess(j,2,i) = hess(j,2,i) + 
     1              dreal(mptemp2(l))*dreal(ympolexy(j,l)) +
     1              dimag(mptemp2(l))*dimag(ympolexy(j,l))
               hess(j,3,i) = hess(j,3,i) + 
     1              dreal(mptemp1(l))*dreal(mbhmpoleyy(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhmpoleyy(j,l))
               hess(j,3,i) = hess(j,3,i) + 
     1              dreal(mptemp2(l))*dreal(ympoleyy(j,l)) +
     1              dimag(mptemp2(l))*dimag(ympoleyy(j,l))
            enddo
         enddo
      enddo

      return 
      end
      

      
      subroutine mbh2dtaevalp_vec(nd,beta,rscale,center,mbhloc,lloc,
     1   nterms,ztarg,ntarg,
     2   pot)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor
c     center(2)          : expansion center
c     mbhloc(0:nterms)   : coefficients for difference-type functions
c     lloc(0:nterms)     : coefficients for Taylor series
c     nterms             : number of terms in expansions
c     ztarg(2,ntarg)     : coordinates of targets
c     ntarg              : number of targets
c
c     OUTPUT:
c      
c     pot(ntarg)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhloc(nd,0:nterms), lloc(nd,0:nterms)
      real *8 rscale, center(2), ztarg(2,ntarg), beta
      real *8 pot(nd,ntarg), grad(nd,2,ntarg), hess(nd,3,ntarg)
      integer nterms, ntarg, nd
c     local
      real *8 zdiff(2), r, theta, ivec(0:200)
      real *8 pih, ders(0:1)
      real *8, allocatable :: diffs(:), pow(:), dpow(:)
      complex *16 z, eye, ztemp1, ztemp2
      complex *16, allocatable, dimension (:,:) :: mbhlocx,mbhlocy,
     1     mbhlocxx,mbhlocxy,mbhlocyy,llocx,llocy,llocxx,llocxy,llocyy
      complex *16, allocatable :: mptemp1(:), mptemp2(:)
      real *8 dtemp
      integer ifder, ifders, l, j, i
      data eye /(0.0d0,1.0d0)/

      allocate(pow(0:nterms+6),dpow(0:nterms+6))
      allocate(mptemp1(0:nterms+6),mptemp2(0:nterms+6))
      allocate(diffs(0:nterms+6))

c     evaluate

      do i = 1,ntarg

         zdiff(1) = ztarg(1,i)-center(1)
         zdiff(2) = ztarg(2,i)-center(2)
         call h2cart2polar(zdiff,r,theta)
         
c     get values of (beta*r)^k
         call mbh2d_rk(pow,dpow,r,beta,rscale,nterms+2)

c     get values of difference functions
         ifders = 0
         call diffszkik_fast(r,beta,rscale,diffs,ifders,ders,ivec,
     1        nterms+2)
         
         mptemp1(0)=diffs(0)
         mptemp2(0)=pow(0)
         ztemp2=exp(eye*theta)
         ztemp1=ztemp2
         do j=1,nterms+2
            mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1           diffs(j)*dimag(ztemp1))
            mptemp2(j)=dcmplx(pow(j)*dreal(ztemp1),pow(j)*dimag(ztemp1))
            ztemp1 = ztemp1*ztemp2
         enddo
         do j = 1,nd
            pot(j,i) = pot(j,i) + dreal(mptemp1(0))*dreal(mbhloc(j,0))
     1           + dreal(mptemp2(0))*dreal(lloc(j,0))
         enddo
         do l = 1,nterms
            do j = 1,nd
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp1(l))*dreal(mbhloc(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhloc(j,l))
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp2(l))*dreal(lloc(j,l)) +
     1              dimag(mptemp2(l))*dimag(lloc(j,l))
            enddo
         enddo
      enddo
      
      return 
      end

      
      subroutine mbh2dtaevalg_vec(nd,beta,rscale,center,mbhloc,lloc,
     1   nterms,ztarg,ntarg,
     2   pot,grad)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor
c     center(2)          : expansion center
c     mbhloc(0:nterms)   : coefficients for difference-type functions
c     lloc(0:nterms)     : coefficients for Taylor series
c     nterms             : number of terms in expansions
c     ztarg(2,ntarg)     : coordinates of targets
c     ntarg              : number of targets
c
c     OUTPUT:
c      
c     pot(ntarg)      : value of potential at targets
c     grad(2,ntarg)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhloc(nd,0:nterms), lloc(nd,0:nterms)
      real *8 rscale, center(2), ztarg(2,ntarg), beta
      real *8 pot(nd,ntarg), grad(nd,2,ntarg), hess(nd,3,ntarg)
      integer nterms, ntarg, nd
c     local
      real *8 zdiff(2), r, theta, ivec(0:200)
      real *8 pih, ders(0:1)
      real *8, allocatable :: diffs(:), pow(:), dpow(:)
      complex *16 z, eye, ztemp1, ztemp2
      complex *16, allocatable, dimension (:,:) :: mbhlocx,mbhlocy,
     1     mbhlocxx,mbhlocxy,mbhlocyy,llocx,llocy,llocxx,llocxy,llocyy
      complex *16, allocatable :: mptemp1(:), mptemp2(:)
      real *8 dtemp
      integer ifder, ifders, l, j, i
      data eye /(0.0d0,1.0d0)/
      allocate(mbhlocx(nd,0:nterms+1),mbhlocy(nd,0:nterms+1),
     1     llocx(nd,0:nterms+1),llocy(nd,0:nterms+1))
      do i=0,nterms+1
         do j = 1,nd
            mbhlocx(j,i)=0
            mbhlocy(j,i)=0
            llocx(j,i)=0
            llocy(j,i)=0
         enddo
      enddo
      
      dtemp = 1.0d0/2.0d0

      do j = 1,nd
         mbhlocx(j,1) = dreal(mbhloc(j,0))*beta*rscale
         mbhlocy(j,1) = dreal(mbhloc(j,0))*beta*rscale*eye
         llocx(j,1) = dreal(mbhloc(j,0))*beta*dtemp*rscale
         llocy(j,1) = dreal(mbhloc(j,0))*beta*dtemp*rscale*eye
      enddo
      do i=1,nterms
         dtemp = dtemp/(2*(i+1))
         do j = 1,nd
            mbhlocx(j,i-1)=mbhlocx(j,i-1)+beta/2/rscale*mbhloc(j,i)
            mbhlocx(j,i+1)=mbhlocx(j,i+1)+beta/2*rscale*mbhloc(j,i)
            mbhlocy(j,i-1)=mbhlocy(j,i-1) -
     1           beta/2*(eye)/rscale*mbhloc(j,i)
            mbhlocy(j,i+1)=mbhlocy(j,i+1) +
     1           beta/2*(eye)*rscale*mbhloc(j,i)
            
            llocx(j,i-1)=llocx(j,i-1) + beta*i/rscale*lloc(j,i)
            llocy(j,i-1)=llocy(j,i-1) - beta*i*(eye)/rscale*lloc(j,i)
            
            llocx(j,i+1)=llocx(j,i+1) + beta/2*rscale*mbhloc(j,i)*dtemp
            llocy(j,i+1)=llocy(j,i+1) +
     1           beta/2*(eye)*rscale*mbhloc(j,i)*dtemp
         enddo
      enddo

      allocate(pow(0:nterms+6),dpow(0:nterms+6))
      allocate(mptemp1(0:nterms+6),mptemp2(0:nterms+6))
      allocate(diffs(0:nterms+6))

c     evaluate

      do i = 1,ntarg

         zdiff(1) = ztarg(1,i)-center(1)
         zdiff(2) = ztarg(2,i)-center(2)
         call h2cart2polar(zdiff,r,theta)
         
c     get values of (beta*r)^k
         call mbh2d_rk(pow,dpow,r,beta,rscale,nterms+2)

c     get values of difference functions
         ifders = 0
         call diffszkik_fast(r,beta,rscale,diffs,ifders,ders,ivec,
     1        nterms+2)
         
         mptemp1(0)=diffs(0)
         mptemp2(0)=pow(0)
         ztemp2=exp(eye*theta)
         ztemp1=ztemp2
         do j=1,nterms+2
            mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1           diffs(j)*dimag(ztemp1))
            mptemp2(j)=dcmplx(pow(j)*dreal(ztemp1),pow(j)*dimag(ztemp1))
            ztemp1 = ztemp1*ztemp2
         enddo
         do j = 1,nd
            pot(j,i) = pot(j,i) + dreal(mptemp1(0))*dreal(mbhloc(j,0))
     1           + dreal(mptemp2(0))*dreal(lloc(j,0))
         enddo
         do l = 1,nterms
            do j = 1,nd
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp1(l))*dreal(mbhloc(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhloc(j,l))
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp2(l))*dreal(lloc(j,l)) +
     1              dimag(mptemp2(l))*dimag(lloc(j,l))
            enddo
         enddo
         do j=1,nd
            grad(j,1,i) = grad(j,1,i) +
     1           dreal(mptemp1(0))*dreal(mbhlocx(j,0)) +
     1           dreal(mptemp2(0))*dreal(llocx(j,0))
            grad(j,2,i) = grad(j,2,i) +
     1           dreal(mptemp1(0))*dreal(mbhlocy(j,0)) +
     1           dreal(mptemp2(0))*dreal(llocy(j,0))
         enddo

         do l = 1,nterms+1
            do j = 1,nd
               grad(j,1,i) = grad(j,1,i) +
     1              dreal(mptemp1(l))*dreal(mbhlocx(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhlocx(j,l))
               grad(j,1,i) = grad(j,1,i) +
     1              dreal(mptemp2(l))*dreal(llocx(j,l)) +
     1              dimag(mptemp2(l))*dimag(llocx(j,l))
               grad(j,2,i) = grad(j,2,i) +
     1              dreal(mptemp1(l))*dreal(mbhlocy(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhlocy(j,l))
               grad(j,2,i) = grad(j,2,i) +
     1              dreal(mptemp2(l))*dreal(llocy(j,l)) +
     1              dimag(mptemp2(l))*dimag(llocy(j,l))
            enddo
         enddo
      enddo
      
      return 
      end

      
      subroutine mbh2dtaevalh_vec(nd,beta,rscale,center,mbhloc,lloc,
     1   nterms,ztarg,ntarg,
     2   pot,grad,hess)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     beta               : Yukawa parameter
c     rscale             : the scaling factor
c     center(2)          : expansion center
c     mbhloc(0:nterms)   : coefficients for difference-type functions
c     lloc(0:nterms)     : coefficients for Taylor series
c     nterms             : number of terms in expansions
c     ztarg(2,ntarg)     : coordinates of targets
c     ntarg              : number of targets
c
c     OUTPUT:
c      
c     pot(ntarg)      : value of potential at targets
c     grad(2,ntarg)   : value of gradient at targets
c     hess(3,ntarg)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhloc(nd,0:nterms), lloc(nd,0:nterms)
      real *8 rscale, center(2), ztarg(2,ntarg), beta
      real *8 pot(nd,ntarg), grad(nd,2,ntarg), hess(nd,3,ntarg)
      integer nterms, ntarg, nd
c     local
      real *8 zdiff(2), r, theta, ivec(0:200)
      real *8 pih, ders(0:1)
      real *8, allocatable :: diffs(:), pow(:), dpow(:)
      complex *16 z, eye, ztemp1, ztemp2
      complex *16, allocatable, dimension (:,:) :: mbhlocx,mbhlocy,
     1     mbhlocxx,mbhlocxy,mbhlocyy,llocx,llocy,llocxx,llocxy,llocyy
      complex *16, allocatable :: mptemp1(:), mptemp2(:)
      real *8 dtemp
      integer ifder, ifders, l, j, i
      data eye /(0.0d0,1.0d0)/
      allocate(mbhlocx(nd,0:nterms+1),mbhlocy(nd,0:nterms+1),
     1     llocx(nd,0:nterms+1),llocy(nd,0:nterms+1))
      do i=0,nterms+1
         do j = 1,nd
            mbhlocx(j,i)=0
            mbhlocy(j,i)=0
            llocx(j,i)=0
            llocy(j,i)=0
         enddo
      enddo
      
      dtemp = 1.0d0/2.0d0

      do j = 1,nd
         mbhlocx(j,1) = dreal(mbhloc(j,0))*beta*rscale
         mbhlocy(j,1) = dreal(mbhloc(j,0))*beta*rscale*eye
         llocx(j,1) = dreal(mbhloc(j,0))*beta*dtemp*rscale
         llocy(j,1) = dreal(mbhloc(j,0))*beta*dtemp*rscale*eye
      enddo
      do i=1,nterms
         dtemp = dtemp/(2*(i+1))
         do j = 1,nd
            mbhlocx(j,i-1)=mbhlocx(j,i-1)+beta/2/rscale*mbhloc(j,i)
            mbhlocx(j,i+1)=mbhlocx(j,i+1)+beta/2*rscale*mbhloc(j,i)
            mbhlocy(j,i-1)=mbhlocy(j,i-1) -
     1           beta/2*(eye)/rscale*mbhloc(j,i)
            mbhlocy(j,i+1)=mbhlocy(j,i+1) +
     1           beta/2*(eye)*rscale*mbhloc(j,i)
            
            llocx(j,i-1)=llocx(j,i-1) + beta*i/rscale*lloc(j,i)
            llocy(j,i-1)=llocy(j,i-1) - beta*i*(eye)/rscale*lloc(j,i)
            
            llocx(j,i+1)=llocx(j,i+1) + beta/2*rscale*mbhloc(j,i)*dtemp
            llocy(j,i+1)=llocy(j,i+1) +
     1           beta/2*(eye)*rscale*mbhloc(j,i)*dtemp
         enddo
      enddo
c
      allocate(mbhlocxx(nd,0:nterms+2),mbhlocxy(nd,0:nterms+2),
     1     mbhlocyy(nd,0:nterms+2),
     1     llocxx(nd,0:nterms+2),llocxy(nd,0:nterms+2),
     1     llocyy(nd,0:nterms+2))
      
      do i=0,nterms+2
         do j = 1,nd
            mbhlocxx(j,i)=0
            mbhlocxy(j,i)=0
            mbhlocyy(j,i)=0
            llocxx(j,i)=0
            llocxy(j,i)=0
            llocyy(j,i)=0
         enddo
      enddo

      dtemp = 1.0d0/2.0d0

      do j = 1,nd
         mbhlocxx(j,1) = dreal(mbhlocx(j,0))*beta*rscale
         mbhlocxy(j,1) = dreal(mbhlocx(j,0))*beta*rscale*eye
         mbhlocyy(j,1) = dreal(mbhlocy(j,0))*beta*rscale*eye
         llocxx(j,1) = dreal(mbhlocx(j,0))*beta*dtemp*rscale
         llocxy(j,1) = dreal(mbhlocx(j,0))*beta*dtemp*rscale*eye
         llocyy(j,1) = dreal(mbhlocy(j,0))*beta*dtemp*rscale*eye
      enddo
      
      do i=1,nterms+1
         dtemp = dtemp/(2*(i+1))
         do j = 1,nd
            mbhlocxx(j,i-1)=mbhlocxx(j,i-1)+beta/2/rscale*mbhlocx(j,i)
            mbhlocxx(j,i+1)=mbhlocxx(j,i+1)+beta/2*rscale*mbhlocx(j,i)
            mbhlocxy(j,i-1)=mbhlocxy(j,i-1)-
     1           beta/2*(eye)/rscale*mbhlocx(j,i)
            mbhlocxy(j,i+1)=mbhlocxy(j,i+1)+
     1           beta/2*(eye)*rscale*mbhlocx(j,i)
            mbhlocyy(j,i-1)=mbhlocyy(j,i-1)-
     1           beta/2*(eye)/rscale*mbhlocy(j,i)
            mbhlocyy(j,i+1)=mbhlocyy(j,i+1)+
     1           beta/2*(eye)*rscale*mbhlocy(j,i)
            
            llocxx(j,i-1)=llocxx(j,i-1)+beta*i/rscale*llocx(j,i)
            llocxy(j,i-1)=llocxy(j,i-1)-beta*i*(eye)/rscale*llocx(j,i)
            llocyy(j,i-1)=llocyy(j,i-1)-beta*i*(eye)/rscale*llocy(j,i)

            llocxx(j,i+1)=llocxx(j,i+1)+beta/2*rscale*mbhlocx(j,i)*dtemp
            llocxy(j,i+1)=llocxy(j,i+1)+
     1           beta/2*(eye)*rscale*mbhlocx(j,i)*dtemp
            llocyy(j,i+1)=llocyy(j,i+1)+
     1           beta/2*(eye)*rscale*mbhlocy(j,i)*dtemp
         enddo
      enddo

      allocate(pow(0:nterms+6),dpow(0:nterms+6))
      allocate(mptemp1(0:nterms+6),mptemp2(0:nterms+6))
      allocate(diffs(0:nterms+6))

c     evaluate

      do i = 1,ntarg

         zdiff(1) = ztarg(1,i)-center(1)
         zdiff(2) = ztarg(2,i)-center(2)
         call h2cart2polar(zdiff,r,theta)
         
c     get values of (beta*r)^k
         call mbh2d_rk(pow,dpow,r,beta,rscale,nterms+2)

c     get values of difference functions
         ifders = 0
         call diffszkik_fast(r,beta,rscale,diffs,ifders,ders,ivec,
     1        nterms+2)
         
         mptemp1(0)=diffs(0)
         mptemp2(0)=pow(0)
         ztemp2=exp(eye*theta)
         ztemp1=ztemp2
         do j=1,nterms+2
            mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1           diffs(j)*dimag(ztemp1))
            mptemp2(j)=dcmplx(pow(j)*dreal(ztemp1),pow(j)*dimag(ztemp1))
            ztemp1 = ztemp1*ztemp2
         enddo
         do j = 1,nd
            pot(j,i) = pot(j,i) + dreal(mptemp1(0))*dreal(mbhloc(j,0))
     1           + dreal(mptemp2(0))*dreal(lloc(j,0))
         enddo
         do l = 1,nterms
            do j = 1,nd
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp1(l))*dreal(mbhloc(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhloc(j,l))
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp2(l))*dreal(lloc(j,l)) +
     1              dimag(mptemp2(l))*dimag(lloc(j,l))
            enddo
         enddo
         do j=1,nd
            grad(j,1,i) = grad(j,1,i) +
     1           dreal(mptemp1(0))*dreal(mbhlocx(j,0)) +
     1           dreal(mptemp2(0))*dreal(llocx(j,0))
            grad(j,2,i) = grad(j,2,i) +
     1           dreal(mptemp1(0))*dreal(mbhlocy(j,0)) +
     1           dreal(mptemp2(0))*dreal(llocy(j,0))
         enddo

         do l = 1,nterms+1
            do j = 1,nd
               grad(j,1,i) = grad(j,1,i) +
     1              dreal(mptemp1(l))*dreal(mbhlocx(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhlocx(j,l))
               grad(j,1,i) = grad(j,1,i) +
     1              dreal(mptemp2(l))*dreal(llocx(j,l)) +
     1              dimag(mptemp2(l))*dimag(llocx(j,l))
               grad(j,2,i) = grad(j,2,i) +
     1              dreal(mptemp1(l))*dreal(mbhlocy(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhlocy(j,l))
               grad(j,2,i) = grad(j,2,i) +
     1              dreal(mptemp2(l))*dreal(llocy(j,l)) +
     1              dimag(mptemp2(l))*dimag(llocy(j,l))
            enddo
         enddo
         do j = 1,nd
            hess(j,1,i) = hess(j,1,i) +
     1           dreal(mptemp1(0))*dreal(mbhlocxx(j,0))
     1           + dreal(mptemp2(0))*dreal(llocxx(j,0))
            hess(j,2,i) = hess(j,2,i) +
     1           dreal(mptemp1(0))*dreal(mbhlocxy(j,0))
     1           + dreal(mptemp2(0))*dreal(llocxy(j,0))
            hess(j,3,i) = hess(j,3,i) +
     1           dreal(mptemp1(0))*dreal(mbhlocyy(j,0))
     1           + dreal(mptemp2(0))*dreal(llocyy(j,0))
         enddo
         do l = 1,nterms+1
            do j = 1,nd
               hess(j,1,i) = hess(j,1,i) +
     1              dreal(mptemp1(l))*dreal(mbhlocxx(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhlocxx(j,l))
               hess(j,1,i) = hess(j,1,i) +
     1              dreal(mptemp2(l))*dreal(llocxx(j,l)) +
     1              dimag(mptemp2(l))*dimag(llocxx(j,l))
               hess(j,2,i) = hess(j,2,i) +
     1              dreal(mptemp1(l))*dreal(mbhlocxy(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhlocxy(j,l))
               hess(j,2,i) = hess(j,2,i) +
     1              dreal(mptemp2(l))*dreal(llocxy(j,l)) +
     1              dimag(mptemp2(l))*dimag(llocxy(j,l))
               hess(j,3,i) = hess(j,3,i) +
     1              dreal(mptemp1(l))*dreal(mbhlocyy(j,l)) +
     1              dimag(mptemp1(l))*dimag(mbhlocyy(j,l))
               hess(j,3,i) = hess(j,3,i) +
     1              dreal(mptemp2(l))*dreal(llocyy(j,l)) +
     1              dimag(mptemp2(l))*dimag(llocyy(j,l))
            enddo
         enddo
      enddo
      
      return 
      end


      
      subroutine mbh2dformmpmp_vec(nd,beta,rscale,source,ns,
     1     mbhmpolesrc,ympolesrc,ntermsrc,center,nterms,mbhmpole,ympole)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     form a multipole expansion based on multipolar sources
c
c     nd - integer, number of expansions per source
c     beta - real *8, modified biharmonic parameter
c     rscale - real *8, scaling factor for output expansions
c     source - real *8 (2,ns) source locations
c     ns - integer, number of sources
c     mbhmpolesrc - complex *16 (nd,0:ntermsrc,ns) unscaled multipolar
c                   charge of difference kernel type at each source
c     ympolesrc - complex *16 (nd,0:ntermsrc,ns) unscaled multipolar
c                   charge of Yukawa type at each source
c     center - real *8 (2), center about which to form the expansion
c     nterms - integer, number of terms in output expansion
c
c     output:
c      
c     mbhmpole - complex *16 (nd,0:nterms), difference kernel type
c                expansion for source charges
c     ympole - complex *16 (nd,0:nterms), Yukawa type
c                expansion for source charges
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ns, ntermsrc, nterms, nd
      real *8 beta, rscale, source(2,*)
      real *8 center(2)
      complex *16 mbhmpolesrc(nd,0:ntermsrc,ns),
     1     ympolesrc(nd,0:ntermsrc,ns)
      complex *16 mbhmpole(nd,0:nterms), ympole(nd,0:nterms)
c     local variables
      real *8 rscale1,rsc
      integer i,j,l
c     
      rscale1=rscale

      do i = 1,ns

         rsc=1
         do j = 0,ntermsrc
            do l = 1,nd
               mbhmpolesrc(l,j,i)=mbhmpolesrc(l,j,i)/rsc
               ympolesrc(l,j,i)=ympolesrc(l,j,i)/rsc               
            enddo
            rsc=rsc*rscale1
         enddo
         call mbh2dmpmp_vec(nd,beta,rscale1,source(1,i),
     1        mbhmpolesrc(1,0,i),ympolesrc(1,0,i),ntermsrc,
     2        rscale,center,mbhmpole,ympole,nterms)
         
         rsc=1
         do j = 0,ntermsrc
            do l = 1,nd
               mbhmpolesrc(l,j,i)=mbhmpolesrc(l,j,i)/rsc
               ympolesrc(l,j,i)=ympolesrc(l,j,i)/rsc               
            enddo
            rsc=rsc/rscale1
         enddo
         
      enddo
      

      return
      end

      subroutine mbh2dformtamp_vec(nd,beta,rscale,source,ns,
     1     mbhmpolesrc,ympolesrc,ntermsrc,center,nterms,
     2     mbhloc,lloc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     form a local expansion based on multipolar sources
c
c     nd - integer, number of expansions per source
c     beta - real *8, modified biharmonic parameter
c     rscale - real *8, scaling factor for output expansions
c     source - real *8 (2,ns) source locations
c     ns - integer, number of sources
c     mbhmpolesrc - complex *16 (nd,0:ntermsrc,ns) unscaled multipolar
c                   charge of difference kernel type at each source
c     ympolesrc - complex *16 (nd,0:ntermsrc,ns) unscaled multipolar
c                   charge of Yukawa type at each source
c     center - real *8 (2), center about which to form the expansion
c     nterms - integer, number of terms in output expansion
c
c     output:
c      
c     mbhloc - complex *16 (nd,0:nterms), difference kernel type
c                expansion for source charges
c     lloc - complex *16 (nd,0:nterms), power series type
c                expansion for source charges
c     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      integer ns, ntermsrc, nterms, nd
      real *8 beta, rscale, source(2,*)
      real *8 center(2)
      complex *16 mbhmpolesrc(nd,0:ntermsrc,ns),
     1     ympolesrc(nd,0:ntermsrc,ns)
      complex *16 mbhloc(nd,0:nterms), lloc(nd,0:nterms)
c     local variables
      integer i, j, l
      real *8 rscale1,rsc

      rscale1=rscale
      
      do i = 1,ns

         rsc=1
         do j = 0,ntermsrc
            do l = 1,nd
               mbhmpolesrc(l,j,i)=mbhmpolesrc(l,j,i)/rsc
               ympolesrc(l,j,i)=ympolesrc(l,j,i)/rsc               
            enddo
            rsc=rsc*rscale1
         enddo
         
         call mbh2dmploc_vec(nd,beta,rscale1,source(1,i),
     1        mbhmpolesrc(1,0,i),ympolesrc(1,0,i),ntermsrc,
     2        rscale,center,mbhloc,lloc,nterms)

         rsc=1
         do j = 0,ntermsrc
            do l = 1,nd
               mbhmpolesrc(l,j,i)=mbhmpolesrc(l,j,i)/rsc
               ympolesrc(l,j,i)=ympolesrc(l,j,i)/rsc               
            enddo
            rsc=rsc/rscale1
         enddo
         
      enddo

      return
      end


      subroutine mbh2dmpmp_vec(nd,beta,rscale1,center1,
     1     mbhmpole1,ympole1,nterms1,rscale2,center2,mbhmpole2,ympole2,
     2     nterms2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     translate multipole expansions
c
c     INPUT:
c
c     nd                : number of densities
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     mbhmpole1(nd,0:nterms1): original coeffs, difference-type functions
c     ympole1(nd,0:nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     nterms2           : number of terms in original expansion
c     beta              : the modified biharmonic parameter
c
c     OUTPUT:
c
c     mbhmpole2(nd,0:nterms2): new coeffs, difference-type functions
c     ympole2(nd,0:nterms2): new coeffs, modified Bessel functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 rscale1, center1(2), beta
      real *8 rscale2, center2(2)
      complex *16 mbhmpole1(nd,0:nterms1), ympole1(nd,0:nterms1)
      complex *16 mbhmpole2(nd,0:nterms2), ympole2(nd,0:nterms2)
      integer nterms1, nterms2, nd
c     local
      integer i, j, ifder, nterms, m, l
      real *8 zdiff(2), r, theta, pi, ders(0:1)
      real *8, allocatable :: diffs(:), ival(:), dfac(:), dfac2(:)
      real *8, allocatable :: pow(:), dpow(:)
      real *8 twojm1
      real *8 rsj, rsi5, rsi7, rsi, fs2, rtemp
      complex *16, allocatable :: jtemp(:), difftemp(:), powtemp(:)
      complex *16 zk, z, ima, ztemp1, zmul, ztemp
      data ima /(0.0d0,1.0d0)/

      pi = 4.0d0*datan(1.0d0)
      zk = ima*beta

      nterms = nterms1+nterms2

      allocate(jtemp(0:nterms+6),difftemp(0:nterms+6),
     1     powtemp(0:nterms+6),
     2     diffs(0:nterms+6),ival(0:nterms+6),pow(0:nterms+6),
     3     dpow(0:nterms+6))


c     get local difference and BesselI vals
      
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)
      theta=theta-pi
      z=zk*r
      ifder=0

      call diffszkik_fast(r,beta,rscale1,diffs,ifder,ders,ival,
     1     nterms)

      rtemp = r
      call mbh2d_rksc(pow,dpow,rtemp,beta,rscale1,nterms)

c     
      jtemp(0) = ival(0)
      difftemp(0) = diffs(0)
      powtemp(0) = pow(0)
      zmul=exp(ima*theta)
      ztemp1= zmul

      do j = 1,nterms
         jtemp( j) = ztemp1*ival(j)
         difftemp(j) = ztemp1*diffs(j)
         powtemp(j) = ztemp1*pow(j)
         ztemp1= ztemp1*zmul
      enddo

c     shift the expansion

      rsj = rscale1

      do l=1,nd
         ympole2(l,0) = ympole2(l,0) + dreal(ympole1(l,0))*jtemp(0)
         ympole2(l,0) = ympole2(l,0) + dreal(mbhmpole1(l,0))*difftemp(0)
         mbhmpole2(l,0) = mbhmpole2(l,0) + dreal(mbhmpole1(l,0))
      enddo

      do j = 1,nterms1
         do l = 1,nd
            ympole2(l,0) = ympole2(l,0)
     1           +(dreal(ympole1(l,j))*dreal(jtemp(j))
     1           +dimag(ympole1(l,j))*dimag(jtemp(j)))*rsj**2
            ympole2(l,0) = ympole2(l,0)
     1           +(dreal(mbhmpole1(l,j))*dreal(jtemp(j))
     1           +dimag(mbhmpole1(l,j))*dimag(jtemp(j)))*rsj**2
         enddo
         rsj=rsj*rscale1
      enddo

c     
      rsi=rscale1

      rsi5=rscale1/rscale2
      do i = 1,nterms2
         do l = 1,nd
            ympole2(l,i) = ympole2(l,i) 
     1           + 2.0d0*(dreal(ympole1(l,0))*dreal(jtemp(i))
     2           +ima*dreal(ympole1(l,0))*dimag(jtemp(i)))*rsi5
            ympole2(l,i) = ympole2(l,i) 
     1           + 2.0d0*(dreal(mbhmpole1(l,0))*dreal(difftemp(i))
     2           +ima*dreal(mbhmpole1(l,0))*dimag(difftemp(i)))*rsi5
            mbhmpole2(l,i) = mbhmpole2(l,i) 
     1           + 2.0d0*(dreal(mbhmpole1(l,0))*dreal(powtemp(i))
     2           +ima*dreal(mbhmpole1(l,0))*dimag(powtemp(i)))*rsi5
         enddo
         rsj=rscale1
         twojm1 = 1.0d0
         do j = 1,min(nterms1,i)
            do l = 1,nd
               ympole2(l,i) = ympole2(l,i)
     1              +(dreal(ympole1(l,j))*dreal(jtemp(i-j))
     1              -dimag(ympole1(l,j))*dimag(jtemp(i-j))
     2              +ima*(dreal(ympole1(l,j))*dimag(jtemp(i-j))
     3              +dimag(ympole1(l,j))*dreal(jtemp(i-j))))*rsi5
               ympole2(l,i) = ympole2(l,i)
     1              +(dreal(ympole1(l,j))*dreal(jtemp(i+j))
     1              +dimag(ympole1(l,j))*dimag(jtemp(i+j))
     2              +ima*(dreal(ympole1(l,j))*dimag(jtemp(i+j))
     3              -dimag(ympole1(l,j))*dreal(jtemp(i+j))))*rsj**2*rsi5
               mbhmpole2(l,i) = mbhmpole2(l,i)
     1              +(dreal(mbhmpole1(l,j))*dreal(powtemp(i-j))
     1              -dimag(mbhmpole1(l,j))*dimag(powtemp(i-j))
     2              +ima*(dreal(mbhmpole1(l,j))*dimag(powtemp(i-j))
     3              +dimag(mbhmpole1(l,j))*dreal(powtemp(i-j))))
     5              *rsi5
               ztemp = difftemp(i-j)
               ympole2(l,i) = ympole2(l,i)
     1              +(dreal(mbhmpole1(l,j))*dreal(ztemp)
     1              -dimag(mbhmpole1(l,j))*dimag(ztemp)
     2              +ima*(dreal(mbhmpole1(l,j))*dimag(ztemp)
     3              +dimag(mbhmpole1(l,j))*dreal(ztemp)))*rsi5
               ympole2(l,i) = ympole2(l,i)
     1              +(dreal(mbhmpole1(l,j))*dreal(jtemp(i+j))
     1              +dimag(mbhmpole1(l,j))*dimag(jtemp(i+j))
     2              +ima*(dreal(mbhmpole1(l,j))*dimag(jtemp(i+j))
     3              -dimag(mbhmpole1(l,j))*dreal(jtemp(i+j))))
     5              *rsj**2*rsi5
            enddo
            rsj=rsj*rscale1
            twojm1 = twojm1*2.0d0
         enddo
         rsj=rscale1**(i+1)
         fs2=rsi5*rscale1**2
         do j = i+1,nterms1
            do l=1,nd
               ympole2(l,i) = ympole2(l,i)
     1              +(dreal(mbhmpole1(l,j))*dreal(jtemp(j-i))
     1              +dimag(mbhmpole1(l,j))*dimag(jtemp(j-i))
     2              -ima*(dreal(mbhmpole1(l,j))*dimag(jtemp(j-i))
     3              -dimag(mbhmpole1(l,j))*dreal(jtemp(j-i))))*fs2
               ympole2(l,i) = ympole2(l,i)
     1              +(dreal(mbhmpole1(l,j))*dreal(jtemp(i+j))
     1              +dimag(mbhmpole1(l,j))*dimag(jtemp(i+j))
     2              +ima*(dreal(mbhmpole1(l,j))*dimag(jtemp(i+j))
     3              -dimag(mbhmpole1(l,j))*dreal(jtemp(i+j))))
     5              *rsj**2*rsi5
               ympole2(l,i) = ympole2(l,i)
     1              +(dreal(ympole1(l,j))*dreal(jtemp(j-i))
     1              +dimag(ympole1(l,j))*dimag(jtemp(j-i))
     2              -ima*(dreal(ympole1(l,j))*dimag(jtemp(j-i))
     3              -dimag(ympole1(l,j))*dreal(jtemp(j-i))))*fs2
               ympole2(l,i) = ympole2(l,i)
     1              +(dreal(ympole1(l,j))*dreal(jtemp(i+j))
     1              +dimag(ympole1(l,j))*dimag(jtemp(i+j))
     2              +ima*(dreal(ympole1(l,j))*dimag(jtemp(i+j))
     3              -dimag(ympole1(l,j))*dreal(jtemp(i+j))))
     5              *rsj**2*rsi5
            enddo
            rsj=rsj*rscale1
            fs2=fs2*rscale1**2
         enddo
         rsi=rsi*rscale1
         rsi5=rsi5*rscale1/rscale2
      enddo

      return
      end


      subroutine mbh2dmploc_vec(nd,beta,rscale1,center1,mbhmpole,ympole,
     1     nterms1,rscale2,center2,mbhloc,lloc,nterms2)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT: 
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     mbhmpole(nd,0:nterms1): original coeffs, difference-type functions
c     ympole(nd,0:nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     nterms2           : number of terms in original expansion
c     zcirc2(2,nterms2) : precomp'd evenly-spaced points on unit circle
c     beta              : the modified biharmonic parameter

c     work(*)           : work array. recommended length 16*nterms+500
c
c     OUTPUT:
c
c     mbhloc(nd,0:nterms2): new coeffs, difference-type functions
c     lloc(nd,0:nterms2): new coeffs, power series type
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      complex *16 mbhmpole(nd,0:nterms1), ympole(nd,0:nterms1)
      complex *16 mbhloc(nd,0:nterms2), lloc(nd,0:nterms2)
      real *8 rscale1, rscale2, center1(2), center2(2), beta
      integer nterms1, nterms2, nd
c     local
      real *8 zdiff(2), r, theta, pi, rsi, rsj, rsi5, rsi7, rsj2
      real *8 rsi52, rsi2, ders(0:1), dtmp1, dtmp2, dtmp3, dtmp4
      real *8, allocatable :: diffs(:), kvec(:), pow(:), dpow(:),
     1     dfac2(:)
      complex *16, allocatable :: ktemp(:), difftemp(:), powtemp(:)
      complex *16 ztemp1, zmul, ima
      integer nterms, i, j, l, isign, ifders
      data ima /(0.0d0,1.0d0)/

      
      pi = 4.0d0*datan(1.0d0)

      nterms = nterms1+nterms2

      allocate(diffs(0:nterms+6),kvec(0:nterms+6),pow(0:nterms+6),
     1     dpow(0:nterms+6),dfac2(0:nterms+6),ktemp(0:nterms+6),
     2     difftemp(0:nterms+6),powtemp(0:nterms+6))
      

      dfac2(0) = 1.0d0
      do i = 1,nterms
         dfac2(i) = dfac2(i-1)/(2.0d0*i)
      enddo

      zdiff(1) = center2(1)-center1(1)
      zdiff(2) = center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)

      theta = theta-pi

c     get values of (beta*r)^-k
      call mbh2d_rmk(pow,dpow,r,beta,rscale1,nterms)

c     get values of difference functions
      ifders = 0
      call diffslogbk_fast(r,beta,rscale1,diffs,ifders,ders,kvec,
     1     nterms)

c     form terms for shifting 

      ktemp(0) = kvec(0)
      difftemp(0) = diffs(0)
      powtemp(0) = pow(0)
      zmul=-exp(ima*theta)
      ztemp1= zmul
      do j = 1,nterms
         ktemp( j) = ztemp1*kvec(j)
         difftemp(j) = ztemp1*diffs(j)
         powtemp(j) = ztemp1*pow(j)
         ztemp1= ztemp1*zmul
      enddo

c     shift

      do l = 1,nd
         mbhloc(l,0) = mbhloc(l,0) + mbhmpole(l,0)*ktemp(0)
         mbhloc(l,0) = mbhloc(l,0) + ympole(l,0)*ktemp(0)
         lloc(l,0) = lloc(l,0) + mbhmpole(l,0)*(ktemp(0)+pow(0))
         lloc(l,0) = lloc(l,0) + ympole(l,0)*ktemp(0)
      enddo
      do j = 1,nterms1
         do l = 1,nd
            mbhloc(l,0) = mbhloc(l,0) +
     1           (dreal(mbhmpole(l,j))*dreal(ktemp(j))
     1           +dimag(mbhmpole(l,j))*dimag(ktemp(j)))
            mbhloc(l,0) = mbhloc(l,0) +
     1           (dreal(ympole(l,j))*dreal(ktemp(j))
     1           +dimag(ympole(l,j))*dimag(ktemp(j)))
            lloc(l,0) = lloc(l,0) +
     1           (dreal(mbhmpole(l,j))*dreal(difftemp(j))
     1           +dimag(mbhmpole(l,j))*dimag(difftemp(j)))
            lloc(l,0) = lloc(l,0) +
     1           (dreal(ympole(l,j))*dreal(ktemp(j))
     1           +dimag(ympole(l,j))*dimag(ktemp(j)))
         enddo
      enddo
c     
      rsi=rscale1
      rsi2=rscale1**2
      rsi5=rscale2/rscale1
      rsi52=rsi5*rscale2/rscale1
      isign = -1
      do i = 1,nterms2
         dtmp1 = isign*rsi5*2
         dtmp2 = dtmp1*dfac2(i)
         do l = 1,nd
            mbhloc(l,i) = mbhloc(l,i) + dtmp1*
     1           dreal(mbhmpole(l,0))*ktemp(i)
            mbhloc(l,i) = mbhloc(l,i) + dtmp1*
     1           dreal(ympole(l,0))*ktemp(i)
            lloc(l,i) = lloc(l,i) + dtmp2*
     1           dreal(mbhmpole(l,0))*difftemp(i)
            lloc(l,i) = lloc(l,i) + dtmp2*
     1           dreal(ympole(l,0))*ktemp(i)
         enddo
         rsj=rscale1
         rsj2=rscale1**2
         dtmp1=isign*rsi5         
         dtmp2=dtmp1*dfac2(i)
         do j = 1,min(nterms1,i)

            dtmp3=dtmp1*rsj2
            dtmp4=dtmp2*rsj2

            do l = 1,nd
               mbhloc(l,i) = mbhloc(l,i) + dtmp3*
     1              (mbhmpole(l,j)+ympole(l,j))*ktemp(i-j)
               mbhloc(l,i) = mbhloc(l,i) + dtmp1*
     1              dconjg(mbhmpole(l,j)+ympole(l,j))*ktemp(i+j)
               lloc(l,i) = lloc(l,i) + dtmp4*
     1              (mbhmpole(l,j)+ympole(l,j))*ktemp(i-j)
               lloc(l,i) = lloc(l,i) + dtmp2*
     1              dconjg(mbhmpole(l,j))*difftemp(i+j)
               lloc(l,i) = lloc(l,i) + dtmp2*
     1              dconjg(ympole(l,j))*ktemp(i+j)
            enddo
            rsj=rsj*rscale1
            rsj2=rsj2*rscale1**2
         enddo
         dtmp3=dtmp1*rsi2
         dtmp4=dtmp2*rsi2
         do j = i+1,nterms1
            do l = 1,nd
               mbhloc(l,i) = mbhloc(l,i) + dtmp3*
     1              (mbhmpole(l,j)+ympole(l,j))*dconjg(ktemp(j-i))
               mbhloc(l,i) = mbhloc(l,i) + dtmp1*
     1              dconjg(mbhmpole(l,j)+ympole(l,j))*ktemp(i+j)
               lloc(l,i) = lloc(l,i) + dtmp4*
     1              (mbhmpole(l,j)+ympole(l,j))*dconjg(ktemp(j-i))
               lloc(l,i) = lloc(l,i) + dtmp2*
     1              dconjg(mbhmpole(l,j))*difftemp(i+j)
               lloc(l,i) = lloc(l,i) + dtmp2*
     1              dconjg(ympole(l,j))*ktemp(i+j)
            enddo
         enddo
         isign = -isign
         rsi=rsi*rscale1
         rsi2=rsi2*rscale1**2
         rsi5=rsi5*rscale2/rscale1
      enddo

      return
      end


      subroutine mbh2dlocloc_vec(nd,beta,rscale1,center1,mbhloc1,lloc1,
     1     nterms1,rscale2,center2,mbhloc2,lloc2,nterms2,carray,ldc)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT: 
c
c     rscale1           : scaling for original expansion
c     center1           : center for original expansion
c     mbhloc1(nd,0:nterms1): original coeffs, difference-type functions
c     lloc1(nd,0:nterms1): original coeffs, modified Bessel functions
c     nterms1           : number of terms in original expansion
c     rscale2           : scaling for new expansion
c     center2           : center for new expansion
c     nterms2           : number of terms in original expansion
c     beta              : the modified biharmonic parameter
c     carray(0:ldc,0:ldc) : precomputed array of binomial coefficients
c     ldc               : dimension of carray
c
c     OUTPUT:
c
c     mbhloc2(nd,0:nterms2): new coeffs, difference-type functions
c     lloc2(nd,0:nterms2): new coeffs, power series-like functions
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 rscale1, center1(2), beta, rscale2, center2(2)
      complex *16 mbhloc1(nd,0:nterms1), lloc1(nd,0:nterms1)
      complex *16 mbhloc2(nd,0:nterms2), lloc2(nd,0:nterms2)
      real *8 carray(0:ldc,0:ldc)
      integer ldc
      integer nterms1, nterms2, nd
c     local
      complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
      complex *16, allocatable :: jtemp(:), difftemp(:), powtemp(:)
      real *8 zdiff(2), r, theta, pi, done
      real *8 rsj, rsi, rsi7, rsi5, fs2, ders(0:1)
      real *8, allocatable :: ival(:), diffs(:), pow(:), dpow(:),
     1     pow2(:), dfac2(:)
      integer i, j, m, l, nterms, ifder
c     
      data ima/(0.0d0,1.0d0)/
c     
      done=1
      pi=4*atan(done)
c     
      nterms = nterms1+nterms2
c

      allocate(jtemp(0:nterms+6),difftemp(0:nterms+6),
     1     powtemp(0:nterms+6),ival(0:nterms+6),diffs(0:nterms+6),
     1     pow(0:nterms+6),dpow(0:nterms+6),pow2(0:nterms+6),
     1     dfac2(0:nterms+6))
      
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)
      theta=theta-pi

      ifder = 0
      call diffszkik_fast(r,beta,rscale1,diffs,ifder,ders,ival,
     1     nterms)

      call mbh2d_rk(pow,dpow,r,beta,rscale1,nterms)
      call mbh2d_rksc(pow2,dpow,r,beta,rscale1,nterms)

      dfac2(0) = 1.0d0
      do i = 1,nterms
         dfac2(i) = dfac2(i-1)/(2.0d0*i)
      enddo
c     
      jtemp(0) = ival(0)
      powtemp(0) = pow(0)
      difftemp(0) = diffs(0)
      zmul=-exp(ima*theta)
      ztemp1= zmul
      do j = 1,nterms
         jtemp( j) = ztemp1*ival(j)
         difftemp(j) = ztemp1*diffs(j)
         powtemp(j) = ztemp1*pow(j)
         ztemp1= ztemp1*zmul
      enddo
c
      do l = 1,nd
         mbhloc2(l,0) = mbhloc2(l,0) + mbhloc1(l,0)*jtemp(0)
         lloc2(l,0) = lloc2(l,0) + mbhloc1(l,0)*difftemp(0)
         lloc2(l,0) = lloc2(l,0) + lloc1(l,0)*powtemp(0)
      enddo
      
      rsj=rscale1
      do j = 1,nterms1
         do l = 1,nd
            mbhloc2(l,0) = mbhloc2(l,0)+
     1           (dreal(mbhloc1(l,j))*dreal(jtemp(j))
     1           +dimag(mbhloc1(l,j))*dimag(jtemp(j)))
            lloc2(l,0) = lloc2(l,0)+
     1           (dreal(mbhloc1(l,j))*dreal(difftemp(j))
     1           +dimag(mbhloc1(l,j))*dimag(difftemp(j)))
            lloc2(l,0) = lloc2(l,0)+
     1           (dreal(lloc1(l,j))*dreal(powtemp(j))
     1           +dimag(lloc1(l,j))*dimag(powtemp(j)))
         enddo
      enddo
c     
      rsi=rscale1
      rsi7=rscale2
      rsi5=rscale2/rscale1
      do i = 1,nterms2
         do l = 1,nd
            mbhloc2(l,i) = mbhloc2(l,i) 
     1           + 2.0d0*(dreal(mbhloc1(l,0))*dreal(jtemp(i))
     2           +ima*dreal(mbhloc1(l,0))*dimag(jtemp(i)))*rsi7*rsi
            lloc2(l,i) = lloc2(l,i) 
     1           + 2.0d0*(dreal(mbhloc1(l,0))*dreal(jtemp(i))
     2           +ima*dreal(mbhloc1(l,0))*dimag(jtemp(i)))*rsi7*rsi*
     1           dfac2(i)
         enddo
         fs2=rsi5*rscale1**2
         if( nterms1 .le. i-1 ) fs2=fs2*rscale1**(2*(i-1-nterms1))
         do j = min(nterms1,i-1),1,-1
            do l = 1,nd
               mbhloc2(l,i) = mbhloc2(l,i)+
     1              (dreal(mbhloc1(l,j))*dreal(jtemp(i-j))
     1              -dimag(mbhloc1(l,j))*dimag(jtemp(i-j))
     2              +ima*(dreal(mbhloc1(l,j))*dimag(jtemp(i-j))
     3              +dimag(mbhloc1(l,j))*dreal(jtemp(i-j))))*fs2
               mbhloc2(l,i) = mbhloc2(l,i)+
     1              (dreal(mbhloc1(l,j))*dreal(jtemp(i+j))
     1              +dimag(mbhloc1(l,j))*dimag(jtemp(i+j))
     2              +ima*(dreal(mbhloc1(l,j))*dimag(jtemp(i+j))
     3              -dimag(mbhloc1(l,j))*dreal(jtemp(i+j))))*rsi*rsi7
               lloc2(l,i) = lloc2(l,i)+
     1              (dreal(mbhloc1(l,j))*dreal(jtemp(i-j))
     1              -dimag(mbhloc1(l,j))*dimag(jtemp(i-j))
     2              +ima*(dreal(mbhloc1(l,j))*dimag(jtemp(i-j))
     3              +dimag(mbhloc1(l,j))*dreal(jtemp(i-j))))*fs2*
     5              dfac2(i)
               lloc2(l,i) = lloc2(l,i)+
     1              (dreal(mbhloc1(l,j))*dreal(jtemp(i+j))
     1              +dimag(mbhloc1(l,j))*dimag(jtemp(i+j))
     2              +ima*(dreal(mbhloc1(l,j))*dimag(jtemp(i+j))
     3              -dimag(mbhloc1(l,j))*dreal(jtemp(i+j))))
     4              *rsi*rsi7*dfac2(i)
            enddo
            fs2=fs2*rscale1**2
         enddo
         do j = i,nterms1
            do l = 1,nd
               mbhloc2(l,i) = mbhloc2(l,i)+
     1              (dreal(mbhloc1(l,j))*dreal(jtemp(j-i))
     1              +dimag(mbhloc1(l,j))*dimag(jtemp(j-i))
     2              -ima*(dreal(mbhloc1(l,j))*dimag(jtemp(j-i))
     3              -dimag(mbhloc1(l,j))*dreal(jtemp(j-i))))*rsi5
               lloc2(l,i) = lloc2(l,i)+
     1              (dreal(mbhloc1(l,j))*dreal(difftemp(j-i))
     1              +dimag(mbhloc1(l,j))*dimag(difftemp(j-i))
     2              -ima*(dreal(mbhloc1(l,j))*dimag(difftemp(j-i))
     3              -dimag(mbhloc1(l,j))*dreal(difftemp(j-i))))*rsi5
     1              *dfac2(i)
               lloc2(l,i) = lloc2(l,i)+
     1              (dreal(lloc1(l,j))*dreal(powtemp(j-i))
     1              +dimag(lloc1(l,j))*dimag(powtemp(j-i))
     2              -ima*(dreal(lloc1(l,j))*dimag(powtemp(j-i))
     3              -dimag(lloc1(l,j))*dreal(powtemp(j-i))))*carray(j,i)
     1              *rsi5
               mbhloc2(l,i) = mbhloc2(l,i)+
     1              (dreal(mbhloc1(l,j))*dreal(jtemp(i+j))
     1              +dimag(mbhloc1(l,j))*dimag(jtemp(i+j))
     2              +ima*(dreal(mbhloc1(l,j))*dimag(jtemp(i+j))
     3              -dimag(mbhloc1(l,j))*dreal(jtemp(i+j))))*rsi*rsi7
               lloc2(l,i) = lloc2(l,i)+
     1              (dreal(mbhloc1(l,j))*dreal(jtemp(i+j))
     1              +dimag(mbhloc1(l,j))*dimag(jtemp(i+j))
     2              +ima*(dreal(mbhloc1(l,j))*dimag(jtemp(i+j))
     3              -dimag(mbhloc1(l,j))*dreal(jtemp(i+j))))
     4              *rsi*rsi7*dfac2(i)
               
            enddo
         enddo
         rsi=rsi*rscale1
         rsi7=rsi7*rscale2
         rsi5=rsi5*rscale2/rscale1
      enddo
      return
      end

      subroutine mbh2d_rk(pow,dpow,r,beta,rscale,nterms)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     r               : radius
c     rscale          : scaling factor
c     nterms          : number of terms in the expansion
c     
c     OUTPUT:
c
c     pow(0:nterms)   : pow(i) = (r/rscale)^i
c     dpow(0:nterms)  : dpow(i) = i*(r/rscale)^(i-1)/rscale
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 pow(0:nterms), dpow(0:nterms), r, rscale, beta
      integer nterms
c     local
      real *8 dtemp1, dtemp2
      integer i

      dtemp2 = (r*beta)/rscale
      dtemp1 = 1.0d0

      pow(0) = 1.0d0
      dpow(0) = 0.0d0

      do i = 1,nterms
         dpow(i) = i*beta*dtemp1/rscale
         dtemp1 = dtemp1*dtemp2
         if (dtemp1 .gt. 1.0d250) dtemp1 = 0.0d0
         pow(i) = dtemp1
      enddo

      return
      end

      subroutine mbh2d_rksc(pow,dpow,r,beta,rscale,nterms)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     r               : radius
c     rscale          : scaling factor
c     nterms          : number of terms in the expansion
c     
c     OUTPUT:
c
c     pow(0:nterms)   : pow(i) = (r*beta/(2*rscale))^i/(i!)
c     dpow(0:nterms)  : dpow(i) = i*(r/rscale)^(i-1)/rscale
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 pow(0:nterms), dpow(0:nterms), r, rscale, beta
      integer nterms
c     local
      real *8 dtemp1, dtemp2
      integer i

      dtemp2 = (r*beta)/(2.0d0*rscale)
      dtemp1 = 1.0d0

      pow(0) = 1.0d0
      dpow(0) = 0.0d0

      do i = 1,nterms
         dpow(i) = beta*dtemp1/(2.0d0*rscale)
         dtemp1 = dtemp1*dtemp2/i
         if (dtemp1 .gt. 1.0d250) dtemp1 = 0.0d0
         pow(i) = dtemp1
      enddo

      return
      end

      subroutine mbh2d_rmk(pow,dpow,r,beta,rscale,nterms)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     INPUT:
c
c     r               : radius
c     rscale          : scaling factor
c     nterms          : number of terms in the expansion
c     
c     OUTPUT:
c
c     pow(0:nterms)   : pow(i) = (rscale/r)^(-i)
c     dpow(0:nterms)  : dpow(i) = -i*(rscale/r)^(-i-1)*rscale
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global
      real *8 pow(0:nterms), dpow(0:nterms), r, rscale, beta
      integer nterms
c     local
      real *8 dtemp1, dtemp2
      integer i

      dtemp2 = rscale/(r*beta)
      dtemp1 = dtemp2

      pow(0) = dlog(r)
      dpow(0) = 1.0d0/r

      do i = 1,nterms
         pow(i) = dtemp1
         if (dtemp1 .lt. 1.0d-200) dtemp1 = 0.0d0
         dtemp1 = dtemp1*dtemp2
         dpow(i) = -i*dtemp1*rscale/beta
      enddo

      return
      end

c
c


      subroutine mbh2d_init_carray(carray,ldc)
      implicit real *8 (a-h,o-z)
      real *8 carray(0:ldc,0:ldc)
      
      do l = 0,ldc
         carray(l,0) = 1.0d0
      enddo
      do m=1,ldc
         carray(m,m) = 1.0d0
         do l=m+1,ldc
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
         enddo
      enddo
c
      return
      end
        
      
