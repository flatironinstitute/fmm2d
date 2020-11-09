cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
cc
c
c
c Modified: (C) 2020, Manas Rachh
c
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c    $Date$
c    $Revision$
c
c
c      This file contains the basic subroutines for 
c      forming and evaluating multipole (partial wave) expansions
c      in two dimensions.
c
c      Since log(z) is a multivalued complex function, we use
c      the real part Re(log(z)) = log(abs(z)) in all computations.
c
c      All multipole and local expansions are properly scaled 
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
c      On the other hand, the direct evaluation routines from particles in 
c      l2d_DIRECT or l2d_DIRECT1 return the result of the requested
c      computation by overwriting the input pot/grad/hess variables.
c-----------------------------------------------------------------------
c
c
c      L2DFORMMPC_VEC:  creates multipole expansions (outgoing) due to 
c                       a collection of charge sources.
c      L2DFORMMPD_VEC:  creates multipole expansions (outgoing) due to 
c                       a collection of dipole sources.
c      L2DFORMMPCD_VEC: creates multipole expansions (outgoing) due to 
c                       a collection of charge and dipole sources.
c      L2DMPEVALP_VEC:  computes potentials from multipole expansions
c                       at a collection of targets
c      L2DMPEVALG_VEC:  computes potentials/gradients from multipole expansions
c                       at a collection of targets
c      L2DMPEVALH_VEC:  computes potentials/gradients/Hessians 
c                       from multipole expansions at a collection of targets
c
c      L2DFORMTAC_VEC:  creates local expansions due to 
c                       a collection of charge sources
c      L2DFORMTAD_VEC:  creates local expansions due to 
c                       a collection of dipole sources
c      L2DFORMTACD_VEC: creates local expansions due to 
c                       a collection of charge and dipole sources
c      L2DTAEVALP_VEC:  computes potentials from local expansions
c                       at a collection of targets
c      L2DTAEVALG_VEC:  computes potentials/gradients from local expansions
c                       at a collection of targets
c      L2DTAEVALH_VEC:  computes potentials/gradients/Hessians 
c                       from local expansions at a collection of targets
c 
c      L2DMPMP_VEC:     Translates center of multipole expansions 
c      L2DLOCLOC_VEC:   Translates center of local expansions 
c      L2DMPLOC_VEC:    Converts multipole expansions to local expansions
c
c      L2DMPZERO_VEC:   utility to initialize multipole coefficients to zero.
C
c-----------------------------------------------------------------------
c
c
C***********************************************************************
      subroutine l2dformmpc_vec(nd,rscale,source,ns,charge,
     1                      center,nterms,mpole)
      implicit none
C***********************************************************************
c
c     This subroutine INCREMENTS multipole (outgoing) expansions about 
c     CENTER due to (ND,NS) charges located at SOURCES(2,NS).
c
c     mpole_0(q) = sum  charge_j(q)
c                   j
c
c     mpole_n(q)=mpole_n(q)-
c
c                sum charge_j(q) 1/n z0^n/rscale^n
c                 j  
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of mpole expansions)
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     ns              : number of sources
c     charge(nd,ns)      : source strengths
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mpole(nd,0:*)     : coeffs for the outgoing expansion are 
c                         incremented.
c-----------------------------------------------------------------------
      integer ns,j,ifder,n,nterms,ii,nd
      real *8 center(2),source(2,ns),zdiff(2),rscale
      complex *16 mpole(nd,0:nterms),charge(nd,ns)
      complex *16 ztemp1,ztemp2,z0,zpow(0:nterms)
      complex *16 z
      
      
      
      
      
c
c
      do j=1,ns
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)

        z0 = dcmplx(zdiff(1),zdiff(2))
        ztemp1 = z0/rscale

        zpow(0) = -1
        do n=1,nterms
          zpow(n) = zpow(n-1)*ztemp1
        enddo
        do n=1,nterms
          zpow(n) = zpow(n)/(n+0.0d0)
        enddo
        zpow(0) = 1
        do n=0,nterms
          do ii=1,nd
            mpole(ii,n) = mpole(ii,n) + charge(ii,j)*zpow(n)
          enddo
        enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l2dformmpd_vec(nd,rscale,source,ns,
     1                      dipstr,center,nterms,mpole)
      implicit none
C***********************************************************************
c
c     This subroutine INCREMENTS multipole (outgoing) expansions about 
c     CENTER due to NS dipoles located at SOURCES(2,*).
c
c     for dipoles
c      
c     mpole_0  =  mpole_0
c                 
c
c     mpole_n  =  mpole_n + sum dipstr_j z0^(n-1) /rscale^n
c                            j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              : vector length (number of mpole expansions)
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     ns              : number of sources
c     dipstr(nd,ns)   : dipole strengths
c
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mpole(nd,*)     : coeffs for the multipole-expansion 
c                       are incremented.
c-----------------------------------------------------------------------
      integer ns,j,ifder,n,nterms,ii,nd
      real *8 center(2),source(2,ns),zdiff(2),rscale
      real *8 r,theta
      complex *16 mpole(nd,0:nterms),dipstr(nd,ns)
      complex *16 zmul,zinv,ztemp1,ztemp2,zpow(nterms)
      complex *16 z,z0
c
c
      do j=1,ns
         zdiff(1)=source(1,j)-center(1)
         zdiff(2)=source(2,j)-center(2)
         z0=dcmplx(zdiff(1),zdiff(2))
         ztemp1 = z0/rscale
         zpow(1) = 1/rscale
         do n=2,nterms
            zpow(n) = zpow(n-1)*ztemp1
         enddo
         do n=1,nterms
            do ii=1,nd
               mpole(ii,n) = mpole(ii,n) + dipstr(ii,j)*zpow(n)
            enddo
         enddo
      enddo
      return
      end
c
c
c
c
C***********************************************************************
      subroutine l2dformmpcd_vec(nd,rscale,source,ns,charge,
     1                       dipstr,center,nterms,mpole)
      implicit none
C***********************************************************************
c
c     This subroutine INCREMENTS multipole (outgoing) expansions about 
c     CENTER due to NS sources located at SOURCES(2,*).
c
c     mpole_0  =  mpole_0 + sum  charge_j
c                            j
c
c     mpole_n  =  mpole_n + sum (dipstr_j z0^(n-1) - charge_j/n z0^n)/rscale^n
c                            j  
c
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of mpole expansions)
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     ns              : number of sources
c     charge(nd,ns)   : source strengths
c     dipstr(nd,ns)   : dipole strengths
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mpole(nd,*)     : coeffs for the multipole-expansion are incremented.
c-----------------------------------------------------------------------
      integer ns,j,ifder,n,nterms,ii,nd
      real *8 center(2),source(2,ns),zdiff(2),rscale
      real *8 r,theta
      complex *16 mpole(nd,0:nterms),dipstr(nd,ns),charge(nd,ns)
      complex *16 ztemp1,ztemp2,zpowd(0:nterms),zpowc(0:nterms)
      complex *16 z,z0
c
c
c
      do j=1,ns
         zdiff(1)=source(1,j)-center(1)
         zdiff(2)=source(2,j)-center(2)
         z0=dcmplx(zdiff(1),zdiff(2))
         ztemp1 = z0/rscale
         zpowd(0) = 0
         zpowd(1) = 1/rscale
         zpowc(0) = 1
         zpowc(1) = -ztemp1
         do n=2,nterms
            zpowd(n) = zpowd(n-1)*ztemp1
            zpowc(n) = zpowc(n-1)*ztemp1
         enddo
         do n=1,nterms
           zpowc(n) = zpowc(n)/(n+0.0d0)
         enddo
         do n=0,nterms
            do ii=1,nd
               mpole(ii,n) = mpole(ii,n) + dipstr(ii,j)*zpowd(n) + 
     1            charge(ii,j)*zpowc(n)
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
c**********************************************************************
      subroutine l2dmpevalp_vec(nd,rscale,center,mpole,nterms,
     1                     ztarg,ntarg,pot1)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials due to
c     outgoing multipole expansions.
c                    +nterms
c     pot1  = pot1 +     sum   mpole_n (rscale/z)^n + mpole_0 log(abs(z))
c                      n=1  
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of mpole expansions)
c     zk           :    Helmholtz parameter
c     rscale       :    scaling parameter 
c     center       :    expansion center
c     mpole(nd,*)  :    multipole expansion 
c     nterms       :    order of the multipole expansion
c     ztarg        :    target locations
c     ntarg        :    number of targets
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot1(nd,ntarg) :   potentials at ztarg locations are incremented
c-----------------------------------------------------------------------
      integer k,ntarg,j,n,nterms,ifder,ii,nd
      real *8 rscale,r,theta,center(2),ztarg(2,ntarg),zdiff(2)
      complex *16 pot,mpole(nd,0:nterms),zpow(0:nterms)
      complex *16 pot1(nd,ntarg)
      complex *16 z,ztemp1,ztemp2,z0
c
c
c
cc      call prin2('pot1=*',pot1,2*nd*ntarg)
cc      call prin2('mpole=*',mpole,2*nd*(nterms+1))
      do k=1,ntarg
         zdiff(1)=ztarg(1,k)-center(1)
         zdiff(2)=ztarg(2,k)-center(2)
         z0 = dcmplx(zdiff(1),zdiff(2))
         zpow(0) = log(abs(z0))
         ztemp1 = rscale/z0
         zpow(1) = ztemp1
         do n=2,nterms
           zpow(n) = zpow(n-1)*ztemp1
         enddo
cc         call prin2('zpow=*',zpow,2*nterms+2)
         do n=0,nterms
           do ii=1,nd
             pot1(ii,k) = pot1(ii,k) + mpole(ii,n)*zpow(n)
           enddo
         enddo
      enddo

      return
      end
c
c
c
c
c**********************************************************************
      subroutine l2dmpevalg_vec(nd,rscale,center,mpole,nterms,
     1                     ztarg,ntarg,pot1,grad1)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials and gradients due to
c     outgoing multipole expansions.
c                    +nterms
c     pot1  = pot1 +     sum   mpole_n (rscale/z)^n + mpole_0 log(abs(z))
c                      n=1  
c     grad = d/dz(pot1) 
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of mpole expansions)
c     rscale       :    scaling parameter 
c     center       :    expansion center
c     mpole(nd,*)  :    multipole expansion 
c     nterms       :    order of the multipole expansion
c     ztarg        :    target locations
c     ntarg        :    number of targets
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot1(nd,*)    :   potentials at ztarg locations are incremented
c     grad1(nd,*) :   gradients at ztarg locations are incremented
c-----------------------------------------------------------------------
      integer k,ntarg,j,n,nterms,i,ifder,ii,nd
      real *8 rscale,center(2),ztarg(2,ntarg),zdiff(2),rinv
      complex *16 pot,grad,mpole(nd,0:nterms)
      complex *16 pot1(nd,ntarg),grad1(nd,ntarg)
      complex *16 zpow(0:nterms+3)
      complex *16 zpowg(0:nterms+3)
c
c
      complex *16 z,z1scale,z2scale,z3scale,z4scale
      complex *16 ztemp1,ztemp2
c
      rinv = 1/rscale
      do k=1,ntarg
         zdiff(1)=ztarg(1,k)-center(1)
         zdiff(2)=ztarg(2,k)-center(2)
         z = dcmplx(zdiff(1),zdiff(2))
         zpow(0) = log(abs(z))
         ztemp1 = rscale/z
         zpow(1) = ztemp1
         do n=2,nterms+1
           zpow(n) = zpow(n-1)*ztemp1
         enddo
         zpowg(0) = zpow(1)*rinv
         do n=1,nterms
           zpowg(n) = -zpow(n+1)*rinv*n
         enddo
         do n=0,nterms
           do ii=1,nd
             pot1(ii,k) = pot1(ii,k) + mpole(ii,n)*zpow(n)
             grad1(ii,k) = grad1(ii,k) + mpole(ii,n)*zpowg(n)
           enddo
         enddo
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine l2dmpevalh_vec(nd,rscale,center,mpole,nterms,
     1                     ztarg,ntarg,pot1,grad1,hess1)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials, gradients, and 
c     hessians due to  outgoing multipole expansions.
c                    +nterms
c     pot1  = pot1 +     sum   mpole_n (rscale/z)^n + mpole_0 log(abs(z))
c                      n=1  
c     grad1 = d/dz(pot1)
c     hess1 = d/dz(grad1)
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of mpole expansions)
c     rscale       :    scaling parameter 
c     center       :    expansion center
c     mpole(nd,*)  :    multipole expansion 
c     nterms       :    order of the multipole expansion
c     ztarg        :    target locations
c     ntarg        :    number of targets
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot1(nd,*)   :   potentials at ztarg locations are incremented
c     grad1(nd,*)  :   gradients at ztarg locations are incremented
c     hess1(nd,*)  :   hessians at ztarg locations are incremented
c-----------------------------------------------------------------------
      integer k,ntarg,j,n,nterms,i,ifder,ii,nd
      real *8 rscale,center(2),ztarg(2,ntarg),zdiff(2),rinv,rinv2
      complex *16 pot,grad,mpole(nd,0:nterms)
      complex *16 pot1(nd,ntarg),grad1(nd,ntarg),hess1(nd,ntarg)
      complex *16 zpow(0:nterms+3)
      complex *16 zpowg(0:nterms+3)
      complex *16 zpowh(0:nterms+3)
c
c
      complex *16 z,z1scale,z2scale,z3scale,z4scale
      complex *16 ztemp1,ztemp2
c
      rinv = 1/rscale
      rinv2 = rinv**2
      zpowg(0) = 0
      zpowh(0) = 0
      zpowh(1) = 0

      do k=1,ntarg
         zdiff(1)=ztarg(1,k)-center(1)
         zdiff(2)=ztarg(2,k)-center(2)
         z = dcmplx(zdiff(1),zdiff(2))
         zpow(0) = log(abs(z))
         ztemp1 = rscale/z
         zpow(1) = ztemp1
         do n=2,nterms+2
           zpow(n) = zpow(n-1)*ztemp1
         enddo

         zpowg(0) = zpow(1)*rinv
         do n=1,nterms
           zpowg(n) = -zpow(n+1)*rinv*n
         enddo

         zpowh(0) = -zpow(2)*rinv2
         do n=1,nterms
           zpowh(n) = zpow(n+2)*n*(n+1)*rinv2
         enddo

         do n=0,nterms
           do ii=1,nd
             pot1(ii,k) = pot1(ii,k) + mpole(ii,n)*zpow(n)
             grad1(ii,k) = grad1(ii,k) + mpole(ii,n)*zpowg(n)
             hess1(ii,k) = hess1(ii,k) + mpole(ii,n)*zpowh(n)
           enddo
         enddo
      enddo
      
      return
      end
c
c
c
c
C***********************************************************************
      subroutine l2dformtac_vec(nd,rscale,source,ns,charge,
     1                center,nterms,local)
      implicit none
C***********************************************************************
c
c     This subroutine INCREMENTS local expansions about CENTER due
c     to (ND,NS) charges located at SOURCES(2,*).
c
c     local_0  = local_0 + sum charge_j  log(abs(z0))
c                           j  
c               
c     local_n  = local_n - sum charge_j 1/n (1/z0)^n rscale^n
c                           j  
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              : vector length (number of mpole expansions)
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     ns              : number of sources
c     charge(nd,ns)   : source strengths
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     local(nd,*)     : coeffs for the local-expansion are incremented
c-----------------------------------------------------------------------
      integer ns,j,ifder,n,nterms,ii,nd
      real *8 center(2),source(2,ns),zdiff(2),rscale
      complex *16 local(nd,0:nterms),charge(nd,ns)
      complex *16 ztemp1,ztemp2,z0,zpow(0:nterms)
      complex *16 z
c
c
      do j=1,ns
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)

        z0 = dcmplx(zdiff(1),zdiff(2))
        ztemp1 = rscale/z0
        zpow(0) = -1
        do n=1,nterms
          zpow(n) = zpow(n-1)*ztemp1
        enddo
        do n=1,nterms
          zpow(n) = zpow(n)/(n+0.0d0)
        enddo
        zpow(0) = log(abs(z0))
        do n=0,nterms
          do ii=1,nd
            local(ii,n) = local(ii,n) + charge(ii,j)*zpow(n)
          enddo
        enddo
      enddo

      return
      end
c
c
C***********************************************************************
      subroutine l2dformtad_vec(nd,rscale,source,ns,
     1           dipstr,center,nterms,local)
      implicit none
C***********************************************************************
c
c     This subroutine INCREMENTS local expansions about CENTER due
c     to (ND,NS) dipoles located at SOURCES(2,*).
c
c     for dipoles
c      
c     local_n  =  -sum dipstr_j rscale^(n)/z0^(n+1)
c                   j  
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              : vector length (number of mpole expansions)
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     ns              : number of sources
c     dipstr(nd,ns)   : dipole strengths
c     center(2)       : expansion center
c     nterms          : order of local expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     local(nd,*)     : coeffs for the local-expansion are incremented
c-----------------------------------------------------------------------
c
      integer ns,j,ifder,n,nterms,ii,nd
      real *8 center(2),source(2,ns),zdiff(2),rscale
      real *8 r,theta
      complex *16 local(nd,0:nterms),dipstr(nd,ns)
      complex *16 zmul,zinv,ztemp1,ztemp2,zpow(0:nterms)
      complex *16 z,z0
c
c
      do j=1,ns
         zdiff(1)=source(1,j)-center(1)
         zdiff(2)=source(2,j)-center(2)
         z0=dcmplx(zdiff(1),zdiff(2))
         ztemp1 = rscale/z0
         zpow(0) = -1/z0
         do n=1,nterms
            zpow(n) = zpow(n-1)*ztemp1
         enddo
         do n=0,nterms
            do ii=1,nd
               local(ii,n) = local(ii,n) + dipstr(ii,j)*zpow(n)
            enddo
         enddo
      enddo

      return
      end
c
c
c
c
C***********************************************************************
      subroutine l2dformtacd_vec(nd,rscale,source,ns,charge,
     1                dipstr,center,nterms,local)
      implicit none
C***********************************************************************
c
c     This subroutine INCREMENTS local expansions about CENTER due
c     to (ND,NS) sources located at SOURCES(2,*).
c
c     local_0  = local_0 + sum charge_j  log(abs(z0))
c                           j
c     local_n  =  local_n -sum (dipstr_j/z0^(n+1) + charge_j/n/z0^n)*rscale^n 
c                           j  
c               
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              : vector length (number of mpole expansions)
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     ns              : number of sources
c     charge(ns)      : source strengths
c     dipstr(ns)      : dipole strengths
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     local           : coeffs for the incoming-expansion are incremented
c-----------------------------------------------------------------------

      integer ns,j,ifder,n,nterms,ii,nd
      real *8 center(2),source(2,ns),zdiff(2),rscale
      complex *16 local(nd,0:nterms),charge(nd,ns),dipstr(nd,ns)
      complex *16 ztemp1,ztemp2,z0,zpow(0:nterms),zpowd(0:nterms)
      complex *16 z
c
c
      do j=1,ns
        zdiff(1)=source(1,j)-center(1)
        zdiff(2)=source(2,j)-center(2)

        z0 = dcmplx(zdiff(1),zdiff(2))
        ztemp1 = rscale/z0
        zpow(0) = -1
        do n=1,nterms
          zpow(n) = zpow(n-1)*ztemp1
        enddo
        do n=1,nterms
          zpow(n) = zpow(n)/(n+0.0d0)

        enddo
        zpow(0) = log(abs(z0))
        
        ztemp1 = rscale/z0
        zpowd(0) = -1/z0
        do n=1,nterms
          zpowd(n) = zpowd(n-1)*ztemp1
        enddo

        do n=0,nterms
          do ii=1,nd
            local(ii,n) = local(ii,n) + charge(ii,j)*zpow(n)+
     1        dipstr(ii,j)*zpowd(n)
          enddo
        enddo
      enddo


      return
      end
c
c
c
c
c**********************************************************************
      subroutine l2dtaevalp_vec(nd,rscale,center,local,nterms,
     1           ztarg,ntarg,pot1)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials 
c     due to incoming local expansions.
c
c                     +nterms
c     POT1  = POT1 +    sum   local_n (z/rscale)^n
c                       n=0  
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of local expansions)
c     rscale       :    scaling parameter 
c     center       :    expansion center
c     local(nd,*)  :    local expansion 
c     nterms       :    order of the multipole expansion
c     ztarg        :    target location
c     ntarg        :    number of targets
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot1(nd,ntarg)  :    potentials at ztarg are incremented
c-----------------------------------------------------------------------
c
      integer j,k,n,ntarg,nterms,ifder,ii,nd,i
      real *8 rscale,r,theta,center(2),ztarg(2,*),zdiff(2)
      complex *16 local(nd,0:nterms)
      complex *16 pot1(nd,ntarg),zpow(0:nterms)
c
c
      complex *16 zmul,zinv,ztemp1,ztemp2,z
c
c

      do k=1,ntarg
         zdiff(1)=ztarg(1,k)-center(1)
         zdiff(2)=ztarg(2,k)-center(2)
         z = dcmplx(zdiff(1),zdiff(2))/rscale
         zpow(0) = 1
         do i=1,nterms
           zpow(i) = zpow(i-1)*z
         enddo
        
         do i=0,nterms
           do ii=1,nd
             pot1(ii,k) = pot1(ii,k) + local(ii,i)*zpow(i)
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
c
c**********************************************************************
      subroutine l2dtaevalg_vec(nd,rscale,center,local,nterms,
     1           ztarg,ntarg,pot1,grad1)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials and gradients 
c     due to incoming local expansions.
c
c                     +nterms
c     POT1  = POT1 +    sum   local_n (z/rscale)^n
c                       n=0  
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of local expansions)
c     rscale       :    scaling parameter 
c     center       :    expansion center
c     local(nd,*)  :    local expansion 
c     nterms       :    order of the multipole expansion
c     ztarg        :    target location
c     ntarg        :    number of targets
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot1(nd,ntarg)   :    potentials at ztarg are incremented
c     grad1(nd,ntarg)  :    gradients at ztarg are incremented
c-----------------------------------------------------------------------
      integer i,j,k,n,nterms,ntarg,ifder,ii,nd,idim
      real *8 rscale,r,theta,center(2),ztarg(2,*),zdiff(2),rinv
      complex *16 local(nd,0:nterms)
      complex *16 pot1(nd,ntarg),grad1(nd,ntarg)
      complex *16 zpow(0:nterms),zpowg(0:nterms)
c
c
      complex *16 zmul,zinv,ztemp1,ztemp2,z
c
c
c
      rinv = 1/rscale
      do k=1,ntarg
         zdiff(1)=ztarg(1,k)-center(1)
         zdiff(2)=ztarg(2,k)-center(2)
         z = dcmplx(zdiff(1),zdiff(2))/rscale
         zpow(0) = 1
         zpowg(0) = 0
         do n=1,nterms
           zpow(n) = zpow(n-1)*z
         enddo
         do n=1,nterms
           zpowg(n) = zpow(n-1)*n*rinv
         enddo

         do n=0,nterms
           do ii=1,nd
             pot1(ii,k) = pot1(ii,k) + local(ii,n)*zpow(n)
             grad1(ii,k) = grad1(ii,k) + local(ii,n)*zpowg(n)
           enddo
         enddo
c
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
      subroutine l2dtaevalh_vec(nd,rscale,center,local,nterms,
     1           ztarg,ntarg,pot1,grad1,hess1)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials, gradients, and hessians 
c     due to incoming local expansions.
c
c                     +nterms
c     POT1  = POT1 +    sum   local_n (z/rscale)^n
c                       n=0  
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of local expansions)
c     rscale       :    scaling parameter 
c     center       :    expansion center
c     local(nd,*)  :    local expansion 
c     nterms       :    order of the multipole expansion
c     ztarg        :    target location
c     ntarg        :    number of targets
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     pot1(nd,ntarg)   :    potentials at ztarg are incremented
c     grad1(nd,ntarg)  :    gradients at ztarg are incremented
c     hess1(nd,ntarg)  :    hessians at ztarg are incremented
c-----------------------------------------------------------------------
      integer i,j,k,n,nterms,ntarg,ifder,ii,nd
      real *8 rscale,r,theta,center(2),ztarg(2,*),zdiff(2),rinv,rinv2
      complex *16 local(nd,0:nterms)
      complex *16 pot1(nd,ntarg),grad1(nd,ntarg),hess1(nd,ntarg)
      complex *16 zpow(0:nterms),zpowg(0:nterms),zpowh(0:nterms)
c
c
      complex *16 zmul,zinv,ztemp1,ztemp2,z
c
c
c
      rinv = 1/rscale
      rinv2 = rinv**2
      do k=1,ntarg
         zdiff(1)=ztarg(1,k)-center(1)
         zdiff(2)=ztarg(2,k)-center(2)
         z = dcmplx(zdiff(1),zdiff(2))/rscale
         zpow(0) = 1
         zpowg(0) = 0
         zpowh(0) = 0
         zpowh(1) = 0
         do n=1,nterms
           zpow(n) = zpow(n-1)*z
         enddo
         do n=1,nterms
           zpowg(n) = zpow(n-1)*n*rinv
         enddo
         do n=2,nterms
           zpowh(n) = zpow(n-2)*n*(n-1)*rinv2
         enddo

         do n=0,nterms
           do ii=1,nd
             pot1(ii,k) = pot1(ii,k) + local(ii,n)*zpow(n)
             grad1(ii,k) = grad1(ii,k) + local(ii,n)*zpowg(n)
             hess1(ii,k) = hess1(ii,k) + local(ii,n)*zpowh(n)
           enddo
         enddo
c
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
      subroutine l2dmpmp_vec(nd,rscale1,center1,hexp1,nterms1,
     $                      rscale2,center2,hexp2,nterms2,carray,ldc)
      implicit none
C**********************************************************************
C
C     This routine shifts multipole expansions HEXP1 to a new center and
C     INCREMENTS the multipole expansions HEXP2 about that point 
C     accordingly.
C
C---------------------------------------------------------------------
C     INPUT:
C
c     nd     :    vector length (number of mpole expansions)
C     rscale1:   scaling parameter for original multipole expansion
C     center1:   center of original multiple expansion
C     hexp1  :   coefficients of original multiple expansions
C     nterms1:   order of original multipole expansion
C     rscale2:   scaling parameter for shifted multipole expansion
C     center2:   center of shifted multipole expansion
C     nterms2:   order of shifted multipole expansion
C     carray :   array of binomial coefficients
C     ldc    :   leading dimension of carray array
C---------------------------------------------------------------------
C     OUTPUT:
C
C     hexp2  :  coefficients of shifted multipole expansions
C---------------------------------------------------------------------
      integer nterms1,nterms2,nterms,i,j,ifder,ii,nd,ldc,nmax
      real *8 rscale1,rscale2,center1(2),center2(2),zdiff(2)
      real *8 carray(0:ldc,0:ldc)
      complex *16 hexp1(nd,0:nterms1),hexp2(nd,0:nterms2)
      complex *16, allocatable :: z0pow1(:),z0pow2(:),hexp1tmp(:,:)
      complex *16, allocatable :: hexp2tmp(:,:)
      real *8 rtmp
      complex *16 z,z0,zmul,zinv,ztemp1,ztemp2
c
c
      nterms = nterms1+nterms2
      nmax = max(nterms1,nterms2)
c
c
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      z0 = -dcmplx(zdiff(1),zdiff(2))

      allocate(z0pow1(0:nmax))
      ztemp1 = 1/(z0/rscale1)
      ztemp2 = ztemp1
      z0pow1(0) = 1
      do i=1,nmax
        z0pow1(i) = ztemp1
        ztemp1 = ztemp1*ztemp2
      enddo

      allocate(z0pow2(0:nmax))
      ztemp1 = z0/rscale2
      ztemp2 = ztemp1
      z0pow2(0) = 1
      do i=1,nmax
        z0pow2(i) = ztemp1
        ztemp1 = ztemp1*ztemp2
      enddo


      allocate(hexp1tmp(nd,0:nterms1),hexp2tmp(nd,0:nterms2))
      do i=0,nterms1
        do ii=1,nd
          hexp1tmp(ii,i) = 0
        enddo
      enddo

      do i=0,nterms2
        do ii=1,nd
          hexp2tmp(ii,i) = 0
        enddo
      enddo

      do i=0,nterms1
        do ii=1,nd
          hexp1tmp(ii,i) = hexp1(ii,i)*z0pow1(i)
        enddo
      enddo

      do ii=1,nd
        hexp2tmp(ii,0) = hexp2tmp(ii,0) + hexp1(ii,0)
      enddo

      do i=1,nterms2
        do ii=1,nd
          hexp2tmp(ii,i) = hexp2tmp(ii,i)-hexp1tmp(ii,0)/i
        enddo
        do j=1,min(i,nterms1)
          do ii=1,nd
            hexp2tmp(ii,i) = hexp2tmp(ii,i) +
     1          hexp1tmp(ii,j)*carray(i-1,j-1)
          enddo
        enddo

        do ii=1,nd
          hexp2tmp(ii,i) = hexp2tmp(ii,i)*z0pow2(i)
        enddo
      enddo

      do i=0,nterms2
        do ii=1,nd
          hexp2(ii,i) = hexp2(ii,i) + hexp2tmp(ii,i)
        enddo
      enddo

      return
      end
c
c
c
c**********************************************************************
      subroutine l2dlocloc_vec(nd,rscale1,center1,jexp1,nterms1,
     $                      rscale2,center2,jexp2,nterms2,carray,ldc)
      implicit none
C**********************************************************************
C
C     This routine shifts local expansions jEXP1 to a new center and
C     INCREMENTS the local expansions jEXP2 about that point 
C     accordingly.
C
C---------------------------------------------------------------------
C     INPUT:
C
c     nd     :    vector length (number of mpole expansions)
C     rscale1:   scaling parameter for original multipole expansion
C     center1:   center of original multiple expansion
C     jexp1  :   coefficients of original local expansions
C     nterms1:   order of original multipole expansion
C     rscale2:   scaling parameter for shifted multipole expansion
C     center2:   center of shifted multipole expansion
C     nterms2:   order of shifted multipole expansion
C     carray :   array of binomial coefficients
C     ldc    :   leading dimension of carray array
C---------------------------------------------------------------------
C     OUTPUT:
C
C     jexp2  :  coefficients of shifted local expansions
C---------------------------------------------------------------------
      integer nterms1,nterms2,nterms,i,j,ifder,ii,nd,ldc,nmax
      real *8 rscale1,rscale2,center1(2),center2(2),zdiff(2)
      real *8 carray(0:ldc,0:ldc)
      complex *16 jexp1(nd,0:nterms1),jexp2(nd,0:nterms2)
      complex *16, allocatable :: z0pow1(:),z0pow2(:),jexp1tmp(:,:)
      complex *16, allocatable :: jexp2tmp(:,:)
      real *8 rtmp
      complex *16 z,z0,zmul,zinv,ztemp1,ztemp2
c
c
      nterms = nterms1+nterms2
      nmax = max(nterms1,nterms2)
c
c
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      z0 = dcmplx(zdiff(1),zdiff(2))

      allocate(z0pow1(0:nmax))
      ztemp1 = (z0/rscale1)
      ztemp2 = ztemp1
      z0pow1(0) = 1
      do i=1,nmax
        z0pow1(i) = ztemp1
        ztemp1 = ztemp1*ztemp2
      enddo

      allocate(z0pow2(0:nmax))
      ztemp1 = 1/(z0/rscale2)
      ztemp2 = ztemp1
      z0pow2(0) = 1
      do i=1,nmax
        z0pow2(i) = ztemp1
        ztemp1 = ztemp1*ztemp2
      enddo
      

      allocate(jexp1tmp(nd,0:nterms1),jexp2tmp(nd,0:nterms2))
      do i=0,nterms1
        do ii=1,nd
          jexp1tmp(ii,i) = 0
        enddo
      enddo

      do i=0,nterms2
        do ii=1,nd
          jexp2tmp(ii,i) = 0
        enddo
      enddo

      do i=0,nterms1
        do ii=1,nd
          jexp1tmp(ii,i) = jexp1(ii,i)*z0pow1(i)
        enddo
      enddo

      do i=0,nterms2
        do j=i,nterms1
          do ii=1,nd
            jexp2tmp(ii,i) = jexp2tmp(ii,i) +
     1          jexp1tmp(ii,j)*carray(j,i)
          enddo
        enddo

        do ii=1,nd
          jexp2tmp(ii,i) = jexp2tmp(ii,i)*z0pow2(i)
        enddo
      enddo

      do i=0,nterms2
        do ii=1,nd
          jexp2(ii,i) = jexp2(ii,i) + jexp2tmp(ii,i)
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
      subroutine l2dmploc_vec(nd,rscale1,center1,hexp1,nterms1,
     $                      rscale2,center2,jexp2,nterms2,carray,ldc)
      implicit none
C**********************************************************************
C
C     This routine shifts local expansions HEXP1 to a new center and
C     INCREMENTS the local expansions JEXP2 about that point 
C     accordingly.
C
C---------------------------------------------------------------------
C     INPUT:
C
c     nd     :    vector length (number of mpole expansions)
C     rscale1:   scaling parameter for original multipole expansion
C     center1:   center of original multiple expansion
C     hexp1  :   coefficients of original local expansions
C     nterms1:   order of original multipole expansion
C     rscale2:   scaling parameter for shifted multipole expansion
C     center2:   center of shifted multipole expansion
C     nterms2:   order of shifted multipole expansion
C     carray :   array of binomial coefficients
C     ldc    :   leading dimension of carray array
C---------------------------------------------------------------------
C     OUTPUT:
C
C     jexp2  :  coefficients of shifted local expansions
C---------------------------------------------------------------------
      integer nterms1,nterms2,nterms,i,j,ifder,ii,nd,ldc,nmax
      real *8 rscale1,rscale2,center1(2),center2(2),zdiff(2)
      real *8 carray(0:ldc,0:ldc)
      complex *16 hexp1(nd,0:nterms1),jexp2(nd,0:nterms2)
      complex *16, allocatable :: z0pow1(:),z0pow2(:),hexp1tmp(:,:)
      complex *16, allocatable :: jexp2tmp(:,:)
      real *8 rtmp
      complex *16 z,z0,zmul,zinv,ztemp1,ztemp2,ztemp3
c
c
      nterms = nterms1+nterms2
      nmax = max(nterms1,nterms2)
c
c
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      z0 = -dcmplx(zdiff(1),zdiff(2))

      allocate(z0pow1(0:nmax))
      allocate(z0pow2(0:nmax))
      ztemp1 = 1/z0
      z0pow1(0) = 1
      z0pow2(0) = 1
      ztemp2 = ztemp1*rscale2
      ztemp3 = -ztemp1*rscale1
      z0pow1(0) = 1
      do i=1,nmax
        z0pow2(i) = ztemp2
        z0pow1(i) = ztemp3
        ztemp2 = ztemp2*ztemp1*rscale2
        ztemp3 = -ztemp3*ztemp1*rscale1
      enddo


      allocate(hexp1tmp(nd,0:nterms1),jexp2tmp(nd,0:nterms2))
      do i=0,nterms1
        do ii=1,nd
          hexp1tmp(ii,i) = 0
        enddo
      enddo

      do i=0,nterms2
        do ii=1,nd
          jexp2tmp(ii,i) = 0
        enddo
      enddo

      do i=0,nterms1
        do ii=1,nd
          hexp1tmp(ii,i) = hexp1(ii,i)*z0pow1(i)
        enddo
      enddo

      rtmp = log(abs(z0))
      do ii=1,nd
        jexp2tmp(ii,0) = hexp1tmp(ii,0)*rtmp
      enddo

      do j=1,nterms1
        do ii=1,nd
          jexp2tmp(ii,0) = jexp2tmp(ii,0)+hexp1tmp(ii,j)
        enddo
      enddo

      do i=1,nterms2
        do ii=1,nd
          jexp2tmp(ii,i) = jexp2tmp(ii,i) - hexp1tmp(ii,0)/i
        enddo
        do j=1,nterms1
          do ii=1,nd
            jexp2tmp(ii,i) = jexp2tmp(ii,i) +
     1          hexp1tmp(ii,j)*carray(i+j-1,j-1)
          enddo
        enddo

        do ii=1,nd
          jexp2tmp(ii,i) = jexp2tmp(ii,i)*z0pow2(i)
        enddo
      enddo

      do i=0,nterms2
        do ii=1,nd
          jexp2(ii,i) = jexp2(ii,i) + jexp2tmp(ii,i)
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
C***********************************************************************
      subroutine l2dmpzero_vec(nd,mpole,nterms)
      implicit none
C***********************************************************************
c
c     This subroutine sets a vector multipole expansion to zero.
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd     :   vector length (number of mpole expansions)
c     nterms :   order of multipole expansion
C---------------------------------------------------------------------
c     OUTPUT:
c
c     mpole  :   coeffs for the expansion set to zero.
C---------------------------------------------------------------------
      integer n,nterms,nd,ii
      complex *16 mpole(nd,0:nterms)
c
      do n=0,nterms
      do ii=1,nd
         mpole(ii,n)=0.0d0
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
        subroutine l2d_init_carray(carray,ldc)
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
        
