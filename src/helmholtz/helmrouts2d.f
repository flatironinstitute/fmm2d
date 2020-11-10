cc Copyright (C) 2010-2011: Leslie Greengard and Zydrunas Gimbutas
cc Contact: greengard@cims.nyu.edu
c
c  Modified: (C) 2019, Leslie Greengard
cc 
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
c      This file contains the basic subroutines for 
c      forming and evaluating multipole (partial wave) expansions
c      in two dimensions.
c
c      Remarks on scaling conventions.
c
c      1)  Hankel and Bessel functions are consistently scaled as
c       	hvec(n)= H_n(z)*rscale^(n)
c       	jvec(n)= J_n(z)/rscale^(n)
c
c          rscale should be of the order of |z| if |z| < 1. Otherwise,
c          rscale should be set to 1.
c
c      2) The potential scaling is that required of the delta function
c      response:  pot = H_0(k*r)*(eye/4)
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
c      H2D_DIRECT or H2D_DIRECT1 return the result of the requested
c      computation by overwriting the input pot/grad/hess variables.
c-----------------------------------------------------------------------
c
c      H2DFORMMPC_VEC:  creates multipole expansions (outgoing) due to 
c                       a collection of charge sources.
c      H2DFORMMPD_VEC:  creates multipole expansions (outgoing) due to 
c                       a collection of dipole sources.
c      H2DFORMMPCD_VEC: creates multipole expansions (outgoing) due to 
c                       a collection of charge and dipole sources.
c      H2DMPEVALP_VEC:  computes potentials from multipole expansions
c                       at a collection of targets
c      H2DMPEVALG_VEC:  computes potentials/gradients from multipole expansions
c                       at a collection of targets
c      H2DMPEVALH_VEC:  computes potentials/gradients/Hessians 
c                       from multipole expansions at a collection of targets
c
c      H2DFORMTAC_VEC:  creates local expansions due to 
c                       a collection of charge sources
c      H2DFORMTAD_VEC:  creates local expansions due to 
c                       a collection of dipole sources
c      H2DFORMTACD_VEC: creates local expansions due to 
c                       a collection of charge and dipole sources
c      H2DTAEVALP_VEC:  computes potentials from local expansions
c                       at a collection of targets
c      H2DTAEVALG_VEC:  computes potentials/gradients from local expansions
c                       at a collection of targets
c      H2DTAEVALH_VEC:  computes potentials/gradients/Hessians 
c                       from local expansions at a collection of targets
c 
c      H2DMPMP_VEC:     Translates center of multipole expansions 
c      H2DLOCLOC_VEC:   Translates center of local expansions 
c      H2DMPLOC_VEC:    Converts multipole expansions to local expansions
c
c      H2DMPZERO_VEC:   utility to initialize multipole coefficients to zero.
C
c-----------------------------------------------------------------------
c
c
C***********************************************************************
      subroutine h2dformmpc_vec(nd,zk,rscale,source,ns,charge,
     1                      center,nterms,mpole)
      implicit none
C***********************************************************************
c
c     This subroutine INCREMENTS multipole (outgoing) expansions about 
c     CENTER due to (ND,NS) charges located at SOURCES(2,NS).
c
c     mpole_n(q)=mpole_n(q)+
c
c                sum charge_j(q) J_n(k r)exp(-i n theta_j)/rscale_j^n
c                 j  
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of mpole expansions)
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     ns              : number of sources
c     charge(nd,ns)      : source strengths
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mpole(nd,*)     : coeffs for the h-expansion are incremented.
c-----------------------------------------------------------------------
      integer ns,j,ifder,n,nterms,ii,nd
      real *8 center(2),source(2,ns),zdiff(2),rscale
      real *8 r,theta
      complex *16 zk,mpole(nd,-nterms:nterms),charge(nd,ns)
      complex *16 zmul,zinv,ztemp1,ztemp2
      complex *16 ima,z
      complex *16, allocatable :: jval(:)
      complex *16, allocatable :: jder(:)
c
      data ima/(0.0d0,1.0d0)/
c
      allocate(jval(0:nterms+5))
      allocate(jder(0:nterms+5))
c
      do j=1,ns
         zdiff(1)=source(1,j)-center(1)
         zdiff(2)=source(2,j)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call jbessel2d(nterms+1,z,rscale,jval,ifder,jder)
         zmul=exp(-ima*theta)
         zinv=conjg(zmul)
         call ctompole_vec(nd,zmul,zinv,mpole,jval,charge(1,j),nterms) 
      enddo
      return
      end
c
c
c
c
C***********************************************************************
      subroutine h2dformmpd_vec(nd,zk,rscale,source,ns,
     1                      dipstr,dipvec,center,nterms,mpole)
      implicit none
C***********************************************************************
c
c     This subroutine INCREMENTS multipole (outgoing) expansions about 
c     CENTER due to NS dipoles located at SOURCES(2,*).
c
c     For charges,               
c     mpole_n = mpole_n +  sum charge_j  J_n(k r) exp(-i n theta_j) /rscale_j^n
c                           j  
c     For dipoles, the formula is more involved.              
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of mpole expansions)
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     ns              : number of sources
c     dipstr(nd,ns)   : dipole strengths
c     dipvec(nd,2,ns) : dipole orientation vectors
c
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mpole(nd,*)     : coeffs for the h-expansion are incremented.
c-----------------------------------------------------------------------
      integer ns,j,ifder,n,nterms,ii,nd
      real *8 center(2),source(2,ns),zdiff(2),dipvec(nd,2,ns),rscale
      real *8 r,theta
      complex *16 zk,mpole(nd,-nterms:nterms),dipstr(nd,ns)
      complex *16 zmul,zinv,ztemp1,ztemp2
      complex *16 ima,ima4,z
      complex *16, allocatable :: jval(:)
      complex *16, allocatable :: jder(:)
c
      data ima/(0.0d0,1.0d0)/
c
      allocate(jval(0:nterms+5))
      allocate(jder(0:nterms+5))
c
      do j=1,ns
         zdiff(1)=source(1,j)-center(1)
         zdiff(2)=source(2,j)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call jbessel2d(nterms+1,z,rscale,jval,ifder,jder)
         zmul=exp(-ima*theta)
         zinv=conjg(zmul)
         call dtompole_vec(nd,zk,rscale,zmul,zinv,
     1           mpole,jval,dipstr(1,j),dipvec(1,1,j),nterms) 
      enddo
      return
      end
c
c
c
c
C***********************************************************************
      subroutine h2dformmpcd_vec(nd,zk,rscale,source,ns,charge,
     1                       dipstr,dipvec,center,nterms,mpole)
      implicit none
C***********************************************************************
c
c     This subroutine INCREMENTS multipole (outgoing) expansions about 
c     CENTER due to NS sources located at SOURCES(2,*).
c
c     For charges,               
c     mpole_n = mpole_n +  sum charge_j  J_n(k r) exp(-i n theta_j) /rscale_j^n
c                           j  
c     For dipoles, the formula is more involved.              
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of mpole expansions)
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     ns              : number of sources
c     charge(nd,ns)   : source strengths
c     dipstr(nd,ns)   : dipole strengths
c     dipvec(nd,2,ns) : dipole orientation vectors
c
c     center(2)       : expansion center
c     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mpole(nd,*)     : coeffs for the h-expansion are incremented.
c-----------------------------------------------------------------------
      integer ns,j,ifder,n,nterms,ii,nd
      real *8 center(2),source(2,ns),zdiff(2),dipvec(nd,2,ns),rscale
      real *8 r,theta
      complex *16 mpole(nd,-nterms:nterms),charge(nd,ns),dipstr(nd,ns)
      complex *16 zmul,zinv,ztemp1,ztemp2
      complex *16 zk,ima,ima4,z
      complex *16, allocatable :: jval(:)
      complex *16, allocatable :: jder(:)
c
      data ima/(0.0d0,1.0d0)/
c
      allocate(jval(0:nterms+5))
      allocate(jder(0:nterms+5))
c
      do j=1,ns
         zdiff(1)=source(1,j)-center(1)
         zdiff(2)=source(2,j)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call jbessel2d(nterms+1,z,rscale,jval,ifder,jder)
         zmul=exp(-ima*theta)
         zinv=conjg(zmul)
c
         call ctompole_vec(nd,zmul,zinv,mpole,jval,charge(1,j),nterms) 
         call dtompole_vec(nd,zk,rscale,zmul,zinv,
     1           mpole,jval,dipstr(1,j),dipvec(1,1,j),nterms) 
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine h2dmpevalp_vec(nd,zk,rscale,center,mpole,nterms,
     1                     ztarg,ntarg,pot1)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials due to
c     outgoing partial wave expansions.
c                    +nterms
c     pot1  = pot1 +     sum   mpole_n H_n(k r) exp(i n theta)
c                  n=-nterms  
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
      complex *16 zk,pot,mpole(nd,-nterms:nterms)
      complex *16 pot1(nd,ntarg)
c
      complex *16, allocatable :: hval(:)
      complex *16, allocatable :: hder(:)
      complex *16, allocatable :: mptemp(:)
c
      complex *16 ima,ima4inv,z
      complex *16 zmul,zinv,ztemp1,ztemp2
      data ima/(0.0d0,1.0d0)/
c
c
      ima4inv=ima/4
c
      allocate(hval(0:nterms+5))
      allocate(hder(0:nterms+5))
      allocate(mptemp(-nterms-2:nterms+2))
c
      do k=1,ntarg
         zdiff(1)=ztarg(1,k)-center(1)
         zdiff(2)=ztarg(2,k)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call h2dall(nterms+3,z,rscale,hval,ifder,hder)
c
         zmul=exp(ima*theta)
         zinv=conjg(zmul)
         call mpole_evalp_vec(nd,zmul,zinv,mpole,mptemp,hval,
     $        nterms,pot1(1,k)) 

      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine h2dmpevalg_vec(nd,zk,rscale,center,mpole,nterms,
     1                     ztarg,ntarg,pot1,grad1)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials and gradients of the 
c     potentials due to outgoing partial wave expansions.
c                    +nterms
c     pot1  = pot1 +     sum   mpole_n H_n(k r) exp(i n theta)
c                  n=-nterms  
c
c     grad = gradient(pot1) 
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
c     pot1(nd,*)    :   potentials at ztarg locations are incremented
c     grad1(nd,2,*) :   gradients at ztarg locations are incremented
c-----------------------------------------------------------------------
      integer k,ntarg,j,n,nterms,i,ifder,ii,nd
      real *8 rscale,r,theta,center(2),ztarg(2,ntarg),zdiff(2)
      complex *16 zk,pot,grad(2),mpole(nd,-nterms:nterms)
      complex *16 pot1(nd,ntarg),grad1(nd,2,ntarg)
c
      complex *16, allocatable :: hval(:)
      complex *16, allocatable :: hder(:)
c
      complex *16, allocatable :: mpolex(:,:)
      complex *16, allocatable :: mpoley(:,:)
c
      complex *16, allocatable :: mptemp(:)
c
      complex *16 ima,ima4,ima4inv,z,z1scale,z2scale,z3scale,z4scale
      complex *16 zmul,zinv,ztemp1,ztemp2
      data ima/(0.0d0,1.0d0)/
c
c
      ima4=-4*ima
      ima4inv=ima/4
c
      allocate(hval(0:nterms+5))
      allocate(hder(0:nterms+5))
c
      allocate(mpolex(nd,-nterms-1:nterms+1))
      allocate(mpoley(nd,-nterms-1:nterms+1))
      allocate(mptemp(-nterms-2:nterms+2))
c
      call mk_mpoleg_vec(nd,zk,rscale,mpole,mpolex,mpoley,nterms)
c
      do k=1,ntarg
         zdiff(1)=ztarg(1,k)-center(1)
         zdiff(2)=ztarg(2,k)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call h2dall(nterms+3,z,rscale,hval,ifder,hder)
c
         zmul=exp(ima*theta)
         zinv=conjg(zmul)
         call mpole_evalp_vec(nd,zmul,zinv,mpole,mptemp,hval,
     $        nterms,pot1(1,k)) 
         call mpole_evalg_vec(nd,mpolex,mpoley,mptemp,nterms,
     $        grad1(1,1,k)) 
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine h2dmpevalh_vec(nd,zk,rscale,center,mpole,nterms,
     1                     ztarg,ntarg,pot1,grad1,hess1)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials and derivatives of the 
c     potential due to outgoing partial wave expansions.
c                      +nterms
c     pot1  = pot1 +     sum   mpole_n H_n(k r) exp(i n theta)
c                     n=-nterms  
c
c     grad = gradient(pot1) 
c     hess = hessian  (pxx,pxy,pyy)
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
c     pot1(nd,ntarg)      :   potentials at ztarg locations are incremented
c     grad1(nd,2,ntarg)   :   gradients at ztarg locations are incremented
c     hess1(nd,3,ntarg)   :   hessians at ztarg locations are incremented
c-----------------------------------------------------------------------
      integer k,ntarg,j,n,nterms,i,ifder,ii,nd
      real *8 rscale,r,theta,center(2),ztarg(2,ntarg),zdiff(2)
      complex *16 zk,pot,grad(2),hess(3),mpole(nd,-nterms:nterms)
      complex *16 pot1(nd,ntarg),grad1(nd,2,ntarg),hess1(nd,3,ntarg)
c
      complex *16, allocatable :: hval(:)
      complex *16, allocatable :: hder(:)
c
      complex *16, allocatable :: mpolex(:,:)
      complex *16, allocatable :: mpoley(:,:)
c
      complex *16, allocatable :: mpolexx(:,:)
      complex *16, allocatable :: mpolexy(:,:)
      complex *16, allocatable :: mpoleyy(:,:)
c
      complex *16, allocatable :: mptemp(:)
c
      complex *16 ima,ima4,ima4inv,z,z1scale,z2scale,z3scale,z4scale
      complex *16 zmul,zinv,ztemp1,ztemp2
      data ima/(0.0d0,1.0d0)/
c
c
      ima4=-4*ima
      ima4inv=ima/4
c
      allocate(hval(0:nterms+5))
      allocate(hder(0:nterms+5))
      allocate(mpolex(nd,-nterms-1:nterms+1))
      allocate(mpoley(nd,-nterms-1:nterms+1))
      allocate(mpolexx(nd,-nterms-2:nterms+2))
      allocate(mpolexy(nd,-nterms-2:nterms+2))
      allocate(mpoleyy(nd,-nterms-2:nterms+2))
      allocate(mptemp(-nterms-2:nterms+2))
c
      call mk_mpoleg_vec(nd,zk,rscale,mpole,mpolex,mpoley,nterms)
      call mk_mpoleh_vec(nd,zk,rscale,mpolex,mpoley,
     $     mpolexx,mpolexy,mpoleyy,nterms)

c
      do k=1,ntarg
         zdiff(1)=ztarg(1,k)-center(1)
         zdiff(2)=ztarg(2,k)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call h2dall(nterms+3,z,rscale,hval,ifder,hder)
c
         zmul=exp(ima*theta)
         zinv=conjg(zmul)
         call mpole_evalp_vec(nd,zmul,zinv,mpole,mptemp,hval,
     $        nterms,pot1(1,k)) 
         call mpole_evalg_vec(nd,mpolex,mpoley,mptemp,nterms,
     $        grad1(1,1,k)) 
         call mpole_evalh_vec(nd,mpolexx,mpolexy,mpoleyy,mptemp,nterms,
     $        hess1(1,1,k)) 
      enddo
      return
      end
c
c
c
c
C***********************************************************************
      subroutine h2dformtac_vec(nd,zk,rscale,source,ns,charge,
     1                center,nterms,local)
      implicit none
C***********************************************************************
c
c     This subroutine INCREMENTS local (j) expansions about CENTER due
c     to (ND,NS) charges located at SOURCES(2,*).
c
c               
c     local_n  = local_n + sum charge_j  H_n(k r) exp(-i n theta_j) *rscale_j^n
c                           j  
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              :    vector length (number of mpole expansions)
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     ns              : number of sources
c     charge(nd,ns)      : source strengths
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     local(nd,*)     : coeffs for the j-expansion are incremented
c-----------------------------------------------------------------------
      integer j,n,ns,nterms,ifder,ii,nd
      real *8 rscale,r,theta,center(2),source(2,ns),zdiff(2)
      complex *16 zk,local(nd,-nterms:nterms),charge(nd,ns)
c
      complex *16, allocatable :: hval(:)
      complex *16, allocatable :: hder(:)
      complex *16 zmul,zinv,ztemp1,ztemp2,ztemp3,ztemp4,ztemp5,ztemp6
c
      complex *16 ima,ima4,z
      data ima/(0.0d0,1.0d0)/
c
c
      allocate(hval(0:nterms+5))
      allocate(hder(0:nterms+5))
c
      do j=1,ns
         zdiff(1)=source(1,j)-center(1)
         zdiff(2)=source(2,j)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call h2dall(nterms+2,z,rscale,hval,ifder,hder)
         zmul=exp(-ima*theta)
         zinv=conjg(zmul)
c
         call ctompole_vec(nd,zmul,zinv,local,hval,charge(1,j),nterms) 
      enddo
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine h2dformtad_vec(nd,zk,rscale,source,ns,
     1           dipstr,dipvec,center,nterms,local)
      implicit none
C***********************************************************************
c
c     This subroutine INCREMENTS local (j) expansions about CENTER due
c     to (ND,NS) dipoles located at SOURCES(2,*).
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              :    vector length (number of mpole expansions)
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     ns              : number of sources
c     dipstr(nd,ns)   : dipole strengths
c     dipvec(nd,2,ns) : dipole vectors
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     local(nd,*)     : coeffs for the j-expansion are incremented
c-----------------------------------------------------------------------
      integer j,n,ns,nterms,ifder,ii,nd
      real *8 rscale,r,theta,rsinv
      real *8 center(2),source(2,ns),zdiff(2),dipvec(nd,2,ns)
      complex *16 zk,local(nd,-nterms:nterms),dipstr(nd,ns)
c
      complex *16, allocatable :: hval(:)
      complex *16, allocatable :: hder(:)
      complex *16 zmul,zinv,ztemp1,ztemp2,ztemp3,ztemp4,ztemp5,ztemp6
c
      complex *16 ima,ima4,z
      data ima/(0.0d0,1.0d0)/
c
c
      allocate(hval(0:nterms+5))
      allocate(hder(0:nterms+5))
c
      rsinv = 1.0d0/rscale
      do j=1,ns
         zdiff(1)=source(1,j)-center(1)
         zdiff(2)=source(2,j)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call h2dall(nterms+2,z,rscale,hval,ifder,hder)
         zmul=exp(-ima*theta)
         zinv=conjg(zmul)
c
         call dtompole_vec(nd,zk,rsinv,zmul,zinv,
     1           local,hval,dipstr(1,j),dipvec(1,1,j),nterms) 
      enddo
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine h2dformtacd_vec(nd,zk,rscale,source,ns,charge,
     1                dipstr,dipvec,center,nterms,local)
      implicit none
C***********************************************************************
c
c     This subroutine INCREMENTS local (j) expansions about CENTER due
c     to (ND,NS) sources located at SOURCES(2,*).
c
c               
c     local_n  = local_n + sum charge_j  H_n(k r) exp(-i n theta_j) *rscale_j^n
c                           j  
c
c                  + more complicated formula for dipole sources
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              :    vector length (number of mpole expansions)
c     zk              : Helmholtz parameter 
c     rscale          : the scaling factor.
c     source(2,ns)    : coordinates of sources
c     ns              : number of sources
c     charge(ns)      : source strengths
c     dipstr(ns)      : dipole strengths
c     dipvec(nd,2,ns) : dipole vectors
c     center(2)       : epxansion center
c     nterms          : order of local expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     local           : coeffs for the j-expansion are incremented
c-----------------------------------------------------------------------
      integer j,n,ns,nterms,ifder,ii,nd
      real *8 rscale,r,theta,rsinv
      real *8 center(2),source(2,ns),zdiff(2),dipvec(nd,2,ns)
      complex *16 local(nd,-nterms:nterms),charge(nd,ns),dipstr(nd,ns)
      complex *16 zk,zmul,zinv,ztemp1,ztemp2,ztemp3,ztemp4,ztemp5,ztemp6
c
      complex *16, allocatable :: hval(:)
      complex *16, allocatable :: hder(:)
c
      complex *16 ima,ima4,z
      data ima/(0.0d0,1.0d0)/
c
c
      allocate(hval(0:nterms+5))
      allocate(hder(0:nterms+5))
c
      rsinv = 1.0d0/rscale
      do j=1,ns
         zdiff(1)=source(1,j)-center(1)
         zdiff(2)=source(2,j)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call h2dall(nterms+2,z,rscale,hval,ifder,hder)
         zmul=exp(-ima*theta)
         zinv=conjg(zmul)
c
         call ctompole_vec(nd,zmul,zinv,local,hval,charge(1,j),nterms) 
c
         call dtompole_vec(nd,zk,rsinv,zmul,zinv,
     1           local,hval,dipstr(1,j),dipvec(1,1,j),nterms) 
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine h2dtaevalp_vec(nd,zk,rscale,center,local,nterms,
     1           ztarg,ntarg,pot1)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials 
c     due to incoming partial wave expansions.
c
c                     +nterms
c     POT1  = POT1 +    sum   local_n J_n(k r) exp(i n theta)
c                     n=-nterms  
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of local expansions)
c     zk           :    Helmholtz parameter
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
      integer j,k,n,ntarg,nterms,ifder,ii,nd
      real *8 rscale,r,theta,center(2),ztarg(2,*),zdiff(2)
      complex *16 zk,local(nd,-nterms:nterms)
      complex *16 pot1(nd,ntarg)
c
      complex *16, allocatable :: jval(:)
      complex *16, allocatable :: jder(:)
      complex *16, allocatable :: mptemp(:)
c
      complex *16 ima,ima4,ima4inv,z,z1scale,z2scale,z3scale,z4scale
      complex *16 zmul,zinv,ztemp1,ztemp2
      data ima/(0.0d0,1.0d0)/
c
c
      ima4=-4*ima
      ima4inv=ima/4
c
      allocate(jval(0:nterms+10))
      allocate(jder(0:nterms+10))
c
      allocate(mptemp(-nterms-2:nterms+2))

      do k=1,ntarg
         zdiff(1)=ztarg(1,k)-center(1)
         zdiff(2)=ztarg(2,k)-center(2)
         call h2cart2polar(zdiff,r,theta)
         z=zk*r
         ifder=0
         call jbessel2d(nterms+3,z,rscale,jval,ifder,jder)
c
         zmul=exp(ima*theta)
         zinv=conjg(zmul)
c
         call mpole_evalp_vec(nd,zmul,zinv,local,mptemp,jval,
     $        nterms,pot1(1,k)) 
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
      subroutine h2dtaevalg_vec(nd,zk,rscale,center,local,nterms,
     1           ztarg,ntarg,pot1,grad1)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials and gradients
c     due to incoming partial wave expansions.
c
c                     +nterms
c     POT1  = POT1 +    sum   local_n J_n(k r) exp(i n theta)
c                     n=-nterms  
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of local expansions)
c     zk           :    Helmholtz parameter
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
      integer i,j,k,n,nterms,ntarg,ifder,ii,nd
      real *8 rscale,r,theta,center(2),ztarg(2,*),zdiff(2),rsinv
      complex *16 zk,local(nd,-nterms:nterms)
      complex *16 pot1(nd,ntarg),grad1(nd,2,ntarg)
c
      complex *16, allocatable :: jval(:)
      complex *16, allocatable :: jder(:)
      complex *16, allocatable :: localx(:,:)
      complex *16, allocatable :: localy(:,:)
      complex *16, allocatable :: mptemp(:)
c
      complex *16 ima,ima4,ima4inv,z,z1scale,z2scale,z3scale,z4scale
      complex *16 zmul,zinv,ztemp1,ztemp2
      data ima/(0.0d0,1.0d0)/
c
c
      ima4=-4*ima
      ima4inv=ima/4
c
      allocate(jval(0:nterms+10))
      allocate(jder(0:nterms+10))
      allocate(localx(nd,-nterms-1:nterms+1))
      allocate(localy(nd,-nterms-1:nterms+1))
      allocate(mptemp(-nterms-2:nterms+2))
c
      rsinv = 1.0d0/rscale
      call mk_mpoleg_vec(nd,zk,rsinv,local,localx,localy,nterms)
c
      do k=1,ntarg
         zdiff(1)=ztarg(1,k)-center(1)
         zdiff(2)=ztarg(2,k)-center(2)
         call h2cart2polar(zdiff,r,theta)
c
         z=zk*r
         ifder=0
         call jbessel2d(nterms+3,z,rscale,jval,ifder,jder)
c
         zmul=exp(ima*theta)
         zinv=conjg(zmul)
c
         call mpole_evalp_vec(nd,zmul,zinv,local,mptemp,jval,
     $        nterms,pot1(1,k)) 
c
         call mpole_evalg_vec(nd,localx,localy,mptemp,nterms,
     $        grad1(1,1,k)) 
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
      subroutine h2dtaevalh_vec(nd,zk,rscale,center,local,nterms,
     1           ztarg,ntarg,pot1,grad1,hess1)
      implicit none
c**********************************************************************
c
c     This subroutine INCREMENTS the potentials and derivatives
c     due to incoming partial wave expansions.
c
c                     +nterms
c     POT1  = POT1 +    sum   local_n J_n(k r) exp(i n theta)
c                     n=-nterms  
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd           :    vector length (number of local expansions)
c     zk           :    Helmholtz parameter
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
c     hess1(nd,ntarg)  :    hessians at ztarg  (pxx,pxy,pyy) are incremented
c-----------------------------------------------------------------------
      integer i,j,k,n,nterms,ntarg,ifder,ii,nd
      real *8 rscale,r,theta,center(2),ztarg(2,*),zdiff(2),rsinv
      complex *16 zk,local(nd,-nterms:nterms)
      complex *16 pot1(nd,ntarg),grad1(nd,2,ntarg),hess1(nd,3,ntarg)
c
      complex *16, allocatable :: jval(:)
      complex *16, allocatable :: jder(:)
      complex *16, allocatable :: localx(:,:)
      complex *16, allocatable :: localy(:,:)
      complex *16, allocatable :: localxx(:,:)
      complex *16, allocatable :: localxy(:,:)
      complex *16, allocatable :: localyy(:,:)
      complex *16, allocatable :: mptemp(:)
c
      complex *16 ima,ima4,ima4inv,z,z1scale,z2scale,z3scale,z4scale
      complex *16 zmul,zinv,ztemp1,ztemp2
      data ima/(0.0d0,1.0d0)/
c
c
      ima4=-4*ima
      ima4inv=ima/4
c
      allocate(jval(0:nterms+10))
      allocate(jder(0:nterms+10))
      allocate(localx(nd,-nterms-1:nterms+1))
      allocate(localy(nd,-nterms-1:nterms+1))
      allocate(localxx(nd,-nterms-2:nterms+2))
      allocate(localxy(nd,-nterms-2:nterms+2))
      allocate(localyy(nd,-nterms-2:nterms+2))
      allocate(mptemp(-nterms-2:nterms+2))
c
      rsinv = 1.0d0/rscale
      call mk_mpoleg_vec(nd,zk,rsinv,local,localx,localy,nterms)
c
      call mk_mpoleh_vec(nd,zk,rsinv,localx,localy,
     $     localxx,localxy,localyy,nterms)
c
c
c
      do k=1,ntarg
         zdiff(1)=ztarg(1,k)-center(1)
         zdiff(2)=ztarg(2,k)-center(2)
         call h2cart2polar(zdiff,r,theta)
c
         z=zk*r
         ifder=0
         call jbessel2d(nterms+3,z,rscale,jval,ifder,jder)
c
         zmul=exp(ima*theta)
         zinv=conjg(zmul)

         call mpole_evalp_vec(nd,zmul,zinv,local,mptemp,jval,
     $        nterms,pot1(1,k)) 
c
         call mpole_evalg_vec(nd,localx,localy,mptemp,nterms,
     $        grad1(1,1,k)) 
c
         call mpole_evalh_vec(nd,localxx,localxy,localyy,mptemp,nterms,
     $        hess1(1,1,k)) 
c
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine h2dmpmp_vec(nd,zk,rscale1,center1,hexp1,nterms1,
     $                      rscale2,center2,hexp2,nterms2)
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
C     zk     :   Helmholtz parameter
C     rscale1:   scaling parameter for original multipole expansion
C     center1:   center of original multiple expansion
C     hexp1  :   coefficients of original multiple expansions
C     nterms1:   order of original multipole expansion
C     rscale2:   scaling parameter for shifted multipole expansion
C     center2:   center of shifted multipole expansion
C     nterms2:   order of shifted multipole expansion
C---------------------------------------------------------------------
C     OUTPUT:
C
C     hexp2  :  coefficients of shifted multipole expansions
C---------------------------------------------------------------------
      integer nterms1,nterms2,nterms,i,j,ifder,ii,nd
      real *8 rscale1,rscale2,r,theta,center1(2),center2(2),zdiff(2)
      real *8 rsj,rsi5,fs2,pi,done,rsj2,fs3
      complex *16 hexp1(nd,-nterms1:nterms1),hexp2(nd,-nterms2:nterms2)
      complex *16 zk,z,ima, zmul,zinv,ztemp1,ztemp2
      complex *16, allocatable :: jval(:), jder(:), jtemp(:)
c
      data ima/(0.0d0,1.0d0)/
c
      done=1
      pi=4*atan(done)
c
      nterms = nterms1+nterms2
c
      allocate(jval(0:nterms+10))
      allocate(jder(0:nterms+10))
      allocate(jtemp(-nterms-5:nterms+5))
c
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)
      theta=theta-pi
      z=zk*r
      ifder=0
      call jbessel2d(nterms+3,z,rscale1,jval,ifder,jder)
c
      jtemp(0) = jval(0)
      zmul=exp(-ima*theta)
      zinv=conjg(zmul)
      ztemp1= zmul
      ztemp2=-zinv
      do j = 1,nterms
         jtemp( j) = ztemp1*jval(j)
         jtemp(-j) = ztemp2*jval(j)
         ztemp1= ztemp1*zmul
         ztemp2=-ztemp2*zinv
      enddo
c
      do ii = 1,nd
         hexp2(ii,0) = hexp2(ii,0) + hexp1(ii,0)*jtemp(0)
      enddo
      rsj=rscale1**2
      rsj2=rsj
      do j = 1,nterms1
         do ii = 1,nd
           hexp2(ii,0) = hexp2(ii,0)+(hexp1(ii,+j)*jtemp(-j))*rsj
           hexp2(ii,0) = hexp2(ii,0)+(hexp1(ii,-j)*jtemp(+j))*rsj
         enddo
         rsj=rsj*rsj2
      enddo
c
      rsi5=rscale1/rscale2
      do i = 1,nterms2
         do ii = 1,nd
           hexp2(ii,i) = hexp2(ii,i) + hexp1(ii,0)*jtemp(i)*rsi5
           hexp2(ii,-i) = hexp2(ii,-i) + hexp1(ii,0)*jtemp(-i)*rsi5
         enddo
         rsi5=rsi5*rscale1/rscale2
      enddo
      rsi5=rscale1/rscale2
      do i = 1,nterms2
         rsj=rsj2
         do j = 1,min(nterms1,i)
            fs3=rsi5*rsj
            do ii = 1,nd
               hexp2(ii,i) = hexp2(ii,i)+
     $         (hexp1(ii,+j)*jtemp(i-j))*rsi5
               hexp2(ii,i) = hexp2(ii,i)+
     $         (hexp1(ii,-j)*jtemp(i+j))*fs3
               hexp2(ii,-i) = hexp2(ii,-i)+
     $         (hexp1(ii,+j)*jtemp(-i-j))*fs3
               hexp2(ii,-i) = hexp2(ii,-i)+
     $         (hexp1(ii,-j)*jtemp(-i+j))*rsi5
            enddo
            rsj=rsj*rsj2
         enddo
         rsj=rsj2**(i+1)
         fs2=rsi5*rsj2
         fs3=rsi5*rsj
         do j = i+1,nterms1
            do ii = 1,nd
               hexp2(ii,i) = hexp2(ii,i)+
     $           (hexp1(ii,+j)*jtemp(i-j))*fs2
               hexp2(ii,i) = hexp2(ii,i)+
     $           (hexp1(ii,-j)*jtemp(i+j))*fs3
               hexp2(ii,-i) = hexp2(ii,-i)+
     $           (hexp1(ii,+j)*jtemp(-i-j))*fs3
               hexp2(ii,-i) = hexp2(ii,-i)+
     $           (hexp1(ii,-j)*jtemp(-i+j))*fs2
            enddo
            fs3=fs3*rsj2
            fs2=fs2*rsj2
         enddo
         rsi5=rsi5*rscale1/rscale2
      enddo
      return
      end
c
c
c
c
c
c**********************************************************************
      subroutine h2dlocloc_vec(nd,zk,rscale1,center1,jexp1,nterms1,
     $                          rscale2,center2,jexp2,nterms2)
      implicit none
c**********************************************************************
C
C     This subroutine shifts local expansions JEXP1 to a new center and
C     INCREMENTS the local expansions JEXP2 about that point accordingly.
C
C---------------------------------------------------------------------
C     INPUT:
C
c     nd     :    vector length (number of mpole expansions)
C     zk     :   Helmholtz parameter
C     rscale1:   scaling parameter for local expansion
C     center1:   center of original local expansion
C     jexp1  :   coefficients of original local expansion
C     nterms1:   order of original local expansion
C     rscale2:   scaling parameter for shifted local expansion
C     center2:   center of shifted local expansion
C     nterms2:   order of shifted local expansion
C---------------------------------------------------------------------
C     OUTPUT:
C
C     jexp2  :   coefficients of shifted local expansion
C---------------------------------------------------------------------
      integer nterms1,nterms2,nterms,i,j,ifder,ii,nd
      real *8 rscale1,rscale2,r,theta,center1(2),center2(2),zdiff(2)
      real *8 rsi,rsj,rsi7,rsi5,fs2,pi,done
      complex *16 jexp1(nd,-nterms1:nterms1),jexp2(nd,-nterms2:nterms2)
      complex *16 zk,z,ima, zmul,zinv,ztemp1,ztemp2
      complex *16, allocatable :: jval(:), jder(:), jtemp(:)
c
      data ima/(0.0d0,1.0d0)/
c
      done=1
      pi=4*atan(done)
c
      nterms = nterms1+nterms2
c
      allocate(jval(0:nterms+10))
      allocate(jder(0:nterms+10))
      allocate(jtemp(-nterms-5:nterms+5))
c
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)
c
      theta=theta-pi
c
      z=zk*r
      ifder=0
      call jbessel2d(nterms+3,z,rscale1,jval,ifder,jder)
c
      jtemp(0) = jval(0)
      zmul=exp(-ima*theta)
      zinv=conjg(zmul)
      ztemp1= zmul
      ztemp2=-zinv
      do j = 1,nterms
         jtemp( j) = ztemp1*jval(j)
         jtemp(-j) = ztemp2*jval(j)
         ztemp1= ztemp1*zmul
         ztemp2=-ztemp2*zinv
      enddo
c
      do ii=1,nd
         jexp2(ii,0) = jexp2(ii,0) + jexp1(ii,0)*jtemp(0)
      enddo
      do j = 1,nterms1
      do ii=1,nd
         jexp2(ii,0) = jexp2(ii,0)+(jexp1(ii,+j)*jtemp(-j))
         jexp2(ii,0) = jexp2(ii,0)+(jexp1(ii,-j)*jtemp(+j))
      enddo
      enddo
c
      rsi=rscale1
      rsi7=rscale2
      rsi5=rscale2/rscale1
      do i = 1,nterms2
         do ii=1,nd
            jexp2(ii,i) = jexp2(ii,i) + jexp1(ii,0)*jtemp(i)*rsi7*rsi
            jexp2(ii,-i)= jexp2(ii,-i)+ jexp1(ii,0)*jtemp(-i)*rsi7*rsi
         enddo
         rsi=rsi*rscale1
         rsi7=rsi7*rscale2
      enddo
      rsi=rscale1
      rsi7=rscale2
      rsi5=rscale2/rscale1
      do i = 1,nterms2
         fs2=rsi5
         if( nterms1 .le. i ) fs2=fs2*rscale1**(2*(i-nterms1))
         do j = min(nterms1,i),1,-1
            do ii=1,nd
              jexp2(ii,i)=jexp2(ii,i)+(jexp1(ii,+j)*jtemp(i-j))*fs2
              jexp2(ii,i)=jexp2(ii,i)+(jexp1(ii,-j)*jtemp(i+j))*rsi7*rsi
              jexp2(ii,-i)=jexp2(ii,-i)+
     $                     (jexp1(ii,+j)*jtemp(-i-j))*rsi7*rsi
              jexp2(ii,-i)=jexp2(ii,-i)+(jexp1(ii,-j)*jtemp(-i+j))*fs2
            enddo
            fs2=fs2*rscale1**2
         enddo
         rsj=rscale1**(i+1)
         do j = i+1,nterms1
            do ii=1,nd
               jexp2(ii,i)=jexp2(ii,i)+(jexp1(ii,+j)*jtemp(i-j))*rsi5
               jexp2(ii,i)=jexp2(ii,i)+
     $                     (jexp1(ii,-j)*jtemp(i+j))*rsi7*rsi
               jexp2(ii,-i)=jexp2(ii,-i)+
     $                     (jexp1(ii,+j)*jtemp(-i-j))*rsi7*rsi
               jexp2(ii,-i)=jexp2(ii,-i)+(jexp1(ii,-j)*jtemp(-i+j))*rsi5
            enddo
            rsj=rsj*rscale1
         enddo
         rsi=rsi*rscale1
         rsi7=rsi7*rscale2
         rsi5=rsi5*rscale2/rscale1
      enddo
      return
      end
c
c
c
c
c**********************************************************************
      subroutine h2dmploc_vec(nd,zk,rscale1,center1,hexp,nterms1,
     $                         rscale2,center2,jexp,nterms2)
      implicit none
c**********************************************************************
C
C     This routine maps multipole expansion hexp to local expansion and
C     INCREMENTS jexp accordingly.
C
C---------------------------------------------------------------------
C     INPUT:
C
c     nd     :    vector length (number of mpole expansions)
C     zk     :   Helmholtz parameter
C     rscale1:   scaling parameter for original multipole expansion
C     center1:   center of original multiple expansion
C     hexp   :   coefficients of original multiple expansion
C     nterms1:   order of original multipole expansion
C     rscale2:   scaling parameter for shifted local expansion
C     center2:   center of shifted local expansion
C     nterms2:   order of shifted local expansion
C---------------------------------------------------------------------
C     OUTPUT:
C
C     jexp   :   coefficients of shifted local expansion are incremented
C---------------------------------------------------------------------
      integer nterms1,nterms2,nterms,i,j,ifder,ii,nd
      real *8 rscale1,rscale2,r,theta,center1(2),center2(2),zdiff(2)
      real *8 rsi,rsj,rsi2,rsj2,pi,done
      complex *16 zk,hexp(nd,-nterms1:nterms1),jexp(nd,-nterms2:nterms2)
      complex *16 z,ima, zmul,zinv,ztemp1,ztemp2
      complex *16, allocatable :: hval(:), hder(:), htemp(:)
      data ima/(0.0d0,1.0d0)/
c
      done=1
      pi=4*atan(done)
c
      nterms = nterms1+nterms2
c
      allocate(hval(0:nterms+5))
      allocate(hder(0:nterms+5))
      allocate(htemp(-nterms-5:nterms+5))
c
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)
c
      theta=theta-pi
c
      z=zk*r
      ifder=0
      call h2dall(nterms+1,z,rscale1,hval,ifder,hder)
c        
      htemp(0) = hval(0)
      zmul=exp(-ima*theta)
      zinv=conjg(zmul)
      ztemp1= zmul
      ztemp2=-zinv
      do j = 1,nterms
         htemp( j) = ztemp1*hval(j)
         htemp(-j) = ztemp2*hval(j)
         ztemp1= ztemp1*zmul
         ztemp2=-ztemp2*zinv
      enddo
c
      do ii=1,nd
         jexp(ii,0) = jexp(ii,0) + hexp(ii,0)*htemp(0)
      enddo
      do j = 1,nterms1
        do ii=1,nd
         jexp(ii,0) = jexp(ii,0)+(hexp(ii,+j)*htemp(-j))
         jexp(ii,0) = jexp(ii,0)+(hexp(ii,-j)*htemp(+j))
        enddo
      enddo
c
      rsi=rscale1
      rsi2=rscale1**2
      do i = 1,nterms2
         do ii=1,nd
            jexp(ii,i) = jexp(ii,i) + hexp(ii,0)*htemp(i)
            jexp(ii,-i) = jexp(ii,-i) + hexp(ii,0)*htemp(-i)
         enddo
         rsj=rscale1
         rsj2=rscale1**2
         do j = 1,min(nterms1,i)
            do ii=1,nd
               jexp(ii,i) = jexp(ii,i)+(hexp(ii,+j)*htemp(i-j))*rsj2
               jexp(ii,i) = jexp(ii,i)+(hexp(ii,-j)*htemp(i+j))
               jexp(ii,-i) = jexp(ii,-i)+(hexp(ii,+j)*htemp(-i-j))
               jexp(ii,-i) = jexp(ii,-i)+(hexp(ii,-j)*htemp(-i+j))*rsj2
            enddo
            rsj=rsj*rscale1
            rsj2=rsj2*rscale1**2
         enddo
         do j = i+1,nterms1
            do ii=1,nd
               jexp(ii,i) = jexp(ii,i)+(hexp(ii,+j)*htemp(i-j))*rsi2
               jexp(ii,i) = jexp(ii,i)+(hexp(ii,-j)*htemp(i+j))
               jexp(ii,-i) = jexp(ii,-i)+(hexp(ii,+j)*htemp(-i-j))
               jexp(ii,-i) = jexp(ii,-i)+(hexp(ii,-j)*htemp(-i+j))*rsi2
            enddo
         enddo
         rsi=rsi*rscale1
         rsi2=rsi2*rscale1**2
      enddo
c
      rsi=rscale2/rscale1
      do i = 1,nterms2
         do ii=1,nd
            jexp(ii,+i) = jexp(ii,+i)*rsi
            jexp(ii,-i) = jexp(ii,-i)*rsi
         enddo
         rsi=rsi*rscale2/rscale1
      enddo
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine h2dmpzero_vec(nd,mpole,nterms)
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
      complex *16 mpole(nd,-nterms:nterms)
c
      do n=-nterms,nterms
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
C***********************************************************************
      subroutine ctompole_vec(nd,zmul,zinv,mpole,jval,charge,nterms) 
      implicit none
C***********************************************************************
c
c     This subroutine is called by formmp and formta routines to get 
c     contribution to expansion from a charge source
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              : vector length (number of mpole expansions)
c     zmul, zinv      : parameters set in calling routine from 
c                       source location
c     jval            : radial expansion (computed J_n or H_n values)
c                       When forming h-expansion, should be J series
c                       When forming j-expansion, should be H series
c     charge          : charge value of source
c     nterms          : order of multipole/local expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mpole           : coeffs for the h- or j-expansion are incremented.
c-----------------------------------------------------------------------
      integer n,nterms,nd,ii
      complex *16 mpole(nd,-nterms:nterms),charge(nd)
      complex *16 zmul,zinv,ztemp1,ztemp2
      complex *16 jval(0:nterms)
c      
      do ii = 1,nd
         mpole(ii,0)=mpole(ii,0)+charge(ii)*jval(0)
      enddo
      ztemp1= zmul
      ztemp2=-zinv
      do n=1,nterms
         do ii = 1,nd
            mpole(ii,n)=mpole(ii,n)+jval(n)*ztemp1*charge(ii)
            mpole(ii,-n)=mpole(ii,-n)+jval(n)*ztemp2*charge(ii)
         enddo
         ztemp1= ztemp1*zmul
         ztemp2=-ztemp2*zinv
      enddo
      return
      end

C***********************************************************************
      subroutine dtompole_vec(nd,zk,rscale,zmul,zinv,
     1           mpole,jval,dipstr,dipvec,nterms) 
      implicit none
C***********************************************************************
c
c     This subroutine is called by formmp and formta routines to get 
c     contribution to expansion from a dipole source
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              : vector length (number of mpole expansions)
c     zk              : Helmholtz parameter
c     rscale          : expansion scaling parameter
c                       In calling routine, when forming h-expansion, 
c                       this should be original rscale
c                       In calling routine, when forming j-expansion, 
c                       this should be 1/rscale
c     zmul, zinv      : parameters set in calling routine from 
c                       source location
c     jval            : radial expansion (computed J_n or H_n values)
c                       When forming h-expansion, should be J series
c                       When forming j-expansion, should be H series
c     dipstr          : source strength
c     dipvec          : dipole orientation
c     nterms          : order of multipole/local expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mpole           : coeffs for the h- or j-expansion are incremented.
c-----------------------------------------------------------------------
      integer n,nterms,nd,ii
      real *8 dipvec(nd,2),rscale
      complex *16 mpole(nd,-nterms:nterms),dipstr(nd),ima
      complex *16 zk,zmul,zinv,ztemp1,ztemp2,ztemp3,ztemp4,ztemp5,ztemp6
      complex *16 ztemp7,ztemp8,ztemp9,ztemp10
      complex *16 jval(0:nterms+1)
      data ima/(0.0d0,1.0d0)/
c
      ztemp3= zinv/rscale
      ztemp4= zmul*rscale
      ztemp5= zmul/rscale
      ztemp6= zinv*rscale
      do ii = 1,nd
         mpole(ii,0)=mpole(ii,0)-dipstr(ii)*jval(1)*zk/2*rscale*
     $         ( (zmul+zinv)*dipvec(ii,1)+
     $         (zmul-zinv)*ima*dipvec(ii,2) )
      enddo
      ztemp1=-zmul*zk/2
      ztemp2=+zinv*zk/2
      do n=1,nterms
         do ii = 1,nd
            ztemp7 = ztemp3*(-dipvec(ii,1)+ima*dipvec(ii,2))
            ztemp8 = ztemp4*(+dipvec(ii,1)+ima*dipvec(ii,2))
            ztemp9 = ztemp5*(-dipvec(ii,1)-ima*dipvec(ii,2))
            ztemp10 = ztemp6*(+dipvec(ii,1)-ima*dipvec(ii,2))
            mpole(ii,n)=mpole(ii,n)+
     $         (jval(n-1)*ztemp7+jval(n+1)*ztemp8)*ztemp1*dipstr(ii)
            mpole(ii,-n)=mpole(ii,-n)+
     $         (jval(n-1)*ztemp9+jval(n+1)*ztemp10)*ztemp2*dipstr(ii)
         enddo
         ztemp1= ztemp1*zmul
         ztemp2=-ztemp2*zinv
      enddo
      return
      end

C***********************************************************************
      subroutine mpole_evalp_vec(nd,zmul,zinv,mpole,mptemp,hval,
     $           nterms,pot1) 
      implicit none
C***********************************************************************
c
c     This subroutine is called by mpeval and taeval routines to get 
c     the potential
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              : vector length (number of mpole expansions)
c     zmul, zinv      : parameters set in calling routine from 
c                       source location
c     mpole           : coeffs for the h- or j-expansion
c     hval            : radial expansion (computed H_n or J_n values)
c                       When evaluating h-expansion, should be H series
c                       When evaluating j-expansion, should be J series
c     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mptemp          : output used also by gradient and Hessian eval codes.
c     pot1            : potential is incremented
c-----------------------------------------------------------------------
      integer n,nterms,j,nd,ii,k
      complex *16 mpole(nd,-nterms:nterms),ima,ima4inv
      complex *16 zmul,zinv,ztemp1,ztemp2
      complex *16 hval(0:nterms+2),pot1(nd)
      complex *16 mptemp(-nterms-2:nterms+2)
      data ima/(0.0d0,1.0d0)/
c      
      ima4inv=ima/4
      mptemp(0)=hval(0)
      ztemp1= zmul*ima4inv
      ztemp2=-zinv*ima4inv
      do j = 1,nterms+2
         mptemp( j) = ztemp1*hval(j)
         mptemp(-j) = ztemp2*hval(j)
         ztemp1= ztemp1*zmul
         ztemp2=-ztemp2*zinv
      enddo
      do ii=1,nd
         pot1(ii)=pot1(ii)+hval(0)*mpole(ii,0)*ima4inv
      enddo
      do n=1,nterms
         do ii = 1,nd
            pot1(ii)=pot1(ii)+
     $        mpole(ii,n)*mptemp(n)+mpole(ii,-n)*mptemp(-n)
         enddo
      enddo
      return
      end

C***********************************************************************
      subroutine mpole_evalg_vec(nd,mpolex,mpoley,mptemp,nterms,grad1) 
      implicit none
C***********************************************************************
c
c     This subroutine is called by mpeval and taeval routines to get 
c     the gradient
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              : vector length (number of mpole expansions)
c     mpolex,mpoley   : coeffs for the expansions of x,y derivatives
c     mptemp          : array computed in mpole_evalp_vec
c     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     grad1            : gradient is incremented
c-----------------------------------------------------------------------
c-----------------------------------------------------------------------
      integer n,nterms,ii,nd
      complex *16 ima,ima4inv,grad1(nd,2)
      complex *16 mpolex(nd,-nterms-1:nterms+1)
      complex *16 mpoley(nd,-nterms-1:nterms+1)
      complex *16 mptemp(-nterms-2:nterms+2)
      data ima/(0.0d0,1.0d0)/
c      
      ima4inv=ima/4
      do ii=1,nd
         grad1(ii,1)=grad1(ii,1)+mptemp(0)*mpolex(ii,0)*ima4inv
         grad1(ii,2)=grad1(ii,2)+mptemp(0)*mpoley(ii,0)*ima4inv
      enddo
      do n=1,nterms+1
         do ii = 1,nd
            grad1(ii,1)=grad1(ii,1)+
     $           mpolex(ii,n)*mptemp(n)+mpolex(ii,-n)*mptemp(-n)
            grad1(ii,2)=grad1(ii,2)+
     $           mpoley(ii,n)*mptemp(n)+mpoley(ii,-n)*mptemp(-n)
         enddo
      enddo
      return
      end

C***********************************************************************
      subroutine mpole_evalh_vec(nd,mpolexx,mpolexy,mpoleyy,mptemp,
     $           nterms,hess1) 
      implicit none
C***********************************************************************
c
c     This subroutine is called by mpeval and taeval routines to get 
c     the Hessian
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              : vector length (number of mpole expansions)
c     mpolexx         
c     mpolexy         : coeffs for the expansions of x,y derivatives
c     mpoleyy
c     mptemp          : array computed in mpole_evalp_vec
c     nterms          : order of multipole expansion
c     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     hess1           : Hessian is incremented
c-----------------------------------------------------------------------
      integer n,nterms,ii,nd
      complex *16 ima,ima4inv,hess1(nd,3)
      complex *16 mpolexx(nd,-nterms-2:nterms+2)
      complex *16 mpolexy(nd,-nterms-2:nterms+2)
      complex *16 mpoleyy(nd,-nterms-2:nterms+2)
      complex *16 mptemp(-nterms-2:nterms+2)
      data ima/(0.0d0,1.0d0)/
c      
      ima4inv=ima/4
      do ii=1,nd
         hess1(ii,1)= hess1(ii,1)+mptemp(0)*mpolexx(ii,0)*ima4inv
         hess1(ii,2)= hess1(ii,2)+mptemp(0)*mpolexy(ii,0)*ima4inv
         hess1(ii,3)= hess1(ii,3)+mptemp(0)*mpoleyy(ii,0)*ima4inv
      enddo
      do n=1,nterms+1
         do ii = 1,nd
            hess1(ii,1)=hess1(ii,1)+
     $              mpolexx(ii,n)*mptemp(n)+mpolexx(ii,-n)*mptemp(-n)
            hess1(ii,2)=hess1(ii,2)+
     $              mpolexy(ii,n)*mptemp(n)+mpolexy(ii,-n)*mptemp(-n)
            hess1(ii,3)=hess1(ii,3)+
     $              mpoleyy(ii,n)*mptemp(n)+mpoleyy(ii,-n)*mptemp(-n)
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
C***********************************************************************
      subroutine mk_mpoleg_vec(nd,zk,rscale,mpole,mpolex,mpoley,nterms)
      implicit none
C***********************************************************************
c
c     This subroutine is called by mpeval and taeval routines to create 
c     coefficient expansions for first derivatives
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd              : vector length (number of mpole expansions)
c     zk              : Helmholtz parameter
c     rscale          : expansion scaling parameter
c                       In calling routine, when forming h-expansion, 
c                       this should be original rscale
c                       In calling routine, when forming j-expansion, 
c                       this should be 1/rscale
c     mpole           : coeffs for the h- or j- expansion
c     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mpolex
c     mpoley          : coeffs for the expansions of second derivatives
c-----------------------------------------------------------------------
      integer n,nterms,i,ii,nd
      real *8 rscale
      complex *16 mpole(nd,-nterms:nterms)
      complex *16 zk,z1scale,z2scale,z3scale,z4scale,ima
      complex *16 mpolex(nd,-nterms-1:nterms+1)
      complex *16 mpoley(nd,-nterms-1:nterms+1)
      data ima/(0.0d0,1.0d0)/

      do i=-nterms-1,nterms+1
         do ii=1,nd
            mpolex(ii,i)=0
            mpoley(ii,i)=0
         enddo
      enddo
      z1scale = zk/2/rscale
      z2scale = zk/2*rscale
      z3scale = zk/2/rscale * ima
      z4scale = zk/2*rscale * ima

      do ii=1,nd
         mpolex(ii,-1)=mpolex(ii,-1)+z1scale*mpole(ii,0)
         mpoley(ii,-1)=mpoley(ii,-1)+z3scale*mpole(ii,0)
         mpolex(ii,1)=mpolex(ii,1)-z1scale*mpole(ii,0)
         mpoley(ii,1)=mpoley(ii,1)+z3scale*mpole(ii,0)
      enddo
c
      do i=1,nterms
      do ii=1,nd
         mpolex(ii,i-1)=mpolex(ii,i-1)+z2scale*mpole(ii,i)
         mpoley(ii,i-1)=mpoley(ii,i-1)+z4scale*mpole(ii,i)
         mpolex(ii,i+1)=mpolex(ii,i+1)-z1scale*mpole(ii,i)
         mpoley(ii,i+1)=mpoley(ii,i+1)+z3scale*mpole(ii,i)
      enddo
      enddo
c
      do i=-nterms,-1,1
      do ii=1,nd
         mpolex(ii,i-1)=mpolex(ii,i-1)+z1scale*mpole(ii,i)
         mpoley(ii,i-1)=mpoley(ii,i-1)+z3scale*mpole(ii,i)
         mpolex(ii,i+1)=mpolex(ii,i+1)-z2scale*mpole(ii,i)
         mpoley(ii,i+1)=mpoley(ii,i+1)+z4scale*mpole(ii,i)
      enddo
      enddo
c
      return
      end
c
c
c
c
C***********************************************************************
      subroutine mk_mpoleh_vec(nd,zk,rscale,mpolex,mpoley,
     $           mpolexx,mpolexy,mpoleyy,nterms)
      implicit none
C***********************************************************************
c
c     This subroutine is called by mpeval and taeval routines to create 
c     coefficient expansions for second derivatives
c-----------------------------------------------------------------------
c     INPUT:
c
c     nd          c   : vector length (number of mpole expansions)
c     zk              : Helmholtz parameter
c     rscale          : expansion scaling parameter
c                       In calling routine, when forming h-expansion, 
c                       this should be original rscale
c                       In calling routine, when forming j-expansion, 
c                       this should be 1/rscale
c     mpolex,mpoley   : coeffs for the expansions of x,y derivatives
c     nterms          : order of multipole expansion
c-----------------------------------------------------------------------
c     OUTPUT:
c
c     mpolexx
c     mpolexy          : coeffs for the expansions of second derivatives
c     mpoleyy
c-----------------------------------------------------------------------
      integer n,nterms,i,ii,nd
      real *8 rscale
      complex *16 mpole(nd,-nterms:nterms)
      complex *16 zk,z1scale,z2scale,z3scale,z4scale,ima
      complex *16 mpolex(nd,-nterms-1:nterms+1)
      complex *16 mpoley(nd,-nterms-1:nterms+1)
      complex *16 mpolexx(nd,-nterms-2:nterms+2)
      complex *16 mpolexy(nd,-nterms-2:nterms+2)
      complex *16 mpoleyy(nd,-nterms-2:nterms+2)
      data ima/(0.0d0,1.0d0)/

      do i=-nterms-2,nterms+2
      do ii=1,nd
         mpolexx(ii,i)=0
         mpolexy(ii,i)=0
         mpoleyy(ii,i)=0
      enddo
      enddo
      z1scale = zk/2/rscale
      z2scale = zk/2*rscale
      z3scale = zk/2/rscale * ima
      z4scale = zk/2*rscale * ima
c
      do ii=1,nd
         mpolexx(ii,-1)=mpolexx(ii,-1)+z1scale*mpolex(ii,0)
         mpolexy(ii,-1)=mpolexy(ii,-1)+z3scale*mpolex(ii,0)
         mpoleyy(ii,-1)=mpoleyy(ii,-1)+z3scale*mpoley(ii,0)
         mpolexx(ii,1)=mpolexx(ii,1)-z1scale*mpolex(ii,0)
         mpolexy(ii,1)=mpolexy(ii,1)+z3scale*mpolex(ii,0)
         mpoleyy(ii,1)=mpoleyy(ii,1)+z3scale*mpoley(ii,0)
      enddo
      do i=1,nterms+1
      do ii=1,nd
         mpolexx(ii,i-1)=mpolexx(ii,i-1)+z2scale*mpolex(ii,i)
         mpolexy(ii,i-1)=mpolexy(ii,i-1)+z4scale*mpolex(ii,i)
         mpoleyy(ii,i-1)=mpoleyy(ii,i-1)+z4scale*mpoley(ii,i)
         mpolexx(ii,i+1)=mpolexx(ii,i+1)-z1scale*mpolex(ii,i)
         mpolexy(ii,i+1)=mpolexy(ii,i+1)+z3scale*mpolex(ii,i)
         mpoleyy(ii,i+1)=mpoleyy(ii,i+1)+z3scale*mpoley(ii,i)
      enddo
      enddo
c
      do i=-nterms+1,-1,1
      do ii=1,nd
         mpolexx(ii,i-1)=mpolexx(ii,i-1)+z1scale*mpolex(ii,i)
         mpolexy(ii,i-1)=mpolexy(ii,i-1)+z3scale*mpolex(ii,i)
         mpoleyy(ii,i-1)=mpoleyy(ii,i-1)+z3scale*mpoley(ii,i)
         mpolexx(ii,i+1)=mpolexx(ii,i+1)-z2scale*mpolex(ii,i)
         mpolexy(ii,i+1)=mpolexy(ii,i+1)+z4scale*mpolex(ii,i)
         mpoleyy(ii,i+1)=mpoleyy(ii,i+1)+z4scale*mpoley(ii,i)
      enddo
      enddo
      return
      end
c
c

