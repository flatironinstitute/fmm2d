cc Copyright (C) 2018-2019: Leslie Greengard, Zydrunas Gimbutas, 
cc and Manas Rachh
cc Contact: greengard@cims.nyu.edu
cc
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
cc 2021/07/07: convert to proper Laplace FMM, calls Cauchy FMM
cc                   Travis Askham
c
c
c     This is the Laplace FMM with real-valued charges and
c     dipole strengths. 
c
c      
      
      subroutine rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1     ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2     nt,targ,ifpghtarg,pottarg,gradtarg,
     3     hesstarg,ier)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   ifcharge      : flag for including charge interactions
c                   charge interactions included if ifcharge =1
c                   not included otherwise
c   charge(nd,ns)    : charge strengths
c   ifdipole      : flag for including dipole interactions
c                   dipole interactions included if ifcharge =1
c                   not included otherwise
c   dipstr(nd,ns) : dipole strengths
c   dipvec(nd,2,ns) : dipole strengths      
c   iper          : flag for periodic implmentations. Currently unused
c   ifpgh         : flag for computing pot/grad/hess
c                   ifpgh = 1, only potential is computed
c                   ifpgh = 2, potential and gradient are computed
c                   ifpgh = 3, potential, gradient, and hessian 
c                   are computed
c   nt            : number of targets
c   targ(2,nt)    : target locations
c   ifpghtarg     : flag for computing pottarg/gradtarg/hesstarg
c                   ifpghtarg = 1, only potential is computed at targets
c                   ifpghtarg = 2, potential and gradient are 
c                   computed at targets
c                   ifpghtarg = 3, potential, gradient, and hessian are 
c                   computed at targets
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : potential at the source locations
c   grad(nd,2,*)    : gradients at the source locations
c   hess(nd,3,*)    : hessian at the source locations
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,2,*): gradient at the target locations
c   hesstarg(nd,3,*): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 omp_get_wtime
      real *8 eps
      integer ns,nt,iper
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(nd,*),dipstr(nd,*)
      real *8 dipvec(nd,2,*)

      real *8 pot(nd,*),grad(nd,2,*),hess(nd,3,*)
      real *8 pottarg(nd,*),gradtarg(nd,2,*),hesstarg(nd,3,*)

      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg,ifprint,ier
      real *8 time1,time2,pi,done

      complex *16, allocatable, dimension(:,:) :: pot1,grad1,hess1,
     1     pottarg1,gradtarg1,hesstarg1,charge1,dipstr1

      integer i,j,nd2
      
      complex *16 :: eye, ztmp
      data eye /(0.0d0,1.0d0)/
      
c     allocate variables for cfmm call
      
      if (ifcharge .eq. 1) then
         allocate(charge1(nd,ns))
      else
         allocate(charge1(nd,1))
      endif

      if (ifdipole .eq. 1) then
         allocate(dipstr1(nd,ns))
      else
         allocate(dipstr1(nd,1))
      endif

      if (ifpgh .eq. 1) then
         allocate(pot1(nd,ns),grad1(nd,1),hess1(nd,1))
      endif
      if (ifpgh .eq. 2) then
         allocate(pot1(nd,ns),grad1(nd,ns),hess1(nd,1))
      endif
      if (ifpgh .eq. 3) then
         allocate(pot1(nd,ns),grad1(nd,ns),hess1(nd,ns))
      endif

      if (ifpghtarg .eq. 1) then
         allocate(pottarg1(nd,nt),gradtarg1(nd,1),hesstarg1(nd,1))
      endif
      if (ifpghtarg .eq. 2) then
         allocate(pottarg1(nd,nt),gradtarg1(nd,nt),
     1        hesstarg1(nd,1))
      endif
      if (ifpghtarg .eq. 3) then
         allocate(pottarg1(nd,nt),gradtarg1(nd,nt),
     1        hesstarg1(nd,nt))
      endif

c     real and complex parts must be split

      if (ifcharge .eq. 1) then
         do i = 1,ns
            do j = 1,nd
               charge1(j,i) = charge(j,i)
            enddo
         enddo
      endif
      
      if (ifdipole .eq. 1) then
         do i = 1,ns
            do j = 1,nd
               ztmp = -(dipvec(j,1,i)+eye*dipvec(j,2,i))
               dipstr1(j,i) = dipstr(j,i)*ztmp
            enddo
         enddo
      endif

c     cfmm does the work

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge1,
     1     ifdipole,dipstr1,iper,ifpgh,pot1,grad1,hess1,
     2     nt,targ,ifpghtarg,pottarg1,gradtarg1,
     3     hesstarg1,ier)

c     unpack the d/dz, d^2/dz^2 as grad/hess and
c     combine real and imaginary parts

      if (ifpgh .eq. 1 .or. ifpgh .eq. 2 .or. ifpgh .eq. 3) then
         do i = 1,ns
            do j = 1,nd
               pot(j,i) = dble(pot1(j,i))
            enddo
         enddo
      endif
      if (ifpgh .eq. 2 .or. ifpgh .eq. 3) then
         do i = 1,ns
            do j = 1,nd
               grad(j,1,i) = dble(grad1(j,i))
               grad(j,2,i) = -imag(grad1(j,i))
            enddo
         enddo
      endif
      if (ifpgh .eq. 3) then
         do i = 1,ns
            do j = 1,nd
               hess(j,1,i) = dble(hess1(j,i))
               hess(j,2,i) = -imag(hess1(j,i))
               hess(j,3,i) = -hess(j,1,i)
            enddo
         enddo
      endif
      
      if (ifpghtarg .eq. 1 .or. ifpghtarg .eq. 2
     1     .or. ifpghtarg .eq. 3) then
         do i = 1,nt
            do j = 1,nd
               pottarg(j,i) = dble(pottarg1(j,i))
            enddo
         enddo
      endif
      if (ifpghtarg .eq. 2 .or. ifpghtarg .eq. 3) then
         do i = 1,nt
            do j = 1,nd
               gradtarg(j,1,i) = dble(gradtarg1(j,i))
               gradtarg(j,2,i) = -imag(gradtarg1(j,i))
            enddo
         enddo
      endif
      if (ifpghtarg .eq. 3) then
         do i = 1,nt
            do j = 1,nd
               hesstarg(j,1,i) = dble(hesstarg1(j,i))
               hesstarg(j,2,i) = -imag(hesstarg1(j,i))
               hesstarg(j,3,i) = -hesstarg(j,1,i)
            enddo
         enddo
      endif
      
      return
      end
c
c
c
c------------------------------------------------------------------
      subroutine rfmm2dpart_direct(nd,istart,iend,jstart,jend,
     $     source,ifcharge,charge,ifdipole,dipstr,dipvec,
     $     targ,ifpgh,pot,grad,hess,thresh)
c--------------------------------------------------------------------
c     This subroutine adds the contribuition due to sources
c     istart to iend in the source array at the expansion centers
c     jstart to jend in the target array to the computed velocities
c     and gradients. Note that contributions for sources
c     within thresh of the targets are not added to the potential
c     
c     
c     INPUT arguments
c-------------------------------------------------------------------
c     nd           in: integer
c     number of charge densities
c     
c     istart       in:Integer
c     Starting index in source array whose expansions
c     we wish to add
c     
c     iend         in:Integer
c     Last index in source array whose expansions
c     we wish to add
c     
c     jstart       in: Integer
c     First index in target array at which we
c     wish to update the potential and gradients
c     
c     jend         in:Integer
c     Last index in target array at which we wish
c     to update the potential and gradients
c     
c     source       in: real *8(2,ns)
c     Source locations
c     
c     ifcharge     in: Integer
c     flag for including expansions due to charges
c     The expansion due to charges will be included
c     if ifcharge == 1
c     
c     charge       in: real *8
c     Charge at the source locations
c     
c     ifdipole     in: Integer
c     flag for including expansions due to dipoles
c     The expansion due to dipoles will be included
c     if ifdipole == 1
c     
c     dipstr        in: real *8(ns)
c     dipole strengths at the source locations
c
c     dipvec        in: real *8 (nd,2,ns)
c     
c     targ        in: real *8(2,nt)
c     target locations
c     
c     ifpgh        in: Integer
c     Flag for computing the potential/gradient/hessian.
c     ifpgh = 1, only potential is computed
c     ifpgh = 2, potential/gradient are computed
c     ifpgh = 3, potential/gradient/hessian are computed
c     
c     thresh       in: real *8
c     threshold for computing interactions
c     if |r| < threshold, then interactions are
c     not included
c     
c     
c------------------------------------------------------------
c     OUTPUT
c     
c     Updated velocity and gradients at the targets
c     pot : potential at the targets
c     grad: gradient at the targets
c     hess: Hessian at the targets
c-------------------------------------------------------
      implicit none
c     
      integer istart,iend,jstart,jend,ns,j,i,nt
      integer ifcharge,ifdipole

      integer nd



      real *8 source(2,*)
      real *8 charge(nd,*),dipstr(nd,*)
      real *8 dipvec(nd,2,*)

      integer ifpgh
      real *8 targ(2,*),thresh
      
c     
      real *8 pot(nd,*)
      real *8 grad(nd,2,*)
      real *8 hess(nd,3,*)

c     
      ns = iend - istart + 1
      nt = jend - jstart + 1
      if(ifcharge.eq.1.and.ifdipole.eq.0) then
         if(ifpgh.eq.1) then
            call r2d_directcp(nd,source(1,istart),ns,
     1           charge(1,istart),targ(1,jstart),nt,
     2           pot(1,jstart),thresh)
         endif

         if(ifpgh.eq.2) then
            call r2d_directcg(nd,source(1,istart),ns,
     1           charge(1,istart),targ(1,jstart),nt,
     2           pot(1,jstart),grad(1,1,jstart),thresh)
         endif
         if(ifpgh.eq.3) then
            call r2d_directch(nd,source(1,istart),ns,
     1           charge(1,istart),targ(1,jstart),nt,
     2           pot(1,jstart),grad(1,1,jstart),hess(1,1,jstart),
     3           thresh)
         endif
      endif

      if(ifcharge.eq.0.and.ifdipole.eq.1) then
         if(ifpgh.eq.1) then
            call r2d_directdp(nd,source(1,istart),ns,
     1           dipstr(1,istart),dipvec(1,1,istart),
     2           targ(1,jstart),nt,pot(1,jstart),thresh)
         endif

         if(ifpgh.eq.2) then
            call r2d_directdg(nd,source(1,istart),ns,
     1           dipstr(1,istart),dipvec(1,1,istart),
     2           targ(1,jstart),nt,pot(1,jstart),grad(1,1,jstart),
     2           thresh)
         endif
         if(ifpgh.eq.3) then
            call r2d_directdh(nd,source(1,istart),ns,
     1           dipstr(1,istart),dipvec(1,1,istart),targ(1,jstart),nt,
     2           pot(1,jstart),grad(1,1,jstart),
     2           hess(1,1,jstart),thresh)
         endif
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.1) then
         if(ifpgh.eq.1) then
            call r2d_directcdp(nd,source(1,istart),ns,
     1           charge(1,istart),dipstr(1,istart),
     2           dipvec(1,1,istart),targ(1,jstart),nt,
     3           pot(1,jstart),thresh)
         endif

         if(ifpgh.eq.2) then
            call r2d_directcdg(nd,source(1,istart),ns,
     1           charge(1,istart),dipstr(1,istart),
     2           dipvec(1,1,istart),targ(1,jstart),nt,
     3           pot(1,jstart),grad(1,1,jstart),
     4           thresh)
         endif
         if(ifpgh.eq.3) then
            call r2d_directcdh(nd,source(1,istart),ns,
     1           charge(1,istart),dipstr(1,istart),
     2           dipvec(1,1,istart),targ(1,jstart),nt,
     3           pot(1,jstart),grad(1,1,jstart),
     2           hess(1,1,jstart),thresh)
         endif
      endif


c     
      return
      end
