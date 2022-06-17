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
c     This is the Laplace FMM with complex-valued charges and
c     dipole strengths. It is kept mostly for historical purposes
c     as its function can also be accomplished with a vectorized
c     rfmm2d call
c
c      
      
      subroutine lfmm2d_ndiv(nd,eps,ns,sources,ifcharge,charge,
     1     ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2     nt,targ,ifpghtarg,pottarg,gradtarg,
     3     hesstarg,ndiv,idivflag,ifnear,timeinfo,ier)
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
c   ndiv          : subdivision criterion, number of points of box
c   idivflag      : subdivision criterion, whether to use sources, targets
c
c   OUTPUT PARAMETERS
c   pot(nd,*)       : potential at the source locations
c   grad(nd,2,*)    : gradients at the source locations
c   hess(nd,3,*)    : hessian at the source locations
c   pottarg(nd,*)   : potential at the target locations
c   gradtarg(nd,2,*): gradient at the target locations
c   hesstarg(nd,3,*): hessian at the target locations
c   timeinfo        : time distribution
c   ier             : error code
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
      complex *16 charge(nd,*),dipstr(nd,*)
      real *8 dipvec(nd,2,*)

      complex *16 pot(nd,*),grad(nd,2,*),hess(nd,3,*)
      complex *16 pottarg(nd,*),gradtarg(nd,2,*),hesstarg(nd,3,*)

      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg,ifprint,ier
      integer ndiv,idivflag,ifnear
      real *8 timeinfo(8)
      real *8 time1,time2,pi,done

      complex *16, allocatable, dimension(:,:,:) :: pot1,grad1,hess1,
     1     pottarg1,gradtarg1,hesstarg1,charge1,dipstr1

      integer i,j,nd2
      
      complex *16 :: eye, ztmp
      data eye /(0.0d0,1.0d0)/
      
c     allocate variables for cfmm call
      
      if (ifcharge .eq. 1) then
         allocate(charge1(2,nd,ns))
      else
         allocate(charge1(2,nd,1))
      endif

      if (ifdipole .eq. 1) then
         allocate(dipstr1(2,nd,ns))
      else
         allocate(dipstr1(2,nd,1))
      endif

      if (ifpgh .eq. 1) then
         allocate(pot1(2,nd,ns),grad1(2,nd,1),hess1(2,nd,1))
      endif
      if (ifpgh .eq. 2) then
         allocate(pot1(2,nd,ns),grad1(2,nd,ns),hess1(2,nd,1))
      endif
      if (ifpgh .eq. 3) then
         allocate(pot1(2,nd,ns),grad1(2,nd,ns),hess1(2,nd,ns))
      endif

      if (ifpghtarg .eq. 1) then
         allocate(pottarg1(2,nd,nt),gradtarg1(2,nd,1),hesstarg1(2,nd,1))
      endif
      if (ifpghtarg .eq. 2) then
         allocate(pottarg1(2,nd,nt),gradtarg1(2,nd,nt),
     1        hesstarg1(2,nd,1))
      endif
      if (ifpghtarg .eq. 3) then
         allocate(pottarg1(2,nd,nt),gradtarg1(2,nd,nt),
     1        hesstarg1(2,nd,nt))
      endif

c     real and complex parts must be split

      if (ifcharge .eq. 1) then
         do i = 1,ns
            do j = 1,nd
               charge1(1,j,i) = dble(charge(j,i))
               charge1(2,j,i) = imag(charge(j,i))
            enddo
         enddo
      endif
      
      if (ifdipole .eq. 1) then
         do i = 1,ns
            do j = 1,nd
               ztmp = -(dipvec(j,1,i)+eye*dipvec(j,2,i))
               dipstr1(1,j,i) = dble(dipstr(j,i))*ztmp
               dipstr1(2,j,i) = imag(dipstr(j,i))*ztmp
            enddo
         enddo
      endif

c     cfmm does the work

      nd2 = 2*nd
      call cfmm2d_ndiv(nd2,eps,ns,sources,ifcharge,charge1,
     1     ifdipole,dipstr1,iper,ifpgh,pot1,grad1,hess1,
     2     nt,targ,ifpghtarg,pottarg1,gradtarg1,
     3     hesstarg1,ndiv,idivflag,ifnear,timeinfo,ier)

c     unpack the d/dz, d^2/dz^2 as grad/hess and
c     combine real and imaginary parts

      if (ifpgh .eq. 1 .or. ifpgh .eq. 2 .or. ifpgh .eq. 3) then
         do i = 1,ns
            do j = 1,nd
               pot(j,i) = dble(pot1(1,j,i))+eye*dble(pot1(2,j,i))
            enddo
         enddo
      endif
      if (ifpgh .eq. 2 .or. ifpgh .eq. 3) then
         do i = 1,ns
            do j = 1,nd
               grad(j,1,i) = dble(grad1(1,j,i))+eye*dble(grad1(2,j,i))
               grad(j,2,i) = -imag(grad1(1,j,i))-eye*imag(grad1(2,j,i))
            enddo
         enddo
      endif
      if (ifpgh .eq. 3) then
         do i = 1,ns
            do j = 1,nd
               hess(j,1,i) = dble(hess1(1,j,i))+eye*dble(hess1(2,j,i))
               hess(j,2,i) = -imag(hess1(1,j,i))-eye*imag(hess1(2,j,i))
               hess(j,3,i) = -hess(j,1,i)
            enddo
         enddo
      endif
      
      if (ifpghtarg .eq. 1 .or. ifpghtarg .eq. 2
     1     .or. ifpghtarg .eq. 3) then
         do i = 1,nt
            do j = 1,nd
               pottarg(j,i) = dble(pottarg1(1,j,i))
     1              +eye*dble(pottarg1(2,j,i))
            enddo
         enddo
      endif
      if (ifpghtarg .eq. 2 .or. ifpghtarg .eq. 3) then
         do i = 1,nt
            do j = 1,nd
               gradtarg(j,1,i) = dble(gradtarg1(1,j,i))
     1              +eye*dble(gradtarg1(2,j,i))
               gradtarg(j,2,i) = -imag(gradtarg1(1,j,i))
     1              -eye*imag(gradtarg1(2,j,i))
            enddo
         enddo
      endif
      if (ifpghtarg .eq. 3) then
         do i = 1,nt
            do j = 1,nd
               hesstarg(j,1,i) = dble(hesstarg1(1,j,i))
     1              +eye*dble(hesstarg1(2,j,i))
               hesstarg(j,2,i) = -imag(hesstarg1(1,j,i))
     1              -eye*imag(hesstarg1(2,j,i))
               hesstarg(j,3,i) = -hesstarg(j,1,i)
            enddo
         enddo
      endif
      
      return
      end
c
c
