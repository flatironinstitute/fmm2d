cc Copyright (C) 2009-2011: Leslie Greengard and Zydrunas Gimbutas
cc Modified 2019: Leslie Greengard 
cc Contact: greengard@cims.nyu.edu
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
c    $Date: 2011-02-22 17:34:23 -0500 (Tue, 22 Feb 2011) $
c    $Revision: 1670 $
c
c     Computation of  Bessel functions via recurrence
c
c**********************************************************************
      subroutine jbessel2d(nterms,z,rscale,fjs,ifder,fjder)
      implicit real *8 (a-h,o-z)
c**********************************************************************
c
c     PURPOSE:
c
c     This subroutine evaluates the first NTERMS  Bessel 
c     functions and if required, their derivatives.
c     It incorporates a scaling parameter SCALE so that
c       
c	fjs_n(z)= j_n(z)/RSCALE^n
c	fjder_n(z)= \frac{\partial fjs_n(z)}{\partial z}
c
c     NOTE: The scaling parameter RSCALE is meant to be used when
c           abs(z) < 1, in which case we recommend setting
c	    RSCALE = abs(z). This prevents the fjs_n from 
c           underflowing too rapidly. Otherwise, set RSCALE=1.
c	    Do not set RSCALE = abs(z) if z could take on the 
c           value zero. 
c           In an FMM, when forming an expansion from a collection of
c           sources, set RSCALE = min( abs(k*r), 1)
c           where k is the Helmholtz parameter and r is the box dimension
c           at the relevant level.
c
c    INPUT:
c
c    nterms    (integer): order of expansion of output array fjs 
c                         nterms should be greater than 1.
c    z     (complex *16): argument of the  Bessel functions
c    rscale     (real *8): scaling factor (discussed above)
c    ifder     (integer): flag indicating whether to calculate "fjder"
c		          0	NO
c		          1	YES
c    OUTPUT:
c
c    fjs   (complex *16): array of scaled Bessel functions.
c    fjder (complex *16): array of derivs of scaled Bessel functions.
c
c
c    INTERNAL VARIABLES
c    lextra   lextra is buffer large enough to run the 
c             recurrence from nterms upward until
c             the value of fjs has blown up by a factor
c             of 1.0d40 (a parameter hardwired as
c             upbound2 below). lextra should be much larger
c             than nterms. It is set below to 1000
c             If that is not sufficient, an error code is returned (ier=8).
c
c    iscale   integer workspace used to keep track of 
c             internal scaling
c
c
      integer, allocatable :: iscale(:)
      complex *16, allocatable :: fjtemp(:)
      complex *16 wavek,fjs(0:nterms),fjder(0:nterms)
      complex *16 z,zinv,com,zscale,ztmp
      complex *16 psi,zmul,zsn,zmulinv
      complex *16 ima,ffn,ffnm1,ffnp1
      data ima/(0.0d0,1.0d0)/
c
      data upbound/1.0d+32/, upbound2/1.0d+40/, upbound2inv/1.0d-40/
      data tiny/1.0d-200/,done/1.0d0/,zero/0.0d0/
c
c ... Initializing ...
c
      ier=0
ccc      call prinf(' ifder is *',ifder,1)
c
c     set to asymptotic values if argument is sufficiently small
c
      if (abs(z).lt.tiny) then
         fjs(0) = done
         do i = 1, nterms
            fjs(i) = zero
	 enddo
c
	 if (ifder.eq.1) then
	    do i=0,nterms
	       fjder(i)=zero
	    enddo
	    fjder(1)=done/(2*rscale)
	 endif
c
         RETURN
      endif
c
c ... Step 1: carry out upward recursion starting at nterms until the 
c     magnitude has grown by uppound2 = 10^40. The top value of
c     the index is stored as NTOP. That becomes the point from which 
c     the downward recurrence begins.
c
      ntop=nterms+nextra
      zinv=done/z
      ffn = done
      ffnm1 = zero
c
      nextra = 10000
      do 1200 i=nterms,nterms+nextra
         dcoef=2*i
         ztmp=dcoef*zinv*ffn-ffnm1
         ffnp1=ztmp
c
         dd = dreal(ztmp)**2 + dimag(ztmp)**2
         if (dd .gt. upbound2) then
            ntop=i+1
            goto 1300
         endif
         ffnm1 = ffn
         ffn = ffnp1
 1200 continue
 1300 continue
      allocate(iscale(0:ntop))
      allocate(fjtemp(0:ntop))
c
c ... Step 2: Recursion back down to generate the unscaled jfuns:
c             This is the stable but growing direction for the recurrence.
c             To avoid overflow, if the magnitude exceeds UPBOUND2, 
c             we rescale and continue the recursion (saving the order 
c             at which rescaling occurred in array iscale for subsequent 
c             correction).
c
      do i=0,ntop
         iscale(i)=0
      enddo
c
      fjtemp(ntop)=zero
      fjtemp(ntop-1)=done
      do i=ntop-1,1,-1
	 dcoef=2*i
         ztmp=dcoef*zinv*fjtemp(i)-fjtemp(i+1)
         fjtemp(i-1)=ztmp
c
         dd = dreal(ztmp)**2 + dimag(ztmp)**2
         if (dd.gt.upbound2) then
            fjtemp(i) = fjtemp(i)*upbound2inv
            fjtemp(i-1) = fjtemp(i-1)*upbound2inv
            iscale(i) = 1
         endif
      enddo
c
c ...  Step 3: go back up to the top and make sure that all
c              Bessel functions are scaled by the same factor
c              (i.e. the net total of times rescaling was invoked
c              on the way down in the previous loop).
c              At the same time, add scaling to fjs array.
c
      ncntr=0
      scalinv=done/rscale
      sctot = 1.0d0
      do i=1,ntop
         sctot = sctot*scalinv
         if(iscale(i-1).eq.1) sctot=sctot*upbound2inv
         fjtemp(i)=fjtemp(i)*sctot
      enddo
c
c ... Determine the normalization parameter:
c
c     From Abramowitz and Stegun (9.1.47) and (9.1.48), Euler's identity
c
        psi = 0d0
c
        if (dimag(z) .lt. 0) zmul = +ima
        if (dimag(z) .ge. 0) zmul = -ima
        zsn = zmul**(mod(ntop,4))
c
        zmulinv=1/zmul
        do i = ntop,1,-1
           psi = rscale*psi+fjtemp(i)*zsn
           zsn = zsn*zmulinv
        enddo
        psi = 2*psi*rscale+fjtemp(0)
c
        if (dimag(z) .lt. 0) zscale = cdexp(+ima*z) / psi
        if (dimag(z) .ge. 0) zscale = cdexp(-ima*z) / psi
c
c
c ... Scale the jfuns by zscale:
c
      ztmp=zscale
      do i=0,nterms
         fjs(i)=fjtemp(i)*ztmp
      enddo
c
c ... Finally, calculate the derivatives if desired:
c
      if (ifder.eq.1) then
         fjs(nterms+1)=fjtemp(nterms+1)*ztmp
c
         fjder(0)=-fjs(1)*rscale
         dc1=0.5d0
         dc2=done-dc1
         dc1=dc1*scalinv
         dc2=dc2*rscale
         do i=1,nterms
            fjder(i)=dc1*fjs(i-1)-dc2*fjs(i+1)
         enddo
      endif
      return
      end
c
