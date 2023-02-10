c
c  This file has the following user callable subroutines:
c
c  
c
c
c  NOTES: even though as operators the mpmphf and loclochf 
c  should be identical. Owing to openmp considerations 
c  and minimizing the number of mp/loc_to_sig or sig_to_exp
C  operations, we have split out certain parts of the routines
C  differently.
C
C  In particular, for the mpmphf, we act on multipole expansions
C  and return the signature at the parent, 
C  while for the loclochf, we act on the signature of the parent
C  and return the local expansion at the children.
c
c
c
c
C***********************************************************************
      subroutine h2dloclochf(nd,zk,rscale1,center1,sig,nterms1,nsig,
     1                rscale2,center2,hexp2,
     2                nterms2,transvec,wsave)
      implicit none
C***********************************************************************
c      
c     This routine shifts a signature expansion SIG to a new center,
C     generating the multipole expansions HEXP2, using FFT convolution.
C     THIS IS AN INCREMENTING ROUTINE - hexp2 is incremented, not
C     overwritten. In the upward pass, this allows accumulation of data 
C     from children. In the downward pass, this allows hexp2 to be 
C     incremented, assuming other data has come in through 
C     <<list 3/4>> considerations.
C
C---------------------------------------------------------------------
c      INPUT:
c      
c      nd      : vector length (number of expansions)
c      zk      : Helmholtz parameter (ignored in opt. version here)
c      rscale1 : scaling parameter (IGNORED and assumed equal to 1)
c      center1 : center of original multiple expansion
c      sig     : coefficients of original signature expansion
c      nterms1 : order of original multipole expansion
c      rscale2 : scaling parameter (IGNORED and assumed equal to 1)
c      center2 : center of shifted multipole expansion
c      nterms2 : order of shifted multipole expansion
c      transvec: precomputed translation operator
c      wsave   : fftpack precomputed array - needs threadsafe dfft if
c                using openmp.
C---------------------------------------------------------------------
c      OUTPUT:
c      
c      hexp2   = coefficients of shifted multipole expansion
c      
c      Note:
c      The FFT accelerated convolution is only valid at high
c      frequencies, and therefore rscale1 and rscale2 parameters
c      are assumed to be 1 (i.e. ignored) inside this routine.
C---------------------------------------------------------------------
      integer nterms1,nterms2,nsig,j,next235,nd,ii
      double precision center1(2),center2(2),dn,rscale1,rscale2
      double complex zk
      double complex hexp2(nd,-nterms2:nterms2)
      double complex sig(nd,nsig)
      double complex, allocatable :: sig2(:,:)
      double complex transvec(*)
      double complex wsave(*)
c     
c     do the convolution via fft
c     
      allocate(sig2(nd,nsig))
c
      do ii = 1,nd
      do j = 1,nsig
         sig2(ii,j) = 0.0d0
      enddo
      enddo
c
      call h2d_diagtrans(nd,nsig,sig,transvec,sig2)
      call h2d_sig2exp(nd,nsig, sig2, wsave, nterms2, hexp2)
      return
      end 
c
c
c      
C***********************************************************************
      subroutine h2dmpmphf(nd,zk,rscale1,center1,hexp1,nterms1,rscale2,
     1                center2,sig2,nterms2,nsig,wsave,transvec)
      implicit none
C***********************************************************************
c      
c     This routine shifts multipole expansions HEXP1 to a new center,
C     generating the multipole expansions HEXP2, using FFT convolution.
C
C     The subroutine h2d_diagtrans(nsig,sig,transvec,sig2)
C     IS INCREMENTAL - that is sig2 is incremented by the translation of
C     sig.
C     
C     Note that this subroutine returns the signature instead of
C     updating the multipole expansion. 
C
C     The accumulated signature is then converted to a multipole
C     expansion in the outer loop
C---------------------------------------------------------------------
c      INPUT:
c      
c      nd      : vector length (number of expansions)
c      zk      : Helmholtz parameter
c      rscale1 : scaling parameter (IGNORED and assumed equal to 1)
c      center1 : center of original multiple expansion
c      hexp1   : coefficients of original multiple expansion
c      nterms1 : order of original multipole expansion
c      rscale2 : scaling parameter (IGNORED and assumed equal to 1)
c      center2 : center of shifted multipole expansion
c      nterms2 : order of shifted multipole expansion
C---------------------------------------------------------------------
c      OUTPUT:
c      
c      sig2   = coefficients of shifted signature expansion
c      
c      Note:
c      The FFT accelerated convolution is only valid at high
c      frequencies, and therefore rscale1 and rscale2 parameters
c      are assumed to be 1 (i.e. ignored) inside this routine.
C---------------------------------------------------------------------
      integer nterms1,nterms2,nsig,j,next235,nd,ii
      double precision center1(2),center2(2),dn,rscale1,rscale2
      double complex zk,hexp1(nd,-nterms1:nterms1) 
      double complex wsave(4*nsig+100)
      double complex sig2(nd,nsig)
      double complex, allocatable :: sig(:,:)
      double complex transvec(nsig)
c     
c     do the convolution via fft
c     
      allocate(sig(nd,nsig))
c
      call h2d_mptosig(nd,nterms1,nsig,hexp1,sig,wsave)
      call h2d_diagtrans(nd,nsig,sig,transvec,sig2)

      return
      end 
c
c
c
c
C***********************************************************************
      subroutine h2d_mptosig(nd,nterms1,nsig,hexp,sig,wsave)
      implicit none 
C***********************************************************************
c
c     Converts h-expansion HEXP to far-field  high-frequency
c     signature function SIG (using FFT).
C---------------------------------------------------------------------
c     INPUT:
c
c     nd      : vector length (number of expansions)
c     nterms1 : order of original local expansion
c     nsig    : order of the HF signature
c     hexp    : coefficients of original local expansion
c     wsave   : FFT work array from call to ZFFTI with size nsig
c     
C---------------------------------------------------------------------
c     OUTPUT:
c     
c     sig     : far-field signature function, stored FFT style.
C---------------------------------------------------------------------
      integer :: nterms1, nsig,j,nd,ii
      double complex :: hexp(nd,-nterms1:nterms1),sig(nd,nsig),wsave(*)
      double complex, allocatable :: sigtmp(:)
c      
c     insert hexp to sig in FFT ordering
c
      allocate(sigtmp(nsig))
c
      do ii = 1,nd
         do j = 1,nsig
            sigtmp(j) = 0
         end do
         do j = 0,nterms1
            sigtmp(j+1) = hexp(ii,j)
         end do
         do j = 1,nterms1
            sigtmp(nsig-j+1) = hexp(ii,-j)
         end do
c
         call zfftf(nsig, sigtmp, wsave)
         do j = 1,nsig
            sig(ii,j) = sigtmp(j)
         end do
      end do

      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine h2d_mkmpshift(zk, center1, nterms1,
     1           center2, nterms2, nsig, wsave, transvec)
      implicit none
C***********************************************************************
c
c     Create diagonal translation operator for shifting center of either
c     h- or j- expansion.
C---------------------------------------------------------------------
c     INPUT:
c
c     zk       : Helmholtz parameter
c     center1  : center of original local expansion
c     nterms1  : order of original local expansion
c     center2  : center of shifted local expansion
c     nterms2  : order of the eventual shifted local expansion
c     nsig     : order of the HF signature, should be good for FFT and
c                greater than 2*(nterm1+nterms2)+1
c     wsave    : FFT work array from call to ZFFTI with size nsig
C---------------------------------------------------------------------
c     OUTPUT:
c     
c     transvec : translation operator in diagonal form
C---------------------------------------------------------------------
      integer :: nterms1,nterms2,nterms,nsig,ifder,j
      double precision :: center1(2),center2(2)
      double precision :: zdiff(2),done,pi,theta,r,rscale1
      double complex :: zk,transvec(nsig),wsave(*)
      double complex :: z,ima,zmul,zinv,ztemp1,ztemp2
      double complex, allocatable :: jval(:),jder(:),jtemp(:)
      data ima/(0.0d0,1.0d0)/
c
      done=1
      pi=4*atan(done)
      nterms = nterms1+nterms2
c
      allocate(jval(0:nterms+3))
      allocate(jder(0:nterms+3))
      allocate(jtemp(-nterms-3:nterms+3))
c
      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)
      theta=theta-pi
      z=zk*r
c
      ifder=0
      rscale1 = 1
      call jbessel2d(nterms,z,rscale1,jval,ifder,jder)
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
c     compute translation operator by FFT
c
      do j = 1,nsig
      transvec(j) = 0
      end do
      do j = 0,nterms
      transvec(j+1) = jtemp(j)/nsig
      end do
      do j = 1,nterms
      transvec(nsig-j+1) = jtemp(-j)/nsig
      end do
      call zfftf(nsig, transvec, wsave)
c
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine h2d_diagtrans(nd,nsig,sig,transvec,sig2)
      implicit none
C***********************************************************************
c
c     Carry out vector translation operator in diagonal form and 
c     INCREMENT SIG2.
c     
C---------------------------------------------------------------------
c     INPUT:
c
c     nd       : vector length (number of expansions)
c     nsig     : order of the HF signature,
c     sig      : original far-field signature function
c     transvec : diagonal translation operator
C---------------------------------------------------------------------
C     OUTPUT:
c     
c     sig2     : translated expansion in diagonal form
C---------------------------------------------------------------------
      integer :: nsig,i,nd,ii
      double complex :: sig(nd,nsig), transvec(nsig), sig2(nd,nsig)
c      
      do i = 1,nsig
         do ii = 1,nd
            sig2(ii,i) = sig2(ii,i) + transvec(i)*sig(ii,i)
         end do
      end do
c
      return
      end
c
c
c
c
c
C***********************************************************************
      subroutine h2dmplochf(nd,zk,rscale1,center1,sig,nterms1,
     1                rscale2,center2,sig2,nterms2,nsig,wsave,tvec)
      implicit none
C***********************************************************************
c      
c     This routine maps multipole expansions HEXP to local expansions
C     using FFT convolution.
C     THIS IS NOT AN INCREMENTING ROUTINE - ONLY A TESTING ROUTINE.
C
C     The subroutine h2d_mploctrans(nsig,sig,transvec,sig2)
C     IS INCREMENTAL - that is sig2 is incremented by the translation of
C     sig. Inside the FMM, sig2 accumulates the contributions from all
C     incoming expansions before transformation back to the "physical" jexp
C     form.
C---------------------------------------------------------------------
c      INPUT:
c      
c      nd      : vector length (number of expansions)
c      zk      : Helmholtz parameter
c      rscale1 : scaling parameter (IGNORED and assumed equal to 1)
c      center1 : center of original multiple expansion
c      sig     : outgoing expansion in diag form
c      nterms1 : order of original multipole expansion
c      rscale2 : scaling parameter (IGNORED and assumed equal to 1)
c      center2 : center of shifted multipole expansion
c      nterms2 : order of shifted multipole expansion
C---------------------------------------------------------------------
c      OUTPUT:
c      
c      sig2    : increments coefficients of incoming signature expansions
c      
c      Note:
c      The FFT accelerated convolution is only valid at high
c      frequencies, and therefore rscale1 and rscale2 parameters
c      are assumed to be 1 (i.e. ignored) inside this routine.
C---------------------------------------------------------------------
      integer nterms1,nterms2,nsig,j,next235,nd,ii,npts
      double precision center1(2),center2(2),dn,rscale1,rscale2
      double complex zk
ccc      double complex wsave(4*nsig+100)
      double complex wsave(*)
      double complex sig(nd,nsig)
      double complex sig2(nd,nsig)
ccc      double complex, allocatable :: sig(:,:)
ccc      double complex, allocatable :: transvec(:)
      double complex tvec(nsig)
ccc      double complex, allocatable :: sig2(:,:)
c     
c     do the convolution via fft
c     
ccc      dn = 2*(nterms1+nterms2) + 1
ccc      nsig = next235(dn)
c
ccc      allocate(sig(nd,nsig))
ccc      allocate(transvec(nsig))
ccc      allocate(sig2(nd,nsig))
ccc      allocate(wsave(4*nsig+100))
c
ccc      call zffti(nsig,wsave)
ccc      call h2d_mptosig(nd,nterms1,nsig,hexp,sig,wsave)
ccc      do ii = 1,nsig
ccc         transvec(ii) = 0.1d0
ccc      enddo
ccc      call h2d_mkm2ltrans(zk, center1, nterms1,
ccc     1           center2, nterms2, nsig, wsave, transvec)
c
ccc      do ii = 1,nd
ccc      do j = 1,nsig
ccc         sig2(ii,j) = 0.0d0
ccc      enddo
ccc      enddo
c
ccc      call h2d_diagtrans(nd,nsig,sig,transvec,sig2)
      call h2d_diagtrans(nd,nsig,sig,tvec,sig2)
c
ccc      call h2d_sig2exp(nd,nsig, sig2, wsave, nterms2, jexp)
      return
      end 
c
c
c
c
c
C***********************************************************************
      subroutine h2d_sig2exp(nd,nsig,sig,wsave,nterms,expans)
      implicit none
C***********************************************************************
c     
c     This routine converts high-frequency signature functions to
c     either h-expansions or j-expansions (since it's the same formula
c     for both).
c     
C---------------------------------------------------------------------
c     INPUT:
c
c     nd     : vector length (number of expansions)
c     nsig   : order of signature function
c     sig    : values of signature function, ordered FFT style
c     wsave  : work array constructed via a previous call to ZFFTI
c     
C---------------------------------------------------------------------
c     OUTPUT:
c
c     nterms : desired order of expansion
c     expans : expansion coefficients, ordered multipole style
C---------------------------------------------------------------------
        integer :: nd,nsig, nterms,ii,j
        double complex :: sig(nd,nsig), wsave(nsig)
        double complex :: expans(nd,-nterms:nterms)
        double complex, allocatable :: sig2(:)
c     
      allocate(sig2(nsig))
c
      do ii = 1,nd
         do j = 1,nsig
           sig2(j) = sig(ii,j)
         end do
         call zfftb(nsig, sig2, wsave)

         do j = 0,nterms
            expans(ii,j) = expans(ii,j) + sig2(j+1)
         end do
         do j = 1,nterms
            expans(ii,-j) = expans(ii,-j) + sig2(nsig-j+1)
         end do
      enddo
  
      return
      end 
c
c
c
c
c
C***********************************************************************
      subroutine h2d_mkm2ltrans(zk, center1, nterms1,
     1           center2, nterms2, nsig, wsave, transvec)
C***********************************************************************
      implicit none
c      
c
c     Create diagonal translation operator for mapping h-expansion to
c     j-expansion.
C---------------------------------------------------------------------
c      INPUT:
c
c      zk      : Helmholtz parameter
c      center1 : center of original multipole expansion
c      nterms1 : order of original multipole expansion
c      center2 : center of shifted local expansion
c      nterms2 : order of the eventual shifted local expansion
c      nsig    : order of the HF signature, should be good for FFT and
c                greater than 2*(nterm1+nterms2)+1
c      wsave   : FFT work array from call to ZFFTI with size nsig
C---------------------------------------------------------------------
c      OUTPUT:
c      
c     transvec : translation operator in diagonal form
C---------------------------------------------------------------------
      integer :: nterms1,nterms2,nsig,nterms,ifder,j
      double complex :: zk, transvec(nsig)
      double complex :: wsave(*)
      double precision :: center1(2), center2(2), zdiff(2)
      double precision :: done,pi,r,theta,rscale1
      double complex :: z,ima, zmul,zinv,ztemp1,ztemp2
      double complex, allocatable :: hval(:), hder(:), htemp(:)

      data ima/(0.0d0,1.0d0)/

      done=1
      pi=4*atan(done)

      nterms = nterms1+nterms2
c
      allocate(hval(0:nterms+5))
      allocate(hder(0:nterms+5))
      allocate(htemp(-nterms-5:nterms+5))

      zdiff(1)=center2(1)-center1(1)
      zdiff(2)=center2(2)-center1(2)
      call h2cart2polar(zdiff,r,theta)

      theta=theta-pi
      z=zk*r
c
      ifder=0
      rscale1 = 1
      call h2dall(nterms+1,z,rscale1,hval,ifder,hder)

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
c     compute transfer function by FFT
c     
      do j = 1,nsig
        transvec(j) = 0
      end do
      do j = 0,nterms
        transvec(j+1) = htemp(j)/nsig
      end do
      do j = 1,nterms
        transvec(nsig-j+1) = htemp(-j)/nsig
      end do
c
      call zfftf(nsig, transvec, wsave)
c
      return
      end 

