      implicit real *8 (a-h,o-z)
      real *8, allocatable :: sources(:,:),targ(:,:)
      complex *16, allocatable :: charges(:,:),dipstr(:,:)
      complex *16, allocatable :: pot(:,:),grad(:,:),hess(:,:)
      complex *16, allocatable :: pottarg(:,:),gradtarg(:,:),
     1    hesstarg(:,:)
      real *8, allocatable :: potr(:,:),pottargr(:,:)
      real *8, allocatable :: potexr(:,:),pottargexr(:,:)

      complex *16, allocatable :: potex(:,:),gradex(:,:),
     1   hessex(:,:)
      complex *16, allocatable :: pottargex(:,:),gradtargex(:,:),
     1                             hesstargex(:,:)

      real *8 expc(100),texps(100),scj(100)

      integer ipass(27)
      
      complex *16 ima
      integer nd,idim
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

      ntests=9
      do i = 1,ntests
         ipass(i) = 0
      enddo
      
      done = 1
      pi = atan(done)*4


      nd = 2


      nsrc = 10000
      ntarg = nsrc

      allocate(sources(2,nsrc),charges(nd,nsrc),dipstr(nd,nsrc))
      allocate(targ(2,ntarg))
      allocate(pot(nd,nsrc),grad(nd,nsrc),hess(nd,nsrc))
      allocate(pottarg(nd,ntarg),gradtarg(nd,ntarg),
     1   hesstarg(nd,ntarg))

      rin = 1.0d0
      rwig = 0.1d0
      nwig = 10
      do i=1,nsrc

         thet = hkrand(0)*2*pi

cc         sources(1,i) = (rin + rwig*cos(thet))*cos(thet)
cc         sources(2,i) = (rin + rwig*cos(thet))*sin(thet)

         sources(1,i) = hkrand(0)
         sources(2,i) = hkrand(0)
         do idim=1,nd

           charges(idim,i) = hkrand(0) 
           dipstr(idim,i) = hkrand(0) + ima*hkrand(0)
        enddo
      enddo


      do i=1,ntarg
         targ(1,i) = hkrand(0)
         targ(2,i) = hkrand(0)
      enddo

      nts = min(20,nsrc)
      ntt = min(20,ntarg)

      allocate(potex(nd,nts),gradex(nd,nts),hessex(nd,nts))
      allocate(pottargex(nd,ntt),gradtargex(nd,ntt),
     1   hesstargex(nd,ntt))
      allocate(potr(nd,nts),potexr(nd,nts))
      allocate(pottargr(nd,ntt),pottargexr(nd,ntt))


c
cc      test low frequency mode
c
      eps = 0.6d-5
c
c
cc      now test source to source  + target, charge
c       with potentials
c
      write(*,*) 'testing stost, charge, potentials'

      call cfmm2d_st_c_p_vec(nd,eps,nsrc,sources,charges,
     1        pot,ntarg,targ,pottarg,ier)


c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(pottargex,2*ntt*nd)

      ifcharge = 1
      ifdipole = 0

      ifpgh = 1
      thresh = 1.0d-14

      call cfmm2dpart_direct(nd,1,nsrc,1,nts,sources,ifcharge,
     1       charges,ifdipole,dipstr,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call cfmm2dpart_direct(nd,1,nsrc,1,ntt,sources,ifcharge,
     1       charges,ifdipole,dipstr,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      do i=1,nts
        do idim=1,nd
          potexr(idim,i) = real(potex(idim,i))
          potr(idim,i) = real(pot(idim,i))
        enddo
      enddo

      do i=1,ntt
        do idim=1,nd
          pottargexr(idim,i) = real(pottargex(idim,i))
          pottargr(idim,i) = real(pottarg(idim,i))
        enddo
      enddo

      call derr(potexr,potr,nts*nd,errps)
      call derr(pottargexr,pottargr,ntt*nd,errpt)

      errgs = 0
      errgt = 0

      errhs = 0
      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht)
      if (errps .lt. eps .and. errgs .lt. eps .and. errhs .lt. eps .and.
     1     errpt .lt. eps .and. errgt .lt. eps .and. errht .lt. eps)
     1     ipass(1) = 1


c
c
cc      now test source to source  + target, charge
c       with gradients
c
      write(*,*) 'testing stost, charge, gradients'

      call cfmm2d_st_c_g_vec(nd,eps,nsrc,sources,charges,
     1        pot,grad,ntarg,targ,pottarg,gradtarg,ier)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(gradex,2*nts*nd)
      call dzero(pottargex,2*ntt*nd)
      call dzero(gradtargex,2*ntt*nd)

      ifcharge = 1
      ifdipole = 0

      ifpgh = 2
      thresh = 1.0d-14

      call cfmm2dpart_direct(nd,1,nsrc,1,nts,sources,ifcharge,
     1       charges,ifdipole,dipstr,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call cfmm2dpart_direct(nd,1,nsrc,1,ntt,sources,ifcharge,
     1       charges,ifdipole,dipstr,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      do i=1,nts
        do idim=1,nd
          potexr(idim,i) = real(potex(idim,i))
          potr(idim,i) = real(pot(idim,i))
        enddo
      enddo

      do i=1,ntt
        do idim=1,nd
          pottargexr(idim,i) = real(pottargex(idim,i))
          pottargr(idim,i) = real(pottarg(idim,i))
        enddo
      enddo

      call derr(potexr,potr,nts*nd,errps)
      call derr(pottargexr,pottargr,ntt*nd,errps)
      
      call derr(gradex,grad,2*nts*nd,errgs)
      call derr(gradtargex,gradtarg,2*ntt*nd,errgt)

      errhs = 0
      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht)
      if (errps .lt. eps .and. errgs .lt. eps .and. errhs .lt. eps .and.
     1     errpt .lt. eps .and. errgt .lt. eps .and. errht .lt. eps)
     1     ipass(2) = 1


c
c
cc      now test source to source  + target, charge
c       with hessians
c



      write(*,*) 'testing stost, charge, hessians'
      call cfmm2d_st_c_h_vec(nd,eps,nsrc,sources,charges,
     1        pot,grad,hess,ntarg,targ,pottarg,gradtarg,
     2        hesstarg,ier)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(gradex,2*nts*nd)
      call dzero(hessex,2*nts*nd)
      call dzero(pottargex,2*ntt*nd)
      call dzero(gradtargex,2*ntt*nd)
      call dzero(hesstargex,2*ntt*nd)

      ifcharge = 1
      ifdipole = 0

      ifpgh = 3
      thresh = 1.0d-14

      call cfmm2dpart_direct(nd,1,nsrc,1,nts,sources,ifcharge,
     1       charges,ifdipole,dipstr,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call cfmm2dpart_direct(nd,1,nsrc,1,ntt,sources,ifcharge,
     1       charges,ifdipole,dipstr,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      do i=1,nts
        do idim=1,nd
          potexr(idim,i) = real(potex(idim,i))
          potr(idim,i) = real(pot(idim,i))
        enddo
      enddo

      do i=1,ntt
        do idim=1,nd
          pottargexr(idim,i) = real(pottargex(idim,i))
          pottargr(idim,i) = real(pottarg(idim,i))
        enddo
      enddo

      call derr(potexr,potr,nts*nd,errps)
      call derr(pottargexr,pottargr,ntt*nd,errps)
      
      call derr(gradex,grad,2*nts*nd,errgs)
      call derr(hessex,hess,2*nts*nd,errhs)
      call derr(gradtargex,gradtarg,2*ntt*nd,errgt)
      call derr(hesstargex,hesstarg,2*nts*nd,errht)

      call errprint(errps,errgs,errhs,errpt,errgt,errht)
      if (errps .lt. eps .and. errgs .lt. eps .and. errhs .lt. eps .and.
     1     errpt .lt. eps .and. errgt .lt. eps .and. errht .lt. eps)
     1     ipass(3) = 1



c
c
c
c
cc      now test source to source  + target, dipole
c       with potentials
c
      write(*,*) 'testing stost, dipole, potentials'

      call cfmm2d_st_d_p_vec(nd,eps,nsrc,sources,dipstr,
     1        pot,ntarg,targ,pottarg,ier)


c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(pottargex,2*ntt*nd)

      ifcharge = 0
      ifdipole = 1

      ifpgh = 1
      thresh = 1.0d-14

      call cfmm2dpart_direct(nd,1,nsrc,1,nts,sources,ifcharge,
     1       charges,ifdipole,dipstr,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call cfmm2dpart_direct(nd,1,nsrc,1,ntt,sources,ifcharge,
     1       charges,ifdipole,dipstr,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      do i=1,nts
        do idim=1,nd
          potexr(idim,i) = real(potex(idim,i))
          potr(idim,i) = real(pot(idim,i))
        enddo
      enddo

      do i=1,ntt
        do idim=1,nd
          pottargexr(idim,i) = real(pottargex(idim,i))
          pottargr(idim,i) = real(pottarg(idim,i))
        enddo
      enddo

      call derr(potexr,potr,nts*nd,errps)
      call derr(pottargexr,pottargr,ntt*nd,errps)
      
  
      errgs = 0
      errgt = 0

      errhs = 0
      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht)
      if (errps .lt. eps .and. errgs .lt. eps .and. errhs .lt. eps .and.
     1     errpt .lt. eps .and. errgt .lt. eps .and. errht .lt. eps)
     1     ipass(4) = 1


c
c
cc      now test source to source  + target, dipole
c       with gradients
c
      write(*,*) 'testing stost, dipole, gradients'

      call cfmm2d_st_d_g_vec(nd,eps,nsrc,sources,dipstr,
     1        pot,grad,ntarg,targ,pottarg,gradtarg,ier)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(gradex,2*nts*nd)
      call dzero(pottargex,2*ntt*nd)
      call dzero(gradtargex,2*ntt*nd)

      ifcharge = 0
      ifdipole = 1

      ifpgh = 2
      thresh = 1.0d-14

      call cfmm2dpart_direct(nd,1,nsrc,1,nts,sources,ifcharge,
     1       charges,ifdipole,dipstr,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call cfmm2dpart_direct(nd,1,nsrc,1,ntt,sources,ifcharge,
     1       charges,ifdipole,dipstr,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      do i=1,nts
        do idim=1,nd
          potexr(idim,i) = real(potex(idim,i))
          potr(idim,i) = real(pot(idim,i))
        enddo
      enddo

      do i=1,ntt
        do idim=1,nd
          pottargexr(idim,i) = real(pottargex(idim,i))
          pottargr(idim,i) = real(pottarg(idim,i))
        enddo
      enddo

      call derr(potexr,potr,nts*nd,errps)
      call derr(pottargexr,pottargr,ntt*nd,errps)
      
      call derr(gradex,grad,2*nts*nd,errgs)
      call derr(gradtargex,gradtarg,2*ntt*nd,errgt)

      errhs = 0
      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht)
      if (errps .lt. eps .and. errgs .lt. eps .and. errhs .lt. eps .and.
     1     errpt .lt. eps .and. errgt .lt. eps .and. errht .lt. eps)
     1     ipass(5) = 1


c
c
cc      now test source to source  + target, dipole
c       with hessians
c



      write(*,*) 'testing stost, dipole, hessians'
      call cfmm2d_st_d_h_vec(nd,eps,nsrc,sources,dipstr,
     1        pot,grad,hess,ntarg,targ,pottarg,gradtarg,
     2        hesstarg,ier)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(gradex,2*nts*nd)
      call dzero(hessex,2*nts*nd)
      call dzero(pottargex,2*ntt*nd)
      call dzero(gradtargex,2*ntt*nd)
      call dzero(hesstargex,2*ntt*nd)

      ifcharge = 0
      ifdipole = 1

      ifpgh = 3
      thresh = 1.0d-14

      call cfmm2dpart_direct(nd,1,nsrc,1,nts,sources,ifcharge,
     1       charges,ifdipole,dipstr,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call cfmm2dpart_direct(nd,1,nsrc,1,ntt,sources,ifcharge,
     1       charges,ifdipole,dipstr,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      do i=1,nts
        do idim=1,nd
          potexr(idim,i) = real(potex(idim,i))
          potr(idim,i) = real(pot(idim,i))
        enddo
      enddo

      do i=1,ntt
        do idim=1,nd
          pottargexr(idim,i) = real(pottargex(idim,i))
          pottargr(idim,i) = real(pottarg(idim,i))
        enddo
      enddo

      call derr(potexr,potr,nts*nd,errps)
      call derr(pottargexr,pottargr,ntt*nd,errps)
      
      call derr(gradex,grad,2*nts*nd,errgs)
      call derr(hessex,hess,2*nts*nd,errhs)
      call derr(gradtargex,gradtarg,2*ntt*nd,errgt)
      call derr(hesstargex,hesstarg,2*nts*nd,errht)

      call errprint(errps,errgs,errhs,errpt,errgt,errht)
      if (errps .lt. eps .and. errgs .lt. eps .and. errhs .lt. eps .and.
     1     errpt .lt. eps .and. errgt .lt. eps .and. errht .lt. eps)
     1     ipass(6) = 1

c
c
c
c
cc      now test source to source  + target, charge + dipole
c       with potentials
c
      write(*,*) 'testing stost, charge + dipole, potentials'

      call cfmm2d_st_cd_p_vec(nd,eps,nsrc,sources,charges,
     1     dipstr,pot,ntarg,targ,pottarg,ier)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(pottargex,2*ntt*nd)

      ifcharge = 1
      ifdipole = 1

      ifpgh = 1
      thresh = 1.0d-14

      call cfmm2dpart_direct(nd,1,nsrc,1,nts,sources,ifcharge,
     1       charges,ifdipole,dipstr,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call cfmm2dpart_direct(nd,1,nsrc,1,ntt,sources,ifcharge,
     1       charges,ifdipole,dipstr,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      do i=1,nts
        do idim=1,nd
          potexr(idim,i) = real(potex(idim,i))
          potr(idim,i) = real(pot(idim,i))
        enddo
      enddo

      do i=1,ntt
        do idim=1,nd
          pottargexr(idim,i) = real(pottargex(idim,i))
          pottargr(idim,i) = real(pottarg(idim,i))
        enddo
      enddo

      call derr(potexr,potr,nts*nd,errps)
      call derr(pottargexr,pottargr,ntt*nd,errps)
      
  
      errgs = 0
      errgt = 0

      errhs = 0
      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht)
      if (errps .lt. eps .and. errgs .lt. eps .and. errhs .lt. eps .and.
     1     errpt .lt. eps .and. errgt .lt. eps .and. errht .lt. eps)
     1     ipass(7) = 1


c
c
cc      now test source to source  + target, charge + dipole
c       with gradients
c
      write(*,*) 'testing stost, charge + dipole, gradients'

      call cfmm2d_st_cd_g_vec(nd,eps,nsrc,sources,charges,
     1    dipstr,pot,grad,ntarg,targ,pottarg,gradtarg,ier)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(gradex,2*nts*nd)
      call dzero(pottargex,2*ntt*nd)
      call dzero(gradtargex,2*ntt*nd)

      ifcharge = 1
      ifdipole = 1

      ifpgh = 2
      thresh = 1.0d-14

      call cfmm2dpart_direct(nd,1,nsrc,1,nts,sources,ifcharge,
     1       charges,ifdipole,dipstr,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call cfmm2dpart_direct(nd,1,nsrc,1,ntt,sources,ifcharge,
     1       charges,ifdipole,dipstr,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      do i=1,nts
        do idim=1,nd
          potexr(idim,i) = real(potex(idim,i))
          potr(idim,i) = real(pot(idim,i))
        enddo
      enddo

      do i=1,ntt
        do idim=1,nd
          pottargexr(idim,i) = real(pottargex(idim,i))
          pottargr(idim,i) = real(pottarg(idim,i))
        enddo
      enddo

      call derr(potexr,potr,nts*nd,errps)
      call derr(pottargexr,pottargr,ntt*nd,errps)
      
      call derr(gradex,grad,2*nts*nd,errgs)
      call derr(gradtargex,gradtarg,2*ntt*nd,errgt)

      errhs = 0
      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht)
      if (errps .lt. eps .and. errgs .lt. eps .and. errhs .lt. eps .and.
     1     errpt .lt. eps .and. errgt .lt. eps .and. errht .lt. eps)
     1     ipass(8) = 1


c
c
cc      now test source to source  + target, charge + dipole
c       with hessians
c



      write(*,*) 'testing stost, charge + dipole, hessians'
      call cfmm2d_st_cd_h_vec(nd,eps,nsrc,sources,charges,
     1    dipstr,pot,grad,hess,ntarg,targ,pottarg,gradtarg,
     2        hesstarg,ier)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(gradex,2*nts*nd)
      call dzero(hessex,2*nts*nd)
      call dzero(pottargex,2*ntt*nd)
      call dzero(gradtargex,2*ntt*nd)
      call dzero(hesstargex,2*ntt*nd)

      ifcharge = 1
      ifdipole = 1

      ifpgh = 3
      thresh = 1.0d-14

      call cfmm2dpart_direct(nd,1,nsrc,1,nts,sources,ifcharge,
     1       charges,ifdipole,dipstr,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call cfmm2dpart_direct(nd,1,nsrc,1,ntt,sources,ifcharge,
     1       charges,ifdipole,dipstr,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      do i=1,nts
        do idim=1,nd
          potexr(idim,i) = real(potex(idim,i))
          potr(idim,i) = real(pot(idim,i))
        enddo
      enddo

      do i=1,ntt
        do idim=1,nd
          pottargexr(idim,i) = real(pottargex(idim,i))
          pottargr(idim,i) = real(pottarg(idim,i))
        enddo
      enddo

      call derr(potexr,potr,nts*nd,errps)
      call derr(pottargexr,pottargr,ntt*nd,errps)
      
      call derr(gradex,grad,2*nts*nd,errgs)
      call derr(hessex,hess,2*nts*nd,errhs)
      call derr(gradtargex,gradtarg,2*ntt*nd,errgt)
      call derr(hesstargex,hesstarg,2*nts*nd,errht)

      call errprint(errps,errgs,errhs,errpt,errgt,errht)
      if (errps .lt. eps .and. errgs .lt. eps .and. errhs .lt. eps .and.
     1     errpt .lt. eps .and. errgt .lt. eps .and. errht .lt. eps)
     1     ipass(9) = 1

      npass = 0
      do i = 1,ntests
         npass = npass + ipass(i)
      enddo

      open(unit=33,file='print_testres.txt',access='append')
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',npass,
     1     ' out of ',ntests,' tests in cfmm2d vec testing suite'
      close(33)
      
      stop
      end
c-----------------------------------------------------     
      subroutine dzero(vec,n)
      implicit real *8 (a-h,o-z)
      real *8 vec(*)

      do i=1,n
         vec(i) = 0
      enddo

      return
      end
c------------------------------------
      subroutine derr(vec1,vec2,n,erra)
      implicit real *8 (a-h,o-z)
      real *8 vec1(*),vec2(*)

      ra = 0
      erra = 0
      do i=1,n
         ra = ra + vec1(i)**2
         erra = erra + (vec1(i)-vec2(i))**2
      enddo

      erra = sqrt(erra/ra)

      return
      end
c----------------------------------
      subroutine errprint(errps,errgs,errhs,errpt,errgt,errht)
      implicit real *8 (a-h,o-z)
 1100 format(3(2x,e11.5))


      write(*,*) 'error in sources'
      write(*,*) 'pot err, grad err, hess err' 
      write(*,1100) errps,errgs,errhs
      write(*,*) 
      write(*,*)
      write(*,* ) 'error in targets'
      write(*,*) 'pot err, grad err, hess err' 
      write(*,1100) errpt,errgt,errht
      write(*,*)
      write(*,*)
      write(*,*)'==================='

      return
      end
      
