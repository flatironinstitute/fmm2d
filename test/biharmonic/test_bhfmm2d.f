      implicit real *8 (a-h,o-z)
      real *8, allocatable :: sources(:,:),targ(:,:)
      complex *16, allocatable :: charges(:,:),dip(:,:)
      complex *16, allocatable :: pot(:),grad(:,:),hess(:,:)
      complex *16, allocatable :: pottarg(:),gradtarg(:,:)
      complex *16, allocatable :: hesstarg(:,:)

      complex *16, allocatable :: potex(:),gradex(:,:)
      complex *16, allocatable :: hessex(:,:)
      complex *16, allocatable :: pottargex(:),gradtargex(:,:),
     1                             hesstargex(:,:)

      real *8 expc(100),texps(100),scj(100)

      integer ipass(27)
      
      complex *16 ima
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

      ntests=1
      do i = 1,ntests
         ipass(i) = 0
      enddo

      done = 1
      pi = atan(done)*4

      nsrc = 9998
      ntarg = 9999

      allocate(sources(2,nsrc),charges(2,nsrc),dip(3,nsrc))
      allocate(targ(2,ntarg))
      allocate(pot(nsrc),grad(3,nsrc),hess(3,nsrc))
      allocate(pottarg(ntarg),gradtarg(3,ntarg),hesstarg(3,ntarg))

      do i=1,nsrc

         sources(1,i) = hkrand(0)
         sources(2,i) = hkrand(0)

         charges(1,i) = hkrand(0) + ima*hkrand(0) 
         charges(2,i) = hkrand(0) + ima*hkrand(0) 
         dip(1,i) = hkrand(0) + ima*hkrand(0)
         dip(2,i) = hkrand(0) + ima*hkrand(0)
         dip(3,i) = hkrand(0) + ima*hkrand(0)
      enddo


      do i=1,ntarg
         targ(1,i) = hkrand(0)
         targ(2,i) = hkrand(0)
      enddo

      nts = min(20,nsrc)
      ntt = min(20,ntarg)

      allocate(potex(nts),gradex(3,nts),hessex(3,nts))
      allocate(pottargex(ntt),gradtargex(3,ntt),hesstargex(3,ntt))

      eps = 0.5d-6

      ifcharge = 0
      ifdipole = 1
      iper = 0
      ifpgh = 2
      ifpghtarg = 2
      nd = 1
      call bhfmm2d(nd,eps,nsrc,sources,ifcharge,charges,
     1  ifdipole,dip,iper,ifpgh,pot,grad,hess,
     2  ntarg,targ,ifpghtarg,pottarg,gradtarg,
     3  hesstarg,ier)


c
cc       test against exact potential
c
      call dzero(potex,2*nts)
      call dzero(pottargex,2*ntt)
      call dzero(gradex,6*nts)
      call dzero(gradtargex,6*ntt)

      thresh = 1.0d-14

      call bhfmm2dpart_direct(1,1,nsrc,1,nts,sources,ifcharge,
     1       charges,ifdipole,dip,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call bhfmm2dpart_direct(1,1,nsrc,1,ntt,sources,ifcharge,
     1       charges,ifdipole,dip,targ,ifpghtarg,pottargex,
     2       gradtargex,hesstargex,thresh)

     
      errps = 0
      errpt = 0

      errgs = 0
      errgt = 0

      errhs = 0
      errht = 0
      if(ifpgh.ge.1) call derr(potex,pot,2*nts,errps)
      if(ifpghtarg.ge.1) call derr(pottargex,pottarg,2*ntt,errpt)


      if(ifpgh.ge.2) call derr(gradex,grad,6*nts,errgs)
      if(ifpghtarg.ge.2) call derr(gradtargex,gradtarg,6*ntt,errgt)

      call errprintbh(errps,errgs,errhs,errpt,errgt,
     1  errht)
      if (errps .lt. eps .and. errgs .lt. eps 
     1     .and. errhs .lt. eps .and.
     1     errpt .lt. eps .and. errgt .lt. eps 
     1     .and. errht .lt. eps)
     1     ipass(1) = 1
      
      print *, "ipass=",ipass(1)
      ntests = 1

      npass = 0
      do i = 1,ntests
         npass = npass + ipass(i)
      enddo



      open(unit=33,file='print_testres.txt',access='append')
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',npass,
     1   ' out of ',ntests,' tests in bhfmm2d testing suite'
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
      subroutine errprintbh(errps,errgs,errhs,errpt,errgt,errht)
      implicit real *8 (a-h,o-z)
 1100 format(3(2x,e11.5))


      write(6,*) 'error in sources'
      write(6,*) 'pot err, grad err, hess err' 
      write(6,1100) errps,errgs,errhs
      write(6,*) 
      write(6,*)
      write(6,* ) 'error in targets'
      write(6,*) 'pot err, grad err, hess err' 
      write(6,1100) errpt,errgt,errht
      write(6,*)
      write(6,*)
      write(6,*)'==================='

      return
      end
      
