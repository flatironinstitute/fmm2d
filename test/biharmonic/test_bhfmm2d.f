      implicit real *8 (a-h,o-z)
      real *8, allocatable :: sources(:,:),targ(:,:)
      complex *16, allocatable :: charges(:),dip1(:),dip2(:)
      complex *16, allocatable :: pot(:),grada(:),gradaa(:),hess(:,:)
      complex *16, allocatable :: pottarg(:),gradatarg(:),gradaatarg(:)
      complex *16, allocatable :: hesstarg(:,:)

      complex *16, allocatable :: potex(:),gradaex(:),gradaaex(:)
      complex *16, allocatable :: hessex(:,:)
      complex *16, allocatable :: pottargex(:),gradatargex(:),
     1                             gradaatargex(:),hesstargex(:,:)

      real *8 expc(100),texps(100),scj(100)

      integer ipass(27)
      
      complex *16 ima
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

      ntests=9
      do i = 1,ntests
         ipass(i) = 0
      enddo

      done = 1
      pi = atan(done)*4

      nsrc = 9998
      ntarg = 9999

      allocate(sources(2,nsrc),charges(nsrc),dip1(nsrc),dip2(nsrc))
      allocate(targ(2,ntarg))
      allocate(pot(nsrc),grada(nsrc),gradaa(nsrc),hess(3,nsrc))
      allocate(pottarg(ntarg),gradatarg(ntarg),hesstarg(3,ntarg))
      allocate(gradaatarg(ntarg))

      do i=1,nsrc

         sources(1,i) = hkrand(0)
         sources(2,i) = hkrand(0)

         charges(i) = hkrand(0) + ima*hkrand(0) 
         dip1(i) = hkrand(0) + ima*hkrand(0)
         dip2(i) = hkrand(0) + ima*hkrand(0)
      enddo


      do i=1,ntarg
         targ(1,i) = hkrand(0)
         targ(2,i) = hkrand(0)
      enddo

      nts = min(20,nsrc)
      ntt = min(20,ntarg)

      allocate(potex(nts),gradaex(nts),gradaaex(nts),hessex(3,nts))
      allocate(pottargex(ntt),gradatargex(ntt),hesstargex(3,ntt))
      allocate(gradaatargex(ntt))

      eps = 0.5d-6

      ifcharge = 1
      ifdipole = 1
      iper = 0
      ifpgh = 2
      ifpghtarg = 2
      nd = 1
      call bhfmm2d(nd,eps,nsrc,sources,ifcharge,charges,
     1  ifdipole,dip1,dip2,iper,ifpgh,pot,grada,gradaa,hess,
     2  ntarg,targ,ifpghtarg,pottarg,gradatarg,gradaatarg,
     3  hesstarg)


c
cc       test against exact potential
c
      call dzero(potex,2*nts)
      call dzero(pottargex,2*ntt)
      call dzero(gradaex,2*nts)
      call dzero(gradatargex,2*ntt)
      call dzero(gradaaex,2*nts)
      call dzero(gradaatargex,2*ntt)

      thresh = 1.0d-14

      call bhfmm2dpart_direct_vec(1,1,nsrc,1,nts,sources,ifcharge,
     1       charges,ifdipole,dip1,dip2,sources,ifpgh,potex,
     2       gradaex,gradaaex,hessex,thresh)

      call bhfmm2dpart_direct_vec(1,1,nsrc,1,ntt,sources,ifcharge,
     1       charges,ifdipole,dip1,dip2,targ,ifpghtarg,pottargex,
     2       gradatargex,gradaatargex,hesstargex,thresh)

     
      errps = 0
      errpt = 0

      errgas = 0
      errgat = 0

      errgaas = 0
      errgaat = 0

      errhs = 0
      errht = 0
      if(ifpgh.ge.1) call derr(potex,pot,2*nts,errps)
      if(ifpghtarg.ge.1) call derr(pottargex,pottarg,2*ntt,errpt)


      if(ifpgh.ge.2) call derr(gradaex,grada,2*nts,errgas)
      if(ifpghtarg.ge.2) call derr(gradatargex,gradatarg,2*ntt,errgat)

      if(ifpgh.ge.2) call derr(gradaaex,gradaa,2*nts,errgaas)
      if(ifpghtarg.ge.2) call derr(gradaatargex,gradaatarg,2*ntt,
     1   errgaat)
      call errprintbh(errps,errgas,errgaas,errhs,errpt,errgat,errgaat,
     1  errht)
      if (errps .lt. eps .and. errgas .lt. eps .and. errgaas .lt. eps 
     1     .and. errhs .lt. eps .and.
     1     errpt .lt. eps .and. errgat .lt. eps .and. errgaat .lt. eps
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
      subroutine errprintbh(errps,errgas,errgaas,errhs,errpt,errgat,
     1  errgaat,errht)
      implicit real *8 (a-h,o-z)
 1100 format(4(2x,e11.5))


      write(6,*) 'error in sources'
      write(6,*) 'pot err, grada err, gradaa err, hess err' 
      write(6,1100) errps,errgas,errgaas,errhs
      write(6,*) 
      write(6,*)
      write(6,* ) 'error in targets'
      write(6,*) 'pot err, grada err, gradaa err, hess err' 
      write(6,1100) errpt,errgat,errgaat,errht
      write(6,*)
      write(6,*)
      write(6,*)'==================='

      return
      end
      
