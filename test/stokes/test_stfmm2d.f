      implicit real *8 (a-h,o-z)
      real *8, allocatable :: sources(:,:),targ(:,:)
      real *8, allocatable :: stoklet(:,:),strslet(:,:),strsvec(:,:)
      real *8, allocatable :: pot(:,:),grad(:,:,:),pre(:)
      real *8, allocatable :: potex(:,:),gradex(:,:,:),preex(:)
      
      real *8, allocatable :: pottarg(:,:),gradtarg(:,:,:),pretarg(:)
      real *8, allocatable :: pottargex(:,:),gradtargex(:,:,:),
     1   pretargex(:)

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

      allocate(sources(2,nsrc),stoklet(2,nsrc))
      allocate(strslet(2,nsrc),strsvec(2,nsrc))
      allocate(targ(2,ntarg))
      allocate(pot(2,nsrc),grad(2,2,nsrc),pre(nsrc))
      allocate(pottarg(2,ntarg),gradtarg(2,2,ntarg),pretarg(ntarg))

      do i=1,nsrc

         sources(1,i) = hkrand(0)
         sources(2,i) = hkrand(0)

         stoklet(1,i) = hkrand(0)
         stoklet(2,i) = hkrand(0)
         
         strslet(1,i) = hkrand(0)
         strslet(2,i) = hkrand(0)

         strsvec(1,i) = hkrand(0)
         strsvec(2,i) = hkrand(0)
      enddo



      do i=1,ntarg
         targ(1,i) = hkrand(0)
         targ(2,i) = hkrand(0)
      enddo

      nts = min(20,nsrc)
      ntt = min(20,ntarg)

      allocate(potex(2,nts),gradex(2,2,nts),preex(nts))
      allocate(pottargex(2,ntt),gradtargex(2,2,ntt),pretargex(ntt))

      eps = 0.5d-6

      ifstoklet = 1
      ifstrslet = 1
      iper = 0
      ifppreg = 3
      ifppregtarg = 3
      nd = 1
      call stfmm2d(nd,eps,nsrc,sources,ifstoklet,stoklet,
     1  ifstrslet,strslet,strsvec,ifppreg,pot,pre,grad,
     2  ntarg,targ,ifppregtarg,pottarg,pretarg,gradtarg,ier)
      
c      call prin2('pot=*',pot,2*nsrc)
c      call prin2('pre=*',pre,nsrc)
c      call prin2('grad=*',grad,4*nsrc)


c
cc       test against exact potential
c
      call dzero(potex,2*nts)
      call dzero(pottargex,2*ntt)
      call dzero(gradex,4*nts)
      call dzero(gradtargex,4*ntt)
      call dzero(preex,nts)
      call dzero(pretargex,ntt)

      thresh = 1.0d-14

      call st2ddirectstokstrsg(nd,sources,ifstoklet,stoklet,ifstrslet,
     1     strslet,strsvec,nsrc,sources,nts,potex,preex,gradex,thresh)

      call st2ddirectstokstrsg(nd,sources,ifstoklet,stoklet,ifstrslet,
     1    strslet,strsvec,nsrc,targ,ntt,pottargex,pretargex,gradtargex,
     2    thresh)
      
c      call prin2('potex=*',potex,2*nts)
c      call prin2('preex=*',preex,nts)
c      call prin2('gradex=*',gradex,4*nts)


     
      errps = 0
      errpt = 0

      errgs = 0
      errgt = 0

      errpres = 0
      errpret = 0
      if(ifppreg.ge.1) call derr(potex,pot,2*nts,errps)
      if(ifppregtarg.ge.1) call derr(pottargex,pottarg,2*ntt,errpt)


      if(ifppreg.ge.2) call derr(preex,pre,nts,errpres)
      if(ifppregtarg.ge.2) call derr(pretargex,pretarg,ntt,errpret)

      if(ifppreg.ge.3) call derr(gradex,grad,4*nts,errgs)
      if(ifppregtarg.ge.3) call derr(gradtargex,gradtarg,4*ntt,errgt)

      call errprintbh(errps,errpres,errgs,errpt,errpret,
     1  errgt)
      if (errps .lt. eps .and. errgs .lt. eps 
     1     .and. errpres .lt. eps .and.
     1     errpt .lt. eps .and. errgt .lt. eps 
     1     .and. errpret .lt. eps)
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
      
