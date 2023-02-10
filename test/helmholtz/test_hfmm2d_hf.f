      implicit real *8 (a-h,o-z)
      real *8, allocatable :: sources(:,:),targ(:,:)
      real *8, allocatable :: dipvec(:,:)
      complex *16, allocatable :: charges(:),dipstr(:)
      complex *16, allocatable :: pot(:),grad(:,:),hess(:,:)
      complex *16, allocatable :: pottarg(:),gradtarg(:,:),hesstarg(:,:)
      character(len=100) str1

      integer ipass(27)
      
      complex *16 ima,zk
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = atan(done)*4

      nsrc = 1000
      ntarg = 1001
cc      ntarg = nsrc+10

      allocate(sources(2,nsrc),charges(nsrc),dipstr(nsrc))
      allocate(dipvec(2,nsrc))
      allocate(targ(2,ntarg))
      allocate(pot(nsrc),grad(2,nsrc),hess(3,nsrc))
      allocate(pottarg(ntarg),gradtarg(2,ntarg),hesstarg(3,ntarg))

      do i=1,nsrc
        sources(1,i) = hkrand(0)
        sources(2,i) = hkrand(0)

        charges(i) = hkrand(0) + ima*hkrand(0)
        dipstr(i) = hkrand(0) + ima*hkrand(0)

        dipvec(1,i) = hkrand(0)
        dipvec(2,i) = hkrand(0)
      enddo

      do i=1,ntarg
         targ(1,i) = hkrand(0)
         targ(2,i) = hkrand(0)
      enddo

c
cc      test low frequency mode
c
      zk = 0.5d0 + 0.1d0*ima
      zk = 10
      zk = 500.0d0
c      zk = 10000.0d0
      eps = 0.51d-8

      write(*,*) "=========================================="
      write(*,*) "Testing suite for hfmm2d"
      write(*,'(a,e11.4)') "Requested precision = ",eps
      write(6,*)
      write(6,*)

      open(unit=33,file='print_testres.txt',access='append')

      ntests = 27
      do i=1,ntests
        ipass(i) = 0
      enddo

      ntest = 20
      nts = min(ntest,nsrc)
      ntt = min(ntest,ntarg)
      thresh = 2.0d0**(-51)

c
c
cc      now test source to source, charge
c       with potentials
c
      write(6,*) 'testing source to sources'
      write(6,*) 'interaction: charges'
      write(6,*) 'output: potentials'
      write(6,*)
      write(6,*)



      call hfmm2d_s_c_p(eps,zk,nsrc,sources,charges,
     1        pot,ier)
      

      ifcharge = 1
      ifdipole = 0
      ifpgh = 1
      ifpghtarg = 0
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(1) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)


c
c
cc      now test source to source, charge
c       with gradients
c
      write(6,*) 'testing source to sources'
      write(6,*) 'interaction: charges'
      write(6,*) 'output: potentials + gradients'
      write(6,*)
      write(6,*)


      call hfmm2d_s_c_g(eps,zk,nsrc,sources,charges,
     1        pot,grad,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 0
      ifpgh = 2
      ifpghtarg = 0
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(2) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)


c
c
cc      now test source to source, charge
c       with hessians
c
      write(6,*) 'testing source to sources'
      write(6,*) 'interaction: charges'
      write(6,*) 'output: potentials + gradients + hessians'
      write(6,*)
      write(6,*)

      call hfmm2d_s_c_h(eps,zk,nsrc,sources,charges,
     1        pot,grad,hess,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 0
      ifpgh = 3
      ifpghtarg = 0
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(3) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)

c
c
c
c
cc      now test source to source, dipole
c       with potentials
c
      write(6,*) 'testing source to sources'
      write(6,*) 'interaction: dipoles'
      write(6,*) 'output: potentials'
      write(6,*)
      write(6,*)


      call hfmm2d_s_d_p(eps,zk,nsrc,sources,dipstr,
     1        dipvec,pot,ier)


c
cc       test against exact potential
c
      ifcharge = 0
      ifdipole = 1
      ifpgh = 1
      ifpghtarg = 0
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(4) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)


c
c
cc      now test source to source, dipole
c       with gradients
c
      write(6,*) 'testing source to sources'
      write(6,*) 'interaction: dipoles'
      write(6,*) 'output: potentials + gradients'
      write(6,*)
      write(6,*)


      call hfmm2d_s_d_g(eps,zk,nsrc,sources,dipstr,
     1        dipvec,pot,grad,ier)
c
cc       test against exact potential
c

      ifcharge = 0
      ifdipole = 1
      ifpgh = 2
      ifpghtarg = 0
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(5) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)



c
c
cc      now test source to source, dipole
c       with hessians
c
      write(6,*) 'testing source to sources'
      write(6,*) 'interaction: dipoles'
      write(6,*) 'output: potentials + gradients + hessians'
      write(6,*)
      write(6,*)

      call hfmm2d_s_d_h(eps,zk,nsrc,sources,dipstr,
     1        dipvec,pot,grad,hess,ier)
c
cc       test against exact potential
c

      ifcharge = 0
      ifdipole = 1
      ifpgh = 3
      ifpghtarg = 0
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(6) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)



c
c
c
c
cc      now test source to source, charge + dipole
c       with potentials
c
      write(6,*) 'testing source to sources'
      write(6,*) 'interaction: charges + dipoles'
      write(6,*) 'output: potentials'
      write(6,*)
      write(6,*)


      call hfmm2d_s_cd_p(eps,zk,nsrc,sources,charges,dipstr,
     1        dipvec,pot,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 1
      ifpgh = 1
      ifpghtarg = 0
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(7) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)

c
c
cc      now test source to source, charge + dipole
c       with gradients
c
      write(6,*) 'testing source to sources'
      write(6,*) 'interaction: charges + dipoles'
      write(6,*) 'output: potentials + gradients'
      write(6,*)
      write(6,*)


      call hfmm2d_s_cd_g(eps,zk,nsrc,sources,charges,dipstr,
     1        dipvec,pot,grad,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 1
      ifpgh = 2
      ifpghtarg = 0
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(8) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)

c
c
cc      now test source to source, charge + dipole
c       with hessians
c
      write(6,*) 'testing source to sources'
      write(6,*) 'interaction: charges + dipoles'
      write(6,*) 'output: potentials + gradients + hessians'
      write(6,*)
      write(6,*)

      call hfmm2d_s_cd_h(eps,zk,nsrc,sources,charges,dipstr,
     1        dipvec,pot,grad,hess,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 1
      ifpgh = 3
      ifpghtarg = 0
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(9) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)

c
cc      now test source to target, charge
c       with potentials
c
      write(6,*) 'testing source to targets'
      write(6,*) 'interaction: charges'
      write(6,*) 'output: potentials'
      write(6,*)
      write(6,*)

      call hfmm2d_t_c_p(eps,zk,nsrc,sources,charges,
     1        ntarg,targ,pottarg,ier)
      

      ifcharge = 1
      ifdipole = 0
      ifpgh = 0
      ifpghtarg = 1
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(10) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)


c
c
cc      now test source to target, charge
c       with gradients
c
      write(6,*) 'testing source to targets'
      write(6,*) 'interaction: charges'
      write(6,*) 'output: potentials + gradients'
      write(6,*)
      write(6,*)


      call hfmm2d_t_c_g(eps,zk,nsrc,sources,charges,
     1        ntarg,targ,pottarg,gradtarg,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 0
      ifpgh = 0
      ifpghtarg = 2
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(11) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)


c
c
cc      now test source to target, charge
c       with hessians
c
      write(6,*) 'testing source to targets'
      write(6,*) 'interaction: charges'
      write(6,*) 'output: potentials + gradients + hessians'
      write(6,*)
      write(6,*)

      call hfmm2d_t_c_h(eps,zk,nsrc,sources,charges,
     1        ntarg,targ,pottarg,gradtarg,
     2        hesstarg,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 0
      ifpgh = 0
      ifpghtarg = 3
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(12) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)

c
c
c
c
cc      now test source to target, dipole
c       with potentials
c
      write(6,*) 'testing source to targets'
      write(6,*) 'interaction: dipoles'
      write(6,*) 'output: potentials'
      write(6,*)
      write(6,*)


      call hfmm2d_t_d_p(eps,zk,nsrc,sources,dipstr,
     1        dipvec,ntarg,targ,pottarg,ier)


c
cc       test against exact potential
c
      ifcharge = 0
      ifdipole = 1
      ifpgh = 0
      ifpghtarg = 1
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(13) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)


c
c
cc      now test source to target, dipole
c       with gradients
c
      write(6,*) 'testing source to targets'
      write(6,*) 'interaction: dipoles'
      write(6,*) 'output: potentials + gradients'
      write(6,*)
      write(6,*)


      call hfmm2d_t_d_g(eps,zk,nsrc,sources,dipstr,
     1        dipvec,ntarg,targ,pottarg,gradtarg,ier)
c
cc       test against exact potential
c

      ifcharge = 0
      ifdipole = 1
      ifpgh = 0
      ifpghtarg = 2
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(14) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)



c
c
cc      now test source to target, dipole
c       with hessians
c
      write(6,*) 'testing source to targets'
      write(6,*) 'interaction: dipoles'
      write(6,*) 'output: potentials + gradients + hessians'
      write(6,*)
      write(6,*)

      call hfmm2d_t_d_h(eps,zk,nsrc,sources,dipstr,
     1        dipvec,ntarg,targ,pottarg,gradtarg,
     2        hesstarg,ier)
c
cc       test against exact potential
c

      ifcharge = 0
      ifdipole = 1
      ifpgh = 0
      ifpghtarg = 3
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(15) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)




cccc      nsrc = n
cccc      ntarg = nsrc

c
c
c
c
cc      now test source to target, charge + dipole
c       with potentials
c
      write(6,*) 'testing source to targets'
      write(6,*) 'interaction: charges + dipoles'
      write(6,*) 'output: potentials'
      write(6,*)
      write(6,*)



      call hfmm2d_t_cd_p(eps,zk,nsrc,sources,charges,dipstr,
     1        dipvec,ntarg,targ,pottarg,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 1
      ifpgh = 0
      ifpghtarg = 1
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(16) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)

c
c
cc      now test source to target, charge + dipole
c       with gradients
c
      write(6,*) 'testing source to targets'
      write(6,*) 'interaction: charges + dipoles'
      write(6,*) 'output: potentials + gradients'
      write(6,*)
      write(6,*)


      call hfmm2d_t_cd_g(eps,zk,nsrc,sources,charges,dipstr,
     1        dipvec,ntarg,targ,pottarg,gradtarg,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 1
      ifpgh = 0 
      ifpghtarg = 2
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(17) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)

c
c
cc      now test source to target, charge + dipole
c       with hessians
c
      write(6,*) 'testing source to targets'
      write(6,*) 'interaction: charges + dipoles'
      write(6,*) 'output: potentials + gradients + hessians'
      write(6,*)
      write(6,*)

      call hfmm2d_t_cd_h(eps,zk,nsrc,sources,charges,dipstr,
     1        dipvec,ntarg,targ,pottarg,gradtarg,
     2        hesstarg,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 1
      ifpgh = 0
      ifpghtarg = 3
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(18) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)

c
c
cc      now test source to source  + target, charge
c       with potentials
c
      write(6,*) 'testing source to source and targets'
      write(6,*) 'interaction: charges'
      write(6,*) 'output: potentials'
      write(6,*)
      write(6,*)

      call hfmm2d_st_c_p(eps,zk,nsrc,sources,charges,
     1        pot,ntarg,targ,pottarg,ier)
      

      ifcharge = 1
      ifdipole = 0
      ifpgh = 1
      ifpghtarg = 1
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(19) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)


c
c
cc      now test source to source  + target, charge
c       with gradients
c
      write(6,*) 'testing source to source and targets'
      write(6,*) 'interaction: charges'
      write(6,*) 'output: potentials + gradients'
      write(6,*)
      write(6,*)


      call hfmm2d_st_c_g(eps,zk,nsrc,sources,charges,
     1        pot,grad,ntarg,targ,pottarg,gradtarg,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 0
      ifpgh = 2
      ifpghtarg = 2
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(20) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)


c
c
cc      now test source to source  + target, charge
c       with hessians
c
      write(6,*) 'testing source to source and targets'
      write(6,*) 'interaction: charges'
      write(6,*) 'output: potentials + gradients + hessians'
      write(6,*)
      write(6,*)

      call hfmm2d_st_c_h(eps,zk,nsrc,sources,charges,
     1        pot,grad,hess,ntarg,targ,pottarg,gradtarg,
     2        hesstarg,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 0
      ifpgh = 3
      ifpghtarg = 3
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(21) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)

c
c
c
c
cc      now test source to source  + target, dipole
c       with potentials
c
      write(6,*) 'testing source to source and targets'
      write(6,*) 'interaction: dipoles'
      write(6,*) 'output: potentials'
      write(6,*)
      write(6,*)


      call hfmm2d_st_d_p(eps,zk,nsrc,sources,dipstr,
     1        dipvec,pot,ntarg,targ,pottarg,ier)


c
cc       test against exact potential
c
      ifcharge = 0
      ifdipole = 1
      ifpgh = 1
      ifpghtarg = 1
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(22) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)


c
c
cc      now test source to source  + target, dipole
c       with gradients
c
      write(6,*) 'testing source to source and targets'
      write(6,*) 'interaction: dipoles'
      write(6,*) 'output: potentials + gradients'
      write(6,*)
      write(6,*)


      call hfmm2d_st_d_g(eps,zk,nsrc,sources,dipstr,
     1        dipvec,pot,grad,ntarg,targ,pottarg,gradtarg,ier)
c
cc       test against exact potential
c

      ifcharge = 0
      ifdipole = 1
      ifpgh = 2
      ifpghtarg = 2
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(23) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)



c
c
cc      now test source to source  + target, dipole
c       with hessians
c
      write(6,*) 'testing source to source and targets'
      write(6,*) 'interaction: dipoles'
      write(6,*) 'output: potentials + gradients + hessians'
      write(6,*)
      write(6,*)

      call hfmm2d_st_d_h(eps,zk,nsrc,sources,dipstr,
     1        dipvec,pot,grad,hess,ntarg,targ,pottarg,gradtarg,
     2        hesstarg,ier)
c
cc       test against exact potential
c

      ifcharge = 0
      ifdipole = 1
      ifpgh = 3
      ifpghtarg = 3
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(24) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)



c
c
c
c
cc      now test source to source  + target, charge + dipole
c       with potentials
c
      write(6,*) 'testing source to source and targets'
      write(6,*) 'interaction: charges + dipoles'
      write(6,*) 'output: potentials'
      write(6,*)
      write(6,*)


      call hfmm2d_st_cd_p(eps,zk,nsrc,sources,charges,dipstr,
     1        dipvec,pot,ntarg,targ,pottarg,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 1
      ifpgh = 1
      ifpghtarg = 1
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(25) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)

c
c
cc      now test source to source  + target, charge + dipole
c       with gradients
c
      write(6,*) 'testing source to source and targets'
      write(6,*) 'interaction: charges + dipoles'
      write(6,*) 'output: potentials + gradients'
      write(6,*)
      write(6,*)


      call hfmm2d_st_cd_g(eps,zk,nsrc,sources,charges,dipstr,
     1        dipvec,pot,grad,ntarg,targ,pottarg,gradtarg,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 1
      ifpgh = 2
      ifpghtarg = 2
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(26) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)

c
c
cc      now test source to source  + target, charge + dipole
c       with hessians
c
      write(6,*) 'testing source to source and targets'
      write(6,*) 'interaction: charges + dipoles'
      write(6,*) 'output: potentials + gradients + hessians'
      write(6,*)
      write(6,*)

      call hfmm2d_st_cd_h(eps,zk,nsrc,sources,charges,dipstr,
     1        dipvec,pot,grad,hess,ntarg,targ,pottarg,gradtarg,
     2        hesstarg,ier)
c
cc       test against exact potential
c

      ifcharge = 1
      ifdipole = 1
      ifpgh = 3
      ifpghtarg = 3
c      
c   compute error 
c

      call comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1  dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2  pottarg,gradtarg,hesstarg,ntest,erra)
      

      call prin2('l2 rel error=*',erra,1)
      write(6,*)
      write(6,*)
      write(6,*) '====================='

      if(erra.lt.eps) ipass(27) = 1
      call geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      if(erra.ge.eps) write(33,*) str1(1:len1)


      isum = 0
      do i=1,ntests
        isum = isum+ipass(i)
      enddo

      write(*,'(a,i2,a,i2,a)') 'Successfully completed ',isum,
     1   ' out of ',ntests,' tests in hf-hfmm2d testing suite'
      write(33,'(a,i2,a,i2,a)') 'Successfully completed ',isum,
     1   ' out of ',ntests,' tests in hf-hfmm2d testing suite'
      close(33)
      


      stop
      end
c      
c
c
c
c
      subroutine comperr(zk,nsrc,sources,ifcharge,charges,ifdipole,
     1   dipstr,dipvec,ifpgh,pot,grad,hess,ntarg,targ,ifpghtarg,
     2   pottarg,gradtarg,hesstarg,ntest,erra)
      implicit real *8 (a-h,o-z)
      complex *16 zk
      real *8 sources(2,*),targ(2,*),dipvec(2,*)
      complex *16 charges(*),dipstr(*)
      integer ifcharge,ifdipole,ifpgh,ifpghtarg
      complex *16 pot(*),grad(2,*),hess(3,*),pottarg(*)
      complex *16 gradtarg(2,*),hesstarg(3,*)
      real *8 erra

      complex *16, allocatable :: potex(:),gradex(:,:),hessex(:,:)
      complex *16, allocatable :: pottargex(:),gradtargex(:,:),
     1                             hesstargex(:,:)


      allocate(potex(ntest),gradex(2,ntest),hessex(3,ntest))
      allocate(pottargex(ntest),gradtargex(2,ntest),hesstargex(3,ntest))


      nts = 0
      ntt = 0
      if(ifpgh.ge.1) nts = min(ntest,nsrc)
      if(ifpghtarg.ge.1) ntt = min(ntest,ntarg)
      
      do i=1,nts
        if(ifpgh.ge.1) potex(i) = 0
        if(ifpgh.ge.2) then
          gradex(1,i) = 0
          gradex(2,i) = 0
        endif
        if(ifpgh.ge.3) then
          hessex(1,i) = 0
          hessex(2,i) = 0
          hessex(3,i) = 0
        endif
      enddo

      
      do i=1,ntt
        if(ifpghtarg.ge.1) pottargex(i) = 0
        if(ifpghtarg.ge.2) then
          gradtargex(1,i) = 0
          gradtargex(2,i) = 0
        endif
        if(ifpghtarg.ge.3) then
          hesstargex(1,i) = 0
          hesstargex(2,i) = 0
          hesstargex(3,i) = 0
        endif
      enddo




      thresh = 2.0d0**(-51)
      
      call hfmm2dpart_direct(1,1,nsrc,1,nts,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call hfmm2dpart_direct(1,1,nsrc,1,ntt,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,targ,ifpghtarg,pottargex,
     2       gradtargex,hesstargex,thresh)
      erra = 0
      ra = 0


      do i=1,nts
        if(ifpgh.ge.1) then
          erra = erra + abs(pot(i) - potex(i))**2
          ra = ra + abs(potex(i))**2
        endif
        if(ifpgh.ge.2) then
          erra = erra + abs(grad(1,i)-gradex(1,i))**2
          erra = erra + abs(grad(2,i)-gradex(2,i))**2
          ra = ra + abs(gradex(1,i))**2
          ra = ra + abs(gradex(2,i))**2
        endif
        if(ifpgh.ge.3) then
          erra = erra + abs(hess(1,i)-hessex(1,i))**2
          erra = erra + abs(hess(2,i)-hessex(2,i))**2
          erra = erra + abs(hess(3,i)-hessex(3,i))**2
          ra = ra + abs(hessex(1,i))**2
          ra = ra + abs(hessex(2,i))**2
          ra = ra + abs(hessex(3,i))**2
        endif
      enddo


      do i=1,ntt
        if(ifpghtarg.ge.1) then
          erra = erra + abs(pottarg(i) - pottargex(i))**2
          ra = ra + abs(pottargex(i))**2
        endif
        if(ifpghtarg.ge.2) then
          erra = erra + abs(gradtarg(1,i)-gradtargex(1,i))**2
          erra = erra + abs(gradtarg(2,i)-gradtargex(2,i))**2
          ra = ra + abs(gradtargex(1,i))**2
          ra = ra + abs(gradtargex(2,i))**2
        endif
        if(ifpghtarg.ge.3) then
          erra = erra + abs(hesstarg(1,i)-hesstargex(1,i))**2
          erra = erra + abs(hesstarg(2,i)-hesstargex(2,i))**2
          erra = erra + abs(hesstarg(3,i)-hesstargex(3,i))**2
          ra = ra + abs(hesstargex(1,i))**2
          ra = ra + abs(hesstargex(2,i))**2
          ra = ra + abs(hesstargex(3,i))**2
        endif
      enddo

      erra = sqrt(erra/ra)

      return
      end
