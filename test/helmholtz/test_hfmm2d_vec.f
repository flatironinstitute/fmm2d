      implicit real *8 (a-h,o-z)
      real *8, allocatable :: sources(:,:),targ(:,:)
      real *8, allocatable :: dipvec(:,:,:)
      complex *16, allocatable :: charges(:,:),dipstr(:,:)
      complex *16, allocatable :: pot(:,:),grad(:,:,:),hess(:,:,:)
      complex *16, allocatable :: pottarg(:,:),gradtarg(:,:,:),
     1    hesstarg(:,:,:)

      complex *16, allocatable :: potex(:,:),gradex(:,:,:),
     1   hessex(:,:,:)
      complex *16, allocatable :: pottargex(:,:),gradtargex(:,:,:),
     1                             hesstargex(:,:,:)

      real *8 expc(100),texps(100),scj(100)
      
      complex *16 ima,zk
      integer nd,idim
      data ima/(0.0d0,1.0d0)/

      call prini(6,13)

      done = 1
      pi = atan(done)*4

      call prin2('Enter n*',n,0)
      read *, n

      nd = 1


      nsrc = 100000
      ntarg = nsrc

      allocate(sources(2,nsrc),charges(nd,nsrc),dipstr(nd,nsrc))
      allocate(targ(2,ntarg),dipvec(nd,2,nsrc))
      allocate(pot(nd,nsrc),grad(nd,2,nsrc),hess(nd,3,nsrc))
      allocate(pottarg(nd,ntarg),gradtarg(nd,2,ntarg),
     1   hesstarg(nd,3,ntarg))

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

           charges(idim,i) = hkrand(0) + ima*hkrand(0)
           dipstr(idim,i) = hkrand(0) + ima*hkrand(0)
           dipvec(idim,1,i) = hkrand(0)
           dipvec(idim,2,i) = hkrand(0)
        enddo
      enddo



      do i=1,ntarg
         targ(1,i) = hkrand(0)
         targ(2,i) = hkrand(0)
      enddo

      nts = min(20,nsrc)
      ntt = min(20,ntarg)

      allocate(potex(nd,nts),gradex(nd,2,nts),hessex(nd,3,nts))
      allocate(pottargex(nd,ntt),gradtargex(nd,2,ntt),
     1   hesstargex(nd,3,ntt))


c
cc      test low frequency mode
c
      zk = 600.5d0 + 0.1d0*ima
      zk = 300.0d0*ima
      eps = 0.5d-14
c
c
cc      now test source to source  + target, charge
c       with potentials
c
      write(*,*) 'testing stost, charge, potentials'

      call hfmm2dpartstostcp_vec(nd,eps,zk,nsrc,sources,charges,
     1        pot,ntarg,targ,pottarg)


c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(pottargex,2*ntt*nd)

      ifcharge = 1
      ifdipole = 0

      ifpgh = 1
      thresh = 1.0d-16

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,nts,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,ntt,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      call derr(potex,pot,2*nts*nd,errps)
      call derr(pottargex,pottarg,2*ntt*nd,errpt)

      errgs = 0
      errgt = 0

      errhs = 0
      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht)


c
c
cc      now test source to source  + target, charge
c       with gradients
c
      write(*,*) 'testing stost, charge, gradients'

      call hfmm2dpartstostcg_vec(nd,eps,zk,nsrc,sources,charges,
     1        pot,grad,ntarg,targ,pottarg,gradtarg)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(gradex,4*nts*nd)
      call dzero(pottargex,2*ntt*nd)
      call dzero(gradtargex,4*ntt*nd)

      ifcharge = 1
      ifdipole = 0

      ifpgh = 2
      thresh = 1.0d-16

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,nts,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,ntt,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      call derr(potex,pot,2*nts*nd,errps)
      call derr(gradex,grad,4*nts*nd,errgs)
      call derr(pottargex,pottarg,2*ntt*nd,errpt)
      call derr(gradtargex,gradtarg,4*ntt*nd,errgt)

      errhs = 0
      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht)


c
c
cc      now test source to source  + target, charge
c       with hessians
c



      write(*,*) 'testing stost, charge, hessians'
      call hfmm2dpartstostch_vec(nd,eps,zk,nsrc,sources,charges,
     1        pot,grad,hess,ntarg,targ,pottarg,gradtarg,
     2        hesstarg)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(gradex,4*nts*nd)
      call dzero(hessex,6*nts*nd)
      call dzero(pottargex,2*ntt*nd)
      call dzero(gradtargex,4*ntt*nd)
      call dzero(hesstargex,6*ntt*nd)

      ifcharge = 1
      ifdipole = 0

      ifpgh = 3
      thresh = 1.0d-16

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,nts,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,ntt,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      call derr(potex,pot,2*nts*nd,errps)
      call derr(gradex,grad,4*nts*nd,errgs)
      call derr(hessex,hess,6*nts*nd,errhs)
      call derr(pottargex,pottarg,2*ntt*nd,errpt)
      call derr(gradtargex,gradtarg,4*ntt*nd,errgt)
      call derr(hesstargex,hesstarg,6*nts*nd,errht)

      call errprint(errps,errgs,errhs,errpt,errgt,errht)


c
c
c
c
cc      now test source to source  + target, dipole
c       with potentials
c
      write(*,*) 'testing stost, dipole, potentials'

      call hfmm2dpartstostdp_vec(nd,eps,zk,nsrc,sources,dipstr,
     1        dipvec,pot,ntarg,targ,pottarg)


c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(pottargex,2*ntt*nd)

      ifcharge = 0
      ifdipole = 1

      ifpgh = 1
      thresh = 1.0d-16

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,nts,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,ntt,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      call derr(potex,pot,2*nts*nd,errps)
      call derr(pottargex,pottarg,2*ntt*nd,errpt)

  
      errgs = 0
      errgt = 0

      errhs = 0
      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht)


c
c
cc      now test source to source  + target, dipole
c       with gradients
c
      write(*,*) 'testing stost, dipole, gradients'

      call hfmm2dpartstostdg_vec(nd,eps,zk,nsrc,sources,dipstr,
     1        dipvec,pot,grad,ntarg,targ,pottarg,gradtarg)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(gradex,4*nts*nd)
      call dzero(pottargex,2*ntt*nd)
      call dzero(gradtargex,4*ntt*nd)

      ifcharge = 0
      ifdipole = 1

      ifpgh = 2
      thresh = 1.0d-16

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,nts,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,ntt,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      call derr(potex,pot,2*nts*nd,errps)
      call derr(gradex,grad,4*nts*nd,errgs)
      call derr(pottargex,pottarg,2*ntt*nd,errpt)
      call derr(gradtargex,gradtarg,4*ntt*nd,errgt)

      errhs = 0
      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht)


c
c
cc      now test source to source  + target, dipole
c       with hessians
c



      write(*,*) 'testing stost, dipole, hessians'
      call hfmm2dpartstostdh_vec(nd,eps,zk,nsrc,sources,dipstr,
     1        dipvec,pot,grad,hess,ntarg,targ,pottarg,gradtarg,
     2        hesstarg)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(gradex,4*nts*nd)
      call dzero(hessex,6*nts*nd)
      call dzero(pottargex,2*ntt*nd)
      call dzero(gradtargex,4*ntt*nd)
      call dzero(hesstargex,6*ntt*nd)

      ifcharge = 0
      ifdipole = 1

      ifpgh = 3
      thresh = 1.0d-16

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,nts,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,ntt,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      call derr(potex,pot,2*nts*nd,errps)
      call derr(gradex,grad,4*nts*nd,errgs)
      call derr(hessex,hess,6*nts*nd,errhs)
      call derr(pottargex,pottarg,2*ntt*nd,errpt)
      call derr(gradtargex,gradtarg,4*ntt*nd,errgt)
      call derr(hesstargex,hesstarg,6*nts*nd,errht)

      call errprint(errps,errgs,errhs,errpt,errgt,errht)

c
c
c
c
cc      now test source to source  + target, charge + dipole
c       with potentials
c
      write(*,*) 'testing stost, charge + dipole, potentials'

      call hfmm2dpartstostcdp_vec(nd,eps,zk,nsrc,sources,charges,
     1     dipstr,dipvec,pot,ntarg,targ,pottarg)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(pottargex,2*ntt*nd)

      ifcharge = 1
      ifdipole = 1

      ifpgh = 1
      thresh = 1.0d-16

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,nts,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,ntt,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      call derr(potex,pot,2*nts*nd,errps)
      call derr(pottargex,pottarg,2*ntt*nd,errpt)

  
      errgs = 0
      errgt = 0

      errhs = 0
      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht)


c
c
cc      now test source to source  + target, charge + dipole
c       with gradients
c
      write(*,*) 'testing stost, charge + dipole, gradients'

      call hfmm2dpartstostcdg_vec(nd,eps,zk,nsrc,sources,charges,
     1    dipstr,dipvec,pot,grad,ntarg,targ,pottarg,gradtarg)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(gradex,4*nts*nd)
      call dzero(pottargex,2*ntt*nd)
      call dzero(gradtargex,4*ntt*nd)

      ifcharge = 1
      ifdipole = 1

      ifpgh = 2
      thresh = 1.0d-16

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,nts,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,ntt,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      call derr(potex,pot,2*nts*nd,errps)
      call derr(gradex,grad,4*nts*nd,errgs)
      call derr(pottargex,pottarg,2*ntt*nd,errpt)
      call derr(gradtargex,gradtarg,4*ntt*nd,errgt)

      errhs = 0
      errht = 0
      call errprint(errps,errgs,errhs,errpt,errgt,errht)


c
c
cc      now test source to source  + target, charge + dipole
c       with hessians
c



      write(*,*) 'testing stost, charge + dipole, hessians'
      call hfmm2dpartstostcdh_vec(nd,eps,zk,nsrc,sources,charges,
     1    dipstr,dipvec,pot,grad,hess,ntarg,targ,pottarg,gradtarg,
     2        hesstarg)
c
cc       test against exact potential
c
      call dzero(potex,2*nts*nd)
      call dzero(gradex,4*nts*nd)
      call dzero(hessex,6*nts*nd)
      call dzero(pottargex,2*ntt*nd)
      call dzero(gradtargex,4*ntt*nd)
      call dzero(hesstargex,6*ntt*nd)

      ifcharge = 1
      ifdipole = 1

      ifpgh = 3
      thresh = 1.0d-16

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,nts,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,sources,ifpgh,potex,
     2       gradex,hessex,thresh)

      call hfmm2dpart_direct_vec(nd,1,nsrc,1,ntt,zk,sources,ifcharge,
     1       charges,ifdipole,dipstr,dipvec,targ,ifpgh,pottargex,
     2       gradtargex,hesstargex,thresh)

      call derr(potex,pot,2*nts*nd,errps)
      call derr(gradex,grad,4*nts*nd,errgs)
      call derr(hessex,hess,6*nts*nd,errhs)
      call derr(pottargex,pottarg,2*ntt*nd,errpt)
      call derr(gradtargex,gradtarg,4*ntt*nd,errgt)
      call derr(hesstargex,hesstarg,6*nts*nd,errht)

      call errprint(errps,errgs,errhs,errpt,errgt,errht)

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
      
