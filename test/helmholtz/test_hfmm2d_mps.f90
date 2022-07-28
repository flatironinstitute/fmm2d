program test_hfmm3d_mp2loc
  implicit double precision (a-h,o-z)
  
  character(len=72) str1
  
  integer :: ns, nt, nc
  integer :: i,j,k,ntest,nd,idim
  integer :: ifcharge,ifdipole,ifpgh,ifpghtarg
  integer :: ipass(18),len1,ntests,isum
  integer, allocatable :: nterms(:), impole(:)
  
  double precision :: eps, err, hkrand, dnorms(1000), force(10)
  double precision, allocatable :: source(:,:), targ(:,:)
  double precision, allocatable :: centers(:,:)
  double precision, allocatable :: wlege(:), rscales(:)
  
  double complex :: eye, zk, ima
  double complex, allocatable :: charge(:,:),dipstr(:,:)
  double precision, allocatable :: dipvec(:,:,:)
  double complex, allocatable :: pot(:,:), pot2(:,:), pottarg(:,:)
  double complex, allocatable :: grad(:,:,:),gradtarg(:,:,:)
  double complex, allocatable :: hess(:,:,:),hesstarg(:,:,:)
  double complex, allocatable :: mpole(:), local(:)


  data eye/(0.0d0,1.0d0)/
  ima = (0,1)
  done = 1
  pi = 4*atan(done)

  !
  ! initialize printing routine
  !
  call prini(6,13)

  nd = 1


  n1 = 10
  ns = n1**2
  nc = ns

  
  nt = 22

  allocate(source(2,ns),targ(2,nt), centers(2,nc))
  allocate(charge(nd,ns),dipvec(nd,2,ns),dipstr(nd,ns))
  allocate(pot(nd,ns), pot2(nd,ns))
  allocate(grad(nd,2,ns))
  allocate(hess(nd,3,ns))

  allocate(pottarg(nd,nt))
  allocate(gradtarg(nd,2,nt))
  allocate(hesstarg(nd,3,nt))
  eps = 0.5d-9

  write(*,*) "=========================================="
  write(*,*) "Testing suite for hfmm2d_mps"
  write(*,'(a,e12.5)') "Requested precision = ",eps

  boxsize = 1
  dlam = 1/boxsize
  zk = 2*pi/dlam 

  call prin2('boxsize in wavelengths = *', boxsize, 1)
  call prin2('dlam = *', dlam, 1)
  call prin2('zk = *', zk, 2)

  !
  ! generate sources uniformly in the unit cube 
  !
  h = 1.0d0/(n1+1)
  ijk = 0
  do i = 1,n1
    do j = 1,n1
      ijk = ijk + 1
      source(1,ijk) = h*i
      source(2,ijk) = h*j
    end do
  end do
  
  
  dnorm = 0
  dnormd = 0
  do i=1,ns
    !source(1,i) = hkrand(0)**2
    !source(2,i) = hkrand(0)**2
    !source(3,i) = hkrand(0)**2

    do idim=1,nd

      charge(idim,i) = hkrand(0) + eye*hkrand(0)
      dipstr(idim,i) = hkrand(0) + eye*hkrand(0)
      dnorm = dnorm + abs(charge(idim,i))**2
      dnormd = dnormd + abs(dipstr(idim,i))**2
      
      dipvec(idim,1,i) = hkrand(0) 
      dipvec(idim,2,i) = hkrand(0) 

      pot(idim,i) = 0
      grad(idim,1,i) = 0
      grad(idim,2,i) = 0

      hess(idim,1,i) = 0
      hess(idim,2,i) = 0
      hess(idim,3,i) = 0
    enddo
  enddo

  dnorm = sqrt(dnorm)
  dnormd = sqrt(dnormd)
  do i=1,ns
    do idim = 1,nd
      charge(idim,i) = charge(idim,i)/dnorm
      charge(idim,i) = 0
      charge(idim,1) = 1.0d0
      dipstr(idim,i) = dipstr(idim,i)/dnormd
      dipstr(idim,i) = 0
    end do
  end do
  

  ! call prin2('min source separation = *', ssep, 1)
  
  shift = h/100
  do i = 1,ns
    centers(1,i) = source(1,i) + shift
    centers(2,i) = source(2,i)
  end do

  !call prin2('centers = *', centers, 3*nc)

  !
  ! now form a multipole expansion at each center
  !
  allocate(nterms(nc), impole(nc))

  ntm = 30
  ntot1 = 0
  do i = 1,nc
    nterms(i) = ntm + 2*cos(pi/2*i)
    nterms(i) = ntm
    ntot1 = ntot1 + (2*nterms(i)+1)

  end do

  ntot = nd*ntot1


  allocate(mpole(ntot))

  impole(1) = 1
  do i = 1,nc-1
    ilen = (2*nterms(i)+1)
    impole(i+1) = impole(i) + nd*ilen
  end do
  
  call zinitialize(ntot, mpole)
  
  ns1 = 1
  rscale = 1
  sc = abs(zk)*shift
  if (sc .lt. 1) rscale = sc

  
  allocate(rscales(nc))
  do i = 1,nc
    rscales(i) = rscale
    call h2dformmpc(nd,zk,rscales(i),source(1,i),ns1,charge(1,i), &
       centers(1,i),nterms(i),mpole(impole(i)))
  end do

  
  !
  ! do the direct calculation
  !
  thresh = 1.0d-15
  ifcharge = 1
  ifdipole = 0
  ifpgh = 1
  ntarg = 0
  ifpghtarg = 0
  call hfmm2d(nd, eps, zk, ns, source, ifcharge, &
      charge, ifdipole,dipstr,dipvec, iper, ifpgh, pot, grad, &
      hess, ntarg, targ, ifpghtarg, pottarg, gradtarg, hesstarg, ier)
  

  
  allocate(local(ntot))
  
  !
  ! now test source to source, charge, 
  ! with potentials
  !
  print *
  print *
  print *
  write(6,*) 'testing multipoles to locals'
  write(6,*) 'input: multipole expansions'
  write(6,*) 'output: local expansions'
  write(6,*) 
  write(6,*) 


  call hfmm2d_mps(nd, eps, zk, &
      nc, centers, rscales, nterms, ntot,mpole, impole, local,ier)
  

  call zinitialize(nd*nc, pot2)
  npts = 1
  do i = 1,nc
    call h2dtaevalp(nd, zk, rscales(i), &
        centers(1,i), local(impole(i)), &
        nterms(i), source(1,i), npts, pot2(1,i))
  end do
  
  nprin = min(2*nd*nc,24)
  call prin2('from hfmm3d_mps, potential = *', pot2, nprin)
  call prin2('via fmm, potential = *', pot, nprin)

  erra = 0
  dnorm = 0
  do j = 1,nc
    do i = 1,nd
      erra = erra + abs(pot(i,j)-pot2(i,j))**2
      dnorm = dnorm + abs(pot(i,j))**2
      if(abs(pot(i,j)-pot2(i,j)).ge.1.0d-6) print *, real(pot(1,j)),real(pot2(1,j))

    end do
  end do

  
  erra = sqrt(erra/dnorm)
  call prin2('l2 rel err=*',erra,1)

  open(unit=33,file='print_testres.txt',access='append')
  isuccess = 0
  ntest = 1
  if(erra.lt.eps) isuccess = 1

  write(33,'(a,i1,a,i1,a)') 'Successfully completed ', &
    isuccess,' out of ',ntest,' tests in helm2d_mps testing suite'
  close(33)

  
  

  stop
end program



! ----------------------------------------------------------
! 
! This is the end of the debugging code.
!
! ----------------------------------------------------------
subroutine prinmp(str, mpole, nterms)
  implicit double precision (a-h,o-z)
  character (len=*) :: str
  double complex :: mpole(0:nterms, -nterms:nterms)
  double complex :: tmp(-10000:10000)

  print *, trim(str)

  do n = 0,nterms
    print *
    write(6,'(a,i2,a,i2,a,i2)') 'n = ', n, '  m = ', -n, '...', n
    do m = -nterms,nterms
      tmp(m) = mpole(n,m)
    end do
    call prin2('coefs = *', tmp(-n), 2*n+1)
  end do
  
  return
end subroutine prinmp



subroutine zinitialize(len, zs)
  implicit double precision (a-h,o-z)
  double complex :: zs(len)

  do i = 1,len
    zs(i) = 0
  end do
  return
end subroutine zinitialize




