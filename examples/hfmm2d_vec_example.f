      implicit none
      integer ns,nt,nd
      double precision, allocatable :: source(:,:),targ(:,:)
      double complex, allocatable :: charge(:,:),dipstr(:,:)
      real *8, allocatable :: dipvec(:,:,:)
      double complex, allocatable :: pot(:,:)
      double complex, allocatable :: pottarg(:,:)

      double precision eps
      double complex eye,zk
      integer i,j,k,idim,ier
      double precision hkrand,pi,thet,phi,t1,t2,omp_get_wtime
      

      data eye/(0.0d0,1.0d0)/

      pi = 4.0d0*atan(1.0d0)

c
cc      initialize printing routine
c
      call prini(6,13)
      write(*,*)
      write(*,*)
      write(*,*) "================================="
      call prin2("This code is an example fortran driver*",i,0)
      call prin2("On output, the code prints sample pot,pottarg*",i,0)
      write(*,*)
      write(*,*)

      zk = 2.2d0

      ns = 1600
      nt = 62500 
      
      nd = 1

      allocate(source(2,ns))
      allocate(targ(2,nt))
      allocate(charge(nd,ns))
      allocate(dipstr(nd,ns))
      allocate(dipvec(nd,2,ns))
      allocate(pot(nd,ns))
      allocate(pottarg(nd,nt))


      eps = 0.51d-6

      write(*,*) "=========================================="

c
c   
c       example demonstrating use of 
c        source to source+targ, charges+dipoles, pot
c



c
cc      generate sources uniformly on the sphere 
c
c
      do i=1,ns
        thet = hkrand(0)*pi
        phi = hkrand(0)*2*pi
        source(1,i) = sin(thet)*cos(phi) 
        source(2,i) = sin(thet)*sin(phi)

        do idim=1,nd
          charge(idim,i) = hkrand(0) + eye*hkrand(0)
          dipstr(idim,i) = hkrand(0) + eye*hkrand(0)
          dipvec(idim,1,i) = hkrand(0) 
          dipvec(idim,2,i) = hkrand(0) 
          pot(idim,i) = 0
        enddo
      enddo

      do i=1,nt
        targ(1,i) = -3.0d0 + hkrand(0)*6.0d0
        targ(2,i) = -3.0d0 + hkrand(0)*6.0d0 

        do idim=1,nd
          pottarg(idim,i) = 0
        enddo
      enddo

       call cpu_time(t1)
C$       t1 = omp_get_wtime()       

       call hfmm2d_t_cd_p(eps,zk,ns,source,charge,
     1      dipstr,dipvec,nt,targ,pottarg,ier)

       call cpu_time(t2)
C$       t2 = omp_get_wtime()
       
       call prin2('time taken=*',t2-t1,1)
       call prin2("potential at targets=*",pottarg,12)


      stop
      end
c----------------------------------------------------------
c
cc
c
c
