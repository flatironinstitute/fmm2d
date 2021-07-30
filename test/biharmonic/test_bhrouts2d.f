c
c      TESTING SUITE FOR BIHARMONIC SUBROUTINE LIBRARY:  Jul 29, 2021
c
c
c
c         Sources at S, Targets at T.
c
c         bh2dformmpcd_vec      forms expansion about x0
c         bh2dmpevalg_vec       evaluates expansion at T
c         bh2dmpmp_vec          shifts expansion to x1
c         bh2dmploc_vec         converts to local about y1
c         bh2dtaevalg_vec       evaluates expansion at T
c         bh2dlocloc_vec        shifts to local about y0
c         bh2dtaevalg_vec       evaluates expansion at T
c         bh2dformtacd_vec      forms local expansion about y0 
c                               directly from S
c         bh2dtaevalg_vec       evaluates expansion at T
c
c         -------------------------------------------------
c         |       |   S   |       |       |T      |       |
c         |       |   c0  |       |       |  c3   |       |
c         |       |       |       |       |       |       |
c         -------c1-------------------------------c2------
c         |       |       |       |       |       |       |
c         |       |       |       |       |       |       |
c         |       |       |       |       |       |       |
c         -------------------------------------------------
c
      implicit real *8 (a-h,o-z)
c
      integer lw,nd
      parameter (lw=100000)
      parameter (nd=2)
      integer iw(0:lw)
      real *8 ztrg(2),source(2,10)
      real *8 ztrg2(2)
      real *8 c0(2),c1(2),c2(2),c3(2)
      real *8, allocatable :: carray(:,:)
c
      complex *16 w(lw)
      complex *16 w2(-lw:lw)
c
      complex *16 pot(nd),grada(nd),gradaa(nd)
      complex *16 opot(nd),ograda(nd),ogradaa(nd)
      complex *16, allocatable :: mpolecd(:,:,:)
      complex *16, allocatable :: mpole2cd(:,:,:)
      complex *16, allocatable :: localcd(:,:,:)
      complex *16, allocatable :: local2cd(:,:,:)
      complex *16, allocatable :: charge(:,:),dip1(:,:),dip2(:,:)
      complex *16, allocatable :: dummy(:,:)
      complex *16 eye
c
      data eye/(0.0d0,1.0d0)/
c
      done=1
      pi=4.0*atan(done)
      call prini(6,13)
c
c  get carray
c
c
      ldc = 120
      allocate(carray(0:ldc,0:ldc))
      call init_carray(carray,ldc)
c
c       set scale parameter neq 1 to make sure it is used correctly.
c
c        rscale = 1.0d0
        rscale = 1.1d0 
        size = 1.0d0
c
c       set center for original multipole expansion
c
      c0(1)=0.248d0
      c0(2)=0.251d0
c
c       create two charge sources.
c
      ns = 2
      allocate(charge(nd,ns),dip1(nd,ns),dip2(nd,ns),dummy(nd,ns))
      call mksource(nd,ns,source,charge,dip1,dip2,dummy,c0)
      call prin2(' charge is *',charge,2*nd*ns)
      call prin2(' dip1 is *',dip1,2*nd*ns)
      call prin2(' dip2 is *',dip2,2*ns*nd)
c
c
c       create center for shifted expansions
c       (location jiggled for good measure)
c
      c1(1)= 0.015d0
      c1(2)= 0.012d0
      c2(1)= 2.015d0
      c2(2)= 0.012d0
      c3(1)= c2(1)-0.25d0
      c3(2)= c2(2)-0.249d0
c
c       create target
c
      ztrg(1)=c3(1)-0.245d0
      ztrg(2)=c3(2)-0.25d0
      call prin2(' xdiff from center is *',ztrg(1)-c3(1),1)
      call prin2(' ydiff from center is *',ztrg(2)-c3(2),1)
c
c       direct calculation:
c
      np = 1
      ng = 1
      call zero_out(nd,np,opot)
      call zero_out(nd,ng,ograda)
      call zero_out(nd,ng,ogradaa)

      thresh = 1.0d-16
      call bh2d_directcdg_vec(nd,source,ns,charge,dip1,dip2,
     1       ztrg,opot,ograda,ogradaa,thresh)
        call prin2('Via direct calculation, potential is*',opot,2*nd)
c
      call prin2('Via direct calculation, potential is*',opot,2*nd)
      call prin2('Via direct calculation, grada is*',ograda,2*nd)
      call prin2('Via direct calculation, grada is*',ogradaa,2*nd)
c
c
c       create h-expansion:
c
      rscale = 1.2d0
      rscale2 = 1.2d0
      rscale3 = 1.3d0

ccc        do 8000 iii = 0,4
      do 8000 iii = 4,4
          
        if (iii.eq.0) then 
          eps = 1.0d-3
        else if (iii.eq.1) then 
          eps = 1.0d-6
        else if (iii.eq.2) then 
          eps = 1.0d-9
        else if (iii.eq.3) then 
          eps = 1.0d-12
        else
          eps = 1.0d-14
        endif

        ier = 0
        call bh2dterms(eps,nterms,ier)
        allocate(mpolecd(nd,5,0:nterms))
        allocate(mpole2cd(nd,5,0:nterms))
        allocate(localcd(nd,5,0:nterms))
        allocate(local2cd(nd,5,0:nterms))

        call prinf('calling formmp and mpeval =*', lused,0)
        call bh2dmpzero_vec(nd,mpolecd,nterms)
        call bh2dformmpcd_vec(nd,rscale,source,ns,charge,
     1          dip1,dip2,c0,nterms,mpolecd)
c
c ... evaluate the h-expansion at the target point:
c
        ntrg = 1
c
        call zero_out(nd,np,pot)
        call zero_out(nd,ng,grada)
        call zero_out(nd,ng,gradaa)
        call bh2dmpevalg_vec(nd,rscale,c0,mpolecd,nterms,ztrg,ntrg,pot,
     1      grada,gradaa)
        call bherrprintvec(nd,pot,opot,grada,ograda,gradaa,ogradaa)
        call prin2('Via mpole calculation, potential is*',pot,2*nd)
        call prin2('Via mpole calculation, grada is*',grada,2*nd*ng)
        call prin2('Via mpole calculation, gradaa is*',gradaa,2*nd*ng)

c
ccc	stop
ccc	goto 8000
c    mpmp shift
c
        call prinf('calling mpmp and mpeval *',nterms,0)
        call bh2dmpzero_vec(nd,mpole2cd,nterms)
        call bh2dmpmp_vec(nd,rscale,c0,mpolecd,nterms,
     1       rscale2,c1,mpole2cd,nterms,carray,ldc)
        call zero_out(nd,np,pot)
        call zero_out(nd,ng,grada)
        call zero_out(nd,ng,gradaa)
        call bh2dmpevalg_vec(nd,rscale2,c1,mpole2cd,nterms,ztrg,ntrg,
     1      pot,grada,gradaa)
        call bherrprintvec(nd,pot,opot,grada,ograda,gradaa,ogradaa)
        call prin2('Via mpole calculation, potential is*',pot,2*nd)
        call prin2('Via mpole calculation, grada is*',grada,2*nd*ng)
        call prin2('Via mpole calculation, gradaa is*',gradaa,2*nd*ng)


c
ccc	   stop
c
c    convert to local
c
        call prin2('calling mploc and taeval *',wavek,0)
        call bh2dmpzero_vec(nd,localcd,nterms)
        call bh2dmploc_vec(nd,rscale2,c1,mpole2cd,nterms,
     1         rscale,c2,localcd,nterms,carray,ldc)
        call zero_out(nd,np,pot)
        call zero_out(nd,ng,grada)
        call zero_out(nd,ng,gradaa)
        call bh2dtaevalg_vec(nd,rscale,c2,localcd,nterms,ztrg,ntrg,
     1      pot,grada,gradaa)
        call bherrprintvec(nd,pot,opot,grada,ograda,gradaa,ogradaa)
        call prin2('Via mpole calculation, potential is*',pot,2*nd)
        call prin2('Via mpole calculation, grada is*',grada,2*nd*ng)
        call prin2('Via mpole calculation, gradaa is*',gradaa,2*nd*ng)

c
ccc	   stop
c
c    shift local and change scaling and nterms.
c
      call prinf('calling locloc and taeval *',nterms3,0)
      call bh2dmpzero_vec(nd,local2cd,nterms)
      call bh2dlocloc_vec(nd,rscale,c2,localcd,nterms,
     1       rscale3,c3,local2cd,nterms,carray,ldc)
        call zero_out(nd,np,pot)
        call zero_out(nd,ng,grada)
        call zero_out(nd,ng,gradaa)
        call bh2dtaevalg_vec(nd,rscale3,c3,local2cd,nterms,ztrg,ntrg,
     1      pot,grada,gradaa)
        call bherrprintvec(nd,pot,opot,grada,ograda,gradaa,ogradaa)
        call prin2('Via mpole calculation, potential is*',pot,2*nd)
        call prin2('Via mpole calculation, grada is*',grada,2*nd*ng)
        call prin2('Via mpole calculation, gradaa is*',gradaa,2*nd*ng)
c
c
c    create local exp from sources
c
      call prin2('calling bh2dformta and taeval *',nterms,0)
      call bh2dmpzero_vec(nd,local2cd,nterms)
      call bh2dformtacd_vec(nd,rscale2,source,ns,charge,dip1,dip2,
     1          c3,nterms,local2cd)
        call zero_out(nd,np,pot)
        call zero_out(nd,ng,grada)
        call zero_out(nd,ng,gradaa)
        call bh2dtaevalg_vec(nd,rscale2,c3,local2cd,nterms,ztrg,ntrg,
     1      pot,grada,gradaa)
        call bherrprintvec(nd,pot,opot,grada,ograda,gradaa,ogradaa)
        call prin2('Via mpole calculation, potential is*',pot,2*nd)
        call prin2('Via mpole calculation, grada is*',grada,2*nd*ng)
        call prin2('Via mpole calculation, gradaa is*',gradaa,2*nd*ng)
        deallocate(mpolecd)
        deallocate(mpole2cd)
        deallocate(localcd)
        deallocate(local2cd)


 8000 continue

      stop
      end
c
c
c

      subroutine bherrprintvec(nd,pot,opot,grada,ograda,gradaa,ogradaa)
      implicit real *8 (a-h,o-z)
      complex *16 pot(nd),opot(nd)
      complex *16 grada(nd),ograda(nd),gradaa(nd),ogradaa(nd)
c
cc	call prin2('error in potential is*',
cc     2 		abs(pot-opot),1)
cc	call prin2('rel error in potential is*',
cc     2 		abs(pot-opot)/abs(opot),1)
      do ii = 1,nd
        ferr = abs(grada(ii)-ograda(ii))
        fddd = abs(ograda(ii))
        herr = abs(gradaa(ii)-ogradaa(ii))
        hddd = abs(ogradaa(ii))
        write(6,*) 
     1     'POT err','     rel err','      GRADA err','     rel err',
     1     '    GRADAA err','   rel err'
        write(6,1000) abs(pot(ii)-opot(ii)),
     1 abs(pot(ii)-opot(ii))/abs(opot(ii)),ferr,ferr/fddd,herr,herr/hddd
        write(13,*) 
     1     'POT err','     rel err','      GRADA err','     rel err',
     1     '    GRADAA err','   rel err'
        write(13,1000) abs(pot(ii)-opot(ii)),
     1 abs(pot(ii)-opot(ii))/abs(opot(ii)),ferr,ferr/fddd,herr,herr/hddd
1000  format(6D12.5) 
      enddo
      return
      end
c
C
C
      subroutine mksource(nd,ns,source,charge,dip1,dip2,dummy,c0)
      implicit none
      integer nd,ns
      real* 8 dd,source(2,ns),c0(2)
      complex* 16 dip1(nd,ns),charge(nd,ns),dummy(nd,ns),dip2(nd,ns)
      complex *16 eye
c
      data eye/(0.0d0,1.0d0)/
c
      dd = 1.0d-6
      source(1,1)=c0(1)+0.24d0
      source(2,1)=c0(2)+0.25d0
      charge(1,1)= (1.0d0 + eye*0.2d0)
      dummy(1,1)= 0.0d0
      dip1(1,1)= charge(1,1)*dd
      dip1(1,1)= charge(1,1)
      dip2(1,1)= charge(1,1)*dd*eye
      dip2(1,1)= charge(1,1)*eye

      charge(2,1)= (1.0d0 + eye*0.2d0)
      dummy(2,1)= 0.0d0
      dip1(2,1)= charge(1,1)*dd
      dip1(2,1)= charge(1,1)
      dip2(2,1) = charge(1,1)*dd
      dip2(2,1) = charge(1,1)
      source(1,2)=c0(1)+0.25d0
      source(2,2)=c0(2)-0.25d0
      source(1,2)=source(1,1)+dd
      source(2,2)=source(2,1)
c
      charge(1,2)= -1.0d0 + eye*0.1d0
      charge(2,2)= -1.0d0 + eye*0.1d0
      dummy(1,2)= 0.0d0
      dummy(2,2)= 0.0d0
      dip1(1,2)= 2.0d0
      dip1(2,2)= 2.0d0
      dip2(1,2)= 2.0d0 + 3.0d0*eye
      dip2(2,2)= 2.0d0 + 3.0d0*eye

      return
      end
           
      subroutine zero_out(nd,ns,zarray)
      implicit none
      integer nd,ns,i,j
      complex* 16 zarray(nd,ns)
c
      do i = 1,nd
      do j = 1,ns
          zarray(i,j) = 0.0d0
      enddo
      enddo
      return
      end
           

