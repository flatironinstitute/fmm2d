


c-------------------------------      

      subroutine hfmm2d_s_c_p_vec(nd,eps,zk,ns,sources,
     1            charge,pot,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge
cf2py  intent(out) pot,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns)

c
cc     temporary variables
c
      complex *16 dipstr(nd,1)
      real *8 dipvec(nd,2,1)
      integer nt
      real *8 targ(2)
      complex *16 pottarg(nd,1)
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 0

      nt = 0

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_s_c_g_vec(nd,eps,zk,ns,sources,
     1            charge,pot,grad,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge
cf2py  intent(out) pot,grad,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns)

c
cc     temporary variables
c
      complex *16 dipstr(nd,1)
      real *8 dipvec(nd,2,1)
      integer nt
      real *8 targ(2)
      complex *16 pottarg(nd,1)
      complex *16 gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 0

      nt = 0

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine hfmm2d_s_c_h_vec(nd,eps,zk,ns,sources,
     1            charge,pot,grad,hess,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge
cf2py  intent(out) pot,grad,hess,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   hess(nd,3,ns)    : hessian at the source locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns),hess(nd,3,ns)

c
cc     temporary variables
c
      complex *16 dipstr(nd,1)
      real *8 dipvec(nd,2,1)
      integer nt
      real *8 targ(2)
      complex *16 pottarg(nd,1)
      complex *16 gradtarg(nd,2)
      complex *16 hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 3
      ifpghtarg = 0
      
      nt = 0

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_s_d_p_vec(nd,eps,zk,ns,sources,
     1            dipstr,dipvec,pot,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec
cf2py  intent(out) pot,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipvec(nd,2,ns)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns)

c
cc     temporary variables
c
      complex *16 charge(nd,1)
      integer nt
      real *8 targ(2)
      complex *16 pottarg(nd,1)
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      nt = 0

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_s_d_g_vec(nd,eps,zk,ns,sources,
     1            dipstr,dipvec,pot,grad,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec
cf2py  intent(out) pot,grad,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipvec(nd,2,ns)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns)

c
cc     temporary variables
c
      complex *16 charge(nd,1)
      integer nt
      real *8 targ(2)
      complex *16 pottarg(nd,1)
      complex *16 gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      nt = 0

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine hfmm2d_s_d_h_vec(nd,eps,zk,ns,sources,
     1            dipstr,dipvec,pot,grad,hess,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec
cf2py  intent(out) pot,grad,hess,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   hess(nd,3,ns)    : hessian at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipvec(nd,2,ns)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns),hess(nd,3,ns)

c
cc     temporary variables
c
      complex *16 charge(nd,1)
      integer nt
      real *8 targ(2)
      complex *16 pottarg(nd,1)
      complex *16 gradtarg(nd,2)
      complex *16 hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 0

      nt = 0

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_s_cd_p_vec(nd,eps,zk,ns,sources,charge,
     1            dipstr,dipvec,pot,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,charge
cf2py  intent(out) pot,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipvec(nd,2,ns)
      complex *16 charge(nd,ns),dipstr(nd,ns)

      complex *16 pot(nd,ns)

c
cc     temporary variables
c
      integer nt
      real *8 targ(2)
      complex *16 pottarg(nd,1)
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      nt = 0

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_s_cd_g_vec(nd,eps,zk,ns,sources,charge,
     1            dipstr,dipvec,pot,grad,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,charge
cf2py  intent(out) pot,grad,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c

      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipvec(nd,2,ns)
      complex *16 charge(nd,ns),dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns)

c
cc     temporary variables
c
      integer nt
      real *8 targ(2)
      complex *16 pottarg(nd,1)
      complex *16 gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      nt = 0

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine hfmm2d_s_cd_h_vec(nd,eps,zk,ns,sources,charge,
     1            dipstr,dipvec,pot,grad,hess,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,charge
cf2py  intent(out) pot,grad,hess,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   hess(nd,3,ns)    : hessian at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipvec(nd,2,ns)
      complex *16 charge(nd,ns),dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns),hess(nd,3,ns)

c
cc     temporary variables
c
      integer nt
      real *8 targ(2)
      complex *16 pottarg(nd,1)
      complex *16 gradtarg(nd,2)
      complex *16 hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 0

      nt = 0

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      


      subroutine hfmm2d_t_c_p_vec(nd,eps,zk,ns,sources,
     1            charge,nt,targ,pottarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pottarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd,1)
      real *8 dipvec(nd,2,1)
      complex *16 pot(nd)
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_t_c_g_vec(nd,eps,zk,ns,sources,
     1            charge,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pottarg,gradtarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd,1)
      real *8 dipvec(nd,2,1)
      complex *16 pot(nd,1)
      complex *16 grad(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 2

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine hfmm2d_t_c_h_vec(nd,eps,zk,ns,sources,
     1            charge,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c   hesstarg(nd,3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt),hesstarg(nd,3,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd,1)
      real *8 dipvec(nd,2,1)
      complex *16 pot(nd,1)
      complex *16 grad(nd,2)
      complex *16 hess(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 3

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_t_d_p_vec(nd,eps,zk,ns,sources,
     1            dipstr,dipvec,nt,targ,pottarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,ns)
      complex *16 dipstr(nd,ns)

      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 charge(nd,1)
      complex *16 pot(nd,1)
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_t_d_g_vec(nd,eps,zk,ns,sources,
     1            dipstr,dipvec,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,gradtarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,ns)
      complex *16 dipstr(nd,ns)

      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt)

c
cc     temporary variables
c
      complex *16 charge(nd,1)
      complex *16 pot(nd,1)
      complex *16 grad(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine hfmm2d_t_d_h_vec(nd,eps,zk,ns,sources,
     1            dipstr,dipvec,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c   hesstarg(nd,3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,ns)
      complex *16 dipstr(nd,ns)

      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt),hesstarg(nd,3,nt)

c
cc     temporary variables
c
      complex *16 charge(nd,1)
      complex *16 pot(nd,1)
      complex *16 grad(nd,2)
      complex *16 hess(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 3

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_t_cd_p_vec(nd,eps,zk,ns,sources,charge,
     1            dipstr,dipvec,nt,targ,pottarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,ns)
      complex *16 charge(nd,ns),dipstr(nd,ns)

      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 pot(nd,1)
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_t_cd_g_vec(nd,eps,zk,ns,sources,charge,
     1            dipstr,dipvec,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,gradtarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c

      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,ns)
      complex *16 charge(nd,ns),dipstr(nd,ns)

      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt)

c
cc     temporary variables
c
      complex *16 pot(nd,1)
      complex *16 grad(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine hfmm2d_t_cd_h_vec(nd,eps,zk,ns,sources,charge,
     1            dipstr,dipvec,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c   hesstarg(nd,3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,ns)
      complex *16 charge(nd,ns),dipstr(nd,ns)

      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt),hesstarg(nd,3,nt)

c
cc     temporary variables
c
      complex *16 pot(nd,1)
      complex *16 grad(nd,1)
      complex *16 hess(nd,1)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 3

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_st_c_p_vec(nd,eps,zk,ns,sources,
     1            charge,pot,nt,targ,pottarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pot,pottarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns)
      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd,1)
      real *8 dipvec(nd,2,1)
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_st_c_g_vec(nd,eps,zk,ns,sources,
     1            charge,pot,grad,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pot,pottarg,grad,gradtarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd,1)
      real *8 dipvec(nd,2,1)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine hfmm2d_st_c_h_vec(nd,eps,zk,ns,sources,
     1            charge,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pot,pottarg,grad,gradtarg,hess,hesstarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   hess(nd,3,ns)    : hessian at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c   hesstarg(nd,3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns),hess(nd,3,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt),hesstarg(nd,3,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd,1)
      real *8 dipvec(nd,2,1)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 3
      ifpghtarg = 3

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_st_d_p_vec(nd,eps,zk,ns,sources,
     1            dipstr,dipvec,pot,nt,targ,pottarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,ns)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns)
      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 charge(nd,1)
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_st_d_g_vec(nd,eps,zk,ns,sources,
     1            dipstr,dipvec,pot,grad,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,grad,gradtarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,ns)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt)

c
cc     temporary variables
c
      complex *16 charge(nd,1)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine hfmm2d_st_d_h_vec(nd,eps,zk,ns,sources,
     1            dipstr,dipvec,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,grad,gradtarg,hess,hesstarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   hess(nd,3,ns)    : hessian at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c   hesstarg(nd,3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,ns)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns),hess(nd,3,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt),hesstarg(nd,3,nt)

c
cc     temporary variables
c
      complex *16 charge(nd,1)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_st_cd_p_vec(nd,eps,zk,ns,sources,charge,
     1            dipstr,dipvec,pot,nt,targ,pottarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,ns)
      complex *16 charge(nd,ns),dipstr(nd,ns)

      complex *16 pot(nd,ns)
      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_st_cd_g_vec(nd,eps,zk,ns,sources,charge,
     1            dipstr,dipvec,pot,grad,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,grad,gradtarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c

      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,ns)
      complex *16 charge(nd,ns),dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt)

c
cc     temporary variables
c
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine hfmm2d_st_cd_h_vec(nd,eps,zk,ns,sources,charge,
     1            dipstr,dipvec,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,grad,gradtarg,hess,hesstarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   dipvec(nd,2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,2,ns)    : gradients at the source locations
c   hess(nd,3,ns)    : hessian at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,2,nt): gradient at the target locations
c   hesstarg(nd,3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(nd,2,ns)
      complex *16 charge(nd,ns),dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns),hess(nd,3,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt),hesstarg(nd,3,nt)

c
cc     temporary variables
c
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      
