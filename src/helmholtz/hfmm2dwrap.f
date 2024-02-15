cc Copyright (C) 2010-2011: Leslie Greengard, Zydrunas Gimbustas 
cc and Manas Rachh
cc Contact: greengard@cims.nyu.edu
cc 
cc This program is free software; you can redistribute it and/or modify 
cc it under the terms of the GNU General Public License as published by 
cc the Free Software Foundation; either version 2 of the License, or 
cc (at your option) any later version.  This program is distributed in 
cc the hope that it will be useful, but WITHOUT ANY WARRANTY; without 
cc even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
cc PARTICULAR PURPOSE.  See the GNU General Public License for more 
cc details. You should have received a copy of the GNU General Public 
cc License along with this program; 
cc if not, see <http://www.gnu.org/licenses/>.



c       
c   Helmholtz FMM in R^2: evaluate all pairwise particle
c   interactions (ignoring self-interaction) 
c   and interactions with targets.
c
c   We use H_0(kr)*(i/4) for the Green's function.
c   Self-interactions are not included
c
c   h2d: charge and dipstr are complex valued, x in \R^2
c
c   \phi(x_i) = (i/4)\sum_{j\ne i} charge_j H^{(1)}_0(k |x_i - x_j|)
c   + dipstr_j (dipvec_j \dot (x_i - x_j)) H^{(1)}_1(k |x_i - x_j|*
c                                          k/|x_i-x_j|
c
c

c
c
c
c
c
c

      subroutine hfmm2d_s_c_p(eps,zk,ns,sources,
     1            charge,pot,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge
cf2py  intent(out) pot,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 charge(ns)

      complex *16 pot(ns)

c
cc     temporary variables
c
      complex *16 dipstr(1)
      real *8 dipvec(2,1)
      integer nt
      real *8 targ(2)
      complex *16 pottarg(1)
      complex *16 grad(2),gradtarg(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg


      integer nd

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 0

      nt = 0

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_s_c_g(eps,zk,ns,sources,
     1            charge,pot,grad,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge
cf2py  intent(out) pot,grad,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)    : gradients at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 charge(ns)

      complex *16 pot(ns),grad(2,ns)

c
cc     temporary variables
c
      complex *16 dipstr(1)
      real *8 dipvec(2,1)
      integer nt
      real *8 targ(2)
      complex *16 pottarg(1)
      complex *16 gradtarg(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 0

      nt = 0
      
      ifpgh = 2
      ifpghtarg = 0

      nd = 1

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
      subroutine hfmm2d_s_c_h(eps,zk,ns,sources,
     1            charge,pot,grad,hess,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge
cf2py  intent(out) pot,grad,hess,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)    : gradients at the source locations
c   hess(3,ns)    : hessian at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 charge(ns)

      complex *16 pot(ns),grad(2,ns),hess(3,ns)

c
cc     temporary variables
c
      complex *16 dipstr(1)
      real *8 dipvec(2,1)
      integer nt
      real *8 targ(2)
      complex *16 pottarg(1)
      complex *16 gradtarg(2)
      complex *16 hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 3
      ifpghtarg = 0

      nd = 1

      nt = 0

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_s_d_p(eps,zk,ns,sources,
     1            dipstr,dipvec,pot,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec
cf2py  intent(out) pot,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipvec(2,ns)
      complex *16 dipstr(ns)

      complex *16 pot(ns)

c
cc     temporary variables
c
      complex *16 charge(1)
      integer nt
      real *8 targ(2)
      complex *16 pottarg(1)
      complex *16 grad(2),gradtarg(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      integer nd

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0
      
      nt = 0

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_s_d_g(eps,zk,ns,sources,
     1            dipstr,dipvec,pot,grad,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec
cf2py  intent(out) pot,grad,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)    : gradients at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipvec(2,ns)
      complex *16 dipstr(ns)

      complex *16 pot(ns),grad(2,ns)

c
cc     temporary variables
c
      complex *16 charge(1)
      integer nt
      real *8 targ(2)
      complex *16 pottarg(1)
      complex *16 gradtarg(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      nt = 0

      nd = 1

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
      subroutine hfmm2d_s_d_h(eps,zk,ns,sources,
     1            dipstr,dipvec,pot,grad,hess,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec
cf2py  intent(out) pot,grad,hess,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)    : gradients at the source locations
c   hess(3,ns)    : hessian at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipvec(2,ns)
      complex *16 dipstr(ns)

      complex *16 pot(ns),grad(2,ns),hess(3,ns)

c
cc     temporary variables
c
      complex *16 charge(1)
      integer nt
      real *8 targ(2)
      complex *16 pottarg(1)
      complex *16 gradtarg(2)
      complex *16 hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 0

      nt = 0

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_s_cd_p(eps,zk,ns,sources,charge,
     1            dipstr,dipvec,pot,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,charge
cf2py  intent(out) pot,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipvec(2,ns)
      complex *16 charge(ns),dipstr(ns)

      complex *16 pot(ns)

c
cc     temporary variables
c
      integer nt
      real *8 targ(2)
      complex *16 pottarg(1)
      complex *16 grad(2),gradtarg(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      integer nd

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      nt = 0

      nd = 1

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

      subroutine hfmm2d_s_cd_g(eps,zk,ns,sources,charge,
     1            dipstr,dipvec,pot,grad,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,charge
cf2py  intent(out) pot,grad,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)    : gradients at the source locations
c

      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipvec(2,ns)
      complex *16 charge(ns),dipstr(ns)

      complex *16 pot(ns),grad(2,ns)

c
cc     temporary variables
c
      integer nt
      real *8 targ(2)
      complex *16 pottarg(1)
      complex *16 gradtarg(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      nt = 0

      nd = 1

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
      subroutine hfmm2d_s_cd_h(eps,zk,ns,sources,charge,
     1            dipstr,dipvec,pot,grad,hess,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,charge
cf2py  intent(out) pot,grad,hess,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)    : gradients at the source locations
c   hess(3,ns)    : hessian at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipvec(2,ns)
      complex *16 charge(ns),dipstr(ns)

      complex *16 pot(ns),grad(2,ns),hess(3,ns)

c
cc     temporary variables
c
      integer nt
      real *8 targ(2)
      complex *16 pottarg(1)
      complex *16 gradtarg(2)
      complex *16 hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 0

      nt = 0

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


      subroutine hfmm2d_t_c_p(eps,zk,ns,sources,
     1            charge,nt,targ,pottarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pottarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(ns)

      complex *16 pottarg(nt)

c
cc     temporary variables
c
      complex *16 dipstr(1)
      real *8 dipvec(2,1)
      complex *16 pot(1)
      complex *16 grad(2),gradtarg(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      integer nd

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 1

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_t_c_g(eps,zk,ns,sources,
     1            charge,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pottarg,gradtarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(ns)

      complex *16 pottarg(nt),gradtarg(2,nt)

c
cc     temporary variables
c
      complex *16 dipstr(1)
      real *8 dipvec(2,1)
      complex *16 pot(1)
      complex *16 grad(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 2

      nd = 1

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
      subroutine hfmm2d_t_c_h(eps,zk,ns,sources,
     1            charge,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt): gradient at the target locations
c   hesstarg(3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(ns)

      complex *16 pottarg(nt),gradtarg(2,nt),hesstarg(3,nt)

c
cc     temporary variables
c
      complex *16 dipstr(1)
      real *8 dipvec(2,1)
      complex *16 pot(1)
      complex *16 grad(2)
      complex *16 hess(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 3

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_t_d_p(eps,zk,ns,sources,
     1            dipstr,dipvec,nt,targ,pottarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(2,ns)
      complex *16 dipstr(ns)

      complex *16 pottarg(nt)

c
cc     temporary variables
c
      complex *16 charge(1)
      complex *16 pot(1)
      complex *16 grad(2),gradtarg(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      integer nd

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_t_d_g(eps,zk,ns,sources,
     1            dipstr,dipvec,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,gradtarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(2,ns)
      complex *16 dipstr(ns)

      complex *16 pottarg(nt),gradtarg(2,nt)

c
cc     temporary variables
c
      complex *16 charge(1)
      complex *16 pot(1)
      complex *16 grad(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      nd = 1

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
      subroutine hfmm2d_t_d_h(eps,zk,ns,sources,
     1            dipstr,dipvec,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt): gradient at the target locations
c   hesstarg(3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(2,ns)
      complex *16 dipstr(ns)

      complex *16 pottarg(nt),gradtarg(2,nt),hesstarg(3,nt)

c
cc     temporary variables
c
      complex *16 charge(1)
      complex *16 pot(1)
      complex *16 grad(2)
      complex *16 hess(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 3

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_t_cd_p(eps,zk,ns,sources,charge,
     1            dipstr,dipvec,nt,targ,pottarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(2,ns)
      complex *16 charge(ns),dipstr(ns)

      complex *16 pottarg(nt)

c
cc     temporary variables
c
      complex *16 pot(1)
      complex *16 grad(2),gradtarg(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      integer nd

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      nd = 1

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

      subroutine hfmm2d_t_cd_g(eps,zk,ns,sources,charge,
     1            dipstr,dipvec,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,gradtarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt): gradient at the target locations
c

      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(2,ns)
      complex *16 charge(ns),dipstr(ns)

      complex *16 pottarg(nt),gradtarg(2,nt)

c
cc     temporary variables
c
      complex *16 pot(1)
      complex *16 grad(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      nd = 1

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
      subroutine hfmm2d_t_cd_h(eps,zk,ns,sources,charge,
     1            dipstr,dipvec,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt): gradient at the target locations
c   hesstarg(3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(2,ns)
      complex *16 charge(ns),dipstr(ns)

      complex *16 pottarg(nt),gradtarg(2,nt),hesstarg(3,nt)

c
cc     temporary variables
c
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd
      complex *16 pot(1)
      complex *16 grad(2)
      complex *16 hess(3)

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 3

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c
c
c
c
c
c

      subroutine hfmm2d_st_c_p(eps,zk,ns,sources,
     1            charge,pot,nt,targ,pottarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pot,pottarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   pottarg(nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(ns)

      complex *16 pot(ns)
      complex *16 pottarg(nt)

c
cc     temporary variables
c
      complex *16 dipstr(1)
      real *8 dipvec(2,1)
      complex *16 grad(2),gradtarg(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      integer nd

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_st_c_g(eps,zk,ns,sources,
     1            charge,pot,grad,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pot,pottarg,grad,gradtarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)    : gradients at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(ns)

      complex *16 pot(ns),grad(2,ns)
      complex *16 pottarg(nt),gradtarg(2,nt)

c
cc     temporary variables
c
      complex *16 dipstr(1)
      real *8 dipvec(2,1)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      nd = 1

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
      subroutine hfmm2d_st_c_h(eps,zk,ns,sources,
     1            charge,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pot,pottarg,grad,gradtarg,hess,hesstarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)    : gradients at the source locations
c   hess(3,ns)    : hessian at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt): gradient at the target locations
c   hesstarg(3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(ns)

      complex *16 pot(ns),grad(2,ns),hess(3,ns)
      complex *16 pottarg(nt),gradtarg(2,nt),hesstarg(3,nt)

c
cc     temporary variables
c
      complex *16 dipstr(1)
      real *8 dipvec(2,1)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 3
      ifpghtarg = 3

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_st_d_p(eps,zk,ns,sources,
     1            dipstr,dipvec,pot,nt,targ,pottarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   pottarg(nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(2,ns)
      complex *16 dipstr(ns)

      complex *16 pot(ns)
      complex *16 pottarg(nt)

c
cc     temporary variables
c
      complex *16 charge(1)
      complex *16 grad(2),gradtarg(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      integer nd

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine hfmm2d_st_d_g(eps,zk,ns,sources,
     1            dipstr,dipvec,pot,grad,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,grad,gradtarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)    : gradients at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(2,ns)
      complex *16 dipstr(ns)

      complex *16 pot(ns),grad(2,ns)
      complex *16 pottarg(nt),gradtarg(2,nt)

c
cc     temporary variables
c
      complex *16 charge(1)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      nd = 1

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
      subroutine hfmm2d_st_d_h(eps,zk,ns,sources,
     1            dipstr,dipvec,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,grad,gradtarg,hess,hesstarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)    : gradients at the source locations
c   hess(3,ns)    : hessian at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt): gradient at the target locations
c   hesstarg(3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(2,ns)
      complex *16 dipstr(ns)

      complex *16 pot(ns),grad(2,ns),hess(3,ns)
      complex *16 pottarg(nt),gradtarg(2,nt),hesstarg(3,nt)

c
cc     temporary variables
c
      complex *16 charge(1)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine hfmm2d_st_cd_p(eps,zk,ns,sources,charge,
     1            dipstr,dipvec,pot,nt,targ,pottarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   pottarg(nt)   : potential at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(2,ns)
      complex *16 charge(ns),dipstr(ns)

      complex *16 pot(ns)
      complex *16 pottarg(nt)

c
cc     temporary variables
c
      complex *16 grad(2),gradtarg(2)
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      integer nd

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      nd = 1

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

      subroutine hfmm2d_st_cd_g(eps,zk,ns,sources,charge,
     1            dipstr,dipvec,pot,grad,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,grad,gradtarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)    : gradients at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt): gradient at the target locations
c

      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(2,ns)
      complex *16 charge(ns),dipstr(ns)

      complex *16 pot(ns),grad(2,ns)
      complex *16 pottarg(nt),gradtarg(2,nt)

c
cc     temporary variables
c
      complex *16 hess(3),hesstarg(3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      nd = 1

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
      subroutine hfmm2d_st_cd_h(eps,zk,ns,sources,charge,
     1            dipstr,dipvec,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) zk
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,grad,gradtarg,hess,hesstarg,ier

c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   zk            : Helmholtz parameter
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : dipole orientation vectors
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)    : gradients at the source locations
c   hess(3,ns)    : hessian at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt): gradient at the target locations
c   hesstarg(3,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      complex *16 zk
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipvec(2,ns)
      complex *16 charge(ns),dipstr(ns)

      complex *16 pot(ns),grad(2,ns),hess(3,ns)
      complex *16 pottarg(nt),gradtarg(2,nt),hesstarg(3,nt)

c
cc     temporary variables
c
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg
      integer nd

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      nd = 1

      call hfmm2d(nd,eps,zk,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end

