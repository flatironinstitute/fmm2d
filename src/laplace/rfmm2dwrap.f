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
c     2021/07/08 : search/replace for real valued version
c                         (Travis Askham)


c       
c   Laplace FMM in R^2: evaluate all pairwise particle
c   interactions (ignoring self-interaction) 
c   and interactions with targets.
c
c   We use log(r) for the Green's function.
c   Self-interactions are not included
c
c   l2d: charge and dipstr are real valued, x in \R^2
c
c   \phi(x_i) = \sum_{j\ne i} charge_j log|x_i-x_j|
c   + dipstr_j/x_i - x_j
c
c

      subroutine rfmm2d_s_c_p(eps,ns,sources,
     1            charge,pot,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge
cf2py  intent(out) pot,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
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
      integer ns,ier
      real *8 sources(2,ns)
      real *8 charge(ns)

      real *8 pot(ns)

c
cc     temporary variables
c
      real *8 dipstr
      real *8 grad,gradtarg
      real *8 hess,hesstarg
      real *8 dipvec(2)
      real *8 targ(2)
      real *8 pottarg(1)
      integer nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nd,iper

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 0

      nt = 0

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1     ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2     nt,targ,ifpghtarg,pottarg,gradtarg,
     3     hesstarg,ier)
      return
      end
c------------------------------


      subroutine rfmm2d_s_c_g(eps,ns,sources,
     1            charge,pot,grad,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge
cf2py  intent(out) pot,grad,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)      : gradients at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,ier
      real *8 sources(2,ns)
      real *8 charge(ns)
      real *8 pot(ns),grad(2,ns)

c
cc     temporary variables
c
      real *8 dipstr
      real *8 hess,hesstarg
      real *8 dipvec(2)
      real *8 targ(2)
      real *8 pottarg(1),gradtarg(2)
      integer nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 0

      nt = 0

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
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
      subroutine rfmm2d_s_c_h(eps,ns,sources,
     1            charge,pot,grad,hess,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge
cf2py  intent(out) pot,grad,hess,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)      : gradients at the source locations
c   hess(3,ns)      : hessian at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,ier
      real *8 sources(2,ns)
      real *8 charge(ns)
      real *8 pot(ns),grad(2,ns),hess(3,ns)

c
cc     temporary variables
c
      real *8 dipstr
      real *8 targ(2)
      real *8 pottarg(1),gradtarg(2),hesstarg(3)
      real *8 dipvec(2)
      integer nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 3
      ifpghtarg = 0

      nt = 0

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine rfmm2d_s_d_p(eps,ns,sources,
     1            dipstr,dipvec,pot,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,dipstr,dipvec
cf2py  intent(out) pot,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : real dipole orientations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipstr(ns)
      real *8 dipvec(2,ns)

      real *8 pot(ns)

c
cc     temporary variables
c
      real *8 charge
      real *8 grad,gradtarg
      real *8 hess,hesstarg
      real *8 targ(2)
      real *8 pottarg(1)
      integer nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nd,iper

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      nt = 0

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine rfmm2d_s_d_g(eps,ns,sources,
     1            dipstr,dipvec,pot,grad,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,dipstr,dipvec
cf2py  intent(out) pot,grad,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : real dipole orientations
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
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pot(ns),grad(2,ns)

c
cc     temporary variables
c
      real *8 charge
      real *8 hess,hesstarg
      real *8 targ(2)
      real *8 pottarg(1),gradtarg(2)
      integer nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      nt = 0

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
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
      subroutine rfmm2d_s_d_h(eps,ns,sources,
     1            dipstr,dipvec,pot,grad,hess,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,dipstr,dipvec
cf2py  intent(out) pot,grad,hess,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : real dipole orientations      
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
      integer ns,ier
      real *8 sources(2,ns)
      real *8 dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pot(ns),grad(2,ns),hess(3,ns)

c
cc     temporary variables
c
      real *8 charge
      real *8 targ(2)
      real *8 pottarg(1),gradtarg(2),hesstarg(3)
      integer nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 0

      nt = 0

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine rfmm2d_s_cd_p(eps,ns,sources,charge,
     1            dipstr,dipvec,pot,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,dipstr,dipvec
cf2py  intent(out) pot,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : real dipole orientations      
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,ier
      real *8 sources(2,ns)
      real *8 charge(ns),dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pot(ns)

c
cc     temporary variables
c
      real *8 grad,gradtarg
      real *8 hess,hesstarg
      real *8 targ(2)
      real *8 pottarg(1)
      integer nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nd,iper

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      nt = 0

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine rfmm2d_s_cd_g(eps,ns,sources,charge,
     1            dipstr,dipvec,pot,grad,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,dipstr,dipvec
cf2py  intent(out) pot,grad,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : real dipole orientations      
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
      integer ns,ier
      real *8 sources(2,ns)
      real *8 charge(ns),dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pot(ns),grad(2,ns)

c
cc     temporary variables
c
      real *8 hess,hesstarg
      real *8 targ(2)
      real *8 pottarg(1),gradtarg(2)
      integer nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      nt = 0

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
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
      subroutine rfmm2d_s_cd_h(eps,ns,sources,charge,
     1            dipstr,dipvec,pot,grad,hess,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,dipstr,dipvec
cf2py  intent(out) pot,grad,hess,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : real dipole orientations      
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
      integer ns,ier
      real *8 sources(2,ns)
      real *8 charge(ns),dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pot(ns),grad(2,ns),hess(3,ns)

c
cc     temporary variables
c
      real *8 targ(2)
      real *8 pottarg(1),gradtarg(2),hesstarg(3)
      integer nt
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 0

      nt = 0

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end

c----------------------------------------------
c
c
c
      subroutine rfmm2d_t_c_p(eps,ns,sources,
     1            charge,nt,targ,pottarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pottarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
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
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(ns)

      real *8 pottarg(nt)

c
cc     temporary variables
c
      real *8 dipstr
      real *8 grad,gradtarg
      real *8 hess,hesstarg
      real *8 dipvec(2)
      real *8 pot(1)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nd,iper

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 1

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1     ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2     nt,targ,ifpghtarg,pottarg,gradtarg,
     3     hesstarg,ier)
      return
      end
c------------------------------


      subroutine rfmm2d_t_c_g(eps,ns,sources,
     1            charge,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pottarg,gradtarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)     : potential at the target locations
c   gradtarg(2,nt)  : gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(ns)
      real *8 pottarg(nt),gradtarg(2,nt)

c
cc     temporary variables
c
      real *8 dipstr
      real *8 hess,hesstarg
      real *8 pot(1),grad(2)
      real *8 dipvec(2)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 2

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
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
      subroutine rfmm2d_t_c_h(eps,ns,sources,
     1            charge,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)     : potential at the target locations
c   gradtarg(2,nt)  : gradient at the target locations
c   hesstarg(3,nt)  : hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(ns)
      real *8 pottarg(nt),gradtarg(2,nt),hesstarg(3,nt)

c
cc     temporary variables
c
      real *8 dipstr
      real *8 dipvec(2)
      real *8 pot(1),grad(2),hess(3)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 3

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine rfmm2d_t_d_p(eps,ns,sources,
     1            dipstr,dipvec,nt,targ,pottarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : real dipole orientations
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
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipstr(ns)
      real *8 dipvec(2,ns)

      real *8 pottarg(nt)

c
cc     temporary variables
c
      real *8 charge
      real *8 grad,gradtarg
      real *8 hess,hesstarg
      real *8 pot(1)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nd,iper

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine rfmm2d_t_d_g(eps,ns,sources,
     1            dipstr,dipvec,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,gradtarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : real dipole orientations
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt)  : gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pottarg(nt),gradtarg(2,nt)

c
cc     temporary variables
c
      real *8 charge
      real *8 hess,hesstarg
      real *8 pot(1),grad(2)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
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
      subroutine rfmm2d_t_d_h(eps,ns,sources,
     1            dipstr,dipvec,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : real dipole orientations      
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt)  : gradient at the target locations
c   hesstarg(3,nt)  : hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pottarg(nt),gradtarg(2,nt),hesstarg(3,nt)

c
cc     temporary variables
c
      real *8 charge
      real *8 pot(1),grad(2),hess(3)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 3

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine rfmm2d_t_cd_p(eps,ns,sources,charge,
     1            dipstr,dipvec,nt,targ,pottarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : real dipole orientations      
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
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(ns),dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pottarg(nt)

c
cc     temporary variables
c
      real *8 grad,gradtarg
      real *8 hess,hesstarg
      real *8 pot(1)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nd,iper

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine rfmm2d_t_cd_g(eps,ns,sources,charge,
     1            dipstr,dipvec,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,gradtarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : real dipole orientations      
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt)  : gradient at the target locations
c

      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(ns),dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pottarg(nt),gradtarg(2,nt)

c
cc     temporary variables
c
      real *8 hess,hesstarg
      real *8 pot(1),grad(2)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
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
      subroutine rfmm2d_t_cd_h(eps,ns,sources,charge,
     1            dipstr,dipvec,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)  : real dipole orientations      
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt)  : gradient at the target locations
c   hesstarg(3,nt)  : hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(ns),dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pottarg(nt),gradtarg(2,nt),hesstarg(3,nt)

c
cc     temporary variables
c
      real *8 pot(1),grad(2),hess(3)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 3

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c----------------------------------------------
c
c
c
      subroutine rfmm2d_st_c_p(eps,ns,sources,
     1            charge,pot,nt,targ,pottarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pot,pottarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
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
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(ns)

      real *8 pot(ns)
      real *8 pottarg(nt)

c
cc     temporary variables
c
      real *8 dipstr
      real *8 grad,gradtarg
      real *8 hess,hesstarg
      real *8 dipvec(2)
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nd,iper

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1     ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2     nt,targ,ifpghtarg,pottarg,gradtarg,
     3     hesstarg,ier)
      return
      end
c------------------------------


      subroutine rfmm2d_st_c_g(eps,ns,sources,
     1            charge,pot,grad,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pot,grad,pottarg,gradtarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(ns)      : gradients at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt)  : gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(ns)
      real *8 dipvec(2)
      real *8 pot(ns),grad(2,ns)
      real *8 pottarg(nt),gradtarg(2,nt)

c
cc     temporary variables
c
      real *8 dipstr
      real *8 hess,hesstarg
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
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
      subroutine rfmm2d_st_c_h(eps,ns,sources,
     1            charge,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pot,grad,hess,pottarg,gradtarg,hesstarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)      : gradients at the source locations
c   hess(3,ns)      : hessian at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt)  : gradient at the target locations
c   hesstarg(3,nt)  : hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(ns)
      real *8 dipvec(2)
      real *8 pot(ns),grad(2,ns),hess(3,ns)
      real *8 pottarg(nt),gradtarg(2,nt),hesstarg(3,nt)

c
cc     temporary variables
c
      real *8 dipstr
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 3
      ifpghtarg = 3

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine rfmm2d_st_d_p(eps,ns,sources,
     1            dipstr,dipvec,pot,nt,targ,pottarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c   dipvec(2,ns)    : real dipole orientations
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
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipstr(ns)
      real *8 dipvec(2,ns)

      real *8 pot(ns)
      real *8 pottarg(nt)

c
cc     temporary variables
c
      real *8 charge
      real *8 grad,gradtarg
      real *8 hess,hesstarg
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nd,iper

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine rfmm2d_st_d_g(eps,ns,sources,
     1            dipstr,dipvec,pot,grad,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,grad,pottarg,gradtarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(ns)    : dipole strengths
c     dipvec(2,ns)    : real dipole orientations
c     nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)      : gradients at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt)  : gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pot(ns),grad(2,ns)
      real *8 pottarg(nt),gradtarg(2,nt)

c
cc     temporary variables
c
      real *8 charge
      real *8 hess,hesstarg
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
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
      subroutine rfmm2d_st_d_h(eps,ns,sources,
     1            dipstr,dipvec,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,grad,hess,pottarg,gradtarg,hesstarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c     dipstr(ns)    : dipole strengths
c     dipvec(2,ns)    : real dipole orientations      
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)      : gradients at the source locations
c   hess(3,ns)      : hessian at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt)  : gradient at the target locations
c   hesstarg(3,nt)  : hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pot(ns),grad(2,ns),hess(3,ns)
      real *8 pottarg(nt),gradtarg(2,nt),hesstarg(3,nt)

c
cc     temporary variables
c
      real *8 charge
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine rfmm2d_st_cd_p(eps,ns,sources,charge,
     1            dipstr,dipvec,pot,nt,targ,pottarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,pottarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c     dipstr(ns)    : dipole strengths
c     dipvec(2,ns)    : real dipole orientations      
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
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(ns),dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pot(ns)
      real *8 pottarg(nt)

c
cc     temporary variables
c
      real *8 grad,gradtarg
      real *8 hess,hesstarg
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg

      integer nd,iper

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine rfmm2d_st_cd_g(eps,ns,sources,charge,
     1            dipstr,dipvec,pot,grad,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,grad,pottarg,gradtarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c     dipstr(ns)    : dipole strengths
c     dipvec(2,ns)    : real dipole orientations      
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)      : gradients at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt)  : gradient at the target locations
c

      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(ns),dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pot(ns),grad(2,ns)
      real *8 pottarg(nt),gradtarg(2,nt)

c
cc     temporary variables
c
      real *8 hess,hesstarg
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
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
      subroutine rfmm2d_st_cd_h(eps,ns,sources,charge,
     1            dipstr,dipvec,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps
cf2py  intent(in) ns,sources,charge,dipstr,dipvec,nt,targ
cf2py  intent(out) pot,grad,hess,pottarg,gradtarg,hesstarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(ns)    : charge strengths
c     dipstr(ns)    : dipole strengths
c     dipvec(2,ns)    : real dipole orientations      
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(ns)       : potential at the source locations
c   grad(2,ns)      : gradients at the source locations
c   hess(3,ns)      : hessian at the source locations
c   pottarg(nt)   : potential at the target locations
c   gradtarg(2,nt)  : gradient at the target locations
c   hesstarg(3,nt)  : hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      real *8 charge(ns),dipstr(ns)
      real *8 dipvec(2,ns)
      real *8 pot(ns),grad(2,ns),hess(3,ns)
      real *8 pottarg(nt),gradtarg(2,nt),hesstarg(3,nt)

c
cc     temporary variables
c
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      integer nd,iper

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      nd = 1

      call rfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,dipvec,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end

