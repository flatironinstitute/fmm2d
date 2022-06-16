cc Copyright (C) 2010-2011: Leslie Greengard, Zydrunas Gimbustas 
cc and Manas Rachh
cc Contact: greengard@cims.nyu.edu
cc
cc      
cc 2021/07/07: create Cauchy wrappers (search and replace)
cc      Travis Askham     
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
c   Cauchy FMM in R^2: evaluate all pairwise particle
c   interactions (ignoring self-interaction) 
c   and interactions with targets.
c
c   We use log(r) for the Green's function.
c   Self-interactions are not included
c
c   l2d: charge and dipstr are complex valued, x in \R^2
c
c   \phi(x_i) = \sum_{j\ne i} charge_j log|x_i-x_j}
c   + dipstr_j/x_i - x_j
c
c

      subroutine cfmm2d_s_c_p_vec(nd,eps,ns,sources,
     1            charge,pot,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge
cf2py  intent(out) pot,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
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
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      complex *16 grad(nd),gradtarg(nd)
      complex *16 hess(nd),hesstarg(nd)
      real *8 targ(2)
      complex *16 pottarg(nd,1)
      integer nt
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 0

      nt = 0

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine cfmm2d_s_c_g_vec(nd,eps,ns,sources,
     1            charge,pot,grad,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge
cf2py  intent(out) pot,grad,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,ns)    : gradients at the source locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns),grad(nd,ns)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      complex *16 hess(nd),hesstarg(nd)
      real *8 targ(2)
      complex *16 pottarg(nd,1),gradtarg(nd,1)
      integer nt
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 0
      
      nt = 0

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
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
      subroutine cfmm2d_s_c_h_vec(nd,eps,ns,sources,
     1            charge,pot,grad,hess,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge
cf2py  intent(out) pot,grad,hess,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,ns)    : gradients at the source locations
c   hess(nd,ns)    : hessian at the source locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns),grad(nd,ns),hess(nd,ns)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      real *8 targ(2)
      complex *16 pottarg(nd,1),gradtarg(nd,1),hesstarg(nd,1)
      integer nt
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 3
      ifpghtarg = 0

      nt = 0

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine cfmm2d_s_d_p_vec(nd,eps,ns,sources,
     1            dipstr,pot,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,dipstr
cf2py  intent(out) pot,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns) : dipole strengths
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
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns)

c
cc     temporary variables
c
      complex *16 charge(nd)
      complex *16 grad(nd),gradtarg(nd)
      complex *16 hess(nd),hesstarg(nd)
      complex *16 pottarg(nd,1)
      real *8 targ(2)
      integer nt
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      nt = 0

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine cfmm2d_s_d_g_vec(nd,eps,ns,sources,
     1            dipstr,pot,grad,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,dipstr
cf2py  intent(out) pot,grad,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns) : dipole strengths
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,ns)    : gradients at the source locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,ns)

c
cc     temporary variables
c
      complex *16 charge(nd)
      complex *16 hess(nd),hesstarg(nd)
      complex *16 pottarg(nd,1),gradtarg(nd,1)
      real *8 targ(2)
      integer nt
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      nt = 0

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
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
      subroutine cfmm2d_s_d_h_vec(nd,eps,ns,sources,
     1            dipstr,pot,grad,hess,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,dipstr
cf2py  intent(out) pot,grad,hess,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns) : dipole strengths
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,ns)    : gradients at the source locations
c   hess(nd,ns)    : hessian at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 dipstr(nd,ns)
      complex *16 pot(nd,ns),grad(nd,ns),hess(nd,ns)

c
cc     temporary variables
c
      complex *16 charge(nd)
      complex *16 pottarg(nd,1),gradtarg(nd,1),hesstarg(nd,1)
      real *8 targ(2)
      integer nt
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 0

      nt = 0

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine cfmm2d_s_cd_p_vec(nd,eps,ns,sources,charge,
     1            dipstr,pot,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,dipstr
cf2py  intent(out) pot,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
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
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 charge(nd,ns),dipstr(nd,ns)
      complex *16 pot(nd,ns)

c
cc     temporary variables
c
      complex *16 grad(nd),gradtarg(nd)
      complex *16 hess(nd),hesstarg(nd)
      complex *16 pottarg(nd,1)
      real *8 targ(2)
      integer nt
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 0

      nt = 0

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine cfmm2d_s_cd_g_vec(nd,eps,ns,sources,charge,
     1            dipstr,pot,grad,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,dipstr
cf2py  intent(out) pot,grad,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,ns)    : gradients at the source locations
c

      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 charge(nd,ns),dipstr(nd,ns)
      complex *16 pot(nd,ns),grad(nd,ns)

c
cc     temporary variables
c
      complex *16 hess(nd),hesstarg(nd)
      complex *16 pottarg(nd,1),gradtarg(nd,1)
      real *8 targ(2)
      integer nt
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 0

      nt = 0

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
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
      subroutine cfmm2d_s_cd_h_vec(nd,eps,ns,sources,charge,
     1            dipstr,pot,grad,hess,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,dipstr
cf2py  intent(out) pot,grad,hess,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,ns)    : gradients at the source locations
c   hess(nd,ns)    : hessian at the source locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      integer ns,ier
      real *8 sources(2,ns)
      complex *16 charge(nd,ns),dipstr(nd,ns)
      complex *16 pot(nd,ns),grad(nd,ns),hess(nd,ns)

c
cc     temporary variables
c
      real *8 targ(2)
      complex *16 pottarg(nd,1),gradtarg(nd,1),hesstarg(nd,1)
      integer nt
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 0

      nt = 0

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c-------------------------------      
c
c
c
c      
      subroutine cfmm2d_t_c_p_vec(nd,eps,ns,sources,
     1            charge,nt,targ,pottarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pottarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
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
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      complex *16 pot(nd)
      complex *16 grad(nd),gradtarg(nd)
      complex *16 hess(nd),hesstarg(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 1

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine cfmm2d_t_c_g_vec(nd,eps,ns,sources,
     1            charge,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pottarg,gradtarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pottarg(nd,nt),gradtarg(nd,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      complex *16 hess(nd),hesstarg(nd)
      complex *16 pot(nd),grad(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 2

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
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
      subroutine cfmm2d_t_c_h_vec(nd,eps,ns,sources,
     1            charge,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,nt): gradient at the target locations
c   hesstarg(nd,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,nt),hesstarg(nd,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      complex *16 pot(nd),grad(nd),hess(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 0
      ifpghtarg = 3

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine cfmm2d_t_d_p_vec(nd,eps,ns,sources,
     1            dipstr,nt,targ,pottarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,dipstr,nt,targ
cf2py  intent(out) pottarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns) : dipole strengths
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
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(nd,ns)

      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 charge(nd)
      complex *16 grad(nd),gradtarg(nd)
      complex *16 hess(nd),hesstarg(nd)
      complex *16 pot(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine cfmm2d_t_d_g_vec(nd,eps,ns,sources,
     1            dipstr,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,dipstr,nt,targ
cf2py  intent(out) pottarg,gradtarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns) : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(nd,ns)

      complex *16 pottarg(nd,nt),gradtarg(nd,nt)

c
cc     temporary variables
c
      complex *16 charge(nd)
      complex *16 hess(nd),hesstarg(nd)
      complex *16 pot(nd),grad(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
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
      subroutine cfmm2d_t_d_h_vec(nd,eps,ns,sources,
     1            dipstr,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,dipstr,nt,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns) : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,nt): gradient at the target locations
c   hesstarg(nd,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(nd,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,nt),hesstarg(nd,nt)

c
cc     temporary variables
c
      complex *16 charge(nd)
      complex *16 pot(nd),grad(nd),hess(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 3

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine cfmm2d_t_cd_p_vec(nd,eps,ns,sources,charge,
     1            dipstr,nt,targ,pottarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,dipstr,nt,targ
cf2py  intent(out) pottarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
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
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns),dipstr(nd,ns)
      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 grad(nd),gradtarg(nd)
      complex *16 hess(nd),hesstarg(nd)
      complex *16 pot(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 1

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine cfmm2d_t_cd_g_vec(nd,eps,ns,sources,charge,
     1            dipstr,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,dipstr,nt,targ
cf2py  intent(out) pottarg,gradtarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,nt): gradient at the target locations
c

      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns),dipstr(nd,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,nt)

c
cc     temporary variables
c
      complex *16 hess(nd),hesstarg(nd)
      complex *16 pot(nd),grad(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 2

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
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
      subroutine cfmm2d_t_cd_h_vec(nd,eps,ns,sources,charge,
     1            dipstr,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,dipstr,nt,targ
cf2py  intent(out) pottarg,gradtarg,hesstarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,nt): gradient at the target locations
c   hesstarg(nd,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns),dipstr(nd,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,nt),hesstarg(nd,nt)

c
cc     temporary variables
c
      complex *16 pot(nd),grad(nd),hess(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 0
      ifpghtarg = 3

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      
c--------------------------------------------
c
c
c
c      
      subroutine cfmm2d_st_c_p_vec(nd,eps,ns,sources,
     1            charge,pot,nt,targ,pottarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pot,pottarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
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
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns)
      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      complex *16 grad(nd),gradtarg(nd)
      complex *16 hess(nd),hesstarg(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine cfmm2d_st_c_g_vec(nd,eps,ns,sources,
     1            charge,pot,grad,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pot,grad,pottarg,gradtarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,ns)    : gradients at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns),grad(nd,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      complex *16 hess(nd),hesstarg(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
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
      subroutine cfmm2d_st_c_h_vec(nd,eps,ns,sources,
     1            charge,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,nt,targ
cf2py  intent(out) pot,grad,hess,pottarg,gradtarg,hesstarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,ns)    : gradients at the source locations
c   hess(nd,ns)    : hessian at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,nt): gradient at the target locations
c   hesstarg(nd,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns),grad(nd,ns),hess(nd,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,nt),hesstarg(nd,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 3
      ifpghtarg = 3

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine cfmm2d_st_d_p_vec(nd,eps,ns,sources,
     1            dipstr,pot,nt,targ,pottarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,dipstr,nt,targ
cf2py  intent(out) pot,pottarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns) : dipole strengths
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
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns)
      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 charge(nd)
      complex *16 grad(nd),gradtarg(nd)
      complex *16 hess(nd),hesstarg(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine cfmm2d_st_d_g_vec(nd,eps,ns,sources,
     1            dipstr,pot,grad,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,dipstr,nt,targ
cf2py  intent(out) pot,grad,pottarg,gradtarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns) : dipole strengths
      
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,ns)    : gradients at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,nt): gradient at the target locations
c


      implicit none
c
cc      calling sequence variables
c
      integer nd
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,nt)

c
cc     temporary variables
c
      complex *16 charge(nd)
      complex *16 hess(nd),hesstarg(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
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
      subroutine cfmm2d_st_d_h_vec(nd,eps,ns,sources,
     1            dipstr,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,dipstr,nt,targ
cf2py  intent(out) pot,grad,hess,pottarg,gradtarg,hesstarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns) : dipole strengths
      
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,ns)    : gradients at the source locations
c   hess(nd,ns)    : hessian at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,nt): gradient at the target locations
c   hesstarg(nd,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(nd,ns)
      complex *16 pot(nd,ns),grad(nd,ns),hess(nd,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,nt),hesstarg(nd,nt)

c
cc     temporary variables
c
      complex *16 charge(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      

      subroutine cfmm2d_st_cd_p_vec(nd,eps,ns,sources,charge,
     1            dipstr,pot,nt,targ,pottarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,dipstr,nt,targ
cf2py  intent(out) pot,pottarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
      
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
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns),dipstr(nd,ns)
      complex *16 pot(nd,ns)
      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 grad(nd),gradtarg(nd)
      complex *16 hess(nd),hesstarg(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end
c------------------------------


      subroutine cfmm2d_st_cd_g_vec(nd,eps,ns,sources,charge,
     1            dipstr,pot,grad,nt,targ,pottarg,gradtarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,dipstr,nt,targ
cf2py  intent(out) pot,grad,pottarg,gradtarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
      
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,ns)    : gradients at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,nt): gradient at the target locations
c

      implicit none
c
cc      calling sequence variables
c 
      integer nd
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns),dipstr(nd,ns)
      complex *16 pot(nd,ns),grad(nd,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,nt)

c
cc     temporary variables
c
      complex *16 hess(nd),hesstarg(nd)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
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
      subroutine cfmm2d_st_cd_h_vec(nd,eps,ns,sources,charge,
     1            dipstr,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg,ier)
cf2py  intent(in) eps,nd
cf2py  intent(in) ns,sources,charge,dipstr,nt,targ
cf2py  intent(out) pot,grad,hess,pottarg,gradtarg,hesstarg,ier
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   charge(nd,ns)    : charge strengths
c   dipstr(nd,ns)    : dipole strengths
      
c   nt            : number of targets
c   targ(2,nt)    : target locations
c
c   OUTPUT PARAMETERS
c   pot(nd,ns)       : potential at the source locations
c   grad(nd,ns)    : gradients at the source locations
c   hess(nd,ns)    : hessian at the source locations
c   pottarg(nd,nt)   : potential at the target locations
c   gradtarg(nd,nt): gradient at the target locations
c   hesstarg(nd,nt): hessian at the target locations
c


      implicit none
c
cc      calling sequence variables
c  
      integer nd
      real *8 eps
      integer ns,nt,ier
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns),dipstr(nd,ns)
      complex *16 pot(nd,ns),grad(nd,ns),hess(nd,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,nt),hesstarg(nd,nt)

c
cc     temporary variables
c
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 1
      
      ifpgh = 3
      ifpghtarg = 3

      call cfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg,ier)
      return
      end


c-------------------------------      
