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
c   Laplace FMM in R^2: evaluate all pairwise particle
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

      subroutine lfmm2d_st_c_p_vec(nd,eps,ns,sources,
     1            charge,pot,nt,targ,pottarg)
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
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns)
      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 1
      ifpghtarg = 1

      call lfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------


      subroutine lfmm2d_st_c_g_vec(nd,eps,ns,sources,
     1            charge,pot,grad,nt,targ,pottarg,gradtarg)
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
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt)

c
cc     temporary variables
c
      complex *16 dipstr(nd)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 1
      ifdipole = 0
      
      ifpgh = 2
      ifpghtarg = 2

      call lfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine lfmm2d_st_c_h_vec(nd,eps,ns,sources,
     1            charge,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg)
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
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 charge(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns),hess(nd,3,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt),hesstarg(nd,3,nt)

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

      call lfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end


c-------------------------------      

      subroutine lfmm2d_st_d_p_vec(nd,eps,ns,sources,
     1            dipstr,pot,nt,targ,pottarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
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
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns)
      complex *16 pottarg(nd,nt)

c
cc     temporary variables
c
      complex *16 charge(nd)
      complex *16 grad(nd,2),gradtarg(nd,2)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 1
      ifpghtarg = 1

      call lfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------


      subroutine lfmm2d_st_d_g_vec(nd,eps,ns,sources,
     1            dipstr,pot,grad,nt,targ,pottarg,gradtarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
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
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt)

c
cc     temporary variables
c
      complex *16 charge(nd)
      complex *16 hess(nd,3),hesstarg(nd,3)
      integer ifcharge,ifdipole,iper
      integer ifpgh,ifpghtarg

      ifcharge = 0
      ifdipole = 1
      
      ifpgh = 2
      ifpghtarg = 2

      call lfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine lfmm2d_st_d_h_vec(nd,eps,ns,sources,
     1            dipstr,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg)
c----------------------------------------------
c   INPUT PARAMETERS:
c   nd            : number of expansions
c   eps           : FMM precision requested
c   ns            : number of sources
c   sources(2,ns) : source locations
c   dipstr(nd,ns)    : dipole strengths
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
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
      complex *16 dipstr(nd,ns)

      complex *16 pot(nd,ns),grad(nd,2,ns),hess(nd,3,ns)
      complex *16 pottarg(nd,nt),gradtarg(nd,2,nt),hesstarg(nd,3,nt)

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

      call lfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end


c-------------------------------      

      subroutine lfmm2d_st_cd_p_vec(nd,eps,ns,sources,charge,
     1            dipstr,pot,nt,targ,pottarg)
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
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
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

      call lfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------


      subroutine lfmm2d_st_cd_g_vec(nd,eps,ns,sources,charge,
     1            dipstr,pot,grad,nt,targ,pottarg,gradtarg)
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
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
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

      call lfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end
c------------------------------
c
c
c
c
c
      subroutine lfmm2d_st_cd_h_vec(nd,eps,ns,sources,charge,
     1            dipstr,pot,grad,hess,nt,targ,pottarg,
     2            gradtarg,hesstarg)
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
      integer ns,nt
      real *8 sources(2,ns),targ(2,nt)
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

      call lfmm2d(nd,eps,ns,sources,ifcharge,charge,
     1            ifdipole,dipstr,iper,ifpgh,pot,grad,hess,
     2            nt,targ,ifpghtarg,pottarg,gradtarg,
     3            hesstarg)
      return
      end


c-------------------------------      
