      subroutine lndiv2d(eps,ns,nt,ifcharge,ifdipole,ifpgh,
     1   ifpghtarg,ndiv,idivflag)
c
c
c       this subroutine estimates ndiv and idivflag 
c       based on geometry parameters
c       
c
      implicit none
      real *8 eps
      integer ns,nt,ifcharge,ifdipole,ifpgh,ifpghtarg,ndiv
      integer idivflag

      idivflag = 0


       if(eps.ge.0.5d-0) then
         ndiv = 3
       else if(eps.ge.0.5d-1) then
         ndiv = 5
       else if(eps.ge.0.5d-2) then
         ndiv = 8
       else if(eps.ge.0.5d-3) then
         ndiv = 10
       else if(eps.ge.0.5d-6) then
         ndiv = 15
       else if(eps.ge.0.5d-9) then
         ndiv = 20
       else if(eps.ge.0.5d-12) then
         ndiv = 25
       else if(eps.ge.0.5d-15) then
         ndiv = 45
       else
         ndiv = ns+nt
       endif


      return
      end
