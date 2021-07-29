c  The subroutine names take the following form:
c    mbh2d_direct<int-ker><out>_vec
c
c      <int-ker>: kernel of interaction
c           (charges/dipoles/quadrupoles/octopoles, any combo thereof)
c        c: charges
c        d: dipoles
c	 q: quadrupoles
c	 o: octopoles
c        cd: charges+dipoles
c	 cdqo: charges+dipoles+quadrupoles+octopoles
c	 ... etc
c 
c      <out>: flag for potential/potential+gradient
c        p: potentials
c        g: potentials+gradients
c        h: potentials+gradients+hessians
c



      subroutine mbh2d_directcp_vec(nd,beta,source,ns,
     1     charge,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdp_vec(nd,beta,source,ns,
     1     dipstr,dipvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdp_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directqp_vec(nd,beta,source,ns,
     1     quastr,quavec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqp_vec(nd,beta,source,ns,
     1     charge,quastr,quavec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqp_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quastr,quavec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqp_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quastr,quavec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directop_vec(nd,beta,source,ns,
     1     octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcop_vec(nd,beta,source,ns,
     1     charge,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdop_vec(nd,beta,source,ns,
     1     dipstr,dipvec,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdop_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directqop_vec(nd,beta,source,ns,
     1     quastr,quavec,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqop_vec(nd,beta,source,ns,
     1     charge,quastr,quavec,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqop_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quastr,quavec,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqop_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quastr,quavec,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcp_vec(nd,beta,source,ns,
     1     charge,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdp_vec(nd,beta,source,ns,
     1     dipstr,dipvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdp_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directqp_vec(nd,beta,source,ns,
     1     quastr,quavec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqp_vec(nd,beta,source,ns,
     1     charge,quastr,quavec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqp_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quastr,quavec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqp_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quastr,quavec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directop_vec(nd,beta,source,ns,
     1     octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcop_vec(nd,beta,source,ns,
     1     charge,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdop_vec(nd,beta,source,ns,
     1     dipstr,dipvec,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdop_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directqop_vec(nd,beta,source,ns,
     1     quastr,quavec,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqop_vec(nd,beta,source,ns,
     1     charge,quastr,quavec,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqop_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quastr,quavec,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqop_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quastr,quavec,octstr,octvec,
     2     target,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      real *8 pot
      
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcg_vec(nd,beta,source,ns,
     1     charge,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdg_vec(nd,beta,source,ns,
     1     dipstr,dipvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdg_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directqg_vec(nd,beta,source,ns,
     1     quastr,quavec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqg_vec(nd,beta,source,ns,
     1     charge,quastr,quavec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqg_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quastr,quavec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqg_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quastr,quavec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directog_vec(nd,beta,source,ns,
     1     octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcog_vec(nd,beta,source,ns,
     1     charge,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdog_vec(nd,beta,source,ns,
     1     dipstr,dipvec,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdog_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directqog_vec(nd,beta,source,ns,
     1     quastr,quavec,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqog_vec(nd,beta,source,ns,
     1     charge,quastr,quavec,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqog_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quastr,quavec,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqog_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quastr,quavec,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcg_vec(nd,beta,source,ns,
     1     charge,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdg_vec(nd,beta,source,ns,
     1     dipstr,dipvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdg_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directqg_vec(nd,beta,source,ns,
     1     quastr,quavec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqg_vec(nd,beta,source,ns,
     1     charge,quastr,quavec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqg_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quastr,quavec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqg_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quastr,quavec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directog_vec(nd,beta,source,ns,
     1     octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcog_vec(nd,beta,source,ns,
     1     charge,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdog_vec(nd,beta,source,ns,
     1     dipstr,dipvec,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdog_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directqog_vec(nd,beta,source,ns,
     1     quastr,quavec,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqog_vec(nd,beta,source,ns,
     1     charge,quastr,quavec,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqog_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quastr,quavec,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqog_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quastr,quavec,octstr,octvec,
     2     target,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      real *8 pot,grad(2)
            
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directch_vec(nd,beta,source,ns,
     1     charge,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdh_vec(nd,beta,source,ns,
     1     dipstr,dipvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdh_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directqh_vec(nd,beta,source,ns,
     1     quastr,quavec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqh_vec(nd,beta,source,ns,
     1     charge,quastr,quavec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqh_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quastr,quavec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqh_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quastr,quavec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directoh_vec(nd,beta,source,ns,
     1     octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcoh_vec(nd,beta,source,ns,
     1     charge,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdoh_vec(nd,beta,source,ns,
     1     dipstr,dipvec,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdoh_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directqoh_vec(nd,beta,source,ns,
     1     quastr,quavec,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqoh_vec(nd,beta,source,ns,
     1     charge,quastr,quavec,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqoh_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quastr,quavec,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqoh_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quastr,quavec,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directch_vec(nd,beta,source,ns,
     1     charge,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdh_vec(nd,beta,source,ns,
     1     dipstr,dipvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdh_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directqh_vec(nd,beta,source,ns,
     1     quastr,quavec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqh_vec(nd,beta,source,ns,
     1     charge,quastr,quavec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqh_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quastr,quavec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqh_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quastr,quavec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
            
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directoh_vec(nd,beta,source,ns,
     1     octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcoh_vec(nd,beta,source,ns,
     1     charge,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdoh_vec(nd,beta,source,ns,
     1     dipstr,dipvec,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdoh_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directqoh_vec(nd,beta,source,ns,
     1     quastr,quavec,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqoh_vec(nd,beta,source,ns,
     1     charge,quastr,quavec,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqoh_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quastr,quavec,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqoh_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quastr,quavec,octstr,octvec,
     2     target,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     dipstr - real *8 (nd,ns) dipole strengths
c     dipvec - real *8 (nd,2,ns) dipole orientations
c     quastr - real *8 (nd,ns) quadrupole strengths
c     quavec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     target - real *8 (2) target location
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c
c     pot - real *8 (nd) potentials at target
c     grad - real *8 (nd,2) gradients at target
c     hess - real *8 (nd,3) Hessians at target
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(*)
      real *8 dipstr(*),dipvec(2,*)
      real *8 quadstr(*),quadvec(3,*)
      real *8 octstr(*),octvec(4,*)      
      real *8 target(2)
      
      
      real *8 pot,grad(2),hess(3)      
      real *8 thresh
      integer ns
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      
      do i = 1,ns

         xdiff=targ(1)-sources(1,i)
         ydiff=targ(2)-sources(2,i)
         rr=xdiff*xdiff+ydiff*ydiff
         if(rr.lt.thresh2) goto 1000
         
         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2        der4,ifder5,der5)

         do ii = 1,nd	 
c     charge contrib
            pot(ii) = pot(ii) + potloc*charge(i)
            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
c     dipole contrib
            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
     1           + gradloc(2)*dipvec(2,i))
            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
     1           + hessloc(2)*dipvec(2,i))
            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
     1           + hessloc(3)*dipvec(2,i))
            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
     1           + der3(2)*dipvec(2,i))
            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
     1           + der3(3)*dipvec(2,i))
            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
     1           + der3(4)*dipvec(2,i))

c     quadrupole contrib
            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
c     octopole contrib            
            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
     1           + der3(4)*octvec(4,i))
            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
     1           + der4(4)*octvec(4,i))
            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
     1           + der4(5)*octvec(4,i))
            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
     1           + der5(4)*octvec(4,i))
            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
     1           + der5(5)*octvec(4,i))
            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
     1           + der5(6)*octvec(4,i))
         enddo
	 
 1000    continue
      enddo
      
      
      return
      end




cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c Below is a basic template of what this routine should do
c
c
c
c      subroutine mbh2d_directcdqoh_vec(nd,beta,source,ns,
c     1     charge,dipstr,dipvec,quastr,quavec,octstr,octvec,
c     2     target,pot,grad,hess,thresh)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc
cc     input:
cc
cc     nd - integer, number of vectors per source
cc     beta - real *8, modified biharmonic parameter
cc     ns - integer, number of sources
cc     source - real *8 (2,ns) source locations
cc     charge - real *8 (nd,ns) charge strengths
cc     dipstr - real *8 (nd,ns) dipole strengths
cc     dipvec - real *8 (nd,2,ns) dipole orientations
cc     quastr - real *8 (nd,ns) quadrupole strengths
cc     quavec - real *8 (nd,3,ns) quadrupole orientations      
cc     octstr - real *8 (nd,ns) octopole strengths
cc     octvec - real *8 (nd,4,ns) octopole orientations
cc     target - real *8 (2) target location
cc     thresh - real *8 threshold, don't compute contribution of
cc                     charge if distance from charge to target is
cc                     less than this threshold
cc
cc     output:
cc
cc     pot - real *8 (nd) potentials at target
cc     grad - real *8 (nd,2) gradients at target
cc     hess - real *8 (nd,3) Hessians at target
cc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c      implicit none
cc     global variables
c      real *8 beta, source(2,*), charge(*), dipstr(*)
c      real *8 dipvec(2,*), quastr(*), quavec(3,*)
c      real *8 octstr(*), octvec(4,*)
c      real *8 target(2)
c      real *8 pot, grad(2), hess(3)
c      real *8 thresh
c      integer ns
cc     local variables
c      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
c      real *8 thresh2
c      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
c      integer i
c
c      thresh2=thresh*thresh
c
c      ifpotloc = 1
c      ifgradloc = 1
c      ifhessloc = 1
c      ifder3 = 1
c      ifder4 = 1
c      ifder5 = 1
c
c      do i = 1,ns
c
c         xdiff=targ(1)-sources(1,i)
c         ydiff=targ(2)-sources(2,i)
c         rr=xdiff*xdiff+ydiff*ydiff
c         if(rr.lt.thresh2) goto 1000
c         
c         call modbhgreen_all(beta,target,source(1,i),ifpotloc,potloc,
c     1        ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
c     2        der4,ifder5,der5)
c
c         do ii = 1,nd
cc     charge contrib
c            
c            pot(ii) = pot(ii) + potloc*charge(i)
c            grad(ii,1) = grad(ii,1) + gradloc(1)*charge(i)
c            grad(ii,2) = grad(ii,2) + gradloc(2)*charge(i)
c            hess(ii,1) = hess(ii,1) + hessloc(1)*charge(i)
c            hess(ii,2) = hess(ii,2) + hessloc(2)*charge(i)
c            hess(ii,3) = hess(ii,3) + hessloc(3)*charge(i)
c
cc     dipole contrib         
c            
c            pot(ii) = pot(ii) - dipstr(i)*(gradloc(1)*dipvec(1,i)
c     1           + gradloc(2)*dipvec(2,i))
c            grad(ii,1) = grad(ii,1) - dipstr(i)*(hessloc(1)*dipvec(1,i)
c     1           + hessloc(2)*dipvec(2,i))
c            grad(ii,2) = grad(ii,2) - dipstr(i)*(hessloc(2)*dipvec(1,i)
c     1           + hessloc(3)*dipvec(2,i))
c
c            hess(ii,1) = hess(ii,1) - dipstr(i)*(der3(1)*dipvec(1,i)
c     1           + der3(2)*dipvec(2,i))
c            hess(ii,2) = hess(ii,2) - dipstr(i)*(der3(2)*dipvec(1,i)
c     1           + der3(3)*dipvec(2,i))
c            hess(ii,3) = hess(ii,3) - dipstr(i)*(der3(3)*dipvec(1,i)
c     1           + der3(4)*dipvec(2,i))
c
c
cc     quadrupole contrib
c
c            pot(ii) = pot(ii) + quastr(i)*(hessloc(1)*quavec(1,i)
c     1           + hessloc(2)*quavec(2,i) + hessloc(3)*quavec(3,i))
c
c            grad(ii,1) = grad(ii,1) + quastr(i)*(der3(1)*quavec(1,i)
c     1           + der3(2)*quavec(2,i) + der3(3)*quavec(3,i))
c            grad(ii,2) = grad(ii,2) + quastr(i)*(der3(2)*quavec(1,i)
c     1           + der3(3)*quavec(2,i) + der3(4)*quavec(3,i))
c
c            hess(ii,1) = hess(ii,1) + quastr(i)*(der4(1)*quavec(1,i)
c     1           + der4(2)*quavec(2,i) + der4(3)*quavec(3,i))
c            hess(ii,2) = hess(ii,2) + quastr(i)*(der4(2)*quavec(1,i)
c     1           + der4(3)*quavec(2,i) + der4(4)*quavec(3,i))
c            hess(ii,3) = hess(ii,3) + quastr(i)*(der4(3)*quavec(1,i)
c     1           + der4(4)*quavec(2,i) + der4(5)*quavec(3,i))
c
cc     octopole contrib
c            
c            pot(ii) = pot(ii) - octstr(i)*(der3(1)*octvec(1,i)
c     1           + der3(2)*octvec(2,i) + der3(3)*octvec(3,i)
c     1           + der3(4)*octvec(4,i))
c
c            grad(ii,1) = grad(ii,1) - octstr(i)*(der4(1)*octvec(1,i)
c     1           + der4(2)*octvec(2,i) + der4(3)*octvec(3,i)
c     1           + der4(4)*octvec(4,i))
c            grad(ii,2) = grad(ii,2) - octstr(i)*(der4(2)*octvec(1,i)
c     1           + der4(3)*octvec(2,i) + der4(4)*octvec(3,i)
c     1           + der4(5)*octvec(4,i))
c
c            hess(ii,1) = hess(ii,1) - octstr(i)*(der5(1)*octvec(1,i)
c     1           + der5(2)*octvec(2,i) + der5(3)*octvec(3,i)
c     1           + der5(4)*octvec(4,i))
c            hess(ii,2) = hess(ii,2) - octstr(i)*(der5(2)*octvec(1,i)
c     1           + der5(3)*octvec(2,i) + der5(4)*octvec(3,i)
c     1           + der5(5)*octvec(4,i))
c            hess(ii,3) = hess(ii,3) - octstr(i)*(der5(3)*octvec(1,i)
c     1           + der5(4)*octvec(2,i) + der5(5)*octvec(3,i)
c     1           + der5(6)*octvec(4,i))
c            
c         enddo
c 1000    continue
c      enddo
c      
c      
c      return
c      end
