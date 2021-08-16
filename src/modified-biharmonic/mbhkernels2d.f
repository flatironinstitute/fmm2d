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
     2     targ,nt,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      
      
            
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directdp_vec(nd,beta,source,ns,
     1     dipstr,dipvec,
     2     targ,nt,pot,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      
            
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdp_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,
     2     targ,nt,pot,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      
            
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directqp_vec(nd,beta,source,ns,
     1     quadstr,quadvec,
     2     targ,nt,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(nd,*),quadvec(nd,3,*)
            
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqp_vec(nd,beta,source,ns,
     1     charge,quadstr,quadvec,
     2     targ,nt,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      
      real *8 quadstr(nd,*),quadvec(nd,3,*)
            
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqp_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quadstr,quadvec,
     2     targ,nt,pot,thresh)
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
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      real *8 quadstr(nd,*),quadvec(nd,3,*)
            
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqp_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quadstr,quadvec,
     2     targ,nt,pot,thresh)
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
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      real *8 quadstr(nd,*),quadvec(nd,3,*)
            
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directop_vec(nd,beta,source,ns,
     1     octstr,octvec,
     2     targ,nt,pot,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcop_vec(nd,beta,source,ns,
     1     charge,octstr,octvec,
     2     targ,nt,pot,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      
      
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directdop_vec(nd,beta,source,ns,
     1     dipstr,dipvec,octstr,octvec,
     2     targ,nt,pot,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdop_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,octstr,octvec,
     2     targ,nt,pot,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directqop_vec(nd,beta,source,ns,
     1     quadstr,quadvec,octstr,octvec,
     2     targ,nt,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(nd,*),quadvec(nd,3,*)
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqop_vec(nd,beta,source,ns,
     1     charge,quadstr,quadvec,octstr,octvec,
     2     targ,nt,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      
      real *8 quadstr(nd,*),quadvec(nd,3,*)
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqop_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quadstr,quadvec,octstr,octvec,
     2     targ,nt,pot,thresh)
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
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      real *8 quadstr(nd,*),quadvec(nd,3,*)
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqop_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quadstr,quadvec,octstr,octvec,
     2     targ,nt,pot,thresh)
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
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      real *8 quadstr(nd,*),quadvec(nd,3,*)
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcg_vec(nd,beta,source,ns,
     1     charge,
     2     targ,nt,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      
      
            
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directdg_vec(nd,beta,source,ns,
     1     dipstr,dipvec,
     2     targ,nt,pot,grad,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      
            
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdg_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,
     2     targ,nt,pot,grad,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      
            
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directqg_vec(nd,beta,source,ns,
     1     quadstr,quadvec,
     2     targ,nt,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(nd,*),quadvec(nd,3,*)
            
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqg_vec(nd,beta,source,ns,
     1     charge,quadstr,quadvec,
     2     targ,nt,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      
      real *8 quadstr(nd,*),quadvec(nd,3,*)
            
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqg_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quadstr,quadvec,
     2     targ,nt,pot,grad,thresh)
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
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      real *8 quadstr(nd,*),quadvec(nd,3,*)
            
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqg_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quadstr,quadvec,
     2     targ,nt,pot,grad,thresh)
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
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      real *8 quadstr(nd,*),quadvec(nd,3,*)
            
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directog_vec(nd,beta,source,ns,
     1     octstr,octvec,
     2     targ,nt,pot,grad,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcog_vec(nd,beta,source,ns,
     1     charge,octstr,octvec,
     2     targ,nt,pot,grad,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      
      
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directdog_vec(nd,beta,source,ns,
     1     dipstr,dipvec,octstr,octvec,
     2     targ,nt,pot,grad,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdog_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,octstr,octvec,
     2     targ,nt,pot,grad,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directqog_vec(nd,beta,source,ns,
     1     quadstr,quadvec,octstr,octvec,
     2     targ,nt,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(nd,*),quadvec(nd,3,*)
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqog_vec(nd,beta,source,ns,
     1     charge,quadstr,quadvec,octstr,octvec,
     2     targ,nt,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      
      real *8 quadstr(nd,*),quadvec(nd,3,*)
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqog_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quadstr,quadvec,octstr,octvec,
     2     targ,nt,pot,grad,thresh)
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
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      real *8 quadstr(nd,*),quadvec(nd,3,*)
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqog_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quadstr,quadvec,octstr,octvec,
     2     targ,nt,pot,grad,thresh)
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
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      real *8 quadstr(nd,*),quadvec(nd,3,*)
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directch_vec(nd,beta,source,ns,
     1     charge,
     2     targ,nt,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      
      
            
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
               hess(ii,1,j) = hess(ii,1,j) + hessloc(1)*charge(ii,i)
               hess(ii,2,j) = hess(ii,2,j) + hessloc(2)*charge(ii,i)
               hess(ii,3,j) = hess(ii,3,j) + hessloc(3)*charge(ii,i)
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directdh_vec(nd,beta,source,ns,
     1     dipstr,dipvec,
     2     targ,nt,pot,grad,hess,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      
            
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))
               hess(ii,1,j) = hess(ii,1,j) -
     1              dipstr(ii,i)*(der3(1)*dipvec(ii,1,i)
     1              + der3(2)*dipvec(ii,2,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              dipstr(ii,i)*(der3(2)*dipvec(ii,1,i)
     1              + der3(3)*dipvec(ii,2,i))
               hess(ii,3,j) = hess(ii,3,j) -
     1              dipstr(ii,i)*(der3(3)*dipvec(ii,1,i)
     1              + der3(4)*dipvec(ii,2,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdh_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,
     2     targ,nt,pot,grad,hess,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      
            
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
               hess(ii,1,j) = hess(ii,1,j) + hessloc(1)*charge(ii,i)
               hess(ii,2,j) = hess(ii,2,j) + hessloc(2)*charge(ii,i)
               hess(ii,3,j) = hess(ii,3,j) + hessloc(3)*charge(ii,i)
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))
               hess(ii,1,j) = hess(ii,1,j) -
     1              dipstr(ii,i)*(der3(1)*dipvec(ii,1,i)
     1              + der3(2)*dipvec(ii,2,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              dipstr(ii,i)*(der3(2)*dipvec(ii,1,i)
     1              + der3(3)*dipvec(ii,2,i))
               hess(ii,3,j) = hess(ii,3,j) -
     1              dipstr(ii,i)*(der3(3)*dipvec(ii,1,i)
     1              + der3(4)*dipvec(ii,2,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directqh_vec(nd,beta,source,ns,
     1     quadstr,quadvec,
     2     targ,nt,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(nd,*),quadvec(nd,3,*)
            
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
               hess(ii,1,j) = hess(ii,1,j) +
     1              quadstr(ii,i)*(der4(1)*quadvec(ii,1,i)
     1              + der4(2)*quadvec(ii,2,i) + der4(3)*quadvec(ii,3,i))
               hess(ii,2,j) = hess(ii,2,j) +
     1              quadstr(ii,i)*(der4(2)*quadvec(ii,1,i)
     1              + der4(3)*quadvec(ii,2,i) + der4(4)*quadvec(ii,3,i))
               hess(ii,3,j) = hess(ii,3,j) +
     1              quadstr(ii,i)*(der4(3)*quadvec(ii,1,i)
     1              + der4(4)*quadvec(ii,2,i) + der4(5)*quadvec(ii,3,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqh_vec(nd,beta,source,ns,
     1     charge,quadstr,quadvec,
     2     targ,nt,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      
      real *8 quadstr(nd,*),quadvec(nd,3,*)
            
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
               hess(ii,1,j) = hess(ii,1,j) + hessloc(1)*charge(ii,i)
               hess(ii,2,j) = hess(ii,2,j) + hessloc(2)*charge(ii,i)
               hess(ii,3,j) = hess(ii,3,j) + hessloc(3)*charge(ii,i)

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
               hess(ii,1,j) = hess(ii,1,j) +
     1              quadstr(ii,i)*(der4(1)*quadvec(ii,1,i)
     1              + der4(2)*quadvec(ii,2,i) + der4(3)*quadvec(ii,3,i))
               hess(ii,2,j) = hess(ii,2,j) +
     1              quadstr(ii,i)*(der4(2)*quadvec(ii,1,i)
     1              + der4(3)*quadvec(ii,2,i) + der4(4)*quadvec(ii,3,i))
               hess(ii,3,j) = hess(ii,3,j) +
     1              quadstr(ii,i)*(der4(3)*quadvec(ii,1,i)
     1              + der4(4)*quadvec(ii,2,i) + der4(5)*quadvec(ii,3,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqh_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quadstr,quadvec,
     2     targ,nt,pot,grad,hess,thresh)
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
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      real *8 quadstr(nd,*),quadvec(nd,3,*)
            
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))
               hess(ii,1,j) = hess(ii,1,j) -
     1              dipstr(ii,i)*(der3(1)*dipvec(ii,1,i)
     1              + der3(2)*dipvec(ii,2,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              dipstr(ii,i)*(der3(2)*dipvec(ii,1,i)
     1              + der3(3)*dipvec(ii,2,i))
               hess(ii,3,j) = hess(ii,3,j) -
     1              dipstr(ii,i)*(der3(3)*dipvec(ii,1,i)
     1              + der3(4)*dipvec(ii,2,i))

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
               hess(ii,1,j) = hess(ii,1,j) +
     1              quadstr(ii,i)*(der4(1)*quadvec(ii,1,i)
     1              + der4(2)*quadvec(ii,2,i) + der4(3)*quadvec(ii,3,i))
               hess(ii,2,j) = hess(ii,2,j) +
     1              quadstr(ii,i)*(der4(2)*quadvec(ii,1,i)
     1              + der4(3)*quadvec(ii,2,i) + der4(4)*quadvec(ii,3,i))
               hess(ii,3,j) = hess(ii,3,j) +
     1              quadstr(ii,i)*(der4(3)*quadvec(ii,1,i)
     1              + der4(4)*quadvec(ii,2,i) + der4(5)*quadvec(ii,3,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqh_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quadstr,quadvec,
     2     targ,nt,pot,grad,hess,thresh)
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
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      real *8 quadstr(nd,*),quadvec(nd,3,*)
            
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
               hess(ii,1,j) = hess(ii,1,j) + hessloc(1)*charge(ii,i)
               hess(ii,2,j) = hess(ii,2,j) + hessloc(2)*charge(ii,i)
               hess(ii,3,j) = hess(ii,3,j) + hessloc(3)*charge(ii,i)
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))
               hess(ii,1,j) = hess(ii,1,j) -
     1              dipstr(ii,i)*(der3(1)*dipvec(ii,1,i)
     1              + der3(2)*dipvec(ii,2,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              dipstr(ii,i)*(der3(2)*dipvec(ii,1,i)
     1              + der3(3)*dipvec(ii,2,i))
               hess(ii,3,j) = hess(ii,3,j) -
     1              dipstr(ii,i)*(der3(3)*dipvec(ii,1,i)
     1              + der3(4)*dipvec(ii,2,i))

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
               hess(ii,1,j) = hess(ii,1,j) +
     1              quadstr(ii,i)*(der4(1)*quadvec(ii,1,i)
     1              + der4(2)*quadvec(ii,2,i) + der4(3)*quadvec(ii,3,i))
               hess(ii,2,j) = hess(ii,2,j) +
     1              quadstr(ii,i)*(der4(2)*quadvec(ii,1,i)
     1              + der4(3)*quadvec(ii,2,i) + der4(4)*quadvec(ii,3,i))
               hess(ii,3,j) = hess(ii,3,j) +
     1              quadstr(ii,i)*(der4(3)*quadvec(ii,1,i)
     1              + der4(4)*quadvec(ii,2,i) + der4(5)*quadvec(ii,3,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directoh_vec(nd,beta,source,ns,
     1     octstr,octvec,
     2     targ,nt,pot,grad,hess,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
               hess(ii,1,j) = hess(ii,1,j) - 
     1              octstr(ii,i)*(der5(1)*octvec(ii,1,i)
     1              + der5(2)*octvec(ii,2,i) + der5(3)*octvec(ii,3,i)
     1              + der5(4)*octvec(ii,4,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              octstr(ii,i)*(der5(2)*octvec(ii,1,i)
     1              + der5(3)*octvec(ii,2,i) + der5(4)*octvec(ii,3,i)
     1              + der5(5)*octvec(ii,4,i))
               hess(ii,3,j) = hess(ii,3,j) - 
     1              octstr(ii,i)*(der5(3)*octvec(ii,1,i)
     1              + der5(4)*octvec(ii,2,i) + der5(5)*octvec(ii,3,i)
     1              + der5(6)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcoh_vec(nd,beta,source,ns,
     1     charge,octstr,octvec,
     2     targ,nt,pot,grad,hess,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      
      
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
               hess(ii,1,j) = hess(ii,1,j) + hessloc(1)*charge(ii,i)
               hess(ii,2,j) = hess(ii,2,j) + hessloc(2)*charge(ii,i)
               hess(ii,3,j) = hess(ii,3,j) + hessloc(3)*charge(ii,i)
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
               hess(ii,1,j) = hess(ii,1,j) - 
     1              octstr(ii,i)*(der5(1)*octvec(ii,1,i)
     1              + der5(2)*octvec(ii,2,i) + der5(3)*octvec(ii,3,i)
     1              + der5(4)*octvec(ii,4,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              octstr(ii,i)*(der5(2)*octvec(ii,1,i)
     1              + der5(3)*octvec(ii,2,i) + der5(4)*octvec(ii,3,i)
     1              + der5(5)*octvec(ii,4,i))
               hess(ii,3,j) = hess(ii,3,j) - 
     1              octstr(ii,i)*(der5(3)*octvec(ii,1,i)
     1              + der5(4)*octvec(ii,2,i) + der5(5)*octvec(ii,3,i)
     1              + der5(6)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directdoh_vec(nd,beta,source,ns,
     1     dipstr,dipvec,octstr,octvec,
     2     targ,nt,pot,grad,hess,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))
               hess(ii,1,j) = hess(ii,1,j) -
     1              dipstr(ii,i)*(der3(1)*dipvec(ii,1,i)
     1              + der3(2)*dipvec(ii,2,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              dipstr(ii,i)*(der3(2)*dipvec(ii,1,i)
     1              + der3(3)*dipvec(ii,2,i))
               hess(ii,3,j) = hess(ii,3,j) -
     1              dipstr(ii,i)*(der3(3)*dipvec(ii,1,i)
     1              + der3(4)*dipvec(ii,2,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
               hess(ii,1,j) = hess(ii,1,j) - 
     1              octstr(ii,i)*(der5(1)*octvec(ii,1,i)
     1              + der5(2)*octvec(ii,2,i) + der5(3)*octvec(ii,3,i)
     1              + der5(4)*octvec(ii,4,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              octstr(ii,i)*(der5(2)*octvec(ii,1,i)
     1              + der5(3)*octvec(ii,2,i) + der5(4)*octvec(ii,3,i)
     1              + der5(5)*octvec(ii,4,i))
               hess(ii,3,j) = hess(ii,3,j) - 
     1              octstr(ii,i)*(der5(3)*octvec(ii,1,i)
     1              + der5(4)*octvec(ii,2,i) + der5(5)*octvec(ii,3,i)
     1              + der5(6)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdoh_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,octstr,octvec,
     2     targ,nt,pot,grad,hess,thresh)
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
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
               hess(ii,1,j) = hess(ii,1,j) + hessloc(1)*charge(ii,i)
               hess(ii,2,j) = hess(ii,2,j) + hessloc(2)*charge(ii,i)
               hess(ii,3,j) = hess(ii,3,j) + hessloc(3)*charge(ii,i)
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))
               hess(ii,1,j) = hess(ii,1,j) -
     1              dipstr(ii,i)*(der3(1)*dipvec(ii,1,i)
     1              + der3(2)*dipvec(ii,2,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              dipstr(ii,i)*(der3(2)*dipvec(ii,1,i)
     1              + der3(3)*dipvec(ii,2,i))
               hess(ii,3,j) = hess(ii,3,j) -
     1              dipstr(ii,i)*(der3(3)*dipvec(ii,1,i)
     1              + der3(4)*dipvec(ii,2,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
               hess(ii,1,j) = hess(ii,1,j) - 
     1              octstr(ii,i)*(der5(1)*octvec(ii,1,i)
     1              + der5(2)*octvec(ii,2,i) + der5(3)*octvec(ii,3,i)
     1              + der5(4)*octvec(ii,4,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              octstr(ii,i)*(der5(2)*octvec(ii,1,i)
     1              + der5(3)*octvec(ii,2,i) + der5(4)*octvec(ii,3,i)
     1              + der5(5)*octvec(ii,4,i))
               hess(ii,3,j) = hess(ii,3,j) - 
     1              octstr(ii,i)*(der5(3)*octvec(ii,1,i)
     1              + der5(4)*octvec(ii,2,i) + der5(5)*octvec(ii,3,i)
     1              + der5(6)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directqoh_vec(nd,beta,source,ns,
     1     quadstr,quadvec,octstr,octvec,
     2     targ,nt,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      
      real *8 quadstr(nd,*),quadvec(nd,3,*)
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
               hess(ii,1,j) = hess(ii,1,j) +
     1              quadstr(ii,i)*(der4(1)*quadvec(ii,1,i)
     1              + der4(2)*quadvec(ii,2,i) + der4(3)*quadvec(ii,3,i))
               hess(ii,2,j) = hess(ii,2,j) +
     1              quadstr(ii,i)*(der4(2)*quadvec(ii,1,i)
     1              + der4(3)*quadvec(ii,2,i) + der4(4)*quadvec(ii,3,i))
               hess(ii,3,j) = hess(ii,3,j) +
     1              quadstr(ii,i)*(der4(3)*quadvec(ii,1,i)
     1              + der4(4)*quadvec(ii,2,i) + der4(5)*quadvec(ii,3,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
               hess(ii,1,j) = hess(ii,1,j) - 
     1              octstr(ii,i)*(der5(1)*octvec(ii,1,i)
     1              + der5(2)*octvec(ii,2,i) + der5(3)*octvec(ii,3,i)
     1              + der5(4)*octvec(ii,4,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              octstr(ii,i)*(der5(2)*octvec(ii,1,i)
     1              + der5(3)*octvec(ii,2,i) + der5(4)*octvec(ii,3,i)
     1              + der5(5)*octvec(ii,4,i))
               hess(ii,3,j) = hess(ii,3,j) - 
     1              octstr(ii,i)*(der5(3)*octvec(ii,1,i)
     1              + der5(4)*octvec(ii,2,i) + der5(5)*octvec(ii,3,i)
     1              + der5(6)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcqoh_vec(nd,beta,source,ns,
     1     charge,quadstr,quadvec,octstr,octvec,
     2     targ,nt,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     charge - real *8 (nd,ns) charge strengths
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      
      real *8 quadstr(nd,*),quadvec(nd,3,*)
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
               hess(ii,1,j) = hess(ii,1,j) + hessloc(1)*charge(ii,i)
               hess(ii,2,j) = hess(ii,2,j) + hessloc(2)*charge(ii,i)
               hess(ii,3,j) = hess(ii,3,j) + hessloc(3)*charge(ii,i)

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
               hess(ii,1,j) = hess(ii,1,j) +
     1              quadstr(ii,i)*(der4(1)*quadvec(ii,1,i)
     1              + der4(2)*quadvec(ii,2,i) + der4(3)*quadvec(ii,3,i))
               hess(ii,2,j) = hess(ii,2,j) +
     1              quadstr(ii,i)*(der4(2)*quadvec(ii,1,i)
     1              + der4(3)*quadvec(ii,2,i) + der4(4)*quadvec(ii,3,i))
               hess(ii,3,j) = hess(ii,3,j) +
     1              quadstr(ii,i)*(der4(3)*quadvec(ii,1,i)
     1              + der4(4)*quadvec(ii,2,i) + der4(5)*quadvec(ii,3,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
               hess(ii,1,j) = hess(ii,1,j) - 
     1              octstr(ii,i)*(der5(1)*octvec(ii,1,i)
     1              + der5(2)*octvec(ii,2,i) + der5(3)*octvec(ii,3,i)
     1              + der5(4)*octvec(ii,4,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              octstr(ii,i)*(der5(2)*octvec(ii,1,i)
     1              + der5(3)*octvec(ii,2,i) + der5(4)*octvec(ii,3,i)
     1              + der5(5)*octvec(ii,4,i))
               hess(ii,3,j) = hess(ii,3,j) - 
     1              octstr(ii,i)*(der5(3)*octvec(ii,1,i)
     1              + der5(4)*octvec(ii,2,i) + der5(5)*octvec(ii,3,i)
     1              + der5(6)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directdqoh_vec(nd,beta,source,ns,
     1     dipstr,dipvec,quadstr,quadvec,octstr,octvec,
     2     targ,nt,pot,grad,hess,thresh)
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
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      real *8 quadstr(nd,*),quadvec(nd,3,*)
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))
               hess(ii,1,j) = hess(ii,1,j) -
     1              dipstr(ii,i)*(der3(1)*dipvec(ii,1,i)
     1              + der3(2)*dipvec(ii,2,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              dipstr(ii,i)*(der3(2)*dipvec(ii,1,i)
     1              + der3(3)*dipvec(ii,2,i))
               hess(ii,3,j) = hess(ii,3,j) -
     1              dipstr(ii,i)*(der3(3)*dipvec(ii,1,i)
     1              + der3(4)*dipvec(ii,2,i))

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
               hess(ii,1,j) = hess(ii,1,j) +
     1              quadstr(ii,i)*(der4(1)*quadvec(ii,1,i)
     1              + der4(2)*quadvec(ii,2,i) + der4(3)*quadvec(ii,3,i))
               hess(ii,2,j) = hess(ii,2,j) +
     1              quadstr(ii,i)*(der4(2)*quadvec(ii,1,i)
     1              + der4(3)*quadvec(ii,2,i) + der4(4)*quadvec(ii,3,i))
               hess(ii,3,j) = hess(ii,3,j) +
     1              quadstr(ii,i)*(der4(3)*quadvec(ii,1,i)
     1              + der4(4)*quadvec(ii,2,i) + der4(5)*quadvec(ii,3,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
               hess(ii,1,j) = hess(ii,1,j) - 
     1              octstr(ii,i)*(der5(1)*octvec(ii,1,i)
     1              + der5(2)*octvec(ii,2,i) + der5(3)*octvec(ii,3,i)
     1              + der5(4)*octvec(ii,4,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              octstr(ii,i)*(der5(2)*octvec(ii,1,i)
     1              + der5(3)*octvec(ii,2,i) + der5(4)*octvec(ii,3,i)
     1              + der5(5)*octvec(ii,4,i))
               hess(ii,3,j) = hess(ii,3,j) - 
     1              octstr(ii,i)*(der5(3)*octvec(ii,1,i)
     1              + der5(4)*octvec(ii,2,i) + der5(5)*octvec(ii,3,i)
     1              + der5(6)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
      enddo
      
      
      return
      end


      subroutine mbh2d_directcdqoh_vec(nd,beta,source,ns,
     1     charge,dipstr,dipvec,quadstr,quadvec,octstr,octvec,
     2     targ,nt,pot,grad,hess,thresh)
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
c     quadstr - real *8 (nd,ns) quadrupole strengths
c     quadvec - real *8 (nd,3,ns) quadrupole orientations
c     octstr - real *8 (nd,ns) octopole strengths
c     octvec - real *8 (nd,4,ns) octopole orientations
c     targ - real *8 (2) target locations
c     nt - integer, number of targets
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      real *8 charge(nd,*)
      real *8 dipstr(nd,*),dipvec(nd,2,*)
      real *8 quadstr(nd,*),quadvec(nd,3,*)
      real *8 octstr(nd,*),octvec(nd,4,*)      
      real *8 targ(2,nt)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,nt
c     local variables
      real *8 potloc,gradloc(2),hessloc(3),der3(4),der4(5),der5(6)
      real *8 thresh2, xdiff, ydiff, rr
      integer ifpotloc, ifgradloc, ifhessloc, ifder3, ifder4, ifder5
      integer i, ii, j

      thresh2=thresh*thresh

      ifpotloc=1
      ifgradloc=1
      ifhessloc=1
      ifder3=1
      ifder4=1
      ifder5=1

      do i = 1,ns

         do j = 1,nt
         
            xdiff=targ(1,j)-source(1,i)
            ydiff=targ(2,j)-source(2,i)
            rr=xdiff*xdiff+ydiff*ydiff
            if(rr.lt.thresh2) goto 1000
         
            call modbhgreen_all(beta,targ(1,j),source(1,i),ifpotloc,
     1           potloc,
     1           ifgradloc,gradloc,ifhessloc,hessloc,ifder3,der3,ifder4,
     2           der4,ifder5,der5)

            do ii = 1,nd	 
c     charge contrib
               pot(ii,j) = pot(ii,j) + potloc*charge(ii,i)
               grad(ii,1,j) = grad(ii,1,j) + gradloc(1)*charge(ii,i)
               grad(ii,2,j) = grad(ii,2,j) + gradloc(2)*charge(ii,i)
               hess(ii,1,j) = hess(ii,1,j) + hessloc(1)*charge(ii,i)
               hess(ii,2,j) = hess(ii,2,j) + hessloc(2)*charge(ii,i)
               hess(ii,3,j) = hess(ii,3,j) + hessloc(3)*charge(ii,i)
c     dipole contrib
               pot(ii,j) = pot(ii,j) -
     1              dipstr(ii,i)*(gradloc(1)*dipvec(ii,1,i)
     1              + gradloc(2)*dipvec(ii,2,i))
               grad(ii,1,j) = grad(ii,1,j) -
     1              dipstr(ii,i)*(hessloc(1)*dipvec(ii,1,i)
     1              + hessloc(2)*dipvec(ii,2,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              dipstr(ii,i)*(hessloc(2)*dipvec(ii,1,i)
     1              + hessloc(3)*dipvec(ii,2,i))
               hess(ii,1,j) = hess(ii,1,j) -
     1              dipstr(ii,i)*(der3(1)*dipvec(ii,1,i)
     1              + der3(2)*dipvec(ii,2,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              dipstr(ii,i)*(der3(2)*dipvec(ii,1,i)
     1              + der3(3)*dipvec(ii,2,i))
               hess(ii,3,j) = hess(ii,3,j) -
     1              dipstr(ii,i)*(der3(3)*dipvec(ii,1,i)
     1              + der3(4)*dipvec(ii,2,i))

c     quadrupole contrib
               pot(ii,j) = pot(ii,j) + 
     1              quadstr(ii,i)*(hessloc(1)*quadvec(ii,1,i)
     1              + hessloc(2)*quadvec(ii,2,i) + 
     1              hessloc(3)*quadvec(ii,3,i))
               grad(ii,1,j) = grad(ii,1,j) + 
     1              quadstr(ii,i)*(der3(1)*quadvec(ii,1,i)
     1              + der3(2)*quadvec(ii,2,i) 
     1              + der3(3)*quadvec(ii,3,i))
               grad(ii,2,j) = grad(ii,2,j) + 
     1              quadstr(ii,i)*(der3(2)*quadvec(ii,1,i)
     1              + der3(3)*quadvec(ii,2,i) +
     1              der3(4)*quadvec(ii,3,i))
               hess(ii,1,j) = hess(ii,1,j) +
     1              quadstr(ii,i)*(der4(1)*quadvec(ii,1,i)
     1              + der4(2)*quadvec(ii,2,i) + der4(3)*quadvec(ii,3,i))
               hess(ii,2,j) = hess(ii,2,j) +
     1              quadstr(ii,i)*(der4(2)*quadvec(ii,1,i)
     1              + der4(3)*quadvec(ii,2,i) + der4(4)*quadvec(ii,3,i))
               hess(ii,3,j) = hess(ii,3,j) +
     1              quadstr(ii,i)*(der4(3)*quadvec(ii,1,i)
     1              + der4(4)*quadvec(ii,2,i) + der4(5)*quadvec(ii,3,i))
c     octopole contrib            
               pot(ii,j) = pot(ii,j) - 
     1              octstr(ii,i)*(der3(1)*octvec(ii,1,i)
     1              + der3(2)*octvec(ii,2,i) + der3(3)*octvec(ii,3,i)
     1              + der3(4)*octvec(ii,4,i))
               grad(ii,1,j) = grad(ii,1,j) - 
     1              octstr(ii,i)*(der4(1)*octvec(ii,1,i)
     1              + der4(2)*octvec(ii,2,i) + der4(3)*octvec(ii,3,i)
     1              + der4(4)*octvec(ii,4,i))
               grad(ii,2,j) = grad(ii,2,j) - 
     1              octstr(ii,i)*(der4(2)*octvec(ii,1,i)
     1              + der4(3)*octvec(ii,2,i) + der4(4)*octvec(ii,3,i)
     1              + der4(5)*octvec(ii,4,i))
               hess(ii,1,j) = hess(ii,1,j) - 
     1              octstr(ii,i)*(der5(1)*octvec(ii,1,i)
     1              + der5(2)*octvec(ii,2,i) + der5(3)*octvec(ii,3,i)
     1              + der5(4)*octvec(ii,4,i))
               hess(ii,2,j) = hess(ii,2,j) - 
     1              octstr(ii,i)*(der5(2)*octvec(ii,1,i)
     1              + der5(3)*octvec(ii,2,i) + der5(4)*octvec(ii,3,i)
     1              + der5(5)*octvec(ii,4,i))
               hess(ii,3,j) = hess(ii,3,j) - 
     1              octstr(ii,i)*(der5(3)*octvec(ii,1,i)
     1              + der5(4)*octvec(ii,2,i) + der5(5)*octvec(ii,3,i)
     1              + der5(6)*octvec(ii,4,i))
            enddo
 1000       continue
            
         enddo	 
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
c     1     charge,dipstr,dipvec,quadstr,quadvec,octstr,octvec,
c     2     targ,pot,grad,hess,thresh)
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
cc     quadstr - real *8 (nd,ns) quadrupole strengths
cc     quadvec - real *8 (nd,3,ns) quadrupole orientations      
cc     octstr - real *8 (nd,ns) octopole strengths
cc     octvec - real *8 (nd,4,ns) octopole orientations
cc     targ - real *8 (2) target location
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
c      real *8 dipvec(2,*), quadstr(*), quadvec(3,*)
c      real *8 octstr(*), octvec(4,*)
c      real *8 targ(2)
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
c         xdiff=targ(1)-source(1,i)
c         ydiff=targ(2)-source(2,i)
c         rr=xdiff*xdiff+ydiff*ydiff
c         if(rr.lt.thresh2) goto 1000
c         
c         call modbhgreen_all(beta,targ,source(1,i),ifpotloc,potloc,
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
c            pot(ii) = pot(ii) + quadstr(i)*(hessloc(1)*quadvec(1,i)
c     1           + hessloc(2)*quadvec(2,i) + hessloc(3)*quadvec(3,i))
c
c            grad(ii,1) = grad(ii,1) + quadstr(i)*(der3(1)*quadvec(1,i)
c     1           + der3(2)*quadvec(2,i) + der3(3)*quadvec(3,i))
c            grad(ii,2) = grad(ii,2) + quadstr(i)*(der3(2)*quadvec(1,i)
c     1           + der3(3)*quadvec(2,i) + der3(4)*quadvec(3,i))
c
c            hess(ii,1) = hess(ii,1) + quadstr(i)*(der4(1)*quadvec(1,i)
c     1           + der4(2)*quadvec(2,i) + der4(3)*quadvec(3,i))
c            hess(ii,2) = hess(ii,2) + quadstr(i)*(der4(2)*quadvec(1,i)
c     1           + der4(3)*quadvec(2,i) + der4(4)*quadvec(3,i))
c            hess(ii,3) = hess(ii,3) + quadstr(i)*(der4(3)*quadvec(1,i)
c     1           + der4(4)*quadvec(2,i) + der4(5)*quadvec(3,i))
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



      subroutine mbh2d_directmpsp_vec(nd,beta,source,ns,
     1     mbhmps,ymps,ntermsmps,
     2     targ,nt,pot,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     multipolar sources direct evaluation routine
c      
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     mbhmps       in: complex *16 (nd,0:ntermsmps,ns)
c                  difference kernel type expansion at each source 
c
c     ymps         in: complex *16 (nd,0:ntermsmps,ns)
c                  Yukawa (Bessel K) type expansion at each source 
c
c     ntermsmps     in: integer, order of the multipolar charges at
c                   each source
c      
c     targ - real *8 (2,nt) target location
c     nt - integer, number of targets 
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      complex *16 :: mbhmps(nd,0:ntermsmps,ns),ymps(nd,0:ntermsmps,ns)
      real *8 targ(2,*)
      real *8 pot(nd,nt)
      
            
      real *8 thresh
      integer ns,nd,ntermsmps,nt
c     local variables
      real *8 zdiff(2), r, theta
      real *8 ders(0:1), pih
      real *8, allocatable :: diffs(:), kvec(:)
      complex *16 z, eye, ztemp1, ztemp2
      complex *16, allocatable :: mptemp1(:),mptemp2(:)      
      real *8 dc, ds, dc2, ds2, dcx, dsx, dcy, dsy
      real *8 rscale
      integer ifder, ifders, l, j, i, nterms, jj
      data eye /(0.0d0,1.0d0)/

      rscale=1
      nterms = ntermsmps

      allocate(mptemp1(0:nterms+6),mptemp2(0:nterms+6),
     1     diffs(0:nterms+6),kvec(0:nterms+6))

      do jj = 1,ns


         do i = 1,nt
            zdiff(1) = targ(1,i)-source(1,jj)
            zdiff(2) = targ(2,i)-source(2,jj)
            call h2cart2polar(zdiff,r,theta)

            if (r .lt. thresh) cycle
            
c     get values of difference functions
            ifders = 0
            call diffslogbk_fast(r,beta,rscale,diffs,ifders,ders,kvec,
     1           nterms+2)       
            
            mptemp1(0)=diffs(0)
            mptemp2(0)=kvec(0)
            ztemp2=exp(eye*theta)
            ztemp1=ztemp2
            do j=1,nterms+2
               mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1              diffs(j)*dimag(ztemp1))
               mptemp2(j)=dcmplx(kvec(j)*dreal(ztemp1),
     1              kvec(j)*dimag(ztemp1))
               ztemp1 = ztemp1*ztemp2
            enddo
         
c     evaluate
            do j = 1,nd
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp1(0))*dreal(mbhmps(j,0,jj)) +
     1              dreal(mptemp2(0))*dreal(ymps(j,0,jj))
            enddo
            
            do l = 1,nterms
               do j=1,nd
                  pot(j,i) = pot(j,i) +
     1                 dreal(mptemp1(l))*dreal(mbhmps(j,l,jj)) +
     1                 dimag(mptemp1(l))*dimag(mbhmps(j,l,jj))
                  pot(j,i) = pot(j,i) +
     1                 dreal(mptemp2(l))*dreal(ymps(j,l,jj)) +
     1                 dimag(mptemp2(l))*dimag(ymps(j,l,jj))
               enddo
            enddo

         enddo
      enddo
      return 
      end


      subroutine mbh2d_directmpsg_vec(nd,beta,source,ns,
     1     mbhmps,ymps,ntermsmps,
     2     targ,nt,pot,grad,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     multipolar sources direct evaluation routine
c      
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     mbhmps       in: complex *16 (nd,0:ntermsmps,ns)
c                  difference kernel type expansion at each source 
c
c     ymps         in: complex *16 (nd,0:ntermsmps,ns)
c                  Yukawa (Bessel K) type expansion at each source 
c
c     ntermsmps     in: integer, order of the multipolar charges at
c                   each source
c      
c     targ - real *8 (2,nt) target location
c     nt - integer, number of targets 
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      complex *16 :: mbhmps(nd,0:ntermsmps,ns),ymps(nd,0:ntermsmps,ns)
      real *8 targ(2,*)
      
      real *8 pot(nd,nt),grad(nd,2,nt)
            
      real *8 thresh
      integer ns,nd,ntermsmps,nt
c     local variables
      real *8 zdiff(2), r, theta
      real *8 ders(0:1), pih
      real *8, allocatable :: diffs(:), kvec(:)
      complex *16 z, eye, ztemp1, ztemp2
      complex *16, allocatable :: mptemp1(:),mptemp2(:)
      complex *16, allocatable :: ympolex(:,:), ympoley(:,:)
      complex *16, allocatable :: mbhmpolex(:,:), mbhmpoley(:,:)      
      real *8 dc, ds, dc2, ds2, dcx, dsx, dcy, dsy
      real *8 rscale
      integer ifder, ifders, l, j, i, nterms, jj
      data eye /(0.0d0,1.0d0)/

      rscale=1
      nterms = ntermsmps
c     
      allocate(ympolex(nd,0:nterms+1),ympoley(nd,0:nterms+1),
     1     mbhmpolex(nd,0:nterms+1),mbhmpoley(nd,0:nterms+1))

      allocate(mptemp1(0:nterms+6),mptemp2(0:nterms+6),
     1     diffs(0:nterms+6),kvec(0:nterms+6))

      do jj = 1,ns
         
         do i=0,nterms+1
            do j = 1,nd
               ympolex(j,i)=0
               ympoley(j,i)=0
               mbhmpolex(j,i)=0
               mbhmpoley(j,i)=0
            enddo
         enddo

         do j=1,nd
            ympolex(j,1) = -beta/rscale*ymps(j,0,jj)
            ympoley(j,1) = -beta/rscale*ymps(j,0,jj)*eye
            
            mbhmpolex(j,1) = -beta/rscale*mbhmps(j,0,jj)
            mbhmpoley(j,1) = -beta/rscale*mbhmps(j,0,jj)*eye
         enddo
         
         do i=1,nterms
            do j = 1,nd
               ympolex(j,i-1) = ympolex(j,i-1) 
     1              -beta/2.0d0*ymps(j,i,jj)*rscale
               ympolex(j,i+1) = ympolex(j,i+1) 
     1              -beta/2.0d0*ymps(j,i,jj)/rscale
               ympoley(j,i-1) = ympoley(j,i-1) 
     1              +beta/2.0d0*ymps(j,i,jj)*rscale*eye
               ympoley(j,i+1) = ympoley(j,i+1)
     1              -beta/2.0d0*ymps(j,i,jj)/rscale*eye

               ympolex(j,i-1) = ympolex(j,i-1) 
     1              -beta/2.0d0*mbhmps(j,i,jj)*rscale
               
               mbhmpolex(j,i+1) = mbhmpolex(j,i+1) 
     1              -beta/2.0d0*mbhmps(j,i,jj)/rscale
               ympoley(j,i-1) = ympoley(j,i-1) 
     1              +beta/2.0d0*mbhmps(j,i,jj)*rscale*eye
               mbhmpoley(j,i+1) = mbhmpoley(j,i+1) 
     1              -beta/2.0d0*mbhmps(j,i,jj)/rscale*eye
            enddo
         enddo

         


         do i = 1,nt
            zdiff(1) = targ(1,i)-source(1,jj)
            zdiff(2) = targ(2,i)-source(2,jj)
            call h2cart2polar(zdiff,r,theta)

            if (r .lt. thresh) cycle
            
c     get values of difference functions
            ifders = 0
            call diffslogbk_fast(r,beta,rscale,diffs,ifders,ders,kvec,
     1           nterms+2)       
            
            mptemp1(0)=diffs(0)
            mptemp2(0)=kvec(0)
            ztemp2=exp(eye*theta)
            ztemp1=ztemp2
            do j=1,nterms+2
               mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1              diffs(j)*dimag(ztemp1))
               mptemp2(j)=dcmplx(kvec(j)*dreal(ztemp1),
     1              kvec(j)*dimag(ztemp1))
               ztemp1 = ztemp1*ztemp2
            enddo
         
c     evaluate
            do j = 1,nd
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp1(0))*dreal(mbhmps(j,0,jj)) +
     1              dreal(mptemp2(0))*dreal(ymps(j,0,jj))
            enddo
            
            do l = 1,nterms
               do j=1,nd
                  pot(j,i) = pot(j,i) +
     1                 dreal(mptemp1(l))*dreal(mbhmps(j,l,jj)) +
     1                 dimag(mptemp1(l))*dimag(mbhmps(j,l,jj))
                  pot(j,i) = pot(j,i) +
     1                 dreal(mptemp2(l))*dreal(ymps(j,l,jj)) +
     1                 dimag(mptemp2(l))*dimag(ymps(j,l,jj))
               enddo
            enddo
            do j = 1,nd
               grad(j,1,i) = grad(j,1,i) + 
     1              dreal(mptemp1(0))*dreal(mbhmpolex(j,0)) +
     1              dreal(mptemp2(0))*dreal(ympolex(j,0))
               grad(j,2,i) = grad(j,2,i) + 
     1              dreal(mptemp1(0))*dreal(mbhmpoley(j,0)) +
     1              dreal(mptemp2(0))*dreal(ympoley(j,0))
            enddo
            do l = 1,nterms+1
               do j = 1,nd
                  grad(j,1,i) = grad(j,1,i) +
     1                 dreal(mptemp1(l))*dreal(mbhmpolex(j,l)) +
     1                 dimag(mptemp1(l))*dimag(mbhmpolex(j,l))
                  grad(j,1,i) = grad(j,1,i) +
     1                 dreal(mptemp2(l))*dreal(ympolex(j,l)) +
     1                 dimag(mptemp2(l))*dimag(ympolex(j,l))
                  grad(j,2,i) = grad(j,2,i) +
     1                 dreal(mptemp1(l))*dreal(mbhmpoley(j,l)) +
     1                 dimag(mptemp1(l))*dimag(mbhmpoley(j,l))
                  grad(j,2,i) = grad(j,2,i) +
     1                 dreal(mptemp2(l))*dreal(ympoley(j,l)) +
     1                 dimag(mptemp2(l))*dimag(ympoley(j,l))
               enddo
            enddo

         enddo
      enddo
      return 
      end


      subroutine mbh2d_directmpsh_vec(nd,beta,source,ns,
     1     mbhmps,ymps,ntermsmps,
     2     targ,nt,pot,grad,hess,thresh)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     multipolar sources direct evaluation routine
c      
c     input:
c
c     nd - integer, number of vectors per source
c     beta - real *8, modified biharmonic parameter
c     ns - integer, number of sources
c     source - real *8 (2,ns) source locations
c     mbhmps       in: complex *16 (nd,0:ntermsmps,ns)
c                  difference kernel type expansion at each source 
c
c     ymps         in: complex *16 (nd,0:ntermsmps,ns)
c                  Yukawa (Bessel K) type expansion at each source 
c
c     ntermsmps     in: integer, order of the multipolar charges at
c                   each source
c      
c     targ - real *8 (2,nt) target location
c     nt - integer, number of targets 
c     thresh - real *8 threshold, don't compute contribution of
c                     charge if distance from charge to target is
c                     less than this threshold
c
c     output:
c      
c     pot(nd,nt)      : value of potential at targets
c     grad(nd,2,nt)   : value of gradient at targets
c     hess(nd,3,nt)   : value of Hessian at targets      
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c     global variables
      real *8 beta, source(2,*)
      complex *16 :: mbhmps(nd,0:ntermsmps,ns),ymps(nd,0:ntermsmps,ns)
      real *8 targ(2,*)
      
      
      real *8 pot(nd,nt),grad(nd,2,nt),hess(nd,3,nt)      
      real *8 thresh
      integer ns,nd,ntermsmps,nt
c     local variables
      real *8 zdiff(2), r, theta
      real *8 ders(0:1), pih
      real *8, allocatable :: diffs(:), kvec(:)
      complex *16 z, eye, ztemp1, ztemp2
      complex *16, allocatable :: mptemp1(:),mptemp2(:)
      complex *16, allocatable :: ympolex(:,:), ympoley(:,:)
      complex *16, allocatable :: mbhmpolex(:,:), mbhmpoley(:,:)
      complex *16, allocatable :: ympolexx(:,:), ympolexy(:,:),
     1     ympoleyy(:,:)
      complex *16, allocatable :: mbhmpolexx(:,:),mbhmpolexy(:,:),
     1     mbhmpoleyy(:,:)      
      real *8 dc, ds, dc2, ds2, dcx, dsx, dcy, dsy
      real *8 rscale
      integer ifder, ifders, l, j, i, nterms, jj
      data eye /(0.0d0,1.0d0)/

      rscale=1
      nterms = ntermsmps
c     
      allocate(ympolex(nd,0:nterms+1),ympoley(nd,0:nterms+1),
     1     mbhmpolex(nd,0:nterms+1),mbhmpoley(nd,0:nterms+1))
c
      allocate(ympolexx(nd,0:nterms+2),ympolexy(nd,0:nterms+2),
     1     ympoleyy(nd,0:nterms+2),
     1     mbhmpolexx(nd,0:nterms+2),mbhmpolexy(nd,0:nterms+2),
     1     mbhmpoleyy(nd,0:nterms+2))

      allocate(mptemp1(0:nterms+6),mptemp2(0:nterms+6),
     1     diffs(0:nterms+6),kvec(0:nterms+6))

      do jj = 1,ns
         
         do i=0,nterms+1
            do j = 1,nd
               ympolex(j,i)=0
               ympoley(j,i)=0
               mbhmpolex(j,i)=0
               mbhmpoley(j,i)=0
            enddo
         enddo

         do j=1,nd
            ympolex(j,1) = -beta/rscale*ymps(j,0,jj)
            ympoley(j,1) = -beta/rscale*ymps(j,0,jj)*eye
            
            mbhmpolex(j,1) = -beta/rscale*mbhmps(j,0,jj)
            mbhmpoley(j,1) = -beta/rscale*mbhmps(j,0,jj)*eye
         enddo
         
         do i=1,nterms
            do j = 1,nd
               ympolex(j,i-1) = ympolex(j,i-1) 
     1              -beta/2.0d0*ymps(j,i,jj)*rscale
               ympolex(j,i+1) = ympolex(j,i+1) 
     1              -beta/2.0d0*ymps(j,i,jj)/rscale
               ympoley(j,i-1) = ympoley(j,i-1) 
     1              +beta/2.0d0*ymps(j,i,jj)*rscale*eye
               ympoley(j,i+1) = ympoley(j,i+1)
     1              -beta/2.0d0*ymps(j,i,jj)/rscale*eye

               ympolex(j,i-1) = ympolex(j,i-1) 
     1              -beta/2.0d0*mbhmps(j,i,jj)*rscale
               
               mbhmpolex(j,i+1) = mbhmpolex(j,i+1) 
     1              -beta/2.0d0*mbhmps(j,i,jj)/rscale
               ympoley(j,i-1) = ympoley(j,i-1) 
     1              +beta/2.0d0*mbhmps(j,i,jj)*rscale*eye
               mbhmpoley(j,i+1) = mbhmpoley(j,i+1) 
     1              -beta/2.0d0*mbhmps(j,i,jj)/rscale*eye
            enddo
         enddo

         
c

         do i=0,nterms+2
            do j=1,nd
               ympolexx(j,i)=0
               ympolexy(j,i)=0
               ympoleyy(j,i)=0
               mbhmpolexx(j,i)=0
               mbhmpolexy(j,i)=0
               mbhmpoleyy(j,i)=0
            enddo
         enddo

         do j = 1,nd
            ympolexx(j,1) = -beta/1.0d0/rscale*dreal(ympolex(j,0))
            ympolexy(j,1) = -beta/1.0d0/rscale*dreal(ympolex(j,0))*eye
            ympoleyy(j,1) = -beta/1.0d0/rscale*dreal(ympoley(j,0))*eye

            mbhmpolexx(j,1) = -beta/1.0d0/rscale*dreal(mbhmpolex(j,0))
            mbhmpolexy(j,1) = 
     1           -beta/1.0d0/rscale*dreal(mbhmpolex(j,0))*eye
            mbhmpoleyy(j,1) = 
     1           -beta/1.0d0/rscale*dreal(mbhmpoley(j,0))*eye
         enddo
         
         do i=1,nterms+1
            do j=1,nd
               ympolexx(j,i-1) = ympolexx(j,i-1)
     1              -beta/2.0d0*ympolex(j,i)*rscale
               ympolexx(j,i+1) = ympolexx(j,i+1)
     1              -beta/2.0d0*ympolex(j,i)/rscale
               ympolexy(j,i-1) = ympolexy(j,i-1) 
     1              +beta/2.0d0*ympolex(j,i)*rscale*eye
               ympolexy(j,i+1) = ympolexy(j,i+1) 
     1              -beta/2.0d0*ympolex(j,i)/rscale*eye
               ympoleyy(j,i-1) = ympoleyy(j,i-1) 
     1              +beta/2.0d0*ympoley(j,i)*rscale*eye
               ympoleyy(j,i+1) = ympoleyy(j,i+1) 
     1              -beta/2.0d0*ympoley(j,i)/rscale*eye

               ympolexx(j,i-1) = ympolexx(j,i-1) 
     1              -beta/2.0d0*mbhmpolex(j,i)*rscale
               mbhmpolexx(j,i+1) = mbhmpolexx(j,i+1) 
     1              -beta/2.0d0*mbhmpolex(j,i)/rscale
               ympolexy(j,i-1) = ympolexy(j,i-1) 
     1              +beta/2.0d0*mbhmpolex(j,i)*rscale*eye
               mbhmpolexy(j,i+1) = mbhmpolexy(j,i+1) 
     1              -beta/2.0d0*mbhmpolex(j,i)/rscale*eye
               ympoleyy(j,i-1) = ympoleyy(j,i-1) 
     1              +beta/2.0d0*mbhmpoley(j,i)*rscale*eye
               mbhmpoleyy(j,i+1) = mbhmpoleyy(j,i+1) 
     1              -beta/2.0d0*mbhmpoley(j,i)/rscale*eye
            enddo
         enddo


         do i = 1,nt
            zdiff(1) = targ(1,i)-source(1,jj)
            zdiff(2) = targ(2,i)-source(2,jj)
            call h2cart2polar(zdiff,r,theta)

            if (r .lt. thresh) cycle
            
c     get values of difference functions
            ifders = 0
            call diffslogbk_fast(r,beta,rscale,diffs,ifders,ders,kvec,
     1           nterms+2)       
            
            mptemp1(0)=diffs(0)
            mptemp2(0)=kvec(0)
            ztemp2=exp(eye*theta)
            ztemp1=ztemp2
            do j=1,nterms+2
               mptemp1(j)=dcmplx(diffs(j)*dreal(ztemp1),
     1              diffs(j)*dimag(ztemp1))
               mptemp2(j)=dcmplx(kvec(j)*dreal(ztemp1),
     1              kvec(j)*dimag(ztemp1))
               ztemp1 = ztemp1*ztemp2
            enddo
         
c     evaluate
            do j = 1,nd
               pot(j,i) = pot(j,i) +
     1              dreal(mptemp1(0))*dreal(mbhmps(j,0,jj)) +
     1              dreal(mptemp2(0))*dreal(ymps(j,0,jj))
            enddo
            
            do l = 1,nterms
               do j=1,nd
                  pot(j,i) = pot(j,i) +
     1                 dreal(mptemp1(l))*dreal(mbhmps(j,l,jj)) +
     1                 dimag(mptemp1(l))*dimag(mbhmps(j,l,jj))
                  pot(j,i) = pot(j,i) +
     1                 dreal(mptemp2(l))*dreal(ymps(j,l,jj)) +
     1                 dimag(mptemp2(l))*dimag(ymps(j,l,jj))
               enddo
            enddo
            do j = 1,nd
               grad(j,1,i) = grad(j,1,i) + 
     1              dreal(mptemp1(0))*dreal(mbhmpolex(j,0)) +
     1              dreal(mptemp2(0))*dreal(ympolex(j,0))
               grad(j,2,i) = grad(j,2,i) + 
     1              dreal(mptemp1(0))*dreal(mbhmpoley(j,0)) +
     1              dreal(mptemp2(0))*dreal(ympoley(j,0))
            enddo
            do l = 1,nterms+1
               do j = 1,nd
                  grad(j,1,i) = grad(j,1,i) +
     1                 dreal(mptemp1(l))*dreal(mbhmpolex(j,l)) +
     1                 dimag(mptemp1(l))*dimag(mbhmpolex(j,l))
                  grad(j,1,i) = grad(j,1,i) +
     1                 dreal(mptemp2(l))*dreal(ympolex(j,l)) +
     1                 dimag(mptemp2(l))*dimag(ympolex(j,l))
                  grad(j,2,i) = grad(j,2,i) +
     1                 dreal(mptemp1(l))*dreal(mbhmpoley(j,l)) +
     1                 dimag(mptemp1(l))*dimag(mbhmpoley(j,l))
                  grad(j,2,i) = grad(j,2,i) +
     1                 dreal(mptemp2(l))*dreal(ympoley(j,l)) +
     1                 dimag(mptemp2(l))*dimag(ympoley(j,l))
               enddo
            enddo
            do j=1,nd
               hess(j,1,i) = hess(j,1,i) +
     1              dreal(mptemp1(0))*dreal(mbhmpolexx(j,0)) +
     1              dreal(mptemp2(0))*dreal(ympolexx(j,0))
               hess(j,2,i) = hess(j,2,i) +
     1              dreal(mptemp1(0))*dreal(mbhmpolexy(j,0)) +
     1              dreal(mptemp2(0))*dreal(ympolexy(j,0))
               hess(j,3,i) = hess(j,3,i) +
     1              dreal(mptemp1(0))*dreal(mbhmpoleyy(j,0)) +
     1              dreal(mptemp2(0))*dreal(ympoleyy(j,0))
            enddo
            do l = 1,nterms+2
               do j = 1,nd
                  hess(j,1,i) = hess(j,1,i) +
     1                 dreal(mptemp1(l))*dreal(mbhmpolexx(j,l)) +
     1                 dimag(mptemp1(l))*dimag(mbhmpolexx(j,l))
                  hess(j,1,i) = hess(j,1,i) + 
     1                 dreal(mptemp2(l))*dreal(ympolexx(j,l)) +
     1                 dimag(mptemp2(l))*dimag(ympolexx(j,l))
                  hess(j,2,i) = hess(j,2,i) + 
     1                 dreal(mptemp1(l))*dreal(mbhmpolexy(j,l)) +
     1                 dimag(mptemp1(l))*dimag(mbhmpolexy(j,l))
                  hess(j,2,i) = hess(j,2,i) + 
     1                 dreal(mptemp2(l))*dreal(ympolexy(j,l)) +
     1                 dimag(mptemp2(l))*dimag(ympolexy(j,l))
                  hess(j,3,i) = hess(j,3,i) + 
     1                 dreal(mptemp1(l))*dreal(mbhmpoleyy(j,l)) +
     1                 dimag(mptemp1(l))*dimag(mbhmpoleyy(j,l))
                  hess(j,3,i) = hess(j,3,i) + 
     1                 dreal(mptemp2(l))*dreal(ympoleyy(j,l)) +
     1                 dimag(mptemp2(l))*dimag(ympoleyy(j,l))
               enddo
            enddo

         enddo
      enddo
      return 
      end



      
