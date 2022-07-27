c----------------------------------------------------------------     

      subroutine ireorderf(ndim,n,arr,arrsort,iarr)
c
cc       this subroutine sorts the array arr and stores
c        it in arrsort using the sorting order defined by
c        iarr
c
c        arrsort(j,i) = arr(j,iarr(i)), j =1,2,\ldots ndim
c                                       i=1,2,\ldots n
c

      implicit none
      integer ndim,idim,i,n
      integer arr(ndim,n),arrsort(ndim,n)
      integer iarr(n)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)      
      do i=1,n
         do idim=1,ndim
            arrsort(idim,i) = arr(idim,iarr(i))
         enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
c----------------------------------------------------------------     

      subroutine dreorderf(ndim,n,arr,arrsort,iarr)
c
cc       this subroutine sorts the array arr and stores
c        it in arrsort using the sorting order defined by
c        iarr
c
c        arrsort(j,i) = arr(j,iarr(i)), j =1,2,\ldots ndim
c                                       i=1,2,\ldots n
c

      implicit none
      integer ndim,idim,i,n
      double precision arr(ndim,n),arrsort(ndim,n)
      integer iarr(n)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)      
      do i=1,n
         do idim=1,ndim
            arrsort(idim,i) = arr(idim,iarr(i))
         enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
c--------------------------------------------------
      subroutine dreorderi(ndim,n,arr,arrsort,iarr)
c
cc       this subroutine sorts the array arr and stores
c        it in arrsort using the inverse of the
c        sorting order defined by
c        iarr.
c
c        Note that this subroutine is the inverse of 
c        dreorderf
c
c        arrsort(j,iarr(i)) = arr(j,i), j =1,2,\ldots ndim
c                                       i=1,2,\ldots n
c

      implicit none
      integer i,idim,ndim,n
      double precision arr(ndim,1),arrsort(ndim,1)
      integer iarr(1)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,idim)      
      do i=1,n
         do idim=1,ndim
            arrsort(idim,iarr(i)) = arr(idim,i)
         enddo
      enddo
C$OMP END PARALLEL DO      

      return
      end
c----------------------------------------------------------     
      
c
c-------------------------------------------------------
      subroutine geterrstr(ifcharge,ifdipole,ifpgh,ifpghtarg,str1,len1)
      implicit real *8 (a-h,o-z)
      character(len=*) str1
      character(len=13) str2
      character(len=14) str3
      character(len=19) str4
      character(len=30) str5

      str2 = "Failed src to"
      len1 = 13
      if(ifpgh.gt.0.and.ifpghtarg.eq.0) then
        str3 = " src,"
        len1 = len1+5  
      endif
      if(ifpgh.eq.0.and.ifpghtarg.gt.0) then
        str3 = " targ,"
        len1 = len1+6
      endif
      if(ifpgh.gt.0.and.ifpghtarg.gt.0) then
        str3 = " src and targ,"
        len1 = len1+14
      endif

      if(ifcharge.eq.1.and.ifdipole.eq.0) then
        str4=" charge,"
        len1 = len1+8
      endif
      
      if(ifcharge.eq.0.and.ifdipole.eq.1) then
        str4=" dipole,"
        len1 = len1+8
      endif
      
      if(ifcharge.eq.1.and.ifdipole.eq.1) then
        str4=" charge and dipole,"
        len1 = len1+19
      endif

      if(ifpgh.eq.1.or.ifpghtarg.eq.1) then
        str5=" pot test"
        len1 = len1 + 9
      endif
      
      if(ifpgh.eq.2.or.ifpghtarg.eq.2) then
        str5=" pot and grad test"
        len1 = len1 + 18
      endif

      if(ifpgh.eq.3.or.ifpghtarg.eq.3) then
        str5=" pot, grad, and hess test"
        len1 = len1+25
      endif

      str1 = str2//trim(str3)//trim(str4)//trim(str5)

      return
      end

c
c
c
        subroutine init_carray(carray,ldc)
        implicit real *8 (a-h,o-z)
        real *8 carray(0:ldc,0:ldc)

        do l = 0,ldc
        carray(l,0) = 1.0d0
        enddo
        do m=1,ldc
        carray(m,m) = 1.0d0
        do l=m+1,ldc
            carray(l,m)=carray(l-1,m)+carray(l-1,m-1)
        enddo
        enddo
c
        return
        end
        
