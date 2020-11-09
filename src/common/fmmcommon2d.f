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
      double precision arr(ndim,1),arrsort(ndim,1)
      integer iarr(1)

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
      
