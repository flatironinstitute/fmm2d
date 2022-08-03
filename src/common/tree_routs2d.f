c
c
c    common routines for generating and processing
c     a level restricted quad tree in 2D
c   
c
c
c
      subroutine tree_refine_boxes(irefinebox,nboxes,
     1  ifirstbox,nbloc,centers,bs,nbctr,nlctr,
     2  ilevel,iparent,nchild,ichild)
      implicit none
      integer nboxes,nbloc,nbctr,nlctr
      real *8 centers(2,nboxes),bs
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(4,nboxes),nchild(nboxes)
      integer irefinebox(nbloc)
      integer ifirstbox
      integer, allocatable :: isum(:)
      integer ii

      integer i,ibox,nel0,j,l,jbox,nel1,nbl

      allocate(isum(nbloc))
      if(nbloc.gt.0) call cumsum(nbloc,irefinebox,isum)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox,ii,l)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(irefinebox(i).eq.1) then
          nbl = nbctr + (isum(i)-1)*4
          
          nchild(ibox) = 4
          do j=1,4
            ii = 2
            if(j.le.2) ii = 1
            jbox = nbl+j
            centers(1,jbox) = centers(1,ibox)+(-1)**j*bs/2
            centers(2,jbox) = centers(2,ibox)+(-1)**ii*bs/2 
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,4
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr 
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*4


      return
      end
c
c
c
c
c
c
       subroutine tree_copy(nb,centers,ilevel,iparent,nchild,ichild,
     1              centers2,ilevel2,iparent2,nchild2,ichild2)

       implicit none
       integer nd,nb,npb
       real *8 centers(2,nb),centers2(2,nb)
       integer ilevel(nb),ilevel2(nb)
       integer iparent(nb),iparent2(nb)
       integer nchild(nb),nchild2(nb)
       integer ichild(4,nb),ichild2(4,nb)

       integer i,j,nel


C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
       do i=1,nb
         centers2(1,i) = centers(1,i)
         centers2(2,i) = centers(2,i)
         ilevel2(i) = ilevel(i)
         iparent2(i) = iparent(i)
         nchild2(i) = nchild(i)
         ichild2(1,i) = ichild(1,i)
         ichild2(2,i) = ichild(2,i)
         ichild2(3,i) = ichild(3,i)
         ichild2(4,i) = ichild(4,i)
       enddo
C$OMP END PARALLEL DO       
       

       return
       end
c
c
c
c
c

      subroutine computecoll(nlevels,nboxes,laddr,boxsize,
     1                       centers,iparent,nchild,ichild,iper,
     2                       nnbors,nbors)

c     This subroutine computes the colleagues for an adaptive
c     pruned tree. box j is a colleague of box i, if they share a
c     vertex or an edge and the two boxes are at the same
c     level in the tree
c
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     laddr       in: integer(2,0:nlevels)
c                 indexing array providing access to boxes at
c                 each level. 
c                 the first box on level i is laddr(1,i)
c                 the last box on level i is laddr(2,i)
c
c     boxsize     in: double precision(0:nlevels)
c                 Array of boxsizes
c 
c     centers     in: double precision(2,nboxes)
c                 array of centers of boxes
c   
c     iparent     in: integer(nboxes)
c                 iparent(i) is the box number of the parent of
c                 box i
c
c     nchild      in: integer(nboxes)
c                 nchild(i) is the number of children of box i
c
c     ichild      in: integer(4,nboxes)
c                 ichild(j,i) is the box id of the jth child of
c                 box i
c
c     iper        in: integer
c                 flag for periodic implementations. 
c                 Currently not used. Feature under construction.
c
c----------------------------------------------------------------
c     OUTPUT
c     nnbors      out: integer(nboxes)
c                 nnbors(i) is the number of colleague boxes of
c                 box i
c
c     nbors       out: integer(9,nboxes)
c                 nbors(j,i) is the box id of the jth colleague
c                 box of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes
      integer iper
      integer laddr(2,0:nlevels)
      double precision boxsize(0:nlevels)
      double precision centers(2,nboxes)
      integer iparent(nboxes), nchild(nboxes), ichild(4,nboxes)
      integer nnbors(nboxes)
      integer nbors(9,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox


c     Setting parameters for level = 0
      nnbors(1) = 1
      nbors(1,1) = 1
      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = laddr(1,ilev)
         ilastbox = laddr(2,ilev)
c        Loop over all boxes to evaluate neighbors, list1 and updating
c        hunglists of targets

C$OMP PARALLEL DO DEFAULT(SHARED)
C$OMP$PRIVATE(ibox,dad,i,jbox,j,kbox)
         do ibox = ifirstbox,ilastbox
c           Find the parent of the current box         
            dad = iparent(ibox)
c           Loop over the neighbors of the parent box
c           to find out list 1 and list 2
            do i=1,nnbors(dad)
                jbox = nbors(i,dad)
                do j=1,4
c               ichild(j,jbox) is one of the children of the
c               neighbors of the parent of the current
c               box
                   kbox = ichild(j,jbox)
                   if(kbox.gt.0) then
c               Check if kbox is a nearest neighbor or in list 2
                      if((abs(centers(1,kbox)-centers(1,ibox)).le.
     1                   1.05*boxsize(ilev)).and.
     2                   (abs(centers(2,kbox)-centers(2,ibox)).le.
     3                   1.05*boxsize(ilev))) then
                     
                         nnbors(ibox) = nnbors(ibox)+1
                         nbors(nnbors(ibox),ibox) = kbox
                      endif
                   endif
                enddo
            enddo
c           End of computing colleagues of box i
         enddo
C$OMP END PARALLEL DO         
      enddo

      return
      end
c
c
c
c
c--------------------------------------------------------------------      
      subroutine updateflags(curlev,nboxes,nlevels,laddr,nchild,ichild,
     1                    nnbors,nbors,centers,boxsize,iflag)

c      This subroutine is to check the boxes flagged as flag++
c      and determine which of the boxes need refinement. The flag
c      of the box which need refinement is updated to iflag(box)=1
c      and that of the boxes which do not need refinement is
c      updated to iflag(box) = 0
c
c      INPUT arguments
c      curlev         in: integer
c                     the level for which boxes need to be processed
c
c      nboxes         in: integer
c                     total number of boxes
c
c      nlevels        in: integer
c                     total number of levels
c
c      laddr          in: integer(2,0:nlevels)
c                     boxes from laddr(1,ilev) to laddr(2,ilev)
c                     are at level ilev
c
c      nchild         in: integer(nboxes)
c                     nchild(ibox) is the number of children
c                     of box ibox
c
c      ichild         in: integer(4,nboxes)
c                     ichild(j,ibox) is the box id of the jth
c                     child of box ibox
c
c      nnbors         in: integer(nboxes)
c                     nnbors(ibox) is the number of colleagues
c                     of box ibox
c
c      nbors          in: integer(9,nboxes)
c                     nbors(j,ibox) is the jth colleague of box
c                     ibox
c
c      centers        in: double precision(2,nboxes)
c                     x and y coordinates of the box centers
c
c      boxsize        in: double precision(0:nlevels)
c                     boxsize(i) is the size of the box at level i
c
c      iflag          in/out: integer(nboxes)
c                     iflag(ibox)=3 if it is flag++. iflag(ibox) =1
c                     or 0 at the end of routine depending on
c                     whether box needs to be subdivided or not
c
      implicit none
c     Calling sequence variables
      integer curlev, nboxes, nlevels
      integer laddr(2,0:nlevels),nchild(nboxes),ichild(4,nboxes)
      integer nnbors(nboxes), nbors(9,nboxes)
      integer iflag(nboxes)
      double precision centers(2,nboxes),boxsize(0:nlevels)

c     Temporary variables
      integer i,j,k,l,ibox,jbox,kbox,lbox, ict
      double precision distest,xdis,ydis,zdis

      distest = 1.05d0*(boxsize(curlev) + boxsize(curlev+1))/2.0d0
c     Loop over all boxes at the current level     

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(ibox,i,jbox,j,kbox,xdis,ydis)
C$OMP$PRIVATE(ict)
      do ibox = laddr(1,curlev),laddr(2,curlev)
         if(iflag(ibox).eq.3) then
            iflag(ibox) = 0
c           Loop over colleagues of the current box      
            do i=1,nnbors(ibox)
c              Loop over colleagues of flag++ box        
               jbox = nbors(i,ibox)
              
c              Loop over the children of the colleague box
c              Note we do not need to exclude self from
c              the list of colleagues as a self box which
c              is flag++ does not have any children 
c              and will not enter the next loop
               do j=1,4
                  kbox = ichild(j,jbox)
                  if(kbox.gt.0) then
                     if(nchild(kbox).gt.0) then
                        xdis = centers(1,kbox) - centers(1,ibox)
                        ydis = centers(2,kbox) - centers(2,ibox)
                        ict = 0
                        if(abs(xdis).le.distest) ict = ict + 1
                        if(abs(ydis).le.distest) ict = ict + 1
                        if(ict.eq.2) then
                           iflag(ibox) = 1
                           goto 1111
                        endif
                     endif
                  endif
c                 End of looping over the children of the child
c                 of the colleague box
               enddo
c              End of looping over the children of the colleague box       
            enddo
c           End of looping over colleagues            
 1111       continue        
         endif
c        End of testing if the current box needs to checked for         
      enddo
c     End of looping over boxes at the current level      
C$OMP END PARALLEL DO      

      return
      end
c
c
c
c
c
c

      subroutine tree_refine_boxes_flag(iflag,nboxes,
     1  ifirstbox,nbloc,centers,bs,nbctr,nlctr,
     2  ilevel,iparent,nchild,ichild)
      implicit none
      integer nboxes,nbloc,nbctr,nlctr
      real *8 centers(2,nboxes),bs
      integer ilevel(nboxes),iparent(nboxes)
      integer ichild(4,nboxes),nchild(nboxes)
      integer iflag(nboxes)
      integer ifirstbox
      integer, allocatable :: isum(:),itmp(:)

      integer i,ibox,nel0,j,l,jbox,nel1,nbl
      integer ii

      allocate(isum(nbloc),itmp(nbloc))
      do i=1,nbloc
        ibox = ifirstbox+i-1
        itmp(i) = 0
        if(iflag(ibox).gt.0) itmp(i) = 1
      enddo
      if(nbloc.gt.0) call cumsum(nbloc,itmp,isum)
      
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,ibox,nbl,j,jbox,l,ii)
      do i = 1,nbloc
        ibox = ifirstbox + i-1
        if(iflag(ibox).gt.0) then
          nbl = nbctr + (isum(i)-1)*4
          
          nchild(ibox) = 4
          do j=1,4
            ii = 2
            if(j.le.2) ii = 1
            jbox = nbl+j
            centers(1,jbox) = centers(1,ibox)+(-1)**j*bs/2
            centers(2,jbox) = centers(2,ibox)+(-1)**ii*bs/2 
            iparent(jbox) = ibox
            nchild(jbox) = 0
            do l=1,4
              ichild(l,jbox) = -1
            enddo
            ichild(j,ibox) = jbox
            ilevel(jbox) = nlctr+1 
            if(iflag(ibox).eq.1) iflag(jbox) = 3
            if(iflag(ibox).eq.2) iflag(jbox) = 0
          enddo
        endif
      enddo
C$OMP END PARALLEL DO      

      if(nbloc.gt.0) nbctr = nbctr + isum(nbloc)*4


      return
      end
c
c
c
c
c
c
      subroutine print_tree(itree,ltree,nboxes,centers,boxsize,nlevels,
     1   iptr,fname)
c
c        this subroutine writes the tree info to a file
c
c        input arguments:
c          itree - integer (ltree)
c             packed array containing tree info
c          ltree - integer
c            length of itree
c          nboxes - integer
c             number of boxes
c          centers - real *8 (2,nboxes)
c             xy coordinates of box centers in tree hierarchy
c          boxsize - real *8 (0:nlevels)
c             size of box at various levels
c          nlevels - integer
c             number of levels
c          iptr - integer(8)
c            pointer to various arrays inside itree
c          fname - character *
c            file name to which tree info is to be written
c 
c          output
c            file with name fname, which contains the tree info
c            file can be plotted using the python script
c              tree_plot.py containted in src/common
c


      implicit real *8 (a-h,o-z)
      integer itree(ltree),ltree,nboxes,nlevels,iptr(8)
      real *8 centers(2,nboxes),boxsize(0:nlevels)
      character (len=*) fname

      open(unit=33,file=trim(fname))
      

 1111 format(10(2x,e11.5))      

      do ibox=1,nboxes
         if(itree(iptr(4)+ibox-1).eq.0) then
           ilev = itree(iptr(2)+ibox-1)
           bs = boxsize(ilev)
           x1 = centers(1,ibox) - bs/2
           x2 = centers(1,ibox) + bs/2

           y1 = centers(2,ibox) - bs/2
           y2 = centers(2,ibox) + bs/2
           
           write(33,1111) x1,x2,x2,x1,x1,y1,y1,y2,y2,y1
         endif
      enddo

      close(33)

      return
      end
c
c
c
c
c

      subroutine computemnlists(nlevels,nboxes,itree,ltree,
     1   iptr,centers,
     1   boxsize,iper,mnlist1,mnlist2,mnlist3,mnlist4)
c
c        determine maximum number of elements in list1,list2,list3,list4
c
c        NOTE in 2D: we use max values
c
      implicit real *8 (a-h,o-z)
      integer ltree
      integer nlevels,nboxes,itree(ltree),iptr(8)
      real *8 centers(2,nboxes),boxsize(0:nlevels)
      integer mnlist1,mnlist2,mnlist3,mnlist4

      mnlist1 = 13
      mnlist2 = 27
      mnlist3 = 20
      mnlist4 = 5

      return
      end
c
c
c
c
c

      subroutine computelists(nlevels,nboxes,itree,ltree,
     1            iptr,centers,
     2            boxsize,iper,mnlist1,nlist1,list1,
     3            mnlist2,nlist2,list2,
     4            mnlist3,nlist3,list3,
     5            mnlist4,nlist4,list4)
c
ck
c     This subroutine computes the various fmm lists of a given tree
c     structure
c     
c     INPUT arguments
c     nlevels     in: integer
c                 Number of levels
c
c     nboxes      in: integer
c                 Total number of boxes
c
c     itree       in: integer(ltree)
c                   array containing tree info - see start of file
c                   for documentation
c     ltree       in: integer
c                   length of itree array
c 
c     iptr        in: integer(8)
c                   pointer for various arrays in itree
c
c     centers     in: real *8(2,nboxes)
c                 xy coordinates of centers of boxes
c   
c     boxsize     in: real *8(0:nlevels)
c                 Array of boxsizes
c   
c     iper        in: integer
c                 flag for periodic implementations. Currently not used.
c                 Feature under construction
c 
c     mnlist1     in: integer
c                 max number of boxes in list 1 of a box
c
c     mnlist2     in: integer
c                 max number of boxes in list 2 of a box
c
c     mnlist3     in: integer
c                 max number of boxes in list 3 of a box
c
c     mnlist4     in: integer
c                 max number of boxes in list 4 of a box 
c  
c--------------------------------------------------------------
c     OUTPUT arguments:
c     nlist1      out: integer(nboxes)
c                 nlist1(i) is the number of boxes in list 1 
c                 of box i
c
c     list1       out: integer(mnlist1,nboxes)
c                 list1(j,i) is the box id of the jth box in 
c                 list1 of box i
c                      
c     nlist2      out: integer(nboxes)
c                 nlist2(i) is the number of boxes in the list 2
c                 of box i
c 
c     list2       out: integer(mnlist2,nboxes)
c                 list2(j,i) is the box id of the jth box in 
c                 list2 of box i
c
c     nlist3      out: integer(nboxes)
c                 nlist3(i) is the number of boxes in list 3
c                 of box i
c
c     list3       out: integer(mnlist3,nboxes)
c                 list3(j,i) is the box id of the jth box in 
c                 list3 of box i
c                      
c     nlist4      out: integer(nboxes)
c                 nlist4(i) is the number of boxes in the list 2
c                 of box i
c 
c     list4       out: integer(mnlist4,nboxes)
c                 list4(j,i) is the box id of the jth box in 
c                 list4 of box i
c---------------------------------------------------------------
      implicit none
      integer nlevels,nboxes,ltree
      integer iper
      integer itree(ltree),iptr(8)
      real *8 boxsize(0:nlevels)
      real *8 centers(2,nboxes)
      integer mnlist1,mnlist2,mnlist3,mnlist4
      integer nlist1(nboxes), list1(mnlist1,nboxes)
      integer nlist2(nboxes), list2(mnlist2,nboxes)
      integer nlist3(nboxes), list3(mnlist3,nboxes)
      integer nlist4(nboxes), list4(mnlist4,nboxes)

c     Temp variables
      integer ilev,ibox,jbox,kbox,dad
      integer i,j,ifirstbox,ilastbox,ii
      real *8 distest,xdis,ydis

      do i=1,nboxes
        nlist1(i) = 0
        nlist2(i) = 0
        nlist3(i) = 0
        nlist4(i) = 0
      enddo


c     Setting parameters for level = 0
      if(itree(iptr(4)).eq.0) then
         nlist1(1) = 1 
         list1(1,1) = 1
      else
         nlist1(1) = 0
      endif
      nlist2(1) = 0
      nlist3(1) = 0
      nlist4(1) = 0

      





      do ilev = 1,nlevels
c        Find the first and the last box at level ilev      
         ifirstbox = itree(2*ilev+1)
         ilastbox = itree(2*ilev+2)

         do ibox = ifirstbox,ilastbox
            dad = itree(iptr(3)+ibox-1)

c           Loop over the neighbors of the parent box
c           to find out list 2

            do i=1,itree(iptr(6)+dad-1)

c
c                 jbox is a colleague of ibox's parent
c
               jbox = itree(iptr(7) + 9*(dad-1)+i-1)
               do j=1,itree(iptr(4)+jbox-1)

c
c                  kbox is a child of a colleague of ibox's parent
c
                  kbox = itree(iptr(5)+4*(jbox-1)+j-1)
                   
                  if((abs(centers(1,kbox)-centers(1,ibox)).ge.
     1                1.05*boxsize(ilev)).or.
     2                (abs(centers(2,kbox)-centers(2,ibox)).ge.
     3                1.05*boxsize(ilev))) then
                     
                     nlist2(ibox) = nlist2(ibox)+1
                     list2(nlist2(ibox),ibox) = kbox
                  endif
               enddo
            enddo

cc            goto 1120

c           Compute list 1 and list 3 of ibox if ibox is childless
            if(itree(iptr(4)+ibox-1).eq.0) then
c              Loop over all colleagues of ibox              
               do i=1,itree(iptr(6)+ibox-1)
                  jbox = itree(iptr(7)+9*(ibox-1)+i-1)

c                 If the colleague box is childless, then
c                 colleague box is in list 1
                  if(itree(iptr(4)+jbox-1).eq.0) then
                     nlist1(ibox) = nlist1(ibox)+1
                     list1(nlist1(ibox),ibox) = jbox

c                 If colleague box is not childless, then
c                 test to see if children of colleague
c                 box are in list1 or list 3. 

                  else
                     distest = 1.05d0*(boxsize(ilev) + 
     1                                 boxsize(ilev+1))/2.0d0
c                    Loop over children of colleague box              
                     do j=1,itree(iptr(4)+jbox-1)
                        kbox = itree(iptr(5)+4*(jbox-1)+j-1)
                        xdis = dabs(centers(1,kbox)-centers(1,ibox))
                        ydis = dabs(centers(2,kbox)-centers(2,ibox))
c                       Test to see if child of colleague box
c                       is in list1
                        if(xdis.lt.distest.and.ydis.lt.distest) then
                           nlist1(ibox) = nlist1(ibox)+1
                           list1(nlist1(ibox),ibox)=kbox

                           nlist1(kbox) = nlist1(kbox)+1
                           list1(nlist1(kbox),kbox) = ibox

c                       If it is not in list 1 of ibox then it 
c                       is in list3
                        else
                           nlist3(ibox) = nlist3(ibox)+1
                           list3(nlist3(ibox),ibox)=kbox

                           nlist4(kbox) = nlist4(kbox)+1
                           list4(nlist4(kbox),kbox)=ibox
                        endif
c                       End of figuring out whether child 
c                       of colleague box is in list 1 or list3
                     enddo
c                    End of looping over of children of colleague
c                    box
                  endif
c                 End of checking if colleague box has children
               enddo
c              End of looping over colleague boxes
            endif 
c           End of checking of current box is childless


 1120      continue
         enddo
c        End of looping over boxes at level ilev         
      enddo
c     End of looping over levels      

      return
      end

