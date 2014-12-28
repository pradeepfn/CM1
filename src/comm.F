!-----------------------------------------------------------------------
!  message passing routines
!-----------------------------------------------------------------------


      function nabor(i,j,nx,ny)
      implicit none
      integer i,j,nx,ny,nabor
      integer newi,newj

      newi=i
      newj=j

      if ( newi .lt.  1 ) newi = nx
      if ( newi .gt.  nx) newi = 1

      if ( newj .lt.  1 ) newj = ny
      if ( newj .gt.  ny) newj = 1

      nabor = (newi-1) + (newj-1)*nx

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

#ifdef MPI
      subroutine getcorner(s,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke)
      real, intent(inout), dimension(nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
 
      integer k,nn,nr,index
      integer reqs(8),status(MPI_STATUS_SIZE)
      integer tag,count

!-----

      nr = 0

      tag=5001
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          se1(k)=s(ni,1,k)
        enddo
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif
 
!-----

      tag=5002
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          ne1(k)=s(ni,nj,k)
        enddo
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif
 
!-----

      tag=5003
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        nr = nr + 1
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          sw1(k)=s(1,1,k)
        enddo
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif
 
!-----

      tag=5004
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        nr = nr + 1
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          nw1(k)=s(1,nj,k)
        enddo
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif
 
!-----

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1),index,status,ierr)
        nn = nn + 1
      enddo

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          s(0,nj+1,k)=nw2(k)
        enddo
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          s(0,0,k)=sw2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          s(ni+1,nj+1,k)=ne2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          s(ni+1,0,k)=se2(k)
        enddo
      endif

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(8),status,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcornert(t,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: t
      real, intent(inout), dimension(nk+1) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer k,nr,nn,index
      integer :: index_nw,index_sw,index_ne,index_se
      integer reqs(8),status(MPI_STATUS_SIZE)
      integer tag,count

!-----

      nr = 0
      index_nw = -1
      index_sw = -1
      index_ne = -1
      index_se = -1

      tag=5001
      count=nkp1

      if(ibw.eq.0 .and. ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_nw = nr
      endif

      if(ibe.eq.0 .and. ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          se1(k)=t(ni,1,k)
        enddo
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif

!-----

      tag=5002
      count=nkp1

      if(ibw.eq.0 .and. ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_sw = nr
      endif

      if(ibe.eq.0 .and. ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          ne1(k)=t(ni,nj,k)
        enddo
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif

!-----

      tag=5003
      count=nkp1

      if(ibe.eq.0 .and. ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_ne = nr
      endif

      if(ibw.eq.0 .and. ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          sw1(k)=t(1,1,k)
        enddo
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif

!-----

      tag=5004
      count=nkp1

      if(ibe.eq.0 .and. ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(nr),ierr)
        index_se = nr
      endif

      if(ibw.eq.0 .and. ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          nw1(k)=t(1,nj,k)
        enddo
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif

!-----

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1),index,status,ierr)
        nn = nn + 1

      if(index.eq.index_nw)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          t(0,nj+1,k)=nw2(k)
        enddo
      elseif(index.eq.index_sw)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          t(0,0,k)=sw2(k)
        enddo
      elseif(index.eq.index_ne)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          t(ni+1,nj+1,k)=ne2(k)
        enddo
      elseif(index.eq.index_se)then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          t(ni+1,0,k)=se2(k)
        enddo
      endif

      enddo

!-----

      if(ibe.eq.0 .and. ibs.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibe.eq.0 .and. ibn.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
      endif

      if(ibw.eq.0 .and. ibs.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(ibw.eq.0 .and. ibn.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcorneru(u,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real u(ib:ie+1,jb:je,kb:ke)
      real, intent(inout), dimension(nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer k
      integer reqs(8),status(MPI_STATUS_SIZE)
      integer tag,count

!-----

      tag=5001
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          se1(k)=u(ni,1,k)
        enddo
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(1),ierr)
      endif

!-----

      tag=5002
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(4),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          ne1(k)=u(ni,nj,k)
        enddo
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(3),ierr)
      endif

!-----

      tag=5003
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          sw1(k)=u(2,1,k)
        enddo
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif

!-----

      tag=5004
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          nw1(k)=u(2,nj,k)
        enddo
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif

!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          u(0,nj+1,k)=nw2(k)
        enddo
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          u(0,0,k)=sw2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          u(ni+2,nj+1,k)=ne2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          u(ni+2,0,k)=se2(k)
        enddo
      endif

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcornerv(v,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real v(ib:ie,jb:je+1,kb:ke)
      real, intent(inout), dimension(nk) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer k
      integer reqs(8),status(MPI_STATUS_SIZE)
      integer tag,count

!-----

      tag=5001
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          se1(k)=v(ni,2,k)
        enddo
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(1),ierr)
      endif

!-----

      tag=5002
      count=nk

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(4),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          ne1(k)=v(ni,nj,k)
        enddo
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(3),ierr)
      endif

!-----

      tag=5003
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          sw1(k)=v(1,2,k)
        enddo
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif

!-----

      tag=5004
      count=nk

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          nw1(k)=v(1,nj,k)
        enddo
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif

!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          v(0,nj+2,k)=nw2(k)
        enddo
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          v(0,0,k)=sw2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          v(ni+1,nj+2,k)=ne2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nk
          v(ni+1,0,k)=se2(k)
        enddo
      endif

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcornerw(w,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real w(ib:ie,jb:je,kb:ke+1)
      real, intent(inout), dimension(nk+1) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer k
      integer reqs(8),status(MPI_STATUS_SIZE)
      integer tag,count

!-----

      tag=5001
      count=nkp1

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          se1(k)=w(ni,1,k)
        enddo
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(1),ierr)
      endif

!-----

      tag=5002
      count=nkp1

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(4),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          ne1(k)=w(ni,nj,k)
        enddo
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(3),ierr)
      endif

!-----

      tag=5003
      count=nkp1

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          sw1(k)=w(1,1,k)
        enddo
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif

!-----

      tag=5004
      count=nkp1

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          nw1(k)=w(1,nj,k)
        enddo
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif

!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          w(0,nj+1,k)=nw2(k)
        enddo
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          w(0,0,k)=sw2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          w(ni+1,nj+1,k)=ne2(k)
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(k)
        do k=1,nkp1
          w(ni+1,0,k)=se2(k)
        enddo
      endif

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2we_start(s,west,newwest,east,neweast,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real s(ib:ie,jb:je,kb:ke)
      real west(2,nj,nk),newwest(2,nj,nk)
      real east(2,nj,nk),neweast(2,nj,nk)
      integer reqs(4)

      integer i,j,k,nr
      integer tag,count

!------------------------------------------------

      nr = 0

      count=cs2we
      nf=nf+1
      tag=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,2
          west(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(3),ierr)
      endif

!----------

      count=cs2we
      nf=nf+1
      tag=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,2
          east(i,j,k)=s(ni-2+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(4),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2we_end(s,west,newwest,east,neweast,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real s(ib:ie,jb:je,kb:ke)
      real west(2,nj,nk),newwest(2,nj,nk)
      real east(2,nj,nk),neweast(2,nj,nk)
      integer reqs(4)

      integer i,j,k,nn,nr,index
      integer :: index_east,index_west
      integer status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------

      index_east = -1
      index_west = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1),index,status,ierr)
        nn = nn + 1

      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,2
          s(ni+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,2
          s(i-2,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2sn_start(s,south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real s(ib:ie,jb:je,kb:ke)
      real south(ni,2,nk),newsouth(ni,2,nk)
      real north(ni,2,nk),newnorth(ni,2,nk)
      integer reqs(4)

      integer i,j,k,nr
      integer tag,count

!----------

      nr = 0

      count=cs2sn
      nf=nf+1
      tag=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,2
        do i=1,ni
          north(i,j,k)=s(i,nj-2+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(3),ierr)
      endif

!----------

      count=cs2sn
      nf=nf+1
      tag=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,2
        do i=1,ni
          south(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(4),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2sn_end(s,south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real s(ib:ie,jb:je,kb:ke)
      real south(ni,2,nk),newsouth(ni,2,nk)
      real north(ni,2,nk),newnorth(ni,2,nk)
      integer reqs(4)

      integer i,j,k,nn,nr,index
      integer :: index_south,index_north
      integer status(MPI_STATUS_SIZE)

!----------

      index_south = -1
      index_north = -1

      nr = 0
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1),index,status,ierr)
        nn = nn + 1

      if(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,2
        do i=1,ni
          s(i,j-2,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,2
        do i=1,ni
          s(i,nj+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

!----------

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_3s_start(s,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke)
      real west(3,nj,nk),newwest(3,nj,nk)
      real east(3,nj,nk),neweast(3,nj,nk)
      real south(ni,3,nk),newsouth(ni,3,nk)
      real north(ni,3,nk),newnorth(ni,3,nk)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag,count
 
!------------------------------------------------

      nr = 0

      count=cs3we
      nf=nf+1
      tag=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          west(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif
 
!----------
 
      count=cs3we
      nf=nf+1
      tag=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          east(i,j,k)=s(ni-3+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(6),ierr)
      endif
 
!----------
 
      count=cs3sn
      nf=nf+1
      tag=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          north(i,j,k)=s(i,nj-3+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif
 
!----------
 
      count=cs3sn
      nf=nf+1
      tag=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          south(i,j,k)=s(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(8),ierr)
      endif
 
!----------
 
      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3s_end(s,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke)
      real west(3,nj,nk),newwest(3,nj,nk)
      real east(3,nj,nk),neweast(3,nj,nk)
      real south(ni,3,nk),newsouth(ni,3,nk)
      real north(ni,3,nk),newnorth(ni,3,nk)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

    nn = 1
    do while( nn .le. nr )
      call MPI_WAITANY(nr,reqs(1),index,status,ierr)
      nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          s(ni+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          s(i-3,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          s(i,j-3,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          s(i,nj+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      endif

    enddo

      if(timestats.ge.1) time_mps2=time_mps2+mytime()
 
!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then

          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(1,0,k)
            enddo
          endif

          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(1,nj+1,k)
            enddo
          endif

        endif

        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni+1,nj,k)
            enddo
          endif
 
        endif
 
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!-----------------------------------------------------------
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_3t_start(t,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real t(ib:ie,jb:je,kb:ke+1)
      real west(3,nj,nk+1),newwest(3,nj,nk+1)
      real east(3,nj,nk+1),neweast(3,nj,nk+1)
      real south(ni,3,nk+1),newsouth(ni,3,nk+1)
      real north(ni,3,nk+1),newnorth(ni,3,nk+1)
      integer reqs(8)

      integer i,j,k,nr
      integer tag,count

!------------------------------------------------

      nr = 0

      count=ct3we
      nf=nf+1
      tag=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,nj
        do i=1,3
          west(i,j,k)=t(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif

!----------

      count=ct3we
      nf=nf+1
      tag=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,nj
        do i=1,3
          east(i,j,k)=t(ni-3+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(6),ierr)
      endif

!----------

      count=ct3sn
      nf=nf+1
      tag=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,3
        do i=1,ni
          north(i,j,k)=t(i,nj-3+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif

!----------

      count=ct3sn
      nf=nf+1
      tag=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,3
        do i=1,ni
          south(i,j,k)=t(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(8),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_3t_end(t,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real t(ib:ie,jb:je,kb:ke+1)
      real west(3,nj,nk+1),newwest(3,nj,nk+1)
      real east(3,nj,nk+1),neweast(3,nj,nk+1)
      real south(ni,3,nk+1),newsouth(ni,3,nk+1)
      real north(ni,3,nk+1),newnorth(ni,3,nk+1)
      integer reqs(8)

      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1),index,status,ierr)
        nn = nn + 1

      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,nj
        do i=1,3
          t(ni+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,nj
        do i=1,3
          t(i-3,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,3
        do i=1,ni
          t(i,j-3,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nkp1
        do j=1,3
        do i=1,ni
          t(i,nj+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner

     if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then

        if(ibw.eq.1)then

          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,0,k)=t(1,0,k)
            enddo
          endif

          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,nj+1,k)=t(1,nj+1,k)
            enddo
          endif

        endif

        if(ibe.eq.1)then

          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,0,k)=t(ni,0,k)
            enddo
          endif

          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,nj+1,k)=t(ni,nj+1,k)
            enddo
          endif

        endif

      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then

        if(ibs.eq.1)then

          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,0,k)=t(0,1,k)
            enddo
          endif

          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,0,k)=t(ni+1,1,k)
            enddo
          endif

        endif

        if(ibn.eq.1)then

          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,nj+1,k)=t(0,nj,k)
            enddo
          endif

          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,nj+1,k)=t(ni+1,nj,k)
            enddo
          endif

        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3u_start(u,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real u(ib:ie+1,jb:je,kb:ke)
      real west(3,nj,nk),newwest(3,nj,nk)
      real east(3,nj,nk),neweast(3,nj,nk)
      real south(ni+1,3,nk),newsouth(ni+1,3,nk)
      real north(ni+1,3,nk),newnorth(ni+1,3,nk)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag,count
 
!------------------------------------------------

      nr = 0

      nu=nu+1
      tag=1000+nu
      count=cs3we

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          west(i,j,k)=u(i+1,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif
 
!----------
 
      nu=nu+1
      tag=1000+nu
      count=cs3we
 
      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          east(i,j,k)=u(ni-3+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(6),ierr)
      endif

!----------
 
      nu=nu+1
      tag=1000+nu
      count=cu3sn

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni+1
          south(i,j,k)=u(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif
 
!----------
 
      nu=nu+1
      tag=1000+nu
      count=cu3sn

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni+1
          north(i,j,k)=u(i,nj-3+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(8),ierr)
      endif
 
!----------

      if(timestats.ge.1) time_mpu1=time_mpu1+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3u_end(u,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real u(ib:ie+1,jb:je,kb:ke)
      real west(3,nj,nk),newwest(3,nj,nk)
      real east(3,nj,nk),neweast(3,nj,nk)
      real south(ni+1,3,nk),newsouth(ni+1,3,nk)
      real north(ni+1,3,nk),newnorth(ni+1,3,nk)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer status(MPI_STATUS_SIZE)

!----------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1),index,status,ierr)
        nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          u(ni+1+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,3
          u(i-3,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni+1
          u(i,nj+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni+1
          u(i,j-3,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mpu2=time_mpu2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,0,k)=u(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,nj+1,k)=u(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,0,k)=u(ni+1,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,nj+1,k)=u(ni+1,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,0,k)=u(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,0,k)=u(ni+2,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(0,nj+1,k)=u(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              u(ni+2,nj+1,k)=u(ni+2,nj,k)
            enddo
          endif
 
        endif
 
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
      endif
 
      if(timestats.ge.1) time_mpu2=time_mpu2+mytime()
 
!----------
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3v_start(v,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real v(ib:ie,jb:je+1,kb:ke)
      real west(3,nj+1,nk),newwest(3,nj+1,nk)
      real east(3,nj+1,nk),neweast(3,nj+1,nk)
      real south(ni,3,nk),newsouth(ni,3,nk)
      real north(ni,3,nk),newnorth(ni,3,nk)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag,count
 
!------------------------------------------------

      nr = 0

      nv=nv+1
      tag=2000+nv
      count=cv3we

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,3
          west(i,j,k)=v(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif
 
!----------
 
      nv=nv+1
      tag=2000+nv
      count=cv3we
 
      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,3
          east(i,j,k)=v(ni-3+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(6),ierr)
      endif

!----------
 
      nv=nv+1
      tag=2000+nv
      count=cs3sn
 
      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          south(i,j,k)=v(i,j+1,k)
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif
 
!----------
 
      nv=nv+1
      tag=2000+nv
      count=cs3sn
 
      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          north(i,j,k)=v(i,nj-3+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(8),ierr)
      endif
 
!----------

      if(timestats.ge.1) time_mpv1=time_mpv1+mytime()
 
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3v_end(v,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real v(ib:ie,jb:je+1,kb:ke)
      real west(3,nj+1,nk),newwest(3,nj+1,nk)
      real east(3,nj+1,nk),neweast(3,nj+1,nk)
      real south(ni,3,nk),newsouth(ni,3,nk)
      real north(ni,3,nk),newnorth(ni,3,nk)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer status(MPI_STATUS_SIZE)
 
!--------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1),index,status,ierr)
        nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,3
          v(ni+i,j,k)=neweast(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,3
          v(i-3,j,k)=newwest(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          v(i,nj+1+j,k)=newnorth(i,j,k)
        enddo
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,3
        do i=1,ni
          v(i,j-3,k)=newsouth(i,j,k)
        enddo
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mpv2=time_mpv2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,0,k)=v(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,nj+2,k)=v(1,nj+2,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,0,k)=v(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,nj+2,k)=v(ni,nj+2,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,0,k)=v(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,0,k)=v(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(0,nj+2,k)=v(0,nj+1,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              v(ni+1,nj+2,k)=v(ni+1,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!--------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
      endif
 
      if(timestats.ge.1) time_mpv2=time_mpv2+mytime()
 
!----------
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3w_start(w,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real w(ib:ie,jb:je,kb:ke+1)
      real west(3,nj,nk-1),newwest(3,nj,nk-1)
      real east(3,nj,nk-1),neweast(3,nj,nk-1)
      real south(ni,3,nk-1),newsouth(ni,3,nk-1)
      real north(ni,3,nk-1),newnorth(ni,3,nk-1)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag,count
 
!------------------------------------------------

      nr = 0

      nw=nw+1
      tag=3000+nw
      count=cw3we

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,3
          west(i,j,k-1)=w(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,    &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif
 
!----------
 
      nw=nw+1
      tag=3000+nw
      count=cw3we
 
      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,3
          east(i,j,k-1)=w(ni-3+i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(6),ierr)
      endif

!----------

      nw=nw+1
      tag=3000+nw
      count=cw3sn
 
      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,3
        do i=1,ni
          south(i,j,k-1)=w(i,j,k)
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif
 
!----------
 
      nw=nw+1
      tag=3000+nw
      count=cw3sn
 
      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,3
        do i=1,ni
          north(i,j,k-1)=w(i,nj-3+j,k)
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(8),ierr)
      endif
 
!----------

      if(timestats.ge.1) time_mpw1=time_mpw1+mytime()
 
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3w_end(w,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real w(ib:ie,jb:je,kb:ke+1)
      real west(3,nj,nk-1),newwest(3,nj,nk-1)
      real east(3,nj,nk-1),neweast(3,nj,nk-1)
      real south(ni,3,nk-1),newsouth(ni,3,nk-1)
      real north(ni,3,nk-1),newnorth(ni,3,nk-1)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer status(MPI_STATUS_SIZE)

!--------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1),index,status,ierr)
        nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,3
          w(ni+i,j,k)=neweast(i,j,k-1)
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,3
          w(i-3,j,k)=newwest(i,j,k-1)
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,3
        do i=1,ni
          w(i,nj+j,k)=newnorth(i,j,k-1)
        enddo
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,3
        do i=1,ni
          w(i,j-3,k)=newsouth(i,j,k-1)
        enddo
        enddo
        enddo
      endif

      enddo
 
      if(timestats.ge.1) time_mpw2=time_mpw2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,0,k)=w(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,nj+1,k)=w(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,0,k)=w(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,nj+1,k)=w(ni,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,0,k)=w(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,0,k)=w(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,nj+1,k)=w(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,nj+1,k)=w(ni+1,nj,k)
            enddo
          endif
 
        endif
 
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!--------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
      endif

      if(timestats.ge.1) time_mpw2=time_mpw2+mytime()
 
!----------
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1s_start(s,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke)
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)
 
      integer i,j,k,nr
      integer tag,count
 
!------------------------------------------------

      nr = 0

      nf=nf+1
      tag=nf
      count=cs1we

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          west(j,k)=s(1,j,k)
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif
 
!----------

      nf=nf+1
      tag=nf
      count=cs1we

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          east(j,k)=s(ni,j,k)
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(6),ierr)
      endif
 
!----------

      nf=nf+1
      tag=nf
      count=cs1sn

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          north(i,k)=s(i,nj,k)
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif
 
!----------

      nf=nf+1
      tag=nf
      count=cs1sn

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          south(i,k)=s(i,1,k)
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(8),ierr)
      endif

!----------
 
      if(timestats.ge.1) time_mps1=time_mps1+mytime()
 
      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1s_end(s,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke)
      real west(nj,nk),newwest(nj,nk)
      real east(nj,nk),neweast(nj,nk)
      real south(ni,nk),newsouth(ni,nk)
      real north(ni,nk),newnorth(ni,nk)
      integer reqs(8)
 
      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer status(MPI_STATUS_SIZE)

!---------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1),index,status,ierr)
        nn = nn + 1
 
      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          s(ni+1,j,k)=neweast(j,k)
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nk
        do j=1,nj
          s(0,j,k)=newwest(j,k)
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          s(i,0,k)=newsouth(i,k)
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nk
        do i=1,ni
          s(i,nj+1,k)=newnorth(i,k)
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,0,k)=s(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,0,k)=s(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(0,nj+1,k)=s(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              s(ni+1,nj+1,k)=s(ni+1,nj,k)
            enddo
          endif
 
        endif
 
      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()
 
!----------
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1t_start(t,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real t(ib:ie,jb:je,kb:ke+1)
      real west(nj,nk+1),newwest(nj,nk+1)
      real east(nj,nk+1),neweast(nj,nk+1)
      real south(ni,nk+1),newsouth(ni,nk+1)
      real north(ni,nk+1),newnorth(ni,nk+1)
      integer reqs(8)

      integer i,j,k,nr
      integer tag,count

!------------------------------------------------

      nr = 0

      nf=nf+1
      tag=nf
      count=ct1we

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nkp1
        do j=1,nj
          west(j,k)=t(1,j,k)
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif

!----------

      nf=nf+1
      tag=nf
      count=ct1we

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nkp1
        do j=1,nj
          east(j,k)=t(ni,j,k)
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(6),ierr)
      endif

!----------

      nf=nf+1
      tag=nf
      count=ct1sn

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nkp1
        do i=1,ni
          north(i,k)=t(i,nj,k)
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif

!----------

      nf=nf+1
      tag=nf
      count=ct1sn

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nkp1
        do i=1,ni
          south(i,k)=t(i,1,k)
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(8),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_1t_end(t,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real t(ib:ie,jb:je,kb:ke+1)
      real west(nj,nk+1),newwest(nj,nk+1)
      real east(nj,nk+1),neweast(nj,nk+1)
      real south(ni,nk+1),newsouth(ni,nk+1)
      real north(ni,nk+1),newnorth(ni,nk+1)
      integer reqs(8)

      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer status(MPI_STATUS_SIZE)

!---------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1),index,status,ierr)
        nn = nn + 1

      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nkp1
        do j=1,nj
          t(ni+1,j,k)=neweast(j,k)
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=1,nkp1
        do j=1,nj
          t(0,j,k)=newwest(j,k)
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nkp1
        do i=1,ni
          t(i,0,k)=newsouth(i,k)
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=1,nkp1
        do i=1,ni
          t(i,nj+1,k)=newnorth(i,k)
        enddo
        enddo
      endif

      enddo

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner

      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then

        if(ibw.eq.1)then

          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,0,k)=t(1,0,k)
            enddo
          endif

          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,nj+1,k)=t(1,nj+1,k)
            enddo
          endif

        endif

        if(ibe.eq.1)then

          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,0,k)=t(ni,0,k)
            enddo
          endif

          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,nj+1,k)=t(ni,nj+1,k)
            enddo
          endif

        endif

      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then

        if(ibs.eq.1)then

          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,0,k)=t(0,1,k)
            enddo
          endif

          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,0,k)=t(ni+1,1,k)
            enddo
          endif

        endif

        if(ibn.eq.1)then

          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(0,nj+1,k)=t(0,nj,k)
            enddo
          endif

          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nkp1
              t(ni+1,nj+1,k)=t(ni+1,nj,k)
            enddo
          endif

        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1w_start(w,ww1,ww2,we1,we2,   &
                                 ws1,ws2,wn1,wn2,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real w(ib:ie,jb:je,kb:ke+1)
      real ww1(nj,nk-1),ww2(nj,nk-1)
      real we1(nj,nk-1),we2(nj,nk-1)
      real ws1(ni,nk-1),ws2(ni,nk-1)
      real wn1(ni,nk-1),wn2(ni,nk-1)
      integer reqs(8)
 
      integer i,j,k
      integer tag,count

!-----

      nw=nw+1
      tag=3000+nw
      count=cw1we
 
        ! receive east
        if(ibe.eq.0)then
          call mpi_irecv(we2,count,MPI_REAL,myeast,tag,   &
                        MPI_COMM_WORLD,reqs(2),ierr)
        endif

        ! send west
        if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=2,nk
          do j=1,nj
            ww1(j,k-1)=w(1,j,k)
          enddo
          enddo
          call mpi_isend(ww1,count,MPI_REAL,mywest,tag,   &
                         MPI_COMM_WORLD,reqs(1),ierr)
        endif

!-----

      nw=nw+1
      tag=3000+nw
      count=cw1we
 
        ! receive west
        if(ibw.eq.0)then
          call mpi_irecv(ww2,count,MPI_REAL,mywest,tag,   &
                        MPI_COMM_WORLD,reqs(4),ierr)
        endif

        ! send east
        if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=2,nk
          do j=1,nj
            we1(j,k-1)=w(ni,j,k)
          enddo
          enddo
          call mpi_isend(we1,count,MPI_REAL,myeast,tag,   &
                         MPI_COMM_WORLD,reqs(3),ierr)
        endif

!-----

      nw=nw+1
      tag=3000+nw
      count=cw1sn
 
        ! receive north
        if(ibn.eq.0)then
          call mpi_irecv(wn2,count,MPI_REAL,mynorth,tag,   &
                        MPI_COMM_WORLD,reqs(6),ierr)
        endif

        ! send south
        if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=2,nk
          do i=1,ni
            ws1(i,k-1)=w(i,1,k)
          enddo
          enddo
          call mpi_isend(ws1,count,MPI_REAL,mysouth,tag,   &
                         MPI_COMM_WORLD,reqs(5),ierr)
        endif

!-----

      nw=nw+1
      tag=3000+nw
      count=cw1sn
 
        ! receive south
        if(ibs.eq.0)then
          call mpi_irecv(ws2,count,MPI_REAL,mysouth,tag,   &
                        MPI_COMM_WORLD,reqs(8),ierr)
        endif

        ! send north
        if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=2,nk
          do i=1,ni
            wn1(i,k-1)=w(i,nj,k)
          enddo
          enddo
          call mpi_isend(wn1,count,MPI_REAL,mynorth,tag,   &
                         MPI_COMM_WORLD,reqs(7),ierr)
        endif

!-----

      if(timestats.ge.1) time_mpw1=time_mpw1+mytime()
 
      return
      end

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_1w_end(w,ww1,ww2,we1,we2,   &
                               ws1,ws2,wn1,wn2,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real w(ib:ie,jb:je,kb:ke+1)
      real ww1(nj,nk-1),ww2(nj,nk-1)
      real we1(nj,nk-1),we2(nj,nk-1)
      real ws1(ni,nk-1),ws2(ni,nk-1)
      real wn1(ni,nk-1),wn2(ni,nk-1)
      integer reqs(8)
 
      integer i,j,k
      integer status(MPI_STATUS_SIZE)
 
!-----

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(2),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=2,nk
        do j=1,nj
          w(ni+1,j,k)=we2(j,k-1)
        enddo
        enddo
      endif
 
      if(ibw.eq.0)then
        call MPI_WAIT (reqs(4),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(j,k)
        do k=2,nk
        do j=1,nj
          w(0,j,k)=ww2(j,k-1)
        enddo
        enddo
      endif
 
!-----
 
      if(ibn.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=2,nk
        do i=1,ni
          w(i,nj+1,k)=wn2(i,k-1)
        enddo
        enddo
      endif
 
      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,k)
        do k=2,nk
        do i=1,ni
          w(i,0,k)=ws2(i,k-1)
        enddo
        enddo
      endif

!-----

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(timestats.ge.1) time_mpw2=time_mpw2+mytime()

!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then
 
          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,0,k)=w(1,0,k)
            enddo
          endif
 
          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,nj+1,k)=w(1,nj+1,k)
            enddo
          endif
 
        endif
 
        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,0,k)=w(ni,0,k)
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,nj+1,k)=w(ni,nj+1,k)
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,0,k)=w(0,1,k)
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,0,k)=w(ni+1,1,k)
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(0,nj+1,k)=w(0,nj,k)
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=2,nk
              w(ni+1,nj+1,k)=w(ni+1,nj,k)
            enddo
          endif
 
        endif
 
      endif

!-----------------------------------------------------------
!  Mirror b.c. patch
!
!      if(ibn.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=2,nk
!        do i=0,ni+1
!          w(i,nj+1,k)=w(i,nj  ,k)
!          w(i,nj+2,k)=w(i,nj-1,k)
!          w(i,nj+3,k)=w(i,nj-2,k)
!        enddo
!        enddo
!      endif
!
!      if(ibw.eq.1)then
!!$omp parallel do default(shared)   &
!!$omp private(i,k)
!        do k=2,nk
!        do j=0,nj+1
!          w(-2,j,k)=w(3,j,k)
!          w(-1,j,k)=w(2,j,k)
!          w( 0,j,k)=w(1,j,k)
!        enddo
!        enddo
!      endif
!
!-----------------------------------------------------------
 
      if(timestats.ge.1) time_bc=time_bc+mytime()
 
!----------

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2d_start(s,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real, dimension(ib:ie,jb:je) :: s
      real, dimension(3,nj) :: west,newwest,east,neweast
      real, dimension(ni,3) :: south,newsouth,north,newnorth
      integer reqs(8)

      integer i,j,nr
      integer tag,count

!------------------------------------------------

      nr = 0

      count=3*nj
      nf=nf+1
      tag=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,3
          west(i,j)=s(i,j)
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif

!----------

      count=3*nj
      nf=nf+1
      tag=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,3
          east(i,j)=s(ni-3+i,j)
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(6),ierr)
      endif

!----------

      count=ni*3
      nf=nf+1
      tag=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,3
        do i=1,ni
          north(i,j)=s(i,nj-3+j)
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif

!----------

      count=ni*3
      nf=nf+1
      tag=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,3
        do i=1,ni
          south(i,j)=s(i,j)
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(8),ierr)
      endif

!----------

      if(timestats.ge.1) time_mpq1=time_mpq1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2dew_end(s,west,newwest,east,neweast,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real, dimension(ib:ie,jb:je) :: s
      real, dimension(3,nj) :: west,newwest,east,neweast
      integer reqs(8)

      integer i,j,nn,nr,index
      integer status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
      endif
      if(ibw.eq.0)then
        nr = nr + 1
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(1),index,status,ierr)
        nn = nn + 1
      enddo

      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,3
          s(ni+i,j)=neweast(i,j)
        enddo
        enddo
      endif

      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,3
          s(i-3,j)=newwest(i,j)
        enddo
        enddo
      endif

!----------

      if(timestats.ge.1) time_mpq2=time_mpq2+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2dns_end(s,south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real, dimension(ib:ie,jb:je) :: s
      real, dimension(ni,3) :: south,newsouth,north,newnorth
      integer reqs(8)

      integer i,j,nn,nr1,nr,index
      integer status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------

      nr1 = 0
      if(ibe.eq.0)then
        nr1 = nr1 + 1
      endif
      if(ibw.eq.0)then
        nr1 = nr1 + 1
      endif

      nr = 0
      if(ibs.eq.0)then
        nr = nr + 1
      endif
      if(ibn.eq.0)then
        nr = nr + 1
      endif

      nn = 1
      do while( nn .le. nr )
        call MPI_WAITANY(nr,reqs(nr1+1),index,status,ierr)
        nn = nn + 1
      enddo

      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,3
        do i=1,ni
          s(i,j-3)=newsouth(i,j)
        enddo
        enddo
      endif

      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,3
        do i=1,ni
          s(i,nj+j)=newnorth(i,j)
        enddo
        enddo
      endif

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
      endif

!----------

      if(timestats.ge.1) time_mpq2=time_mpq2+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_2d_corner(s)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real, dimension(ib:ie,jb:je) :: s

      integer reqs(8),status(MPI_STATUS_SIZE)
      real, dimension(3,3) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      integer :: i,j,tag,count

!-----

      tag=5001
      count=9

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(nw2,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(2),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        do j=1,3
        do i=1,3
          se1(i,j)=s(ni-3+i,j)
        enddo
        enddo
        call mpi_isend(se1,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(1),ierr)
      endif

!-----

      tag=5002
      count=9

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(sw2,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(4),ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        do j=1,3
        do i=1,3
          ne1(i,j)=s(ni-3+i,nj-3+j)
        enddo
        enddo
        call mpi_isend(ne1,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(3),ierr)
      endif

!-----

      tag=5003
      count=9

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call mpi_irecv(ne2,count,MPI_REAL,myne,tag,MPI_COMM_WORLD,   &
                       reqs(6),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        do j=1,3
        do i=1,3
          sw1(i,j)=s(i,j)
        enddo
        enddo
        call mpi_isend(sw1,count,MPI_REAL,mysw,tag,MPI_COMM_WORLD,   &
                       reqs(5),ierr)
      endif

!-----

      tag=5004
      count=9

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call mpi_irecv(se2,count,MPI_REAL,myse,tag,MPI_COMM_WORLD,   &
                       reqs(8),ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        do j=1,3
        do i=1,3
          nw1(i,j)=s(i,nj-3+j)
        enddo
        enddo
        call mpi_isend(nw1,count,MPI_REAL,mynw,tag,MPI_COMM_WORLD,   &
                       reqs(7),ierr)
      endif

!-----

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(2),status,ierr)
        do j=1,3
        do i=1,3
          s(-3+i,nj+j)=nw2(i,j)
        enddo
        enddo
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(4),status,ierr)
        do j=1,3
        do i=1,3
          s(-3+i,-3+j)=sw2(i,j)
        enddo
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(6),status,ierr)
        do j=1,3
        do i=1,3
          s(ni+i,nj+j)=ne2(i,j)
        enddo
        enddo
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(8),status,ierr)
        do j=1,3
        do i=1,3
          s(ni+i,-3+j)=se2(i,j)
        enddo
        enddo
      endif

!-----

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(1),status,ierr)
      endif

      if((ibe.eq.0.or.ibe.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(3),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibs.eq.0.or.ibs.eq.2))then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if((ibw.eq.0.or.ibw.eq.2) .and. (ibn.eq.0.or.ibn.eq.2))then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

!-----

      if(timestats.ge.1) time_mptk1=time_mptk1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_3r_start(th,pp,west,newwest,east,neweast,   &
                                     south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real th(ib:ie,jb:je,kb:ke)
      real pp(ib:ie,jb:je,kb:ke)
      real west(3,nj,nk,2),newwest(3,nj,nk,2)
      real east(3,nj,nk,2),neweast(3,nj,nk,2)
      real south(ni,3,nk,2),newsouth(ni,3,nk,2)
      real north(ni,3,nk,2),newnorth(ni,3,nk,2)
      integer reqs(8)

      integer i,j,k,nr
      integer tag,count

!------------------------------------------------

      nr = 0

      count=cs3we*2
      nf=nf+1
      tag=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,3
          west(i,j,k,1)=th(i,j,k)
        enddo
        enddo
        do j=1,nj
        do i=1,3
          west(i,j,k,2)=pp(i,j,k)
        enddo
        enddo
      enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif

!----------

      count=cs3we*2
      nf=nf+1
      tag=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)

      endif

      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,3
          east(i,j,k,1)=th(ni-3+i,j,k)
        enddo
        enddo
        do j=1,nj
        do i=1,3
          east(i,j,k,2)=pp(ni-3+i,j,k)
        enddo
        enddo
      enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(6),ierr)
      endif

!----------

      count=cs3sn*2
      nf=nf+1
      tag=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,3
        do i=1,ni
          north(i,j,k,1)=th(i,nj-3+j,k)
        enddo
        enddo
        do j=1,3
        do i=1,ni
          north(i,j,k,2)=pp(i,nj-3+j,k)
        enddo
        enddo
      enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif

!----------

      count=cs3sn*2
      nf=nf+1
      tag=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,3
        do i=1,ni
          south(i,j,k,1)=th(i,j,k)
        enddo
        enddo
        do j=1,3
        do i=1,ni
          south(i,j,k,2)=pp(i,j,k)
        enddo
        enddo
      enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(8),ierr)
      endif

!----------

      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_3r_end(th,pp,west,newwest,east,neweast,   &
                                   south,newsouth,north,newnorth,reqs)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'

      real th(ib:ie,jb:je,kb:ke)
      real pp(ib:ie,jb:je,kb:ke)
      real west(3,nj,nk,2),newwest(3,nj,nk,2)
      real east(3,nj,nk,2),neweast(3,nj,nk,2)
      real south(ni,3,nk,2),newsouth(ni,3,nk,2)
      real north(ni,3,nk,2),newnorth(ni,3,nk,2)
      integer reqs(8)

      integer i,j,k,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

    nn = 1
    do while( nn .le. nr )
      call MPI_WAITANY(nr,reqs(1),index,status,ierr)
      nn = nn + 1

    if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,3
          th(ni+i,j,k)=neweast(i,j,k,1)
        enddo
        enddo
        do j=1,nj
        do i=1,3
          pp(ni+i,j,k)=neweast(i,j,k,2)
        enddo
        enddo
      enddo
    elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,3
          th(i-3,j,k)=newwest(i,j,k,1)
        enddo
        enddo
        do j=1,nj
        do i=1,3
          pp(i-3,j,k)=newwest(i,j,k,2)
        enddo
        enddo
      enddo
    elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,3
        do i=1,ni
          th(i,j-3,k)=newsouth(i,j,k,1)
        enddo
        enddo
        do j=1,3
        do i=1,ni
          pp(i,j-3,k)=newsouth(i,j,k,2)
        enddo
        enddo
      enddo
    elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,3
        do i=1,ni
          th(i,nj+j,k)=newnorth(i,j,k,1)
        enddo
        enddo
        do j=1,3
        do i=1,ni
          pp(i,nj+j,k)=newnorth(i,j,k,2)
        enddo
        enddo
      enddo
    else
      print *,'  unknown index '
      print *,'  myid,index = ',myid,index
      call stopcm1
    endif

    enddo

    if(timestats.ge.1) time_mps2=time_mps2+mytime()

!----------
!  patch for corner

      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then

        if(ibw.eq.1)then

          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              th(0,0,k)=th(1,0,k)
              pp(0,0,k)=pp(1,0,k)
            enddo
          endif

          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              th(0,nj+1,k)=th(1,nj+1,k)
              pp(0,nj+1,k)=pp(1,nj+1,k)
            enddo
          endif

        endif

        if(ibe.eq.1)then

          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              th(ni+1,0,k)=th(ni,0,k)
              pp(ni+1,0,k)=pp(ni,0,k)
            enddo
          endif

          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              th(ni+1,nj+1,k)=th(ni,nj+1,k)
              pp(ni+1,nj+1,k)=pp(ni,nj+1,k)
            enddo
          endif

        endif

      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then

        if(ibs.eq.1)then

          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              th(0,0,k)=th(0,1,k)
              pp(0,0,k)=pp(0,1,k)
            enddo
          endif

          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              th(ni+1,0,k)=th(ni+1,1,k)
              pp(ni+1,0,k)=pp(ni+1,1,k)
            enddo
          endif

        endif

        if(ibn.eq.1)then

          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              th(0,nj+1,k)=th(0,nj,k)
              pp(0,nj+1,k)=pp(0,nj,k)
            enddo
          endif

          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k)
            do k=1,nk
              th(ni+1,nj+1,k)=th(ni+1,nj,k)
              pp(ni+1,nj+1,k)=pp(ni+1,nj,k)
            enddo
          endif

        endif

      endif

      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!-----------------------------------------------------------

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine comm_3q_start(s,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke,numq)
      real west(3,nj,nk,numq),newwest(3,nj,nk,numq)
      real east(3,nj,nk,numq),neweast(3,nj,nk,numq)
      real south(ni,3,nk,numq),newsouth(ni,3,nk,numq)
      real north(ni,3,nk,numq),newnorth(ni,3,nk,numq)
      integer reqs(8)
 
      integer i,j,k,n,nr
      integer tag,count
 
!------------------------------------------------

      nr = 0

      count=cs3we*numq
      nf=nf+1
      tag=nf

      ! receive east
      if(ibe.eq.0)then
        nr = nr + 1
        call mpi_irecv(neweast,count,MPI_REAL,myeast,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif

      ! send west
      if(ibw.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,n)
        do k=1,nk
        do n=1,numq
        do j=1,nj
        do i=1,3
          west(i,j,k,n)=s(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
        call mpi_isend(west,count,MPI_REAL,mywest,tag,   &
                       MPI_COMM_WORLD,reqs(5),ierr)
      endif
 
!----------
 
      count=cs3we*numq
      nf=nf+1
      tag=nf

      ! receive west
      if(ibw.eq.0)then
        nr = nr + 1
        call mpi_irecv(newwest,count,MPI_REAL,mywest,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send east
      if(ibe.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,n)
        do k=1,nk
        do n=1,numq
        do j=1,nj
        do i=1,3
          east(i,j,k,n)=s(ni-3+i,j,k,n)
        enddo
        enddo
        enddo
        enddo
        call mpi_isend(east,count,MPI_REAL,myeast,tag,   &
                       MPI_COMM_WORLD,reqs(6),ierr)
      endif
 
!----------
 
      count=cs3sn*numq
      nf=nf+1
      tag=nf

      ! receive south
      if(ibs.eq.0)then
        nr = nr + 1
        call mpi_irecv(newsouth,count,MPI_REAL,mysouth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send north
      if(ibn.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,n)
        do k=1,nk
        do n=1,numq
        do j=1,3
        do i=1,ni
          north(i,j,k,n)=s(i,nj-3+j,k,n)
        enddo
        enddo
        enddo
        enddo
        call mpi_isend(north,count,MPI_REAL,mynorth,tag,   &
                       MPI_COMM_WORLD,reqs(7),ierr)
      endif
 
!----------
 
      count=cs3sn*numq
      nf=nf+1
      tag=nf

      ! receive north
      if(ibn.eq.0)then
        nr = nr + 1
        call mpi_irecv(newnorth,count,MPI_REAL,mynorth,tag,   &
                      MPI_COMM_WORLD,reqs(nr),ierr)
      endif
 
      ! send south
      if(ibs.eq.0)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,n)
        do k=1,nk
        do n=1,numq
        do j=1,3
        do i=1,ni
          south(i,j,k,n)=s(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
        call mpi_isend(south,count,MPI_REAL,mysouth,tag,   &
                       MPI_COMM_WORLD,reqs(8),ierr)
      endif
 
!----------
 
      if(timestats.ge.1) time_mps1=time_mps1+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 
      subroutine comm_3q_end(s,west,newwest,east,neweast,   &
                               south,newsouth,north,newnorth,reqs)
      implicit none
 
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'mpif.h'
 
      real s(ib:ie,jb:je,kb:ke,numq)
      real west(3,nj,nk,numq),newwest(3,nj,nk,numq)
      real east(3,nj,nk,numq),neweast(3,nj,nk,numq)
      real south(ni,3,nk,numq),newsouth(ni,3,nk,numq)
      real north(ni,3,nk,numq),newnorth(ni,3,nk,numq)
      integer reqs(8)
 
      integer i,j,k,n,nn,nr,index
      integer :: index_east,index_west,index_south,index_north
      integer status(MPI_STATUS_SIZE)

!-------------------------------------------------------------------

      index_east = -1
      index_west = -1
      index_south = -1
      index_north = -1

      nr = 0
      if(ibe.eq.0)then
        nr = nr + 1
        index_east = nr
      endif
      if(ibw.eq.0)then
        nr = nr + 1
        index_west = nr
      endif
      if(ibs.eq.0)then
        nr = nr + 1
        index_south = nr
      endif
      if(ibn.eq.0)then
        nr = nr + 1
        index_north = nr
      endif

    nn = 1
    do while( nn .le. nr )
      call MPI_WAITANY(nr,reqs(1),index,status,ierr)
      nn = nn + 1

      if(index.eq.index_east)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,n)
        do k=1,nk
        do n=1,numq
        do j=1,nj
        do i=1,3
          s(ni+i,j,k,n)=neweast(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
      elseif(index.eq.index_west)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,n)
        do k=1,nk
        do n=1,numq
        do j=1,nj
        do i=1,3
          s(i-3,j,k,n)=newwest(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
      elseif(index.eq.index_south)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,n)
        do k=1,nk
        do n=1,numq
        do j=1,3
        do i=1,ni
          s(i,j-3,k,n)=newsouth(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
      elseif(index.eq.index_north)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k,n)
        do k=1,nk
        do n=1,numq
        do j=1,3
        do i=1,ni
          s(i,nj+j,k,n)=newnorth(i,j,k,n)
        enddo
        enddo
        enddo
        enddo
      else
        print *,'  unknown index '
        print *,'  myid,index = ',myid,index
        call stopcm1
      endif

    enddo

    if(timestats.ge.1) time_mps2=time_mps2+mytime()
 
!----------
!  patch for corner
 
      if( (ebc.eq.2.or.wbc.eq.2).and.(sbc.eq.1.or.nbc.eq.1) )then
 
        if(ibw.eq.1)then

          if(p2tchsww)then
!$omp parallel do default(shared)   &
!$omp private(k,n)
            do k=1,nk
            do n=1,numq
              s(0,0,k,n)=s(1,0,k,n)
            enddo
            enddo
          endif

          if(p2tchnww)then
!$omp parallel do default(shared)   &
!$omp private(k,n)
            do k=1,nk
            do n=1,numq
              s(0,nj+1,k,n)=s(1,nj+1,k,n)
            enddo
            enddo
          endif

        endif

        if(ibe.eq.1)then
 
          if(p2tchsee)then
!$omp parallel do default(shared)   &
!$omp private(k,n)
            do k=1,nk
            do n=1,numq
              s(ni+1,0,k,n)=s(ni,0,k,n)
            enddo
            enddo
          endif
 
          if(p2tchnee)then
!$omp parallel do default(shared)   &
!$omp private(k,n)
            do k=1,nk
            do n=1,numq
              s(ni+1,nj+1,k,n)=s(ni,nj+1,k,n)
            enddo
            enddo
          endif
 
        endif
 
      endif

      if( (ebc.eq.1.or.wbc.eq.1).and.(sbc.eq.2.or.nbc.eq.2) )then
 
        if(ibs.eq.1)then
 
          if(p2tchsws)then
!$omp parallel do default(shared)   &
!$omp private(k,n)
            do k=1,nk
            do n=1,numq
              s(0,0,k,n)=s(0,1,k,n)
            enddo
            enddo
          endif
 
          if(p2tchses)then
!$omp parallel do default(shared)   &
!$omp private(k,n)
            do k=1,nk
            do n=1,numq
              s(ni+1,0,k,n)=s(ni+1,1,k,n)
            enddo
            enddo
          endif
 
        endif
 
        if(ibn.eq.1)then
 
          if(p2tchnwn)then
!$omp parallel do default(shared)   &
!$omp private(k,n)
            do k=1,nk
            do n=1,numq
              s(0,nj+1,k,n)=s(0,nj,k,n)
            enddo
            enddo
          endif
 
          if(p2tchnen)then
!$omp parallel do default(shared)   &
!$omp private(k,n)
            do k=1,nk
            do n=1,numq
              s(ni+1,nj+1,k,n)=s(ni+1,nj,k,n)
            enddo
            enddo
          endif
 
        endif
 
      endif
 
      if(timestats.ge.1) time_bc=time_bc+mytime()

!----------

      if(ibw.eq.0)then
        call MPI_WAIT (reqs(5),status,ierr)
      endif

      if(ibe.eq.0)then
        call MPI_WAIT (reqs(6),status,ierr)
      endif

      if(ibn.eq.0)then
        call MPI_WAIT (reqs(7),status,ierr)
      endif

      if(ibs.eq.0)then
        call MPI_WAIT (reqs(8),status,ierr)
      endif

      if(timestats.ge.1) time_mps2=time_mps2+mytime()

!-----------------------------------------------------------
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
#endif


      subroutine prepcorners(s,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      implicit none

      include 'input.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: s
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer :: i,j

#ifdef MPI
      call getcorner(s,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcs2(s)
#endif
      do j=0,nj+1
      do i=0,ni+1
        s(i,j,0)    = s(i,j,1)
        s(i,j,nk+1) = s(i,j,nk)
      enddo
      enddo

      return
      end


