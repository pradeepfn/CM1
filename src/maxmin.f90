/* Copyright (C) 1991-2014 Free Software Foundation, Inc.
   This file is part of the GNU C Library.

   The GNU C Library is free software; you can redistribute it and/or
   modify it under the terms of the GNU Lesser General Public
   License as published by the Free Software Foundation; either
   version 2.1 of the License, or (at your option) any later version.

   The GNU C Library is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
   Lesser General Public License for more details.

   You should have received a copy of the GNU Lesser General Public
   License along with the GNU C Library; if not, see
   <http://www.gnu.org/licenses/>.  */


/* This header is separate from features.h so that the compiler can
   include it implicitly at the start of every compilation.  It must
   not itself include <features.h> or any other header that includes
   <features.h> because the implicit include comes before any feature
   test macros that may be defined in a source file before it first
   explicitly includes a system header.  GCC knows the name of this
   header in order to preinclude it.  */

/* glibc's intent is to support the IEC 559 math functionality, real
   and complex.  If the GCC (4.9 and later) predefined macros
   specifying compiler intent are available, use them to determine
   whether the overall intent is to support these features; otherwise,
   presume an older compiler has intent to support these features and
   define these macros by default.  */



/* wchar_t uses ISO/IEC 10646 (2nd ed., published 2011-03-15) /
   Unicode 6.0.  */

/* We do not support C11 <threads.h>.  */



      subroutine maxmin(izz,jzz,kzz,f,nstat,rstat,amax,amin)
      implicit none

      include 'input.incl'
      include 'timestat.incl'
      include 'mpif.h'

      integer :: izz,jzz,kzz,nstat
      real, dimension(stat_out) :: rstat
      real, dimension(1-ngxy:izz+ngxy,1-ngxy:jzz+ngxy,1-ngz:kzz+ngz) :: f
      character*6 :: amax,amin

!-----------------------------------------------------------------------

      integer :: i,j,k
      integer :: imax,jmax,kmax,imin,jmin,kmin
      integer, dimension(nk+1) :: imaxt,jmaxt,kmaxt,imint,jmint,kmint
      real, dimension(nk+1) :: tmax,tmin
      real :: fmax,fmin,rmax,rmin
      integer :: loc
      real, dimension(2) :: mmax,nmax,mmin,nmin

!-----------------------------------------------------------------------

      imin = 1
      jmin = 1
      kmin = 1
      imax = 1
      jmax = 1
      kmax = 1

!$omp parallel do default(shared)    &
!$omp private(i,j,k)
      do k=1,kzz
        tmax(k)= -1.e30
        tmin(k)=  1.e30
        do j=1,jzz
        do i=1,izz
          if(f(i,j,k).gt.tmax(k))then
            tmax(k)=f(i,j,k)
            imaxt(k)=i
            jmaxt(k)=j
            kmaxt(k)=k
          endif
          if(f(i,j,k).lt.tmin(k))then
            tmin(k)=f(i,j,k)
            imint(k)=i
            jmint(k)=j
            kmint(k)=k
          endif
        enddo
        enddo
      enddo

      fmax= -1.e30
      fmin=  1.e30
      do k=1,kzz
        if(tmax(k).gt.fmax)then
          fmax=tmax(k)
          imax=imaxt(k)
          jmax=jmaxt(k)
          kmax=kmaxt(k)
        endif
        if(tmin(k).lt.fmin)then
          fmin=tmin(k)
          imin=imint(k)
          jmin=jmint(k)
          kmin=kmint(k)
        endif 
      enddo

      mmax(1)=fmax
      mmax(2)=myid
      call MPI_ALLREDUCE(mmax,nmax,1,MPI_2REAL,MPI_MAXLOC,   &
                         MPI_COMM_WORLD,ierr)
      loc=nint(nmax(2))
      imax=imax+(myi-1)*ni
      jmax=jmax+(myj-1)*nj
      call MPI_BCAST(imax,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jmax,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kmax,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)

      mmin(1)=fmin
      mmin(2)=myid
      call MPI_ALLREDUCE(mmin,nmin,1,MPI_2REAL,MPI_MINLOC,   &
                         MPI_COMM_WORLD,ierr)
      loc=nint(nmin(2))
      imin=imin+(myi-1)*ni
      jmin=jmin+(myj-1)*nj
      call MPI_BCAST(imin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jmin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kmin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)

      fmax=nmax(1)
      fmin=nmin(1)

    if(myid.eq.0)then
      write(6,100) amax,fmax,imax,jmax,kmax,    &
                   amin,fmin,imin,jmin,kmin
100   format(2x,a6,':',1x,g13.6,i5,i5,i5,    &
             4x,a6,':',1x,g13.6,i5,i5,i5)

      nstat = nstat + 1
      rstat(nstat) = fmax
      nstat = nstat + 1
      rstat(nstat) = fmin
    endif

      if(timestats.ge.1) time_stat=time_stat+mytime()
 
      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine maxmin2d(izz,jzz,f,nstat,rstat,amax,amin)
      implicit none
        
      include 'input.incl'
      include 'timestat.incl'
      include 'mpif.h'

      integer :: izz,jzz,nstat
      real, dimension(stat_out) :: rstat
      real, dimension(1-ngxy:izz+ngxy,1-ngxy:jzz+ngxy) :: f
      character*6 :: amax,amin
        
!-----------------------------------------------------------------------
          
      integer :: i,j
      integer :: imax,jmax,imin,jmin
      integer, dimension(jzz) :: imaxt,jmaxt,imint,jmint
      real, dimension(jzz) :: tmax,tmin
      real :: fmax,fmin,rmax,rmin
      integer :: loc
      real, dimension(2) :: mmax,nmax,mmin,nmin
          
!-----------------------------------------------------------------------

!$omp parallel do default(shared)    &
!$omp private(i,j)
      do j=1,jzz
        tmax(j)= -1.e30
        tmin(j)=  1.e30
        do i=1,izz
          if(f(i,j).gt.tmax(j))then
            tmax(j)=f(i,j)
            imaxt(j)=i
            jmaxt(j)=j
          endif
          if(f(i,j).lt.tmin(j))then
            tmin(j)=f(i,j)
            imint(j)=i
            jmint(j)=j
          endif
        enddo
      enddo

      fmax= -1.e30
      fmin=  1.e30
      do j=1,jzz
        if(tmax(j).gt.fmax)then
          fmax=tmax(j)
          imax=imaxt(j)
          jmax=jmaxt(j)
        endif
        if(tmin(j).lt.fmin)then
          fmin=tmin(j)
          imin=imint(j)
          jmin=jmint(j)
        endif
      enddo

      mmax(1)=fmax
      mmax(2)=myid
      call MPI_ALLREDUCE(mmax,nmax,1,MPI_2REAL,MPI_MAXLOC,   &
                         MPI_COMM_WORLD,ierr)
      loc=nint(nmax(2))
      imax=imax+(myi-1)*ni
      jmax=jmax+(myj-1)*nj
      call MPI_BCAST(imax,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jmax,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)

      mmin(1)=fmin
      mmin(2)=myid
      call MPI_ALLREDUCE(mmin,nmin,1,MPI_2REAL,MPI_MINLOC,   &
                         MPI_COMM_WORLD,ierr)
      loc=nint(nmin(2))
      imin=imin+(myi-1)*ni
      jmin=jmin+(myj-1)*nj
      call MPI_BCAST(imin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jmin,1,MPI_INTEGER,loc,MPI_COMM_WORLD,ierr)

      fmax=nmax(1)
      fmin=nmin(1)

    if(myid.eq.0)then
      write(6,100) amax,fmax,imax,jmax,1,    &
                   amin,fmin,imin,jmin,1
100   format(2x,a6,':',1x,g13.6,i5,i5,i5,    &
             4x,a6,':',1x,g13.6,i5,i5,i5)

      nstat = nstat + 1
      rstat(nstat) = fmax
      nstat = nstat + 1
      rstat(nstat) = fmin
    endif

      if(timestats.ge.1) time_stat=time_stat+mytime()

      return
      end


