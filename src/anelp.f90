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



      subroutine anelp(xh,uh,xf,uf,yh,vh,yf,vf,                     &
                       zh,mh,rmh,mf,rmf,pi0,thv0,rho0,prs0,rf0,     &
                       radbcw,radbce,radbcs,radbcn,dum1,divx,       &
                       u0,ua,u3d,uten,v0,va,v3d,vten,wa,w3d,wten,   &
                       ppi,pp3d,tha,th3d,thten,thv,cfb,cfa,cfc,     &
                       d1,d2,pdt,deft,rhs,trans,dttmp,nrk,rtime)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie) :: xh,uh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je) :: yh,vh
      real, dimension(jb:je+1) :: yf,vf
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,rho0,prs0,rf0
      real, dimension(jb:je,kb:ke) :: radbcw,radbce
      real, dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, dimension(ib:ie,jb:je,kb:ke) :: dum1,divx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0,ua,u3d,uten
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0,va,v3d,vten
      real, dimension(ib:ie,jb:je,kb:ke+1) :: wa,w3d,wten
      real, dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,tha,th3d,thten,thv
      real, dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: cfb
      real, dimension(kpb:kpe) :: cfa,cfc,d1,d2
      complex, dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: pdt,deft
      complex, dimension(ipb:ipe,jpb:jpe) :: rhs,trans
      real dttmp
      integer nrk
      real, intent(in) :: rtime

!-----

      integer :: i,j,k
      real :: tem
      real*8 :: fluxout,fluxin,u1,v1

!---------------------------------------------------------------------

        if(irbc.eq.2)then
 
          if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,ua)
 
          if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,va)
 
        endif

!---------------------------------------------------------------------

        if(ibw.eq.1.and.wbc.eq.2)then
          do k=1,nk
          do j=1,nj
            uten(1,j,k)=uten(1,j,k)-radbcw(j,k)    &
                      *(ua(2,j,k)-ua(1,j,k))*rdx
          enddo
          enddo

          if(dowr) write(outfile,*) '  This model configuration has not been tested with open boundary conditions '
          call stopcm1
        endif

        if(ibe.eq.1.and.ebc.eq.2)then
          do k=1,nk
          do j=1,nj
            uten(ni+1,j,k)=uten(ni+1,j,k)-radbce(j,k)     &
                         *(ua(ni+1,j,k)-ua(ni  ,j,k))*rdx
          enddo
          enddo

          if(dowr) write(outfile,*) '  This model configuration has not been tested with open boundary conditions '
          call stopcm1
        endif

        if(ibs.eq.1.and.sbc.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            vten(i,1,k)=vten(i,1,k)-radbcs(i,k)    &
                      *(va(i,2,k)-va(i,1,k))*rdy
          enddo
          enddo

          if(dowr) write(outfile,*) '  This model configuration has not been tested with open boundary conditions '
          call stopcm1
        endif

        if(ibn.eq.1.and.nbc.eq.2)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            vten(i,nj+1,k)=vten(i,nj+1,k)-radbcn(i,k)     &
                         *(va(i,nj+1,k)-va(i,nj  ,k))*rdy
          enddo
          enddo

          if(dowr) write(outfile,*) '  This model configuration has not been tested with open boundary conditions '
          call stopcm1
        endif

!---------------------------------------------------------------------
!  convergence forcing:

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. nx.gt.1 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni+1
              u3d(i,j,k)=ua(i,j,k)+dttmp*uten(i,j,k)
            enddo
            enddo
            enddo
            call convinitu(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibw,ibe,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xf,yh,zh,u0,u3d)
            tem=1.0/dttmp
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni+1
              uten(i,j,k) = (u3d(i,j,k)-ua(i,j,k))*tem
            enddo
            enddo
            enddo
          ENDIF
          IF( rtime.le.convtime .and. ny.gt.1 )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj+1
            do i=1,ni
              v3d(i,j,k)=va(i,j,k)+dttmp*vten(i,j,k)
            enddo
            enddo
            enddo
            call convinitv(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xh,yf,zh,v0,v3d)
            tem=1.0/dttmp
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj+1
            do i=1,ni
              vten(i,j,k) = (v3d(i,j,k)-va(i,j,k))*tem
            enddo
            enddo
            enddo
          ENDIF
        ENDIF
        if(timestats.ge.1) time_sound=time_sound+mytime()

!---------------------------------------------------------------------
!  Get pressure

        call poiss(uh,vh,mh,rmh,mf,rmf,pi0,thv0,rho0,rf0,    &
                   dum1,divx,ppi,uten,vten,wten,             &
                   cfb,cfa,cfc,d1,d2,pdt,deft,rhs,trans,dttmp)

!---------------------------------------------------------------------

        tem=dttmp*rdx

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          u3d(i,j,k)=ua(i,j,k)+dttmp*uten(i,j,k)             &
                  -(tem*(ppi(i,j,k)-ppi(i-1,j,k))*uf(i))
        enddo
        enddo
        enddo

        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcu(u3d)

!-----

      IF(axisymm.eq.0)THEN

        tem=dttmp*rdy

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          v3d(i,j,k)=va(i,j,k)+dttmp*vten(i,j,k)             &
                  -(tem*(ppi(i,j,k)-ppi(i,j-1,k))*vf(j))
        enddo
        enddo
        enddo

        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcv(v3d)

      ENDIF

!-----

        tem=dttmp*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=2,nk
        do j=1,nj
        do i=1,ni
          w3d(i,j,k)=wa(i,j,k)+dttmp*wten(i,j,k)             &
                  -(tem*(ppi(i,j,k)-ppi(i,j,k-1))*mf(i,j,k))
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcw(w3d,1)

!-------------------------------------------------------------------- 

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          ppi(i,j,k)=((prs0(1,1,k)+ppi(i,j,k)*rho0(1,1,k))*rp00)**rovcp   &
                    -pi0(1,1,k)
          pp3d(i,j,k)=ppi(i,j,k)
        enddo
        enddo
        do j=1,nj
        do i=1,ni
          th3d(i,j,k)=tha(i,j,k)+dttmp*thten(i,j,k)
        enddo
        enddo
      enddo

        if(timestats.ge.1) time_sound=time_sound+mytime()

        call bcs(ppi)
        call bcs(pp3d)
        call bcs(th3d)

!-------------------------------------------------------------------- 

      return
      end


