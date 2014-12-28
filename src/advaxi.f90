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



      subroutine advsaxi(doweno,bflag,bsq,xh,rxh,arh1,arh2,uh,ruh,xf,vh,rvh,rmh,gz,rgz, &
                  rho0,rr0,rf0,rrf0,advx,dum,mass,rru,s0,s,sten,pdef,dt,weps, &
                  flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      logical, intent(in) :: doweno
      integer bflag
      real*8 bsq
      real, dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: vh,rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh,rho0,rr0,rf0,rrf0
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz
      real, dimension(ib:ie,jb:je,kb:ke) :: advx,dum,mass
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je,kb:ke) :: s0,s,sten
      integer pdef
      real, intent(in) :: dt
      double precision, intent(in) :: weps
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
 
      integer i,j,k,i1,i2,j1,j2
      real :: tem0
      real*8, dimension(nk) :: budx,budy

      real :: dd,rr,phi
      real :: s1,s2,s3,s4,s5
      real :: f1,f2,f3
      real :: b1,b2,b3
      double precision :: bmax
      real :: w1,w2,w3
      double precision :: a1,a2,a3,a4
      logical :: doit

!----------------------------------------------------------------

      i1 = 1
      i2 = ni+1

!-----------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k,s1,s2,s3,s4,s5,f1,f2,f3,b1,b2,b3,w1,w2,w3,a1,a2,a3,a4,bmax,dd,rr,phi,doit)
    DO k=1,nk   ! start of k-loop

      budx(k) = 0.0d0

! Advection in x-direction

    if(doweno)then
      do j=1,nj
      do i=3,i2
        if(rru(i,j,k).ge.0.0)then
          s1=s(i-3,j,k)
          s2=s(i-2,j,k)
          s3=s(i-1,j,k)
          s4=s(i  ,j,k)
          s5=s(i+1,j,k)
        else
          s1=s(i+2,j,k)
          s2=s(i+1,j,k)
          s3=s(i  ,j,k)
          s4=s(i-1,j,k)
          s5=s(i-2,j,k)
        endif

      doit = .true.
      IF( pdef.eq.1 )THEN
        bmax = max(s1,s2,s3,s4,s5)
        if( bmax.lt.Min(1.0d-20,weps) ) doit = .false.
      ENDIF
      IF( i.eq.3 .and. rru(i,j,k).gt.0.0 ) doit = .false.

      IF(doit)THEN

        b1=thdtw*( s1 -2.0*s2 +s3 )**2 + 0.25*(     s1 -4.0*s2 +3.0*s3 )**2
        b2=thdtw*( s2 -2.0*s3 +s4 )**2 + 0.25*(     s2             -s4 )**2
        b3=thdtw*( s3 -2.0*s4 +s5 )**2 + 0.25*( 3.0*s3 -4.0*s4     +s5 )**2

        ! from Jerry Straka (Univ of Oklahoma):
        ! based on Shen and Zha (2010, Int J Num Meth Fluids)
        ! (GHB 120201:  added the "min" part to prevent overflows)
        a1 = 0.10*(1.0+min(1.0d30,abs(b1-b3)/(b1+weps))**2)
        a2 = 0.60*(1.0+min(1.0d30,abs(b1-b3)/(b2+weps))**2)
        a3 = 0.30*(1.0+min(1.0d30,abs(b1-b3)/(b3+weps))**2)

        a4 = 1.0/(a1+a2+a3)
        w1 = a1*a4
        w2 = a2*a4
        w3 = a3*a4

        f1=( f1a*s1 + f1b*s2 + f1c*s3 )
        f2=( f2a*s2 + f2b*s3 + f2c*s4 )
        f3=( f3a*s3 + f3b*s4 + f3c*s5 )

        dum(i,j,k)=rru(i,j,k)*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      ELSE
        dum(i,j,k)=0.0
      ENDIF
      enddo
      enddo
      do j=1,nj
        ! flux at i=3 if u > 0
        i = 3
        if(rru(i,j,k).ge.0.0)then
          dd = s(i-1,j,k)-s(i-2,j,k)
          rr = (s(i,j,k)-s(i-1,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = rru(i,j,k)*( s(i-1,j,k) + 0.5*phi*(s(i-1,j,k)-s(i-2,j,k)) )
        endif
        i = 2
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*0.5*(s(i-1,j,k)+s(i,j,k))
        else
          dd = s(i,j,k)-s(i+1,j,k)
          rr = (s(i-1,j,k)-s(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = rru(i,j,k)*( s(i,j,k) + 0.5*phi*(s(i,j,k)-s(i+1,j,k)) )
        endif
        dum(1,j,k)=0.0
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
      do j=1,nj
        i = ni-1
        if(rru(i,j,k).le.0.0)then
          dd = s(i,j,k)-s(i+1,j,k)
          rr = (s(i-1,j,k)-s(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = rru(i,j,k)*( s(i,j,k) + 0.5*phi*(s(i,j,k)-s(i+1,j,k)) )
        endif
        i = ni
        if(rru(i,j,k).ge.0.0)then
          dd = s(i-1,j,k)-s(i-2,j,k)
          rr = (s(i,j,k)-s(i-1,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = rru(i,j,k)*( s(i-1,j,k) + 0.5*phi*(s(i-1,j,k)-s(i-2,j,k)) )
        else
          dum(i,j,k)=rru(i,j,k)*0.5*(s(i-1,j,k)+s(i,j,k))
        endif
        dum(ni+1,j,k)=0.0
      enddo
      ENDIF
    elseif(hadvordrs.eq.5)then
      do j=1,nj
      do i=4,i2
        if(rru(i,j,k).ge.0.)then
          dum(i,j,k)=rru(i,j,k)*( 2.*s(i-3,j,k)-13.*s(i-2,j,k)   &
                +47.*s(i-1,j,k)+27.*s(i,j,k)-3.*s(i+1,j,k) )*onedsixty
        else
          dum(i,j,k)=rru(i,j,k)*( 2.*s(i+2,j,k)-13.*s(i+1,j,k)   &
                +47.*s(i,j,k)+27.*s(i-1,j,k)-3.*s(i-2,j,k) )*onedsixty
        endif
      enddo
      enddo
      do j=1,nj
        i = 3
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*(-s(i-2,j,k)+5.*s(i-1,j,k)+2.*s(i,j,k))*onedsix
        else
          dum(i,j,k)=rru(i,j,k)*( 2.*s(i+2,j,k)-13.*s(i+1,j,k)   &
                +47.*s(i,j,k)+27.*s(i-1,j,k)-3.*s(i-2,j,k) )*onedsixty
        endif
        i = 2
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*0.5*(s(i-1,j,k)+s(i,j,k))
        else
          dum(i,j,k)=rru(i,j,k)*(-s(i+1,j,k)+5.*s(i,j,k)+2.*s(i-1,j,k))*onedsix
        endif
        dum(1,j,k)=0.0
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
      do j=1,nj
        i = ni-1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*( 2.*s(i-3,j,k)-13.*s(i-2,j,k)   &
                +47.*s(i-1,j,k)+27.*s(i,j,k)-3.*s(i+1,j,k) )*onedsixty
        else
          dum(i,j,k)=rru(i,j,k)*(-s(i+1,j,k)+5.*s(i,j,k)+2.*s(i-1,j,k))*onedsix
        endif
        i = ni
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*(-s(i-2,j,k)+5.*s(i-1,j,k)+2.*s(i,j,k))*onedsix
        else
          dum(i,j,k)=rru(i,j,k)*0.5*(s(i-1,j,k)+s(i,j,k))
        endif
        dum(ni+1,j,k)=0.0
      enddo
      ENDIF
    elseif(hadvordrs.eq.6)then
      do j=1,nj
      do i=4,i2
        dum(i,j,k)=rru(i,j,k)*( 37.0*(s(i  ,j,k)+s(i-1,j,k))     &
                                -8.0*(s(i+1,j,k)+s(i-2,j,k))     &
                                    +(s(i+2,j,k)+s(i-3,j,k)) )*onedsixty
      enddo
      enddo
      do j=1,nj
        i = 3
        dum(i,j,k)=rru(i,j,k)*( 7.0*(s(i  ,j,k)+s(i-1,j,k))          &
                                   -(s(i+1,j,k)+s(i-2,j,k)) )*onedtwelve
        i = 2
        dum(i,j,k)=rru(i,j,k)*0.5*(s(i-1,j,k)+s(i,j,k))
        dum(1,j,k)=0.0
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
      do j=1,nj
        i = ni-1
        dum(i,j,k)=rru(i,j,k)*( 7.0*(s(i  ,j,k)+s(i-1,j,k))          &
                                   -(s(i+1,j,k)+s(i-2,j,k)) )*onedtwelve
        i = ni
        dum(i,j,k)=rru(i,j,k)*0.5*(s(i-1,j,k)+s(i,j,k))
        dum(ni+1,j,k)=0.0
      enddo
      ENDIF
    endif

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
      do j=1,nj
        i=ni+1
        if(rru(i,j,k).le.0.0)then
          dum(i,j,k)=dum(i-1,j,k)*arh1(i-1)/arh2(i-1)
        endif
        budx(k)=budx(k)-dum(ni+1,j,k)*rvh(j)*rmh(ni+1,j,k)
      enddo
      ENDIF

      do j=1,nj
      do i=1,ni
        advx(i,j,k)=-(arh2(i)*dum(i+1,j,k)-arh1(i)*dum(i,j,k))*rdx*uh(i)
      enddo
      enddo

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
        do j=1,nj
          if(rru(ni+1,j,k).le.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-s(i,j,k)*(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)
          endif
        enddo
      ENDIF

    ENDDO   ! end of k-loop

!----------------------------------------------------------------
!  Misc for x-direction

      IF(stat_qsrc.eq.1.and.(wbc.eq.2.or.ebc.eq.2).and.bflag.eq.1)THEN
        tem0=dt*dy*dz
        do k=1,nk
          bsq=bsq+budx(k)*tem0
        enddo
      ENDIF

      IF(pdscheme.eq.1 .and. pdef.eq.1)THEN
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call pdefx(xh,arh1,arh2,uh,rho0,gz,rgz,rru,advx,dum,mass,s0,s,dt,flag,sw31,sw32,se31,se32)
      ENDIF

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advs=time_advs+mytime()
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advuaxi(doweno,arh1,arh2,xf,rxf,arf1,arf2,uf,vh,rho0,rr0,rf0,rrf0,dum,advx,rru,u3d,uten)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      logical, intent(in) :: doweno
      real, dimension(ib:ie) :: arh1,arh2
      real, dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf
      real, dimension(jb:je) :: vh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru,u3d,uten
 
      integer i,j,k,i1,i2,j1,j2,id1,id2
      real :: ubar

      real :: dd,rr,phi
      real :: s1,s2,s3,s4,s5
      real :: f1,f2,f3
      real :: b1,b2,b3
      real :: w1,w2,w3
      double precision :: a1,a2,a3,a4
      double precision :: weps

!------------------------------------------------------------

      weps = 100.0*epsilon

      if(ibw.eq.1)then
        i1=2
      else
        i1=1
      endif
 
      if(ibe.eq.1)then
        i2=ni+1-1
      else
        i2=ni+1
      endif
 
      id1 = i1-1
      id2 = i2

      i1 = max(2,i1)
      id1 = max(1,id1)

!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
    DO k=1,nk

! Advection in x-direction

    if(doweno)then
      do j=1,nj
      do i=2,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          s1=u3d(i-2,j,k)
          s2=u3d(i-1,j,k)
          s3=u3d(i  ,j,k)
          s4=u3d(i+1,j,k)
          s5=u3d(i+2,j,k)
        else
          s1=u3d(i+3,j,k)
          s2=u3d(i+2,j,k)
          s3=u3d(i+1,j,k)
          s4=u3d(i  ,j,k)
          s5=u3d(i-1,j,k)
        endif

        b1=thdtw*( s1 -2.0*s2 +s3 )**2 + 0.25*(     s1 -4.0*s2 +3.0*s3 )**2
        b2=thdtw*( s2 -2.0*s3 +s4 )**2 + 0.25*(     s2             -s4 )**2
        b3=thdtw*( s3 -2.0*s4 +s5 )**2 + 0.25*( 3.0*s3 -4.0*s4     +s5 )**2

        ! from Jerry Straka (Univ of Oklahoma):
        ! based on Shen and Zha (2010, Int J Num Meth Fluids)
        ! (GHB 120201:  added the "min" part to prevent overflows)
        a1 = 0.10*(1.0+min(1.0d30,abs(b1-b3)/(b1+weps))**2)
        a2 = 0.60*(1.0+min(1.0d30,abs(b1-b3)/(b2+weps))**2)
        a3 = 0.30*(1.0+min(1.0d30,abs(b1-b3)/(b3+weps))**2)

        a4 = 1.0/(a1+a2+a3)
        w1 = a1*a4
        w2 = a2*a4
        w3 = a3*a4

        f1=( f1a*s1 + f1b*s2 + f1c*s3 )
        f2=( f2a*s2 + f2b*s3 + f2c*s4 )
        f3=( f3a*s3 + f3b*s4 + f3c*s5 )

        dum(i,j,k)=ubar*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo
      do j=1,nj
        i = 2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dd = u3d(i,j,k)-u3d(i-1,j,k)
          rr = (u3d(i+1,j,k)-u3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = ubar*( u3d(i,j,k) + 0.5*phi*(u3d(i,j,k)-u3d(i-1,j,k)) )
        endif
        i = 1
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k)=ubar*0.5*(arh1(i)*u3d(i,j,k)+arh2(i)*u3d(i+1,j,k))
        else
          dd = u3d(i+1,j,k)-u3d(i+2,j,k)
          rr = (u3d(i,j,k)-u3d(i+1,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = ubar*( u3d(i+1,j,k) + 0.5*phi*(u3d(i+1,j,k)-u3d(i+2,j,k)) )
        endif
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
      do j=1,nj
        i = ni-1
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.le.0.0)then
          dd = u3d(i+1,j,k)-u3d(i+2,j,k)
          rr = (u3d(i,j,k)-u3d(i+1,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = ubar*( u3d(i+1,j,k) + 0.5*phi*(u3d(i+1,j,k)-u3d(i+2,j,k)) )
        endif
        i = ni
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dd = u3d(i,j,k)-u3d(i-1,j,k)
          rr = (u3d(i+1,j,k)-u3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = ubar*( u3d(i,j,k) + 0.5*phi*(u3d(i,j,k)-u3d(i-1,j,k)) )
        else
          dum(i,j,k)=ubar*0.5*(arh1(i)*u3d(i,j,k)+arh2(i)*u3d(i+1,j,k))
        endif
      enddo
      ENDIF
    elseif(hadvordrv.eq.5)then
      do j=1,nj
      do i=3,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.)then
          dum(i,j,k)=ubar*( 2.*u3d(i-2,j,k)-13.*u3d(i-1,j,k)+47.*u3d(i,j,k)   &
                          +27.*u3d(i+1,j,k)-3.*u3d(i+2,j,k) )*onedsixty
        else
          dum(i,j,k)=ubar*( 2.*u3d(i+3,j,k)-13.*u3d(i+2,j,k)+47.*u3d(i+1,j,k)   &
                          +27.*u3d(i,j,k)-3.*u3d(i-1,j,k) )*onedsixty
        endif
      enddo
      enddo
      do j=1,nj
        i = 2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k)=ubar*(-u3d(i-1,j,k)+5.*u3d(i  ,j,k)+2.*u3d(i+1,j,k))*onedsix
        else
          dum(i,j,k)=ubar*( 2.*u3d(i+3,j,k)-13.*u3d(i+2,j,k)+47.*u3d(i+1,j,k)   &
                          +27.*u3d(i,j,k)-3.*u3d(i-1,j,k) )*onedsixty
        endif
        i = 1
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k)=ubar*0.5*(arh1(i)*u3d(i,j,k)+arh2(i)*u3d(i+1,j,k))
        else
          dum(i,j,k)=ubar*(-u3d(i+2,j,k)+5.*u3d(i+1,j,k)+2.*u3d(i  ,j,k))*onedsix
        endif
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
      do j=1,nj
        i = ni-1
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k)=ubar*( 2.*u3d(i-2,j,k)-13.*u3d(i-1,j,k)+47.*u3d(i,j,k)   &
                          +27.*u3d(i+1,j,k)-3.*u3d(i+2,j,k) )*onedsixty
        else
          dum(i,j,k)=ubar*(-u3d(i+2,j,k)+5.*u3d(i+1,j,k)+2.*u3d(i  ,j,k))*onedsix
        endif
        i = ni
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        if(ubar.ge.0.0)then
          dum(i,j,k)=ubar*(-u3d(i-1,j,k)+5.*u3d(i  ,j,k)+2.*u3d(i+1,j,k))*onedsix
        else
          dum(i,j,k)=ubar*0.5*(arh1(i)*u3d(i,j,k)+arh2(i)*u3d(i+1,j,k))
        endif
      enddo
      ENDIF
    elseif(hadvordrv.eq.6)then
      do j=1,nj
      do i=3,id2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k)=ubar*( 37.0*(u3d(i+1,j,k)+u3d(i  ,j,k)) &
                          -8.0*(u3d(i+2,j,k)+u3d(i-1,j,k)) &
                              +(u3d(i+3,j,k)+u3d(i-2,j,k)) )*onedsixty
      enddo
      enddo
      do j=1,nj
        i = 2
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*(  7.0*(u3d(i+1,j,k)+u3d(i  ,j,k)) &
                            -1.0*(u3d(i+2,j,k)+u3d(i-1,j,k)) )*onedtwelve
        i = 1
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*0.5*(u3d(i+1,j,k)+u3d(i,j,k))
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
      do j=1,nj
        i = ni-1
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*(  7.0*(u3d(i+1,j,k)+u3d(i  ,j,k)) &
                            -1.0*(u3d(i+2,j,k)+u3d(i-1,j,k)) )*onedtwelve
        i = ni
        ubar = 0.5*(arh1(i)*rru(i,j,k)+arh2(i)*rru(i+1,j,k))
        dum(i,j,k) = ubar*0.5*(u3d(i+1,j,k)+u3d(i,j,k))
      enddo
      ENDIF
    endif

      do j=1,nj
      advx(1,j,k)=0.0
      do i=2,i2
        advx(i,j,k)=-(arf2(i)*dum(i,j,k)-arf1(i)*dum(i-1,j,k))*rdx*uf(i)
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
        advx(ni+1,j,k)=0.0
      ENDIF
      enddo

    ENDDO

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advu=time_advu+mytime()
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advvaxi(doweno,xh,rxh,arh1,arh2,uh,xf,vf,rho0,rr0,rf0,rrf0,dum,advx,rru,v3d,vten)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      logical, intent(in) :: doweno
      real, dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v3d,vten
 
      integer i,j,k,i1,i2

      real :: dd,rr,phi
      real :: s1,s2,s3,s4,s5
      real :: f1,f2,f3
      real :: b1,b2,b3
      real :: w1,w2,w3
      double precision :: a1,a2,a3,a4
      double precision :: weps
 
!------------------------------------------------------------

      weps = 100.0*epsilon

!-----------------

      i1 = 1
      i2 = ni+1

!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk

! Advection in x-direction

    if(doweno)then
      do j=1,nj
      do i=3,i2
        if(rru(i,j,k).ge.0.0)then
          s1=v3d(i-3,j,k)
          s2=v3d(i-2,j,k)
          s3=v3d(i-1,j,k)
          s4=v3d(i  ,j,k)
          s5=v3d(i+1,j,k)
        else
          s1=v3d(i+2,j,k)
          s2=v3d(i+1,j,k)
          s3=v3d(i  ,j,k)
          s4=v3d(i-1,j,k)
          s5=v3d(i-2,j,k)
        endif

        b1=thdtw*( s1 -2.0*s2 +s3 )**2 + 0.25*(     s1 -4.0*s2 +3.0*s3 )**2
        b2=thdtw*( s2 -2.0*s3 +s4 )**2 + 0.25*(     s2             -s4 )**2
        b3=thdtw*( s3 -2.0*s4 +s5 )**2 + 0.25*( 3.0*s3 -4.0*s4     +s5 )**2

        ! from Jerry Straka (Univ of Oklahoma):
        ! based on Shen and Zha (2010, Int J Num Meth Fluids)
        ! (GHB 120201:  added the "min" part to prevent overflows)
        a1 = 0.10*(1.0+min(1.0d30,abs(b1-b3)/(b1+weps))**2)
        a2 = 0.60*(1.0+min(1.0d30,abs(b1-b3)/(b2+weps))**2)
        a3 = 0.30*(1.0+min(1.0d30,abs(b1-b3)/(b3+weps))**2)

        a4 = 1.0/(a1+a2+a3)
        w1 = a1*a4
        w2 = a2*a4
        w3 = a3*a4

        f1=( f1a*s1 + f1b*s2 + f1c*s3 )
        f2=( f2a*s2 + f2b*s3 + f2c*s4 )
        f3=( f3a*s3 + f3b*s4 + f3c*s5 )

        dum(i,j,k)=rru(i,j,k)*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo
      do j=1,nj
        i = 3
        if(rru(i,j,k).ge.0.0)then
          dd = v3d(i-1,j,k)-v3d(i-2,j,k)
          rr = (v3d(i,j,k)-v3d(i-1,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = rru(i,j,k)*( v3d(i-1,j,k) + 0.5*phi*(v3d(i-1,j,k)-v3d(i-2,j,k)) )
        endif
        i = 2
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*0.5*(v3d(i-1,j,k)+v3d(i,j,k))
        else
          dd = v3d(i,j,k)-v3d(i+1,j,k)
          rr = (v3d(i-1,j,k)-v3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = rru(i,j,k)*( v3d(i,j,k) + 0.5*phi*(v3d(i,j,k)-v3d(i+1,j,k)) )
        endif
        dum(1,j,k)=0.0
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
      do j=1,nj
        i = ni-1
        if(rru(i,j,k).le.0.0)then
          dd = v3d(i,j,k)-v3d(i+1,j,k)
          rr = (v3d(i-1,j,k)-v3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = rru(i,j,k)*( v3d(i,j,k) + 0.5*phi*(v3d(i,j,k)-v3d(i+1,j,k)) )
        endif
        i = ni
        if(rru(i,j,k).ge.0.0)then
          dd = v3d(i-1,j,k)-v3d(i-2,j,k)
          rr = (v3d(i,j,k)-v3d(i-1,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = rru(i,j,k)*( v3d(i-1,j,k) + 0.5*phi*(v3d(i-1,j,k)-v3d(i-2,j,k)) )
        else
          dum(i,j,k)=rru(i,j,k)*0.5*(v3d(i-1,j,k)+v3d(i,j,k))
        endif
        dum(ni+1,j,k)=0.0
      enddo
      ENDIF
    elseif(hadvordrv.eq.5)then
      do j=1,nj
      do i=4,i2
        if(rru(i,j,k).ge.0.)then
          dum(i,j,k)=rru(i,j,k)*( 2.*v3d(i-3,j,k)-13.*v3d(i-2,j,k)+47.*v3d(i-1,j,k)    &
                                +27.*v3d(i,j,k)-3.*v3d(i+1,j,k) )*onedsixty
        else
          dum(i,j,k)=rru(i,j,k)*( 2.*v3d(i+2,j,k)-13.*v3d(i+1,j,k)+47.*v3d(i,j,k)    &
                                +27.*v3d(i-1,j,k)-3.*v3d(i-2,j,k) )*onedsixty
        endif
      enddo
      enddo
      do j=1,nj
        i = 3
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*(-v3d(i-2,j,k)+5.*v3d(i-1,j,k)+2.*v3d(i,j,k))*onedsix
        else
          dum(i,j,k)=rru(i,j,k)*( 2.*v3d(i+2,j,k)-13.*v3d(i+1,j,k)+47.*v3d(i,j,k)    &
                                +27.*v3d(i-1,j,k)-3.*v3d(i-2,j,k) )*onedsixty
        endif
        i = 2
        if(rru(i,j,k).ge.0.)then
          dum(i,j,k)=rru(i,j,k)*0.5*(v3d(i-1,j,k)+v3d(i,j,k))
        else
          dum(i,j,k)=rru(i,j,k)*(-v3d(i+1,j,k)+5.*v3d(i,j,k)+2.*v3d(i-1,j,k))*onedsix
        endif
        dum(1,j,k)=0.0
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
      do j=1,nj
        i = ni-1
        if(rru(i,j,k).ge.0.)then
          dum(i,j,k)=rru(i,j,k)*( 2.*v3d(i-3,j,k)-13.*v3d(i-2,j,k)+47.*v3d(i-1,j,k)    &
                                +27.*v3d(i,j,k)-3.*v3d(i+1,j,k) )*onedsixty
        else
          dum(i,j,k)=rru(i,j,k)*(-v3d(i+1,j,k)+5.*v3d(i,j,k)+2.*v3d(i-1,j,k))*onedsix
        endif
        i = ni
        if(rru(i,j,k).ge.0.)then
          dum(i,j,k)=rru(i,j,k)*(-v3d(i-2,j,k)+5.*v3d(i-1,j,k)+2.*v3d(i,j,k))*onedsix
        else
          dum(i,j,k)=rru(i,j,k)*0.5*(v3d(i-1,j,k)+v3d(i,j,k))
        endif
        dum(ni+1,j,k)=0.0
      enddo
      ENDIF
    elseif(hadvordrv.eq.6)then
      do j=1,nj
      do i=4,i2
        dum(i,j,k)=rru(i,j,k)*( 37.0*(v3d(i  ,j,k)+v3d(i-1,j,k)) &
                                -8.0*(v3d(i+1,j,k)+v3d(i-2,j,k)) &
                                    +(v3d(i+2,j,k)+v3d(i-3,j,k)) )*onedsixty
      enddo
      enddo
      do j=1,nj
        i = 3
        dum(i,j,k)=rru(i,j,k)*( 7.0*(v3d(i  ,j,k)+v3d(i-1,j,k))          &
                                   -(v3d(i+1,j,k)+v3d(i-2,j,k)) )*onedtwelve
        i = 2
        dum(i,j,k)=rru(i,j,k)*0.5*(v3d(i-1,j,k)+v3d(i,j,k))
        dum(1,j,k)=0.0
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
      do j=1,nj
        i = ni-1
        dum(i,j,k)=rru(i,j,k)*( 7.0*(v3d(i  ,j,k)+v3d(i-1,j,k))          &
                                   -(v3d(i+1,j,k)+v3d(i-2,j,k)) )*onedtwelve
        i = ni
        dum(i,j,k)=rru(i,j,k)*0.5*(v3d(i-1,j,k)+v3d(i,j,k))
        dum(ni+1,j,k)=0.0
      enddo
      ENDIF
    endif

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
      do j=1,nj
        i=ni+1
        if(rru(i,j,k).le.0.0)then
          i = ni
          dum(ni+1,j,k)=arh1(i)*arh1(i)*dum(i,j,k)/(arh2(i)*arh2(i))
        endif
      enddo
      ENDIF

      do j=1,nj
      do i=1,ni
        ! cm1r17: include centrifugal accel term here
        advx(i,j,k)=-(arh2(i)*arh2(i)*dum(i+1,j,k)-arh1(i)*arh1(i)*dum(i,j,k))*rdx*uh(i)
      enddo
      enddo

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
        do j=1,nj
          if(rru(ni+1,j,k).le.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-v3d(i,j,k)*(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)
          endif
        enddo
      ENDIF

    ENDDO

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advv=time_advv+mytime()
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advwaxi(doweno,xh,rxh,arh1,arh2,uh,xf,vh,rho0,rr0,rf0,rrf0,dum,advx,rru,w3d,wten,c1,c2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      logical, intent(in) :: doweno
      real, dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: vh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w3d,wten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
 
      integer i,j,k,i1,i2,j1,j2
      real :: ubar

      real :: dd,rr,phi
      real :: s1,s2,s3,s4,s5
      real :: f1,f2,f3
      real :: b1,b2,b3
      real :: w1,w2,w3
      double precision :: a1,a2,a3,a4
      double precision :: weps

!----------------------------------------------------------------

      weps = 100.0*epsilon

      i1 = 1
      i2 = ni+1

!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar)
    DO k=2,nk

! Advection in x-direction

    if(doweno)then
      do j=1,nj
      do i=3,i2
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          s1=w3d(i-3,j,k)
          s2=w3d(i-2,j,k)
          s3=w3d(i-1,j,k)
          s4=w3d(i  ,j,k)
          s5=w3d(i+1,j,k)
        else
          s1=w3d(i+2,j,k)
          s2=w3d(i+1,j,k)
          s3=w3d(i  ,j,k)
          s4=w3d(i-1,j,k)
          s5=w3d(i-2,j,k)
        endif

        b1=thdtw*( s1 -2.0*s2 +s3 )**2 + 0.25*(     s1 -4.0*s2 +3.0*s3 )**2
        b2=thdtw*( s2 -2.0*s3 +s4 )**2 + 0.25*(     s2             -s4 )**2
        b3=thdtw*( s3 -2.0*s4 +s5 )**2 + 0.25*( 3.0*s3 -4.0*s4     +s5 )**2

        ! from Jerry Straka (Univ of Oklahoma):
        ! based on Shen and Zha (2010, Int J Num Meth Fluids)
        ! (GHB 120201:  added the "min" part to prevent overflows)
        a1 = 0.10*(1.0+min(1.0d30,abs(b1-b3)/(b1+weps))**2)
        a2 = 0.60*(1.0+min(1.0d30,abs(b1-b3)/(b2+weps))**2)
        a3 = 0.30*(1.0+min(1.0d30,abs(b1-b3)/(b3+weps))**2)

        a4 = 1.0/(a1+a2+a3)
        w1 = a1*a4
        w2 = a2*a4
        w3 = a3*a4

        f1=( f1a*s1 + f1b*s2 + f1c*s3 )
        f2=( f2a*s2 + f2b*s3 + f2c*s4 )
        f3=( f3a*s3 + f3b*s4 + f3c*s5 )

        dum(i,j,k)=ubar*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo
      do j=1,nj
        i = 3
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dd = w3d(i-1,j,k)-w3d(i-2,j,k)
          rr = (w3d(i,j,k)-w3d(i-1,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = ubar*( w3d(i-1,j,k) + 0.5*phi*(w3d(i-1,j,k)-w3d(i-2,j,k)) )
        endif
        i = 2
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k)=ubar*0.5*(w3d(i-1,j,k)+w3d(i,j,k))
        else
          dd = w3d(i,j,k)-w3d(i+1,j,k)
          rr = (w3d(i-1,j,k)-w3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = ubar*( w3d(i,j,k) + 0.5*phi*(w3d(i,j,k)-w3d(i+1,j,k)) )
        endif
        dum(1,j,k)=0.0
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
      do j=1,nj
        i = ni-1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.le.0.0)then
          dd = w3d(i,j,k)-w3d(i+1,j,k)
          rr = (w3d(i-1,j,k)-w3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = ubar*( w3d(i,j,k) + 0.5*phi*(w3d(i,j,k)-w3d(i+1,j,k)) )
        endif
        i = ni
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dd = w3d(i-1,j,k)-w3d(i-2,j,k)
          rr = (w3d(i,j,k)-w3d(i-1,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = ubar*( w3d(i-1,j,k) + 0.5*phi*(w3d(i-1,j,k)-w3d(i-2,j,k)) )
        else
          dum(i,j,k)=ubar*0.5*(w3d(i-1,j,k)+w3d(i,j,k))
        endif
        dum(ni+1,j,k)=0.0
      enddo
      ENDIF
    elseif(hadvordrv.eq.5)then
      do j=1,nj
      do i=4,i2
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.)then
          dum(i,j,k)=ubar*                                              &
                 ( 2.*w3d(i-3,j,k)-13.*w3d(i-2,j,k)+47.*w3d(i-1,j,k)    &
                  +27.*w3d(i,j,k)-3.*w3d(i+1,j,k) )*onedsixty
        else
          dum(i,j,k)=ubar*                                            &
                 ( 2.*w3d(i+2,j,k)-13.*w3d(i+1,j,k)+47.*w3d(i,j,k)    &
                  +27.*w3d(i-1,j,k)-3.*w3d(i-2,j,k) )*onedsixty
        endif
      enddo
      enddo
      do j=1,nj
        i = 3
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k)=ubar*(-w3d(i-2,j,k)+5.*w3d(i-1,j,k)+2.*w3d(i  ,j,k))*onedsix
        else
          dum(i,j,k)=ubar*                                            &
                 ( 2.*w3d(i+2,j,k)-13.*w3d(i+1,j,k)+47.*w3d(i,j,k)    &
                  +27.*w3d(i-1,j,k)-3.*w3d(i-2,j,k) )*onedsixty
        endif
        i = 2
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.)then
          dum(i,j,k)=ubar*0.5*(w3d(i-1,j,k)+w3d(i,j,k))
        else
          dum(i,j,k)=ubar*(-w3d(i+1,j,k)+5.*w3d(i  ,j,k)+2.*w3d(i-1,j,k))*onedsix
        endif
        dum(1,j,k)=0.0
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
      do j=1,nj
        i = ni-1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.)then
          dum(i,j,k)=ubar*                                              &
                 ( 2.*w3d(i-3,j,k)-13.*w3d(i-2,j,k)+47.*w3d(i-1,j,k)    &
                  +27.*w3d(i,j,k)-3.*w3d(i+1,j,k) )*onedsixty
        else
          dum(i,j,k)=ubar*(-w3d(i+1,j,k)+5.*w3d(i  ,j,k)+2.*w3d(i-1,j,k))*onedsix
        endif
        i = ni
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.ge.0.)then
          dum(i,j,k)=ubar*(-w3d(i-2,j,k)+5.*w3d(i-1,j,k)+2.*w3d(i  ,j,k))*onedsix
        else
          dum(i,j,k)=ubar*0.5*(w3d(i-1,j,k)+w3d(i,j,k))
        endif
        dum(ni+1,j,k)=0.0
      enddo
      ENDIF
    elseif(hadvordrv.eq.6)then
      do j=1,nj
      do i=4,i2
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k)=ubar*                                  &
                   ( 37.0*(w3d(i  ,j,k)+w3d(i-1,j,k))     &
                     -8.0*(w3d(i+1,j,k)+w3d(i-2,j,k))     &
                         +(w3d(i+2,j,k)+w3d(i-3,j,k)) )*onedsixty
      enddo
      enddo
      do j=1,nj
        i = 3
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k)=ubar*( 7.0*(w3d(i  ,j,k)+w3d(i-1,j,k))          &
                             -(w3d(i+1,j,k)+w3d(i-2,j,k)) )*onedtwelve
        i = 2
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k)=ubar*0.5*(w3d(i-1,j,k)+w3d(i,j,k))
        dum(1,j,k)=0.0
      enddo
      IF(ebc.eq.3.or.ebc.eq.4)THEN
      do j=1,nj
        i = ni-1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k)=ubar*( 7.0*(w3d(i  ,j,k)+w3d(i-1,j,k))          &
                             -(w3d(i+1,j,k)+w3d(i-2,j,k)) )*onedtwelve
        i = ni
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        dum(i,j,k)=ubar*0.5*(w3d(i-1,j,k)+w3d(i,j,k))
        dum(ni+1,j,k)=0.0
      enddo
      ENDIF
    endif

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
      do j=1,nj
        i=ni+1
        ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
        if(ubar.le.0.0)then
          i=ni
          dum(i+1,j,k)=dum(i,j,k)*arh1(i)/arh2(i)
        endif
      enddo
      ENDIF

      do j=1,nj
      do i=1,ni
        advx(i,j,k)=-(arh2(i)*dum(i+1,j,k)-arh1(i)*dum(i,j,k))*rdx*uh(i)
      enddo
      enddo

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
        do j=1,nj
          i=ni+1
          ubar = c2(1,1,k)*rru(i,j,k)+c1(1,1,k)*rru(i,j,k-1)
          if(ubar.le.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-w3d(i,j,k)*(                                  &
                    c1(1,1,k)*(arh2(i)*rru(i+1,j,k-1)-arh1(i)*rru(i,j,k-1))         &
                   +c2(1,1,k)*(arh2(i)*rru(i+1,j,k  )-arh1(i)*rru(i,j,k  )) )*rdx*uh(i)
          endif
        enddo
      ENDIF

    ENDDO

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advw=time_advw+mytime()
 
      return
      end


