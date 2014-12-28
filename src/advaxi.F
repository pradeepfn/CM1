

      subroutine advsaxi(doweno,bflag,bsq,xh,rxh,uh,ruh,xf,vh,rvh,rmh,        &
                  rho0,rr0,rf0,rrf0,advx,dum,mass,rru,s0,s,sten,pdef,dt,weps, &
                  flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      logical, intent(in) :: doweno
      integer bflag
      real*8 bsq
      real, dimension(ib:ie) :: xh,rxh,uh,ruh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: vh,rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh,rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: advx,dum,mass
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je,kb:ke) :: s0,s,sten
      integer pdef
      real, intent(in) :: dt,weps
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
 
      integer i,j,k,i1,i2,j1,j2
      real, parameter :: tem = 1.0/60.0
      real, parameter :: tem2 = 1.0/6.0
      logical :: doitw,doite,doits,doitn
      real :: tem0
      real*8, dimension(nk) :: budx,budy

      real s1,s2,s3,s4,s5
      real f1,f2,f3
      real b1,b2,b3,bmax
      real w1,w2,w3
      real a1,a2,a3,a4
      logical :: doit

!----------------------------------------------------------------

      i1 = 1
      i2 = ni+1

      if(wbc.eq.2 .and. ibw.eq.1)then
        i1=4
      endif

      if(ebc.eq.2 .and. ibe.eq.1)then
        i2=ni-2
      endif

      doitw = .false.
      doite = .false.

      IF(wbc.eq.2 .and. ibw.eq.1) doitw = .true.
      IF(ebc.eq.2 .and. ibe.eq.1) doite = .true.

!-----------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k,s1,s2,s3,s4,s5,f1,f2,f3,b1,b2,b3,w1,w2,w3,a1,a2,a3,a4,bmax,doit)
    DO k=1,nk   ! start of k-loop

      budx(k) = 0.0d0

! Advection in x-direction

    if(doweno)then
      do j=1,nj
      do i=i1,i2
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
        if( bmax.lt.1.0e-20 ) doit = .false.
      ENDIF

      IF(doit)THEN

        b1=thdtw*( s1 -2.0*s2 +s3 )**2 + 0.25*(     s1 -4.0*s2 +3.0*s3 )**2
        b2=thdtw*( s2 -2.0*s3 +s4 )**2 + 0.25*(     s2             -s4 )**2
        b3=thdtw*( s3 -2.0*s4 +s5 )**2 + 0.25*( 3.0*s3 -4.0*s4     +s5 )**2

        ! from Jerry Straka (Univ of Oklahoma):
        ! based on Shen and Zha (2010, Int J Num Meth Fluids)
        ! (GHB 120201:  added the "min" part to prevent overflows)
        a1 = 0.10*(1.0+min(1.0e15,abs(b1-b3)/(b1+weps))**2)
        a2 = 0.60*(1.0+min(1.0e15,abs(b1-b3)/(b2+weps))**2)
        a3 = 0.30*(1.0+min(1.0e15,abs(b1-b3)/(b3+weps))**2)

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
    elseif(hadvorder.eq.5)then
      do j=1,nj
      do i=i1,i2
#ifdef SGN
        ! this form of code is used for IBM/xlf compiler:
        dum(i,j,k)=rru(i,j,k)*( 37.0*(s(i  ,j,k)+s(i-1,j,k))     &
                                -8.0*(s(i+1,j,k)+s(i-2,j,k))     &
                                    +(s(i+2,j,k)+s(i-3,j,k))     &
        -sign(1.0,rru(i,j,k))*( 10.0*(s(i  ,j,k)-s(i-1,j,k))     &
                                -5.0*(s(i+1,j,k)-s(i-2,j,k))     &
                                    +(s(i+2,j,k)-s(i-3,j,k)) ) )*tem
#else
        ! this seems to be faster on most other platforms:
        if(rru(i,j,k).ge.0.)then
          dum(i,j,k)=rru(i,j,k)*( 2.*s(i-3,j,k)-13.*s(i-2,j,k)   &
                +47.*s(i-1,j,k)+27.*s(i,j,k)-3.*s(i+1,j,k) )*tem
        else
          dum(i,j,k)=rru(i,j,k)*( 2.*s(i+2,j,k)-13.*s(i+1,j,k)   &
                +47.*s(i,j,k)+27.*s(i-1,j,k)-3.*s(i-2,j,k) )*tem
        endif
#endif
      enddo
      enddo
    elseif(hadvorder.eq.6)then
      do j=1,nj
      do i=i1,i2
        dum(i,j,k)=rru(i,j,k)*( 37.0*(s(i  ,j,k)+s(i-1,j,k))     &
                                -8.0*(s(i+1,j,k)+s(i-2,j,k))     &
                                    +(s(i+2,j,k)+s(i-3,j,k)) )*tem
      enddo
      enddo
    endif

      IF(doitw)THEN
      do j=1,nj
        i=3
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*(-s(i-2,j,k)+5.*s(i-1,j,k)+2.*s(i,j,k))*tem2
        else
          dum(i,j,k)=rru(i,j,k)*(-s(i+1,j,k)+5.*s(i,j,k)+2.*s(i-1,j,k))*tem2
        endif
        i=2
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*s(i-1,j,k)
        else
          dum(i,j,k)=rru(i,j,k)*s(i,j,k)
        endif
        i=1
        if(rru(i,j,k).lt.0.0)then
          dum(i,j,k)=rru(i,j,k)*s(i,j,k)
        else
          dum(i,j,k)=dum(i+1,j,k)
        endif
        budx(k)=budx(k)+dum(1,j,k)*rvh(j)*rmh(1,j,k)
      enddo
      ENDIF

      IF(doite)THEN
      do j=1,nj
        i=ni-1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*(-s(i-2,j,k)+5.*s(i-1,j,k)+2.*s(i,j,k))*tem2
        else
          dum(i,j,k)=rru(i,j,k)*(-s(i+1,j,k)+5.*s(i,j,k)+2.*s(i-1,j,k))*tem2
        endif
        i=ni
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*s(i-1,j,k)
        else
          dum(i,j,k)=rru(i,j,k)*s(i,j,k)
        endif
        i=ni+1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=rru(i,j,k)*s(i-1,j,k)
        else
          dum(i,j,k)=dum(i-1,j,k)*xf(ni)/xf(ni+1)
        endif
        budx(k)=budx(k)-dum(ni+1,j,k)*rvh(j)*rmh(ni+1,j,k)
      enddo
      ENDIF

      do j=1,nj
      do i=1,ni
        advx(i,j,k)=-(xf(i+1)*dum(i+1,j,k)-xf(i)*dum(i,j,k))*rdx*uh(i)*rxh(i)
      enddo
      enddo

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
        call pdefx(xh,rho0,advx,dum,mass,s0,dt,flag,sw31,sw32,se31,se32)
      ENDIF

      IF(doitw)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
          i=1
          if(rru(i,j,k).ge.0.0)then
            advx(i,j,k)=advx(i,j,k)-s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)
          endif
        enddo
        enddo
      ENDIF

      IF(doite)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
          i=ni+1
          if(rru(i,j,k).lt.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-s(i,j,k)*(xf(i+1)*rru(i+1,j,k)-xf(i)*rru(i,j,k))*rdx*uh(i)*rxh(i)
          endif
        enddo
        enddo
      ENDIF

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advs=time_advs+mytime()
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advuaxi(xf,rxf,uf,vh,rho0,rr0,rf0,rrf0,dum,advx,rru,u3d,uten)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(ib:ie+1) :: xf,rxf,uf
      real, dimension(jb:je) :: vh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru,u3d,uten
 
      integer i,j,k,i1,i2,j1,j2,id1,id2
      real, parameter :: tem = 1.0/120.0
      real, parameter :: tem2 = 1.0/12.0

!------------------------------------------------------------

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
 
!-----------------

      id1 = i1-1
      id2 = i2

      if(wbc.eq.2 .and. ibw.eq.1)then
        id1 = 3
      endif

      if(ebc.eq.2 .and. ibe.eq.1)then
        id2 = ni-2
      endif

!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk

! Advection in x-direction

    if(hadvorder.eq.5)then
      do j=1,nj
      do i=id1,id2
#ifdef SGN
        ! this form of code is used for IBM/xlf compiler:
        dum(i,j,k)=(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))           &
                             *( 37.0*(u3d(i+1,j,k)+u3d(i  ,j,k))     &
                                -8.0*(u3d(i+2,j,k)+u3d(i-1,j,k))     &
                                    +(u3d(i+3,j,k)+u3d(i-2,j,k))     &
        -sign(1.0,(rru(i,j,k)+rru(i+1,j,k)))                         &
                             *( 10.0*(u3d(i+1,j,k)-u3d(i  ,j,k))     &
                                -5.0*(u3d(i+2,j,k)-u3d(i-1,j,k))     &
                                    +(u3d(i+3,j,k)-u3d(i-2,j,k)) ) )*tem
#else
        ! this seems to be faster on most other platforms:
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*                     &
                 ( 2.*u3d(i-2,j,k)-13.*u3d(i-1,j,k)+47.*u3d(i,j,k)   &
                  +27.*u3d(i+1,j,k)-3.*u3d(i+2,j,k) )*tem
        else
          dum(i,j,k)=(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*                       &
                 ( 2.*u3d(i+3,j,k)-13.*u3d(i+2,j,k)+47.*u3d(i+1,j,k)   &
                  +27.*u3d(i,j,k)-3.*u3d(i-1,j,k) )*tem
        endif
#endif
      enddo
      enddo
    elseif(hadvorder.eq.6)then
      do j=1,nj
      do i=id1,id2
        dum(i,j,k)=(rru(i,j,k)+rru(i+1,j,k))*             &
                   ( 37.0*(u3d(i+1,j,k)+u3d(i  ,j,k))     &
                     -8.0*(u3d(i+2,j,k)+u3d(i-1,j,k))     &
                         +(u3d(i+3,j,k)+u3d(i-2,j,k)) )*tem
      enddo
      enddo
    endif

      IF(wbc.eq.2 .and. ibw.eq.1)THEN
      do j=1,nj
        i=2
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i+1,j,k))*                     &
                 (-u3d(i-1,j,k)+5.*u3d(i  ,j,k)+2.*u3d(i+1,j,k))*tem2
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i+1,j,k))*                       &
                 (-u3d(i+2,j,k)+5.*u3d(i+1,j,k)+2.*u3d(i  ,j,k))*tem2
        endif
        i=1
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i+1,j,k))*u3d(i  ,j,k)
        else
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i+1,j,k))*u3d(i+1,j,k)
        endif
      enddo
      ENDIF

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
      do j=1,nj
        i=ni-1
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*                     &
                 (-u3d(i-1,j,k)+5.*u3d(i  ,j,k)+2.*u3d(i+1,j,k))*tem
        else
          dum(i,j,k)=(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*                       &
                 (-u3d(i+2,j,k)+5.*u3d(i+1,j,k)+2.*u3d(i  ,j,k))*tem
        endif
        i=ni
        if((rru(i,j,k)+rru(i+1,j,k)).ge.0.)then
          dum(i,j,k)=0.5*(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*u3d(i  ,j,k)
        else
          dum(i,j,k)=0.5*(xf(i)*rru(i,j,k)+xf(i+1)*rru(i+1,j,k))*u3d(i+1,j,k)
        endif
      enddo
      ENDIF

      do j=1,nj
      do i=i1,i2
        advx(i,j,k)=-(dum(i,j,k)-dum(i-1,j,k))*rdx*uf(i)*rxf(i)
      enddo
      enddo

    ENDDO

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advu=time_advu+mytime()
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advvaxi(xh,rxh,uh,xf,vf,rho0,rr0,rf0,rrf0,dum,advx,rru,v3d,vten)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(ib:ie) :: xh,rxh,uh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v3d,vten
 
      integer i,j,k,i1,i2,j1,j2,jd1,jd2
      real, parameter :: tem = 1.0/120.0
      real, parameter :: tem2 = 1.0/12.0
 
!------------------------------------------------------------

      if(ibs.eq.1)then
        j1=2
      else
        j1=1
      endif
 
      if(ibn.eq.1)then
        j2=nj+1-1
      else
        j2=nj+1
      endif
 
!-----------------

      i1 = 1
      i2 = ni+1

      if(wbc.eq.2 .and. ibw.eq.1)then
        i1=4
      endif

      if(ebc.eq.2 .and. ibe.eq.1)then
        i2=ni-2
      endif

!-----------------

      jd1 = j1-1
      jd2 = j2

      if(sbc.eq.2 .and. ibs.eq.1)then
        jd1 = 3
      endif

      if(nbc.eq.2 .and. ibn.eq.1)then
        jd2 = nj-2
      endif

!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk

! Advection in x-direction

    if(hadvorder.eq.5)then
      do j=j1,j2
      do i=i1,i2
#ifdef SGN
        ! this form of code is used for IBM/xlf compiler:
        dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))                         &
                             *( 37.0*(v3d(i  ,j,k)+v3d(i-1,j,k))     &
                                -8.0*(v3d(i+1,j,k)+v3d(i-2,j,k))     &
                                    +(v3d(i+2,j,k)+v3d(i-3,j,k))     &
        -sign(1.0,(rru(i,j,k)+rru(i,j-1,k)))                         &
                             *( 10.0*(v3d(i  ,j,k)-v3d(i-1,j,k))     &
                                -5.0*(v3d(i+1,j,k)-v3d(i-2,j,k))     &
                                    +(v3d(i+2,j,k)-v3d(i-3,j,k)) ) )*tem
#else
        ! this seems to be faster on most other platforms:
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*                        &
                 ( 2.*v3d(i-3,j,k)-13.*v3d(i-2,j,k)+47.*v3d(i-1,j,k)    &
                  +27.*v3d(i,j,k)-3.*v3d(i+1,j,k) )*tem
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*                      &
                 ( 2.*v3d(i+2,j,k)-13.*v3d(i+1,j,k)+47.*v3d(i,j,k)    &
                  +27.*v3d(i-1,j,k)-3.*v3d(i-2,j,k) )*tem
        endif
#endif
      enddo
      enddo
    elseif(hadvorder.eq.6)then
      do j=j1,j2
      do i=i1,i2
        dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*             &
                   ( 37.0*(v3d(i  ,j,k)+v3d(i-1,j,k))     &
                     -8.0*(v3d(i+1,j,k)+v3d(i-2,j,k))     &
                         +(v3d(i+2,j,k)+v3d(i-3,j,k)) )*tem
      enddo
      enddo
    endif

      IF(wbc.eq.2 .and. ibw.eq.1)THEN
      do j=j1,j2
        i=3
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*    &
                 (-v3d(i-2,j,k)+5.*v3d(i-1,j,k)+2.*v3d(i  ,j,k))*tem2
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*    &
                 (-v3d(i+1,j,k)+5.*v3d(i  ,j,k)+2.*v3d(i-1,j,k))*tem2
        endif
        i=2
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i-1,j,k)
        else
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i  ,j,k)
        endif
        i=1
        if((rru(i,j,k)+rru(i,j-1,k)).lt.0.0)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i,j,k)
        else
          dum(i,j,k)=dum(i+1,j,k)
        endif
      enddo
      ENDIF

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
      do j=j1,j2
        i=ni-1
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*    &
                 (-v3d(i-2,j,k)+5.*v3d(i-1,j,k)+2.*v3d(i  ,j,k))*tem2
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i,j-1,k))*    &
                 (-v3d(i+1,j,k)+5.*v3d(i  ,j,k)+2.*v3d(i-1,j,k))*tem2
        endif
        i=ni
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i-1,j,k)
        else
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i  ,j,k)
        endif
        i=ni+1
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.0)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j-1,k))*v3d(i-1,j,k)
        else
          dum(i,j,k)=dum(i-1,j,k)*xf(ni)/xf(ni+1)
        endif
      enddo
      ENDIF

      do j=j1,j2
      do i=1,ni
        advx(i,j,k)=-(xf(i+1)*dum(i+1,j,k)-xf(i)*dum(i,j,k))*rdx*uh(i)*rxh(i)
      enddo
      enddo

      IF(wbc.eq.2 .and. ibw.eq.1)THEN
        do j=j1,j2
          i=1
          if((rru(i,j,k)+rru(i,j-1,k)).ge.0.0)then
            advx(i,j,k)=advx(i,j,k)-v3d(i,j,k)*0.5*(            &
                    (rru(i+1,j-1,k)-rru(i,j-1,k))               &
                   +(rru(i+1,j  ,k)-rru(i,j  ,k)) )*rdx*uh(i)
          endif
        enddo
      ENDIF

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
        do j=j1,j2
          i=ni+1
          if((rru(i,j,k)+rru(i,j-1,k)).lt.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-v3d(i,j,k)*0.5*(            &
                    (xf(i+1)*rru(i+1,j-1,k)-xf(i)*rru(i,j-1,k))               &
                   +(xf(i+1)*rru(i+1,j  ,k)-xf(i)*rru(i,j  ,k)) )*rdx*uh(i)*rxh(i)
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


      subroutine advwaxi(xh,rxh,uh,xf,vh,rho0,rr0,rf0,rrf0,dum,advx,rru,w3d,wten)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, dimension(ib:ie) :: xh,rxh,uh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: vh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w3d,wten
 
      integer i,j,k,i1,i2,j1,j2
      real, parameter :: tem = 1.0/120.0
      real, parameter :: tem2 = 1.0/12.0

!----------------------------------------------------------------

      i1 = 1
      i2 = ni+1

      if(wbc.eq.2 .and. ibw.eq.1)then
        i1=4
      endif

      if(ebc.eq.2 .and. ibe.eq.1)then
        i2=ni-2
      endif

!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=2,nk

! Advection in x-direction

    if(hadvorder.eq.5)then
      do j=1,nj
      do i=i1,i2
#ifdef SGN
        ! this form of code is used for IBM/xlf compiler:
        dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))                         &
                             *( 37.0*(w3d(i  ,j,k)+w3d(i-1,j,k))     &
                                -8.0*(w3d(i+1,j,k)+w3d(i-2,j,k))     &
                                    +(w3d(i+2,j,k)+w3d(i-3,j,k))     &
        -sign(1.0,(rru(i,j,k)+rru(i,j,k-1)))                         &
                             *( 10.0*(w3d(i  ,j,k)-w3d(i-1,j,k))     &
                                -5.0*(w3d(i+1,j,k)-w3d(i-2,j,k))     &
                                    +(w3d(i+2,j,k)-w3d(i-3,j,k)) ) )*tem
#else
        ! this seems to be faster on most other platforms:
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*                        &
                 ( 2.*w3d(i-3,j,k)-13.*w3d(i-2,j,k)+47.*w3d(i-1,j,k)    &
                  +27.*w3d(i,j,k)-3.*w3d(i+1,j,k) )*tem
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*                      &
                 ( 2.*w3d(i+2,j,k)-13.*w3d(i+1,j,k)+47.*w3d(i,j,k)    &
                  +27.*w3d(i-1,j,k)-3.*w3d(i-2,j,k) )*tem
        endif
#endif
      enddo
      enddo
    elseif(hadvorder.eq.6)then
      do j=1,nj
      do i=i1,i2
        dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*             &
                   ( 37.0*(w3d(i  ,j,k)+w3d(i-1,j,k))     &
                     -8.0*(w3d(i+1,j,k)+w3d(i-2,j,k))     &
                         +(w3d(i+2,j,k)+w3d(i-3,j,k)) )*tem
      enddo
      enddo
    endif

      IF(wbc.eq.2 .and. ibw.eq.1)THEN
      do j=1,nj
        i=3
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*     &
                 (-w3d(i-2,j,k)+5.*w3d(i-1,j,k)+2.*w3d(i  ,j,k))*tem2
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*     &
                 (-w3d(i+1,j,k)+5.*w3d(i  ,j,k)+2.*w3d(i-1,j,k))*tem2
        endif
        i=2
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i-1,j,k)
        else
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i  ,j,k)
        endif
        i=1
        if((rru(i,j,k)+rru(i,j,k-1)).lt.0.0)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i,j,k)
        else
          dum(i,j,k)=dum(i+1,j,k)
        endif
      enddo
      ENDIF

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
      do j=1,nj
        i=ni-1
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.)then
          dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*     &
                 (-w3d(i-2,j,k)+5.*w3d(i-1,j,k)+2.*w3d(i  ,j,k))*tem2
        else
          dum(i,j,k)=(rru(i,j,k)+rru(i,j,k-1))*     &
                 (-w3d(i+1,j,k)+5.*w3d(i  ,j,k)+2.*w3d(i-1,j,k))*tem2
        endif
        i=ni
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i-1,j,k)
        else
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i  ,j,k)
        endif
        i=ni+1
        if((rru(i,j,k)+rru(i,j,k-1)).ge.0.0)then
          dum(i,j,k)=0.5*(rru(i,j,k)+rru(i,j,k-1))*w3d(i-1,j,k)
        else
          dum(i,j,k)=dum(i-1,j,k)*xf(ni)/xf(ni+1)
        endif
      enddo
      ENDIF

      do j=1,nj
      do i=1,ni
        advx(i,j,k)=-(xf(i+1)*dum(i+1,j,k)-xf(i)*dum(i,j,k))*rdx*uh(i)*rxh(i)
      enddo
      enddo

      IF(wbc.eq.2 .and. ibw.eq.1)THEN
        do j=1,nj
          i=1
          if((rru(i,j,k)+rru(i,j,k-1)).ge.0.0)then
            advx(i,j,k)=advx(i,j,k)-w3d(i,j,k)*0.5*(      &
                    (rru(i+1,j,k-1)-rru(i,j,k-1))         &
                   +(rru(i+1,j,k  )-rru(i,j,k  )) )*rdx*uh(i)
          endif
        enddo
      ENDIF

      IF(ebc.eq.2 .and. ibe.eq.1)THEN
        do j=1,nj
          i=ni+1
          if((rru(i,j,k)+rru(i,j,k-1)).lt.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-w3d(i,j,k)*0.5*(      &
                    (xf(i+1)*rru(i+1,j,k-1)-xf(i)*rru(i,j,k-1))         &
                   +(xf(i+1)*rru(i+1,j,k  )-xf(i)*rru(i,j,k  )) )*rdx*uh(i)*rxh(i)
          endif
        enddo
      ENDIF

    ENDDO

!----------------------------------------------------------------
 
      if(timestats.ge.1) time_advw=time_advw+mytime()
 
      return
      end


