

      subroutine advs(nrk,wflag,bflag,bsq,xh,rxh,arh1,arh2,uh,ruh,xf,vh,rvh,gz,rgz,mh,rmh,  &
                       rho0,rr0,rf0,rrf0,advx,advy,advz,dum,divx,mass,        &
                       rru,rrv,rrw,s0,s,sten,pdef,dt,weps,                    &
                       flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,rdsf,c1,c2,rho,ri,diffit)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      integer, intent(in) :: nrk
      integer, intent(in) :: wflag,bflag
      real*8 bsq
      real, dimension(ib:ie) :: xh,rxh,arh1,arh2,uh,ruh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: vh,rvh
      real, dimension(itb:ite,jtb:jte) :: gz,rgz
      real, dimension(ib:ie,jb:je,kb:ke) :: mh,rmh,rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: advx,advy,advz,dum,divx,mass
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, dimension(ib:ie,jb:je,kb:ke) :: s0,s,sten
      integer pdef
      real, intent(in) :: dt
      double precision, intent(in) :: weps
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(in), dimension(kb:ke+1) :: rdsf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,ri
      integer, intent(in) :: diffit
 
      integer i,j,k,i1,i2,j1,j2
      logical :: doitw,doite,doits,doitn
      logical :: doweno
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
      real :: coef

!----------------------------------------------------------------

      doweno = .false.
      IF( wflag.eq.1 )THEN
        IF( (advwenos.eq.1) .or. (advwenos.eq.2.and.nrk.eq.3) ) doweno = .true.
      ENDIF

      IF(diffit.eq.1)THEN
        coef = kdiff6/64.0/dt
      ENDIF

!-----------------

      i1 = 1
      i2 = ni+1

      doitw = .false.
      doite = .false.

      IF(wbc.eq.2 .and. ibw.eq.1) doitw = .true.
      IF(ebc.eq.2 .and. ibe.eq.1) doite = .true.

!-----------------

      j1 = 1
      j2 = nj+1

      doits = .false.
      doitn = .false.

      IF(sbc.eq.2 .and. ibs.eq.1) doits = .true.
      IF(nbc.eq.2 .and. ibn.eq.1) doitn = .true.

!-----------------

    hadvsection:  IF(axisymm.eq.1)THEN
      call advsaxi(doweno,bflag,bsq,xh,rxh,arh1,arh2,uh,ruh,xf,vh,rvh,rmh,gz,rgz, &
                   rho0,rr0,rf0,rrf0,advx,dum,mass,rru,s0,s,sten,pdef,dt,weps, &
                   flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        advy(i,j,k)=0.0
      enddo
      enddo
      enddo
    ELSE

!-----------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k,s1,s2,s3,s4,s5,f1,f2,f3,b1,b2,b3,w1,w2,w3,a1,a2,a3,a4,bmax,doit)
    DO k=1,nk   ! start of k-loop

      budx(k) = 0.0d0
      budy(k) = 0.0d0

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
        if( bmax.lt.Min(1.0d-20,weps) ) doit = .false.
      ENDIF

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
    elseif(hadvordrs.eq.5)then
      do j=1,nj
      do i=i1,i2
#ifdef SGN
        ! this form of code is used for IBM/xlf compiler:
        dum(i,j,k)=rru(i,j,k)*( 37.0*(s(i  ,j,k)+s(i-1,j,k))     &
                                -8.0*(s(i+1,j,k)+s(i-2,j,k))     &
                                    +(s(i+2,j,k)+s(i-3,j,k))     &
        -sign(1.0,rru(i,j,k))*( 10.0*(s(i  ,j,k)-s(i-1,j,k))     &
                                -5.0*(s(i+1,j,k)-s(i-2,j,k))     &
                                    +(s(i+2,j,k)-s(i-3,j,k)) ) )*onedsixty
#else
        ! this seems to be faster on most other platforms:
        if(rru(i,j,k).ge.0.)then
          dum(i,j,k)=rru(i,j,k)*( 2.*s(i-3,j,k)-13.*s(i-2,j,k)   &
                +47.*s(i-1,j,k)+27.*s(i,j,k)-3.*s(i+1,j,k) )*onedsixty
        else
          dum(i,j,k)=rru(i,j,k)*( 2.*s(i+2,j,k)-13.*s(i+1,j,k)   &
                +47.*s(i,j,k)+27.*s(i-1,j,k)-3.*s(i-2,j,k) )*onedsixty
        endif
#endif
      enddo
      enddo
    elseif(hadvordrs.eq.6)then
      do j=1,nj
      do i=i1,i2
        dum(i,j,k)=rru(i,j,k)*( 37.0*(s(i  ,j,k)+s(i-1,j,k))     &
                                -8.0*(s(i+1,j,k)+s(i-2,j,k))     &
                                    +(s(i+2,j,k)+s(i-3,j,k)) )*onedsixty
      enddo
      enddo
    endif

      if(doitw)then
      do j=1,nj
        i=1
        if(rru(i,j,k).ge.0.0)then
          dum(i,j,k)=dum(i+1,j,k)
        endif
      enddo
      endif

      if(doite)then
      do j=1,nj
        i=ni+1
        if(rru(i,j,k).le.0.0)then
          dum(i,j,k)=dum(i-1,j,k)
        endif
      enddo
      endif

      do j=1,nj
      do i=1,ni
        advx(i,j,k)=-(dum(i+1,j,k)-dum(i,j,k))*rdx*uh(i)
      enddo
      enddo

      IF(doitw)THEN
        do j=1,nj
          if(rru(1,j,k).ge.0.0)then
            i=1
            advx(i,j,k)=advx(i,j,k)-s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)
          endif
        enddo
      ENDIF

      IF(doite)THEN
        do j=1,nj
          if(rru(ni+1,j,k).le.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)
          endif
        enddo
      ENDIF

    !-------------------------------------------------------
    ! 6th-order diffusion:
    IF(diffit.eq.1)THEN
      do j=1,nj
      do i=1,ni+1
        dum(i,j,k)=( 10.0*(s(i  ,j,k)-s(i-1,j,k))     &
                     -5.0*(s(i+1,j,k)-s(i-2,j,k))     &
                         +(s(i+2,j,k)-s(i-3,j,k)) )   &
                  *0.5*(rho(i-1,j,k)+rho(i,j,k))
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj
        do i=1,ni+1
          if( dum(i,j,k)*(s(i,j,k)-s(i-1,j,k)).le.0.0 )then
            dum(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj
      do i=1,ni
        advx(i,j,k)=advx(i,j,k)+coef*(dum(i+1,j,k)-dum(i,j,k))*ri(i,j,k)*rho0(i,j,k)
      enddo
      enddo
    ENDIF
    !-------------------------------------------------------

! Advection in y-direction

    if(doweno)then
      do j=j1,j2
      do i=1,ni
        if(rrv(i,j,k).ge.0.0)then
          s1=s(i,j-3,k)
          s2=s(i,j-2,k)
          s3=s(i,j-1,k)
          s4=s(i,j  ,k)
          s5=s(i,j+1,k)
        else
          s1=s(i,j+2,k)
          s2=s(i,j+1,k)
          s3=s(i,j  ,k)
          s4=s(i,j-1,k)
          s5=s(i,j-2,k)
        endif

      doit = .true.
      IF( pdef.eq.1 )THEN
        bmax = max(s1,s2,s3,s4,s5)
        if( bmax.lt.Min(1.0d-20,weps)  ) doit = .false.
      ENDIF

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

        dum(i,j,k)=rrv(i,j,k)*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      ELSE
        dum(i,j,k)=0.0
      ENDIF
      enddo
      enddo
    elseif(hadvordrs.eq.5)then
      do j=j1,j2
      do i=1,ni
#ifdef SGN
        ! this form of code is used for IBM/xlf compiler:
        dum(i,j,k)=rrv(i,j,k)*( 37.0*(s(i,j  ,k)+s(i,j-1,k))     &
                                -8.0*(s(i,j+1,k)+s(i,j-2,k))     &
                                    +(s(i,j+2,k)+s(i,j-3,k))     &
        -sign(1.0,rrv(i,j,k))*( 10.0*(s(i,j  ,k)-s(i,j-1,k))     &
                                -5.0*(s(i,j+1,k)-s(i,j-2,k))     &
                                    +(s(i,j+2,k)-s(i,j-3,k)) ) )*onedsixty
#else
        ! this seems to be faster on most other platforms:
        if(rrv(i,j,k).ge.0.)then
          dum(i,j,k)=rrv(i,j,k)*( 2.*s(i,j-3,k)-13.*s(i,j-2,k)   &
                +47.*s(i,j-1,k)+27.*s(i,j,k)-3.*s(i,j+1,k) )*onedsixty
        else
          dum(i,j,k)=rrv(i,j,k)*( 2.*s(i,j+2,k)-13.*s(i,j+1,k)   &
                +47.*s(i,j,k)+27.*s(i,j-1,k)-3.*s(i,j-2,k) )*onedsixty
        endif
#endif
      enddo
      enddo
    elseif(hadvordrs.eq.6)then
      do j=j1,j2
      do i=1,ni
        dum(i,j,k)=rrv(i,j,k)*( 37.0*(s(i,j  ,k)+s(i,j-1,k))     &
                                -8.0*(s(i,j+1,k)+s(i,j-2,k))     &
                                    +(s(i,j+2,k)+s(i,j-3,k)) )*onedsixty
      enddo
      enddo
    endif

      if(doits)then
      do i=1,ni
        j=1
        if(rrv(i,j,k).ge.0.0)then
          dum(i,j,k)=dum(i,j+1,k)
        endif
      enddo
      endif

      if(doitn)then
      do i=1,ni
        j=nj+1
        if(rrv(i,j,k).le.0.0)then
          dum(i,j,k)=dum(i,j-1,k)
        endif
      enddo
      endif

      do j=1,nj
      do i=1,ni
        advy(i,j,k)=-(dum(i,j+1,k)-dum(i,j,k))*rdy*vh(j)
      enddo
      enddo

      IF(doits)THEN
        do i=1,ni
          if(rrv(i,1,k).ge.0.0)then
            j=1
            advy(i,j,k)=advy(i,j,k)-s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)
          endif
        enddo
      ENDIF

      IF(doitn)THEN
        do i=1,ni
          if(rrv(i,nj+1,k).le.0.0)then
            j=nj
            advy(i,j,k)=advy(i,j,k)-s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)
          endif
        enddo
      ENDIF

    !-------------------------------------------------------
    ! 6th-order diffusion:
    IF(diffit.eq.1)THEN
      do j=1,nj+1
      do i=1,ni
        dum(i,j,k)=( 10.0*(s(i,j  ,k)-s(i,j-1,k))     &
                     -5.0*(s(i,j+1,k)-s(i,j-2,k))     &
                         +(s(i,j+2,k)-s(i,j-3,k)) )   &
                  *0.5*(rho(i,j-1,k)+rho(i,j,k))
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj+1
        do i=1,ni
          if( dum(i,j,k)*(s(i,j,k)-s(i,j-1,k)).le.0.0 )then
            dum(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj
      do i=1,ni
        advy(i,j,k)=advy(i,j,k)+coef*(dum(i,j+1,k)-dum(i,j,k))*ri(i,j,k)*rho0(i,j,k)
      enddo
      enddo
    ENDIF
    !-------------------------------------------------------

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
!  Misc for y-direction

      IF(stat_qsrc.eq.1.and.(sbc.eq.2.or.nbc.eq.2).and.bflag.eq.1)THEN
        tem0=dt*dx*dz
        do k=1,nk
          bsq=bsq+budy(k)*tem0
        enddo
      ENDIF

      IF(pdscheme.eq.1 .and. pdef.eq.1)THEN
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call pdefy(vh,rho0,gz,rgz,rrv,advy,dum,mass,s0,s,dt,flag,ss31,ss32,sn31,sn32)
      ENDIF

    ENDIF  hadvsection

!----------------------------------------------------------------
! Advection in z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k,s1,s2,s3,s4,s5,f1,f2,f3,b1,b2,b3,w1,w2,w3,a1,a2,a3,a4,bmax,dd,rr,phi,doit)
  jloops:  DO j=1,nj

    IF(doweno)THEN

      do k=3,nk-1
      do i=1,ni
        if(rrw(i,j,k).ge.0.0)then
          s1=s(i,j,k-3)
          s2=s(i,j,k-2)
          s3=s(i,j,k-1)
          s4=s(i,j,k  )
          s5=s(i,j,k+1)
        else
          s1=s(i,j,k+2)
          s2=s(i,j,k+1)
          s3=s(i,j,k  )
          s4=s(i,j,k-1)
          s5=s(i,j,k-2)
        endif

      doit = .true.
      IF( pdef.eq.1 )THEN
        bmax = max(s1,s2,s3,s4,s5)
        if( bmax.lt.Min(1.0d-20,weps)  ) doit = .false.
      ENDIF
      IF( k.eq.3 .and. rrw(i,j,k).gt.0.0 ) doit = .false.
      IF( k.eq.(nk-1) .and. rrw(i,j,k).lt.0.0 ) doit = .false.

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

        dum(i,j,k)=rrw(i,j,k)*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      ELSE
        dum(i,j,k)=0.0
      ENDIF
      enddo
      enddo

      ! flux at k=3 if w > 0
      k = 3
      do i=1,ni
        if( rrw(i,j,k).ge.0.0 )then
          dd = s(i,j,k-1)-s(i,j,k-2)
          rr = (s(i,j,k)-s(i,j,k-1))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = rrw(i,j,k)*( s(i,j,k-1) + 0.5*phi*(s(i,j,k-1)-s(i,j,k-2)) )
        endif
      enddo

      ! flux at k=(nk-1) if w < 0
      k = nk-1
      do i=1,ni
        if( rrw(i,j,k).le.0.0 )then
          dd = s(i,j,k)-s(i,j,k+1)
          rr = (s(i,j,k-1)-s(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = rrw(i,j,k)*( s(i,j,k) + 0.5*phi*(s(i,j,k)-s(i,j,k+1)) )
        endif
      enddo

      k = 2
      do i=1,ni
        if( rrw(i,j,k).ge.0.0 )then
          dum(i,j,k)=rrw(i,j,k)*(c1(i,j,k)*s(i,j,k-1)+c2(i,j,k)*s(i,j,k))
        else
          dd = s(i,j,k)-s(i,j,k+1)
          rr = (s(i,j,k-1)-s(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = rrw(i,j,k)*( s(i,j,k) + 0.5*phi*(s(i,j,k)-s(i,j,k+1)) )
        endif
      enddo

      k = nk
      do i=1,ni
        if( rrw(i,j,k).gt.0.0 )then
          dd = s(i,j,k-1)-s(i,j,k-2)
          rr = (s(i,j,k)-s(i,j,k-1))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = rrw(i,j,k)*( s(i,j,k-1) + 0.5*phi*(s(i,j,k-1)-s(i,j,k-2)) )
        else
          dum(i,j,k)=rrw(i,j,k)*(c1(i,j,k)*s(i,j,k-1)+c2(i,j,k)*s(i,j,k))
        endif
      enddo

    ELSEIF(vadvordrs.eq.5)THEN

      do k=4,nk-2
      do i=1,ni
        if(rrw(i,j,k).ge.0.)then
          dum(i,j,k)=rrw(i,j,k)*( 2.*s(i,j,k-3)-13.*s(i,j,k-2)      &
                +47.*s(i,j,k-1)+27.*s(i,j,k)-3.*s(i,j,k+1) )*onedsixty
        else
          dum(i,j,k)=rrw(i,j,k)*( 2.*s(i,j,k+2)-13.*s(i,j,k+1)      &
                +47.*s(i,j,k)+27.*s(i,j,k-1)-3.*s(i,j,k-2) )*onedsixty
        endif
      enddo
      enddo

      k = 3
      do i=1,ni
        if(rrw(i,j,k).ge.0.)then
          dum(i,j,k)=rrw(i,j,k)*(-s(i,j,k-2)+5.*s(i,j,k-1)+2.*s(i,j,k))*onedsix
        else
          dum(i,j,k)=rrw(i,j,k)*( 2.*s(i,j,k+2)-13.*s(i,j,k+1)      &
                +47.*s(i,j,k)+27.*s(i,j,k-1)-3.*s(i,j,k-2) )*onedsixty
        endif
      enddo

      k = nk-1
      do i=1,ni
        if(rrw(i,j,k).ge.0.)then
          dum(i,j,k)=rrw(i,j,k)*( 2.*s(i,j,k-3)-13.*s(i,j,k-2)      &
                +47.*s(i,j,k-1)+27.*s(i,j,k)-3.*s(i,j,k+1) )*onedsixty
        else
          dum(i,j,k)=rrw(i,j,k)*(-s(i,j,k+1)+5.*s(i,j,k)+2.*s(i,j,k-1))*onedsix
        endif
      enddo

      k = 2
      do i=1,ni
        if(rrw(i,j,k).ge.0.)then
          dum(i,j,k)=rrw(i,j,k)*(c1(i,j,k)*s(i,j,k-1)+c2(i,j,k)*s(i,j,k))
        else
          dum(i,j,k)=rrw(i,j,k)*(-s(i,j,k+1)+5.*s(i,j,k)+2.*s(i,j,k-1))*onedsix
        endif
      enddo

      k = nk
      do i=1,ni
        if(rrw(i,j,k).ge.0.)then
          dum(i,j,k)=rrw(i,j,k)*(-s(i,j,k-2)+5.*s(i,j,k-1)+2.*s(i,j,k))*onedsix
        else
          dum(i,j,k)=rrw(i,j,k)*(c1(i,j,k)*s(i,j,k-1)+c2(i,j,k)*s(i,j,k))
        endif
      enddo

    ELSEIF(vadvordrs.eq.6)THEN

      do k=4,nk-2
      do i=1,ni
        dum(i,j,k)=rrw(i,j,k)                                 &
                            *( 37.0*(s(i,j,k  )+s(i,j,k-1))          &
                               -8.0*(s(i,j,k+1)+s(i,j,k-2))          &
                                   +(s(i,j,k+2)+s(i,j,k-3)) )*onedsixty
      enddo
      enddo

      do k=3,(nk-1),(nk-4)
      do i=1,ni
        dum(i,j,k)=rrw(i,j,k)                                &
                            *( 7.0*(s(i,j,k  )+s(i,j,k-1))          &
                                  -(s(i,j,k+1)+s(i,j,k-2)) )*onedtwelve
      enddo
      enddo

      do k=2,nk,(nk-2)
      do i=1,ni
        dum(i,j,k)=rrw(i,j,k)*(c1(i,j,k)*s(i,j,k-1)+c2(i,j,k)*s(i,j,k))
      enddo
      enddo

    ENDIF

!------

    IF(terrain_flag)THEN

      k=1
      do i=1,ni
        advz(i,j,k)=-dum(i,j,k+1)*rdsf(k)
      enddo
      do k=2,nk-1
      do i=1,ni
        advz(i,j,k)=-(dum(i,j,k+1)-dum(i,j,k))*rdsf(k)
      enddo
      enddo
      k=nk
      do i=1,ni
        advz(i,j,k)=+dum(i,j,k)*rdsf(k)
      enddo

    ELSE

      k=1
      do i=1,ni
        advz(i,j,k)=-dum(i,j,k+1)*rdz*mh(1,1,k)
      enddo
      do k=2,nk-1
      do i=1,ni
        advz(i,j,k)=-(dum(i,j,k+1)-dum(i,j,k))*rdz*mh(1,1,k)
      enddo
      enddo
      k=nk
      do i=1,ni
        advz(i,j,k)=+dum(i,j,k)*rdz*mh(1,1,k)
      enddo

    ENDIF

    ENDDO  jloops

!----------------------------------------------------------------
!  Misc for z-direction

      IF(pdscheme.eq.1 .and. pdef.eq.1)THEN
        if(timestats.ge.1) time_advs=time_advs+mytime()
        call pdefz(mh,rho0,gz,rgz,rdsf,rrw,advz,dum,mass,s0,s,dt,flag)
      ENDIF

!----------------------------------------------------------------
!  Total advection tendency:

    IF(terrain_flag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                                 +s(i,j,k)*divx(i,j,k) )*rr0(i,j,k)*gz(i,j)
      enddo
      enddo
      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                                 +s(i,j,k)*divx(i,j,k) )*rr0(1,1,k)
      enddo
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


      subroutine advu(nrk,arh1,arh2,xf,rxf,arf1,arf2,uf,vh,gz,rgz,gzu,mh,rho0,rr0,rf0,rrf0,dum,advx,advy,advz,divx, &
                       rru,u3d,uten,rrv,rrw,rdsf,c1,c2,rho,dt)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      integer, intent(in) :: nrk
      real, dimension(ib:ie) :: arh1,arh2
      real, dimension(ib:ie+1) :: xf,rxf,arf1,arf2,uf
      real, dimension(jb:je) :: vh
      real, dimension(itb:ite,jtb:jte) :: gz,rgz,gzu
      real, dimension(ib:ie,jb:je,kb:ke) :: mh,rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advx,advy,advz,divx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru,u3d,uten
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, intent(in), dimension(kb:ke+1) :: rdsf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho
      real, intent(in) :: dt
 
      integer i,j,k,i1,i2,j1,j2,id1,id2
      real :: ubar,vbar,wbar,cc1,cc2
      logical :: doitw,doite,doits,doitn

      logical :: doweno
      real :: dd,rr,phi
      real :: s1,s2,s3,s4,s5
      real :: f1,f2,f3
      real :: b1,b2,b3
      real :: w1,w2,w3
      double precision :: a1,a2,a3,a4
      double precision :: weps
      real :: coef

!------------------------------------------------------------

      doweno = .false.
      IF( (advwenov.eq.1) .or. (advwenov.eq.2.and.nrk.eq.3) ) doweno = .true.
      weps = 100.0*epsilon

      IF( idiff.ge.1 .and. difforder.eq.6 )THEN
        coef = kdiff6/64.0/dt
      ENDIF

!-----------------

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

      doitw = .false.
      doite = .false.

      IF(wbc.eq.2 .and. ibw.eq.1) doitw = .true.
      IF(ebc.eq.2 .and. ibe.eq.1) doite = .true.

!-----------------

      j1 = 1
      j2 = nj+1

      doits = .false.
      doitn = .false.

      IF(sbc.eq.2 .and. ibs.eq.1) doits = .true.
      IF(nbc.eq.2 .and. ibn.eq.1) doitn = .true.

!----------------------------------------------------------------

    hadvsection:  IF(axisymm.eq.1)THEN
      call advuaxi(doweno,arh1,arh2,xf,rxf,arf1,arf2,uf,vh,rho0,rr0,rf0,rrf0,dum,advx,rru,u3d,uten)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=i1,i2
        advy(i,j,k)=0.0
      enddo
      enddo
      enddo
    ELSE

!-----------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk

! Advection in x-direction

    if(doweno)then
      do j=1,nj
      do i=id1,id2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
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
    elseif(hadvordrv.eq.5)then
      do j=1,nj
      do i=id1,id2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        if(ubar.ge.0.)then
          dum(i,j,k)=ubar*( 2.*u3d(i-2,j,k)-13.*u3d(i-1,j,k)+47.*u3d(i,j,k)   &
                          +27.*u3d(i+1,j,k)-3.*u3d(i+2,j,k) )*onedsixty
        else
          dum(i,j,k)=ubar*( 2.*u3d(i+3,j,k)-13.*u3d(i+2,j,k)+47.*u3d(i+1,j,k)   &
                          +27.*u3d(i,j,k)-3.*u3d(i-1,j,k) )*onedsixty
        endif
      enddo
      enddo
    elseif(hadvordrv.eq.6)then
      do j=1,nj
      do i=id1,id2
        ubar = 0.5*(rru(i,j,k)+rru(i+1,j,k))
        dum(i,j,k)=ubar*( 37.0*(u3d(i+1,j,k)+u3d(i  ,j,k)) &
                          -8.0*(u3d(i+2,j,k)+u3d(i-1,j,k)) &
                              +(u3d(i+3,j,k)+u3d(i-2,j,k)) )*onedsixty
      enddo
      enddo
    endif

      do j=1,nj
      do i=i1,i2
        advx(i,j,k)=-(dum(i,j,k)-dum(i-1,j,k))*rdx*uf(i)
      enddo
      enddo

    !-------------------------------------------------------
    ! 6th-order diffusion:
    IF( idiff.ge.1 .and. difforder.eq.6 )THEN
      do j=1,nj
      do i=1,ni+2
        dum(i,j,k)=( 10.0*(u3d(i  ,j,k)-u3d(i-1,j,k))     &
                     -5.0*(u3d(i+1,j,k)-u3d(i-2,j,k))     &
                         +(u3d(i+2,j,k)-u3d(i-3,j,k)) )*rho(i-1,j,k)
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj
        do i=1,ni+2
          if( dum(i,j,k)*(u3d(i,j,k)-u3d(i-1,j,k)).le.0.0 )then
            dum(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj
      do i=1,ni+1
        advx(i,j,k)=advx(i,j,k)+coef*(dum(i+1,j,k)-dum(i,j,k))*(rho0(i-1,j,k)+rho0(i,j,k))/(rho(i-1,j,k)+rho(i,j,k))
      enddo
      enddo
    ENDIF
    !-------------------------------------------------------

! Advection in y-direction

    if(doweno)then
      do j=j1,j2
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        if(vbar.ge.0.0)then
          s1=u3d(i,j-3,k)
          s2=u3d(i,j-2,k)
          s3=u3d(i,j-1,k)
          s4=u3d(i,j  ,k)
          s5=u3d(i,j+1,k)
        else
          s1=u3d(i,j+2,k)
          s2=u3d(i,j+1,k)
          s3=u3d(i,j  ,k)
          s4=u3d(i,j-1,k)
          s5=u3d(i,j-2,k)
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

        dum(i,j,k)=vbar*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo
    elseif(hadvordrv.eq.5)then
      do j=j1,j2
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        if(vbar.ge.0.)then
          dum(i,j,k)=vbar*( 2.*u3d(i,j-3,k)-13.*u3d(i,j-2,k)+47.*u3d(i,j-1,k)   &
                          +27.*u3d(i,j,k)-3.*u3d(i,j+1,k) )*onedsixty
        else
          dum(i,j,k)=vbar*( 2.*u3d(i,j+2,k)-13.*u3d(i,j+1,k)+47.*u3d(i,j,k)   &
                          +27.*u3d(i,j-1,k)-3.*u3d(i,j-2,k) )*onedsixty
        endif
      enddo
      enddo
    elseif(hadvordrv.eq.6)then
      do j=j1,j2
      do i=i1,i2
        vbar = 0.5*(rrv(i,j,k)+rrv(i-1,j,k))
        dum(i,j,k)=vbar*( 37.0*(u3d(i,j  ,k)+u3d(i,j-1,k)) &
                          -8.0*(u3d(i,j+1,k)+u3d(i,j-2,k)) &
                              +(u3d(i,j+2,k)+u3d(i,j-3,k)) )*onedsixty
      enddo
      enddo
    endif

      if(doits)then
      do i=i1,i2
        j=1
        if((rrv(i,j,k)+rrv(i-1,j,k)).ge.0.0)then
          dum(i,j,k)=dum(i,j+1,k)
        endif
      enddo
      endif

      if(doitn)then
      do i=i1,i2
        j=nj+1
        if((rrv(i,j,k)+rrv(i-1,j,k)).le.0.0)then
          dum(i,j,k)=dum(i,j-1,k)
        endif
      enddo
      endif

      do j=1,nj
      do i=i1,i2
        advy(i,j,k)=-(dum(i,j+1,k)-dum(i,j,k))*rdy*vh(j)
      enddo
      enddo

      IF(doits)THEN
        do i=i1,i2
          if((rrv(i,1,k)+rrv(i-1,1,k)).ge.0.0)then
            j=1
            advy(i,j,k)=advy(i,j,k)-u3d(i,j,k)*0.5*(                    &
                            (rrv(i-1,j+1,k)-rrv(i-1,j,k))               &
                           +(rrv(i  ,j+1,k)-rrv(i  ,j,k)) )*rdy*vh(j)
          endif
        enddo
      ENDIF

      IF(doitn)THEN
        do i=i1,i2
          if((rrv(i,nj+1,k)+rrv(i-1,nj+1,k)).le.0.0)then
            j=nj
            advy(i,j,k)=advy(i,j,k)-u3d(i,j,k)*0.5*(                    &
                            (rrv(i-1,j+1,k)-rrv(i-1,j,k))               &
                           +(rrv(i  ,j+1,k)-rrv(i  ,j,k)) )*rdy*vh(j)
          endif
        enddo
      ENDIF

    !-------------------------------------------------------
    ! 6th-order diffusion:
    IF( idiff.ge.1 .and. difforder.eq.6 )THEN
      do j=1,nj+1
      do i=1,ni+1
        dum(i,j,k)=( 10.0*(u3d(i,j  ,k)-u3d(i,j-1,k))     &
                     -5.0*(u3d(i,j+1,k)-u3d(i,j-2,k))     &
                         +(u3d(i,j+2,k)-u3d(i,j-3,k)) )   &
                  *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))     &
                         +(rho(i-1,j,k)+rho(i,j-1,k)) )
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj+1
        do i=1,ni+1
          if( dum(i,j,k)*(u3d(i,j,k)-u3d(i,j-1,k)).le.0.0 )then
            dum(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj
      do i=1,ni+1
        advy(i,j,k)=advy(i,j,k)+coef*(dum(i,j+1,k)-dum(i,j,k))*(rho0(i-1,j,k)+rho0(i,j,k))/(rho(i-1,j,k)+rho(i,j,k))
      enddo
      enddo
    ENDIF
    !-------------------------------------------------------

    ENDDO

    ENDIF  hadvsection

!----------------------------------------------------------------
! Advection in z-direction  (Cartesian grid)

  vadvu:  IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar,cc1,cc2)
  jloopu:  DO j=1,nj

    IF(doweno)THEN

      do k=3,nk-1
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          s1=u3d(i,j,k-3)
          s2=u3d(i,j,k-2)
          s3=u3d(i,j,k-1)
          s4=u3d(i,j,k  )
          s5=u3d(i,j,k+1)
        else
          s1=u3d(i,j,k+2)
          s2=u3d(i,j,k+1)
          s3=u3d(i,j,k  )
          s4=u3d(i,j,k-1)
          s5=u3d(i,j,k-2)
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

        dum(i,j,k)=wbar*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo

      ! flux at k=3 if w > 0
      k = 3
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dd = u3d(i,j,k-1)-u3d(i,j,k-2)
          rr = (u3d(i,j,k)-u3d(i,j,k-1))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( u3d(i,j,k-1) + 0.5*phi*(u3d(i,j,k-1)-u3d(i,j,k-2)) )
        endif
      enddo

      ! flux at k=(nk-1) if w < 0
      k = nk-1
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.lt.0.0 )then
          dd = u3d(i,j,k)-u3d(i,j,k+1)
          rr = (u3d(i,j,k-1)-u3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( u3d(i,j,k) + 0.5*phi*(u3d(i,j,k)-u3d(i,j,k+1)) )
        endif
      enddo

      k = 2
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.ge.0.0 )then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dum(i,j,k)=wbar*(cc1*u3d(i,j,k-1)+cc2*u3d(i,j,k))
        else
          dd = u3d(i,j,k)-u3d(i,j,k+1)
          rr = (u3d(i,j,k-1)-u3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( u3d(i,j,k) + 0.5*phi*(u3d(i,j,k)-u3d(i,j,k+1)) )
        endif
      enddo

      k = nk
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dd = u3d(i,j,k-1)-u3d(i,j,k-2)
          rr = (u3d(i,j,k)-u3d(i,j,k-1))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( u3d(i,j,k-1) + 0.5*phi*(u3d(i,j,k-1)-u3d(i,j,k-2)) )
        else
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dum(i,j,k)=wbar*(cc1*u3d(i,j,k-1)+cc2*u3d(i,j,k))
        endif
      enddo

    ELSEIF(vadvordrv.eq.5)THEN

      do k=4,nk-2
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*( 2.*u3d(i,j,k-3)-13.*u3d(i,j,k-2)+47.*u3d(i,j,k-1)   &
                          +27.*u3d(i,j,k)-3.*u3d(i,j,k+1) )*onedsixty
        else
          dum(i,j,k)=wbar*( 2.*u3d(i,j,k+2)-13.*u3d(i,j,k+1)+47.*u3d(i,j,k)   &
                          +27.*u3d(i,j,k-1)-3.*u3d(i,j,k-2) )*onedsixty
        endif
      enddo
      enddo

      k = 3
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if((rrw(i,j,k)+rrw(i-1,j,k)).ge.0.)then
          dum(i,j,k)=wbar*(-u3d(i,j,k-2)+5.*u3d(i,j,k-1)+2.*u3d(i,j,k  ))*onedsix
        else
          dum(i,j,k)=wbar*( 2.*u3d(i,j,k+2)-13.*u3d(i,j,k+1)+47.*u3d(i,j,k)   &
                          +27.*u3d(i,j,k-1)-3.*u3d(i,j,k-2) )*onedsixty
        endif
      enddo

      k = nk-1
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if((rrw(i,j,k)+rrw(i-1,j,k)).ge.0.)then
          dum(i,j,k)=wbar*( 2.*u3d(i,j,k-3)-13.*u3d(i,j,k-2)+47.*u3d(i,j,k-1)   &
                          +27.*u3d(i,j,k)-3.*u3d(i,j,k+1) )*onedsixty
        else
          dum(i,j,k)=wbar*(-u3d(i,j,k+1)+5.*u3d(i,j,k  )+2.*u3d(i,j,k-1))*onedsix
        endif
      enddo

      k = 2
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if((rrw(i,j,k)+rrw(i-1,j,k)).ge.0.)then
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dum(i,j,k)=wbar*(cc1*u3d(i,j,k-1)+cc2*u3d(i,j,k))
        else
          dum(i,j,k)=wbar*(-u3d(i,j,k+1)+5.*u3d(i,j,k  )+2.*u3d(i,j,k-1))*onedsix
        endif
      enddo

      k = nk
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        if((rrw(i,j,k)+rrw(i-1,j,k)).ge.0.)then
          dum(i,j,k)=wbar*(-u3d(i,j,k-2)+5.*u3d(i,j,k-1)+2.*u3d(i,j,k  ))*onedsix
        else
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dum(i,j,k)=wbar*(cc1*u3d(i,j,k-1)+cc2*u3d(i,j,k))
        endif
      enddo

    ELSEIF(vadvordrv.eq.6)THEN

      do k=4,nk-2
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        dum(i,j,k)=wbar*( 37.0*(u3d(i,j,k  )+u3d(i,j,k-1)) &
                          -8.0*(u3d(i,j,k+1)+u3d(i,j,k-2)) &
                              +(u3d(i,j,k+2)+u3d(i,j,k-3)) )*onedsixty
      enddo
      enddo

      do k=3,(nk-1),(nk-4)
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        dum(i,j,k)=wbar*( 7.0*(u3d(i,j,k  )+u3d(i,j,k-1)) &
                             -(u3d(i,j,k+1)+u3d(i,j,k-2)) )*onedtwelve
      enddo
      enddo

      do k=2,nk,(nk-2)
      do i=i1,i2
        wbar = 0.5*(rrw(i,j,k)+rrw(i-1,j,k))
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        dum(i,j,k)=wbar*(cc1*u3d(i,j,k-1)+cc2*u3d(i,j,k))
      enddo
      enddo

    ENDIF

!------

      IF(terrain_flag)THEN

        k=1
        do i=i1,i2
          advz(i,j,k)=-dum(i,j,k+1)*rdsf(k)
          uten(i,j,k)=uten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(divx(i,j,k)+divx(i-1,j,k)) )         &
                  *gzu(i,j)/(0.5*(rho0(i-1,j,k)+rho0(i,j,k)))
        enddo
        do k=2,nk-1
        do i=i1,i2
          advz(i,j,k)=-(dum(i,j,k+1)-dum(i,j,k))*rdsf(k)
          uten(i,j,k)=uten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(divx(i,j,k)+divx(i-1,j,k)) )         &
                  *gzu(i,j)/(0.5*(rho0(i-1,j,k)+rho0(i,j,k)))
        enddo
        enddo
        k=nk
        do i=i1,i2
          advz(i,j,k)=+dum(i,j,k)*rdsf(k)
          uten(i,j,k)=uten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(divx(i,j,k)+divx(i-1,j,k)) )         &
                  *gzu(i,j)/(0.5*(rho0(i-1,j,k)+rho0(i,j,k)))
        enddo

      ELSE

        k=1
        do i=i1,i2
          advz(i,j,k)=-dum(i,j,k+1)*rdz*mh(1,1,k)
          uten(i,j,k)=uten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(divx(i,j,k)+divx(i-1,j,k)) )*rr0(1,1,k)
        enddo
        do k=2,nk-1
        do i=i1,i2
          advz(i,j,k)=-(dum(i,j,k+1)-dum(i,j,k))*rdz*mh(1,1,k)
          uten(i,j,k)=uten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(divx(i,j,k)+divx(i-1,j,k)) )*rr0(1,1,k)
        enddo
        enddo
        k=nk
        do i=i1,i2
          advz(i,j,k)=+dum(i,j,k)*rdz*mh(1,1,k)
          uten(i,j,k)=uten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(divx(i,j,k)+divx(i-1,j,k)) )*rr0(1,1,k)
        enddo

      ENDIF

    ENDDO    jloopu

!  end vadvu for Cartesian grid
!----------------------------------------------------------------
! Advection in z-direction  (axisymmetric grid)

  ELSEIF(axisymm.eq.1)THEN

    IF(ebc.eq.3.or.ebc.eq.4) i2 = ni

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
  jloopuasymm:  DO j=1,nj

    IF(doweno)THEN

      do k=3,nk-1
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.0)then
          s1=u3d(i,j,k-3)
          s2=u3d(i,j,k-2)
          s3=u3d(i,j,k-1)
          s4=u3d(i,j,k  )
          s5=u3d(i,j,k+1)
        else
          s1=u3d(i,j,k+2)
          s2=u3d(i,j,k+1)
          s3=u3d(i,j,k  )
          s4=u3d(i,j,k-1)
          s5=u3d(i,j,k-2)
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

        dum(i,j,k)=wbar*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo

      ! flux at k=3 if w > 0
      k = 3
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dd = u3d(i,j,k-1)-u3d(i,j,k-2)
          rr = (u3d(i,j,k)-u3d(i,j,k-1))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( u3d(i,j,k-1) + 0.5*phi*(u3d(i,j,k-1)-u3d(i,j,k-2)) )
        endif
      enddo

      ! flux at k=(nk-1) if w < 0
      k = nk-1
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if( wbar.lt.0.0 )then
          dd = u3d(i,j,k)-u3d(i,j,k+1)
          rr = (u3d(i,j,k-1)-u3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( u3d(i,j,k) + 0.5*phi*(u3d(i,j,k)-u3d(i,j,k+1)) )
        endif
      enddo

      k = 2
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if( wbar.ge.0.0 )then
          dum(i,j,k)=wbar*(c1(1,1,k)*u3d(i,j,k-1)+c2(1,1,k)*u3d(i,j,k))
        else
          dd = u3d(i,j,k)-u3d(i,j,k+1)
          rr = (u3d(i,j,k-1)-u3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( u3d(i,j,k) + 0.5*phi*(u3d(i,j,k)-u3d(i,j,k+1)) )
        endif
      enddo

      k = nk
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if( wbar.gt.0.0 )then
          dd = u3d(i,j,k-1)-u3d(i,j,k-2)
          rr = (u3d(i,j,k)-u3d(i,j,k-1))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( u3d(i,j,k-1) + 0.5*phi*(u3d(i,j,k-1)-u3d(i,j,k-2)) )
        else
          dum(i,j,k)=wbar*(c1(1,1,k)*u3d(i,j,k-1)+c2(1,1,k)*u3d(i,j,k))
        endif
      enddo

    ELSEIF(vadvordrv.eq.5)THEN

      do k=4,nk-2
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*( 2.*u3d(i,j,k-3)-13.*u3d(i,j,k-2)+47.*u3d(i,j,k-1)   &
                          +27.*u3d(i,j,k)-3.*u3d(i,j,k+1) )*onedsixty
        else
          dum(i,j,k)=wbar*( 2.*u3d(i,j,k+2)-13.*u3d(i,j,k+1)+47.*u3d(i,j,k)   &
                          +27.*u3d(i,j,k-1)-3.*u3d(i,j,k-2) )*onedsixty
        endif
      enddo
      enddo

      k = 3
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*(-u3d(i,j,k-2)+5.*u3d(i,j,k-1)+2.*u3d(i,j,k  ))*onedsix
        else
          dum(i,j,k)=wbar*( 2.*u3d(i,j,k+2)-13.*u3d(i,j,k+1)+47.*u3d(i,j,k)   &
                          +27.*u3d(i,j,k-1)-3.*u3d(i,j,k-2) )*onedsixty
        endif
      enddo

      k = nk-1
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*( 2.*u3d(i,j,k-3)-13.*u3d(i,j,k-2)+47.*u3d(i,j,k-1)   &
                          +27.*u3d(i,j,k)-3.*u3d(i,j,k+1) )*onedsixty
        else
          dum(i,j,k)=wbar*(-u3d(i,j,k+1)+5.*u3d(i,j,k  )+2.*u3d(i,j,k-1))*onedsix
        endif
      enddo

      k = 2
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*(c1(1,1,k)*u3d(i,j,k-1)+c2(1,1,k)*u3d(i,j,k))
        else
          dum(i,j,k)=wbar*(-u3d(i,j,k+1)+5.*u3d(i,j,k  )+2.*u3d(i,j,k-1))*onedsix
        endif
      enddo

      k = nk
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*(-u3d(i,j,k-2)+5.*u3d(i,j,k-1)+2.*u3d(i,j,k  ))*onedsix
        else
          dum(i,j,k)=wbar*(c1(1,1,k)*u3d(i,j,k-1)+c2(1,1,k)*u3d(i,j,k))
        endif
      enddo

    ELSEIF(vadvordrv.eq.6)THEN

      do k=4,nk-2
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        dum(i,j,k)=wbar*( 37.0*(u3d(i,j,k  )+u3d(i,j,k-1)) &
                          -8.0*(u3d(i,j,k+1)+u3d(i,j,k-2)) &
                              +(u3d(i,j,k+2)+u3d(i,j,k-3)) )*onedsixty
      enddo
      enddo

      do k=3,(nk-1),(nk-4)
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        dum(i,j,k)=wbar*( 7.0*(u3d(i,j,k  )+u3d(i,j,k-1)) &
                             -(u3d(i,j,k+1)+u3d(i,j,k-2)) )*onedtwelve
      enddo
      enddo

      do k=2,nk,(nk-2)
      do i=2,i2
        wbar = 0.5*(arf2(i)*rrw(i,j,k)+arf1(i)*rrw(i-1,j,k))
        dum(i,j,k)=wbar*(c1(1,1,k)*u3d(i,j,k-1)+c2(1,1,k)*u3d(i,j,k))
      enddo
      enddo

    ENDIF

!------

        k=1
        do i=2,i2
          advz(i,j,k)=-dum(i,j,k+1)*rdz*mh(1,1,k)
          uten(i,j,k)=uten(i,j,k)+( advx(i,j,k)+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(arf2(i)*divx(i,j,k)+arf1(i)*divx(i-1,j,k)) )*rr0(1,1,k)
        enddo
        do k=2,nk-1
        do i=2,i2
          advz(i,j,k)=-(dum(i,j,k+1)-dum(i,j,k))*rdz*mh(1,1,k)
          uten(i,j,k)=uten(i,j,k)+( advx(i,j,k)+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(arf2(i)*divx(i,j,k)+arf1(i)*divx(i-1,j,k)) )*rr0(1,1,k)
        enddo
        enddo
        k=nk
        do i=2,i2
          advz(i,j,k)=+dum(i,j,k)*rdz*mh(1,1,k)
          uten(i,j,k)=uten(i,j,k)+( advx(i,j,k)+advz(i,j,k)    &
                     +u3d(i,j,k)*0.5*(arf2(i)*divx(i,j,k)+arf1(i)*divx(i-1,j,k)) )*rr0(1,1,k)
        enddo

    ENDDO    jloopuasymm

  ELSE
    stop 55555
  ENDIF  vadvu

!  end vadvu for axisymmetric grid
!----------------------------------------------------------------

      if(timestats.ge.1) time_advu=time_advu+mytime()
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advv(nrk,xh,rxh,arh1,arh2,uh,xf,vf,gz,rgz,gzv,mh,rho0,rr0,rf0,rrf0,dum,advx,advy,advz,divx, &
                       rru,rrv,v3d,vten,rrw,rdsf,c1,c2,rho,dt)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      integer, intent(in) :: nrk
      real, dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je+1) :: vf
      real, dimension(itb:ite,jtb:jte) :: gz,rgz,gzv
      real, dimension(ib:ie,jb:je,kb:ke) :: mh,rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advx,advy,advz,divx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv,v3d,vten
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, intent(in), dimension(kb:ke+1) :: rdsf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho
      real, intent(in) :: dt
 
      integer i,j,k,i1,i2,j1,j2,jd1,jd2
      real :: ubar,vbar,wbar,cc1,cc2
      logical :: doitw,doite,doits,doitn

      logical :: doweno
      real :: dd,rr,phi
      real :: s1,s2,s3,s4,s5
      real :: f1,f2,f3
      real :: b1,b2,b3
      real :: w1,w2,w3
      double precision :: a1,a2,a3,a4
      double precision :: weps
      real :: coef

!------------------------------------------------------------

      doweno = .false.
      IF( (advwenov.eq.1) .or. (advwenov.eq.2.and.nrk.eq.3) ) doweno = .true.
      weps = 100.0*epsilon

      IF( idiff.ge.1 .and. difforder.eq.6 )THEN
        coef = kdiff6/64.0/dt
      ENDIF

!-----------------

      i1 = 1
      i2 = ni+1

      doitw = .false.
      doite = .false.

      IF(wbc.eq.2 .and. ibw.eq.1) doitw = .true.
      IF(ebc.eq.2 .and. ibe.eq.1) doite = .true.

!-----------------

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

      jd1 = j1-1
      jd2 = j2

      doits = .false.
      doitn = .false.

      IF(sbc.eq.2 .and. ibs.eq.1) doits = .true.
      IF(nbc.eq.2 .and. ibn.eq.1) doitn = .true.

!----------------------------------------------------------------

    hadvsection:  IF(axisymm.eq.1)THEN
      call advvaxi(doweno,xh,rxh,arh1,arh2,uh,xf,vf,rho0,rr0,rf0,rrf0,dum,advx,rru,v3d,vten)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=j1,j2
      do i=1,ni
        advy(i,j,k)=0.0
      enddo
      enddo
      enddo
    ELSE

!-----------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar)
    DO k=1,nk

! Advection in x-direction

    if(doweno)then
      do j=j1,j2
      do i=i1,i2
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        if(ubar.ge.0.0)then
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

        dum(i,j,k)=ubar*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo
    elseif(hadvordrv.eq.5)then
      do j=j1,j2
      do i=i1,i2
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        if(ubar.ge.0.)then
          dum(i,j,k)=ubar*( 2.*v3d(i-3,j,k)-13.*v3d(i-2,j,k)+47.*v3d(i-1,j,k)    &
                          +27.*v3d(i,j,k)-3.*v3d(i+1,j,k) )*onedsixty
        else
          dum(i,j,k)=ubar*( 2.*v3d(i+2,j,k)-13.*v3d(i+1,j,k)+47.*v3d(i,j,k)    &
                          +27.*v3d(i-1,j,k)-3.*v3d(i-2,j,k) )*onedsixty
        endif
      enddo
      enddo
    elseif(hadvordrv.eq.6)then
      do j=j1,j2
      do i=i1,i2
        ubar = 0.5*(rru(i,j,k)+rru(i,j-1,k))
        dum(i,j,k)=ubar*( 37.0*(v3d(i  ,j,k)+v3d(i-1,j,k))     &
                          -8.0*(v3d(i+1,j,k)+v3d(i-2,j,k))     &
                              +(v3d(i+2,j,k)+v3d(i-3,j,k)) )*onedsixty
      enddo
      enddo
    endif

      if(doitw)then
      do j=j1,j2
        i=1
        if((rru(i,j,k)+rru(i,j-1,k)).ge.0.0)then
          dum(i,j,k)=dum(i+1,j,k)
        endif
      enddo
      endif

      if(doite)then
      do j=j1,j2
        i=ni+1
        if((rru(i,j,k)+rru(i,j-1,k)).le.0.0)then
          dum(i,j,k)=dum(i-1,j,k)
        endif
      enddo
      endif

      do j=j1,j2
      do i=1,ni
        advx(i,j,k)=-(dum(i+1,j,k)-dum(i,j,k))*rdx*uh(i)
      enddo
      enddo

      IF(doitw)THEN
        do j=j1,j2
          if((rru(1,j,k)+rru(1,j-1,k)).ge.0.0)then
            i=1
            advx(i,j,k)=advx(i,j,k)-v3d(i,j,k)*0.5*(            &
                    (rru(i+1,j-1,k)-rru(i,j-1,k))               &
                   +(rru(i+1,j  ,k)-rru(i,j  ,k)) )*rdx*uh(i)
          endif
        enddo
      ENDIF

      IF(doite)THEN
        do j=j1,j2
          if((rru(ni+1,j,k)+rru(ni+1,j-1,k)).le.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-v3d(i,j,k)*0.5*(            &
                    (rru(i+1,j-1,k)-rru(i,j-1,k))               &
                   +(rru(i+1,j  ,k)-rru(i,j  ,k)) )*rdx*uh(i)
          endif
        enddo
      ENDIF

    !-------------------------------------------------------
    ! 6th-order diffusion:
    IF( idiff.ge.1 .and. difforder.eq.6 )THEN
      do j=1,nj+1
      do i=1,ni+1
        dum(i,j,k)=( 10.0*(v3d(i  ,j,k)-v3d(i-1,j,k))     &
                     -5.0*(v3d(i+1,j,k)-v3d(i-2,j,k))     &
                         +(v3d(i+2,j,k)-v3d(i-3,j,k)) )   &
                  *0.25*( (rho(i-1,j-1,k)+rho(i,j,k))     &
                         +(rho(i-1,j,k)+rho(i,j-1,k)) )
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj+1
        do i=1,ni+1
          if( dum(i,j,k)*(v3d(i,j,k)-v3d(i-1,j,k)).le.0.0 )then
            dum(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj+1
      do i=1,ni
        advx(i,j,k)=advx(i,j,k)+coef*(dum(i+1,j,k)-dum(i,j,k))*(rho0(i,j-1,k)+rho0(i,j,k))/(rho(i,j-1,k)+rho(i,j,k))
      enddo
      enddo
    ENDIF
    !-------------------------------------------------------

! Advection in y-direction

    if(doweno)then
      do j=jd1,jd2
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        if(vbar.ge.0.0)then
          s1=v3d(i,j-2,k)
          s2=v3d(i,j-1,k)
          s3=v3d(i,j  ,k)
          s4=v3d(i,j+1,k)
          s5=v3d(i,j+2,k)
        else
          s1=v3d(i,j+3,k)
          s2=v3d(i,j+2,k)
          s3=v3d(i,j+1,k)
          s4=v3d(i,j  ,k)
          s5=v3d(i,j-1,k)
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

        dum(i,j,k)=vbar*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo
    elseif(hadvordrv.eq.5)then
      do j=jd1,jd2
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        if(vbar.ge.0.)then
          dum(i,j,k)=vbar*( 2.*v3d(i,j-2,k)-13.*v3d(i,j-1,k)+47.*v3d(i,j,k)    &
                          +27.*v3d(i,j+1,k)-3.*v3d(i,j+2,k) )*onedsixty
        else
          dum(i,j,k)=vbar*( 2.*v3d(i,j+3,k)-13.*v3d(i,j+2,k)+47.*v3d(i,j+1,k)    &
                          +27.*v3d(i,j,k)-3.*v3d(i,j-1,k) )*onedsixty
        endif
      enddo
      enddo
    elseif(hadvordrv.eq.6)then
      do j=jd1,jd2
      do i=1,ni
        vbar = 0.5*(rrv(i,j,k)+rrv(i,j+1,k))
        dum(i,j,k)=vbar*( 37.0*(v3d(i,j+1,k)+v3d(i,j  ,k)) &
                          -8.0*(v3d(i,j+2,k)+v3d(i,j-1,k)) &
                              +(v3d(i,j+3,k)+v3d(i,j-2,k)) )*onedsixty
      enddo
      enddo
    endif

      do j=j1,j2
      do i=1,ni
        advy(i,j,k)=-(dum(i,j,k)-dum(i,j-1,k))*rdy*vf(j)
      enddo
      enddo

    !-------------------------------------------------------
    ! 6th-order diffusion:
    IF( idiff.ge.1 .and. difforder.eq.6 )THEN
      do j=1,nj+2
      do i=1,ni
        dum(i,j,k)=( 10.0*(v3d(i,j  ,k)-v3d(i,j-1,k))     &
                     -5.0*(v3d(i,j+1,k)-v3d(i,j-2,k))     &
                         +(v3d(i,j+2,k)-v3d(i,j-3,k)) )*rho(i,j-1,k)
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj+2
        do i=1,ni
          if( dum(i,j,k)*(v3d(i,j,k)-v3d(i,j-1,k)).le.0.0 )then
            dum(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj+1
      do i=1,ni
        advy(i,j,k)=advy(i,j,k)+coef*(dum(i,j+1,k)-dum(i,j,k))*(rho0(i,j-1,k)+rho0(i,j,k))/(rho(i,j-1,k)+rho(i,j,k))
      enddo
      enddo
    ENDIF
    !-------------------------------------------------------

    ENDDO

    ENDIF  hadvsection

!----------------------------------------------------------------
! Advection in z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar,cc1,cc2)
  jloopv:  DO j=j1,j2

    IF(doweno)THEN

      do k=3,nk-1
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.0)then
          s1=v3d(i,j,k-3)
          s2=v3d(i,j,k-2)
          s3=v3d(i,j,k-1)
          s4=v3d(i,j,k  )
          s5=v3d(i,j,k+1)
        else
          s1=v3d(i,j,k+2)
          s2=v3d(i,j,k+1)
          s3=v3d(i,j,k  )
          s4=v3d(i,j,k-1)
          s5=v3d(i,j,k-2)
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

        dum(i,j,k)=wbar*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo

      ! flux at k=3 if w > 0
      k = 3
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.gt.0.0 )then
          dd = v3d(i,j,k-1)-v3d(i,j,k-2)
          rr = (v3d(i,j,k)-v3d(i,j,k-1))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( v3d(i,j,k-1) + 0.5*phi*(v3d(i,j,k-1)-v3d(i,j,k-2)) )
        endif
      enddo

      ! flux at k=(nk-1) if w < 0
      k = nk-1
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.lt.0.0 )then
          dd = v3d(i,j,k)-v3d(i,j,k+1)
          rr = (v3d(i,j,k-1)-v3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( v3d(i,j,k) + 0.5*phi*(v3d(i,j,k)-v3d(i,j,k+1)) )
        endif
      enddo

      k = 2
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.ge.0.0 )then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dum(i,j,k)=wbar*(cc1*v3d(i,j,k-1)+cc2*v3d(i,j,k))
        else
          dd = v3d(i,j,k)-v3d(i,j,k+1)
          rr = (v3d(i,j,k-1)-v3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( v3d(i,j,k) + 0.5*phi*(v3d(i,j,k)-v3d(i,j,k+1)) )
        endif
      enddo

      k = nk
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if( wbar.gt.0.0 )then
          dd = v3d(i,j,k-1)-v3d(i,j,k-2)
          rr = (v3d(i,j,k)-v3d(i,j,k-1))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( v3d(i,j,k-1) + 0.5*phi*(v3d(i,j,k-1)-v3d(i,j,k-2)) )
        else
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dum(i,j,k)=wbar*(cc1*v3d(i,j,k-1)+cc2*v3d(i,j,k))
        endif
      enddo

    ELSEIF(vadvordrv.eq.5)THEN

      do k=4,nk-2
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*( 2.*v3d(i,j,k-3)-13.*v3d(i,j,k-2)+47.*v3d(i,j,k-1)    &
                          +27.*v3d(i,j,k)-3.*v3d(i,j,k+1) )*onedsixty
        else
          dum(i,j,k)=wbar*( 2.*v3d(i,j,k+2)-13.*v3d(i,j,k+1)+47.*v3d(i,j,k)    &
                          +27.*v3d(i,j,k-1)-3.*v3d(i,j,k-2) )*onedsixty
        endif
      enddo
      enddo

      k = 3
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*(-v3d(i,j,k-2)+5.*v3d(i,j,k-1)+2.*v3d(i,j,k  ))*onedsix
        else
          dum(i,j,k)=wbar*( 2.*v3d(i,j,k+2)-13.*v3d(i,j,k+1)+47.*v3d(i,j,k)    &
                          +27.*v3d(i,j,k-1)-3.*v3d(i,j,k-2) )*onedsixty
        endif
      enddo

      k = nk-1
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*( 2.*v3d(i,j,k-3)-13.*v3d(i,j,k-2)+47.*v3d(i,j,k-1)    &
                          +27.*v3d(i,j,k)-3.*v3d(i,j,k+1) )*onedsixty
        else
          dum(i,j,k)=wbar*(-v3d(i,j,k+1)+5.*v3d(i,j,k  )+2.*v3d(i,j,k-1))*onedsix
        endif
      enddo

      k = 2
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.)then
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dum(i,j,k)=wbar*(cc1*v3d(i,j,k-1)+cc2*v3d(i,j,k))
        else
          dum(i,j,k)=wbar*(-v3d(i,j,k+1)+5.*v3d(i,j,k  )+2.*v3d(i,j,k-1))*onedsix
        endif
      enddo

      k = nk
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*(-v3d(i,j,k-2)+5.*v3d(i,j,k-1)+2.*v3d(i,j,k  ))*onedsix
        else
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          dum(i,j,k)=wbar*(cc1*v3d(i,j,k-1)+cc2*v3d(i,j,k))
        endif
      enddo

    ELSEIF(vadvordrv.eq.6)THEN

      do k=4,nk-2
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        dum(i,j,k)=wbar*( 37.0*(v3d(i,j,k  )+v3d(i,j,k-1)) &
                          -8.0*(v3d(i,j,k+1)+v3d(i,j,k-2)) &
                              +(v3d(i,j,k+2)+v3d(i,j,k-3)) )*onedsixty
      enddo
      enddo

      do k=3,(nk-1),(nk-4)
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        dum(i,j,k)=wbar*( 7.0*(v3d(i,j,k  )+v3d(i,j,k-1)) &
                             -(v3d(i,j,k+1)+v3d(i,j,k-2)) )*onedtwelve
      enddo
      enddo

      do k=2,nk,(nk-2)
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j-1,k))
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        dum(i,j,k)=wbar*(cc1*v3d(i,j,k-1)+cc2*v3d(i,j,k))
      enddo
      enddo

    ENDIF

!------

      IF(terrain_flag)THEN

        k=1
        do i=1,ni
          advz(i,j,k)=-dum(i,j,k+1)*rdsf(k)
          vten(i,j,k)=vten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +v3d(i,j,k)*0.5*(divx(i,j,k)+divx(i,j-1,k)) )         &
                  *gzv(i,j)/(0.5*(rho0(i,j-1,k)+rho0(i,j,k)))
        enddo
        do k=2,nk-1
        do i=1,ni
          advz(i,j,k)=-(dum(i,j,k+1)-dum(i,j,k))*rdsf(k)
          vten(i,j,k)=vten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +v3d(i,j,k)*0.5*(divx(i,j,k)+divx(i,j-1,k)) )         &
                  *gzv(i,j)/(0.5*(rho0(i,j-1,k)+rho0(i,j,k)))
        enddo
        enddo
        k=nk
        do i=1,ni
          advz(i,j,k)=+dum(i,j,k)*rdsf(k)
          vten(i,j,k)=vten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +v3d(i,j,k)*0.5*(divx(i,j,k)+divx(i,j-1,k)) )         &
                  *gzv(i,j)/(0.5*(rho0(i,j-1,k)+rho0(i,j,k)))
        enddo

      ELSE

        k=1
        do i=1,ni
          advz(i,j,k)=-dum(i,j,k+1)*rdz*mh(1,1,k)
          vten(i,j,k)=vten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +v3d(i,j,k)*0.5*(divx(i,j,k)+divx(i,j-1,k)) )*rr0(1,1,k)
        enddo
        do k=2,nk-1
        do i=1,ni
          advz(i,j,k)=-(dum(i,j,k+1)-dum(i,j,k))*rdz*mh(1,1,k)
          vten(i,j,k)=vten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +v3d(i,j,k)*0.5*(divx(i,j,k)+divx(i,j-1,k)) )*rr0(1,1,k)
        enddo
        enddo
        k=nk
        do i=1,ni
          advz(i,j,k)=+dum(i,j,k)*rdz*mh(1,1,k)
          vten(i,j,k)=vten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                     +v3d(i,j,k)*0.5*(divx(i,j,k)+divx(i,j-1,k)) )*rr0(1,1,k)
        enddo

      ENDIF

    ENDDO  jloopv
!----------------------------------------------------------------

      if(timestats.ge.1) time_advv=time_advv+mytime()
 
      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine advw(nrk,xh,rxh,arh1,arh2,uh,xf,vh,gz,rgz,mf,rho0,rr0,rf0,rrf0,dum,advx,advy,advz,divx, &
                       rru,rrv,rrw,w3d,wten,rds,c1,c2,rho,dt)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      integer, intent(in) :: nrk
      real, dimension(ib:ie) :: xh,rxh,arh1,arh2,uh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: vh
      real, dimension(itb:ite,jtb:jte) :: gz,rgz
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: dum,advx,advy,advz,divx
      real, dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, dimension(ib:ie,jb:je,kb:ke+1) :: rrw,w3d,wten
      real, intent(in), dimension(kb:ke+1) :: rds
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho
      real, intent(in) :: dt
 
      integer i,j,k,i1,i2,j1,j2
      real :: ubar,vbar,wbar,cc1,cc2
      logical :: doitw,doite,doits,doitn

      logical :: doweno
      real :: dd,rr,phi
      real :: s1,s2,s3,s4,s5
      real :: f1,f2,f3
      real :: b1,b2,b3
      real :: w1,w2,w3
      double precision :: a1,a2,a3,a4
      double precision :: weps
      real :: coef

!----------------------------------------------------------------

      doweno = .false.
      IF( (advwenov.eq.1) .or. (advwenov.eq.2.and.nrk.eq.3) ) doweno = .true.
      weps = 100.0*epsilon

      IF( idiff.ge.1 .and. difforder.eq.6 )THEN
        coef = kdiff6/64.0/dt
      ENDIF

!-----------------

      i1 = 1
      i2 = ni+1

      doitw = .false.
      doite = .false.

      IF(wbc.eq.2 .and. ibw.eq.1) doitw = .true.
      IF(ebc.eq.2 .and. ibe.eq.1) doite = .true.

!-----------------

      j1 = 1
      j2 = nj+1

      doits = .false.
      doitn = .false.

      IF(sbc.eq.2 .and. ibs.eq.1) doits = .true.
      IF(nbc.eq.2 .and. ibn.eq.1) doitn = .true.

!----------------------------------------------------------------

    hadvsection:  IF(axisymm.eq.1)THEN
      call advwaxi(doweno,xh,rxh,arh1,arh2,uh,xf,vh,rho0,rr0,rf0,rrf0,dum,advx,rru,w3d,wten,c1,c2)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        advy(i,j,k)=0.0
      enddo
      enddo
      enddo
    ELSE

!-----------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k,ubar,vbar,cc1,cc2)
    DO k=2,nk

! Advection in x-direction

    if(doweno)then
      do j=1,nj
      do i=i1,i2
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
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
    elseif(hadvordrv.eq.5)then
      do j=1,nj
      do i=i1,i2
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        if(ubar.ge.0.)then
          dum(i,j,k)=ubar*( 2.*w3d(i-3,j,k)-13.*w3d(i-2,j,k)+47.*w3d(i-1,j,k)    &
                          +27.*w3d(i,j,k)-3.*w3d(i+1,j,k) )*onedsixty
        else
          dum(i,j,k)=ubar*( 2.*w3d(i+2,j,k)-13.*w3d(i+1,j,k)+47.*w3d(i,j,k)    &
                          +27.*w3d(i-1,j,k)-3.*w3d(i-2,j,k) )*onedsixty
        endif
      enddo
      enddo
    elseif(hadvordrv.eq.6)then
      do j=1,nj
      do i=i1,i2
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        dum(i,j,k)=ubar*( 37.0*(w3d(i  ,j,k)+w3d(i-1,j,k)) &
                          -8.0*(w3d(i+1,j,k)+w3d(i-2,j,k)) &
                              +(w3d(i+2,j,k)+w3d(i-3,j,k)) )*onedsixty
      enddo
      enddo
    endif

      if(doitw)then
      do j=1,nj
        i=1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        if(ubar.ge.0.0)then
          dum(i,j,k)=dum(i+1,j,k)
        endif
      enddo
      endif

      if(doite)then
      do j=1,nj
        i=ni+1
        cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
        if(ubar.le.0.0)then
          dum(i,j,k)=dum(i-1,j,k)
        endif
      enddo
      endif

      do j=1,nj
      do i=1,ni
        advx(i,j,k)=-(dum(i+1,j,k)-dum(i,j,k))*rdx*uh(i)
      enddo
      enddo

      IF(doitw)THEN
        do j=1,nj
          i=1
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
          if(ubar.ge.0.0)then
            i=1
            advx(i,j,k)=advx(i,j,k)-w3d(i,j,k)*(                    &
                    c1(i,j,k)*(rru(i+1,j,k-1)-rru(i,j,k-1))         &
                   +c2(i,j,k)*(rru(i+1,j,k  )-rru(i,j,k  )) )*rdx*uh(i)
          endif
        enddo
      ENDIF

      IF(doite)THEN
        do j=1,nj
          i=ni+1
          cc2 = 0.5*(c2(i-1,j,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          ubar = cc2*rru(i,j,k)+cc1*rru(i,j,k-1)
          if(ubar.le.0.0)then
            i=ni
            advx(i,j,k)=advx(i,j,k)-w3d(i,j,k)*(                    &
                    c1(i,j,k)*(rru(i+1,j,k-1)-rru(i,j,k-1))         &
                   +c2(i,j,k)*(rru(i+1,j,k  )-rru(i,j,k  )) )*rdx*uh(i)
          endif
        enddo
      ENDIF

    !-------------------------------------------------------
    ! 6th-order diffusion:
    IF( idiff.ge.1 .and. difforder.eq.6 )THEN
      do j=1,nj
      do i=1,ni+1
        dum(i,j,k)=( 10.0*(w3d(i  ,j,k)-w3d(i-1,j,k))     &
                     -5.0*(w3d(i+1,j,k)-w3d(i-2,j,k))     &
                         +(w3d(i+2,j,k)-w3d(i-3,j,k)) )   &
                  *0.25*( (rho(i-1,j,k  )+rho(i,j,k  ))   &
                         +(rho(i-1,j,k-1)+rho(i,j,k-1)) )
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj
        do i=1,ni+1
          if( dum(i,j,k)*(w3d(i,j,k)-w3d(i-1,j,k)).le.0.0 )then
            dum(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj
      do i=1,ni
        advx(i,j,k)=advx(i,j,k)+coef*(dum(i+1,j,k)-dum(i,j,k))*rf0(i,j,k)/(0.5*(rho(i,j,k-1)+rho(i,j,k)))
      enddo
      enddo
    ENDIF
    !-------------------------------------------------------

! Advection in y-direction

    if(doweno)then
      do j=j1,j2
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        if(vbar.ge.0.0)then
          s1=w3d(i,j-3,k)
          s2=w3d(i,j-2,k)
          s3=w3d(i,j-1,k)
          s4=w3d(i,j  ,k)
          s5=w3d(i,j+1,k)
        else
          s1=w3d(i,j+2,k)
          s2=w3d(i,j+1,k)
          s3=w3d(i,j  ,k)
          s4=w3d(i,j-1,k)
          s5=w3d(i,j-2,k)
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

        dum(i,j,k)=vbar*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo
    elseif(hadvordrv.eq.5)then
      do j=j1,j2
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        if(vbar.ge.0.)then
          dum(i,j,k)=vbar*( 2.*w3d(i,j-3,k)-13.*w3d(i,j-2,k)+47.*w3d(i,j-1,k)    &
                          +27.*w3d(i,j,k)-3.*w3d(i,j+1,k) )*onedsixty
        else
          dum(i,j,k)=vbar*( 2.*w3d(i,j+2,k)-13.*w3d(i,j+1,k)+47.*w3d(i,j,k)    &
                          +27.*w3d(i,j-1,k)-3.*w3d(i,j-2,k) )*onedsixty
        endif
      enddo
      enddo
    elseif(hadvordrv.eq.6)then
      do j=j1,j2
      do i=1,ni
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        dum(i,j,k)=vbar*( 37.0*(w3d(i,j  ,k)+w3d(i,j-1,k)) &
                          -8.0*(w3d(i,j+1,k)+w3d(i,j-2,k)) &
                              +(w3d(i,j+2,k)+w3d(i,j-3,k)) )*onedsixty
      enddo
      enddo
    endif

      if(doits)then
      do i=1,ni
        j=1
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        if(vbar.ge.0.0)then
          dum(i,j,k)=dum(i,j+1,k)
        endif
      enddo
      endif

      if(doitn)then
      do i=1,ni
        j=nj+1
        cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
        cc1 = 1.0-cc2
        vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
        if(vbar.le.0.0)then
          dum(i,j,k)=dum(i,j-1,k)
        endif
      enddo
      endif

      do j=1,nj
      do i=1,ni
        advy(i,j,k)=-(dum(i,j+1,k)-dum(i,j,k))*rdy*vh(j)
      enddo
      enddo

      IF(doits)THEN
        do i=1,ni
          j=1
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
          if(vbar.ge.0.0)then
            j=1
            advy(i,j,k)=advy(i,j,k)-w3d(i,j,k)*(                     &
                           c1(i,j,k)*(rrv(i,j+1,k-1)-rrv(i,j,k-1))   &
                          +c2(i,j,k)*(rrv(i,j+1,k  )-rrv(i,j,k  )) )*rdy*vh(j)
          endif
        enddo
      ENDIF

      IF(doitn)THEN
        do i=1,ni
          j=nj+1
          cc2 = 0.5*(c2(i,j-1,k)+c2(i,j,k))
          cc1 = 1.0-cc2
          vbar = cc2*rrv(i,j,k)+cc1*rrv(i,j,k-1)
          if(vbar.le.0.0)then
            j=nj
            advy(i,j,k)=advy(i,j,k)-w3d(i,j,k)*(                     &
                           c1(i,j,k)*(rrv(i,j+1,k-1)-rrv(i,j,k-1))   &
                          +c2(i,j,k)*(rrv(i,j+1,k  )-rrv(i,j,k  )) )*rdy*vh(j)
          endif
        enddo
      ENDIF

    !-------------------------------------------------------
    ! 6th-order diffusion:
    IF( idiff.ge.1 .and. difforder.eq.6 )THEN
      do j=1,nj+1
      do i=1,ni
        dum(i,j,k)=( 10.0*(w3d(i,j  ,k)-w3d(i,j-1,k))     &
                     -5.0*(w3d(i,j+1,k)-w3d(i,j-2,k))     &
                         +(w3d(i,j+2,k)-w3d(i,j-3,k)) )   &
                  *0.25*( (rho(i,j-1,k  )+rho(i,j,k  ))   &
                         +(rho(i,j-1,k-1)+rho(i,j,k-1)) )
      enddo
      enddo
      if(mdiff.eq.1)then
        do j=1,nj+1
        do i=1,ni
          if( dum(i,j,k)*(w3d(i,j,k)-w3d(i,j-1,k)).le.0.0 )then
            dum(i,j,k)=0.0
          endif
        enddo
        enddo
      endif
      do j=1,nj
      do i=1,ni
        advy(i,j,k)=advy(i,j,k)+coef*(dum(i,j+1,k)-dum(i,j,k))*rf0(i,j,k)/(0.5*(rho(i,j,k-1)+rho(i,j,k)))
      enddo
      enddo
    ENDIF
    !-------------------------------------------------------

    ENDDO

    ENDIF  hadvsection

!----------------------------------------------------------------
! Advection in z-direction

!$omp parallel do default(shared)   &
!$omp private(i,j,k,wbar)
  jloopw:  DO j=1,nj

    IF(doweno)THEN

      do k=2,nk-1
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.0)then
          s1=w3d(i,j,k-2)
          s2=w3d(i,j,k-1)
          s3=w3d(i,j,k  )
          s4=w3d(i,j,k+1)
          s5=w3d(i,j,k+2)
        else
          s1=w3d(i,j,k+3)
          s2=w3d(i,j,k+2)
          s3=w3d(i,j,k+1)
          s4=w3d(i,j,k  )
          s5=w3d(i,j,k-1)
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

        dum(i,j,k)=wbar*((w1*f1)+(w2*f2)+(w3*f3))/(w1+w2+w3)
      enddo
      enddo

      ! flux at k=2 if w > 0
      k = 2
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.gt.0.0 )then
          dd = w3d(i,j,k)-w3d(i,j,k-1)
          rr = (w3d(i,j,k+1)-w3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( w3d(i,j,k) + 0.5*phi*(w3d(i,j,k)-w3d(i,j,k-1)) )
        endif
      enddo

      ! flux at k=(nk-1) if w < 0
      k = nk-1
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.lt.0.0 )then
          dd = w3d(i,j,k+1)-w3d(i,j,k+2)
          rr = (w3d(i,j,k)-w3d(i,j,k+1))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( w3d(i,j,k+1) + 0.5*phi*(w3d(i,j,k+1)-w3d(i,j,k+2)) )
        endif
      enddo

      k = 1
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.ge.0.0 )then
          dum(i,j,k)=wbar*0.5*(w3d(i,j,k)+w3d(i,j,k+1))
        else
          dd = w3d(i,j,k+1)-w3d(i,j,k+2)
          rr = (w3d(i,j,k)-w3d(i,j,k+1))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( w3d(i,j,k+1) + 0.5*phi*(w3d(i,j,k+1)-w3d(i,j,k+2)) )
        endif
      enddo

      k = nk
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if( wbar.gt.0.0 )then
          dd = w3d(i,j,k)-w3d(i,j,k-1)
          rr = (w3d(i,j,k+1)-w3d(i,j,k))/(sign(sngl(weps),dd)+dd)
          phi = max(0.0,min(2.0*rr,min( onedthree+twodthree*rr , 2.0 ) ) )
          dum(i,j,k) = wbar*( w3d(i,j,k) + 0.5*phi*(w3d(i,j,k)-w3d(i,j,k-1)) )
        else
          dum(i,j,k)=wbar*0.5*(w3d(i,j,k)+w3d(i,j,k+1))
        endif
      enddo

    ELSEIF(vadvordrv.eq.5)THEN

      do k=3,nk-2
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*( 2.*w3d(i,j,k-2)-13.*w3d(i,j,k-1)+47.*w3d(i,j,k)    &
                          +27.*w3d(i,j,k+1)-3.*w3d(i,j,k+2) )*onedsixty
        else
          dum(i,j,k)=wbar*( 2.*w3d(i,j,k+3)-13.*w3d(i,j,k+2)+47.*w3d(i,j,k+1)  &
                          +27.*w3d(i,j,k)-3.*w3d(i,j,k-1) )*onedsixty
        endif
      enddo
      enddo

      k = 2
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*(-w3d(i,j,k-1)+5.*w3d(i,j,k  )+2.*w3d(i,j,k+1))*onedsix
        else
          dum(i,j,k)=wbar*( 2.*w3d(i,j,k+3)-13.*w3d(i,j,k+2)+47.*w3d(i,j,k+1)  &
                          +27.*w3d(i,j,k)-3.*w3d(i,j,k-1) )*onedsixty
        endif
      enddo

      k = nk-1
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*( 2.*w3d(i,j,k-2)-13.*w3d(i,j,k-1)+47.*w3d(i,j,k)    &
                          +27.*w3d(i,j,k+1)-3.*w3d(i,j,k+2) )*onedsixty
        else
          dum(i,j,k)=wbar*(-w3d(i,j,k+2)+5.*w3d(i,j,k+1)+2.*w3d(i,j,k  ))*onedsix
        endif
      enddo

      k = 1
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*0.5*(w3d(i,j,k)+w3d(i,j,k+1))
        else
          dum(i,j,k)=wbar*(-w3d(i,j,k+2)+5.*w3d(i,j,k+1)+2.*w3d(i,j,k  ))*onedsix
        endif
      enddo

      k = nk
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        if(wbar.ge.0.)then
          dum(i,j,k)=wbar*(-w3d(i,j,k-1)+5.*w3d(i,j,k  )+2.*w3d(i,j,k+1))*onedsix
        else
          dum(i,j,k)=wbar*0.5*(w3d(i,j,k)+w3d(i,j,k+1))
        endif
      enddo

    ELSEIF(vadvordrv.eq.6)THEN

      do k=3,nk-2
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        dum(i,j,k)=wbar*( 37.0*(w3d(i,j,k+1)+w3d(i,j,k  )) &
                          -8.0*(w3d(i,j,k+2)+w3d(i,j,k-1)) &
                              +(w3d(i,j,k+3)+w3d(i,j,k-2)) )*onedsixty
      enddo
      enddo

      do k=2,(nk-1),(nk-3)
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        dum(i,j,k)=wbar*( 7.0*(w3d(i,j,k+1)+w3d(i,j,k  )) &
                             -(w3d(i,j,k+2)+w3d(i,j,k-1)) )*onedtwelve
      enddo
      enddo

      do k=1,nk,(nk-1)
      do i=1,ni
        wbar = 0.5*(rrw(i,j,k)+rrw(i,j,k+1))
        dum(i,j,k)=wbar*0.5*(w3d(i,j,k)+w3d(i,j,k+1))
      enddo
      enddo

    ENDIF

!------

      IF(terrain_flag)THEN

      do k=2,nk
      do i=1,ni
        advz(i,j,k)=-(dum(i,j,k)-dum(i,j,k-1))*rds(k)
        wten(i,j,k)=wten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                   +w3d(i,j,k)*(c2(i,j,k)*divx(i,j,k)+c1(i,j,k)*divx(i,j,k-1)) )*rrf0(i,j,k)*gz(i,j)
      enddo
      enddo

      ELSE

      do k=2,nk
      do i=1,ni
        advz(i,j,k)=-(dum(i,j,k)-dum(i,j,k-1))*rdz*mf(1,1,k)
        wten(i,j,k)=wten(i,j,k)+( (advx(i,j,k)+advy(i,j,k))+advz(i,j,k)    &
                   +w3d(i,j,k)*(c2(1,1,k)*divx(i,j,k)+c1(1,1,k)*divx(i,j,k-1)) )*rrf0(1,1,k)
      enddo
      enddo

      ENDIF

    ENDDO  jloopw

!----------------------------------------------------------------

      if(timestats.ge.1) time_advw=time_advw+mytime()
 
      return
      end


