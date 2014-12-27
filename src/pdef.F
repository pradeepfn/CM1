!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefx(xh,arh1,arh2,uh,rho0,gz,rgz,rru,advx,dum,mass,s0,s,dt,flag,west,newwest,east,neweast)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(ib:ie) :: xh,arh1,arh2,uh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho0,s0,s
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: rru
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: advx,dum,mass
      real, intent(in) :: dt
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(cmp,jmp,kmp) :: west,newwest,east,neweast

      integer i,j,k
      real foo1,foo2,foo3,rdt

#ifdef MPI
      include 'mpif.h'
      integer, dimension(4) :: reqs_s
#endif

!----------------------------------------------------------------
! cm1r17:  include divx component

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    kloop:  DO k=1,nk

      IF(.not.terrain_flag)THEN
        IF(axisymm.eq.0)THEN
          do j=1,nj
          do i=1,ni
            dum(i,j,k)=rho0(1,1,k)*s0(i,j,k)+dt*( advx(i,j,k)  &
                   +s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i) )
          enddo
          enddo
        ELSE
          do j=1,nj
          do i=1,ni
            dum(i,j,k)=rho0(1,1,k)*s0(i,j,k)+dt*( advx(i,j,k)  &
                   +s(i,j,k)*(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i) )
          enddo
          enddo
        ENDIF
      ELSE
          do j=1,nj
          do i=1,ni
            dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)+dt*( advx(i,j,k)  &
                   +s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i) )*gz(i,j)
          enddo
          enddo
      ENDIF

        if(wbc.eq.2 .and. ibw.eq.1)then
          do j=1,nj
            if(rru(1,j,k).ge.0.0)then
              i=1
              dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)
            endif
          enddo
        endif

        if(ebc.eq.2 .and. ibe.eq.1)then
          do j=1,nj
            if(rru(ni+1,j,k).le.0.0)then
              i=ni
              dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)
            endif
          enddo
        endif

    ENDDO  kloop

        if(timestats.ge.1) time_pdef=time_pdef+mytime()

        call bcs(dum)
#ifdef MPI
        call comm_2we_start(dum,west(1,1,1),newwest(1,1,1),east(1,1,1),neweast(1,1,1),reqs_s)
        call comm_2we_end(dum,west(1,1,1),newwest(1,1,1),east(1,1,1),neweast(1,1,1),reqs_s)
#endif

        rdt=1.0/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k,foo1,foo2,foo3)
    DO k=1,nk

        do j=1,nj
        do i=-1,ni+2
          mass(i,j,k)=0.0
        enddo
        enddo

      IF(axisymm.eq.0)THEN

        jloop:  do j=1,nj
        do i=-1,ni+2
          flag(i,j,k)=.false.
        enddo
        do i=0,ni+1
          if(dum(i,j,k).lt.0.0)then
            foo1=max(0.0,dum(i-1,j,k))
            foo2=max(0.0,dum(i+1,j,k))
            if(foo1+foo2.gt.1.0e-30)then
              foo3=max(dum(i,j,k),-(foo1+foo2))/(foo1+foo2)
              mass(i-1,j,k)=mass(i-1,j,k)+foo1*foo3
              mass(i  ,j,k)=mass(i  ,j,k)-(foo1+foo2)*foo3
              mass(i+1,j,k)=mass(i+1,j,k)+foo2*foo3
              if(dum(i-1,j,k).gt.1.0e-30) flag(i-1,j,k)=.true.
                                          flag(i  ,j,k)=.true.
              if(dum(i+1,j,k).gt.1.0e-30) flag(i+1,j,k)=.true.
            endif
          endif
        enddo
        !-----
        IF(.not.terrain_flag)THEN
        do i=1,ni
        if(flag(i,j,k))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advx(i,j,k)=(dum(i,j,k)-rho0(1,1,k)*s0(i,j,k))*rdt  &
                   -s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)
        endif
        enddo
        !-----
        ELSE
        do i=1,ni
        if(flag(i,j,k))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advx(i,j,k)=(dum(i,j,k)-rho0(i,j,k)*s0(i,j,k))*rdt*rgz(i,j)  &
                   -s(i,j,k)*(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)
        endif
        enddo
        ENDIF
        !-----
        enddo  jloop

      ELSE

        do j=1,nj
        do i=-1,ni+2
          flag(i,j,k)=.false.
        enddo
        do i=0,ni+1
          if(dum(i,j,k).lt.0.0)then
            foo1=max(0.0,dum(i-1,j,k))
            foo2=max(0.0,dum(i+1,j,k))
            if(foo1+foo2.gt.1.0e-30)then
              foo3=max(xh(i)*dum(i,j,k),-(xh(i-1)*foo1+xh(i+1)*foo2))   &
                                        /(xh(i-1)*foo1+xh(i+1)*foo2)
              mass(i-1,j,k)=mass(i-1,j,k)+foo1*foo3
              mass(i  ,j,k)=mass(i  ,j,k)-(foo1+foo2)*foo3
              mass(i+1,j,k)=mass(i+1,j,k)+foo2*foo3
              if(dum(i-1,j,k).gt.1.0e-30) flag(i-1,j,k)=.true.
                                          flag(i  ,j,k)=.true.
              if(dum(i+1,j,k).gt.1.0e-30) flag(i+1,j,k)=.true.
            endif
          endif
        enddo
        do i=1,ni
        if(flag(i,j,k))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advx(i,j,k)=(dum(i,j,k)-rho0(1,1,k)*s0(i,j,k))*rdt  &
                   -s(i,j,k)*(arh2(i)*rru(i+1,j,k)-arh1(i)*rru(i,j,k))*rdx*uh(i)
        endif
        enddo
        enddo

      ENDIF

    ENDDO

!----------------------------------------------------------------

      if(timestats.ge.1) time_pdef=time_pdef+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefy(vh,rho0,gz,rgz,rrv,advy,dum,mass,s0,s,dt,flag,south,newsouth,north,newnorth)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho0,s0,s
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: rrv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: advy,dum,mass
      real, intent(in) :: dt
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(imp,cmp,kmp) :: south,newsouth,north,newnorth

      integer i,j,k
      real foo1,foo2,foo3,rdt

#ifdef MPI
      include 'mpif.h'
      integer, dimension(4) :: reqs_s
#endif

!----------------------------------------------------------------
! cm1r17:  include divx component

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    kloop:  DO k=1,nk

      IF(.not.terrain_flag)THEN
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=rho0(1,1,k)*s0(i,j,k)+dt*( advy(i,j,k)  &
                       +s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j) )
        enddo
        enddo
      ELSE
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)+dt*( advy(i,j,k)  &
                       +s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j) )*gz(i,j)
        enddo
        enddo
      ENDIF

        if(sbc.eq.2 .and. ibs.eq.1)then
          do i=1,ni
            if(rrv(i,1,k).ge.0.0)then
              j=1
              dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)
            endif
          enddo
        endif

        if(nbc.eq.2 .and. ibn.eq.1)then
          do i=1,ni
            if(rrv(i,nj+1,k).le.0.0)then
              j=nj
              dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)
            endif
          enddo
        endif

    ENDDO  kloop

        if(timestats.ge.1) time_pdef=time_pdef+mytime()

        call bcs(dum)
#ifdef MPI
        call comm_2sn_start(dum,south(1,1,1),newsouth(1,1,1),north(1,1,1),newnorth(1,1,1),reqs_s)
        call comm_2sn_end(dum,south(1,1,1),newsouth(1,1,1),north(1,1,1),newnorth(1,1,1),reqs_s)
#endif

        rdt=1.0/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k,foo1,foo2,foo3)
      DO k=1,nk

        do j=-1,nj+2
        do i=1,ni
          mass(i,j,k)=0.0
        enddo
        enddo

        do i=1,ni
        do j=-1,nj+2
          flag(i,j,k)=.false.
        enddo
        do j=0,nj+1
          if(dum(i,j,k).lt.0.0)then
            foo1=max(0.0,dum(i,j-1,k))
            foo2=max(0.0,dum(i,j+1,k))
            if(foo1+foo2.gt.1.0e-30)then
              foo3=max(dum(i,j,k),-(foo1+foo2))/(foo1+foo2)
              mass(i,j-1,k)=mass(i,j-1,k)+foo1*foo3
              mass(i,j  ,k)=mass(i,j  ,k)-(foo1+foo2)*foo3
              mass(i,j+1,k)=mass(i,j+1,k)+foo2*foo3
              if(dum(i,j-1,k).gt.1.0e-30) flag(i,j-1,k)=.true.
                                          flag(i,j  ,k)=.true.
              if(dum(i,j+1,k).gt.1.0e-30) flag(i,j+1,k)=.true.
            endif
          endif
        enddo
        !-----
        IF(.not.terrain_flag)THEN
        do j=1,nj
        if(flag(i,j,k))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advy(i,j,k)=(dum(i,j,k)-rho0(1,1,k)*s0(i,j,k))*rdt  &
                       -s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)
        endif
        enddo
        !-----
        ELSE
        do j=1,nj
        if(flag(i,j,k))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advy(i,j,k)=(dum(i,j,k)-rho0(i,j,k)*s0(i,j,k))*rdt*rgz(i,j)  &
                       -s(i,j,k)*(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)
        endif
        enddo
        ENDIF
        !-----
        enddo

      ENDDO

!----------------------------------------------------------------

      if(timestats.ge.1) time_pdef=time_pdef+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefz(mh,rho0,gz,rgz,rdsf,rrw,advz,dum,mass,s0,s,dt,flag)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rho0,s0,s
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz
      real, intent(in), dimension(kb:ke+1) :: rdsf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rrw
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: advz,dum,mass
      real, intent(in) :: dt
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag

      integer i,j,k
      real foo1,foo2,foo3,rdt

!----------------------------------------------------------------

        rdt=1.0/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k,foo1,foo2,foo3)
    DO j=1,nj

      IF(.not.terrain_flag)THEN
        do k=1,nk
        do i=1,ni
          dum(i,j,k)=rho0(1,1,k)*s0(i,j,k)+dt*( advz(i,j,k)  &
                       +s(i,j,k)*(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k) )
        enddo
        enddo
      ELSE
        do k=1,nk
        do i=1,ni
          dum(i,j,k)=rho0(i,j,k)*s0(i,j,k)+dt*( advz(i,j,k)  &
                       +s(i,j,k)*(rrw(i,j,k+1)-rrw(i,j,k))*rdsf(k) )*gz(i,j)
        enddo
        enddo
      ENDIF

        do i=1,ni
          dum(i,j, 0)=0.0
          dum(i,j,nk+1)=0.0
        enddo

        do k=0,nk+1
        do i=0,ni+1
          mass(i,j,k)=0.0
        enddo
        enddo

      do i=1,ni
        do k=0,nk+1
          flag(i,j,k)=.false.
        enddo
        do k=1,nk
          if(dum(i,j,k).lt.0.0)then
            foo1=max(0.0,dum(i,j,k-1))
            foo2=max(0.0,dum(i,j,k+1))
            if(foo1+foo2.gt.1.0e-30)then
              foo3=max(dum(i,j,k),-(foo1+foo2))/(foo1+foo2)
              mass(i,j,k-1)=mass(i,j,k-1)+foo1*foo3
              mass(i,j,k  )=mass(i,j,k  )-(foo1+foo2)*foo3
              mass(i,j,k+1)=mass(i,j,k+1)+foo2*foo3
              if(dum(i,j,k-1).gt.1.0e-30) flag(i,j,k-1)=.true.
                                          flag(i,j,k  )=.true.
              if(dum(i,j,k+1).gt.1.0e-30) flag(i,j,k+1)=.true.
            endif
          endif
        enddo
        !-----
        IF(.not.terrain_flag)THEN
        do k=1,nk
        if(flag(i,j,k))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advz(i,j,k)=(dum(i,j,k)-rho0(1,1,k)*s0(i,j,k))*rdt  &
                       -s(i,j,k)*(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)
        endif
        enddo
        !-----
        ELSE
        do k=1,nk
        if(flag(i,j,k))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advz(i,j,k)=(dum(i,j,k)-rho0(i,j,k)*s0(i,j,k))*rdt*rgz(i,j)  &
                       -s(i,j,k)*(rrw(i,j,k+1)-rrw(i,j,k))*rdsf(k)
        endif
        enddo
        ENDIF
        !-----
      enddo

    ENDDO

!----------------------------------------------------------------

      if(timestats.ge.1) time_pdef=time_pdef+mytime()

      return
      end


