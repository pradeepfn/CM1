!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefx(xh,rho0,advx,dum,mass,s1,dt,flag,west,newwest,east,neweast)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie) :: xh
      real, dimension(ib:ie,jb:je,kb:ke) :: rho0
      real, dimension(ib:ie,jb:je,kb:ke) :: advx,dum,mass,s1
      real :: dt
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(cmp,jmp,kmp) :: west,newwest,east,neweast

      integer i,j,k
      real foo1,foo2,foo3,rdt

#ifdef MPI
      include 'mpif.h'
      integer, dimension(4) :: reqs_s
#endif

!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=rho0(i,j,k)*s1(i,j,k)+dt*advx(i,j,k)
        enddo
        enddo
        enddo
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

        if(ibw.eq.1.and.wbc.eq.2)then
          do j=1,nj
            dum(-1,j,k)=0.0
            dum( 0,j,k)=0.0
            dum( 1,j,k)=0.0
          enddo
        endif

        if(ibe.eq.1.and.ebc.eq.2)then
          do j=1,nj
            dum(ni  ,j,k)=0.0
            dum(ni+1,j,k)=0.0
            dum(ni+2,j,k)=0.0
          enddo
        endif

      IF(axisymm.eq.0)THEN

        do j=1,nj
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
        do i=1,ni
        if(flag(i,j,k))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advx(i,j,k)=(dum(i,j,k)-rho0(i,j,k)*s1(i,j,k))*rdt
        endif
        enddo
        enddo

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
          advx(i,j,k)=(dum(i,j,k)-rho0(i,j,k)*s1(i,j,k))*rdt
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


      subroutine pdefy(rho0,advy,dum,mass,s1,dt,flag,south,newsouth,north,newnorth)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: rho0
      real, dimension(ib:ie,jb:je,kb:ke) :: advy,dum,mass,s1
      real :: dt
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      real, intent(inout), dimension(imp,cmp,kmp) :: south,newsouth,north,newnorth

      integer i,j,k
      real foo1,foo2,foo3,rdt

#ifdef MPI
      include 'mpif.h'
      integer, dimension(4) :: reqs_s
#endif

!----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum(i,j,k)=rho0(i,j,k)*s1(i,j,k)+dt*advy(i,j,k)
        enddo
        enddo
        enddo
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

        if(ibs.eq.1.and.sbc.eq.2)then
          do i=1,ni
            dum(i,-1,k)=0.0
            dum(i, 0,k)=0.0
            dum(i, 1,k)=0.0
          enddo
        endif

        if(ibn.eq.1.and.nbc.eq.2)then
          do i=1,ni
            dum(i,nj  ,k)=0.0
            dum(i,nj+1,k)=0.0
            dum(i,nj+2,k)=0.0
          enddo
        elseif(ibn.eq.1.and.nbc.eq.3)then
          do i=1,ni
            dum(i,nj+1,k)=dum(i,nj,k)
            dum(i,nj+2,k)=dum(i,nj,k)
          enddo
        endif

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
        do j=1,nj
        if(flag(i,j,k))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advy(i,j,k)=(dum(i,j,k)-rho0(i,j,k)*s1(i,j,k))*rdt
        endif
        enddo
        enddo

      ENDDO

!----------------------------------------------------------------

      if(timestats.ge.1) time_pdef=time_pdef+mytime()

      return
      end


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine pdefz(rho0,advz,dum,mass,s1,dt,flag)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, dimension(ib:ie,jb:je,kb:ke) :: rho0
      real, dimension(ib:ie,jb:je,kb:ke) :: advz,dum,mass,s1
      real :: dt
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag

      integer i,j,k
      real foo1,foo2,foo3,rdt

!----------------------------------------------------------------

        rdt=1.0/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k,foo1,foo2,foo3)
    DO j=1,nj

        do k=1,nk
        do i=1,ni
          dum(i,j,k)=rho0(i,j,k)*s1(i,j,k)+dt*advz(i,j,k)
        enddo
        enddo

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
        do k=1,nk
        if(flag(i,j,k))then
          dum(i,j,k)=dum(i,j,k)+mass(i,j,k)
          advz(i,j,k)=(dum(i,j,k)-rho0(i,j,k)*s1(i,j,k))*rdt
        endif
        enddo
      enddo

    ENDDO

!----------------------------------------------------------------

      if(timestats.ge.1) time_pdef=time_pdef+mytime()

      return
      end


