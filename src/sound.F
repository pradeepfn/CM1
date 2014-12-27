

      subroutine sound(dt,xh,arh1,arh2,uh,ruh,xf,uf,yh,vh,rvh,yf,vf,      &
                        rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,mf,         &
                        pi0,rho0,rr0,rf0,rrf0,th0,th00,rth00,zs,          &
                        gz,rgz,gzu,rgzu,gzv,rgzv,                         &
                        dzdx,dzdy,gx,gxu,gy,gyv,                          &
                        radbcw,radbce,radbcs,radbcn,                      &
                        fwk,ppd,dum1,dum2,pk,mm,qk,fpk,dum3,              &
                        u0,ua,u3d,uten,                                   &
                        v0,va,v3d,vten,                                   &
                        wa,w3d,wten,                                      &
                        ppi,pp3d,ppten,tha,th3d,thten,thterm,             &
                        thv,ppterm,nrk,dttmp,rtime,th00s,thlr,            &
                        pkz,bkz,ftk,bk,sk,                                &
                        reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,               &
                        uw31,uw32,ue31,ue32,us31,us32,un31,un32,          &
                        vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,          &
                        ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,          &
                        sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,          &
                        rw31,rw32,re31,re32,rs31,rs32,rn31,rn32,          &
                        pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,rr0,rf0,rrf0,th0,th00,rth00
      real, intent(in), dimension(itb:ite,jtb:jte) :: zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: fwk,ppd,dum1,dum2,pk,mm,qk,fpk,dum3
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua,u3d,uten
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va,v3d,vten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa,w3d,wten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha,th3d,thten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: thterm,thv,ppterm
      integer, intent(in) :: nrk
      real, intent(in)  :: dttmp,rtime,th00s,thlr
      real, intent(inout), dimension(kb:ke) :: pkz,bkz
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ftk,bk,sk
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, intent(inout), dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
      real, intent(inout), dimension(cmp,jmp,kmp-1) :: ww31,ww32,we31,we32
      real, intent(inout), dimension(imp,cmp,kmp-1) :: ws31,ws32,wn31,wn32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(inout), dimension(cmp,jmp,kmp,2) :: rw31,rw32,re31,re32
      real, intent(inout), dimension(imp,cmp,kmp,2) :: rs31,rs32,rn31,rn32
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2

!-----

      integer :: i,j,k,n,nloop
      real :: tem,tem1,tem2,tem3,tem4,r1,r2,dts
      real :: aa,bb,cc,dd

      real :: beta,asq
      real :: temx,temy,temz,div,w1,w2

#ifdef MPI
      include 'mpif.h'
#endif

!---------------------------------------------------------------------
      ! beta = backward weight
      ! alph = forward weight

      beta=1.0-alph
      asq=alph*alph

      if(nrk.eq.1)then
!!!        nloop=1
!!!        dts=dt/3.
        nloop=nint(float(nsound)/3.0)
        dts=dt/(nloop*3.0)
        if( dts.gt.(dt/nsound) )then
          nloop=nloop+1
          dts=dt/(nloop*3.0)
        endif
      elseif(nrk.eq.2)then
        nloop=0.5*nsound
        dts=dt/nsound
      elseif(nrk.eq.3)then
        nloop=nsound
        dts=dt/nsound
      endif

!!!      print *,'  nloop,dts,dttmp = ',nloop,dts,nloop*dts

!-----------------------------------------------------------------
!  Arrays for implicit solver:

    IF(.not.terrain_flag)THEN

      ! without terrain:

      tem1 = dts*g*0.5/cp
      tem2 = dts*cp*rdz
      tem3 = dts*g
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem,r1)
      do k=1,nk
        tem = tem1*rr0(1,1,k)*rth00(1,1,k)
        pkz(k)=tem
        IF(k.gt.1)THEN
          bkz(k)=tem3/(c1(1,1,k)*th0(1,1,k-1)+c2(1,1,k)*th0(1,1,k))
          tem = tem2*mf(1,1,k)
          do j=1,nj
          do i=1,ni
            mm(i,j,k)= tem*(c1(1,1,k)*thv(i,j,k-1)+c2(1,1,k)*thv(i,j,k))
          enddo
          enddo
        ENDIF
        r1 = dts*rdz*mh(1,1,k)
        do j=1,nj
        do i=1,ni
          qk(i,j,k)=r1*ppterm(i,j,k)
          sk(i,j,k)=r1*thterm(i,j,k)
        enddo
        enddo
      enddo

      bkz(1) = 0.0
      bkz(nk+1) = 0.0

    ELSE

      ! with terrain:

      tem1 = dts*g*0.5/cp
      tem2 = dts*cp
      tem3 = dts*g
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          pk(i,j,k)=tem1*rr0(i,j,k)*rth00(i,j,k)
        enddo
        enddo
        IF(k.gt.1)THEN
          do j=1,nj
          do i=1,ni
            mm(i,j,k)=tem2*(c1(i,j,k)*thv(i,j,k-1)+c2(i,j,k)*thv(i,j,k))  &
                          *gz(i,j)*rds(k)
            bk(i,j,k)=tem3/(c1(i,j,k)*th0(i,j,k-1)+c2(i,j,k)*th0(i,j,k))
          enddo
          enddo
        ENDIF
        r1 = dts*rdsf(k)
        do j=1,nj
        do i=1,ni
          qk(i,j,k)=r1*gz(i,j)*ppterm(i,j,k)
          sk(i,j,k)=r1*gz(i,j)*thterm(i,j,k)
        enddo
        enddo
      enddo

    ENDIF


!---------------------------------------------------------------------
!  Prepare for acoustic steps

      if(nrk.ge.2)then

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1,ni+1
            u3d(i,j,k)=ua(i,j,k)
          enddo
          enddo
          IF(axisymm.eq.0)THEN
            ! Cartesian grid:
            do j=1,nj+1
            do i=1,ni
              v3d(i,j,k)=va(i,j,k)
            enddo
            enddo
          ENDIF
          IF(k.gt.1)THEN
            do j=1,nj
            do i=1,ni
              w3d(i,j,k)=wa(i,j,k)
            enddo
            enddo
          ENDIF
          do j=0,nj+1
          do i=0,ni+1
            pp3d(i,j,k)=ppi(i,j,k)
            ppd(i,j,k)=pp3d(i,j,k)
          enddo
          enddo
          do j=1,nj
          do i=1,ni
            th3d(i,j,k)=tha(i,j,k)
          enddo
          enddo
        enddo

      else

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            ppd(i,j,k)=pp3d(i,j,k)
          enddo
          enddo
        enddo

      endif

      if(timestats.ge.1) time_sound=time_sound+mytime()

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  Begin small steps:

      small_step_loop:  DO N=1,NLOOP

!-----

        if(irbc.eq.2)then
 
          if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,u3d)
 
          if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,v3d)
 
        endif

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  Open boundary conditions:
!  (cm1r17: proper accounting of terrain)

        ! west open bc tendency:
        IF(wbc.eq.2.and.ibw.eq.1)THEN
          IF(.not.terrain_flag)THEN
            ! no terrain:
!$omp parallel do default(shared)   &
!$omp private(j,k)
            do k=1,nk
            do j=1,nj
              u3d(1,j,k)=u3d(1,j,k)+dts*( -radbcw(j,k)       &
                        *(u3d(2,j,k)-u3d(1,j,k))*rdx*uh(1)   &
                           +uten(1,j,k) )
            enddo
            enddo
          ELSE
            ! with terrain:
            ! dum1 stores u at w points:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
            do j=1,nj
              ! lowest model level:
              do i=1,2
                dum1(i,j,1) = cgs1*u3d(i,j,1)+cgs2*u3d(i,j,2)+cgs3*u3d(i,j,3)
              enddo
              ! upper-most model level:
              do i=1,2
                dum1(i,j,nk+1) = cgt1*u3d(i,j,nk)+cgt2*u3d(i,j,nk-1)+cgt3*u3d(i,j,nk-2)
              enddo
              ! interior:
              do k=2,nk
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              do i=1,2
                dum1(i,j,k) = r1*u3d(i,j,k-1)+r2*u3d(i,j,k)
              enddo
              enddo
            enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
              u3d(1,j,k)=u3d(1,j,k)+dts*( -radbcw(j,k)*(                            &
                     gz(1,j)*(u3d(2,j,k)*rgzu(2,j)-u3d(1,j,k)*rgzu(1,j))*rdx*uh(1)  &
                    +0.5*( gx(1,j,k+1)*(dum1(2,j,k+1)+dum1(1,j,k+1))                &
                          -gx(1,j,k  )*(dum1(2,j,k  )+dum1(1,j,k  )) )*rdsf(k)      &
                                                       )+uten(1,j,k) )
            enddo
            enddo
          ENDIF  ! end check for terrain
        ENDIF    ! end open bc tendency

!-----

        ! east open bc tendency:
        IF(ebc.eq.2.and.ibe.eq.1)THEN
          IF(.not.terrain_flag)THEN
            ! no terrain:
!$omp parallel do default(shared)   &
!$omp private(j,k)
            do k=1,nk
            do j=1,nj
              u3d(ni+1,j,k)=u3d(ni+1,j,k)+dts*( -radbce(j,k)           &
                           *(u3d(ni+1,j,k)-u3d(ni  ,j,k))*rdx*uh(ni)   &
                           +uten(ni+1,j,k) )
            enddo
            enddo
          ELSE
            ! with terrain:
            ! dum1 stores u at w points:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
            do j=1,nj
              ! lowest model level:
              do i=ni,ni+1
                dum1(i,j,1) = cgs1*u3d(i,j,1)+cgs2*u3d(i,j,2)+cgs3*u3d(i,j,3)
              enddo
              ! upper-most model level:
              do i=ni,ni+1
                dum1(i,j,nk+1) = cgt1*u3d(i,j,nk)+cgt2*u3d(i,j,nk-1)+cgt3*u3d(i,j,nk-2)
              enddo
              ! interior:
              do k=2,nk
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              do i=ni,ni+1
                dum1(i,j,k) = r1*u3d(i,j,k-1)+r2*u3d(i,j,k)
              enddo
              enddo
            enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
              u3d(ni+1,j,k)=u3d(ni+1,j,k)+dts*( -radbce(j,k)*(                                &
                     gz(ni,j)*(u3d(ni+1,j,k)*rgzu(ni+1,j)-u3d(ni,j,k)*rgzu(ni,j))*rdx*uh(ni)  &
                    +0.5*( gx(ni,j,k+1)*(dum1(ni+1,j,k+1)+dum1(ni,j,k+1))                     &
                          -gx(ni,j,k  )*(dum1(ni+1,j,k  )+dum1(ni,j,k  )) )*rdsf(k)           &
                                                             )+uten(ni+1,j,k) )
            enddo
            enddo
          ENDIF  ! end check for terrain
        ENDIF    ! end open bc tendency

        if(timestats.ge.1) time_sound=time_sound+mytime()

!-----

        IF(roflux.eq.1)THEN
          call restrict_openbc_we(rvh,rmh,rho0,u3d)
        ENDIF

!-----

        ! south open bc tendency:
        IF(sbc.eq.2.and.ibs.eq.1)THEN
          IF(.not.terrain_flag)THEN
            ! no terrain:
!$omp parallel do default(shared)   &
!$omp private(i,k)
            do k=1,nk
            do i=1,ni
              v3d(i,1,k)=v3d(i,1,k)+dts*( -radbcs(i,k)      &
                        *(v3d(i,2,k)-v3d(i,1,k))*rdy*vh(1)  &
                        +vten(i,1,k) )
            enddo
            enddo
          ELSE
            ! with terrain:
            ! dum1 stores v at w points:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
            do j=1,2
              ! lowest model level:
              do i=1,ni
                dum1(i,j,1) = cgs1*v3d(i,j,1)+cgs2*v3d(i,j,2)+cgs3*v3d(i,j,3)
              enddo
              ! upper-most model level:
              do i=1,ni
                dum1(i,j,nk+1) = cgt1*v3d(i,j,nk)+cgt2*v3d(i,j,nk-1)+cgt3*v3d(i,j,nk-2)
              enddo
              ! interior:
              do k=2,nk
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              do i=1,ni
                dum1(i,j,k) = r1*v3d(i,j,k-1)+r2*v3d(i,j,k)
              enddo
              enddo
            enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do i=1,ni
              v3d(i,1,k)=v3d(i,1,k)+dts*( -radbcs(i,k)*(                            &
                     gz(i,1)*(v3d(i,2,k)*rgzv(i,2)-v3d(i,1,k)*rgzv(i,1))*rdy*vh(1)  &
                    +0.5*( gy(i,1,k+1)*(dum1(i,2,k+1)+dum1(i,1,k+1))                &
                          -gy(i,1,k  )*(dum1(i,2,k  )+dum1(i,1,k  )) )*rdsf(k)      &
                                                       )+vten(i,1,k) )
            enddo
            enddo
          ENDIF  ! end check for terrain
        ENDIF    ! end open bc tendency

!-----
 
        ! north open bc tendency:
        IF(nbc.eq.2.and.ibn.eq.1)THEN
          IF(.not.terrain_flag)THEN
            ! no terrain:
!$omp parallel do default(shared)   &
!$omp private(i,k)
            do k=1,nk
            do i=1,ni
              v3d(i,nj+1,k)=v3d(i,nj+1,k)+dts*( -radbcn(i,k)        &
                           *(v3d(i,nj+1,k)-v3d(i,nj,k))*rdy*vh(nj)  &
                           +vten(i,nj+1,k) )
            enddo
            enddo
          ELSE
            ! with terrain:
            ! dum1 stores v at w points:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
            do j=nj,nj+1
              ! lowest model level:
              do i=1,ni
                dum1(i,j,1) = cgs1*v3d(i,j,1)+cgs2*v3d(i,j,2)+cgs3*v3d(i,j,3)
              enddo
              ! upper-most model level:
              do i=1,ni
                dum1(i,j,nk+1) = cgt1*v3d(i,j,nk)+cgt2*v3d(i,j,nk-1)+cgt3*v3d(i,j,nk-2)
              enddo
              ! interior:
              do k=2,nk
              r2 = (sigmaf(k)-sigma(k-1))*rds(k)
              r1 = 1.0-r2
              do i=1,ni
                dum1(i,j,k) = r1*v3d(i,j,k-1)+r2*v3d(i,j,k)
              enddo
              enddo
            enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do i=1,ni
              v3d(i,nj+1,k)=v3d(i,nj+1,k)+dts*( -radbcn(i,k)*(                                &
                     gz(i,nj)*(v3d(i,nj+1,k)*rgzv(i,nj+1)-v3d(i,nj,k)*rgzv(i,nj))*rdy*vh(nj)  &
                    +0.5*( gy(i,nj,k+1)*(dum1(i,nj+1,k+1)+dum1(i,nj,k+1))                     &
                          -gy(i,nj,k  )*(dum1(i,nj+1,k  )+dum1(i,nj,k  )) )*rdsf(k)           &
                                                             )+vten(i,nj+1,k) )
            enddo
            enddo
          ENDIF  ! end check for terrain
        ENDIF    ! end open bc tendency

        if(timestats.ge.1) time_sound=time_sound+mytime()

!-----

        IF(roflux.eq.1)THEN
          call restrict_openbc_sn(ruh,rmh,rho0,v3d)
        ENDIF

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  integrate u,v forward in time:

#ifdef MPI
        IF(n.gt.1)THEN
          call comm_1s_end(ppd,pw1,pw2,pe1,pe2,   &
                               ps1,ps2,pn1,pn2,reqs_p)
        ENDIF
#endif
!-----

    IF(.not.terrain_flag)THEN

      IF(axisymm.eq.0)THEN
        ! Cartesian grid without terrain:

        tem1 = rdx*cp*0.5
        tem2 = rdy*cp*0.5
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)         &
                   -tem1*(ppd(i,j,k)-ppd(i-1,j,k))*uf(i)    &
                        *(thv(i,j,k)+thv(i-1,j,k)) )
          enddo
          enddo
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            v3d(i,j,k)=v3d(i,j,k)+dts*( vten(i,j,k)         &
                   -tem2*(ppd(i,j,k)-ppd(i,j-1,k))*vf(j)    &
                        *(thv(i,j,k)+thv(i,j-1,k)) )
          enddo
          enddo
        enddo

      ELSE
        ! axisymmetric grid:

        tem1 = rdx*cp*0.5
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)         &
                   -tem1*(ppd(i,j,k)-ppd(i-1,j,k))*uf(i)    &
                        *(thv(i,j,k)+thv(i-1,j,k)) )
          enddo
          enddo
        enddo

      ENDIF

    ELSE

        ! Cartesian grid with terrain:

        ! dum1 stores ppd at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
        do j=0,nj+1

          ! lowest model level:
          do i=0,ni+1
            dum1(i,j,1) = cgs1*ppd(i,j,1)+cgs2*ppd(i,j,2)+cgs3*ppd(i,j,3)
          enddo

          ! upper-most model level:
          do i=0,ni+1
            dum1(i,j,nk+1) = cgt1*ppd(i,j,nk)+cgt2*ppd(i,j,nk-1)+cgt3*ppd(i,j,nk-2)
          enddo

          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=0,ni+1
            dum1(i,j,k) = r1*ppd(i,j,k-1)+r2*ppd(i,j,k)
          enddo
          enddo

        enddo

        tem1 = rdx*cp*0.5
        tem2 = rdy*cp*0.5
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
          ! x-dir
          do j=1,nj
          do i=1+ibw,ni+1-ibe
            u3d(i,j,k)=u3d(i,j,k)+dts*( uten(i,j,k)               &
                   -cp*0.5*(thv(i,j,k)+thv(i-1,j,k))*(            &
                  (ppd(i,j,k)*rgz(i,j)-ppd(i-1,j,k)*rgz(i-1,j))   &
                   *gzu(i,j)*rdx*uf(i)                            &
              +0.5*( gxu(i,j,k+1)*(dum1(i,j,k+1)+dum1(i-1,j,k+1)) &
                    -gxu(i,j,k  )*(dum1(i,j,k  )+dum1(i-1,j,k  )) &
                   )*rdsf(k) ) )
          enddo
          enddo
          do j=1+ibs,nj+1-ibn
          do i=1,ni
            v3d(i,j,k)=v3d(i,j,k)+dts*( vten(i,j,k)               &
                   -cp*0.5*(thv(i,j,k)+thv(i,j-1,k))*(            &
                  (ppd(i,j,k)*rgz(i,j)-ppd(i,j-1,k)*rgz(i,j-1))   &
                   *gzv(i,j)*rdy*vf(j)                            &
              +0.5*( gyv(i,j,k+1)*(dum1(i,j,k+1)+dum1(i,j-1,k+1)) &
                    -gyv(i,j,k  )*(dum1(i,j,k  )+dum1(i,j-1,k  )) &
                   )*rdsf(k) ) )
          enddo
          enddo
        enddo

    ENDIF

!----------------------------------------------
!  convergence forcing:

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. nx.gt.1 )THEN
            call convinitu(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibw,ibe,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xf,yh,zh,u0,u3d)
          ENDIF
        ENDIF

!----------------------------------------------

        if(timestats.ge.1) time_sound=time_sound+mytime()

        if(n.eq.nloop)then
          call bcu(u3d)
        endif

#ifdef MPI
        if(n.eq.nloop)then
          call comm_3u_start(u3d,uw31,uw32,ue31,ue32,   &
                                 us31,us32,un31,un32,reqs_u)
        endif
#endif

!----------------------------------------------

      IF(axisymm.eq.0)THEN
        ! Cartesian grid:

!----------------------------------------------
!  convergence forcing:

        IF( convinit.eq.1 )THEN
          IF( rtime.le.convtime .and. ny.gt.1 )THEN
            call convinitv(myid,ib,ie,jb,je,kb,ke,ni,nj,nk,ibs,ibn,   &
                           zdeep,lamx,lamy,xcent,ycent,aconv,    &
                           xh,yf,zh,v0,v3d)
          ENDIF
        ENDIF

!----------------------------------------------

        if(timestats.ge.1) time_sound=time_sound+mytime()

        if(n.eq.nloop)then
          call bcv(v3d)
        endif

#ifdef MPI
        if(n.eq.nloop)then
          call comm_3v_start(v3d,vw31,vw32,ve31,ve32,   &
                                 vs31,vs32,vn31,vn32,reqs_v)
        endif
#endif

      ENDIF

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  get variables for terrain:

    IF(terrain_flag)THEN
      call bcwsfc(gz,dzdx,dzdy,u3d,v3d,w3d)
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk
        do j=1,nj
        do i=1,ni+1
          dum1(i,j,k)=u3d(i,j,k)*rgzu(i,j)
        enddo
        enddo
        do j=1,nj+1
        do i=1,ni
          dum2(i,j,k)=v3d(i,j,k)*rgzv(i,j)
        enddo
        enddo
      ENDDO
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      DO k=1,nk
        IF(k.eq.1)THEN
          do j=1,nj
          do i=1,ni
            dum3(i,j,1)=0.0
            dum3(i,j,nk+1)=0.0
          enddo
          enddo
        ELSE
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do j=1,nj
          do i=1,ni
            dum3(i,j,k)=0.5*( ( r2*(dum1(i,j,k  )+dum1(i+1,j,k  ))               &
                               +r1*(dum1(i,j,k-1)+dum1(i+1,j,k-1)) )*dzdx(i,j)   &
                             +( r2*(dum2(i,j,k  )+dum2(i,j+1,k  ))               &
                               +r1*(dum2(i,j,k-1)+dum2(i,j+1,k-1)) )*dzdy(i,j)   &
                           )*(sigmaf(k)-zt)*gz(i,j)*rzt
          enddo
          enddo
        ENDIF
      ENDDO
    ENDIF

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  update arrays for vertically implicit solver

      tem1 = dts*rdz
      tem2 = -0.5*thlr*dts

      temx = -dts*0.5*rdx
      temy = -dts*0.5*rdy
      temz = -dts*0.5*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k,div,w1,w2)
    DO k=1,nk
    IF(axisymm.eq.0)THEN
      IF(.not.terrain_flag)THEN
        ! Cartesian grid without terrain:
        do j=1,nj
        do i=1,ni
          div = (u3d(i+1,j,k)-u3d(i,j,k))*rdx*uh(i)    &
               +(v3d(i,j+1,k)-v3d(i,j,k))*rdy*vh(j)
          fpk(i,j,k)=dts*( ppten(i,j,k)-ppterm(i,j,k)*div )                 &
              +beta*( w3d(i,j,k  )*(pkz(k)*rf0(1,1,k  )+qk(i,j,k))          &
                     +w3d(i,j,k+1)*(pkz(k)*rf0(1,1,k+1)-qk(i,j,k)) )
          ftk(i,j,k)=dts*( thten(i,j,k)-thterm(i,j,k)*div )                 &
              +beta*( w3d(i,j,k  )*(tem2*rr0(1,1,k)*rf0(1,1,k  )+sk(i,j,k)) &
                     +w3d(i,j,k+1)*(tem2*rr0(1,1,k)*rf0(1,1,k+1)-sk(i,j,k)) )
        enddo
        enddo
      ELSE
        ! Cartesian grid with terrain:
        do j=1,nj
        do i=1,ni
          if(k.eq. 1)then
            w1 = 0.0
          else
            w1 = w3d(i,j,k)
          endif
          if(k.eq.nk)then
            w2 = 0.0
          else
            w2 = w3d(i,j,k+1)
          endif
          div = gz(i,j)*(                             &
                (dum1(i+1,j,k)-dum1(i,j,k))*rdx*uh(i) &
               +(dum2(i,j+1,k)-dum2(i,j,k))*rdy*vh(j) &
               +(dum3(i,j,k+1)-dum3(i,j,k))*rdsf(k) )
          fpk(i,j,k)=dts*( ppten(i,j,k)-ppterm(i,j,k)*div )       &
              +beta*( w1*(pk(i,j,k)*rf0(i,j,k  )+qk(i,j,k))       &
                     +w2*(pk(i,j,k)*rf0(i,j,k+1)-qk(i,j,k)) )
          ftk(i,j,k)=dts*( thten(i,j,k)-thterm(i,j,k)*div )       &
              +beta*( w1*(tem2*rr0(i,j,k)*rf0(i,j,k  )+sk(i,j,k)) &
                     +w2*(tem2*rr0(i,j,k)*rf0(i,j,k+1)-sk(i,j,k)) )
        enddo
        enddo
      ENDIF
    ELSE
        ! axisymmetric grid:
        do j=1,nj
        do i=1,ni
          div = (arh2(i)*u3d(i+1,j,k)-arh1(i)*u3d(i,j,k))*rdx*uh(i)
          fpk(i,j,k)=dts*( ppten(i,j,k)-ppterm(i,j,k)*div )                 &
              +beta*( w3d(i,j,k  )*(pkz(k)*rf0(1,1,k  )+qk(i,j,k))          &
                     +w3d(i,j,k+1)*(pkz(k)*rf0(1,1,k+1)-qk(i,j,k)) )
          ftk(i,j,k)=dts*( thten(i,j,k)-thterm(i,j,k)*div )                 &
              +beta*( w3d(i,j,k  )*(tem2*rr0(1,1,k)*rf0(1,1,k  )+sk(i,j,k)) &
                     +w3d(i,j,k+1)*(tem2*rr0(1,1,k)*rf0(1,1,k+1)-sk(i,j,k)) )
        enddo
        enddo
    ENDIF

    IF(k.gt.1)THEN
      IF(.not.terrain_flag)THEN
        do j=1,nj
        do i=1,ni
          fwk(i,j,k)=dts*wten(i,j,k)                                        &
             +beta*( bkz(k)*(c2(1,1,k)*th3d(i,j,k)+c1(1,1,k)*th3d(i,j,k-1)) &
                    -mm(i,j,k)*(pp3d(i,j,k)-pp3d(i,j,k-1)) )
        enddo
        enddo
      ELSE
        do j=1,nj
        do i=1,ni
          fwk(i,j,k)=dts*wten(i,j,k)                                           &
             +beta*( bk(i,j,k)*(c2(i,j,k)*th3d(i,j,k)+c1(i,j,k)*th3d(i,j,k-1)) &
                    -mm(i,j,k)*(pp3d(i,j,k)-pp3d(i,j,k-1)) )
        enddo
        enddo
      ENDIF
    ENDIF

    ENDDO

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  the vertcally implicit solver .... get new w

        tem1 = -0.5*thlr*dts

!$omp parallel do default(shared)   &
!$omp private(i,j,k,aa,bb,cc,dd,r1)
      DO j=1,nj

      IF(.not.terrain_flag)THEN
        ! without terrain:
        do i=1,ni
          dum1(i,j,1)=0.0
          dum2(i,j,1)=w3d(i,j,1)
        enddo
        do k=2,nk
        do i=1,ni
          aa=asq*( mm(i,j,k)*(qk(i,j,k)-pkz(k)*rf0(1,1,k+1))   &
                  +bkz(k)*(tem1*rf0(1,1,k+1)*c2(1,1,k)*rr0(1,1,k)-sk(i,j,k)) )
          bb=1.0+asq*( mm(i,j,k)*( qk(i,j,k)+qk(i,j,k-1)                    &
                                  +pkz(k)*rf0(1,1,k)-pkz(k-1)*rf0(1,1,k) ) &
                      -bkz(k)*(tem1*rf0(1,1,k)*(c2(1,1,k)*rr0(1,1,k)+c1(1,1,k)*rr0(1,1,k-1))   &
                                  +sk(i,j,k)-sk(i,j,k-1)) )
          cc=asq*( mm(i,j,k)*(pkz(k-1)*rf0(1,1,k-1)+qk(i,j,k-1))   &
                  +bkz(k)*(tem1*rf0(1,1,k-1)*c1(1,1,k)*rr0(1,1,k-1)+sk(i,j,k-1)) )
          dd=w3d(i,j,k)+fwk(i,j,k)                                                    &
              -alph*( mm(i,j,k)*(pp3d(i,j,k)+fpk(i,j,k)-pp3d(i,j,k-1)-fpk(i,j,k-1))   &
                     -bkz(k)*( c2(1,1,k)*(th3d(i,j,k)+ftk(i,j,k))+c1(1,1,k)*(th3d(i,j,k-1)+ftk(i,j,k-1)) ) )
          r1 = 1.0/(bb-cc*dum1(i,j,k-1))
          dum1(i,j,k)=aa*r1
          dum2(i,j,k)=(dd+cc*dum2(i,j,k-1))*r1
        enddo
        enddo
      ELSE
        ! with terrain:
        do i=1,ni
          w3d(i,j,1)=0.0
          dum1(i,j,1)=0.0
          dum2(i,j,1)=w3d(i,j,1)
        enddo
        do k=2,nk
        do i=1,ni
        if(k.eq.nk)then
          aa=0.0
        else
          aa=asq*( mm(i,j,k)*(qk(i,j,k)-pk(i,j,k)*rf0(i,j,k+1))   &
                  +bk(i,j,k)*(tem1*rf0(i,j,k+1)*c2(i,j,k)*rr0(i,j,k)-sk(i,j,k)) )
        endif
          bb=1.0+asq*( mm(i,j,k)*( qk(i,j,k)+qk(i,j,k-1)                    &
                                  +pk(i,j,k)*rf0(i,j,k)-pk(i,j,k-1)*rf0(i,j,k) ) &
                      -bk(i,j,k)*(tem1*rf0(i,j,k)*(c2(i,j,k)*rr0(i,j,k)+c1(i,j,k)*rr0(i,j,k-1))   &
                                  +sk(i,j,k)-sk(i,j,k-1)) )
        if(k.eq.2)then
          cc=0.0
        else
          cc=asq*( mm(i,j,k)*(pk(i,j,k-1)*rf0(i,j,k-1)+qk(i,j,k-1))   &
                  +bk(i,j,k)*(tem1*rf0(i,j,k-1)*c1(i,j,k)*rr0(i,j,k-1)+sk(i,j,k-1)) )
        endif
          dd=w3d(i,j,k)+fwk(i,j,k)                                                    &
              -alph*( mm(i,j,k)*(pp3d(i,j,k)+fpk(i,j,k)-pp3d(i,j,k-1)-fpk(i,j,k-1))   &
                     -bk(i,j,k)*( c2(i,j,k)*(th3d(i,j,k)+ftk(i,j,k))+c1(i,j,k)*(th3d(i,j,k-1)+ftk(i,j,k-1)) ) )
          r1 = 1.0/(bb-cc*dum1(i,j,k-1))
          dum1(i,j,k)=aa*r1
          dum2(i,j,k)=(dd+cc*dum2(i,j,k-1))*r1
        enddo
        enddo
      ENDIF

        do i=1,ni
          w3d(i,j,nk+1)=0.0
        enddo

        do k=nk,2,-1
        do i=1,ni
          w3d(i,j,k)=dum1(i,j,k)*w3d(i,j,k+1)+dum2(i,j,k)
        enddo
        enddo

      ENDDO
      if(timestats.ge.1) time_sound=time_sound+mytime()

        if(n.eq.nloop)then
          call bcw(w3d,1)
          if(terrain_flag)  &
          call bcwsfc(gz,dzdx,dzdy,u3d,v3d,w3d)
        endif

#ifdef MPI
        if(n.eq.nloop)then
          call comm_3w_start(w3d,ww31,ww32,we31,we32,   &
                                 ws31,ws32,wn31,wn32,reqs_w)
        endif
#endif

!-----------------------------------------------------------------------
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------------------------------------------------
!  get new pp,th

      tem2 = -0.5*thlr*dts
 
!!$omp parallel do default(shared)   &
!!$omp private(i,j,k,w1,w2)
      DO k=1,nk

      IF(.not.terrain_flag)THEN
        do j=1,nj
        do i=1,ni
          ! without terrain:
          !-----
          ppd(i,j,k)=pp3d(i,j,k)
          pp3d(i,j,k)=pp3d(i,j,k)+fpk(i,j,k)                       &
              +alph*( w3d(i,j,k  )*(pkz(k)*rf0(1,1,k  )+qk(i,j,k)) &
                     +w3d(i,j,k+1)*(pkz(k)*rf0(1,1,k+1)-qk(i,j,k)) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          ppd(i,j,k)=pp3d(i,j,k)+kdiv*( pp3d(i,j,k)-ppd(i,j,k) )
          !-----
          th3d(i,j,k)=th3d(i,j,k)+ftk(i,j,k)                      &
              +alph*( w3d(i,j,k  )*(tem2*rr0(1,1,k)*rf0(1,1,k  )+sk(i,j,k)) &
                     +w3d(i,j,k+1)*(tem2*rr0(1,1,k)*rf0(1,1,k+1)-sk(i,j,k)) )
          if(abs(th3d(i,j,k)).lt.smeps) th3d(i,j,k)=0.0
          !-----
        enddo
        enddo
      ELSE
        do j=1,nj
        do i=1,ni
          ! with terrain:
          if(k.eq. 1)then
            w1 = 0.0
          else
            w1 = w3d(i,j,k)
          endif
          if(k.eq.nk)then
            w2 = 0.0
          else
            w2 = w3d(i,j,k+1)
          endif
          !-----
          ppd(i,j,k)=pp3d(i,j,k)
          pp3d(i,j,k)=pp3d(i,j,k)+fpk(i,j,k)                          &
              +alph*( w1*(pk(i,j,k)*rf0(i,j,k  )+qk(i,j,k)) &
                     +w2*(pk(i,j,k)*rf0(i,j,k+1)-qk(i,j,k)) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          ppd(i,j,k)=pp3d(i,j,k)+kdiv*( pp3d(i,j,k)-ppd(i,j,k) )
          !-----
          th3d(i,j,k)=th3d(i,j,k)+ftk(i,j,k)                      &
              +alph*( w1*(tem2*rr0(i,j,k)*rf0(i,j,k  )+sk(i,j,k)) &
                     +w2*(tem2*rr0(i,j,k)*rf0(i,j,k+1)-sk(i,j,k)) )
          if(abs(th3d(i,j,k)).lt.smeps) th3d(i,j,k)=0.0
          !-----
        enddo
        enddo
      ENDIF

      ENDDO
      if(timestats.ge.1) time_sound=time_sound+mytime()

        IF(n.lt.nloop)THEN
          call bcs(ppd)
#ifdef MPI
          call comm_1s_start(ppd,pw1,pw2,pe1,pe2,   &
                                 ps1,ps2,pn1,pn2,reqs_p)
#endif
        ENDIF

        if(n.eq.nloop)then
          if(nrk.lt.3.or.imoist.eq.0)then
            call bcs(th3d)
            call bcs(pp3d)
#ifdef MPI
            call comm_3r_start(th3d,pp3d,rw31,rw32,re31,re32,   &
                                         rs31,rs32,rn31,rn32,reqs_p)
#endif
          endif
        endif

!--------------------------------------------------------------------

      ENDDO  small_step_loop

!  end of small steps
!--------------------------------------------------------------------

      return
      end


