

      subroutine sound(dt,th00s,thlr,xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf, &
                        rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,zf,mf,rmf,  &
                        pi0,th00,thv0,rho0,rr0,rf0,th0,zs,                &
                        gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv, &
                        radbcw,radbce,radbcs,radbcn,                    &
                        fwk,ppd,dum1,dum2,pk,mm,qk,fpk,                 &
                        u0,ua,u3d,uten,v0,va,v3d,vten,wa,w3d,wten,      &
                        ppi,pp3d,ppten,tha,th3d,thten,thterm,           &
                        thv,ppterm,ftk,bk,sk,nrk,rtime,                 &
                        reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,             &
                        uw31,uw32,ue31,ue32,us31,us32,un31,un32,        &
                        vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,        &
                        ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,        &
                        sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,        &
                        rw31,rw32,re31,re32,rs31,rs32,rn31,rn32,        &
                        pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt,th00s,thlr
      real, intent(in), dimension(ib:ie) :: xh,rxh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf,rmf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,th00,thv0,rho0,rr0,rf0,th0
      real, intent(in), dimension(itb:ite,jtb:jte) :: zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: fwk,ppd,dum1,dum2,pk,mm,qk,fpk
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua,u3d,uten
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va,v3d,vten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa,w3d,wten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha,th3d,thten,thterm
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: thv,ppterm,ftk,bk,sk
      integer, intent(in) :: nrk
      real, intent(in)  :: rtime
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
      real :: tem,tem1,tem2,tem3,r1,r2,dts
      real :: aa,bb,cc,dd

!!!      real alph,beta,asq
!!!      parameter(alph=0.60)
!!!      parameter(beta=1.0-alph)
!!!      parameter(asq=alph*alph)

      real :: beta,asq

#ifdef MPI
      include 'mpif.h'
#endif

!---------------------------------------------------------------------

      beta=1.0-alph
      asq=alph*alph

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
          do j=1,nj
          do i=1,ni
            w3d(i,j,k)=wa(i,j,k)
          enddo
          enddo
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

!---------------------------------------------------------------------
 
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

    IF(terrain_flag)THEN

      ! with terrain:

      tem1 = dts*g*0.5/cp
      tem2 = dts*cp
      tem3 = dts*g
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          pk(i,j,k)=tem1*rr0(i,j,k)/th00(i,j,k)
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
      enddo

    ELSE

      ! without terrain:

      tem1 = dts*g*0.5/cp
      tem2 = dts*cp*rdz
      tem3 = dts*g
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem)
      do k=1,nk
        tem = tem1*rr0(1,1,k)/th00(1,1,k)
        do j=1,nj
        do i=1,ni
          pk(i,j,k)=tem
        enddo
        enddo
        IF(k.gt.1)THEN
          tem = tem2*mf(1,1,k)
          do j=1,nj
          do i=1,ni
            mm(i,j,k)= tem*(c1(1,1,k)*thv(i,j,k-1)+c2(1,1,k)*thv(i,j,k))
            bk(i,j,k)=tem3/(c1(1,1,k)*th0(i,j,k-1)+c2(1,1,k)*th0(i,j,k))
          enddo
          enddo
        ENDIF
      enddo

    ENDIF

      if(timestats.ge.1) time_sound=time_sound+mytime()

!---------------------------------------------------------------------

      DO N=1,NLOOP

!-----

        if(irbc.eq.2)then
 
          if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,u3d)
 
          if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,v3d)
 
        endif

!-----

        if(wbc.eq.2.and.ibw.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            u3d(1,j,k)=u3d(1,j,k)+dts*( -radbcw(j,k)       &
                      *(u3d(2,j,k)-u3d(1,j,k))*rdx*uh(1)   &
                         +uten(1,j,k) )
          enddo
          enddo
        endif

        if(ebc.eq.2.and.ibe.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(j,k)
          do k=1,nk
          do j=1,nj
            u3d(ni+1,j,k)=u3d(ni+1,j,k)+dts*( -radbce(j,k)           &
                         *(u3d(ni+1,j,k)-u3d(ni  ,j,k))*rdx*uh(ni)   &
                         +uten(ni+1,j,k) )
          enddo
          enddo
        endif

        if(timestats.ge.1) time_sound=time_sound+mytime()

        IF(roflux.eq.1)THEN
          call restrict_openbc_we(rvh,rmh,rho0,u3d)
        ENDIF

!-----

        if(sbc.eq.2.and.ibs.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            v3d(i,1,k)=v3d(i,1,k)+dts*( -radbcs(i,k)       &
                      *(v3d(i,2,k)-v3d(i,1,k))*rdy*vh(1)   &
                      +vten(i,1,k) )
          enddo
          enddo
        endif
 
        if(nbc.eq.2.and.ibn.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,k)
          do k=1,nk
          do i=1,ni
            v3d(i,nj+1,k)=v3d(i,nj+1,k)+dts*( -radbcn(i,k)           &
                         *(v3d(i,nj+1,k)-v3d(i,nj  ,k))*rdy*vh(nj)   &
                         +vten(i,nj+1,k) )
          enddo
          enddo
        endif

        if(timestats.ge.1) time_sound=time_sound+mytime()

!-----

        IF(roflux.eq.1)THEN
          call restrict_openbc_sn(ruh,rmh,rho0,v3d)
        ENDIF

!-----
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
          r1 = (sigma(1)-sigmaf(1))*0.5*rdsf(1)
          do i=0,ni+1
            dum1(i,j,1) = ppd(i,j,1)-r1*(-3.0*ppd(i,j,1)+4.0*ppd(i,j,2)-ppd(i,j,3))
          enddo

          ! upper-most model level:
          r1 = (sigmaf(nk+1)-sigma(nk))*0.5*rdsf(nk)
          do i=0,ni+1
            dum1(i,j,nk+1) = ppd(i,j,nk)+r1*(3.0*ppd(i,j,nk)-4.0*ppd(i,j,nk-1)+ppd(i,j,nk-2))
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

!---------------------------------------------------
!  update implicit solver arrays

      tem1 = dts*rdz
      tem2 = -0.5*thlr*dts

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1)
    DO k=1,nk
      IF(terrain_flag)THEN
        ! with terrain:
        r1 = dts*rdsf(k)
        do j=1,nj
        do i=1,ni
          qk(i,j,k)=r1*gz(i,j)*ppterm(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))
          sk(i,j,k)=r1*gz(i,j)*thterm(i,j,k)*(th0(i,j,k)+th3d(i,j,k))
        enddo
        enddo
      ELSE
        ! without terrain:
        r1 = tem1*mh(1,1,k)
        do j=1,nj
        do i=1,ni
          qk(i,j,k)=r1*ppterm(i,j,k)*(pi0(1,1,k)+pp3d(i,j,k))
          sk(i,j,k)=r1*thterm(i,j,k)*(th0(1,1,k)+th3d(i,j,k))
        enddo
        enddo
      ENDIF
    ! store divh in dum1:
    IF(axisymm.eq.0)THEN
      IF(.not.terrain_flag)THEN
        ! Cartesian grid without terrain:
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = (u3d(i+1,j,k)-u3d(i,j,k))*rdx*uh(i)    &
                       +(v3d(i,j+1,k)-v3d(i,j,k))*rdy*vh(j)
        enddo
        enddo
      ELSE
        ! Cartesian grid with terrain:
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = gz(i,j)*(                                         &
                (u3d(i+1,j,k)*rgzu(i+1,j)-u3d(i,j,k)*rgzu(i,j))*rdx*uh(i) &
               +(v3d(i,j+1,k)*rgzv(i,j+1)-v3d(i,j,k)*rgzv(i,j))*rdy*vh(j) )
        enddo
        enddo
      ENDIF
    ELSE
        ! axisymmetric grid:
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = (xf(i+1)*u3d(i+1,j,k)-xf(i)*u3d(i,j,k))*rdx*uh(i)*rxh(i)
        enddo
        enddo
    ENDIF
        do j=1,nj
        do i=1,ni
          fpk(i,j,k)=dts*( ppten(i,j,k)                                     &
                 -ppterm(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))*dum1(i,j,k) )      &
              +beta*( w3d(i,j,k  )*(pk(i,j,k)*rf0(i,j,k  )+qk(i,j,k))       &
                     +w3d(i,j,k+1)*(pk(i,j,k)*rf0(i,j,k+1)-qk(i,j,k)) )
          ftk(i,j,k)=dts*( thten(i,j,k)                                     &
                 -thterm(i,j,k)*(th0(i,j,k)+th3d(i,j,k))*dum1(i,j,k) )      &
              +beta*( w3d(i,j,k  )*(tem2*rr0(i,j,k)*rf0(i,j,k  )+sk(i,j,k)) &
                     +w3d(i,j,k+1)*(tem2*rr0(i,j,k)*rf0(i,j,k+1)-sk(i,j,k)) )
        enddo
        enddo

      IF(k.gt.1)THEN
        do j=1,nj
        do i=1,ni
          fwk(i,j,k)=dts*wten(i,j,k)                                           &
             +beta*( bk(i,j,k)*(c2(i,j,k)*th3d(i,j,k)+c1(i,j,k)*th3d(i,j,k-1)) &
                    -mm(i,j,k)*(pp3d(i,j,k)-pp3d(i,j,k-1)) )
        enddo
        enddo
      ENDIF

    ENDDO

!---------------------------------------------
!  remaining piece of horizontal
!  divergence for terrain:

      IF(terrain_flag)THEN

        ! dum1 stores u3d at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
        do j=1,nj
          ! lowest model level:
          r1 = (sigma(1)-sigmaf(1))*0.5*rdsf(1)
          do i=1,ni+1
            dum1(i,j,1) = u3d(i,j,1)-r1*(-3.0*u3d(i,j,1)+4.0*u3d(i,j,2)-u3d(i,j,3))
          enddo

          ! upper-most model level:
          r1 = (sigmaf(nk+1)-sigma(nk))*0.5*rdsf(nk)
          do i=1,ni+1
            dum1(i,j,nk+1) = u3d(i,j,nk)+r1*(3.0*u3d(i,j,nk)-4.0*u3d(i,j,nk-1)+u3d(i,j,nk-2))
          enddo

          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=1,ni+1
            dum1(i,j,k) = r1*u3d(i,j,k-1)+r2*u3d(i,j,k)
          enddo
          enddo
        enddo

        ! dum2 stores v3d at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
        do j=1,nj+1
          ! lowest model level:
          r1 = (sigma(1)-sigmaf(1))*0.5*rdsf(1)
          do i=1,ni
            dum2(i,j,1) = v3d(i,j,1)-r1*(-3.0*v3d(i,j,1)+4.0*v3d(i,j,2)-v3d(i,j,3))
          enddo

          ! upper-most model level:
          r1 = (sigmaf(nk+1)-sigma(nk))*0.5*rdsf(nk)
          do i=1,ni
            dum2(i,j,nk+1) = v3d(i,j,nk)+r1*(3.0*v3d(i,j,nk)-4.0*v3d(i,j,nk-1)+v3d(i,j,nk-2))
          enddo

          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=1,ni
            dum2(i,j,k) = r1*v3d(i,j,k-1)+r2*v3d(i,j,k)
          enddo
          enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1)
        do j=1,nj
          do k=1,nk
          do i=1,ni
            r1 = ( ( gx(i,j,k+1)*(dum1(i,j,k+1)+dum1(i+1,j,k+1))      &
                    -gx(i,j,k  )*(dum1(i,j,k  )+dum1(i+1,j,k  )) )    &
                  +( gy(i,j,k+1)*(dum2(i,j,k+1)+dum2(i,j+1,k+1))      &
                    -gy(i,j,k  )*(dum2(i,j,k  )+dum2(i,j+1,k  )) )    &
                 )*0.5*rdsf(k)
            fpk(i,j,k)=fpk(i,j,k)-dts*ppterm(i,j,k)*(pi0(i,j,k)+pp3d(i,j,k))*r1
            ftk(i,j,k)=ftk(i,j,k)-dts*thterm(i,j,k)*(th0(i,j,k)+th3d(i,j,k))*r1
          enddo
          enddo
        enddo

      ENDIF

      if(terrain_flag) call bcwsfc(zh,zf,sigma,sigmaf,dzdx,dzdy,dum1,dum2,u3d,v3d,w3d)

!---------------------------------------------
!  the implicit solver ....

        tem1 = -0.5*thlr*dts

!$omp parallel do default(shared)   &
!$omp private(i,j,k,aa,bb,cc,dd,r1)
      DO j=1,nj

        !  Get boundary conditions for w
        if(terrain_flag)then
          ! with terrain:
          do i=1,ni
            ! w3d at level 1 is set above:
            dum1(i,j,1)=0.0
            dum2(i,j,1)=w3d(i,j,1)
          enddo
        else
          ! without terrain:
          do i=1,ni
            w3d(i,j,1)=0.0
            dum1(i,j,1)=0.0
            dum2(i,j,1)=w3d(i,j,1)
          enddo
        endif

        do k=2,nk
        do i=1,ni
          aa=asq*( mm(i,j,k)*(qk(i,j,k)-pk(i,j,k)*rf0(i,j,k+1))   &
                  +bk(i,j,k)*(tem1*rf0(i,j,k+1)*c2(i,j,k)*rr0(i,j,k)-sk(i,j,k)) )
          bb=1.0+asq*( mm(i,j,k)*( qk(i,j,k)+qk(i,j,k-1)                    &
                                  +pk(i,j,k)*rf0(i,j,k)-pk(i,j,k-1)*rf0(i,j,k) ) &
                      -bk(i,j,k)*(tem1*rf0(i,j,k)*(c2(i,j,k)*rr0(i,j,k)+c1(i,j,k)*rr0(i,j,k-1))   &
                                  +sk(i,j,k)-sk(i,j,k-1)) )
          cc=asq*( mm(i,j,k)*(pk(i,j,k-1)*rf0(i,j,k-1)+qk(i,j,k-1))   &
                  +bk(i,j,k)*(tem1*rf0(i,j,k-1)*c1(i,j,k)*rr0(i,j,k-1)+sk(i,j,k-1)) )
          dd=w3d(i,j,k)+fwk(i,j,k)                                                    &
              -alph*( mm(i,j,k)*(pp3d(i,j,k)+fpk(i,j,k)-pp3d(i,j,k-1)-fpk(i,j,k-1))   &
                     -bk(i,j,k)*( c2(i,j,k)*(th3d(i,j,k)+ftk(i,j,k))+c1(i,j,k)*(th3d(i,j,k-1)+ftk(i,j,k-1)) ) )
          r1 = 1.0/(bb-cc*dum1(i,j,k-1))
          dum1(i,j,k)=aa*r1
          dum2(i,j,k)=(dd+cc*dum2(i,j,k-1))*r1
        enddo
        enddo

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
          call bcwsfc(zh,zf,sigma,sigmaf,dzdx,dzdy,dum1,dum2,u3d,v3d,w3d)
        endif

#ifdef MPI
        if(n.eq.nloop)then
          call comm_3w_start(w3d,ww31,ww32,we31,we32,   &
                                 ws31,ws32,wn31,wn32,reqs_w)
        endif
#endif

!-------------------------

      tem2 = -0.5*thlr*dts
 
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      DO k=1,nk

        do j=1,nj
        do i=1,ni
          ppd(i,j,k)=pp3d(i,j,k)
          pp3d(i,j,k)=pp3d(i,j,k)+fpk(i,j,k)                          &
              +alph*( w3d(i,j,k  )*(pk(i,j,k)*rf0(i,j,k  )+qk(i,j,k)) &
                     +w3d(i,j,k+1)*(pk(i,j,k)*rf0(i,j,k+1)-qk(i,j,k)) )
          if(abs(pp3d(i,j,k)).lt.smeps) pp3d(i,j,k)=0.0
          ppd(i,j,k)=pp3d(i,j,k)+kdiv*( pp3d(i,j,k)-ppd(i,j,k)-dts*ppten(i,j,k) )
        enddo
        enddo

        do j=1,nj
        do i=1,ni
          th3d(i,j,k)=th3d(i,j,k)+ftk(i,j,k)                      &
              +alph*( w3d(i,j,k  )*(tem2*rr0(i,j,k)*rf0(i,j,k  )+sk(i,j,k)) &
                     +w3d(i,j,k+1)*(tem2*rr0(i,j,k)*rf0(i,j,k+1)-sk(i,j,k)) )
          if(abs(th3d(i,j,k)).lt.smeps) th3d(i,j,k)=0.0
        enddo
        enddo

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

      ENDDO

!--------------------------------------------------------------------

      return
      end


