!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcecd(u0,v0,u1,v1,w1,u,v,zh,u10,v10,w10,xland,znt,ust,cd,ch,cq)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(ib:ie,jb:je) :: u1,v1,w1
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,u
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,v
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: u10,v10,w10
      real, intent(in), dimension(ibl:iel,jbl:jel) :: xland
      real, intent(inout), dimension(ib:ie,jb:je) :: znt,ust,cd,ch,cq

      integer i,j,n,nmax
      real wsp,wlast,var,rznt

      real, parameter :: dcd1  =  1.0e-3
      real, parameter :: dcd2  =  2.4e-3
      real, parameter :: dwsp1 =  5.0
      real, parameter :: dwsp2 = 25.0

      real, parameter :: dfac = (dcd2-dcd1)/(dwsp2-dwsp1)

!-----------------------------------------------------------------------
!  Get Cd:

    nmax = 0

!$omp parallel do default(shared)   &
!$omp private(i,j,wlast,n,var,rznt)
    DO j=1,nj
    do i=1,ni
      IF(xland(i,j).gt.1.5)THEN
        ! water:  roughness length (z0) is a function of windspeed
        ! use last known z0 for first guess:
        rznt = 1.0/znt(i,j)
        var = alog(10.0*rznt)/alog(zh(i,j,1)*rznt)
        w10(i,j) = w1(i,j)*var
        wlast = -1.0
        n = 0
        do while( abs(w10(i,j)-wlast).gt.0.001 )
          n = n + 1
          wlast = w10(i,j)
          IF(cecd.eq.1)THEN
            ! constant value:
            cd(i,j) = max(1.0e-4,cnstcd)
          ELSEIF(cecd.eq.2)THEN
            ! Deacon's formula:  see Rotunno and Emanuel (1987, JAS, p. 547)
            cd(i,j) = 1.1e-3+(4.0e-5*w10(i,j))
          ELSEIF(cecd.eq.3)THEN
            ! based on Fairall et al (2003, JClim) at low wind speeds
            ! based on Donelan et al (2004, GRL) at high wind speeds
            cd(i,j) = dcd1+(w10(i,j)-dwsp1)*dfac
            cd(i,j) = min(cd(i,j),dcd2)
            cd(i,j) = max(cd(i,j),dcd1)
          ENDIF
          znt(i,j) = 10.0/exp(karman/sqrt(cd(i,j)))
          rznt = 1.0/znt(i,j)
          var = alog(10.0*rznt)/alog(zh(i,j,1)*rznt)
          w10(i,j) = w1(i,j)*var
          if(n.gt.10) print *,n,w10(i,j)
          if(n.gt.20) stop 12555
        enddo
!!!        nmax = max(nmax,n)
      ELSE
        ! land:  roughness length (z0) is specified
        IF(cecd.eq.1)THEN
          cd(i,j) = max(1.0e-4,cnstcd)
          znt(i,j) = 10.0/exp(karman/sqrt(cd(i,j)))
          rznt = 1.0/znt(i,j)
        ELSE
          rznt = 1.0/znt(i,j)
          cd(i,j) = ( karman/alog(10.0*rznt) )**2
        ENDIF
        var = alog(10.0*rznt)/alog(zh(i,j,1)*rznt)
        w10(i,j) = w1(i,j)*var
      ENDIF
      u10(i,j) = u1(i,j)*var
      v10(i,j) = v1(i,j)*var
!!!      ust(i,j) = sqrt(cd(i,j))*w10(i,j)
      ust(i,j) = w1(i,j)*karman/alog(zh(i,j,1)*rznt)
    enddo
!--------------------------------------------------
!  Get Ce:
    IF(isfcflx.eq.1)THEN
      do i=1,ni
        IF(xland(i,j).gt.1.5)THEN
          ! water:
          IF(cecd.eq.1)THEN
            ! constant value (from namelist.input):
            ch(i,j) = cnstce
            cq(i,j) = cnstce
          ELSEIF(cecd.eq.2)THEN
            ! Deacon's formula:  see Rotunno and Emanuel (1987, JAS, p. 547)
            ch(i,j) = 1.1e-3+(4.0e-5*w10(i,j))
            cq(i,j) = 1.1e-3+(4.0e-5*w10(i,j))
          ELSEIF(cecd.eq.3)THEN
            ! Constant, based on Drennan et al. (2007, JAS, p. 1103)
            ch(i,j) = 1.20e-3
            cq(i,j) = 1.20e-3
          ENDIF
        ELSE
          ! land ... just set Ce to Cd:
          ch(i,j) = cd(i,j)
          cq(i,j) = cd(i,j)
        ENDIF
      enddo
    ENDIF
!--------------------------------------------------

    ENDDO  ! enddo for j-loop

!-----------------------------------------------------------------------

!!!      print *,'  nmax = ',nmax

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine sfcflux(dt,ruh,xf,rvh,pi0s,ch,cq,zh,pi0,thv0,th0,u0,v0,tsk,thflux,qvflux,mavail,   &
                         rho,rf,u1,v1,w1,u,v,ppi,tha,qva,qsfc,psfc,u10,v10,w10,znt)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: rvh
      real, intent(in), dimension(ib:ie,jb:je) :: pi0s,ch,cq
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,pi0,thv0,th0
      real, intent(in), dimension(ib:ie,jb:je) :: tsk
      real, intent(inout), dimension(ib:ie,jb:je) :: thflux,qvflux
      real, intent(in), dimension(ibl:iel,jbl:jel) :: mavail
      real, intent(in), dimension(ib:ie,jb:je) :: u1,v1,w1
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,u
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,v
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rf,ppi,tha
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem) :: qva
      real*8, intent(inout) :: qsfc
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: psfc,u10,v10,w10
      real, intent(inout), dimension(ib:ie,jb:je) :: znt

      integer i,j
      real pisfc,qvsat
      real :: rhosfc,tem,p1,p2,p3
      real rslf
      real*8, dimension(nj) :: bud1

!-----------------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,p1,p2,p3,pisfc,qvsat,rhosfc)
    DO j=1,nj
      bud1(j)=0.0d0

      !  sensible heat flux:
      do i=1,ni
!!!        pisfc = pi0s(i,j) + ppi(i,j,1)
        ! cm1r17, 2nd-order extrapolation:
        p1 = pi0(i,j,1)+ppi(i,j,1)
        p2 = pi0(i,j,2)+ppi(i,j,2)
        p3 = pi0(i,j,3)+ppi(i,j,3)
        pisfc = cgs1*p1+cgs2*p2+cgs3*p3
        psfc(i,j) = p00*(pisfc**cpdrd)
        thflux(i,j)=ch(i,j)*w10(i,j)*(tsk(i,j)/pisfc-th0(i,j,1)-tha(i,j,1))
      enddo

      !  latent heat flux:
      IF(imoist.eq.1)THEN
        do i=1,ni
          qvsat=rslf(psfc(i,j),tsk(i,j))
          qvflux(i,j)=cq(i,j)*w10(i,j)*(qvsat-qva(i,j,1))*mavail(i,j)
          ! some budget calculations (only calculated if imoist=1):
          rhosfc=rf(i,j,1)
          if(axisymm.eq.1) rhosfc=rhosfc*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
          bud1(j)=bud1(j)+qvflux(i,j)*ruh(i)*rvh(j)*rhosfc
        enddo
      ENDIF
    ENDDO

    IF(imoist.eq.1)THEN
      tem = dt*dx*dy
      do j=1,nj
        qsfc=qsfc+bud1(j)*tem
      enddo
    ENDIF

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine sfcdiags(tsk,thflux,qvflux,cd,ch,cq,u1,v1,w1,             &
                          xland,psfc,qsfc,u10,v10,hfx,qfx,cda,znt,gz1oz0,  &
                          psim,psih,br,zol,mol,hpbl,dsxy,th2,t2,q2,        &
                          zs,zh,pi0s,pi0,th0,ppi,tha,rho,rf,qa,ua,va)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(ib:ie,jb:je) :: tsk,thflux,qvflux,   &
                                                  cd,ch,cq,u1,v1,w1
      real, intent(in), dimension(ibl:iel,jbl:jel) :: xland
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: psfc,qsfc,u10,v10,hfx,qfx, &
                                    cda,gz1oz0,psim,psih,br,zol,mol,hpbl,dsxy,th2,t2,q2
      real, intent(inout), dimension(ib:ie,jb:je) :: znt
      real, intent(in), dimension(itb:ite,jtb:jte) :: zs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je) :: pi0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,th0,ppi,tha,rho,rf
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va

      integer :: i,j
      real :: pisfc,thgb,thx,thvx,tskv,za,govrth,dthvdz,vconv,vsgd,dthvm,   &
              val,fluxc,wspd
      real :: rslf

      REAL    , PARAMETER ::  VCONVC=1.
      REAL    , PARAMETER ::  CZO=0.0185
      REAL    , PARAMETER ::  OZO=1.59E-5
      REAL    , PARAMETER ::  EP1 = rv/rd - 1.0
      REAL,     PARAMETER ::  PRT=1.

      ! surface layer diagnostics:

!$omp parallel do default(shared)   &
!$omp private(i,j,pisfc,thgb,thx,thvx,tskv,za,govrth,dthvdz,vconv,vsgd,   &
!$omp dthvm,val,fluxc,wspd)
      do j=1,nj
      do i=1,ni
!!!        pisfc = pi0s(i,j) + ppi(i,j,1)
!!!        psfc(i,j) = p00*(pisfc**cpdrd)
        pisfc = (psfc(i,j)*rp00)**rovcp
        thgb = tsk(i,j)/pisfc
        thx = th0(i,j,1)+tha(i,j,1)
        thvx = thx*(1.+EP1*qa(i,j,1,nqv))
        qsfc(i,j) = rslf(psfc(i,j),tsk(i,j))
        tskv = thgb*(1.0+ep1*qsfc(i,j))
        za = zh(i,j,1)
        if( terrain_flag ) za = za - zs(i,j)
        govrth = g/thx
        ! use old znt (as in sfclay):
        gz1oz0(i,j) = alog(za/znt(i,j))
        DTHVDZ = THVX-TSKV
        if (xland(i,j).lt.1.5) then
          ! land:
          fluxc = max(thflux(i,j) + ep1*tskv*qvflux(i,j),0.)
          VCONV = vconvc*(g/tsk(i,j)*hpbl(i,j)*fluxc)**.33
        else
          ! ocean:
          IF(-DTHVDZ.GE.0)THEN
            DTHVM=-DTHVDZ
          ELSE
            DTHVM=0.
          ENDIF
          VCONV = 2.*SQRT(DTHVM)
        endif
! Mahrt and Sun low-res correction
        VSGD = 0.32 * (max(dsxy(i,j)/5000.-1.,0.))**.33
        wspd = sqrt( w1(i,j)*w1(i,j) + vconv*vconv + vsgd*vsgd )
        wspd = max(0.1,wspd)
        br(i,j) = govrth*za*DTHVDZ/(wspd**2)
        hfx(i,j) = thflux(i,j)*cp*rf(i,j,1)
        qfx(i,j) = qvflux(i,j)*rf(i,j,1)
        cda(i,j) = cd(i,j)
        ! new ust/znt:
        ! TO PREVENT OSCILLATIONS AVERAGE WITH OLD VALUE
!!!        ust = 0.5*( ust + sqrt(cda(i,j)*wspd*wspd) )
!!!        if( (xland(i,j)-1.5).ge.0.0 )then
!!!          znt(i,j) = CZO*ust*ust/G+OZO
!!!        endif
        ! impose neutral sfc layer:
        psim(i,j) = 0.0
        psih(i,j) = 0.0
        zol(i,j) = 0.0
        mol(i,j) = 0.0
!!!        val = alog(10.0/znt(i,j))/alog(za/znt(i,j))
!!!        u10(i,j) = 0.5*(ua(i,j,1)+ua(i+1,j,1))*val
!!!        v10(i,j) = 0.5*(va(i,j,1)+va(i,j+1,1))*val
        val = alog(2.0/znt(i,j))/alog(za/znt(i,j))
        th2(i,j) = thgb+(thx-thgb)*val
        q2(i,j) = qsfc(i,j)+(qa(i,j,1,nqv)-qsfc(i,j))*val
        t2(i,j) = th2(i,j)*pisfc
      enddo
      enddo

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()
      end subroutine sfcdiags


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine gethpbl(zh,th0,tha,qa,hpbl)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,th0,tha
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: hpbl

      integer :: i,j,kk
      real :: thx,thvx,thv,thvlast,thcrit

      REAL    , PARAMETER ::  EP1 = rv/rd - 1.0

      ! (NEEDED BY SFCLAY ... THIS IS A ROUGH ESTIMATE ONLY)
      ! (ONLY NEEDED WHEN IPBL=0)
      ! (USE WITH CAUTION)
      ! extraordinarily simple calculation:  define pbl depth as 
      ! level where thv is first greater than thv at lowest model level
      ! 110104:  add 0.5 K, for the sake of slightly stable PBLs

!$omp parallel do default(shared)   &
!$omp private(i,j,kk,thx,thvx,thv,thvlast,thcrit)
      do j=1,nj
      do i=1,ni
        hpbl(i,j) = 0.0
        kk = 1
        thx = th0(i,j,1)+tha(i,j,1)
        thvx = thx*(1.+EP1*qa(i,j,1,nqv))
        thvlast = thvx
        thcrit = thvx+0.5
        do while( hpbl(i,j).lt.1.0e-12 .and. kk.lt.nk )
          kk = kk + 1
          thv = (th0(i,j,kk)+tha(i,j,kk))*(1.0+EP1*qa(i,j,kk,nqv))
          if( thv.ge.thcrit )then
            hpbl(i,j) = zh(i,j,kk-1)+(zh(i,j,kk)-zh(i,j,kk-1))   &
                                    *(thcrit-thvlast)/(thv-thvlast)
          endif
          thvlast = thv
        enddo
        if( kk.gt.(nk-1) .or. hpbl(i,j).lt.1.0e-12 ) hpbl(i,j) = 0.0
      enddo
      enddo

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()
      end subroutine gethpbl


