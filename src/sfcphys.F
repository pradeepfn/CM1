!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine getcecd(cdu,cdv,ce,u0,v0,u1,v1,w1,u,v,zh,u10,v10,w10,xland,znt,ust,cd, &
                         west,newwest,east,neweast,south,newsouth,north,newnorth,reqs_s)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(inout), dimension(ib:ie,jb:je) :: cdu,cdv,ce
      real, intent(inout), dimension(ib:ie,jb:je) :: u1,v1,w1
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,u
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,v
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: u10,v10,w10,xland,znt,ust,cd
      integer, intent(inout), dimension(rmp) :: reqs_s
      real, intent(inout), dimension(cmp,jmp) :: west,newwest,east,neweast
      real, intent(inout), dimension(imp,cmp) :: south,newsouth,north,newnorth

      integer i,j,n,nmax
      real wsp,wlast,var,rznt

      real, parameter :: dcd1  =  1.0e-3
      real, parameter :: dcd2  =  2.4e-3
      real, parameter :: dwsp1 =  5.0
      real, parameter :: dwsp2 = 25.0

      real, parameter :: dfac = (dcd2-dcd1)/(dwsp2-dwsp1)

!-----------------------------------------------------------------------
! for pertflx=1, the base-state is not included in calculation of wind speed

    ! u1 is u at lowest model level
    ! v1 is v at lowest model level
    ! w1 is total wind speed at lowest model level
    ! (all defined at the scalar point of the staggered grid)

!$omp parallel do default(shared)   &
!$omp private(i,j)
    do j=1,nj
      IF(pertflx.eq.1)THEN
        do i=1,ni
          u1(i,j) = 0.5*( u(i,j,1)-u0(i,j,1) + u(i+1,j,1)-u0(i+1,j,1) )
          v1(i,j) = 0.5*( v(i,j,1)-v0(i,j,1) + v(i,j+1,1)-v0(i,j+1,1) )
        enddo
      ELSE
        do i=1,ni
          u1(i,j) = 0.5*( u(i,j,1) + u(i+1,j,1) )
          v1(i,j) = 0.5*( v(i,j,1) + v(i,j+1,1) )
        enddo
      ENDIF
      do i=1,ni
        w1(i,j)=sqrt(u1(i,j)**2+v1(i,j)**2)
      enddo
    enddo

!-----------------------------------------------------------------------
!  Get Cd:

    nmax = 0

  IF(idrag.eq.1.or.isfcflx.eq.1)THEN

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
      ust(i,j) = sqrt(cd(i,j)*w10(i,j)*w10(i,j))
    enddo
!--------------------------------------------------
!  Get Ce:
    IF(isfcflx.eq.1)THEN
      do i=1,ni
        IF(xland(i,j).gt.1.5)THEN
          ! water:
          IF(cecd.eq.1)THEN
            ! constant value (from namelist.input):
            ce(i,j) = cnstce
          ELSEIF(cecd.eq.2)THEN
            ! Deacon's formula:  see Rotunno and Emanuel (1987, JAS, p. 547)
            ce(i,j) = 1.1e-3+(4.0e-5*w10(i,j))
          ELSEIF(cecd.eq.3)THEN
            ! Constant, based on Drennan et al. (2007, JAS, p. 1103)
            ce(i,j) = 1.20e-3
          ENDIF
        ELSE
          ! land ... just set Ce to Cd:
          ce(i,j) = cd(i,j)
        ENDIF
      enddo
    ENDIF
!--------------------------------------------------

    ENDDO  ! enddo for j-loop

  ENDIF

        call bc2d(cd)
        call bc2d(w10)
        call bc2d(u10)
        call bc2d(v10)
#ifdef MPI
        call comm_2d_start(cd,west,newwest,east,neweast,   &
                              south,newsouth,north,newnorth,reqs_s)
        call comm_2dew_end(cd,west,newwest,east,neweast,reqs_s)
        call comm_2dns_end(cd,south,newsouth,north,newnorth,reqs_s)
        call bcs2_2d(cd)
        call comm_2d_corner(cd)
        call comm_2d_start(w10,west,newwest,east,neweast,   &
                              south,newsouth,north,newnorth,reqs_s)
        call comm_2dew_end(w10,west,newwest,east,neweast,reqs_s)
        call comm_2dns_end(w10,south,newsouth,north,newnorth,reqs_s)
        call bcs2_2d(w10)
        call comm_2d_corner(w10)
        call comm_2d_start(u10,west,newwest,east,neweast,   &
                              south,newsouth,north,newnorth,reqs_s)
        call comm_2dew_end(u10,west,newwest,east,neweast,reqs_s)
        call comm_2dns_end(u10,south,newsouth,north,newnorth,reqs_s)
        call bcs2_2d(u10)
        call comm_2d_corner(u10)
        call comm_2d_start(v10,west,newwest,east,neweast,   &
                              south,newsouth,north,newnorth,reqs_s)
        call comm_2dew_end(v10,west,newwest,east,neweast,reqs_s)
        call comm_2dns_end(v10,south,newsouth,north,newnorth,reqs_s)
        call bcs2_2d(v10)
        call comm_2d_corner(v10)
#endif

!-----------------------------------------------------------------------
!  Cd on u,v pts:

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1
        cdu(i,j) = 0.5*(cd(i-1,j)+cd(i,j))
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        cdv(i,j) = 0.5*(cd(i,j-1)+cd(i,j))
      enddo
      enddo

!-----------------------------------------------------------------------

!!!      print *,'  nmax = ',nmax

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine sfcdrag(cdu,cdv,u0,v0,u1,v1,t13,t23,u,v,u10,v10,wspd)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(ib:ie,jb:je) :: cdu,cdv
      real, intent(in), dimension(ib:ie,jb:je) :: u1,v1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t13,t23
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0,u
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0,v
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: u10,v10,wspd

      integer i,j

!-----------------------------------------------------------------------
!  x component:

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1
        t13(i,j,1) = 0.25*cdu(i,j)*(u10(i-1,j)+u10(i,j))*(wspd(i-1,j)+wspd(i,j))
      enddo
      enddo

!-----------------------------------------------------------------------
!  y component:

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        t23(i,j,1) = 0.25*cdv(i,j)*(v10(i,j-1)+v10(i,j))*(wspd(i,j-1)+wspd(i,j))
      enddo
      enddo

!-----------------------------------------------------------------------

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine sfcflux(dt,ruh,xf,rvh,pi0s,ce,zh,pi0,thv0,th0,u0,v0,tsk,thflux,qvflux,mavail,   &
                         rho,rf,u1,v1,w1,u,v,ppi,tha,qva,qsfc,psfc,u10,v10,wspd,znt)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: rvh
      real, intent(in), dimension(ib:ie,jb:je) :: pi0s,ce
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
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: psfc,u10,v10,wspd,znt

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
        ! 2nd-order extrapolation:
        p1 = pi0(i,j,1)+ppi(i,j,1)
        p2 = pi0(i,j,2)+ppi(i,j,2)
        p3 = pi0(i,j,3)+ppi(i,j,3)
        pisfc=(1.75*p1-p2+0.25*p3)
        psfc(i,j) = p00*(pisfc**cpdrd)
        thflux(i,j)=ce(i,j)*wspd(i,j)*(tsk(i,j)/pisfc-th0(i,j,1)-tha(i,j,1))
      enddo

      !  latent heat flux:
      IF(imoist.eq.1)THEN
        do i=1,ni
          qvsat=rslf(psfc(i,j),tsk(i,j))
          qvflux(i,j)=ce(i,j)*wspd(i,j)*(qvsat-qva(i,j,1))*mavail(i,j)
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


      subroutine sfcdiags(tsk,thflux,qvflux,cdu,cdv,ce,u1,v1,w1,           &
                          xland,psfc,qsfc,u10,v10,hfx,qfx,cda,znt,gz1oz0,  &
                          psim,psih,br,zol,mol,hpbl,dsxy,th2,t2,q2,        &
                          zs,zh,pi0s,pi0,th0,ppi,tha,rho,rf,qa,ua,va)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in), dimension(ib:ie,jb:je) :: tsk,thflux,qvflux,   &
                                                  cdu,cdv,ce,u1,v1,w1
      real, intent(in), dimension(ibl:iel,jbl:jel) :: xland
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: psfc,qsfc,u10,v10,hfx,qfx, &
                                    cda,znt,gz1oz0,psim,psih,br,zol,mol,hpbl,dsxy,th2,t2,q2
      real, intent(in), dimension(itb:ite,jtb:jte) :: zs
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je) :: pi0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,th0,ppi,tha,rho,rf
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va

      integer :: i,j
      real :: pisfc,thgb,thx,thvx,tskv,za,govrth,dthvdz,vconv,vsgd,dthvm,   &
              val,fluxc,wspd,ust
      real :: rslf

      REAL    , PARAMETER ::  VCONVC=1.
      REAL    , PARAMETER ::  CZO=0.0185
      REAL    , PARAMETER ::  OZO=1.59E-5
      REAL    , PARAMETER ::  EP1 = rv/rd - 1.0
      REAL,     PARAMETER ::  PRT=1.

      ! surface layer diagnostics:

!$omp parallel do default(shared)   &
!$omp private(i,j,pisfc,thgb,thx,thvx,tskv,za,govrth,dthvdz,vconv,vsgd,   &
!$omp dthvm,val,fluxc,wspd,ust)
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
        cda(i,j) = 0.25*( (cdu(i,j)+cdu(i+1,j))+(cdv(i,j)+cdv(i,j+1)) )
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


