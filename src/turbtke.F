

      subroutine turbtke(dt,dodrag,dosfcflx,ruh,rvh,rmh,zf,mf,rmf,th0,thflux,qvflux,rth0s,    &
                         nm,defsq,defh,tk,lenscl,lenh,grdscl,rgrdscl,   &
                         kmh,kmv,khh,khv,tkea,tketen,t13,t23,ua,va,dissten, &
                         nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                         khcw1,khcw2,khce1,khce2,khcs1,khcs2,khcn1,khcn2, &
                         khdw1,khdw2,khde1,khde2,khds1,khds2,khdn1,khdn2, &
                         kvcw1,kvcw2,kvce1,kvce2,kvcs1,kvcs2,kvcn1,kvcn2, &
                         kvdw1,kvdw2,kvde1,kvde2,kvds1,kvds2,kvdn1,kvdn2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real :: dt
      logical, intent(in) :: dodrag,dosfcflx
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf,rmf
      real, dimension(ib:ie,jb:je,kb:ke) :: th0
      real, dimension(ib:ie,jb:je) :: thflux,qvflux,rth0s
      real, dimension(ib:ie,jb:je,kb:ke) :: nm,defsq,defh,tk,   &
                                            lenscl,lenh,grdscl,rgrdscl
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tketen
      real, dimension(ib:ie,jb:je,kb:ke) :: t13,t23
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke) :: dissten
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmt) :: khcw1,khcw2,khce1,khce2
      real, intent(inout), dimension(imp,kmt) :: khcs1,khcs2,khcn1,khcn2
      real, intent(inout), dimension(jmp,kmt) :: khdw1,khdw2,khde1,khde2
      real, intent(inout), dimension(imp,kmt) :: khds1,khds2,khdn1,khdn2
      real, intent(inout), dimension(jmp,kmt) :: kvcw1,kvcw2,kvce1,kvce2
      real, intent(inout), dimension(imp,kmt) :: kvcs1,kvcs2,kvcn1,kvcn2
      real, intent(inout), dimension(jmp,kmt) :: kvdw1,kvdw2,kvde1,kvde2
      real, intent(inout), dimension(imp,kmt) :: kvds1,kvds2,kvdn1,kvdn2

!----------------------------------------

      integer i,j,k
      real prinv,tem1,tem2,tem3

#ifdef MPI
      integer reqs_khc(8)
      integer reqs_kvc(8)
      integer reqs_khd(8)
      integer reqs_kvd(8)
#endif

!------------------------------------------------------------------

      tem1 = 0.125*dx*dx/dt
      tem2 = 0.125*dy*dy/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,k,prinv)
    DO k=1,nk+1

    !  get grid scale
    IF(tconfig.eq.1)THEN
      ! single length scale:  appropriate if dx,dy are nearly the same as dz
      do j=1,nj
      do i=1,ni
        grdscl(i,j,k)=(dx*ruh(i)*dy*rvh(j)*dz*rmf(i,j,k))**0.33333333
        rgrdscl(i,j,k)=1.0/grdscl(i,j,k)
      enddo
      enddo
    ELSEIF(tconfig.eq.2)THEN
      ! two length scales:  one for horizontal, one for vertical
      do j=1,nj
      do i=1,ni
        lenh(i,j,k)=sqrt(dx*ruh(i)*dy*rvh(j))
      enddo
      enddo
      do j=1,nj
      do i=1,ni
        grdscl(i,j,k)=dz*rmf(i,j,k)
        rgrdscl(i,j,k)=1.0/grdscl(i,j,k)
      enddo
      enddo
    ENDIF

      ! Get turbulence length scale
      do j=1,nj
      do i=1,ni
        tk(i,j,k)=max(tkea(i,j,k),1.0e-6)
        lenscl(i,j,k)=grdscl(i,j,k)
        if(nm(i,j,k).gt.1.0e-6)then
          lenscl(i,j,k)=0.8165*sqrt(tk(i,j,k)/nm(i,j,k))
          lenscl(i,j,k)=min(lenscl(i,j,k),grdscl(i,j,k))
          lenscl(i,j,k)=max(lenscl(i,j,k),1.0e-6*grdscl(i,j,k))
        endif 
      enddo
      enddo

    !  Get km, kh
    IF(tconfig.eq.1)THEN

      do j=1,nj
      do i=1,ni
        kmh(i,j,k)=0.10*sqrt(tk(i,j,k))*lenscl(i,j,k)
        kmv(i,j,k)=kmh(i,j,k)
        prinv=3.00
        if(nm(i,j,k).gt.1.0e-6)then
          prinv=min(1.0+2.00*lenscl(i,j,k)*rgrdscl(i,j,k),3.00)
        endif
        khh(i,j,k)=kmh(i,j,k)*prinv
        khv(i,j,k)=khh(i,j,k)
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN

      do j=1,nj
      do i=1,ni
        kmh(i,j,k)=0.10*sqrt(tk(i,j,k))*lenh(i,j,k)
        kmv(i,j,k)=0.10*sqrt(tk(i,j,k))*lenscl(i,j,k)
        prinv=3.00
        if(nm(i,j,k).gt.1.0e-6)then
          prinv=min(1.0+2.00*lenscl(i,j,k)*rgrdscl(i,j,k),3.00)
        endif
        khh(i,j,k)=kmh(i,j,k)*prinv
        khv(i,j,k)=kmv(i,j,k)*prinv
      enddo
      enddo

    ENDIF


    IF(k.ge.2)THEN
      !  Buoyancy Term
      do j=1,nj
      do i=1,ni
        tketen(i,j,k)=tketen(i,j,k)-khv(i,j,k)*nm(i,j,k)
      enddo
      enddo
      ! Shear term 
      IF(tconfig.eq.1)THEN
        do j=1,nj
        do i=1,ni
          tketen(i,j,k)=tketen(i,j,k)+kmv(i,j,k)*(defsq(i,j,k)+defh(i,j,k))
        enddo
        enddo
      ELSEIF(tconfig.eq.2)THEN
        do j=1,nj
        do i=1,ni
          tketen(i,j,k)=tketen(i,j,k)+kmv(i,j,k)*defsq(i,j,k)   &
                                     +kmh(i,j,k)*defh(i,j,k)
        enddo
        enddo
      ENDIF
    ELSE
      ! k=1 only:
      IF(dosfcflx)THEN
        do j=1,nj
        do i=1,ni
          tketen(i,j,1)=tketen(i,j,1)   &
                       +g*( thflux(i,j)*rth0s(i,j)+repsm1*qvflux(i,j) )
        enddo
        enddo
      ENDIF
      IF(bcturbu.eq.3)THEN
        do j=1,nj
        do i=1,ni
          tketen(i,j,1)=tketen(i,j,1)+(              &
                        ( t13(i  ,j,1)*ua(i  ,j,1)     &
                         +t13(i+1,j,1)*ua(i+1,j,1) )   &
                      + ( t23(i,j  ,1)*va(i,j  ,1)     &
                         +t23(i,j+1,1)*va(i,j+1,1) )   &
                                      )*rdz*mf(i,j,1)
        enddo
        enddo
      ENDIF
    ENDIF

      !  Dissipation Term
      do j=1,nj
      do i=1,ni
        dissten(i,j,k) = (0.191+0.796*lenscl(i,j,k)*rgrdscl(i,j,k))*   &
                      tk(i,j,k)*sqrt(tk(i,j,k))/lenscl(i,j,k)
        tketen(i,j,k)=tketen(i,j,k) - dissten(i,j,k)
      enddo
      enddo

      !  limit for numerical stability:
      do j=1,nj
      do i=1,ni
        kmh(i,j,k) = min( kmh(i,j,k) , tem1*ruh(i)*ruh(i) , tem2*rvh(j)*rvh(j) )
        khh(i,j,k) = min( khh(i,j,k) , tem1*ruh(i)*ruh(i) , tem2*rvh(j)*rvh(j) )
      enddo
      enddo

    ENDDO

      if(timestats.ge.1) time_turb=time_turb+mytime()

!------------------------------------------------------------
! Set values at boundaries

      call bct(kmh)
      call bct(kmv)
#ifdef MPI
      call comm_1t_start(kmh,khcw1,khcw2,khce1,khce2,   &
                             khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call comm_1t_start(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                             kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
#endif
      call bct(khh)
      call bct(khv)
#ifdef MPI
      call comm_1t_start(khh,khdw1,khdw2,khde1,khde2,   &
                             khds1,khds2,khdn1,khdn2,reqs_khd)
      call comm_1t_start(khv,kvdw1,kvdw2,kvde1,kvde2,   &
                             kvds1,kvds2,kvdn1,kvdn2,reqs_kvd)
!--------------------------------------------------------------
      call comm_1t_end(kmh,khcw1,khcw2,khce1,khce2,   &
                           khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call bct2(kmh)
      call getcornert(kmh,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call comm_1t_end(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                           kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call bct2(kmv)
      call getcornert(kmv,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call comm_1t_end(khh,khdw1,khdw2,khde1,khde2,   &
                           khds1,khds2,khdn1,khdn2,reqs_khd)
      call comm_1t_end(khv,kvdw1,kvdw2,kvde1,kvde2,   &
                           kvds1,kvds2,kvdn1,kvdn2,reqs_kvd)
#endif
!--------------------------------------------------------------
!  finished
      
      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbsmag(dt,dodrag,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,  &
                          nm,defsq,defh,lenscl,grdscl,lenh,                &
                          kmh,kmv,khh,khv,t13,t23,ua,va,dissten,           &
                          nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                          khcw1,khcw2,khce1,khce2,khcs1,khcs2,khcn1,khcn2, &
                          kvcw1,kvcw2,kvce1,kvce2,kvcs1,kvcs2,kvcn1,kvcn2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real :: dt
      logical, intent(in) :: dodrag,dosfcflx
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, dimension(ib:ie,jb:je,kb:ke) :: th0
      real, dimension(ib:ie,jb:je) :: thflux,qvflux,rth0s
      real, dimension(ib:ie,jb:je,kb:ke) :: nm,defsq,defh,lenscl,grdscl,lenh
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ib:ie,jb:je,kb:ke) :: t13,t23
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke) :: dissten
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmt) :: khcw1,khcw2,khce1,khce2
      real, intent(inout), dimension(imp,kmt) :: khcs1,khcs2,khcn1,khcn2
      real, intent(inout), dimension(jmp,kmt) :: kvcw1,kvcw2,kvce1,kvce2
      real, intent(inout), dimension(imp,kmt) :: kvcs1,kvcs2,kvcn1,kvcn2

      integer i,j,k
      real :: tem,tem1,tem2,temx,temy

#ifdef MPI
      integer reqs_khc(8)
      integer reqs_kvc(8)
#endif

      real, parameter :: cs      = 0.18
      real, parameter :: prandtl = 1.0/3.00
      real, parameter :: prinv   = 1.0/prandtl

!!!      real, parameter :: c_m = 0.0856
!!!      real, parameter :: c_h = 0.214
!!!      real, parameter :: c_l = 0.816
!!!      real, parameter :: ce1 = 0.191
!!!      real, parameter :: ce2 = 0.654
!!!      real, parameter :: ric = 0.23

!-----------------------------------------------------------------------

      temx = 0.125*dx*dx/dt
      temy = 0.125*dy*dy/dt

!$omp parallel do default(shared)   &
!$omp private(i,j,tem,tem1,tem2)
    DO j=1,nj

    IF(tconfig.eq.1)THEN
      ! single length scale:  appropriate if dx,dy are nearly the same as dz

      do k=1,nk+1
      do i=1,ni
        grdscl(i,j,k)=(dx*ruh(i)*dy*rvh(j)*dz*rmf(i,j,k))**0.33333333
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN
      ! two length scales:  one for horizontal, one for vertical

      do i=1,ni
        tem=sqrt(dx*ruh(i)*dy*rvh(j))
        do k=1,nk+1
          lenh(i,j,k)=tem
        enddo
      enddo

      do k=1,nk+1
      do i=1,ni
        grdscl(i,j,k)=dz*rmf(i,j,k)
      enddo
      enddo

    ENDIF

!-----------------------------------------------------------------------
!  Interior points:

    IF(tconfig.eq.1)THEN

      do k=2,nk
      do i=1,ni
        kmh(i,j,k)=((cs*grdscl(i,j,k))**2)     &
                 *sqrt( max(defsq(i,j,k)+defh(i,j,k)-nm(i,j,k)*prinv,0.0) )
        kmh(i,j,k) = min( kmh(i,j,k) , temx*ruh(i)*ruh(i)   &
                                     , temy*rvh(j)*rvh(j) )
        kmv(i,j,k)=kmh(i,j,k)
      enddo
      enddo

    ELSEIF(tconfig.eq.2)THEN

      do k=2,nk
      do i=1,ni
        kmh(i,j,k)=((cs*lenh(i,j,k))**2)     &
                 *sqrt( max(defh(i,j,k),0.0) )
        kmh(i,j,k) = min( kmh(i,j,k) , temx*ruh(i)*ruh(i)   &
                                     , temy*rvh(j)*rvh(j) )
        kmv(i,j,k)=((cs*grdscl(i,j,k))**2)     &
                 *sqrt( max(defsq(i,j,k)-nm(i,j,k)*prinv,0.0) )
      enddo
      enddo

    ENDIF

!------------------------------------------------------------
!  Surface:

    IF(bcturbu.eq.1.or.bcturbu.eq.2)THEN

      do i=1,ni
        tem1 = 0.0
      IF(dosfcflx)THEN
        tem2 = g*( thflux(i,j)*rth0s(i,j)+repsm1*qvflux(i,j) )
      ELSE
        tem2 = 0.0
      ENDIF
        kmv(i,j,1) = ( ((cs*grdscl(i,j,1))**4)*max(tem1+tem2,0.0) )**0.33333333
      enddo

    ELSEIF(bcturbu.eq.3)THEN

      do i=1,ni
        tem1 = ( ( t13(i  ,j,1)*ua(i  ,j,1)     &
                  +t13(i+1,j,1)*ua(i+1,j,1) )   &
               + ( t23(i,j  ,1)*va(i,j  ,1)     &
                  +t23(i,j+1,1)*va(i,j+1,1) )   &
               )*rdz*mf(i,j,1)
      IF(dosfcflx)THEN
        tem2 = g*( thflux(i,j)*rth0s(i,j)+repsm1*qvflux(i,j) )
      ELSE
        tem2 = 0.0
      ENDIF
        kmv(i,j,1) = ( ((cs*grdscl(i,j,1))**4)*max(tem1+tem2,0.0) )**0.33333333
      enddo

    ENDIF

!------------------------------------------------------------
!  Top of domain:

      do i=1,ni
        kmh(i,j,1) = kmh(i,j,2)
        kmh(i,j,nk+1) = 0.0
        kmv(i,j,nk+1) = 0.0
      enddo

    ENDDO

!--------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()
      call bct(kmh)
      call bct(kmv)
#ifdef MPI
      call comm_1t_start(kmh,khcw1,khcw2,khce1,khce2,   &
                             khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call comm_1t_start(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                             kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call comm_1t_end(kmh,khcw1,khcw2,khce1,khce2,   &
                           khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call bct2(kmh)
      call getcornert(kmh,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call comm_1t_end(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                           kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call bct2(kmv)
      call getcornert(kmv,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
#endif

!--------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    do j=0,nj+1

      do k=2,nk+2
      do i=0,ni+1
        khh(i,j,k)=kmh(i,j,k)*prinv
        ! limit for numerical stability:
        khh(i,j,k) = min( khh(i,j,k) , temx*ruh(i)*ruh(i)   &
                                     , temy*rvh(j)*rvh(j) )
        khv(i,j,k)=kmv(i,j,k)*prinv
      enddo
      enddo

      do i=0,ni+1
        khh(i,j,0)=kmh(i,j,0)
        khh(i,j,1)=kmh(i,j,1)
        khv(i,j,0)=kmv(i,j,0)
        khv(i,j,1)=kmv(i,j,1)
      enddo

    IF( j.ge.1 .and. j.le.nj )THEN
    IF( tconfig.eq.1 )THEN
      do k=1,nk+1
      do i=1,ni
        dissten(i,j,k) = (kmv(i,j,k)**3)/((cs*grdscl(i,j,k))**4)
      enddo
      enddo
    ELSE
      do k=1,nk+1
      do i=1,ni
        dissten(i,j,k) = (kmv(i,j,k)**3)/((cs*grdscl(i,j,k))**4)    &
                        +(kmh(i,j,k)**3)/((cs*lenh(i,j,k))**4)
      enddo
      enddo
    ENDIF
    ENDIF

    enddo

!--------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbparam(nstep,zf,dt,dodrag,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s, &
                          nm,defsq,defh,lvz,kmh,kmv,khh,khv,t13,t23,ua,va,dissten,znt, &
                          nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                         &
                          khcw1,khcw2,khce1,khce2,khcs1,khcs2,khcn1,khcn2,         &
                          kvcw1,kvcw2,kvce1,kvce2,kvcs1,kvcs2,kvcn1,kvcn2)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      integer, intent(in) :: nstep
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf
      real, intent(in) :: dt
      logical, intent(in) :: dodrag,dosfcflx
      real, dimension(ib:ie) :: ruh
      real, dimension(jb:je) :: rvh
      real, dimension(ib:ie,jb:je,kb:ke) :: rmh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, dimension(ib:ie,jb:je,kb:ke) :: th0
      real, dimension(ib:ie,jb:je) :: thflux,qvflux,rth0s
      real, dimension(ib:ie,jb:je,kb:ke) :: nm,defsq,defh,lvz
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, dimension(ib:ie,jb:je,kb:ke) :: t13,t23
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, dimension(ib:ie,jb:je,kb:ke) :: dissten
      real, intent(in), dimension(ibl:iel,jbl:jel) :: znt
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmt) :: khcw1,khcw2,khce1,khce2
      real, intent(inout), dimension(imp,kmt) :: khcs1,khcs2,khcn1,khcn2
      real, intent(inout), dimension(jmp,kmt) :: kvcw1,kvcw2,kvce1,kvce2
      real, intent(inout), dimension(imp,kmt) :: kvcs1,kvcs2,kvcn1,kvcn2

      integer i,j,k
      real :: tem,tem1,tem2,temx,temy

#ifdef MPI
      integer reqs_khc(8)
      integer reqs_kvc(8)
#endif

      real, parameter :: prandtl = 1.0
      real, parameter :: prinv   = 1.0/prandtl

!--------------------------------------------------------------
!  Smagorinsky-type scheme for parameterized turbulence:
!--------------------------------------------------------------
!  Interior:

    tem = 1.0/(1.0e-6+l_inf)
    if(ny.eq.1)then
      temx =  0.250*dx*dx/dt
      temy = 1000.0*dy*dy/dt
    elseif(nx.eq.1)then
      temx = 1000.0*dx*dx/dt
      temy =  0.250*dy*dy/dt
    else
      temx =  0.125*dx*dx/dt
      temy =  0.125*dy*dy/dt
    endif

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tem1,tem2)
    DO j=1,nj

    do k=2,nk
    do i=1,ni
      kmh(i,j,k)=(l_h**2)*sqrt( defh(i,j,k) )
      kmh(i,j,k) = min( kmh(i,j,k) , temx*ruh(i)*ruh(i) , temy*rvh(j)*rvh(j) )
      lvz(i,j,k)=karman*zf(i,j,k)/(1.0+karman*zf(i,j,k)*tem)
      kmv(i,j,k)=(lvz(i,j,k)**2)*sqrt( max(defsq(i,j,k)-nm(i,j,k)*prinv,0.0) )
    enddo
    enddo

  IF(sfcmodel.eq.0)THEN
    do i=1,ni
      ! temp roughness length at surface:
      lvz(i,j,1) = 0.001
      lvz(i,j,nk+1)=karman*zf(i,j,nk+1)/(1.0+karman*zf(i,j,nk+1)*tem)
    enddo
  ELSE
    do i=1,ni
      ! use roughness length at surface:
      lvz(i,j,1) = znt(i,j)
      lvz(i,j,nk+1)=karman*zf(i,j,nk+1)/(1.0+karman*zf(i,j,nk+1)*tem)
    enddo
  ENDIF

!--------------------------------------------------------------
!  Surface:

    IF(bcturbu.eq.1.or.bcturbu.eq.2)THEN

      do i=1,ni
        tem1 = 0.0
      IF(dosfcflx)THEN
        tem2 = g*( thflux(i,j)*rth0s(i,j)+repsm1*qvflux(i,j) )
      ELSE
        tem2 = 0.0
      ENDIF
        kmv(i,j,1) = ( (lvz(i,j,1)**4)*max(tem1+tem2,0.0) )**0.33333333
      enddo

    ELSEIF(bcturbu.eq.3)THEN

      do i=1,ni
        tem1 = ( ( t13(i  ,j,1)*ua(i  ,j,1)     &
                  +t13(i+1,j,1)*ua(i+1,j,1) )   &
               + ( t23(i,j  ,1)*va(i,j  ,1)     &
                  +t23(i,j+1,1)*va(i,j+1,1) )   &
               )*rdz*mf(i,j,1)
      IF(dosfcflx)THEN
        tem2 = g*( thflux(i,j)*rth0s(i,j)+repsm1*qvflux(i,j) )
      ELSE
        tem2 = 0.0
      ENDIF
        kmv(i,j,1) = ( (lvz(i,j,1)**4)*max(tem1+tem2,0.0) )**0.33333333
      enddo

    ENDIF

!------------------------------------------------------------
!  Top/bottom of domain:

      do i=1,ni
        kmh(i,j,  1) = kmh(i,j,2)
        kmh(i,j,nk+1) = 0.0
        kmv(i,j,nk+1) = 0.0
      enddo

    ENDDO

!--------------------------------------------------------------
! boundary conditions:

      if(timestats.ge.1) time_turb=time_turb+mytime()
      call bct(kmh)
      call bct(kmv)
#ifdef MPI
      call comm_1t_start(kmh,khcw1,khcw2,khce1,khce2,   &
                             khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call comm_1t_start(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                             kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call comm_1t_end(kmh,khcw1,khcw2,khce1,khce2,   &
                           khcs1,khcs2,khcn1,khcn2,reqs_khc)
      call bct2(kmh)
      call getcornert(kmh,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call comm_1t_end(kmv,kvcw1,kvcw2,kvce1,kvce2,   &
                           kvcs1,kvcs2,kvcn1,kvcn2,reqs_kvc)
      call bct2(kmv)
      call getcornert(kmv,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
#endif

!--------------------------------------------------------------
!  calculate Kh
!  and also limit horizontal coeffs for numerical stability:
!  calculate dissipation rate

      IF( l_h .gt. 0.001 )THEN
        tem = 1.0/(l_h**4)
      ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO j=0,nj+1

      do k=1,nk+1
      do i=0,ni+1
        khh(i,j,k)=kmh(i,j,k)*prinv
        khh(i,j,k) = min( khh(i,j,k) , temx*ruh(i)*ruh(i) , temy*rvh(j)*rvh(j) )
        khv(i,j,k)=kmv(i,j,k)*prinv
      enddo
      enddo

    IF( j.ge.1 .and. j.le.nj )THEN
    IF( l_h .gt. 0.001 )THEN
      do k=1,nk+1
      do i=1,ni
        dissten(i,j,k) = (kmv(i,j,k)**3)/(lvz(i,j,k)**4)   &
                        +(kmh(i,j,k)**3)*tem
      enddo
      enddo
    ELSE
      do k=1,nk+1
      do i=1,ni
        dissten(i,j,k) = (kmv(i,j,k)**3)/(lvz(i,j,k)**4)
      enddo
      enddo
    ENDIF
    ENDIF

    ENDDO

!--------------------------------------------------------------

!!!      if(nstep.eq.1)then
!!!        do k=1,nk+1
!!!          print *,k,lvz(1,1,k)
!!!        enddo
!!!      endif

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine gettau(dodrag,xf,rxf,rho,rf,kmh,kmv,t11,t12,t13,t22,t23,t33,ua)
      implicit none
      
      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      logical, intent(in) :: dodrag
      real, dimension(ib:ie+1) :: xf,rxf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho,rf
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv
      real, dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, dimension(ib:ie+1,jb:je,kb:ke) :: ua
        
      integer i,j,k

!-----------------------------------------------------------------------
! Note:  turb coefficients are defined on w points

  IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    do k=1,nk

      do j=0,nj+1
      do i=0,ni+1
        t11(i,j,k)=t11(i,j,k)*(kmh(i,j,k)+kmh(i,j,k+1))*rho(i,j,k)
        t22(i,j,k)=t22(i,j,k)*(kmh(i,j,k)+kmh(i,j,k+1))*rho(i,j,k)
        t33(i,j,k)=t33(i,j,k)*(kmv(i,j,k)+kmv(i,j,k+1))*rho(i,j,k)
      enddo
      enddo

      do j=1,nj+1
      do i=1,ni+1
        t12(i,j,k)=t12(i,j,k)*0.03125                                         &
     *( ( (kmh(i-1,j-1,k  )+kmh(i,j,k  ))+(kmh(i-1,j,k  )+kmh(i,j-1,k  )) )   &
       +( (kmh(i-1,j-1,k+1)+kmh(i,j,k+1))+(kmh(i-1,j,k+1)+kmh(i,j-1,k+1)) ) ) &
           *( (rho(i-1,j-1,k)+rho(i,j,k))+(rho(i-1,j,k)+rho(i,j-1,k)) )
      enddo
      enddo

    IF(k.ge.2)THEN
      do j=1,nj+1
      do i=1,ni+1
        t13(i,j,k)=t13(i,j,k)*0.25        &
           *( kmv(i-1,j,k)+kmv(i,j,k) )   &
           *( rf(i-1,j,k)+rf(i,j,k) )
        t23(i,j,k)=t23(i,j,k)*0.25        &
           *( kmv(i,j-1,k)+kmv(i,j,k) )   &
           *( rf(i,j-1,k)+rf(i,j,k) )
      enddo
      enddo
    ENDIF

    enddo

  ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    do k=1,nk

      do j=1,nj
      do i=0,ni+1
        t11(i,j,k)=t11(i,j,k)*(kmh(i,j,k)+kmh(i,j,k+1))*rho(i,j,k)
        t33(i,j,k)=t33(i,j,k)*(kmv(i,j,k)+kmv(i,j,k+1))*rho(i,j,k)
      enddo
      enddo

      do j=1,nj
      do i=2,ni+1
        t22(i,j,k)=2.0*rho(1,1,k)     &
                  *0.25*(kmh(i-1,j,k)+kmh(i,j,k)+kmh(i-1,j,k+1)+kmh(i,j,k+1))  &
                  *ua(i,j,k)*rxf(i)
      enddo
      enddo

      do j=1,nj
      do i=2,ni+1
        t12(i,j,k)=t12(i,j,k)*rho(1,1,k)   &
                  *0.25*(kmh(i,j,k+1)+kmh(i,j,k)+kmh(i-1,j,k+1)+kmh(i-1,j,k))
      enddo
      enddo

    IF(k.ge.2)THEN
      do j=1,nj+1
      do i=1,ni+1
        t13(i,j,k)=t13(i,j,k)*0.25        &
           *( kmv(i-1,j,k)+kmv(i,j,k) )   &
           *( rf(i-1,j,k)+rf(i,j,k) )
        t23(i,j,k)=t23(i,j,k)*0.25        &
           *( kmv(i,j-1,k)+kmv(i,j,k) )   &
           *( rf(i,j-1,k)+rf(i,j,k) )
      enddo
      enddo
    ENDIF

    enddo

  ENDIF

!--------------------------------------------------------------
!  lateral boundary conditions for axisymmetric simulations

    IF(axisymm.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk
        t22(1,1,k)=0.0
      enddo

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk
        t12(1,1,k)=0.0
      enddo

!$omp parallel do default(shared)   &
!$omp private(k)
      do k=1,nk+1
        t13(1,1,k)=0.0
      enddo

    ENDIF

!--------------------------------------------------------------
!  lower boundary conditions

    IF(dodrag)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1
        t13(i,j,1)=t13(i,j,1)*0.5*(rf(i-1,j,1)+rf(i,j,1))
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        t23(i,j,1)=t23(i,j,1)*0.5*(rf(i,j-1,1)+rf(i,j,1))
      enddo
      enddo

    ELSE

      IF(bcturbu.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni+1
          t13(i,j,1)=0.0
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni
          t23(i,j,1)=0.0
        enddo
        enddo

      ELSEIF(bcturbu.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni+1
          t13(i,j,1)=t13(i,j,2)
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni
          t23(i,j,1)=t23(i,j,2)
        enddo
        enddo

      ENDIF

    ENDIF

!--------------------------------------------------------------
!  upper boundary conditions

      IF(bcturbu.eq.1)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni+1
          t13(i,j,nk+1)=0.0
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni
          t23(i,j,nk+1)=0.0
        enddo
        enddo

      ELSEIF(bcturbu.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni+1
          t13(i,j,nk+1)=t13(i,j,nk)
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj+1
        do i=1,ni
          t23(i,j,nk+1)=t23(i,j,nk)
        enddo
        enddo

      ENDIF

!--------------------------------------------------------------
!  finished

      if(timestats.ge.1) time_turb=time_turb+mytime()
 
      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine calcdef(dodrag,rds,sigma,rdsf,sigmaf,zs,gz,rgz,gzu,rgzu,gzv,rgzv,  &
                         xh,rxh,uh,xf,rxf,uf,vh,vf,mh,c1,c2,mf,defsq,defh,       &
                         dum1,dum2,ua,va,wa,t11,t12,t13,t22,t23,t33,gx,gy)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      logical, intent(in) :: dodrag
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(itb:ite,jtb:jte) :: zs,gz,rgz,gzu,rgzu,gzv,rgzv
      real, intent(in), dimension(ib:ie) :: xh,rxh,uh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: defsq,defh,dum1,dum2
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gy
        
      integer :: i,j,k
      real :: r1,r2
      real :: tmp11,tmp22,tmp33,tmp12,tmp13,tmp23

!----------------------------------------------------------------------
!
!  Reference:  Mason, 1989, JAS, p. 1497
!
!----------------------------------------------------------------------

  IF(.not.terrain_flag)THEN

  IF( axisymm.eq.0 )THEN
    ! Cartesian without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk

      do j=1,nj
      do i=0,ni+1 
        t11(i,j,k)=(ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)
      enddo
      enddo
      do j=1,nj+1 
      do i=1,ni+1
        t12(i,j,k)=(ua(i,j,k)-ua(i,j-1,k))*rdy*vf(j)   &
                  +(va(i,j,k)-va(i-1,j,k))*rdx*uf(i)
      enddo
      enddo       
      do j=0,nj+1
      do i=1,ni
        t22(i,j,k)=(va(i,j+1,k)-va(i,j,k))*rdy*vh(j)
      enddo
      enddo
      do j=1,nj
      do i=1,ni
        t33(i,j,k)=(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(i,j,k)
      enddo
      enddo
    IF(k.ge.2)THEN
      do j=1,nj
      do i=1,ni+1
        t13(i,j,k)=(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)   &
                  +(ua(i,j,k)-ua(i,j,k-1))*rdz*0.5*(mf(i-1,j,k)+mf(i,j,k))
      enddo
      enddo
      do j=1,nj+1   
      do i=1,ni
        t23(i,j,k)=(wa(i,j,k)-wa(i,j-1,k))*rdy*vf(j)   &
                  +(va(i,j,k)-va(i,j,k-1))*rdz*0.5*(mf(i,j-1,k)+mf(i,j,k))
      enddo
      enddo
    ENDIF

    ENDDO

  ELSE
    ! axisymmetric:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk

      do j=1,nj
      do i=0,ni+1 
        t11(i,j,k)=(ua(i+1,j,k)-ua(i,j,k))*rdx*uh(i)
      enddo
      enddo
      do j=1,nj
      do i=2,ni+1
        t12(i,j,k)=xf(i)*(va(i,j,k)*rxh(i)-va(i-1,j,k)*rxh(i-1))*rdx*uf(i)
      enddo
      enddo
      do j=1,nj
      do i=2,ni+1
        t22(i,j,k)=ua(i,j,k)*rxf(i)
      enddo
      enddo
      do j=1,nj
      do i=1,ni
        t33(i,j,k)=(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(i,j,k)
      enddo
      enddo
    IF(k.ge.2)THEN
      do j=1,nj
      do i=1,ni+1
        t13(i,j,k)=(wa(i,j,k)-wa(i-1,j,k))*rdx*uf(i)   &
                  +(ua(i,j,k)-ua(i,j,k-1))*rdz*0.5*(mf(i-1,j,k)+mf(i,j,k))
      enddo
      enddo
      do j=1,nj   
      do i=1,ni
        t23(i,j,k)=(va(i,j,k)-va(i,j,k-1))*rdz*mf(1,1,k)
      enddo
      enddo
    ENDIF

    ENDDO

      ! lateral boundary condition:
!$omp parallel do default(shared)   &
!$omp private(k)
      do k=0,nk+1
        t22(1,1,k)=0.0
        t12(1,1,k)=0.0
        t13(1,1,k)=0.0
      enddo

  ENDIF

!------------------------------------------------------------------
!  Cartesian with terrain:

  ELSE

    ! dum1 stores u at w-pts:
    ! dum2 stores v at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
    do j=0,nj+2
      ! lowest model level:
      r1 = (sigma(1)-sigmaf(1))*0.5*rdsf(1)
      do i=0,ni+2
        dum1(i,j,1) = ua(i,j,1)-r1*(-3.0*ua(i,j,1)+4.0*ua(i,j,2)-ua(i,j,3))
        dum2(i,j,1) = va(i,j,1)-r1*(-3.0*va(i,j,1)+4.0*va(i,j,2)-va(i,j,3))
      enddo

      ! upper-most model level:
      r1 = (sigmaf(nk+1)-sigma(nk))*0.5*rdsf(nk)
      do i=0,ni+2
        dum1(i,j,nk+1) = ua(i,j,nk)+r1*(3.0*ua(i,j,nk)-4.0*ua(i,j,nk-1)+ua(i,j,nk-2))
        dum2(i,j,nk+1) = va(i,j,nk)+r1*(3.0*va(i,j,nk)-4.0*va(i,j,nk-1)+va(i,j,nk-2))
      enddo

      ! interior:
      do k=2,nk
      r2 = (sigmaf(k)-sigma(k-1))*rds(k)
      r1 = 1.0-r2
      do i=0,ni+2
        dum1(i,j,k) = r1*ua(i,j,k-1)+r2*ua(i,j,k)
        dum2(i,j,k) = r1*va(i,j,k-1)+r2*va(i,j,k)
      enddo
      enddo
    enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1)
    DO k=1,nk
      do j=1,nj
      do i=0,ni+1 
        t11(i,j,k)=gz(i,j)*(ua(i+1,j,k)*rgzu(i+1,j)-ua(i,j,k)*rgzu(i,j))*rdx*uh(i) &
                  +( gx(i,j,k+1)*(dum1(i,j,k+1)+dum1(i+1,j,k+1))      &
                    -gx(i,j,k  )*(dum1(i,j,k  )+dum1(i+1,j,k  ))      &
                   )*0.5*rdsf(k)
      enddo
      enddo
      do j=0,nj+1
      do i=1,ni
        t22(i,j,k)=gz(i,j)*(va(i,j+1,k)*rgzv(i,j+1)-va(i,j,k)*rgzv(i,j))*rdy*vh(j) &
                  +( gy(i,j,k+1)*(dum2(i,j,k+1)+dum2(i,j+1,k+1))      &
                    -gy(i,j,k  )*(dum2(i,j,k  )+dum2(i,j+1,k  ))      &
                   )*0.5*rdsf(k)
      enddo
      enddo
      do j=1,nj+1 
      do i=1,ni+1
        r1 = zt/(zt-0.25*((zs(i-1,j-1)+zs(i,j))+(zs(i-1,j)+zs(i,j-1))))
        t12(i,j,k)=( r1*(ua(i,j,k)*rgzu(i,j)-ua(i,j-1,k)*rgzu(i,j-1))*rdy*vf(j)  &
                    +0.5*( (zt-sigmaf(k+1))*(dum1(i,j-1,k+1)+dum1(i,j,k+1))      &
                          -(zt-sigmaf(k  ))*(dum1(i,j-1,k  )+dum1(i,j,k  ))      &
                         )*rdsf(k)*r1*(rgzu(i,j)-rgzu(i,j-1))*rdy*vf(j) )        &
                  +( r1*(va(i,j,k)*rgzv(i,j)-va(i-1,j,k)*rgzv(i-1,j))*rdx*uf(i)  &
                    +0.5*( (zt-sigmaf(k+1))*(dum2(i-1,j,k+1)+dum2(i,j,k+1))      &
                          -(zt-sigmaf(k  ))*(dum2(i-1,j,k  )+dum2(i,j,k  ))      &
                         )*rdsf(k)*r1*(rgzv(i,j)-rgzv(i-1,j))*rdx*uf(i) )
      enddo
      enddo       
      do j=1,nj
      do i=1,ni
        t33(i,j,k)=(wa(i,j,k+1)-wa(i,j,k))*rdz*mh(i,j,k)
      enddo
      enddo
    ENDDO

    ! dum1 stores w at scalar-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=0,nj+1
      do i=0,ni+1
        dum1(i,j,k)=0.5*(wa(i,j,k)+wa(i,j,k+1))
      enddo
      enddo
    ENDDO

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=2,nk
      do j=1,nj
      do i=1,ni+1
        t13(i,j,k)=(ua(i,j,k)-ua(i,j,k-1))*rdz*0.5*(mf(i-1,j,k)+mf(i,j,k))        &
                  +(wa(i,j,k)*rgz(i,j)-wa(i-1,j,k)*rgz(i-1,j))*gzu(i,j)*rdx*uf(i) &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum1(i,j,k  )+dum1(i-1,j,k  ))       &
                               -(zt-sigma(k-1))*(dum1(i,j,k-1)+dum1(i-1,j,k-1)) )     &
                             *gzu(i,j)*(rgz(i,j)-rgz(i-1,j))*rdx*uf(i)
      enddo
      enddo
      do j=1,nj+1   
      do i=1,ni
        t23(i,j,k)=(va(i,j,k)-va(i,j,k-1))*rdz*0.5*(mf(i,j-1,k)+mf(i,j,k))        &
                  +(wa(i,j,k)*rgz(i,j)-wa(i,j-1,k)*rgz(i,j-1))*gzv(i,j)*rdy*vf(j) &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum1(i,j,k  )+dum1(i,j-1,k  ))       &
                               -(zt-sigma(k-1))*(dum1(i,j,k-1)+dum1(i,j-1,k-1)) )     &
                             *gzv(i,j)*(rgz(i,j)-rgz(i,j-1))*rdy*vf(j)
      enddo
      enddo
    ENDDO

  ENDIF

!  end of calculations for terrain
!------------------------------------------------------------------
!  lower boundary conditions

  IF(.not.dodrag)THEN

    IF(bcturbu.eq.1)THEN
      
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1 
        t13(i,j,1)=0.0
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        t23(i,j,1)=0.0
      enddo
      enddo

    ELSEIF(bcturbu.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1 
        t13(i,j,1)=t13(i,j,2)
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        t23(i,j,1)=t23(i,j,2)
      enddo
      enddo

    ELSEIF(bcturbu.eq.3)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1
        t13(i,j,1   )= 2.0*ua(i,j,1 )*rdz*0.5*(mf(i-1,j,1   )+mf(i,j,1   ))
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1   
      do i=1,ni
        t23(i,j,1   )= 2.0*va(i,j,1 )*rdz*0.5*(mf(i,j-1,1   )+mf(i,j,1   ))
      enddo
      enddo

    ENDIF

  ENDIF

!------------------------------------------------------------------
!  upper boundary conditions

    IF(bcturbu.eq.1)THEN
      
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1 
        t13(i,j,nk+1)=0.0
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        t23(i,j,nk+1)=0.0
      enddo
      enddo

    ELSEIF(bcturbu.eq.2)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1 
        t13(i,j,nk+1)=t13(i,j,nk)
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1
      do i=1,ni
        t23(i,j,nk+1)=t23(i,j,nk)
      enddo
      enddo

    ELSEIF(bcturbu.eq.3)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni+1
        t13(i,j,nk+1)=-2.0*ua(i,j,nk)*rdz*0.5*(mf(i-1,j,nk+1)+mf(i,j,nk+1))
      enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj+1   
      do i=1,ni
        t23(i,j,nk+1)=-2.0*va(i,j,nk)*rdz*0.5*(mf(i,j-1,nk+1)+mf(i,j,nk+1))
      enddo
      enddo

    ENDIF

!----------------------------------------------------------------------
!  if l_h or l_v is zero, set appropriate terms to zero:
!    (just to be sure)

    IF( iturb.eq.3 .and. l_h.lt.tsmall )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        t11(i,j,k) = 0.0
        t22(i,j,k) = 0.0
        t12(i,j,k) = 0.0
      enddo
      enddo
      enddo
    ENDIF

    IF( iturb.eq.3 .and. l_inf.lt.tsmall )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        t13(i,j,k) = 0.0
        t23(i,j,k) = 0.0
        t33(i,j,k) = 0.0
      enddo
      enddo
      enddo
    ENDIF

!----------------------------------------------------------------------
!  calculate D term:

!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,nj
      do i=1,ni
        defsq(i,j,   1) = 0.0
        defsq(i,j,nk+1) = 0.0
         defh(i,j,   1) = 0.0
         defh(i,j,nk+1) = 0.0
      enddo
      enddo

    IF(axisymm.eq.0)THEN
      ! Cartesian domain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23)
      do k=2,nk
      do j=1,nj
      do i=1,ni

        tmp11=( c1(i,j,k)*t11(i,j,k-1)**2 + c2(i,j,k)*t11(i,j,k)**2 )
        tmp22=( c1(i,j,k)*t22(i,j,k-1)**2 + c2(i,j,k)*t22(i,j,k)**2 )
        tmp33=( c1(i,j,k)*t33(i,j,k-1)**2 + c2(i,j,k)*t33(i,j,k)**2 )

        tmp12=0.25*( c1(i,j,k)*( ( t12(i,j  ,k-1)**2 + t12(i+1,j+1,k-1)**2 )     &
                               + ( t12(i,j+1,k-1)**2 + t12(i+1,j  ,k-1)**2 ) )   &
                    +c2(i,j,k)*( ( t12(i,j  ,k  )**2 + t12(i+1,j+1,k  )**2 )     &
                               + ( t12(i,j+1,k  )**2 + t12(i+1,j  ,k  )**2 ) ) )

        tmp13=0.5*( t13(i,j,k)**2 + t13(i+1,j,k)**2 )

        tmp23=0.5*( t23(i,j,k)**2 + t23(i,j+1,k)**2 )

        defsq(i,j,k)= 2.0*( tmp33 ) + ( tmp13 + tmp23 )

        defh(i,j,k) = 2.0*( tmp11 + tmp22 ) + tmp12

      enddo
      enddo
      enddo

!--------------------------------------------
    ELSE
      ! axisymmetric domain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k,tmp11,tmp22,tmp33,tmp12,tmp13,tmp23)
      do k=2,nk
      do j=1,nj
      do i=1,ni

        tmp11=( c1(i,j,k)*t11(i,j,k-1)**2 + c2(i,j,k)*t11(i,j,k)**2 )
        tmp22=0.5*( c1(i,j,k)*(t22(i,j,k-1)**2 + t22(i+1,j,k-1)**2)  &
                   +c2(i,j,k)*(t22(i,j,k  )**2 + t22(i+1,j,k  )**2) )
        tmp33=( c1(i,j,k)*t33(i,j,k-1)**2 + c2(i,j,k)*t33(i,j,k)**2 )

        tmp12=0.5*(  c1(i,j,k)*( t12(i,j  ,k-1)**2 + t12(i+1,j  ,k-1)**2 )     &
                   + c2(i,j,k)*( t12(i,j  ,k  )**2 + t12(i+1,j  ,k  )**2 ) )

        tmp13=0.5*( t13(i,j,k)**2 + t13(i+1,j,k)**2 )

        tmp23=      t23(i,j,k)**2

        defsq(i,j,k)= 2.0*( tmp33 ) + tmp13 + tmp23

        defh(i,j,k) = 2.0*( tmp11 + tmp22 ) + tmp12

      enddo
      enddo
      enddo

    ENDIF  ! endif for axisymm

!--------------------------------------------------------------
!  finished

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine calcnm(c1,c2,mf,pi0,thv0,th0,cloudvar,nm,t,qt,thv,cloud,   &
                        prs,pp,th,qa)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
      include 'goddard.incl'

      logical, dimension(maxq) :: cloudvar
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,th0
      real, dimension(ib:ie,jb:je,kb:ke) :: nm,t,qt,thv,cloud,prs
      real, dimension(ib:ie,jb:je,kb:ke) :: pp,th
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa

      integer i,j,k,n
      real pavg,tavg,qtavg,qvs,lhv,cpml,gamma,qiavg,qsavg,qgavg,drdt
      real qlavg,qvl,qvi,fliq,fice,nmtmp
      real rslf,rsif

!----------------------------------------------------------------------
!  Dry nm

    IF(imoist.eq.0)then

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
    do j=1,nj

      do k=2,nk
      do i=1,ni
        nm(i,j,k)=alog( (th0(i,j,k)+th(i,j,k))/(th0(i,j,k-1)+th(i,j,k-1)) ) &
                    *g*rdz*mf(i,j,k)
      enddo
      enddo
      do i=1,ni
        nm(i,j,   1)=0.0
        nm(i,j,nk+1)=0.0
      enddo

    enddo

!-----------------------------------------------------------------------
!  Moist nm

    ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,pavg,tavg,qtavg,qvs,lhv,cpml,drdt,gamma,nmtmp)
    DO j=1,nj

      do k=1,nk
      do i=1,ni
        t(i,j,k)=(th0(i,j,k)+th(i,j,k))*(pi0(i,j,k)+pp(i,j,k))
      enddo
      enddo

      do k=1,nk
      do i=1,ni
        qt(i,j,k)=0.0
      enddo
      enddo

      DO n=1,numq
        IF( (n.eq.nqv) .or.                                 &
            (n.ge.nql1.and.n.le.nql2) .or.                  &
            (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )THEN
          do k=1,nk
          do i=1,ni
            qt(i,j,k)=qt(i,j,k)+qa(i,j,k,n)
          enddo
          enddo
        ENDIF
      ENDDO

      do k=1,nk
      do i=1,ni
        thv(i,j,k)=(th0(i,j,k)+th(i,j,k))*(1.0+reps*qa(i,j,k,nqv))   &
                                         /(1.0+qt(i,j,k))
      enddo
      enddo

      do k=2,nk
      do i=1,ni
        nm(i,j,k)=g*alog(thv(i,j,k)/thv(i,j,k-1))*rdz*mf(i,j,k)
      enddo
      enddo

      do i=1,ni
        nm(i,j,   1)=0.0
        nm(i,j,nk+1)=0.0
      enddo

      do k=1,nk
      do i=1,ni
        cloud(i,j,k)=0.0
      enddo
      enddo
      do n=1,numq
        if(cloudvar(n))then
          do k=1,nk
          do i=1,ni
            cloud(i,j,k)=cloud(i,j,k)+qa(i,j,k,n)
          enddo
          enddo
        endif
      enddo

      do k=2,nk
      do i=1,ni
        IF( (cloud(i,j,k).ge.clwsat) .or. (cloud(i,j,k-1).ge.clwsat) )THEN
          pavg = c1(i,j,k)*prs(i,j,k-1)+c2(i,j,k)*prs(i,j,k)
          tavg =   c1(i,j,k)*t(i,j,k-1)+  c2(i,j,k)*t(i,j,k)
          qtavg=  c1(i,j,k)*qt(i,j,k-1)+ c2(i,j,k)*qt(i,j,k)
          qvs=rslf(pavg,tavg)
          lhv=lv1-lv2*tavg
          cpml=cp+cpv*qvs+cpl*(qtavg-qvs)

          drdt=17.67*(273.15-29.65)*qvs/((tavg-29.65)**2)
          gamma=g*(1.0+qtavg)*(1.0+lhv*qvs/(rd*tavg))/(cpml+lhv*drdt)
          nmtmp=g*( ( alog(t(i,j,k)/t(i,j,k-1))*rdz*mf(i,j,k)      &
                            +gamma/tavg )*(1.0+tavg*drdt/(eps+qvs))   &
                         -alog((1.0+qt(i,j,k))/(1.0+qt(i,j,k-1)))*rdz*mf(i,j,k) )
        IF( (cloud(i,j,k).ge.clwsat) .and. (cloud(i,j,k-1).ge.clwsat) )THEN
          nm(i,j,k)=nmtmp
        ELSE
          nm(i,j,k)=0.5*(nm(i,j,k)+nmtmp)
        ENDIF
        ENDIF
      enddo
      enddo

    ENDDO

    ENDIF    ! endif for imoist

!----------------------------------------------------------------------

      if(timestats.ge.1) time_turb=time_turb+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbs(iflux,dt,dosfcflx,xh,rxh,uh,xf,uf,vh,vf,sflux,   &
                       rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                       turbx,turby,turbz,dumx,dumy,dumz,rho,rr,rf,s,sten,khh,khv,doimpl)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      integer iflux
      real :: dt
      logical, intent(in) :: dosfcflx
      real, dimension(ib:ie) :: xh,rxh,uh
      real, dimension(ib:ie+1) :: xf,uf
      real, dimension(jb:je) :: vh
      real, dimension(jb:je+1) :: vf
      real, dimension(ib:ie,jb:je) :: sflux
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dumx,dumy,dumz,rho,rr,rf,s,sten
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: khh,khv
      logical, dimension(ib:ie,jb:je) :: doimpl

      integer :: i,j,k,k1,k2
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2
      real, dimension(nk) :: cfa,cfb,cfc,cfd,s2
      real, dimension(nk) :: lgbth,lgbph

!---------------------------------------------------------------

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN
      ! Cartesian without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        !  x-direction
        do j=1,nj
        do i=1,ni+1
          dumx(i,j,k)= -0.125*( rho(i,j,k)+rho(i-1,j,k) )           &
                             *(  (khh(i,j,k  )+ khh(i-1,j,k  ))     &
                                +(khh(i,j,k+1)+ khh(i-1,j,k+1)) )   &
                             *(    s(i,j,k)-   s(i-1,j,k) )*rdx*uf(i)
        enddo
        enddo
        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=-(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj+1
        do i=1,ni
          dumy(i,j,k)= -0.125*( rho(i,j,k)+rho(i,j-1,k) )           &
                             *(  (khh(i,j,k  )+ khh(i,j-1,k  ))     &
                                +(khh(i,j,k+1)+ khh(i,j-1,k+1)) )   &
                            *(    s(i,j,k)-   s(i,j-1,k) )*rdy*vf(j)
        enddo
        enddo
        do j=1,nj
        do i=1,ni
          turby(i,j,k)=-(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
        enddo
        enddo

      enddo

    ELSE
      ! axisymmetric:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        do j=1,nj
        do i=1,ni+1
          dumx(i,j,k)= -0.125*( rho(i,j,k)+rho(i-1,j,k) )           &
                             *(  (khh(i,j,k  )+ khh(i-1,j,k  ))     &
                                +(khh(i,j,k+1)+ khh(i-1,j,k+1)) )   &
                             *(    s(i,j,k)-   s(i-1,j,k) )*rdx*uf(i)
        enddo
        enddo
        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=-(xf(i+1)*dumx(i+1,j,k)-xf(i)*dumx(i,j,k))*rdx*uh(i)*rxh(i)
          turby(i,j,k)=0.0
        enddo
        enddo

      enddo

    ENDIF   ! endif for axisymm check

!---------------------------------------------------------------

  ELSE
      ! Cartesian with terrain:

      ! use turbz as a temporary array for s at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do j=0,nj+1

        ! lowest model level:
        r1 = (sigma(1)-sigmaf(1))*0.5*rdsf(1)
        do i=0,ni+1
          turbz(i,j,1) = s(i,j,1)-r1*(-3.0*s(i,j,1)+4.0*s(i,j,2)-s(i,j,3))
        enddo

        ! upper-most model level:
        r1 = (sigmaf(nk+1)-sigma(nk))*0.5*rdsf(nk)
        do i=0,ni+1
          turbz(i,j,nk+1) = s(i,j,nk)+r1*(3.0*s(i,j,nk)-4.0*s(i,j,nk-1)+s(i,j,nk-2))
        enddo

        ! interior:
        do k=2,nk
        r2 = (sigmaf(k)-sigma(k-1))*rds(k)
        r1 = 1.0-r2
        do i=0,ni+1
          turbz(i,j,k) = r1*s(i,j,k-1)+r2*s(i,j,k)
        enddo
        enddo

      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        ! x-flux
        do j=1,nj
        do i=1,ni+1
          dumx(i,j,k)= -0.125*( rho(i,j,k)+rho(i-1,j,k) )           &
                             *(  (khh(i,j,k  )+ khh(i-1,j,k  ))     &
                                +(khh(i,j,k+1)+ khh(i-1,j,k+1)) )*( &
                  (s(i,j,k)*rgz(i,j)-s(i-1,j,k)*rgz(i-1,j))         &
                   *gzu(i,j)*rdx*uf(i)                              &
              +0.5*( gxu(i,j,k+1)*(turbz(i,j,k+1)+turbz(i-1,j,k+1)) &
                    -gxu(i,j,k  )*(turbz(i,j,k  )+turbz(i-1,j,k  )) &
                   )*rdsf(k) )
        enddo
        enddo

        ! y-flux
        do j=1,nj+1
        do i=1,ni
          dumy(i,j,k)= -0.125*( rho(i,j,k)+rho(i,j-1,k) )           &
                             *(  (khh(i,j,k  )+ khh(i,j-1,k  ))     &
                                +(khh(i,j,k+1)+ khh(i,j-1,k+1)) )*( &
                  (s(i,j,k)*rgz(i,j)-s(i,j-1,k)*rgz(i,j-1))         &
                   *gzv(i,j)*rdy*vf(j)                              &
              +0.5*( gyv(i,j,k+1)*(turbz(i,j,k+1)+turbz(i,j-1,k+1)) &
                    -gyv(i,j,k  )*(turbz(i,j,k  )+turbz(i,j-1,k  )) &
                   )*rdsf(k) )
        enddo
        enddo

      enddo

      ! use turbz,dumz as temporary arrays for fluxes at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do j=1,nj+1
        ! lowest model level:
        r1 = (sigma(1)-sigmaf(1))*0.5*rdsf(1)
        do i=1,ni+1
          turbz(i,j,1) = dumx(i,j,1)-r1*(-3.0*dumx(i,j,1)+4.0*dumx(i,j,2)-dumx(i,j,3))
           dumz(i,j,1) = dumy(i,j,1)-r1*(-3.0*dumy(i,j,1)+4.0*dumy(i,j,2)-dumy(i,j,3))
        enddo

        ! upper-most model level:
        r1 = (sigmaf(nk+1)-sigma(nk))*0.5*rdsf(nk)
        do i=1,ni+1
          turbz(i,j,nk+1) = dumx(i,j,nk)+r1*(3.0*dumx(i,j,nk)-4.0*dumx(i,j,nk-1)+dumx(i,j,nk-2))
           dumz(i,j,nk+1) = dumy(i,j,nk)+r1*(3.0*dumy(i,j,nk)-4.0*dumy(i,j,nk-1)+dumy(i,j,nk-2))
        enddo

        ! interior:
        do k=2,nk
        r2 = (sigmaf(k)-sigma(k-1))*rds(k)
        r1 = 1.0-r2
        do i=1,ni+1
          turbz(i,j,k) = r1*dumx(i,j,k-1)+r2*dumx(i,j,k)
           dumz(i,j,k) = r1*dumy(i,j,k-1)+r2*dumy(i,j,k)
        enddo
        enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        ! x-tendency
        do j=1,nj
        do i=1,ni
          turbx(i,j,k) = -gz(i,j)*( dumx(i+1,j,k)*rgzu(i+1,j)             &
                                   -dumx(i  ,j,k)*rgzu(i  ,j) )*rdx*uh(i) &
                -( ( gx(i,j,k+1)*(turbz(i,j,k+1)+turbz(i+1,j,k+1))        &
                    -gx(i,j,k  )*(turbz(i,j,k  )+turbz(i+1,j,k  )) )      &
                 )*0.5*rdsf(k)
        enddo
        enddo

        ! y-tendency
        do j=1,nj
        do i=1,ni
          turby(i,j,k) = -gz(i,j)*( dumy(i,j+1,k)*rgzv(i,j+1)             &
                                   -dumy(i,j  ,k)*rgzv(i,j  ) )*rdy*vh(j) &
                -( ( gy(i,j,k+1)*( dumz(i,j,k+1)+ dumz(i,j+1,k+1))        &
                    -gy(i,j,k  )*( dumz(i,j,k  )+ dumz(i,j+1,k  )) )      &
                 )*0.5*rdsf(k)
        enddo
        enddo

      enddo

  ENDIF  ! endif for terrain check

!---------------------------------------------------------------------
!  z-direction

    IF( iturb.eq.3 .and. l_inf.lt.tsmall )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

    ELSE

      rdt = 1.0/dt
      tema = -1.0*dt*vialpha*rdz*rdz
      temb = dt*vibeta*rdz*rdz
      temc = dt*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,lgbth,lgbph,tem,s2)
    DO j=1,nj
    DO i=1,ni

      IF(doimpl(i,j))THEN
        ! implicit calculation:

        do k=1,nk
          r1 = khv(i,j,k  )*mf(i,j,k  )*rf(i,j,k  )*mh(i,j,k)*rr(i,j,k)
          r2 = khv(i,j,k+1)*mf(i,j,k+1)*rf(i,j,k+1)*mh(i,j,k)*rr(i,j,k)
          cfa(k) = tema*r1
          cfc(k) = tema*r2
          cfd(k) = s(i,j,k)   &
                 + temb*(r1*s(i,j,k-1)-(r1+r2)*s(i,j,k)+r2*s(i,j,k+1) )
        enddo
        IF(bcturbs.eq.1)THEN
          cfa( 1) = 0.0
          cfc(nk) = 0.0
        ELSEIF(bcturbs.eq.2)THEN
          cfa( 1) = 0.0
          cfc( 1) = 0.0
          cfa(nk) = 0.0
          cfc(nk) = 0.0
        ENDIF
        do k=1,nk
          cfb(k) = 1.0 - cfa(k) - cfc(k)
        enddo
        if(iflux.eq.1 .and. dosfcflx)then
          cfd(1) = cfd(1) + temc*sflux(i,j)*rf(i,j,1)*mh(i,j,1)*rr(i,j,1)
        endif

        lgbth(1)=-cfc(1)/cfb(1)
        lgbph(1)= cfd(1)/cfb(1)
        do k=2,nk
          tem = 1.0/(cfa(k)*lgbth(k-1)+cfb(k))
          lgbth(k)=-cfc(k)*tem
          lgbph(k)=(cfd(k)-cfa(k)*lgbph(k-1))*tem
        enddo
        s2(nk)=lgbph(nk)
        do k=nk-1,1,-1
          s2(k)=lgbth(k)*s2(k+1)+lgbph(k)
        enddo

        do k=1,nk
          turbz(i,j,k) = rho(i,j,k)*(s2(k)-s(i,j,k))*rdt
        enddo

      ELSE
        ! explicit calculation:

        do k=2,nk
          dumz(i,j,k)=-khv(i,j,k)*(s(i,j,k)-s(i,j,k-1))*rdz*mf(i,j,k)*rf(i,j,k)
        enddo

        IF(bcturbs.eq.1)THEN
 
          dumz(i,j,1)=0.0
          dumz(i,j,nk+1)=0.0

        ELSEIF(bcturbs.eq.2)THEN
 
          dumz(i,j,1)=dumz(i,j,2)
          dumz(i,j,nk+1)=dumz(i,j,nk)

        ENDIF

        if(iflux.eq.1 .and. dosfcflx)then
          dumz(i,j,1)=sflux(i,j)*rf(i,j,1)
        endif

        do k=1,nk
          turbz(i,j,k)=-(dumz(i,j,k+1)-dumz(i,j,k))*rdz*mh(i,j,k)
        enddo

      ENDIF

    ENDDO
    ENDDO

  ENDIF

!---------------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=sten(i,j,k)+(turbx(i,j,k)+turby(i,j,k)+turbz(i,j,k))*rr(i,j,k)
      enddo
      enddo
      enddo

!---------------------------------------------------------------------

      if(timestats.ge.1) time_tmix=time_tmix+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbt(dt,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho,rr,rf,          &
                       rds,sigma,gz,rgz,gzu,rgzu,gzv,rgzv,                &
                       turbx,turby,turbz,dumx,dumy,dumz,t,tten,kmh,kmv,doimpl)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,rxh,uh
      real, intent(in), dimension(ib:ie+1) :: xf,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rr
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: rf
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,rgzu,gzv,rgzv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dumx,dumy,dumz
      real, intent(in), dimension(ibt:iet,jbt:jet,kbt:ket) :: t
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tten
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv
      logical, intent(in), dimension(ib:ie,jb:je) :: doimpl

      integer :: i,j,k
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,rrf
      real, dimension(nk+1) :: cfa,cfb,cfc,cfd,cfe,cff,t2

!---------------------------------------------------------------

    IF(.not.terrain_flag)THEN
      ! Cartesian without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        !  x-direction
        do j=1,nj
        do i=1,ni+1
          dumx(i,j,k)= -0.25*( rf(i,j,k)+rf(i-1,j,k) )   &
                           *( kmh(i,j,k)+kmh(i-1,j,k) )   &
                           *(   t(i,j,k)-  t(i-1,j,k) )*rdx*uf(i)
        enddo
        enddo
        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=-(dumx(i+1,j,k)-dumx(i,j,k))*rdx*uh(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj+1
        do i=1,ni
          dumy(i,j,k)= -0.25*( rf(i,j,k)+rf(i,j-1,k) )   &
                           *( kmh(i,j,k)+kmh(i,j-1,k) )   &
                           *(   t(i,j,k)-  t(i,j-1,k) )*rdy*vf(j)
        enddo
        enddo
        do j=1,nj
        do i=1,ni
          turby(i,j,k)=-(dumy(i,j+1,k)-dumy(i,j,k))*rdy*vh(j)
        enddo
        enddo

      enddo

!---------------------------------------------------------------
!  Cartesian with terrain:

    ELSE

      ! turbz stores t at s-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=0,nj+1
        do i=0,ni+1
          turbz(i,j,k) = 0.5*(t(i,j,k)+t(i,j,k+1))
        enddo
        enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        ! x-flux:
        do j=1,nj
        do i=1,ni+1
          dumx(i,j,k)= -0.25*( rf(i,j,k)+rf(i-1,j,k) )                                   &
                           *( kmh(i,j,k)+kmh(i-1,j,k) )*(                                &
                            (t(i,j,k)*rgz(i,j)-t(i-1,j,k)*rgz(i-1,j))*gzu(i,j)*rdx*uf(i) &
                     +0.5*( (zt-sigma(k  ))*(turbz(i-1,j,k  )+turbz(i,j,k  ))            &
                           -(zt-sigma(k-1))*(turbz(i-1,j,k-1)+turbz(i,j,k-1))            &
                          )*rds(k)*(rgz(i,j)-rgz(i-1,j))*gzu(i,j)*rdx*uf(i)              &
                                                        )
        enddo
        enddo

        ! y-flux:
        do j=1,nj+1
        do i=1,ni
          dumy(i,j,k)= -0.25*( rf(i,j,k)+rf(i,j-1,k) )                                   &
                           *( kmh(i,j,k)+kmh(i,j-1,k) )*(                                &
                            (t(i,j,k)*rgz(i,j)-t(i,j-1,k)*rgz(i,j-1))*gzv(i,j)*rdy*vf(j) &
                     +0.5*( (zt-sigma(k  ))*(turbz(i,j-1,k  )+turbz(i,j,k  ))            &
                           -(zt-sigma(k-1))*(turbz(i,j-1,k-1)+turbz(i,j,k-1))            &
                          )*rds(k)*(rgz(i,j)-rgz(i,j-1))*gzv(i,j)*rdy*vf(j)              &
                                                        )
        enddo
        enddo

      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do j=1,nj+1
        do i=1,ni+1
          dumx(i,j,   1)=0.0
          dumx(i,j,nk+1)=0.0
          dumy(i,j,   1)=0.0
          dumy(i,j,nk+1)=0.0
        enddo
        enddo

      ! turbz stores dumx at s-pts:
      !  dumz stores dumy at s-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          turbz(i,j,k)=0.5*(dumx(i,j,k)+dumx(i,j,k+1))
           dumz(i,j,k)=0.5*(dumy(i,j,k)+dumy(i,j,k+1))
        enddo
        enddo
      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        ! x-tendency:
        do j=1,nj
        do i=1,ni
          turbx(i,j,k) = -(dumx(i+1,j,k)*rgzu(i+1,j)-dumx(i,j,k)*rgzu(i,j))*gz(i,j)*rdx*uh(i) &
                         -0.5*( (zt-sigma(k  ))*(turbz(i,j,k  )+turbz(i+1,j,k  ))             &
                               -(zt-sigma(k-1))*(turbz(i,j,k-1)+turbz(i+1,j,k-1))             &
                              )*rds(k)*(rgzu(i+1,j)-rgzu(i,j))*gz(i,j)*rdx*uh(i)
        enddo
        enddo

        ! y-tendency:
        do j=1,nj
        do i=1,ni
          turby(i,j,k) = -(dumy(i,j+1,k)*rgzv(i,j+1)-dumy(i,j,k)*rgzv(i,j))*gz(i,j)*rdy*vh(j) &
                         -0.5*( (zt-sigma(k  ))*( dumz(i,j,k  )+ dumz(i,j+1,k  ))             &
                               -(zt-sigma(k-1))*( dumz(i,j,k-1)+ dumz(i,j+1,k-1))             &
                              )*rds(k)*(rgzv(i,j+1)-rgzv(i,j))*gz(i,j)*rdy*vh(j)
        enddo
        enddo

      enddo

    ENDIF

!---------------------------------------------------------------------
!  z-direction

      rdt = 1.0/dt
      tema = dt*vialpha*rdz*rdz
      temb = dt*vibeta*rdz*rdz
      temc = dt*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,cfe,cff,tem,t2,rrf)
    do j=1,nj
    do i=1,ni

    IF(doimpl(i,j))THEN

        do k=2,nk
          rrf = 1.0/(0.5*(rho(i,j,k-1)+rho(i,j,k)))
          r1 = 0.5*(kmv(i,j,k-1)+kmv(i,j,k  ))*mh(i,j,k-1)*rho(i,j,k-1)*mf(i,j,k)*rrf
          r2 = 0.5*(kmv(i,j,k  )+kmv(i,j,k+1))*mh(i,j,k  )*rho(i,j,k  )*mf(i,j,k)*rrf
          cfa(k) = tema*r1
          cfc(k) = tema*r2
          cfb(k) = 1.0 + cfa(k) + cfc(k)
          cfd(k) = t(i,j,k)    &
                 + temb*(r2*t(i,j,k+1)-(r1+r2)*t(i,j,k)+r1*t(i,j,k-1))
        enddo

        cfe(1)=0.0
        cff(1)=t(i,j,1)
        do k=2,nk
          tem = 1.0/(cfb(k)-cfc(k)*cfe(k-1))
          cfe(k)=cfa(k)*tem
          cff(k)=(cfd(k)+cfc(k)*cff(k-1))*tem
        enddo
        t2(nk+1)=0.0
        do k=nk,2,-1
          t2(k)=cfe(k)*t2(k+1)+cff(k)
        enddo
        do k=2,nk
          turbz(i,j,k) = rf(i,j,k)*(t2(k)-t(i,j,k))*rdt
        enddo

    ELSE

      do k=1,nk
        dumz(i,j,k)=-0.5*(kmv(i,j,k)+kmv(i,j,k+1))*(t(i,j,k+1)-t(i,j,k))*rdz*mh(i,j,k)*rho(i,j,k)
      enddo
      do k=2,nk
        turbz(i,j,k)=-(dumz(i,j,k)-dumz(i,j,k-1))*rdz*mf(i,j,k)
      enddo

    ENDIF

    enddo
    enddo

!---------------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k,rrf)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        rrf = 1.0/(0.5*(rho(i,j,k-1)+rho(i,j,k)))
        tten(i,j,k)=tten(i,j,k)+2.0*((turbx(i,j,k)+turby(i,j,k))+turbz(i,j,k))*rrf
      enddo
      enddo
      enddo

!---------------------------------------------------------------------

      if(timestats.ge.1) time_tmix=time_tmix+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbu(dt,dodrag,xh,ruh,xf,rxf,uf,vh,mh,mf,rmf,rho,rf,  &
                       zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gxu,     &
                       turbx,turby,turbz,dum1,dum2,u,uten,w,t11,t12,t13,t22,kmv,doimpl)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      logical, intent(in) :: dodrag
      real, intent(in), dimension(ib:ie) :: xh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,uf
      real, intent(in), dimension(jb:je) :: vh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf,rmf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rf
      real, intent(in), dimension(itb:ite,jtb:jte) :: zs,gz,rgz,gzu,gzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gxu
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum1,dum2
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: uten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv
      logical, intent(in), dimension(ib:ie,jb:je) :: doimpl

      integer :: i,j,k,ip
      real :: rdt,tema,temb,temc,temd,teme
      real :: tem,r1,r2,rru0
      real, dimension(nk) :: cfa,cfb,cfc,cfd,u2
      real, dimension(nk) :: lgbth,lgbph

!---------------------------------------------------------------

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        !  x-direction
        do j=1,nj
        do i=1,ni+1
          turbx(i,j,k)=(t11(i,j,k)-t11(i-1,j,k))*rdx*uf(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj
        do i=1,ni+1
          turby(i,j,k)=(t12(i,j+1,k)-t12(i,j,k))*rdy*vh(j)
        enddo
        enddo

      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(j,k)
      do k=1,nk

        do j=1,nj
        turbx(1,j,k)=0.0
        turby(1,j,k)=0.0
        do i=2,ni+1
          turbx(i,j,k)=( (xh(i)*t11(i,j,k)-xh(i-1)*t11(i-1,j,k))*rdx*uf(i)  &
                        -t22(i,j,k) )*rxf(i)
          turby(i,j,k)=0.0
        enddo
        enddo

      enddo

    ENDIF

!---------------------------------------------------------------
!  Terrain:

  ELSE

      ! dum1 stores t11 at w-pts:
      ! dum2 stores t12 at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do j=1,nj+1

          ! lowest model level:
          r1 = (sigma(1)-sigmaf(1))*0.5*rdsf(1)
          do i=0,ni+1
            dum1(i,j,1) = t11(i,j,1)-r1*(-3.0*t11(i,j,1)+4.0*t11(i,j,2)-t11(i,j,3))
            dum2(i,j,1) = t12(i,j,1)-r1*(-3.0*t12(i,j,1)+4.0*t12(i,j,2)-t12(i,j,3))
          enddo

          ! upper-most model level:
          r1 = (sigmaf(nk+1)-sigma(nk))*0.5*rdsf(nk)
          do i=0,ni+1
            dum1(i,j,nk+1) = t11(i,j,nk)+r1*(3.0*t11(i,j,nk)-4.0*t11(i,j,nk-1)+t11(i,j,nk-2))
            dum2(i,j,nk+1) = t12(i,j,nk)+r1*(3.0*t12(i,j,nk)-4.0*t12(i,j,nk-1)+t12(i,j,nk-2))
          enddo

          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=0,ni+1
            dum1(i,j,k) = r1*t11(i,j,k-1)+r2*t11(i,j,k)
            dum2(i,j,k) = r1*t12(i,j,k-1)+r2*t12(i,j,k)
          enddo
          enddo

      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do k=1,nk

        !  x-direction
        do j=1,nj
        do i=1,ni+1
          turbx(i,j,k)=gzu(i,j)*(t11(i,j,k)*rgz(i,j)-t11(i-1,j,k)*rgz(i-1,j))*rdx*uf(i)  &
                      +0.5*( gxu(i,j,k+1)*(dum1(i-1,j,k+1)+dum1(i,j,k+1))                &
                            -gxu(i,j,k  )*(dum1(i-1,j,k  )+dum1(i,j,k  )) )*rdsf(k)
        enddo
        enddo

        !  y-direction
        do j=1,nj
        do i=1,ni+1
          r1 = (zt-0.25*((zs(i-1,j-1)+zs(i,j))+(zs(i-1,j)+zs(i,j-1))))*rzt
          r2 = (zt-0.25*((zs(i-1,j+1)+zs(i,j))+(zs(i-1,j)+zs(i,j+1))))*rzt
          turby(i,j,k)=gzu(i,j)*(t12(i,j+1,k)*r2-t12(i,j,k)*r1)*rdy*vh(j)      &
                      +0.5*( (zt-sigmaf(k+1))*(dum2(i,j,k+1)+dum2(i,j+1,k+1))  &
                            -(zt-sigmaf(k  ))*(dum2(i,j,k  )+dum2(i,j+1,k  ))  &
                           )*gzu(i,j)*(r2-r1)*rdy*vh(j)
        enddo
        enddo

      enddo

  ENDIF

!-----------------------------------------------------------------
!  z-direction

    IF( iturb.eq.3 .and. l_inf.lt.tsmall )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

    ELSE

      rdt = 0.5/dt
      tema = -0.0625*dt*vialpha*rdz*rdz
      temb =  0.0625*dt*vibeta*rdz*rdz
      temd =  0.25*dt*rdz
      temc =  0.5*dt*rdz
      teme = 0.25*dx*dx/dt

      ip = 0
      if( axisymm.eq.1 ) ip = 1

  IF(terrain_flag)THEN
    ! dum1 stores w at scalar-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=0,nj+1
      do i=0,ni+1
        dum1(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))
      enddo
      enddo
    ENDDO
  ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,lgbth,lgbph,tem,u2,rru0)
    DO j=1,nj
    DO i=1+ip,ni+1

    IF(doimpl(i-1,j).or.doimpl(i,j))THEN

        do k=1,nk
          rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
          tem = (mh(i-1,j,k)+mh(i,j,k))*rru0
          r1 = (kmv(i-1,j,k  )+kmv(i,j,k  ))*(mf(i-1,j,k  )+mf(i,j,k  ))   &
              *(rf(i-1,j,k  )+rf(i,j,k  ))*tem
          r2 = (kmv(i-1,j,k+1)+kmv(i,j,k+1))*(mf(i-1,j,k+1)+mf(i,j,k+1))   &
              *(rf(i-1,j,k+1)+rf(i,j,k+1))*tem
          cfa(k) = tema*r1
          cfc(k) = tema*r2
          cfd(k) = u(i,j,k)    &
                 + temb*( r2*u(i,j,k+1)-(r1+r2)*u(i,j,k)+r1*u(i,j,k-1) )
        enddo
        IF(bcturbu.eq.1)THEN
          cfa( 1) = 0.0
          cfc(nk) = 0.0
        ELSEIF(bcturbu.eq.2)THEN
          cfa( 1) = 0.0
          cfc( 1) = 0.0
          cfa(nk) = 0.0
          cfc(nk) = 0.0
        ELSEIF(bcturbu.eq.3)THEN
          cfa( 1) = 0.0
          cfc(nk) = 0.0
        ENDIF
        do k=1,nk
          cfb(k) = 1.0 - cfa(k) - cfc(k)
        enddo
        if(dodrag)then
          rru0 = 1.0/(0.5*(rho(i-1,j,1)+rho(i,j,1)))
          tem = temc*t13(i,j,1)*(mh(i-1,j,1)+mh(i,j,1))*rru0
          cfd(1) = cfd(1) - tem
        endif

        lgbth(1)=-cfc(1)/cfb(1)
        lgbph(1)= cfd(1)/cfb(1)
        do k=2,nk
          tem = 1.0/(cfa(k)*lgbth(k-1)+cfb(k))
          lgbth(k)=-cfc(k)*tem
          lgbph(k)=(cfd(k)-cfa(k)*lgbph(k-1))*tem
        enddo
        u2(nk)=lgbph(nk)
        do k=nk-1,1,-1
          u2(k)=lgbth(k)*u2(k+1)+lgbph(k)
        enddo

        do k=1,nk
          turbz(i,j,k) = (rho(i-1,j,k)+rho(i,j,k))*(u2(k)-u(i,j,k))*rdt
        enddo

      ! explicit piece ... dwdx term
      IF(.not.terrain_flag)THEN
        dum2(i,j,1) = 0.0
        do k=2,nk
          dum2(i,j,k)=(w(i,j,k)-w(i-1,j,k))*rdx*uf(i)  &
                     *0.25*( kmv(i-1,j,k)+kmv(i,j,k) )   &
                          *( rf(i-1,j,k)+rf(i,j,k) )
        enddo
        dum2(i,j,nk+1) = 0.0
      ELSE
        dum2(i,j,1) = 0.0
        do k=2,nk
          dum2(i,j,k)=(w(i,j,k)*rgz(i,j)-w(i-1,j,k)*rgz(i-1,j))*gzu(i,j)*rdx*uf(i) &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum1(i,j,k  )+dum1(i-1,j,k  ))    &
                               -(zt-sigma(k-1))*(dum1(i,j,k-1)+dum1(i-1,j,k-1)) )  &
                             *gzu(i,j)*(rgz(i,j)-rgz(i-1,j))*rdx*uf(i)             &
                     *0.25*( kmv(i-1,j,k)+kmv(i,j,k) )                             &
                          *( rf(i-1,j,k)+rf(i,j,k) )
        enddo
        dum2(i,j,nk+1) = 0.0
      ENDIF
      do k=1,nk
        turbz(i,j,k)=turbz(i,j,k)+(dum2(i,j,k+1)-dum2(i,j,k))*rdz*0.5*(mh(i-1,j,k)+mh(i,j,k))
      enddo

    ELSE

      do k=1,nk
        turbz(i,j,k)=(t13(i,j,k+1)-t13(i,j,k))*rdz*0.5*(mh(i-1,j,k)+mh(i,j,k))
      enddo

    ENDIF

    ENDDO
    ENDDO

      IF(axisymm.eq.1)THEN
        DO k=1,nk
          turbz(1,1,k) = 0.0
        ENDDO
      ENDIF

  ENDIF

!-----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k,rru0)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        rru0 = 1.0/(0.5*(rho(i-1,j,k)+rho(i,j,k)))
        uten(i,j,k)=uten(i,j,k)+(turbx(i,j,k)+turby(i,j,k)+turbz(i,j,k))*rru0
      enddo
      enddo
      enddo

!-------------------------------------------------------------------
!  All done

      if(timestats.ge.1) time_tmix=time_tmix+mytime()

      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


      subroutine turbv(dt,dodrag,xh,rxh,uh,xf,rvh,vf,mh,mf,rho,rf,   &
                       zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gyv,  &
                       turbx,turby,turbz,dum1,dum2,v,vten,w,t12,t22,t23,kmv,doimpl)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real, intent(in) :: dt
      logical, intent(in) :: dodrag
      real, intent(in), dimension(ib:ie) :: xh,rxh,uh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: rvh
      real, intent(in), dimension(jb:je+1) :: vf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,rf
      real, intent(in), dimension(itb:ite,jtb:jte) :: zs,gz,rgz,gzu,gzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gyv
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum1,dum2
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: vten
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: w
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: t12,t22,t23
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv
      logical, intent(in), dimension(ib:ie,jb:je) :: doimpl
 
      integer :: i,j,k,ip
      real :: rdt,tema,temb,temc,temd,teme
      real :: tem,r1,r2,rrv0
      real, dimension(nk) :: cfa,cfb,cfc,cfd,v2
      real, dimension(nk) :: lgbth,lgbph

!---------------------------------------------------------------

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        !  x-direction
        do j=1,nj+1
        do i=1,ni
          turbx(i,j,k)=(t12(i+1,j,k)-t12(i,j,k))*rdx*uh(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj+1
        do i=1,ni
          turby(i,j,k)=(t22(i,j,k)-t22(i,j-1,k))*rdy*vf(j)
        enddo
        enddo

      enddo

    ELSE

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk

        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=(xf(i+1)*xf(i+1)*t12(i+1,j,k)-xf(i)*xf(i)*t12(i,j,k))*rdx*uh(i)*rxh(i)*rxh(i)
          turby(i,j,k)=0.0
        enddo
        enddo

      enddo

    ENDIF

!---------------------------------------------------------------
!  Terrain:

  ELSE

      ! dum1 stores t12 at w-pts:
      ! dum2 stores t22 at w-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do j=0,nj+1

          ! lowest model level:
          r1 = (sigma(1)-sigmaf(1))*0.5*rdsf(1)
          do i=1,ni+1
            dum1(i,j,1) = t12(i,j,1)-r1*(-3.0*t12(i,j,1)+4.0*t12(i,j,2)-t12(i,j,3))
            dum2(i,j,1) = t22(i,j,1)-r1*(-3.0*t22(i,j,1)+4.0*t22(i,j,2)-t22(i,j,3))
          enddo

          ! upper-most model level:
          r1 = (sigmaf(nk+1)-sigma(nk))*0.5*rdsf(nk)
          do i=1,ni+1
            dum1(i,j,nk+1) = t12(i,j,nk)+r1*(3.0*t12(i,j,nk)-4.0*t12(i,j,nk-1)+t12(i,j,nk-2))
            dum2(i,j,nk+1) = t22(i,j,nk)+r1*(3.0*t22(i,j,nk)-4.0*t22(i,j,nk-1)+t22(i,j,nk-2))
          enddo

          ! interior:
          do k=2,nk
          r2 = (sigmaf(k)-sigma(k-1))*rds(k)
          r1 = 1.0-r2
          do i=1,ni+1
            dum1(i,j,k) = r1*t12(i,j,k-1)+r2*t12(i,j,k)
            dum2(i,j,k) = r1*t22(i,j,k-1)+r2*t22(i,j,k)
          enddo
          enddo

      enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2)
      do k=1,nk

        !  x-direction
        do j=1,nj+1
        do i=1,ni
          r1 = (zt-0.25*((zs(i-1,j-1)+zs(i,j))+(zs(i-1,j)+zs(i,j-1))))*rzt
          r2 = (zt-0.25*((zs(i+1,j-1)+zs(i,j))+(zs(i+1,j)+zs(i,j-1))))*rzt
          turbx(i,j,k)=gzv(i,j)*(t12(i+1,j,k)*r2-t12(i,j,k)*r1)*rdx*uh(i)      &
                      +0.5*( (zt-sigmaf(k+1))*(dum1(i,j,k+1)+dum1(i+1,j,k+1))  &
                            -(zt-sigmaf(k  ))*(dum1(i,j,k  )+dum1(i+1,j,k  ))  &
                           )*gzv(i,j)*(r2-r1)*rdx*uh(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj+1
        do i=1,ni
          turby(i,j,k)=gzv(i,j)*(t22(i,j,k)*rgz(i,j)-t22(i,j-1,k)*rgz(i,j-1))*rdy*vf(j)  &
                      +0.5*( gyv(i,j,k+1)*(dum2(i,j-1,k+1)+dum2(i,j,k+1))                &
                            -gyv(i,j,k  )*(dum2(i,j-1,k  )+dum2(i,j,k  )) )*rdsf(k)
        enddo
        enddo

      enddo

  ENDIF

!-----------------------------------------------------------------
!  z-direction

    IF( iturb.eq.3 .and. l_inf.lt.tsmall )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

    ELSE

      rdt = 0.5/dt
      tema = -0.0625*dt*vialpha*rdz*rdz
      temb =  0.0625*dt*vibeta*rdz*rdz
      temd =  0.25*dt*rdz
      temc =  0.5*dt*rdz
      teme = 0.25*dy*dy/dt

      ip = 1
      if( axisymm.eq.1 ) ip = 0

  IF(terrain_flag)THEN
    ! dum1 stores w at scalar-pts:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    DO k=1,nk
      do j=0,nj+1
      do i=0,ni+1
        dum1(i,j,k)=0.5*(w(i,j,k)+w(i,j,k+1))
      enddo
      enddo
    ENDDO
  ENDIF

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,lgbth,lgbph,tem,v2,rrv0)
    do j=1,nj+ip
    do i=1,ni

    IF(doimpl(i,j-1).or.doimpl(i,j))THEN

        do k=1,nk
          rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
          tem = (mh(i,j-1,k)+mh(i,j,k))*rrv0
          r1 = (kmv(i,j-1,k  )+kmv(i,j,k  ))*(mf(i,j-1,k  )+mf(i,j,k  ))   &
              *(rf(i,j-1,k  )+rf(i,j,k  ))*tem
          r2 = (kmv(i,j-1,k+1)+kmv(i,j,k+1))*(mf(i,j-1,k+1)+mf(i,j,k+1))   &
              *(rf(i,j-1,k+1)+rf(i,j,k+1))*tem
          cfa(k) = tema*r1
          cfc(k) = tema*r2
          cfd(k) = v(i,j,k)    &
                 + temb*( r2*v(i,j,k+1)-(r1+r2)*v(i,j,k)+r1*v(i,j,k-1) )
        enddo
        IF(bcturbu.eq.1)THEN
          cfa( 1) = 0.0
          cfc(nk) = 0.0
        ELSEIF(bcturbu.eq.2)THEN
          cfa( 1) = 0.0
          cfc( 1) = 0.0
          cfa(nk) = 0.0
          cfc(nk) = 0.0
        ELSEIF(bcturbu.eq.3)THEN
          cfa( 1) = 0.0
          cfc(nk) = 0.0
        ENDIF
        do k=1,nk
          cfb(k) = 1.0 - cfa(k) - cfc(k)
        enddo
        if(dodrag)then
          rrv0 = 1.0/(0.5*(rho(i,j-1,1)+rho(i,j,1)))
          tem = temc*t23(i,j,1)*(mh(i,j-1,1)+mh(i,j,1))*rrv0
          cfd(1) = cfd(1) - tem
        endif

        lgbth(1)=-cfc(1)/cfb(1)
        lgbph(1)= cfd(1)/cfb(1)
        do k=2,nk
          tem = 1.0/(cfa(k)*lgbth(k-1)+cfb(k))
          lgbth(k)=-cfc(k)*tem
          lgbph(k)=(cfd(k)-cfa(k)*lgbph(k-1))*tem
        enddo
        v2(nk)=lgbph(nk)
        do k=nk-1,1,-1
          v2(k)=lgbth(k)*v2(k+1)+lgbph(k)
        enddo

        do k=1,nk
          turbz(i,j,k) = (rho(i,j-1,k)+rho(i,j,k))*(v2(k)-v(i,j,k))*rdt
        enddo

      ! explicit piece ... dwdy term
      ! Cartesian grid only:
      IF(axisymm.eq.0)THEN
      IF(.not.terrain_flag)THEN
        dum2(i,j,1) = 0.0
        do k=2,nk
          dum2(i,j,k)=(w(i,j,k)-w(i,j-1,k))*rdy*vf(j)  &
                     *0.25*( kmv(i,j-1,k)+kmv(i,j,k) )   &
                          *( rf(i,j-1,k)+rf(i,j,k) )
        enddo
        dum2(i,j,nk+1) = 0.0
      ELSE
        dum2(i,j,1) = 0.0
        do k=2,nk
          dum2(i,j,k)=(w(i,j,k)*rgz(i,j)-w(i,j-1,k)*rgz(i,j-1))*gzv(i,j)*rdy*vf(j) &
                  +0.5*rds(k)*( (zt-sigma(k  ))*(dum1(i,j,k  )+dum1(i,j-1,k  ))    &
                               -(zt-sigma(k-1))*(dum1(i,j,k-1)+dum1(i,j-1,k-1)) )  &
                             *gzv(i,j)*(rgz(i,j)-rgz(i,j-1))*rdy*vf(j)             &
                     *0.25*( kmv(i,j-1,k)+kmv(i,j,k) )   &
                          *( rf(i,j-1,k)+rf(i,j,k) )
        enddo
        dum2(i,j,nk+1) = 0.0
      ENDIF
      do k=1,nk
        turbz(i,j,k)=turbz(i,j,k)+(dum2(i,j,k+1)-dum2(i,j,k))*rdz*0.5*(mh(i,j-1,k)+mh(i,j,k))
      enddo
      ENDIF

    ELSE

      do k=1,nk
        turbz(i,j,k)=(t23(i,j,k+1)-t23(i,j,k))*rdz*0.5*(mh(i,j-1,k)+mh(i,j,k))
      enddo

    ENDIF

    enddo
    enddo

  ENDIF

!-----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k,rrv0)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        rrv0 = 1.0/(0.5*(rho(i,j-1,k)+rho(i,j,k)))
        vten(i,j,k)=vten(i,j,k)+(turbx(i,j,k)+turby(i,j,k)+turbz(i,j,k))*rrv0
      enddo
      enddo
      enddo

!-------------------------------------------------------------------
!  All done
 
      if(timestats.ge.1) time_tmix=time_tmix+mytime()
 
      return
      end


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

 
      subroutine turbw(dt,xh,rxh,uh,xf,vh,mh,mf,rho,rf,gz,rgzu,rgzv,rds,sigma,   &
                       turbx,turby,turbz,dum1,dum2,w,wten,t13,t23,t33,kmv,doimpl)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'
 
      real :: dt
      real, dimension(ib:ie) :: xh,rxh,uh
      real, dimension(ib:ie+1) :: xf
      real, dimension(jb:je) :: vh
      real, dimension(ib:ie,jb:je,kb:ke) :: mh
      real, dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, dimension(ib:ie,jb:je,kb:ke) :: rho,rf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgzu,rgzv
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, dimension(ib:ie,jb:je,kb:ke) :: turbx,turby,turbz,dum1,dum2
      real, dimension(ib:ie,jb:je,kb:ke+1) :: w,wten
      real, dimension(ib:ie,jb:je,kb:ke) :: t13,t23,t33
      real, dimension(ibc:iec,jbc:jec,kbc:kec) :: kmv
      logical, dimension(ib:ie,jb:je) :: doimpl
 
      integer :: i,j,k
      real :: rdt,tema,temb,temc
      real :: tem,r1,r2,rrf
      real, dimension(nk+1) :: cfa,cfb,cfc,cfd,cfe,cff,w2

!----------------------------------------------------------------

  IF(.not.terrain_flag)THEN

    IF(axisymm.eq.0)THEN
      ! Cartesian without terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        !  x-direction
        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=(t13(i+1,j,k)-t13(i,j,k))*rdx*uh(i)
        enddo
        enddo

        !  y-direction
        do j=1,nj
        do i=1,ni
          turby(i,j,k)=(t23(i,j+1,k)-t23(i,j,k))*rdy*vh(j)
        enddo
        enddo

      enddo

    ELSE
      ! axisymmetric:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=(xf(i+1)*t13(i+1,j,k)-xf(i)*t13(i,j,k))*rdx*uh(i)*rxh(i)
          turby(i,j,k)=0.0
        enddo
        enddo

      enddo

    ENDIF

!----------------------------------------------------------------

  ELSE
      ! Cartesian with terrain:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k) = 0.25*( (t13(i,j,k+1)+t13(i+1,j,k+1)) &
                              +(t13(i,j,k  )+t13(i+1,j,k  )) )
        enddo
        enddo
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = 0.25*( (t23(i,j,k+1)+t23(i,j+1,k+1)) &
                              +(t23(i,j,k  )+t23(i,j+1,k  )) )
        enddo
        enddo
      enddo


!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk

        do j=1,nj
        do i=1,ni
          turbx(i,j,k)=gz(i,j)*( t13(i+1,j,k)*rgzu(i+1,j)             &
                                -t13(i  ,j,k)*rgzu(i  ,j) )*rdx*uh(i) &
              +( (zt-sigma(k  ))*dum1(i,j,k  )                        &
                -(zt-sigma(k-1))*dum1(i,j,k-1)                        &
               )*gz(i,j)*(rgzu(i+1,j)-rgzu(i,j))*rdx*uh(i)*rds(k)
        enddo
        enddo

        do j=1,nj
        do i=1,ni
          turby(i,j,k)=gz(i,j)*( t23(i,j+1,k)*rgzv(i,j+1)             &
                                -t23(i,j  ,k)*rgzv(i,j  ) )*rdy*vh(j) &
              +( (zt-sigma(k  ))*dum2(i,j,k  )                        &
                -(zt-sigma(k-1))*dum2(i,j,k-1)                        &
               )*gz(i,j)*(rgzv(i,j+1)-rgzv(i,j))*rdy*vh(j)*rds(k)
        enddo
        enddo

      enddo

  ENDIF

!-----------------------------------------------------------------
!  z-direction

    IF( iturb.eq.3 .and. l_inf.lt.tsmall )THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        turbz(i,j,k)=0.0
      enddo
      enddo
      enddo

    ELSE

      rdt = 1.0/dt
      tema = dt*vialpha*rdz*rdz
      temb = dt*vibeta*rdz*rdz
      temc = dt*rdz

!$omp parallel do default(shared)   &
!$omp private(i,j,k,r1,r2,cfa,cfb,cfc,cfd,cfe,cff,tem,w2,rrf)
    do j=1,nj
    do i=1,ni

    IF(doimpl(i,j))THEN

        do k=2,nk
          rrf = 1.0/(0.5*(rho(i,j,k-1)+rho(i,j,k)))
          r1 = (kmv(i,j,k-1)+kmv(i,j,k  ))*mh(i,j,k-1)*rho(i,j,k-1)*mf(i,j,k)*rrf
          r2 = (kmv(i,j,k  )+kmv(i,j,k+1))*mh(i,j,k  )*rho(i,j,k  )*mf(i,j,k)*rrf
          cfa(k) = tema*r1
          cfc(k) = tema*r2
          cfb(k) = 1.0 + cfa(k) + cfc(k)
          cfd(k) = w(i,j,k)    &
                 + temb*(r2*w(i,j,k+1)-(r1+r2)*w(i,j,k)+r1*w(i,j,k-1))
        enddo

        cfe(1)=0.0
        cff(1)=w(i,j,1)
        do k=2,nk
          tem = 1.0/(cfb(k)-cfc(k)*cfe(k-1))
          cfe(k)=cfa(k)*tem
          cff(k)=(cfd(k)+cfc(k)*cff(k-1))*tem
        enddo
        w2(nk+1)=0.0
        do k=nk,2,-1
          w2(k)=cfe(k)*w2(k+1)+cff(k)
        enddo
        do k=2,nk
          turbz(i,j,k) = rf(i,j,k)*(w2(k)-w(i,j,k))*rdt
        enddo

    ELSE

      do k=2,nk
        turbz(i,j,k)=(t33(i,j,k)-t33(i,j,k-1))*rdz*mf(i,j,k)
      enddo

    ENDIF

    enddo
    enddo

  ENDIF

!-----------------------------------------------------------------

!$omp parallel do default(shared)   &
!$omp private(i,j,k,rrf)
      do k=2,nk
      do j=1,nj
      do i=1,ni
        rrf = 1.0/(0.5*(rho(i,j,k-1)+rho(i,j,k)))
        wten(i,j,k)=wten(i,j,k)+(turbx(i,j,k)+turby(i,j,k)+turbz(i,j,k))*rrf
      enddo
      enddo
      enddo

!-------------------------------------------------------------------
!  All done

      if(timestats.ge.1) time_tmix=time_tmix+mytime()
 
      return
      end


