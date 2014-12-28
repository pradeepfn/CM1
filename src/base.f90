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



      subroutine base(th00s,thlr,zh,mh,c1,c2,zf,mf,rho0s,pi0s,prs0s,rth0s,  &
                      pi0,prs0,rho0,thv0,th0,th00,rth00,pi00,qv0,u0,v0,     &
                      qc0,ql0,rr0,rf0,rrf0,t0,rh0,                          &
                      reqs_u,reqs_v,reqs_s,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2, &
                      uw31,uw32,ue31,ue32,us31,us32,un31,un32,              &
                      vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,              &
                      sw31,sw32,se31,se32,ss31,ss32,sn31,sn32)
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'goddard.incl'

      real, intent(inout) :: th00s,thlr
      real, dimension(ib:ie,jb:je,kb:ke) :: zh,mh,c1,c2
      real, dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, dimension(ib:ie,jb:je) :: rho0s,pi0s,prs0s,rth0s
      real, dimension(ib:ie,jb:je,kb:ke) :: pi0,prs0,rho0,thv0,th0,th00,rth00,pi00,qv0
      real, dimension(ib:ie,jb:je,kb:ke) :: qc0,ql0,rr0,rf0,rrf0
      real, dimension(ib:ie,jb:je,kb:ke) :: t0,rh0
      real, dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, dimension(ib:ie,jb:je+1,kb:ke) :: v0
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_s
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(cmp,jmp,kmp) :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, intent(inout), dimension(imp,cmp,kmp) :: vs31,vs32,vn31,vn32
      real, intent(inout), dimension(cmp,jmp,kmp) :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp) :: ss31,ss32,sn31,sn32

!-----------------------------------------------------------------------

      integer i,j,k,n,nn,irec,niter,nsnd,kbot,ktop,tflag,nmax
      real zu,zv
      real z_trop,th_trop,th_sfc,t_trop,prs_sfc,qv_pbl,pi_sfc,t_sfc,rh_sfc,rh_pbl
      real qv_sfc,thv_sfc,psurf,tsurf,qsurf,thsurf,thvsurf
      real qv1,th1,tlast,th2,qv2,ql2,thbar,qvbar,tlcl
      real ns,ns1,ns2,ns3,zl1,zl2,zsfc
      real qcval,qtval,p_sfc,ql_sfc,thn,thlast,qtm,qtp,tavg,lhv,     &
           qvavg,qlavg,qtavg,desdt,gamma,tp,delz,pim,thvm,qvm,qlm,   &
           tm,pm,thm,drdt,pavg,qvl,qvi,fliq,fice,cpml,qim
      real thex,es,qvs,aaa
      real udep,uconst1,uconst2
      real udep1,udep2,umax1,umax2,vmax1,angle
      real tmp,ql,t1,t2,pitmp
      real hs,lapse
      real alpha,umax,nm,dudz,dvdz,rinum
      real, dimension(:), allocatable :: zsnd,thsnd,qvsnd,usnd,vsnd,   &
                                         thvsnd,pisnd,psnd,tsnd,rhsnd
      real, dimension(:), allocatable :: thinterp,qvinterp,uinterp,vinterp, &
                                         pinterp,tinterp,rhinterp
      integer :: kk,kup,kdn
      real :: interp_frac
      real*8 dblepi
      real rslf,rsif

      integer :: flag,ttype
      real :: pisfc,the_sfc,thv1,thv2,pi1,pi2,z1,z2,p2,theq,qt_sfc
      real :: getthe

      real :: ztrop,zmix,qv_mix,zmin,dtheta,thv_trop,qv_trop,pi_trop,p_trop,rhexp
      real, dimension(nk) :: thep,the0
      real :: depth_layer_01,depth_layer_02,shear_layer_01,shear_layer_02
      real :: du,dv

      real, dimension(:), allocatable :: pfoo,tfoo,qvfoo
      real :: cape,cin,zlcl,zlfc,zel,psource,tsource,thsource,qvsource

      integer :: k1,k2
      real :: lr1,ths1
      real :: tem

      include 'mpif.h'

!------------------------------------------------------------------

      if(dowr) write(outfile,*) 'Inside BASE'

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc  Start definition of base state sounding (isnd opton)  cccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      psurf = 0.0
      tsurf = 0.0
      qsurf = 0.0

!-----------------------------------------------------------------------
!  isnd = 1
!  Dry adiabatic base state

      IF(isnd.eq.1)THEN

        ! Set these two variables for dry adiabatic sounding

        th_sfc   =   300.0   ! Potential temperature of atmosphere (K)
        pi_sfc   =   1.0     ! Exner function at surface

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          th0(i,j,k)=th_sfc
          pi0(i,j,k)=pi_sfc-g*zh(i,j,k)/(cp*th_sfc)
          prs0(i,j,k)=p00*(pi0(i,j,k)**cpdrd)
        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          qv0(i,j,k)=0.0
          rh0(i,j,k)=0.0
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  isnd = 2
!  Dry isothermal base state
            
      ELSEIF(isnd.eq.2)THEN
          
        ! Set these two variables for dry isothermal sounding

        t_sfc    =   250.0    ! Temperature of atmosphere (K)
        prs_sfc  =   p00      ! Pressure at surface (Pa)


        hs=rd*t_sfc/g        ! scale height of atmosphere

        do k=kb,ke
        do j=jb,je
        do i=ib,ie

          ! calculate pressure field
          prs0(i,j,k)=prs_sfc*EXP(-zh(i,j,k)/hs)

          ! using the pressure field, calculate the exner pressure
          pi0(i,j,k)=(prs0(i,j,k)/p00)**(rd/cp)

          ! using exner pressure, determine the potential temperature
          th0(i,j,k)=t_sfc/pi0(i,j,k)

        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          qv0(i,j,k)=0.0
          rh0(i,j,k)=0.0
        enddo
        enddo
        enddo


!-----------------------------------------------------------------------
!  isnd = 3
!  Dry, constant dT/dz sounding.
!  Lapse rate of 0.0065  =  standard dry atmosphere.

      ELSEIF(isnd.eq.3)THEN

        ! Set these three variables for dry constant lapse rate sounding

        th_sfc   =  300.0     ! theta at surface (K)
        prs_sfc  =  p00       ! pressure at surface (Pa)
        lapse    =  0.0065    ! dT/dz (K m^-1)

        do k=kb,ke
        do j=jb,je
        do i=ib,ie

          ! Calculate the temperature using the specified lapse rate.
          t0(i,j,k)=th_sfc-lapse*zh(i,j,k)

          ! Calculate the pressure from the temperature field.
          prs0(i,j,k)=p00*(t0(i,j,k)/th_sfc)**(g/(lapse*rd))

          ! Calculate the exner pressure from the pressure field.
          pi0(i,j,k)=(prs0(i,j,k)/p00)**(rd/cp)

          ! Calculate the theta field from temperature and the 
          ! specified lapse rate.
          th0(i,j,k)=th_sfc*(t0(i,j,k)/th_sfc)**(1-g/(lapse*cp))

        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          qv0(i,j,k)=0.0
          rh0(i,j,k)=0.0
        enddo
        enddo
        enddo


!------------------------------------------------------------------
!  isnd = 4
!  Saturated, neutrally-stable sounding for moist benchmark simulation.
!  reference:  Bryan and Fritsch, 2002, MWR, 130, 2917-2928.

      ELSEIF(isnd.eq.4)THEN

        ! these two parameters define the sounding

        thec_mb   =   320.0     ! wet equivalent potential temp (K)
        qt_mb     =   0.020     ! total water mixing ratio (unitless)


        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Saturated, neutrally-stable sounding'
        if(dowr) write(outfile,*) '    thec,qt=',thec_mb,qt_mb
        if(dowr) write(outfile,*)

      do j=jb,je
      do i=ib,ie

        !! First guesses at lowest model level
        prs0(i,j,1)=100000.
        qv1=qt_mb
        tmp=thec_mb
        tlast=tmp

        do n=1,20
          thbar=(tmp*(p00/prs0(i,j,1))**(rd/cp))*(1.0+qv1*reps)/(1.0+qt_mb)
          pi0(i,j,1)=1.0-g*zh(i,j,1)/(cp*thbar)
          prs0(i,j,1)=p00*(pi0(i,j,1)**(cpdrd))
          qv1=rslf(prs0(i,j,1),tmp)
          tmp=thec_mb*((prs0(i,j,1)/(1.0+qv1/eps)/p00)**(rd/(cp+cpl*qt_mb)))  &
             /exp((lv1-lv2*tmp)*qv1/(tmp*(cp+cpl*qt_mb)))
          tmp=tlast+0.3*(tmp-tlast)
          tlast=tmp
        enddo
 
        th0(i,j,1)=tmp/pi0(i,j,1)
        rh0(i,j,1)=1.0
        qv0(i,j,1)=rslf(prs0(i,j,1),tmp)
        t1=tmp
        qv1=qv0(i,j,1)
        th1=th0(i,j,1)

        do k=2,nk
          tlast=th1
          t2=t1
          qv2=qv1
          th2=th1
          n=0
          pi0(i,j,k)=pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*th2)
100       continue
            n=n+1
            th2=tlast
            t2=th2*pi0(i,j,k)
            thbar=0.5*( th1*(1.0+qv1*reps)/(1.0+qt_mb)    &
                       +th2*(1.0+qv2*reps)/(1.0+qt_mb) )
            pi0(i,j,k)=pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*thbar)
            prs0(i,j,k)=p00*(pi0(i,j,k)**(cpdrd))
            qv2=rslf(prs0(i,j,k),t2)
            t2=thec_mb*((prs0(i,j,k)/(1.0+qv2/eps)/p00)**(rd/(cp+cpl*qt_mb))) &
             /exp((lv1-lv2*t2)*qv2/(t2*(cp+cpl*qt_mb)))
            th2=t2/pi0(i,j,k)

            if(n.gt.50.and.dowr) write(outfile,*) n,th2
            if(abs(th2-tlast).gt.0.0001 .and. n.lt.100)then
              tlast=tlast+0.3*(th2-tlast)
              go to 100
            elseif(n.ge.100)then
              if(dowr) write(outfile,*) '  stuck in loop!'
              call stopcm1
            endif
 
          th0(i,j,k)=th2
          qv0(i,j,k)=rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k))
          rh0(i,j,k)=1.0
 
          th1=th2
          qv1=qv2
          t1=t2
 
        enddo

        do k=1,nk
          qc0(i,j,k)=qt_mb-qv0(i,j,k)
        enddo

      enddo
      enddo


!-----------------------------------------------------------------------
!  isnd = 5
!  Weisman-Klemp analytic sounding.
!  reference:  Weisman and Klemp, 1982, MWR, 110, 504-520.

      ELSEIF(isnd.eq.5)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  WK sounding'

!  variables related to Weisman-Klemp analytic sounding
        z_trop   = 12000.0      ! height of tropopause (m)
        th_trop  = 343.0        ! theta at tropopause (K)
        t_trop   = 213.0        ! temp at tropopause (K)
        th_sfc   = 300.0        ! theta at surface (K)
        prs_sfc  = 100000.0     ! pressure at surface (Pa)
        qv_pbl   = 0.014        ! constant value of mixing ratio in PBL

!--------------

        pi_sfc  = (prs_sfc/p00)**(rd/cp)
        qv_sfc  = rslf(prs_sfc,th_sfc*pi_sfc)
        thv_sfc = th_sfc*(1.0+qv_sfc*reps)/(1.0+qv_sfc)

      do j=jb,je
      do i=ib,ie

        do k=kb,ke
          rh0(i,j,k)=0.0
        enddo

        do k=1,nk
          if(zh(i,j,k).le.z_trop)then
            th0(i,j,k)=th_sfc+(th_trop-th_sfc)*((zh(i,j,k)/z_trop)**1.25)
            if(imoist.eq.1) rh0(i,j,k)=1.0-0.75*((zh(i,j,k)/z_trop)**1.25)
          else
            th0(i,j,k)=th_trop*exp((g/(t_trop*cp))*(zh(i,j,k)-z_trop))
            if(imoist.eq.1) rh0(i,j,k)=0.25
          endif
        enddo

        th0(i,j,0)=th0(i,j,1)
        th0(i,j,nk+1)=th0(i,j,nk)


!  Get pressure, temperature, and mixing ratio using hydrostatic eqt.

        do k=kb,ke
          qv0(i,j,k)=0.0
        enddo

        do n=1,20
! virtual potential temperature
          do k=kb,ke
            thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k))
          enddo

          pi0(i,j,1)=pi_sfc-g*zh(i,j,1)/(cp*0.5*(thv_sfc+thv0(i,j,1)))
          do k=2,nk
            pi0(i,j,k)=pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))/(cp*0.5*(thv0(i,j,k)+thv0(i,j,k-1)))
          enddo

! pressure
          do k=1,nk
            prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
          enddo

! mixing ratio
          do k=1,nk
            qv0(i,j,k)=rh0(i,j,k)*rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k))
            if(qv0(i,j,k).gt.qv_pbl) qv0(i,j,k)=qv_pbl
          enddo

        enddo

        do k=1,nk
          rh0(i,j,k)=qv0(i,j,k)/(rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k)))
        enddo

      enddo
      enddo

!------------------------------------------------------------------
!  isnd = 6

      ELSEIF(isnd.eq.6)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' isnd = 6 code has been removed (for now)'
        if(dowr) write(outfile,*)
        call stopcm1

!------------------------------------------------------------------
!  isnd = 7
!  Read base-state sounding from an external text file.
!  Assumes file name = input_sounding
!  NOTE:  for isnd=7, iwnd is ignored.
!
!  The format is the same as that for the WRF Model, as well as 
!  for the Klemp-Wilhelmson Model.  Format is:
!
!  1 line header containing: sfc pres (mb)    sfc theta (K)    sfc qv (g/kg)
!  nsnd lines of:  zheight (m)    theta (K)   qv (g/kg)    u (m/s)    v (m/s)
!
!  (Thanks to Leigh Orf, Central Michigan University, for contributing this
!   code.)
!  Note:  This now works with terrain.  (GHB, 061011)
!

      ELSEIF(isnd.eq.7)then

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' Reading sounding from external file, input_sounding'
        if(dowr) write(outfile,*)

        nmax = 100000 !should be enough!

        allocate(   zsnd(nmax) )
        allocate(  thsnd(nmax) )
        allocate(  qvsnd(nmax) )
        allocate(   usnd(nmax) )
        allocate(   vsnd(nmax) )
        allocate( thvsnd(nmax) )
        allocate(  pisnd(nmax) )
        allocate(   psnd(nmax) )
        allocate(   tsnd(nmax) )
        allocate(  rhsnd(nmax) )

        allocate( thinterp(nk) )
        allocate( qvinterp(nk) )
        allocate(  uinterp(nk) )
        allocate(  vinterp(nk) )
        allocate(  pinterp(nk) )
        allocate(  tinterp(nk) )
        allocate( rhinterp(nk) )

        open(unit=40,file='input_sounding',status='old')

        ! read surface parameters:
        read(40,*) p_sfc, th_sfc, qv_sfc
        if(dowr) write(outfile,*) ' p_sfc, th_sfc, qv_sfc = ',p_sfc, th_sfc, qv_sfc

        if(imoist.eq.0) qv_sfc = 0.0

        p_sfc = p_sfc * 100.0
        ! put qv in g/g
        qv_sfc = qv_sfc / 1000.0
        pi_sfc  = (p_sfc/p00)**(rd/cp)
        thv_sfc = th_sfc*(1.0+qv_sfc*reps)/(1.0+qv_sfc)

        psurf  = p_sfc
        thsurf = th_sfc
        tsurf  = th_sfc * pi_sfc
        qsurf  = qv_sfc

        zsnd(1) = 0.0
        thsnd(1) = th_sfc
        qvsnd(1) = qv_sfc

        ! now, read entire sounding until end of file is discovered
        nsnd=1
        do k=1,nmax
          read(40,*,end=445) zsnd(k+1),thsnd(k+1),qvsnd(k+1),usnd(k+1),vsnd(k+1)
          ! put qv in g/g
          qvsnd(k+1) = qvsnd(k+1)/1000.0
          nsnd=nsnd+1
        enddo
445     continue
        if(dowr) write(outfile,*) '  Found ',nsnd,'  levels (including surface)'
        if(dowr) write(outfile,*)
        close(unit=40)

        if(imoist.eq.0) qvsnd = 0.0

        usnd(1) = 1.75*usnd(2)-usnd(3)+0.25*usnd(4)
        vsnd(1) = 1.75*vsnd(2)-vsnd(3)+0.25*vsnd(4)

!--------------------------------------------------------------------
!  Added by GHB, 061021:
!  Get thv and prs ... check if qv is too small.  If so, set rh to 5%
!  (This code has no effect on the sounding if qv > 1e-12 everywere)
!  (It was added to deal with the 0 g/kg qv values in the Trier sounding.)
        do k=1,nsnd
          thvsnd(k)=thsnd(k)*(1.0+reps*qvsnd(k))/(1.0+qvsnd(k))
        enddo
        pisnd(1)=pi_sfc
        do k=2,nsnd
          pisnd(k)=pisnd(k-1)-g*(zsnd(k)-zsnd(k-1))   &
                               /(cp*0.5*(thvsnd(k)+thvsnd(k-1)))
        enddo
        do k=1,nsnd
          psnd(k)=p00*(pisnd(k)**(cp/rd))
          tsnd(k)=thsnd(k)*pisnd(k)
        enddo
      if(imoist.eq.1)then
        do k=1,nsnd
          if(qvsnd(k).lt.1.0e-12)then
            if(dowr) write(outfile,*) '  Qv is too small.  Setting rh to 1%.  k,zsnd=',k,zsnd(k)
            qvsnd(k)=0.01*rslf(psnd(k),thsnd(k)*pisnd(k))
          endif
          rhsnd(k)=qvsnd(k)/rslf(psnd(k),tsnd(k))
        enddo
      endif
!--------------------------------------------------------------------

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '    k,z,th,qv:'
        do k=1,nsnd
          if(dowr) write(outfile,*) k,zsnd(k),thsnd(k),1000.0*qvsnd(k)
        enddo
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '    k,u,v:'
        do k=1,nsnd
          if(dowr) write(outfile,*) k,usnd(k),vsnd(k)
        enddo
        if(dowr) write(outfile,*)

! check to make sure sounding levels span computational grid (WICKER)
!!!        061011, GHB:  Commented this out ... we'll use the surface info.
!!!        if (zsnd(1) .gt. zh(1,1,1)) then
!!!          write(*,*) 'zmin of sounding > zmin of grid!'
!!!          write(*,*) 'zmin of sounding = ',zsnd(1)
!!!          write(*,*) 'zmin of grid = ',zh(1,1,1)
!!!          call stopcm1
!!!        endif

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) 'interpolating sounding to grid: '
        if(dowr) write(outfile,*)

      DO j=jb,je
      DO i=ib,ie

        if (zsnd(nsnd) .lt. zh(i,j,nk)) then
          if(dowr) write(outfile,*) 'zmax of sounding < zmax of grid!'
          if(dowr) write(outfile,*) 'zmax of sounding = ',zsnd(nsnd)
          if(dowr) write(outfile,*) 'zmax of grid  = ',zh(i,j,nk)
          call stopcm1
        endif

        DO k=1,nk

            kk = 1
            do while( zsnd(kk) .lt. zh(i,j,k) )
              kk = kk+1
            enddo
            kdn = kk-1
            kup = kk

            zu=0.5*(zh(max(ib,i-1),j,k)+zh(i,j,k))
            interp_frac = (   zu        - zsnd(kdn) )   &
                        / ( zsnd( kup ) - zsnd(kdn) )
            uinterp(k) =  usnd(kdn) + ( usnd(kup)- usnd(kdn))*interp_frac

            zv=0.5*(zh(i,max(jb,j-1),k)+zh(i,j,k))
            interp_frac = (   zv        - zsnd(kdn) )   &
                        / ( zsnd( kup ) - zsnd(kdn) )
            vinterp(k) =  vsnd(kdn) + ( vsnd(kup)- vsnd(kdn))*interp_frac

!!!            ! if this is first grid point, utilize surface values of th,qv:
!!!            if( k.eq.1 ) kdn = 1

            interp_frac = (   zh(i,j,k) - zsnd(kdn) )   &
                        / ( zsnd( kup ) - zsnd(kdn) )
            thinterp(k) = thsnd(kdn) + (thsnd(kup)-thsnd(kdn))*interp_frac
            qvinterp(k) = qvsnd(kdn) + (qvsnd(kup)-qvsnd(kdn))*interp_frac
             pinterp(k) =  psnd(kdn) + ( psnd(kup)- psnd(kdn))*interp_frac
             tinterp(k) =  tsnd(kdn) + ( tsnd(kup)- tsnd(kdn))*interp_frac
            rhinterp(k) = rhsnd(kdn) + (rhsnd(kup)-rhsnd(kdn))*interp_frac

            if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) '       ',zsnd(kdn),zh(i,j,k),zsnd(kup),interp_frac

        ENDDO
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*)

        do k=1,nk
           u0(i,j,k) =  uinterp(k)
           v0(i,j,k) =  vinterp(k)
          qv0(i,j,k) = qvinterp(k)
          th0(i,j,k) = thinterp(k)
         prs0(i,j,k) =  pinterp(k)
           t0(i,j,k) =  tinterp(k)
          rh0(i,j,k) = rhinterp(k)
        enddo

        ! get pi0 and prs0 from thv0, using hydrostatic equation

        do k=1,nk
          ! get qv from linear interpolation of rh:
          qv0(i,j,k) = rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
        enddo

        do k=1,nk
          thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k))
        enddo

        pi0(i,j,1)=pi_sfc-g*zh(i,j,1)/(cp*0.5*(thv_sfc+thv0(i,j,1)))
        do k=2,nk
          pi0(i,j,k)=pi0(i,j,k-1)-g*(zh(i,j,k)-zh(i,j,k-1))   &
                                   /(cp*0.5*(thv0(i,j,k)+thv0(i,j,k-1)))
        enddo

        do k=1,nk
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
        enddo

        ! rh, just in case we want/need it later

      if(imoist.eq.1)then
        do k=1,nk
          rh0(i,j,k)=qv0(i,j,k)/(rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k)))
        enddo
      endif

      ENDDO    ! enddo for i loop
      ENDDO    ! enddo for j loop


        ! deallocate temporary 1D arrays
        deallocate(   zsnd )
        deallocate(  thsnd )
        deallocate(  qvsnd )
        deallocate(   usnd )
        deallocate(   vsnd )
        deallocate( thvsnd )
        deallocate(  pisnd )
        deallocate(   psnd )
        deallocate(   tsnd )
        deallocate(  rhsnd )
        deallocate( thinterp )
        deallocate( qvinterp )
        deallocate(  uinterp )
        deallocate(  vinterp )
        deallocate(  pinterp )
        deallocate(  tinterp )
        deallocate( rhinterp )


!-----------------------------------------------------------------------
!  isnd = 8
!  Dry, constant d(theta)/dz sounding

      ELSEIF(isnd.eq.8)THEN

        ! Set these three variables for dry, constant d(theta)/dz sounding

        th_sfc   =  300.0     ! theta at surface (K)
        pi_sfc   =    1.0     ! Exner function at surface
        lapse    =  0.0035    ! potential temperature lapse rate (K/m)

        do k=kb,ke
        do j=jb,je
        do i=ib,ie

          ! Calculate theta using the specified lapse rate.
          th0(i,j,k)=th_sfc+lapse*zh(i,j,k)

          ! Calculate pi from theta
          pi0(i,j,k)=pi_sfc-(g/(cp*lapse))*alog(th0(i,j,k)/th_sfc)

          ! Calculate pressure from pi
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))

          ! Calculate temperature from theta and pi
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)

        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          qv0(i,j,k)=0.0
          rh0(i,j,k)=0.0
        enddo
        enddo
        enddo

!------------------------------------------------------------------
!  isnd = 9
!  Constant Brunt-Vaisala frequency

      ELSEIF(isnd.eq.9)then

        ns1 =  0.0001
        ns2 =  0.0000

        zl1 = 40000.0

      do j=jb,je
      do i=ib,ie
        th_sfc   =  288.00
        pi_sfc   =    1.0
        if(zh(i,j,1).lt.zl1)then
          zsfc = 0.0
          ns   = ns1
        else
          if(abs(ns1).lt.1.0e-6)then
            pi_sfc=pi_sfc-g*zl1/(cp*th_sfc)
          else
            pi_sfc=pi_sfc+g*g/(cp*ns1*th_sfc)*(exp(-ns1*zl1/g)-1.0)
            th_sfc=th_sfc*exp(ns1*zl1/g)
          endif
          zsfc = zl1
          ns   = ns2
        endif
        do k=1,nk
          qv0(i,j,k)=0.0
          qc0(i,j,k)=0.0
          rh0(i,j,k)=0.0
          if(abs(ns).lt.1.0e-6)then
            thv0(i,j,k)=th_sfc
            pi0(i,j,k)=pi_sfc-g*(zh(i,j,k)-zsfc)/(cp*th_sfc)
          else
            thv0(i,j,k)=th_sfc*exp(ns*(zh(i,j,k)-zsfc)/g)
            pi0(i,j,k)=pi_sfc+g*g/(cp*ns*th_sfc)   &
                           *(exp(-ns*(zh(i,j,k)-zsfc)/g)-1.0)
          endif
          prs0(i,j,k)=p00*(pi0(i,j,k)**cpdrd)
          th0(i,j,k)=thv0(i,j,k)
          do n=1,20
            t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
            qv0(i,j,k)=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
            th0(i,j,k)=thv0(i,j,k)*(1.0+qv0(i,j,k)+qc0(i,j,k))/(1.0+reps*qv0(i,j,k))
!!!            if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,n,th0(i,j,k),qv0(i,j,k)
          enddo
          if(zh(i,j,k+1).gt.zl1.and.zsfc.lt.1.0)then
            if(abs(ns1).lt.1.0e-6)then
              pi_sfc=pi_sfc-g*zl1/(cp*th_sfc)
            else
              pi_sfc=pi_sfc+g*g/(cp*ns1*th_sfc)*(exp(-ns1*zl1/g)-1.0)
              th_sfc=th_sfc*exp(ns1*zl1/g)
            endif
            zsfc=zl1
            ns=ns2
          endif
        enddo

      enddo
      enddo

!------------------------------------------------------------------
!  isnd = 10
!  Moist, constant Brunt-Vaisala frequency
!    (assuming reversible moist microphysics)

      ELSEIF(isnd.eq.10)then

        if(dowr) write(outfile,*) '  isnd = 10'

        qcval    =   0.005
        qtval    =   0.00000

        ns1 =  0.0
        ns2 =  4.0e-4

        zl1 = 111700.0

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        ! DEFINE RH here:  saturated (RH = 100%) is assumed by default
        rh0(i,j,k)=1.0
        qc0(i,j,k)=qcval
      enddo
      enddo
      enddo

      if(dowr) write(outfile,*)
      do j=jb,je
      do i=ib,ie
        th_sfc   =  288.0
        pi_sfc   =    1.0
        ns=ns1
        p_sfc=p00*(pi_sfc**(cp/rd))
        t_sfc=th_sfc*pi_sfc
        qv_sfc=rh0(i,j,1)*rslf(p_sfc,t_sfc)
        ql_sfc=qcval
!!!        ql_sfc=qtval-qv_sfc
        thv_sfc=th_sfc*(1.0+reps*qv_sfc)/(1.0+qv_sfc+ql_sfc)
        thn=th_sfc
        thlast=thn
      !---------- do k=1 first ----------!
        k=1
        n=1
        qv0(i,j,k)=qv_sfc
        qc0(i,j,k)=ql_sfc
 571      continue
          th0(i,j,k)=thn
          thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))   &
                                /(1.0+qv0(i,j,k)+qc0(i,j,k))
          pi0(i,j,k)=pi_sfc-g*zh(i,j,k)/(cp*0.5*(thv0(i,j,k)+thv_sfc))
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          qvs=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
          qv0(i,j,k)=qvs
          qc0(i,j,k)=qcval
!!!          qc0(i,j,k)=qtval-qv0(i,j,k)
          qtm=qv_sfc+ql_sfc
          qtp=qv0(i,j,k)+qc0(i,j,k)
          tavg=0.5*( t0(i,j,k)+t_sfc )
          lhv=lv1-lv2*tavg
          qvavg=0.5*( qv_sfc + qv0(i,j,k) )
          qlavg=0.5*( ql_sfc + qc0(i,j,k) )
          qtavg=qvavg+qlavg
          desdt=17.67*(273.15-29.65)/((tavg-29.65)**2)
          drdt=17.67*(273.15-29.65)*qvavg/((tavg-29.65)**2)
        if(eqtset.eq.2)then
          gamma=g*(1.0+qtavg)*(1.0+lhv*qvavg/rd/tavg)   &
               /( cp+cpv*qvavg+cpl*qlavg          &
                 +lhv*drdt )
        else
          gamma=g*(1.0+qtavg)*(rd/(rd+rv*qvavg)+lhv*qvavg/rd/tavg)   &
               /( cp+lhv*(1.0+qvavg*reps)*qvavg*desdt )
        endif
          tp=t_sfc*exp( zh(i,j,k)*(                     &
                ((ns/g)+alog((1.0+qtp)/(1.0+qtm))/zh(i,j,k))  &
                   /(1.0+tavg*drdt/(eps+qvavg))-gamma/tavg                  &
                                                      ) )
          thn=tp/pi0(i,j,k)

          n=n+1

          if(n.gt.180.and.dowr) write(outfile,*) n,thn
          if(abs(thn-thlast).gt.0.0001 .and. n.lt.200)then
            thn=thlast+0.3*(thn-thlast)
            thlast=thn
            go to 571
          elseif(n.ge.200)then
            if(dowr) write(outfile,*) '  stuck in loop (111)!'
            call stopcm1
          endif

        th0(i,j,k)=thn
        prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
        t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
        qvs=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
        qv0(i,j,k)=qvs
        qc0(i,j,k)=qcval
!!!        qc0(i,j,k)=qtval-qv0(i,j,k)
        if(qc0(i,j,1).lt.0.0)then
          print *,'  calling stopcm1 (222) '
          call stopcm1
        endif
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,n,th0(i,j,k),prs0(i,j,k)
      !---------- done with k=1, do all other levels ----------!
        qvs=qv0(i,j,1)
        ns=ns1
        do k=2,nk
        tflag=0
        do niter=1,2
        if(tflag.lt.2)then
          if(zh(i,j,k).gt.zl1.and.ns.eq.ns1.and.tflag.eq.0.and.zh(i,j,k-1).lt.zl1)then
            tflag=1
            delz=zl1-zh(i,j,k-1)
            pim=pi0(i,j,k-1)
            thm=th0(i,j,k-1)
            thvm=thv0(i,j,k-1)
            qtm=qv0(i,j,k-1)+qc0(i,j,k-1)+ql0(i,j,k-1)
            qvm=qv0(i,j,k-1)
            qlm=qc0(i,j,k-1)
            qim=ql0(i,j,k-1)
            tm=t0(i,j,k-1)
            pm=prs0(i,j,k-1)
          elseif(tflag.eq.1)then
            tflag=0
            ns=ns2
            delz=zh(i,j,k)-zl1
            pim=pi0(i,j,k)
            thm=th0(i,j,k)
            thvm=thv0(i,j,k)
            qtm=qv0(i,j,k)+qc0(i,j,k)+ql0(i,j,k)
            qvm=qv0(i,j,k)
            qlm=qc0(i,j,k)
            qim=ql0(i,j,k)
            tm=t0(i,j,k)
            pm=prs0(i,j,k)
          else
            tflag=2
            delz=zh(i,j,k)-zh(i,j,k-1)
            pim=pi0(i,j,k-1)
            thm=th0(i,j,k-1)
            thvm=thv0(i,j,k-1)
            qtm=qv0(i,j,k-1)+qc0(i,j,k-1)+ql0(i,j,k-1)
            qvm=qv0(i,j,k-1)
            qlm=qc0(i,j,k-1)
            qim=ql0(i,j,k-1)
            tm=t0(i,j,k-1)
            pm=prs0(i,j,k-1)
          endif
!!!          t0(i,j,k)=tm
          th0(i,j,k)=thm
          qv0(i,j,k)=qvm
          qc0(i,j,k)=qlm
          ql0(i,j,k)=qim
          thv0(i,j,k)=thvm
!!!          thlast=tm*pim
          thlast=thm
          n=0
          thn=thlast
          thlast=thn
 572      continue
            th0(i,j,k)=thn
            thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))   &
                                  /(1.0+qv0(i,j,k)+qc0(i,j,k)+ql0(i,j,k))
            pi0(i,j,k)=pim-g*delz/(cp*0.5*(thv0(i,j,k)+thvm))
            prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
            t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          if(iice.eq.0)then
            qvs=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
            qv0(i,j,k)=qvs
            qc0(i,j,k)=qcval
!!!            qc0(i,j,k)=qtval-qv0(i,j,k)
            qtp=qv0(i,j,k)+qc0(i,j,k)
            pavg=0.5*(prs0(i,j,k-1)+prs0(i,j,k))
            tavg=0.5*( t0(i,j,k)+tm )
            qvavg=rh0(i,j,k)*rslf(pavg,tavg)
            qlavg=0.5*( qlm + qc0(i,j,k) )
            qtavg=0.5*( qtm + qtp )
            drdt=17.67*(273.15-29.65)*qvavg/((tavg-29.65)**2)
            lhv=lv1-lv2*tavg
            cpml=cp+cpv*qvavg+cpl*qlavg
          else
            qvl=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
            qvi=rh0(i,j,k)*rsif(prs0(i,j,k),t0(i,j,k))
            fliq=max(min((t0(i,j,k)-t00k)*rt0,1.0),0.0)
            fice=1.0-fliq
            qvs=fliq*qvl+fice*qvi
            qv0(i,j,k)=qvs
            qc0(i,j,k)=fliq*qcval
            ql0(i,j,k)=fice*qcval
            qtp=qv0(i,j,k)+qc0(i,j,k)+ql0(i,j,k)
            pavg=0.5*(prs0(i,j,k-1)+prs0(i,j,k))
            tavg=0.5*( t0(i,j,k)+tm )
            qvl=rh0(i,j,k)*rslf(pavg,tavg)
            qvi=rh0(i,j,k)*rsif(pavg,tavg)
            fliq=max(min((tavg-t00k)*rt0,1.0),0.0)
            fice=1.0-fliq
            qvavg=fliq*qvl+fice*qvi
            qlavg=0.5*( qlm + qc0(i,j,k) )
            qtavg=0.5*( qtm + qtp )
            drdt=fliq*17.67*(273.15-29.65)*qvl/((tavg-29.65)**2)    &
                +fice*21.8745584*(273.15-7.66)*qvi/((tavg-7.66)**2)
            if(tavg.gt.t00k.and.tavg.lt.t0k)then
              drdt=drdt+(qvl-qvi)*rt0
            endif
            lhv=fliq*(lv1-lv2*tavg)+fice*(ls1-ls2*tavg)
            cpml=cp+cpv*qvavg+cpl*qlavg+cpi*(qtavg-qlavg-qvavg)
          endif
          if(eqtset.eq.2)then
            gamma=g*(1.0+qtavg)*(1.0+lhv*qvavg/rd/tavg)   &
                 /( cpml+lhv*drdt )
          else
            gamma=g*(1.0+qtavg)*(rd/(rd+rv*qvavg)+lhv*qvavg/rd/tavg)   &
                 /( cp+lhv*(1.0+qvavg*reps)*qvavg*desdt )
          endif
            tp=tm*exp( delz*(                                 &
                  ((ns/g)+alog((1.0+qtp)/(1.0+qtm))/delz)     &
                     /(1.0+tavg*drdt/(eps+qvavg))-gamma/tavg  &
                                                        ) )
            thn=tp/pi0(i,j,k)

            n=n+1

            if(n.gt.180.and.dowr) write(outfile,*) n,tp,tm,delz
            if(abs(thn-thlast).gt.0.0001 .and. n.lt.200)then
              thn=thlast+0.3*(thn-thlast)
              thlast=thn
              go to 572
            elseif(n.ge.200)then
              if(dowr) write(outfile,*) '  stuck in loop (333)!'
              call stopcm1
            endif

          t0(i,j,k)=tp
          th0(i,j,k)=thn
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          if(iice.eq.0)then
            qvs=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
            qv0(i,j,k)=qvs
          else
            qvl=rh0(i,j,k)*rslf(prs0(i,j,k),t0(i,j,k))
            qvi=rh0(i,j,k)*rsif(prs0(i,j,k),t0(i,j,k))
            fliq=max(min((t0(i,j,k)-t00k)*rt0,1.0),0.0)
            fice=1.0-fliq
            qvs=fliq*qvl+fice*qvi
            qv0(i,j,k)=qvs
            qc0(i,j,k)=fliq*qcval
            ql0(i,j,k)=fice*qcval
          endif
          if(qc0(i,j,1).lt.0.0)then
            print *,'  calling stopcm1 (444) '
            call stopcm1
          endif
!!!          if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,n,zh(i,j,k),th0(i,j,k),prs0(i,j,k)
          if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,n,t0(i,j,k),fliq,fice
        endif
        enddo  ! enddo for iteration loop
        enddo  ! enddo for k loop

      enddo  ! enddo for i loop
      enddo  ! enddo for j loop
      if(dowr) write(outfile,*)

!------------------------------------------------------------------
!  constant theta-e, saturated
!  Reference:  Bryan and Rotunno, 2009, JAS, v10, pp. 3042-3060

      ELSEIF(isnd.eq.11)THEN

        IF(imoist.eq.0)THEN
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  isnd=11 requires imoist=1'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  stopping model ...'
          if(dowr) write(outfile,*)
          call stopcm1
        ENDIF

        zl1 = 15000.0    ! tropopause height

        ttype =  2    ! 1 = reversible
                      ! 2 = pseudoadiabatic

        psurf = 101510.0    ! surface pressure (Pa)

      DO j=jb,je
      DO i=ib,ie

        ! use psurf:
        pi_sfc   =  (psurf*rp00)**rovcp
        ! use tsk from namelist.input, make 3 K cooler:
        th_sfc   =  ( tsk0 - 3.0 )/pi_sfc

        t_sfc=th_sfc*pi_sfc
        p_sfc=p00*(pi_sfc**(cp/rd))
        qv_sfc=rslf(p_sfc,t_sfc)
        qt_sfc=qv_sfc
        thv_sfc=th_sfc*(1.0+reps*qv_sfc)/(1.0+qv_sfc)

        the_sfc=getthe(ttype,p_sfc,t_sfc,qv_sfc,qt_sfc)
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) '  p_sfc,t_sfc,qv_sfc = ',p_sfc,t_sfc,qv_sfc
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) '  the_sfc,qt_sfc = ',the_sfc,qt_sfc

        th_sfc=0.0
        qv_sfc=0.0
        call revthe(ttype,the_sfc,p_sfc,qt_sfc,t_sfc,qv_sfc)
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) '  p_sfc,t_sfc,qv_sfc = ',p_sfc,t_sfc,qv_sfc

        pi1=pi_sfc
        z1=0.0
        t1=t_sfc
        thv1=thv_sfc
        thv2=thv1
        flag=0

        do k=1,nk

          z2=zh(i,j,k)
          t2=t1
          tlast=0.0
          n=0

        if(z2.le.zl1)then
          ! troposphere
          theq=the_sfc
          ! iterate:
          do while( abs(t2-tlast).gt.0.0001 )
            tlast=t2
            n=n+1
            pi2=pi1-g*(z2-z1)/(cp*0.5*(thv1+thv2))
            p2=p00*(pi2**(cp/rd))
            call revthe(ttype,theq,p2,qt_sfc,t2,qv2)
            th2=t2/pi2
            if(ttype.eq.1)then
              ql2=qt_sfc-qv2
            elseif(ttype.eq.2)then
              ql2=0.0
            endif
            thv2=th2*(1.0+reps*qv2)/(1.0+qv2+ql2)
          enddo
        else
          ! stratosphere
          ns=4.0e-4
          if(flag.eq.0)then
            ! first time in this section ... get sfc params
            !-------------------
            z2=zl1
            t2=t1
            tlast=0.0
            n=0
            theq=the_sfc
            do while( abs(t2-tlast).gt.0.0001 )
              tlast=t2
              n=n+1
              pi2=pi1-g*(z2-z1)/(cp*0.5*(thv1+thv2))
              p2=p00*(pi2**(cp/rd))
              call revthe(ttype,theq,p2,qt_sfc,t2,qv2)
              th2=t2/pi2
              if(ttype.eq.1)then
                ql2=qt_sfc-qv2
              elseif(ttype.eq.2)then
                ql2=0.0
              endif
              thv2=th2*(1.0+reps*qv2)/(1.0+qv2+ql2)
            enddo
            flag=1
            th_sfc=thv2
            pi_sfc=pi2
            zsfc=z2
            !-------------------
            z2=zh(i,j,k)
            t2=t1
            tlast=0.0
            n=0
          endif
          if(abs(ns).lt.1.0e-6)then
            thv2=th_sfc
            pi2=pi_sfc-g*(z2-zsfc)/(cp*th_sfc)
          else
            thv2=th_sfc*exp(ns*(z2-zsfc)/g)
            pi2=pi_sfc+g*g/(cp*ns*th_sfc)   &
                           *(exp(-ns*(z2-zsfc)/g)-1.0)
          endif
          p2=p00*(pi2**cpdrd)
          th2=thv2
          do n=1,20
            t2=th2*pi2
            qv2=rslf(p2,t2)
            if(ttype.eq.1)then
              ql2=qt_sfc-qv2
            elseif(ttype.eq.2)then
              ql2=0.0
            endif
            th2=thv2*(1.0+qv2+ql2)/(1.0+reps*qv2)
          enddo
        endif
!!!          if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) n,p2,th2

          t1=t2
          thv1=thv2
          z1=z2
          pi1=pi2

          pi0(i,j,k)=pi2
          prs0(i,j,k)=p2
          t0(i,j,k)=t2
          th0(i,j,k)=th2
          thv0(i,j,k)=thv2
          rh0(i,j,k)=1.0
          qv0(i,j,k)=qv2
          qc0(i,j,k)=ql2

        enddo

      ENDDO
      ENDDO

!------------------------------------------------------------------
!  PBL simulation:  assumed dry

      ELSEIF(isnd.eq.12)THEN

        pi_sfc = 1.0
        th_sfc = 300.0
        zl1    = 960.0
        lapse  = 0.010

        pisfc = pi_sfc-g*zl1/(cp*th_sfc)

        do k=1,nk
        do j=jb,je
        do i=ib,ie
          IF(zh(i,j,k).le.zl1)THEN
            th0(i,j,k) = th_sfc
            pi0(i,j,k)=pi_sfc-g*zh(i,j,k)/(cp*th_sfc)
          ELSE
            th0(i,j,k) = th_sfc+lapse*(zh(i,j,k)-zl1)
            pi0(i,j,k)=pisfc-(g/(cp*lapse))*alog(th0(i,j,k)/th_sfc)
          ENDIF
          prs0(i,j,k)=p00*(pi0(i,j,k)**(cp/rd))
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          rh0(i,j,k)=0.0
          qv0(i,j,k)=0.0
          if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,zh(i,j,k),th0(i,j,k),prs0(i,j,k)
        enddo
        enddo
        enddo

!------------------------------------------------------------------
!  from Klemp (2011)

      ELSEIF(isnd.eq.13)then

        ns1 =  0.0001
        ns2 =  0.0004
        ns3 =  0.0001

        zl1 =  2000.0
        zl2 =  3000.0

      do j=jb,je
      do i=ib,ie
        th_sfc   =  288.00
        pi_sfc   =    1.0
        if(zh(i,j,1).lt.zl1)then
          zsfc = 0.0
          ns   = ns1
        else
          if(abs(ns1).lt.1.0e-6)then
            pi_sfc=pi_sfc-g*zl1/(cp*th_sfc)
          else
            pi_sfc=pi_sfc+g*g/(cp*ns1*th_sfc)*(exp(-ns1*zl1/g)-1.0)
            th_sfc=th_sfc*exp(ns1*zl1/g)
          endif
          zsfc = zl1
          ns   = ns2
        endif
        do k=1,nk
          if(abs(ns).lt.1.0e-6)then
            th0(i,j,k)=th_sfc
            pi0(i,j,k)=pi_sfc-g*(zh(i,j,k)-zsfc)/(cp*th_sfc)
          else
            th0(i,j,k)=th_sfc*exp(ns*(zh(i,j,k)-zsfc)/g)
            pi0(i,j,k)=pi_sfc+g*g/(cp*ns*th_sfc)   &
                           *(exp(-ns*(zh(i,j,k)-zsfc)/g)-1.0)
          endif
          prs0(i,j,k)=p00*(pi0(i,j,k)**cpdrd)
          thv0(i,j,k)=th0(i,j,k)
          qv0(i,j,k)=0.0
          rh0(i,j,k)=0.0
          t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
          if(zh(i,j,k+1).gt.zl1.and.zsfc.lt.1.0)then
            if(abs(ns1).lt.1.0e-6)then
              pi_sfc=pi_sfc-g*zl1/(cp*th_sfc)
            else
              pi_sfc=pi_sfc+g*g/(cp*ns1*th_sfc)*(exp(-ns1*zl1/g)-1.0)
              th_sfc=th_sfc*exp(ns1*zl1/g)
            endif
            zsfc=zl1
            ns=ns2
          elseif(zh(i,j,k+1).gt.zl2.and.zsfc.lt.zl2)then
            if(abs(ns1).lt.1.0e-6)then
              pi_sfc=pi_sfc-g*(zl2-zl1)/(cp*th_sfc)
            else
              pi_sfc=pi_sfc+g*g/(cp*ns2*th_sfc)*(exp(-ns2*(zl2-zl1)/g)-1.0)
              th_sfc=th_sfc*exp(ns2*(zl2-zl1)/g)
            endif
            zsfc=zl2
            ns=ns3
          endif
        enddo

      enddo
      enddo

!------------------------------------------------------------------

      ENDIF

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc  End definition of base state sounding (isnd opton)  cccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!------------------------------------------------------------------
!  fill in ghost cells

      call bcs(pi0)
      call bcs(prs0)
      call bcs(th0)
      call bcs(qv0)
      call bcs(qc0)
      call bcs(rh0)

      nf=0
      nu=0
      nv=0
      nw=0
      call comm_3s_start(pi0,sw31,sw32,se31,se32,   &
                             ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(pi0,sw31,sw32,se31,se32,   &
                           ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_start(prs0,sw31,sw32,se31,se32,   &
                             ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(prs0,sw31,sw32,se31,se32,   &
                           ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_start(th0,sw31,sw32,se31,se32,   &
                             ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(th0,sw31,sw32,se31,se32,   &
                           ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_start(qv0,sw31,sw32,se31,se32,   &
                             ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(qv0,sw31,sw32,se31,se32,   &
                           ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_start(qc0,sw31,sw32,se31,se32,   &
                             ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(qc0,sw31,sw32,se31,se32,   &
                           ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_start(rh0,sw31,sw32,se31,se32,   &
                             ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(rh0,sw31,sw32,se31,se32,   &
                           ss31,ss32,sn31,sn32,reqs_s)

    do j=jb,je
    do i=ib,ie

      pi0(i,j,0)=pi0(i,j,1)
      pi0(i,j,nk+1)=pi0(i,j,nk)

      prs0(i,j,0)=prs0(i,j,1)
      prs0(i,j,nk+1)=prs0(i,j,nk)

      th0(i,j,0)=th0(i,j,1)
      th0(i,j,nk+1)=th0(i,j,nk)

      qv0(i,j,0)=qv0(i,j,1)
      qv0(i,j,nk+1)=qv0(i,j,nk)

      qc0(i,j,0)=qc0(i,j,1)
      qc0(i,j,nk+1)=qc0(i,j,nk)

      rh0(i,j,0)=rh0(i,j,1)
      rh0(i,j,nk+1)=rh0(i,j,nk)

!------------------------------------------------------------------
!  check thv0
!  Assumes th0, qv0 are accurate

      if(imoist.eq.1)then
        do k=kb,ke
          thv0(i,j,k)=th0(i,j,k)*(1.0+reps*qv0(i,j,k))/(1.0+qv0(i,j,k)+qc0(i,j,k))
        enddo
      else
        do k=kb,ke
          qv0(i,j,k)=0.0
          rh0(i,j,k)=0.0
          thv0(i,j,k)=th0(i,j,k)
        enddo
      endif

!----------------------------
!  calculate pressure, density, and temperature

      do k=kb,ke
        prs0(i,j,k)=p00*(pi0(i,j,k)**cpdrd)
        rho0(i,j,k)=prs0(i,j,k)/(rd*th0(i,j,k)*pi0(i,j,k)*(1.0+qv0(i,j,k)*reps))
        t0(i,j,k)=th0(i,j,k)*pi0(i,j,k)
      enddo

!----------------------------
!  This reduces errors associated with buoyancy term
!  (seems kind of redundant ... but it works)

      !  This qv0 must match specification of qva array in 
      !  the INIT3D subroutine
      !  (i.e., identical bit-for-bit calculation)

!    IF(imoist.eq.1)THEN
!      do k=kb,ke
!        qv0(i,j,k)=rh0(i,j,k)*rslf(prs0(i,j,k),th0(i,j,k)*pi0(i,j,k))
!      enddo
!    ENDIF


      !  This thv0 must exactly match the manner in which thv
      !  is calculated in the SOLVE subroutine
      !  (i.e., identical bit-for-bit calculation)

      do k=kb,ke
        if(imoist.eq.1)then
          thv0(i,j,k)=(th0(i,j,k)+0.0)*(1.0+reps*max(0.0,qv0(i,j,k)))/(1.0+max(0.0,qv0(i,j,k))+max(0.0,qc0(i,j,k)))
        else
          thv0(i,j,k)=th0(i,j,k)
        endif
      enddo

    enddo
    enddo

      call bcs(rho0)
      call comm_3s_start(rho0,sw31,sw32,se31,se32,   &
                              ss31,ss32,sn31,sn32,reqs_s)
      call comm_3s_end(rho0,sw31,sw32,se31,se32,   &
                            ss31,ss32,sn31,sn32,reqs_s)
      call bcs2(rho0)
!!!      call getcorner(rho0,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call getcorner3(rho0)
    do j=jb,je
    do i=ib,ie
      rho0(i,j,0)=rho0(i,j,1)
      rho0(i,j,nk+1)=rho0(i,j,nk)
    enddo
    enddo

      IF(psolver.eq.5)THEN   ! incompressible option:  set rho0 to a constant

        do k=kb,ke
        do j=jb,je
        do i=ib,ie
          rho0(i,j,k) = 1.0
        enddo
        enddo
        enddo

      ENDIF

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        rr0(i,j,k)=1.0/rho0(i,j,k)
      enddo
      enddo
      enddo

      do k=2,nk
      do j=jb,je
      do i=ib,ie
!!!        rf0(i,j,k)=0.5*(rho0(i,j,k-1)+rho0(i,j,k))
        rf0(i,j,k)=c1(i,j,k)*rho0(i,j,k-1)+c2(i,j,k)*rho0(i,j,k)
      enddo
      enddo
      enddo

      do j=jb,je
      do i=ib,ie
        ! cm1r17, 2nd-order extrapolation:
        rf0(i,j,1) = cgs1*rho0(i,j,1)+cgs2*rho0(i,j,2)+cgs3*rho0(i,j,3)
        rf0(i,j,0)=rf0(i,j,1)
        rho0s(i,j) = rf0(i,j,1)
        ! cm1r17, 2nd-order extrapolation:
        rf0(i,j,nk+1) = cgt1*rho0(i,j,nk)+cgt2*rho0(i,j,nk-1)+cgt3*rho0(i,j,nk-2)
      enddo
      enddo

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        rrf0(i,j,k)=1.0/rf0(i,j,k)
      enddo
      enddo
      enddo

!        i = 1
!        j = 1
!        print *
!        print *,'  nk = ',nk
!        k = 1
!        print *,k,rf0(i,j,k)
!      do k=2,nk+1
!        print *,k,rf0(i,j,k),rf0(i,j,k)-rf0(i,j,k-1)
!      enddo
!        print *
!      stop 11111

!-----------------------------------------------------------------------
!  values at surface:

      ! Get surface p/T/q (for surface models and for CAPE calculation):

      IF( (isnd.eq.7) .and. (.not.terrain_flag) )THEN
        ! this section of code only if no terrain
        do j=jb,je
        do i=ib,ie
          prs0s(i,j) = psurf
          pi0s(i,j) = (psurf/p00)**(rd/cp)
          rth0s(i,j) = thsurf**(-1)
        enddo
        enddo
      ELSE
        do j=jb,je
        do i=ib,ie
          thsurf = cgs1*th0(i,j,1)+cgs2*th0(i,j,2)+cgs3*th0(i,j,3)
           qsurf = cgs1*qv0(i,j,1)+cgs2*qv0(i,j,2)+cgs3*qv0(i,j,3)
          thvsurf = thsurf*(1.0+qsurf*reps)/(1.0+qsurf)
          ! use hydrostatic equation:
          pi0s(i,j) = pi0(i,j,1)+zh(i,j,1)*g/(cp*0.5*(thvsurf+thv0(i,j,1)))
          prs0s(i,j) = p00*( pi0s(i,j)**(cp/rd) )
          rth0s(i,j) = thsurf**(-1)
        enddo
        enddo
      ENDIF

      if( psurf.lt.tsmall ) psurf = prs0s(1,1)

      i = 1
      j = 1
      tsurf = cgs1*t0(i,j,1)+cgs2*t0(i,j,2)+cgs3*t0(i,j,3)
      qsurf = cgs1*qv0(i,j,1)+cgs2*qv0(i,j,2)+cgs3*qv0(i,j,3)


!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc  Start definition of base state wind (iwnd option)  ccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!-----------------------------------------------------------------------
!  Get wind profiles  ...  assume zero wind to start

!--------------------------
! 061012:
! Ignore this section if isnd = 7;  in that case, wind profile has
! already been retrieved from the input_sounding file
!--------------------------

    IF(isnd.ne.7)THEN

      do k=1,nk
      do j=1,nj
      do i=1,ni+1
        u0(i,j,k)= 0.0
      enddo
      enddo
      enddo

      do k=1,nk
      do j=1,nj+1
      do i=1,ni
        v0(i,j,k)= 0.0
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------
!  iwnd = 1
!  RKW-type wind profile
!  reference: Rotunno, Klemp, and Weisman, 1988, JAS, 463-485.

      if(iwnd.eq.1)then

        udep1   =     0.0    ! height of bottom of shear layer (m)
        udep2   =  2500.0    ! height of top of shear layer (m)
        uconst1 =     0.0    ! u at bottom of shear layer
        uconst2 =    10.0    ! u at top of shear layer

        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          zu=0.5*(zh(i-1,j,k)+zh(i,j,k))
          if(zu.lt.udep1)then
            u0(i,j,k)=uconst1
          elseif(zu.gt.udep1 .and. zu.lt.udep2)then
            u0(i,j,k)=(uconst2-uconst1)*(zu-udep1)/(udep2-udep1)+uconst1
          else
            u0(i,j,k)=uconst2
          endif
        enddo
        enddo
        enddo

        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
          v0(i,j,k)=0.0
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iwnd = 2
!  Weisman-Klemp type supercell profile

      elseif(iwnd.eq.2)then

        udep1=2000.0
        udep2=6000.0
        umax1=7.0
        umax2=31.0

        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          zu=0.5*(zh(i-1,j,k)+zh(i,j,k))
          if(zu.le.udep1)THEN
            ANGLE=90.0*(zu/udep1)*(pi/180.0)
            u0(i,j,k)=umax1-umax1*cos(ANGLE)
          elseif(zu.gt.udep1 .and. zu.le.udep2)THEN
            u0(i,j,k)=umax1+(zu-udep1)*(umax2-umax1)/(udep2-udep1)
          ELSE
            u0(i,j,k)=umax2
          ENDIF
        enddo
        enddo
        enddo

        vmax1=umax1

        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
          if(zv.le.udep1)THEN
            ANGLE=90.0*(zv/udep1)*(pi/180.0)
            v0(i,j,k)=vmax1*SIN(ANGLE)
          elseif(zv.gt.udep1 .and. zv.le.udep2)THEN
            v0(i,j,k)=vmax1
          ELSE
            v0(i,j,k)=vmax1
          ENDIF
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iwnd = 3
!  Mulit-cell type profile (?)

      elseif(iwnd.eq.3)then

        udep1=0.0
        udep2=7500.0
        umax1=-40.0/pi
        umax2=40.0/pi+40.0

        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          zu=0.5*(zh(i-1,j,k)+zh(i,j,k))
          if(zu.le.udep2)then
            u0(i,j,k)=umax1+(zu-udep1)*(umax2-umax1)/(udep2-udep1)
          else
            u0(i,j,k)=umax2
          endif
        enddo
        enddo
        enddo

        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
          v0(i,j,k)=40.0/pi
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iwnd = 4
!  Multi-cell
!  reference:  Weisman and Klemp, 1982, MWR, 110, 504-520.

      elseif(iwnd.eq.4)then

        umax1=35.0

        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          u0(i,j,k)=umax1*tanh(0.5*(zh(i-1,j,k)+zh(i,j,k))/3000.0)
        enddo
        enddo
        enddo

        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
          v0(i,j,k)=0.
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iwnd = 5
!  reference:  Dornbrack et al., 2005, Atmos. Sci. Let., 6, 118-122

      elseif(iwnd.eq.5)then

        umax  =   15.0
        udep1 = 4000.0
        udep2 = 6000.0

        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          zu=0.5*(zh(i-1,j,k)+zh(i,j,k))
          if(zu.lt.udep1)then
            u0(i,j,k) = umax
          elseif(zu.lt.udep2)then
            alpha=0.25*pi*(1.0+cos(pi*(zu-udep1)/(udep2-udep1)))
            u0(i,j,k) = umax*sin(alpha)
          else
            u0(i,j,k) = 0.0
          endif
        enddo
        enddo
        enddo

        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          zv=0.5*(zh(i,j-1,k)+zh(i,j,k))
          if(zv.lt.udep1)then
            v0(i,j,k) = 0.0
          elseif(zv.lt.udep2)then
            alpha=0.25*pi*(1.0+cos(pi*(zv-udep1)/(udep2-udep1)))
            v0(i,j,k) = umax*cos(alpha)
          else
            v0(i,j,k) = umax
          endif
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------
!  iwnd = 6
!  constant wind

      ELSEIF(iwnd.eq.6)THEN

        do k=kb,ke
        do j=jb,je
        do i=ib,ie+1
          u0(i,j,k) = 10.0
        enddo
        enddo
        enddo

        do k=kb,ke
        do j=jb,je+1
        do i=ib,ie
          v0(i,j,k) =  0.0
        enddo
        enddo
        enddo

!-----------------------------------------------------------------------

      ENDIF    ! endif for iwnd options

    ENDIF   ! endif for isnd=7 check


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!-----------------------------------------------------------------------
!  subtract off umove and vmove (if applicable)

      if(imove.eq.1)then
        do k=1,nk
        do j=1,nj+1
        do i=1,ni+1
          u0(i,j,k)=u0(i,j,k)-umove
          v0(i,j,k)=v0(i,j,k)-vmove
        enddo
        enddo
        enddo
      else
        umove=0.0
        vmove=0.0
      endif

!-----------------------------------------------------------------------
!  Fill in ghost cells

      call bcu(u0)
      call bcv(v0)

      call comm_3u_start(u0,uw31,uw32,ue31,ue32,   &
                            us31,us32,un31,un32,reqs_u)
      call comm_3u_end(u0,uw31,uw32,ue31,ue32,   &
                          us31,us32,un31,un32,reqs_u)
      call comm_3v_start(v0,vw31,vw32,ve31,ve32,   &
                            vs31,vs32,vn31,vn32,reqs_v)
      call comm_3v_end(v0,vw31,vw32,ve31,ve32,   &
                          vs31,vs32,vn31,vn32,reqs_v)

      do j=jb,je
      do i=ib,ie+1
        u0(i,j,   0) = u0(i,j, 1)
        u0(i,j,nk+1) = u0(i,j,nk)
      enddo
      enddo

      do j=jb,je+1
      do i=ib,ie
        v0(i,j,   0) = v0(i,j, 1)
        v0(i,j,nk+1) = v0(i,j,nk)
      enddo
      enddo

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc  End definition of base state wind (iwnd option)  ccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


!---------------------------------------------------------------------

!  Print out base state

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,101)
      if(dowr) write(outfile,102)
      do k=1,nk
        if(dowr) write(outfile,103) k,prs0(1,1,k),pi0(1,1,k),rho0(1,1,k)
      enddo
      if(dowr) write(outfile,*)

101   format(7x,'k    prs0 (Pa)        pi0        rho0 (kg/m^3)')
102   format(4x,'---------------------------------------------------------------------')
103   format(4x,i4,4x,f9.2,4x,f10.7,4x,f10.7)

!-----

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,104)
      if(dowr) write(outfile,102)
      do k=1,nk
        if(dowr) write(outfile,105) k,th0(1,1,k),thv0(1,1,k),t0(1,1,k),   &
                   -1000*(t0(1,1,k)-t0(1,1,k-1))*rdz*mh(1,1,k),  &
                    1.0e4*g*alog(thv0(1,1,k)/thv0(1,1,k-1))*rdz*mh(1,1,k)
      enddo
      if(dowr) write(outfile,*)

104   format(7x,'k    th0 (K)     thv0 (K)     t0 (K)    l.r. (K/km)     N^2')
105   format(4x,i4,4x,f8.4,4x,f8.4,4x,f8.4,4x,f8.4,4x,f8.4)

!-----
! Ri check ... warn user if Ri < 0.25, but do not stop model
! 061021: added N^2 check ... stop model if N^2 < 0

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Richardson number:  NOTE!!!  Ri should be > 0.25 for most applications!'
      if(dowr) write(outfile,*)
      do k=2,nk
        nm = g*alog(thv0(1,1,k)/thv0(1,1,k-1))*rdz*mf(1,1,k)
        if(nm .lt. -1.0e-12)then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  k,z,N^2:',k,0.5*(zh(1,1,k-1)+zh(1,1,k)),nm
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) ' Warning.  N^2 (Brunt-Vaisala frequency squared) is less than zero!'
          if(dowr) write(outfile,*) ' This is really, really not recommended for the base state.'
          if(dowr) write(outfile,*) ' Stopping model ....'
          if(dowr) write(outfile,*)
          call stopcm1
        endif
        dudz = (u0(1,1,k)-u0(1,1,k-1))*rdz*mf(1,1,k)
        dvdz = (v0(1,1,k)-v0(1,1,k-1))*rdz*mf(1,1,k)
        rinum = nm/(1.0e-12+dudz*dudz+dvdz*dvdz)
        if(rinum.gt.0.25)then
          if(dowr) write(outfile,*) '  k,z,Ri:',k,0.5*(zh(1,1,k-1)+zh(1,1,k)),rinum
        else
          if(dowr) write(outfile,*) '  k,z,Ri:',k,0.5*(zh(1,1,k-1)+zh(1,1,k)),rinum,'<---- NOTE!  Ri < 0.25'
        endif
      enddo
      if(dowr) write(outfile,*)

!-----

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,106)
      if(dowr) write(outfile,102)
      do k=1,nk
        if(rh0(1,1,k).gt.0.999 .or. imoist.eq.0)then
          tlcl=t0(1,1,k)
        else
          tlcl=55.0+(2840./(3.5*alog(t0(1,1,k))-    &
                   alog(0.01*prs0(1,1,k)*qv0(1,1,k)/(0.622+qv0(1,1,k)))-4.805))
        endif
        if(dowr) write(outfile,107) k,                                              &
            th0(1,1,k)*exp((3376./tlcl-2.54)*qv0(1,1,k)*(1.0+0.81*qv0(1,1,k))),  &
            rh0(1,1,k),qv0(1,1,k),qc0(1,1,k),ql0(1,1,k)
      enddo
      if(dowr) write(outfile,*)

106   format(7x,'k  theta-e (K)       rh0          qv0          qc0          ql0')
107   format(4x,i4,4x,f8.4,4x,f9.6,4x,f9.6,4x,f9.6,4x,f9.6)

!------------------------------------------------------------------
!  Get CAPE,CIN,etc:

  IF(imoist.eq.1)THEN

    allocate(  pfoo(nk+1) )
    allocate(  tfoo(nk+1) )
    allocate( qvfoo(nk+1) )

    pfoo(1) = psurf * 0.01
    tfoo(1) = tsurf - 273.15
   qvfoo(1) = qsurf

    do k=1,nk
      pfoo(k+1) = 0.01*prs0(1,1,k)
      tfoo(k+1) = t0(1,1,k) - 273.15
      qvfoo(k+1) = qv0(1,1,k)
    enddo

    if(dowr) write(outfile,*)
    if(dowr) write(outfile,*) '  Thermodynamic properties of base-state sounding:'
  IF(terrain_flag)THEN
    if(dowr) write(outfile,*) '    (for lower-left corner of domain:  i=1,j=1)  '
  ENDIF
    if(dowr) write(outfile,*)
    do n=1,3
      call getcape( n , nk+1 , pfoo , tfoo , qvfoo , cape , cin ,   &
                    zlcl, zlfc, zel , psource , tsource , qvsource )
      if(n.eq.1)then
        if(dowr) write(outfile,*) '    for surface parcel:'
      elseif(n.eq.2)then
        if(dowr) write(outfile,*) '    for most-unstable parcel:'
      elseif(n.eq.3)then
        if(dowr) write(outfile,*) '    for mixed-layer parcel:'
      endif
      if(dowr) write(outfile,116) 0.01*psource,tsource,1000.0*qvsource
      if(dowr) write(outfile,118) zlcl,zlfc,zel
      if(dowr) write(outfile,117) cape,cin
116   format('        source p(mb),T(K),qv(g/kg) = ',3(4x,f6.1))
118   format('        LCL,LFC,EL (m AGL)         = ',3(3x,f7.1))
117   format('        CAPE,CIN (J/kg)            = ',2(4x,f6.1))
      if(dowr) write(outfile,*)
    enddo
    if(dowr) write(outfile,*)

    deallocate(  pfoo )
    deallocate(  tfoo )
    deallocate( qvfoo )

  ENDIF

!------------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,108)
      if(dowr) write(outfile,102)
      do k=1,nk
        if(dowr) write(outfile,109) k,u0(1,1,k),v0(1,1,k)
      enddo
      if(dowr) write(outfile,*)

108   format(7x,'k     u0 (m/s)     v0 (m/s)')
109   format(4x,i4,4x,f9.4,4x,f9.4)

!------------------------------------------------------------------
!  Get an average lapse rate:

    IF(myid.eq.0)THEN
      i=1
      j=1
      k=1
      do while( zh(i,j,k).lt.0.25*maxz )
        k=k+1
      enddo
      k1=k
      k=1
      do while( zh(i,j,k).lt.0.75*maxz )
        k=k+1
      enddo
      k2=k
      lr1=(th0(i,j,k2)-th0(i,j,k1))/(zh(i,j,k2)-zh(i,j,k1))
      ths1=th0(i,j,k1)-lr1*zh(i,j,k1)
!!!      print *,'  i,j   = ',i,j
!!!      print *,'  k1,k2 = ',k1,k2
!!!      print *,'  z1,z2 = ',zh(i,j,k1),zh(i,j,k2)
!!!      print *,'  t1,t2 = ',th0(i,j,k1),th0(i,j,k2)
!!!      print *,'  lr,th = ',lr1,ths1
!!!      print *,'  th    = ',ths1,ths1+lr1*zh(i,j,k1),ths1+lr1*zh(i,j,k2),ths1+lr1*maxz
!!!      print *
    ENDIF
    call MPI_BCAST(lr1   ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    call MPI_BCAST(ths1  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
    th00s = ths1
    thlr = lr1
    IF( psolver.eq.1 .or. psolver.eq.4 .or. psolver.eq.5 )THEN
      ! if no time-splitting, just use isentropic reference:
      ! (so base-state theta advection is exactly zero)
      th00s = 300.0
      thlr = 0.0
    ENDIF
    if(dowr) write(outfile,*) '  th00s = ',th00s
    if(dowr) write(outfile,*) '  thlr  = ',thlr
    if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) '  th00,pi00:'
      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        th00(i,j,k)=th00s+thlr*zh(i,j,k)
        rth00(i,j,k)=1.0/th00(i,j,k)
        IF( abs(thlr).lt.tsmall )THEN
          pi00(i,j,k)=1.0-(g/(cp*th00s))*zh(i,j,k)
        ELSE
          pi00(i,j,k)=1.0-(g/(cp*thlr))*alog(th00(i,j,k)/th00s)
        ENDIF
        ! for psolver=1, there is no pi0 advection on small steps:
        IF( psolver.eq.1 .or. psolver.eq.4 .or. psolver.eq.5 ) pi00(i,j,k)=0.0
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,*) k,th00(i,j,k),pi00(i,j,k)
      enddo
      enddo
      enddo
      if(dowr) write(outfile,*)

!------------------------------------------------------------------

      if(dowr) write(outfile,*) 'Leaving BASE'

      return
      end

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    real function getthe(ttype,p,t,qv,qt)
    implicit none
    include 'constants.incl'
    integer ttype
    real p,t,qv,qt

    ! Assumes air is saturated

    real :: pd,cpm,lhv

  IF(ttype.eq.1)THEN
    pd=p/(1.0+reps*qv)
    cpm=cp+cpl*qt
    lhv=lv1-lv2*t
    getthe=t*((p00/pd)**(rd/cpm))*exp(lhv*qv/(cpm*t))
  ELSEIF(ttype.eq.2)THEN
    getthe=t*( (p00/p)**(0.2854*(1.0-0.28*qv)) )   &
            *exp( ((3376.0/t)-2.54)*qv*(1.0+0.81*qv) )
  ENDIF
    end function getthe

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    subroutine revthe(ttype,the,p,qt,t,qv)
    implicit none
    include 'constants.incl'
    integer ttype
    real the,p,qt,rh,t,qv

    ! Assumes air is saturated
    ! Input:   the,p,qt  (t is first guess value upon input)
    ! Output:  t,qv

    integer n
    real tlast,tinc,pd,cpm,lhv,thx,diff
    real rslf

    n=1
    tlast=t
    tinc=0.0
250 continue
      t=tlast+tinc
      qv=rslf(p,t)
    IF(ttype.eq.1)THEN
      pd=p/(1.0+reps*qv)
      cpm=cp+cpl*qt
      lhv=lv1-lv2*t
      thx=t*((p00/pd)**(rd/cpm))*exp(lhv*qv/(cpm*t))
    ELSEIF(ttype.eq.2)THEN
      thx=t*( (p00/p)**(0.2854*(1.0-0.28*qv)) )   &
           *exp( ((3376.0/t)-2.54)*qv*(1.0+0.81*qv) )
    ENDIF

      diff=the-thx
      if(n.ge.40) print *,n,p,the,thx
      if(abs(diff).gt.0.0001 .and. n.lt.50)then
        n=n+1
        tinc=0.30*diff
        tlast=t
        if(abs(tinc).ge.0.0001) go to 250
      elseif(n.ge.50)then
        print *,'n exceeded 50!'
        print *,'n=',n
        stop 1222
      endif

    return
    end subroutine revthe

