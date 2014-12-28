

      subroutine parcel_driver(dt,xh,uh,ruh,yh,vh,rvh,zh,mh,rmh,mf,  &
                               ua,va,wa,pdata,ploc,packet,reqs_p,    &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,      &
                               nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)

!-----------------------------------------------------------------------
!  This subroutine updates the parcel locations
!-----------------------------------------------------------------------

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(npvals,nparcels) :: pdata
      real, intent(inout), dimension(3,nparcels) :: ploc
      real, intent(inout), dimension(4,nparcels) :: packet
      integer, intent(inout), dimension(rmp) :: reqs_p
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer :: n,np,i,j,k,iflag,jflag,kflag
      real :: uval,vval,wval,rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8

#ifdef MPI
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE)
      integer :: proc,ihave
#endif

#ifdef MPI
      ihave = 0
#endif

!----------------------------------------------------------------------
!  get corner info for MPI runs
!  (may not parallelize correctly if this is not done)

#ifdef MPI
      call getcorneru(ua,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcu2(ua)
      call getcornerv(va,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcv2(va)
      call getcornerw(wa,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcw2(wa)
#endif

      do j=0,nj+1
      do i=0,ni+2
        ua(i,j,0)    = ua(i,j,1)
        ua(i,j,nk+1) = ua(i,j,nk)
      enddo
      enddo

      do j=0,nj+2
      do i=0,ni+1
        va(i,j,0)    = va(i,j,1)
        va(i,j,nk+1) = va(i,j,nk)
      enddo
      enddo

      do j=0,nj+1
      do i=0,ni+1
        wa(i,j,0)    = -wa(i,j,2)
        wa(i,j,nk+2) = -wa(i,j,nk)
      enddo
      enddo

      if(timestats.ge.1) time_parcels=time_parcels+mytime()

!----------------------------------------------------------------------
!  Loop through all parcels:  if you have it, update it's location:

    DO np=1,nparcels

      iflag=0
      jflag=0

      do i=1,ni
        if( abs(xh(i)-pdata(1,np)).le.0.5*dx*ruh(i) ) iflag=i
      enddo

    IF(axisymm.eq.1)THEN
      jflag = 1
    ELSE
      do j=1,nj
        if( abs(yh(j)-pdata(2,np)).le.0.5*dy*rvh(j) ) jflag=j
      enddo
    ENDIF

      IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
          (jflag.ge.1.and.jflag.le.nj) )THEN

        i=iflag
        j=jflag

        do k=1,nk
          if( abs(zh(i,j,k)-pdata(3,np)).le.0.5*dz*rmh(i,j,k) ) kflag=k
        enddo

!----------------------------------------------------------------------
!  Data on u points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(2,np).lt.yh(j) )then
          j=j-1
        endif
        if( pdata(3,np).lt.zh(i,j,k) )then
          k=k-1
        endif

        rx = ( pdata(1,np)-xh(i)+0.5*dx*ruh(i) )*rdx*uh(i)
        ry = ( pdata(2,np)-yh(j) )*rdy*vh(j)
        rz = ( pdata(3,np)-zh(iflag,jflag,k) )*rdz*mh(i,j,k)

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        call tri_interp(ni+1,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,ua,uval)

!----------------------------------------------------------------------
!  Data on v points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(1,np).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(3,np).lt.zh(i,j,k) )then
          k=k-1
        endif

        rx = ( pdata(1,np)-xh(i) )*rdx*uh(i)
        ry = ( pdata(2,np)-yh(j)+0.5*dy*rvh(j) )*rdy*vh(j)
        rz = ( pdata(3,np)-zh(iflag,jflag,k) )*rdz*mh(i,j,k)

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        call tri_interp(ni,nj+1,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,va,vval)

!----------------------------------------------------------------------
!  Data on w points


        i=iflag
        j=jflag
        k=kflag

        if( pdata(1,np).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(2,np).lt.yh(j) )then
          j=j-1
        endif

        rx = ( pdata(1,np)-xh(i) )*rdx*uh(i)
        ry = ( pdata(2,np)-yh(j) )*rdy*vh(j)
        rz = ( pdata(3,np)-(zh(iflag,jflag,k)-0.5*dz*rmh(i,j,k)) )*rdz*mh(i,j,k)

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,wa  ,wval)

!-----------------------------------------------------
!  Step locations of parcels forward in time
!-----------------------------------------------------

        pdata(1,np)=pdata(1,np)+dt*uval
      IF(axisymm.eq.1)THEN
        pdata(2,np)=0.0
      ELSE
        pdata(2,np)=pdata(2,np)+dt*vval
      ENDIF
        pdata(3,np)=pdata(3,np)+dt*wval

!-----------------------------------------------------
!  Account for boundary conditions (if necessary)
!-----------------------------------------------------

        if(pdata(1,np).lt. 0.0)then
          if(wbc.eq.1)then
            pdata(1,np)=pdata(1,np)+(maxx-minx)
          else
            pdata(1,np)=minx
          endif
        endif
        if(pdata(1,np).ge.maxx)then
           if(ebc.eq.1)then
             pdata(1,np)=pdata(1,np)-(maxx-minx)
           else
             pdata(1,np)=maxx
           endif
        endif

        if((pdata(2,np).ge.maxy).and.(axisymm.ne.1))then
          if(nbc.eq.1)then
            pdata(2,np)=pdata(2,np)-(maxy-miny)
          else
            pdata(2,np)=maxy
          endif
        endif
        if((pdata(2,np).lt. 0.0).and.(axisymm.ne.1))then
          if(sbc.eq.1)then
            pdata(2,np)=pdata(2,np)+(maxy-miny)
          else
            pdata(2,np)=miny
          endif
        endif

        pdata(3,np)=max(pdata(3,np),0.0)
        pdata(3,np)=min(pdata(3,np),maxz)

#ifdef MPI
        ! keep track of data that needs to be communicated:
        ihave = ihave + 1
        packet(1,ihave) = pdata(1,np)
        packet(2,ihave) = pdata(2,np)
        packet(3,ihave) = pdata(3,np)
        packet(4,ihave) = np
#endif

      ENDIF

    ENDDO

!----------------------------------------------------------------------
!  communicate data  (for MPI runs)

#ifdef MPI
      IF(myid.eq.0)THEN

        ! send parcel locations to processor 0:
        do i=1,ihave
          np=nint(packet(4,i))
          do n=1,3
            pdata(n,np)=packet(n,i)
          enddo
        enddo

        DO proc=1,(numprocs-1)
          call MPI_RECV(ihave,1,MPI_INTEGER,proc,proc,MPI_COMM_WORLD,status,ierr)
          if(ihave.gt.0) call MPI_RECV(packet,4*ihave,MPI_REAL,proc,1000+proc, &
                                       MPI_COMM_WORLD,status,ierr)
          do i=1,ihave
            np=nint(packet(4,i))
            do n=1,3
              pdata(n,np)=packet(n,i)
            enddo
          enddo
        ENDDO

      ELSE

        call MPI_SEND(ihave,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
        if(ihave.gt.0) call MPI_SEND(packet,4*ihave,MPI_REAL,0,1000+myid, &
                                     MPI_COMM_WORLD,ierr)

      ENDIF

      ! send location of all parcels to ALL processors:
      IF(myid.eq.0)THEN
        DO np=1,nparcels
          ploc(1,np)=pdata(1,np)
          ploc(2,np)=pdata(2,np)
          ploc(3,np)=pdata(3,np)
        ENDDO
      ENDIF

      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ploc,3*nparcels,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      IF(myid.ne.0)THEN
        DO np=1,nparcels
          pdata(1,np)=ploc(1,np)
          pdata(2,np)=ploc(2,np)
          pdata(3,np)=ploc(3,np)
        ENDDO
      ENDIF
#endif

!----------------------------------------------------------------------

      if(timestats.ge.1) time_parcels=time_parcels+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine parcel_interp(dt,xh,uh,ruh,yh,vh,rvh,zh,mh,rmh,mf,  &
                               pi0,thv0,th0,the,b,dpdz,thv,qt,prs,   &
                               ua,va,wa,ppi,th,tha,qa,kh,pdata,      &
                               ploc,packet,reqs_p,                   &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,      &
                               nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      implicit none

!-----------------------------------------------------------------------
!  This subroutine interpolates model information to the parcel locations
!  (diagnostic only ... not used for model integration)
!-----------------------------------------------------------------------

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh,mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: mf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,th0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: the,b,dpdz,thv,qt,prs
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: th
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kh
      real, intent(inout), dimension(npvals,nparcels) :: pdata
      real, intent(inout), dimension(3,nparcels) :: ploc
      real, intent(inout), dimension(npvals+1,nparcels) :: packet
      integer, intent(inout), dimension(rmp) :: reqs_p
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer :: n,np,i,j,k,iflag,jflag,kflag
      real :: tx,cpm,qvs,tem
      real :: uval,vval,wval,rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8
      real :: rslf

#ifdef MPI
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE)
      integer :: proc,ihave
#endif

!----------------------------------------------------------------------

    call bcs(prs)
#ifdef MPI
    call comm_1s_start(prs,pw1,pw2,pe1,pe2,   &
                           ps1,ps2,pn1,pn2,reqs_p)
    call comm_1s_end(prs,pw1,pw2,pe1,pe2,   &
                         ps1,ps2,pn1,pn2,reqs_p)
#endif

    IF(imoist.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,tx,cpm,qvs)
    do k=0,nk+1

      do j=0,nj+1
      do i=0,ni+1
        qt(i,j,k)=qa(i,j,k,nqv)
      enddo
      enddo
      do n=nql1,nql2
        do j=0,nj+1
        do i=0,ni+1
          qt(i,j,k)=qt(i,j,k)+qa(i,j,k,n)
        enddo
        enddo
      enddo
      do n=nqs1,nqs2
        do j=0,nj+1
        do i=0,ni+1
          qt(i,j,k)=qt(i,j,k)+qa(i,j,k,n)
        enddo
        enddo
      enddo
      do j=0,nj+1
      do i=0,ni+1
        th(i,j,k)=(th0(i,j,k)+tha(i,j,k))
        tx=th(i,j,k)*(pi0(i,j,k)+ppi(i,j,k))
        cpm=cp+cpl*qa(i,j,k,nqv)
        do n=nql1,nql2
          cpm=cpm+cpl*qa(i,j,k,n)
        enddo
        qvs=rslf( prs(i,j,k) , tx )
        ! reversible theta-e:
        the(i,j,k)=tx*((p00*(1.0+qa(i,j,k,nqv)*reps)      &
                        /prs(i,j,k))**(rd/cpm))       &
           *((qa(i,j,k,nqv)/qvs)**(-qa(i,j,k,nqv)*rv/cpm))    &
           *exp((lv1-lv2*tx)*qa(i,j,k,nqv)/(cpm*tx))
        thv(i,j,k)=(th0(i,j,k)+tha(i,j,k))*(1.0+reps*qa(i,j,k,nqv))/(1.0+qt(i,j,k))
        b(i,j,k)=g*( thv(i,j,k)/thv0(i,j,k) - 1.0 )
      enddo
      enddo

    enddo

    ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+1
      do j=0,nj+1
      do i=0,ni+1
        the(i,j,k)=th0(i,j,k)+tha(i,j,k)
        thv(i,j,k)=th0(i,j,k)+tha(i,j,k)
        b(i,j,k)=g*( thv(i,j,k)/thv0(i,j,k) - 1.0 )
      enddo
      enddo
      enddo

    ENDIF

    if(timestats.ge.1) time_parcels=time_parcels+mytime()

#ifdef MPI
      ihave = 0
#endif

!----------------------------------------------------------------------
!  get corner info for MPI runs
!  (may not parallelize correctly if this is not done)

#ifdef MPI
      call getcorneru(ua,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcu2(ua)
      call getcornerv(va,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcv2(va)
      call getcornerw(wa,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bcw2(wa)
#endif

      do j=0,nj+1
      do i=0,ni+2
        ua(i,j,0)    = ua(i,j,1)
        ua(i,j,nk+1) = ua(i,j,nk)
      enddo
      enddo

      do j=0,nj+2
      do i=0,ni+1
        va(i,j,0)    = va(i,j,1)
        va(i,j,nk+1) = va(i,j,nk)
      enddo
      enddo

      do j=0,nj+1
      do i=0,ni+1
        wa(i,j,0)    = -wa(i,j,2)
        wa(i,j,nk+2) = -wa(i,j,nk)
      enddo
      enddo

      if(timestats.ge.1) time_parcels=time_parcels+mytime()

      call prepcorners(th,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call prepcorners(prs,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      if(imoist.eq.1)then
        if(nqv.ne.0) call prepcorners(qa(ib,jb,kb,nqv),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        if(nqc.ne.0) call prepcorners(qa(ib,jb,kb,nqc),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        if(nqr.ne.0) call prepcorners(qa(ib,jb,kb,nqr),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      endif
      call prepcorners(the,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call prepcorners(b,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
#ifdef MPI
      call getcornert(kh,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bct2(kh)
#endif

!----------------------------------------------------------------------
!
!  Currently, pdata( 1) = x
!                  ( 2) = y
!                  ( 3) = z
!                  ( 4) = qv
!                  ( 5) = qc  ( = qa(2) for kessler and goddard)
!                  ( 6) = qr  ( = qa(3) for kessler and goddard)
!                       (may have to be changed for other microphysics schemes)
!                  ( 7) = theta
!                  ( 8) = u
!                  ( 9) = v
!                  (10) = w
!                  (11) = kh
!                  (12) = theta-e
!                  (13) = b
!                  (14) = prs
!
!    npvals should be 14 in param.F
!
!----------------------------------------------------------------------
!  Loop through all parcels:  if you have it, get interpolated info:

    DO np=1,nparcels

      iflag=0
      jflag=0

      do i=1,ni
        if( abs(xh(i)-pdata(1,np)).le.0.5*dx*ruh(i) ) iflag=i
      enddo

    IF(axisymm.eq.1)THEN
      jflag = 1
    ELSE
      do j=1,nj
        if( abs(yh(j)-pdata(2,np)).le.0.5*dy*rvh(j) ) jflag=j
      enddo
    ENDIF

      IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
          (jflag.ge.1.and.jflag.le.nj) )THEN

        i=iflag
        j=jflag

        do k=1,nk
          if( abs(zh(i,j,k)-pdata(3,np)).le.0.5*dz*rmh(i,j,k) ) kflag=k
        enddo

!----------------------------------------------------------------------
!  Data on u points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(2,np).lt.yh(j) )then
          j=j-1
        endif
        if( pdata(3,np).lt.zh(i,j,k) )then
          k=k-1
        endif

        rx = ( pdata(1,np)-xh(i)+0.5*dx*ruh(i) )*rdx*uh(i)
        ry = ( pdata(2,np)-yh(j) )*rdy*vh(j)
        rz = ( pdata(3,np)-zh(iflag,jflag,k) )*rdz*mh(i,j,k)

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        call tri_interp(ni+1,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,ua,uval)

!----------------------------------------------------------------------
!  Data on v points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(1,np).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(3,np).lt.zh(i,j,k) )then
          k=k-1
        endif

        rx = ( pdata(1,np)-xh(i) )*rdx*uh(i)
        ry = ( pdata(2,np)-yh(j)+0.5*dy*rvh(j) )*rdy*vh(j)
        rz = ( pdata(3,np)-zh(iflag,jflag,k) )*rdz*mh(i,j,k)

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        call tri_interp(ni,nj+1,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,va,vval)

!----------------------------------------------------------------------
!  Data on w points


        i=iflag
        j=jflag
        k=kflag

        if( pdata(1,np).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(2,np).lt.yh(j) )then
          j=j-1
        endif

        rx = ( pdata(1,np)-xh(i) )*rdx*uh(i)
        ry = ( pdata(2,np)-yh(j) )*rdy*vh(j)
        rz = ( pdata(3,np)-(zh(iflag,jflag,k)-0.5*dz*rmh(i,j,k)) )*rdz*mh(i,j,k)

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,wa  ,wval)
      if(iturb.ge.1)then
        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,kh  ,pdata(11,np))
      endif

!----------------------------------------------------------------------
!  Data on scalar points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(1,np).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(2,np).lt.yh(j) )then
          j=j-1
        endif
        if( pdata(3,np).lt.zh(i,j,k) )then
          k=k-1
        endif

        rx = ( pdata(1,np)-xh(i) )*rdx*uh(i)
        ry = ( pdata(2,np)-yh(j) )*rdy*vh(j)
        rz = ( pdata(3,np)-zh(iflag,jflag,k) )*rdz*mh(i,j,k)

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz


      if(imoist.eq.1)then
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qa(ib,jb,kb,nqv),pdata( 4,np))
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qa(ib,jb,kb,nqc),pdata( 5,np))
        if(ptype.ne.6)then
          call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qa(ib,jb,kb,nqr),pdata( 6,np))
        endif
      endif
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,th ,pdata( 7,np))
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,the,pdata(12,np))
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,b  ,pdata(13,np))
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,prs,pdata(14,np))

!----------------------------------------------------------------------

        pdata( 8,np)=uval
        pdata( 9,np)=vval
        pdata(10,np)=wval

#ifdef MPI
        ihave = ihave + 1
        do n=1,npvals
          packet(n,ihave) = pdata(n,np)
        enddo
        packet(npvals+1,ihave) = np
#endif

      ENDIF

    ENDDO

!----------------------------------------------------------------------
!  communicate data

#ifdef MPI
      IF(myid.eq.0)THEN

        do i=1,ihave
          np=nint(packet(npvals+1,i))
          do n=1,npvals
            pdata(n,np)=packet(n,i)
          enddo
        enddo

        DO proc=1,(numprocs-1)
          call MPI_RECV(ihave,1,MPI_INTEGER,proc,proc,MPI_COMM_WORLD,status,ierr)
          if(ihave.gt.0) call MPI_RECV(packet,(npvals+1)*ihave,MPI_REAL,proc,1000+proc,MPI_COMM_WORLD,status,ierr)
          do i=1,ihave
            np=nint(packet(npvals+1,i))
            do n=1,npvals
              pdata(n,np)=packet(n,i)
            enddo
          enddo
        ENDDO

      ELSE

        call MPI_SEND(ihave,1,MPI_INTEGER,0,myid,MPI_COMM_WORLD,ierr)
        if(ihave.gt.0) call MPI_SEND(packet,(npvals+1)*ihave,MPI_REAL,0,1000+myid,MPI_COMM_WORLD,ierr)

      ENDIF
#endif

!----------------------------------------------------------------------

      if(timestats.ge.1) time_parcels=time_parcels+mytime()

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine parcel_write(prec,rtime,pdata)
      implicit none

      include 'input.incl'

      integer, intent(inout) :: prec
      real, intent(in) :: rtime
      real, intent(in), dimension(npvals,nparcels) :: pdata

      integer :: n,np

!----------------------------------------------------------------------
!  write out data

    IF(myid.eq.0)THEN

      if(dowr) write(outfile,*)
      if(dowr) print *,'  pdata nwrite = ',1+rtime/prclfrq 

      IF(output_format.eq.1)THEN
        ! GrADS format:

        string(totlen+1:totlen+22) = '_pdata.dat            '
        if(dowr) write(outfile,*) string
        open(unit=61,file=string,form='unformatted',access='direct',   &
             recl=4,status='unknown')

          do n=1,npvals
          do np=1,nparcels
            write(61,rec=prec) pdata(n,np)
            prec=prec+1
          enddo
          enddo

        close(unit=61)

#ifdef NETCDF
      ELSEIF(output_format.eq.2)THEN

        call writepdata_nc(prec,rtime,pdata)

#endif
#ifdef HDFOUT
      ELSEIF(output_format.ge.3)THEN

        call writepdata_hdf5(prec,rtime,pdata)

#endif
      ENDIF
      if(dowr) write(outfile,*)

    ENDIF   ! endif for myid=0

      return
      end subroutine parcel_write


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine tri_interp(iz,jz,kz,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,s,pdata)
      implicit none

      include 'input.incl'

      integer :: iz,jz,kz,i,j,k
      real :: w1,w2,w3,w4,w5,w6,w7,w8
      real, dimension(1-ngxy:iz+ngxy,1-ngxy:jz+ngxy,1-ngz:kz+ngz) :: s
      real :: pdata

      pdata=s(i  ,j  ,k  )*w1    &
           +s(i+1,j  ,k  )*w2    &
           +s(i  ,j+1,k  )*w3    &
           +s(i  ,j  ,k+1)*w4    &
           +s(i+1,j  ,k+1)*w5    &
           +s(i  ,j+1,k+1)*w6    &
           +s(i+1,j+1,k  )*w7    &
           +s(i+1,j+1,k+1)*w8

      return
      end


