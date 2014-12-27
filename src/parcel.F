

      subroutine parcel_driver(dt,xh,uh,ruh,xf,yh,vh,rvh,yf,zh,mh,rmh,zf,mf,  &
                               znt,ua,va,wa,pdata,ploc,packet,reqs_p,         &
                               pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,               &
                               nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)

!-----------------------------------------------------------------------
!  This subroutine updates the parcel locations
!-----------------------------------------------------------------------

      include 'input.incl'
      include 'constants.incl'
      include 'timestat.incl'

      real, intent(in) :: dt
      real, intent(in), dimension(ib:ie) :: xh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je) :: znt
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
      real :: uval,vval,wval,rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8,wsum
      real :: x3d,y3d,z3d

#ifdef MPI
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE)
      integer :: proc,ihave
#endif

      logical, parameter :: debug = .false.

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

!----------------------------------------------------------------------
!  apply bottom/top boundary conditions:
!  [Note:  for u,v the array index (i,j,0) means the surface, ie z=0]
!     (for the parcel subroutines only!)

!$omp parallel do default(shared)  &
!$omp private(i,j)
  DO j=jb,je+1

    IF(bbc.eq.1)THEN
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,0) = cgs1*ua(i,j,1)+cgs2*ua(i,j,2)+cgs3*ua(i,j,3)
      enddo
      do i=ib,ie
        zh(i,j,0) = cgs1*zh(i,j,1)+cgs2*zh(i,j,2)+cgs3*zh(i,j,3)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,0) = cgs1*va(i,j,1)+cgs2*va(i,j,2)+cgs3*va(i,j,3)
      enddo
    ELSEIF(bbc.eq.2)THEN
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,0) = 0.0
      enddo
      do i=ib,ie
        zh(i,j,0) = cgs1*zh(i,j,1)+cgs2*zh(i,j,2)+cgs3*zh(i,j,3)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,0) = 0.0
      enddo
    ELSEIF(bbc.eq.3)THEN
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,0) = 0.0
      enddo
      do i=ib,ie
        zh(i,j,0) = znt(i,j)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,0) = 0.0
      enddo
    ENDIF

!----------

    IF(tbc.eq.1)THEN
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,nk+1) = cgt1*ua(i,j,nk)+cgt2*ua(i,j,nk-1)+cgt3*ua(i,j,nk-2)
      enddo
      do i=ib,ie
        zh(i,j,nk+1) = cgt1*zh(i,j,nk)+cgt2*zh(i,j,nk-1)+cgt3*zh(i,j,nk-2)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,nk+1) = cgt1*va(i,j,nk)+cgt2*va(i,j,nk-1)+cgt3*va(i,j,nk-2)
      enddo
    ELSEIF(tbc.eq.2)THEN
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,nk+1) = 0.0
      enddo
      do i=ib,ie
        zh(i,j,nk+1) = cgt1*zh(i,j,nk)+cgt2*zh(i,j,nk-1)+cgt3*zh(i,j,nk-2)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,nk+1) = 0.0
      enddo
    ENDIF

!----------

      ! assuming no terrain:
      IF(j.le.je)THEN
      do i=ib,ie
        wa(i,j,1)    = 0.0
        wa(i,j,nk+1) = 0.0
      enddo
      ENDIF

  ENDDO

      if(timestats.ge.1) time_parcels=time_parcels+mytime()

!----------------------------------------------------------------------
!  Loop through all parcels:  if you have it, update it's location:

    nploop:  DO np=1,nparcels

      x3d = pdata(1,np)
      y3d = pdata(2,np)
      z3d = pdata(3,np)

      iflag=0
      jflag=0
      kflag=0

      do i=1,ni
        if( abs(xh(i)-x3d).le.0.5*dx*ruh(i) ) iflag=i
      enddo

    IF(axisymm.eq.1.or.ny.eq.1)THEN
      jflag = 1
    ELSE
      do j=1,nj
        if( abs(yh(j)-y3d).le.0.5*dy*rvh(j) ) jflag=j
      enddo
    ENDIF

      myparcel:  IF( (iflag.ge.1.and.iflag.le.ni) .and.   &
                     (jflag.ge.1.and.jflag.le.nj) )THEN

        i=iflag
        j=jflag

        kflag = 0
        do while( z3d.ge.zf(iflag,jflag,kflag) )
          kflag = kflag+1
        enddo
        kflag = kflag-1

        IF(debug)THEN
        if( kflag.le.0 .or. kflag.ge.(nk+1) )then
          print *,iflag,jflag,pdata(3,np),zf(iflag,jflag,1)
          print *,myid,pdata(1,np),pdata(2,np),pdata(3,np)
          print *,'  16667 '
          call stopcm1
        endif
        ENDIF

!----------------------------------------------------------------------
!  Data on u points

        i=iflag
        j=jflag
        k=kflag

        if( y3d.lt.yh(j) )then
          j=j-1
        endif
        if( z3d.lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

        rx = ( x3d-xf(i) )/( xf(i+1)-xf(i) )
        ry = ( y3d-yh(j) )/( yh(j+1)-yh(j) )
        rz = ( z3d-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.0 .or. i.gt.(ni+1)   .or.        &
            j.lt.-1 .or. j.gt.(nj+1)   .or.       &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  13333a: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xf1,x3d,xf2 = ',xf(i),x3d,xf(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),y3d,yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),z3d,zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni+1,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,ua,uval)

!----------------------------------------------------------------------
!  Data on v points

        i=iflag
        j=jflag
        k=kflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( z3d.lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

        rx = ( x3d-xh(i) )/( xh(i+1)-xh(i) )
        ry = ( y3d-yf(j) )/( yf(j+1)-yf(j) )
        rz = ( z3d-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.(ni+1)   .or.       &
            j.lt.0 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  23333b: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),x3d,xh(i+1)
          print *,'  yf1,y3d,yh2 = ',yf(j),y3d,yf(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),z3d,zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni,nj+1,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,va,vval)

!----------------------------------------------------------------------
!  Data on w points

        i=iflag
        j=jflag
        k=kflag

        if( x3d.lt.xh(i) )then
          i=i-1
        endif
        if( y3d.lt.yh(j) )then
          j=j-1
        endif

        rx = ( x3d-xh(i) )/( xh(i+1)-xh(i) )
        ry = ( y3d-yh(j) )/( yh(j+1)-yh(j) )
        rz = ( z3d-zf(iflag,jflag,k) )/( zf(iflag,jflag,k+1)-zf(iflag,jflag,k) )

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.(ni+1)   .or.       &
            j.lt.-1 .or. j.gt.(nj+1)   .or.       &
            k.lt.1 .or. k.gt.nk                   )then
          print *
          print *,'  43333a: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),x3d,xh(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),y3d,yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zf(iflag,jflag,k),z3d,zf(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,wa,wval)

!-----------------------------------------------------
!  Update parcel positions:
!-----------------------------------------------------

        IF(nx.eq.1)THEN
          x3d=0.0
        ELSE
          x3d=pdata(1,np)+dt*uval
        ENDIF
        IF(axisymm.eq.1.or.ny.eq.1)THEN
          y3d=0.0
        ELSE
          y3d=pdata(2,np)+dt*vval
        ENDIF
        z3d=pdata(3,np)+dt*wval

!-----------------------------------------------------
!  Account for boundary conditions (if necessary)
!-----------------------------------------------------

        ! New for cm1r17:  if parcel exits domain,
        ! just assume periodic lateral boundary conditions
        ! (no matter what actual settings are for wbc,ebc,sbc,nbc)

        if(x3d.lt.minx)then
          x3d=x3d+(maxx-minx)
        endif
        if(x3d.gt.maxx)then
          x3d=x3d-(maxx-minx)
        endif

        if( (y3d.gt.maxy).and.(axisymm.ne.1).and.(ny.ne.1) )then
          y3d=y3d-(maxy-miny)
        endif
        if( (y3d.lt.miny).and.(axisymm.ne.1).and.(ny.ne.1) )then
          y3d=y3d+(maxy-miny)
        endif

        IF( z3d.lt.0.0 )THEN
          print *,'  parcel is below surface:  np,x3d,y3d,z3d = ',np,x3d,y3d,z3d
          z3d=1.0e-6
        ENDIF
        z3d=min(z3d,maxz)

        pdata(1,np)=x3d
        pdata(2,np)=y3d
        pdata(3,np)=z3d

#ifdef MPI
        ! keep track of data that needs to be communicated:
        ihave = ihave + 1
        packet(1,ihave) = pdata(1,np)
        packet(2,ihave) = pdata(2,np)
        packet(3,ihave) = pdata(3,np)
        packet(4,ihave) = np
#endif

      ENDIF  myparcel

    ENDDO  nploop

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


      subroutine parcel_interp(dt,xh,uh,ruh,xf,yh,vh,rvh,yf,zh,mh,rmh,zf,mf,  &
                               znt,pi0,thv0,th0,the,b,dpdz,thv,qt,prs,   &
                               ua,va,wa,ppi,th,tha,qa,kh,pta,pdata,  &
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
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh,rmh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf,mf
      real, intent(in), dimension(ib:ie,jb:je) :: znt
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,thv0,th0
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: the,b,dpdz,thv,qt,prs
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: th
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(in), dimension(ibc:iec,jbc:jec,kbc:kec) :: kh
      real, intent(in), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, intent(inout), dimension(npvals,nparcels) :: pdata
      real, intent(inout), dimension(3,nparcels) :: ploc
      real, intent(inout), dimension(npvals+1,nparcels) :: packet
      integer, intent(inout), dimension(rmp) :: reqs_p
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2

      integer :: n,np,i,j,k,iflag,jflag,kflag
      real :: tx,qvs,tem
      real :: uval,vval,wval,rx,ry,rz,w1,w2,w3,w4,w5,w6,w7,w8,wsum
      real :: rslf

#ifdef MPI
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE)
      integer :: proc,ihave
#endif

      logical, parameter :: debug = .false.

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
!$omp private(i,j,k,n,tx,qvs)
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
      IF(iice.eq.1)THEN
        do n=nqs1,nqs2
        do j=0,nj+1
        do i=0,ni+1
          qt(i,j,k)=qt(i,j,k)+qa(i,j,k,n)
        enddo
        enddo
        enddo
      ENDIF
      do j=0,nj+1
      do i=0,ni+1
        th(i,j,k)=(th0(i,j,k)+tha(i,j,k))
        tx=th(i,j,k)*(pi0(i,j,k)+ppi(i,j,k))
        qvs=rslf( prs(i,j,k) , tx )
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
        th(i,j,k)=th0(i,j,k)+tha(i,j,k)
        thv(i,j,k)=th(i,j,k)
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
    if(prkh.ge.1)then
      call getcornert(kh,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
      call bct2(kh)
    endif
#endif

!----------------------------------------------------------------------
!  apply bottom/top boundary conditions:
!  [Note:  for u,v,s the array index (i,j,0) means the surface, ie z=0]
!     (for the parcel subroutines only!)

!$omp parallel do default(shared)  &
!$omp private(i,j)
  DO j=jb,je+1

    IF(bbc.eq.1)THEN
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,0) = cgs1*ua(i,j,1)+cgs2*ua(i,j,2)+cgs3*ua(i,j,3)
      enddo
      do i=ib,ie
        zh(i,j,0) = cgs1*zh(i,j,1)+cgs2*zh(i,j,2)+cgs3*zh(i,j,3)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,0) = cgs1*va(i,j,1)+cgs2*va(i,j,2)+cgs3*va(i,j,3)
      enddo
    ELSEIF(bbc.eq.2)THEN
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,0) = 0.0
      enddo
      do i=ib,ie
        zh(i,j,0) = cgs1*zh(i,j,1)+cgs2*zh(i,j,2)+cgs3*zh(i,j,3)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,0) = 0.0
      enddo
    ELSEIF(bbc.eq.3)THEN
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,0) = 0.0
      enddo
      do i=ib,ie
        zh(i,j,0) = znt(i,j)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,0) = 0.0
      enddo
    ENDIF

!----------

    IF(tbc.eq.1)THEN
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,nk+1) = cgt1*ua(i,j,nk)+cgt2*ua(i,j,nk-1)+cgt3*ua(i,j,nk-2)
      enddo
      do i=ib,ie
        zh(i,j,nk+1) = cgt1*zh(i,j,nk)+cgt2*zh(i,j,nk-1)+cgt3*zh(i,j,nk-2)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,nk+1) = cgt1*va(i,j,nk)+cgt2*va(i,j,nk-1)+cgt3*va(i,j,nk-2)
      enddo
    ELSEIF(tbc.eq.2)THEN
      IF(j.le.je)THEN
      do i=ib,ie+1
        ua(i,j,nk+1) = 0.0
      enddo
      do i=ib,ie
        zh(i,j,nk+1) = cgt1*zh(i,j,nk)+cgt2*zh(i,j,nk-1)+cgt3*zh(i,j,nk-2)
      enddo
      ENDIF
      do i=ib,ie
        va(i,j,nk+1) = 0.0
      enddo
    ENDIF

!----------

      ! assuming no terrain:
      IF(j.le.je)THEN
      do i=ib,ie
        wa(i,j,1)    = 0.0
        wa(i,j,nk+1) = 0.0
      enddo
      ENDIF

  ENDDO

      if(timestats.ge.1) time_parcels=time_parcels+mytime()

      if(imoist.eq.1)then
        if(prqv.ge.1) call prepcorners(qa(ib,jb,kb,nqv),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        if(prqc.ge.1) call prepcorners(qa(ib,jb,kb,nqc),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        if(prqr.ge.1) call prepcorners(qa(ib,jb,kb,nqr),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        if(prqi.ge.1) call prepcorners(qa(ib,jb,kb,nqi),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        if(prqs.ge.1) call prepcorners(qa(ib,jb,kb,nqs),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        if(prqg.ge.1) call prepcorners(qa(ib,jb,kb,nqg),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      endif
      call prepcorners( th,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call prepcorners(prs,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      call prepcorners(  b,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      if(iptra.eq.1)then
        do n=1,npt
          call prepcorners(pta(ib,jb,kb,n),nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        enddo
      endif

!----------------------------------------------------------------------
!
!      Note:  pdata( 1) = x
!                  ( 2) = y
!                  ( 3) = z
!
!----------------------------------------------------------------------
!  Loop through all parcels:  if you have it, get interpolated info:

    DO np=1,nparcels

      iflag=0
      jflag=0

      do i=1,ni
        if( abs(xh(i)-pdata(1,np)).le.0.5*dx*ruh(i) ) iflag=i
      enddo

    IF(axisymm.eq.1.or.ny.eq.1)THEN
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

        kflag = 0
        do while( pdata(3,np).ge.zf(iflag,jflag,kflag) )
          kflag = kflag+1
        enddo
        kflag = kflag-1

!----------------------------------------------------------------------
!  Data on u points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(2,np).lt.yh(j) )then
          j=j-1
        endif
        if( pdata(3,np).lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

        rx = ( pdata(1,np)-xf(i) )/( xf(i+1)-xf(i) )
        ry = ( pdata(2,np)-yh(j) )/( yh(j+1)-yh(j) )
        rz = ( pdata(3,np)-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.0 .or. i.gt.(ni+1)   .or.        &
            j.lt.-1 .or. j.gt.(nj+1)   .or.       &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  13333b: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xf1,x3d,xf2 = ',xf(i),pdata(1,np),xf(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),pdata(2,np),yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),pdata(3,np),zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni+1,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,ua,uval)

!----------------------------------------------------------------------
!  Data on v points

        i=iflag
        j=jflag
        k=kflag

        if( pdata(1,np).lt.xh(i) )then
          i=i-1
        endif
        if( pdata(3,np).lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

        rx = ( pdata(1,np)-xh(i) )/( xh(i+1)-xh(i) )
        ry = ( pdata(2,np)-yf(j) )/( yf(j+1)-yf(j) )
        rz = ( pdata(3,np)-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.(ni+1)   .or.       &
            j.lt.0 .or. j.gt.(nj+1)   .or.        &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  23333a: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),pdata(1,np),xh(i+1)
          print *,'  yf1,y3d,yh2 = ',yf(j),pdata(2,np),yf(j+1)
          print *,'  zh1,z3d,zh2 = ',zh(iflag,jflag,k),pdata(3,np),zh(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

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

        rx = ( pdata(1,np)-xh(i) )/( xh(i+1)-xh(i) )
        ry = ( pdata(2,np)-yh(j) )/( yh(j+1)-yh(j) )
        rz = ( pdata(3,np)-zf(iflag,jflag,k) )/( zf(iflag,jflag,k+1)-zf(iflag,jflag,k) )

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.ni   .or.           &
            j.lt.-1 .or. j.gt.nj   .or.           &
            k.lt.1 .or. k.gt.nk                   )then
          print *
          print *,'  43333b: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *,'  xh1,x3d,xh2 = ',xh(i),pdata(1,np),xh(i+1)
          print *,'  yh1,y3d,yh2 = ',yh(j),pdata(2,np),yh(j+1)
          print *,'  zh1,z3d,zh2 = ',zf(iflag,jflag,k),pdata(3,np),zf(iflag,jflag,k+1)
          print *
          call stopcm1
        endif
        ENDIF

        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,wa  ,wval)
      if(prkh.ge.1)then
        call tri_interp(ni,nj,nk+1,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,kh  ,pdata(prkh,np))
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
        if( pdata(3,np).lt.zh(iflag,jflag,k) )then
          k=k-1
        endif

        rx = ( pdata(1,np)-xh(i) )/( xh(i+1)-xh(i) )
        ry = ( pdata(2,np)-yh(j) )/( yh(j+1)-yh(j) )
        rz = ( pdata(3,np)-zh(iflag,jflag,k) )/( zh(iflag,jflag,k+1)-zh(iflag,jflag,k) )

        w1=(1.0-rx)*(1.0-ry)*(1.0-rz)
        w2=rx*(1.0-ry)*(1.0-rz)
        w3=(1.0-rx)*ry*(1.0-rz)
        w4=(1.0-rx)*(1.0-ry)*rz
        w5=rx*(1.0-ry)*rz
        w6=(1.0-rx)*ry*rz
        w7=rx*ry*(1.0-rz)
        w8=rx*ry*rz

        IF(debug)THEN
        wsum = w1+w2+w3+w4+w5+w6+w7+w8
        if( rx.lt.-0.0001 .or. rx.gt.1.0001 .or.  &
            ry.lt.-0.0001 .or. ry.gt.1.0001 .or.  &
            rz.lt.-0.0001 .or. rz.gt.1.0001 .or.  &
            wsum.le.0.99999 .or.                  &
            wsum.ge.1.00001 .or.                  &
            i.lt.-1 .or. i.gt.ni   .or.           &
            j.lt.-1 .or. j.gt.nj   .or.           &
            k.lt.0 .or. k.gt.nk                   )then
          print *
          print *,'  15558: '
          print *,'  np          = ',np
          print *,'  myid,i,j,k  = ',myid,i,j,k
          print *,'  rx,ry,rz    = ',rx,ry,rz
          print *,'  wsum        = ',wsum
          print *
          call stopcm1
        endif
        ENDIF

      if(imoist.eq.1)then
        if(prqv.ge.1)  &
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qa(ib,jb,kb,nqv),pdata(prqv,np))
        if(prqc.ge.1)  &
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qa(ib,jb,kb,nqc),pdata(prqc,np))
        if(prqr.ge.1)  &
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qa(ib,jb,kb,nqr),pdata(prqr,np))
        if(prqi.ge.1)  &
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qa(ib,jb,kb,nqi),pdata(prqi,np))
        if(prqs.ge.1)  &
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qa(ib,jb,kb,nqs),pdata(prqs,np))
        if(prqg.ge.1)  &
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,qa(ib,jb,kb,nqg),pdata(prqg,np))
      endif
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,th ,pdata( 7,np))
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,prs,pdata( 8,np))
        call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,b  ,pdata( 9,np))

        if(iptra.eq.1)then
          do n=1,npt
            call tri_interp(ni,nj,nk,i,j,k,w1,w2,w3,w4,w5,w6,w7,w8,pta(ib,jb,kb,n),pdata(9+n,np))
          enddo
        endif

!----------------------------------------------------------------------

        pdata(4,np)=uval
        pdata(5,np)=vval
        pdata(6,np)=wval

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
             recl=4*npvals*nparcels,status='unknown')

        write(61,rec=prec) ((pdata(n,np),np=1,nparcels),n=1,npvals)
        prec=prec+1

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

