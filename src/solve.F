
!-----------------------------------------------------------------------
!  CM1 Numerical Model, Release 16  (cm1r16)
!  6 February 2012
!  http://www.mmm.ucar.edu/people/bryan/cm1/
!-----------------------------------------------------------------------
!  Quick Index:
!    ua/u3d     = velocity in x-direction (m/s)
!    va/v3d     = velocity in y-direction (m/s)
!    wa/w3d     = velocity in z-direction (m/s)
!    tha/th3d   = perturbation potential temperature (K)
!    ppi/pp3d   = perturbation nondimensional pressure ("Exner function")
!    qa/q3d     = mixing ratios of moisture (kg/kg)
!    tkea/tke3d = SUBGRID turbulence kinetic energy (m^2/s^2)
!    kmh/kmv    = turbulent diffusion coefficients for momentum (m^2/s)
!    khh/khv    = turbulent diffusion coefficients for scalars (m^2/s)
!                 (h = horizontal, v = vertical)
!    prs        = pressure (Pa)
!    rho        = density (kg/m^3)
!
!    th0,pi0,prs0,etc = base-state arrays
!
!    xh         = x (m) at scalar points
!    xf         = x (m) at u points
!    yh         = y (m) at scalar points
!    yf         = y (m) at v points
!    zh         = z (m above sea level) of scalar points (aka, "half levels")
!    zf         = z (m above sea level) of w points (aka, "full levels")
!
!    For the axisymmetric model (axisymm=1), xh and xf are radius (m).
!
!  See "The governing equations for CM1" for more details:
!        http://www.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf
!-----------------------------------------------------------------------
!  Some notes:
!
!  - Upon entering solve, the arrays ending in "a" (eg, ua,wa,tha,qa,etc)
!    are equivalent to the arrays ending in "3d" (eg, u3d,w3d,th3d,q3d,etc).
!  - The purpose of solve is to update the variables from time "t" to time
!    "t+dt".  Values at time "t+dt" are stored in the "3d" arrays.
!  - The "ghost zones" (boundaries beyond the computational subdomain) are
!    filled out completely (3 rows/columns) for the "3d" arrays.  To save 
!    unnecessary computations, starting with cm1r15 the "ghost zones" of 
!    the "a" arrays are only filled out to 1 row/column.  Hence, if you 
!    need to do calculations that use a large stencil, you must use the 
!    "3d" arrays (not the "a") arrays.
!  - Arrays named "ten" store tendencies.  Those ending "ten1" store
!    pre-RK tendencies that are calculated once and then held fixed during
!    the RK (Runge-Kutta) sub-steps. 
!  - CM1 uses a low-storage three-step Runge-Kutta scheme.  See Wicker
!    and Skamarock (2002, MWR, p 2088) for more information.
!  - CM1 uses a staggered C grid.  Hence, u arrays have one more grid point
!    in the i direction, v arrays have one more grid point in the j 
!    direction, and w arrays have one more grid point in the k direction
!    (compared to scalar arrays).
!  - CM1 assumes the subgrid turbulence parameters (tke,km,kh) are located
!    at the w points. 
!-----------------------------------------------------------------------

      subroutine solve(nstep,nrec,prec,nwrite,nrst,rbufsz,num_soil_layers,ndt, &
                   dt,dtlast,th00s,thlr,mtime,stattim,taptim,rsttim,radtim,adt,acfl, &
                   dodrag,dosfcflx,cloudvar,rhovar,qname,budname,bud,bud2,qbudget,asq,bsq, &
                   xh,rxh,uh,ruh,xf,rxf,uf,ruf,yh,vh,rvh,yf,vf,rvf,  &
                   xfref,yfref,rds,sigma,rdsf,sigmaf,tauh,taus,zh,mh,rmh,c1,c2,tauf,zf,mf,rmf,  &
                   rstat,rho0s,pi0s,prs0s,rth0s,pi0,rho0,prs0,thv0,th0,th00,pi00,qv0,qc0, &
                   ql0,rr0,rf0,rrf0,                                 &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv, &
                   rain,sws,svs,sps,srs,sgs,sus,shs,                 &
                   doimpl,tsk,thflux,qvflux,cdu,cdv,ce,u1,v1,w1,     &
                   radbcw,radbce,radbcs,radbcn,                      &
                   dum1,dum2,dum3,dum4,divx,rho,rr,rf,prs,           &
                   t11,t12,t13,t22,t23,t33,                          &
                   u0,rru,ua,u3d,uten,uten1,                         &
                   v0,rrv,va,v3d,vten,vten1,                         &
                   rrw,wa,w3d,wten,wten1,ppi,pp3d,ppten,sten,        &
                   tha,th3d,thten,thten1,thterm,                     &
                   qpten,qtten,qvten,qcten,qiten,qa,q3d,qten,zvdarray, &
                   kmh,kmv,khh,khv,tkea,tke3d,tketen,                &
                   dissten,thpten,qvpten,qcpten,qipten,upten,vpten,  &
                   swten,lwten,o30,radsw,rnflx,radswnet,radlwin,rad2d,   &
                   x,y,z,za,zp,                                      &
                   lu_index,kpbl2d,psfc,u10,v10,hfx,qfx,xland,znt,ust,  &
                   hpbl,wspd,psim,psih,gz1oz0,br,                    &
                   CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                    &
                   MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                 &
                   CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD, &
                   f2d,gsw,glw,chklowq,capg,snowc,dsxy,              &
                   slab_zs,slab_dzs,tslb,tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,       &
                   pta,pt3d,ptten,                                   &
                   pdata,cfb,cfa,cfc,ad1,ad2,pdt,deft,rhs,trans,flag, &
                   reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_tk,reqs_q,reqs_t, &
                   nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                  &
                   ww1,ww2,we1,we2,ws1,ws2,wn1,wn2,                  &
                   pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                  &
                   vw1,vw2,ve1,ve2,vs1,vs2,vn1,vn2,                  &
                   uw31,uw32,ue31,ue32,us31,us32,un31,un32,          &
                   vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,          &
                   ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,          &
                   sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,          &
                   rw31,rw32,re31,re32,rs31,rs32,rn31,rn32,          &
                   qw31,qw32,qe31,qe32,qs31,qs32,qn31,qn32,          &
                   tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2,          &
                   kw1,kw2,ke1,ke2,ks1,ks2,kn1,kn2,                  &
                   tw1,tw2,te1,te2,ts1,ts2,tn1,tn2,ploc,packet,dat1,dat2,dat3)
      use module_mp_thompson
      use module_mp_graupel
      use module_sf_sfclay
      use module_bl_ysu
      use module_sf_slab
      use module_sf_oml
      implicit none

      include 'input.incl'
      include 'constants.incl'
      include 'radcst.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif

!-----------------------------------------------------------------------
! Arrays and variables passed into solve
      integer :: notify
      integer, intent(in) :: nstep
      integer, intent(inout) :: nrec,prec,nwrite,nrst
      integer, intent(in) :: rbufsz,num_soil_layers
      integer, intent(inout) :: ndt
      real, intent(inout) :: dt,dtlast
      real, intent(in) :: th00s,thlr
      real*8, intent(inout) :: mtime
      real*8, intent(inout) :: stattim,taptim,rsttim,radtim,adt,acfl
      logical, intent(in) :: dodrag,dosfcflx
      logical, intent(in), dimension(maxq) :: cloudvar,rhovar
      character*3, intent(in), dimension(maxq) :: qname
      character*6, intent(in), dimension(maxq) :: budname
      real*8, intent(inout), dimension(nk) :: bud
      real*8, intent(inout), dimension(nj) :: bud2
      real*8, intent(inout), dimension(nbudget) :: qbudget
      real*8, intent(inout), dimension(numq) :: asq,bsq
      real, intent(in), dimension(ib:ie) :: xh,rxh,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: xf,rxf,uf,ruf
      real, intent(in), dimension(jb:je) :: yh,vh,rvh
      real, intent(in), dimension(jb:je+1) :: yf,vf,rvf
      real, intent(in), dimension(-2:nx+4) :: xfref
      real, intent(in), dimension(-2:ny+4) :: yfref
      real, intent(in), dimension(kb:ke) :: rds,sigma
      real, intent(in), dimension(kb:ke+1) :: rdsf,sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: tauh,taus,zh,mh,rmh,c1,c2
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: tauf,zf,mf,rmf
      real, intent(inout), dimension(stat_out) :: rstat
      real, intent(in), dimension(ib:ie,jb:je) :: rho0s,pi0s,prs0s,rth0s
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: pi0,rho0,prs0,thv0,th0,th00,pi00,qv0,qc0
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: ql0,rr0,rf0,rrf0
      real, intent(in), dimension(itb:ite,jtb:jte) :: zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      real, intent(inout), dimension(ib:ie,jb:je,nrain) :: rain,sws,svs,sps,srs,sgs,sus,shs
      logical, intent(inout), dimension(ib:ie,jb:je) :: doimpl
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,thflux,qvflux,cdu,cdv,ce,u1,v1,w1
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1,dum2,dum3,dum4,divx,rho,rr,rf,prs
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: t11,t12,t13,t22,t23,t33
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: u0
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: rru,ua,u3d,uten,uten1
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: v0
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: rrv,va,v3d,vten,vten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: rrw,wa,w3d,wten,wten1
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,pp3d,ppten,sten
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: tha,th3d,thten,thten1,thterm
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem) :: qpten,qtten,qvten,qcten,qiten
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa,q3d,qten
      real, intent(inout), dimension(ibzvd:iezvd,jbzvd:jezvd,kbzvd:kezvd,nqzvd) :: zvdarray
      real, intent(inout), dimension(ibc:iec,jbc:jec,kbc:kec) :: kmh,kmv,khh,khv
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea,tke3d,tketen
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dissten
      real, intent(inout), dimension(ibb:ieb,jbb:jeb,kbb:keb) :: thpten,qvpten,qcpten,qipten,upten,vpten
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten
      real, intent(in), dimension(ibr:ier,jbr:jer,kbr:ker) :: o30
      real, intent(inout), dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin
      real, intent(inout), dimension(ni,nj,nrad2d) :: rad2d
      real, intent(inout), dimension(ni+1) :: x
      real, intent(inout), dimension(nj+1) :: y
      real, intent(inout), dimension(nk+3) :: z,za
      real, intent(inout), dimension(ni,nj,nk+3) :: zp
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: psfc,u10,v10,hfx,qfx,xland,znt,ust, &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,       &
                                      CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      f2d,gsw,glw,chklowq,capg,snowc,dsxy
      real, intent(in), dimension(num_soil_layers) :: slab_zs,slab_dzs
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta,pt3d,ptten
      real, intent(inout), dimension(npvals,nparcels) :: pdata
      real, intent(in), dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: cfb
      real, intent(in), dimension(kpb:kpe) :: cfa,cfc,ad1,ad2
      complex, intent(inout), dimension(ipb:ipe,jpb:jpe,kpb:kpe) :: pdt,deft
      complex, intent(inout), dimension(ipb:ipe,jpb:jpe) :: rhs,trans
      logical, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: flag
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_tk
      integer, intent(inout), dimension(rmp,numq) :: reqs_q
      integer, intent(inout), dimension(rmp,npt) :: reqs_t
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(jmp,kmp-1) :: ww1,ww2,we1,we2
      real, intent(inout), dimension(imp,kmp-1) :: ws1,ws2,wn1,wn2
      real, intent(inout), dimension(jmp,kmp) :: pw1,pw2,pe1,pe2
      real, intent(inout), dimension(imp,kmp) :: ps1,ps2,pn1,pn2
      real, intent(inout), dimension(jmp,kmp) :: vw1,vw2,ve1,ve2
      real, intent(inout), dimension(imp,kmp) :: vs1,vs2,vn1,vn2
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
      real, intent(inout), dimension(cmp,jmp,kmp,numq) :: qw31,qw32,qe31,qe32
      real, intent(inout), dimension(imp,cmp,kmp,numq) :: qs31,qs32,qn31,qn32
      real, intent(inout), dimension(cmp,jmp,kmt)   :: tkw1,tkw2,tke1,tke2
      real, intent(inout), dimension(imp,cmp,kmt)   :: tks1,tks2,tkn1,tkn2
      real, intent(inout), dimension(jmp,kmt,4)     :: kw1,kw2,ke1,ke2
      real, intent(inout), dimension(imp,kmt,4)     :: ks1,ks2,kn1,kn2
      real, intent(inout), dimension(cmp,jmp,kmp,npt) :: tw1,tw2,te1,te2
      real, intent(inout), dimension(imp,cmp,kmp,npt) :: ts1,ts2,tn1,tn2
      real, intent(inout), dimension(3,nparcels) :: ploc
      real, intent(inout), dimension(npvals+1,nparcels) :: packet
      real, intent(inout), dimension(ni+1,nj+1) :: dat1
      real, intent(inout), dimension(nx+1,ny+1) :: dat2
      real, intent(inout), dimension(ni+1,nj+1,numprocs) :: dat3

!-----------------------------------------------------------------------
! Arrays and variables defined inside solve

      integer i,j,k,n,nrk,bflag,pdef,nn,fnum
      real :: delqv,delpi,delth,delt,fac,weps
      real :: foo1,foo2

      logical :: dorad

      real :: tout,cfl_limit,max_change,dtsm

      real dttmp,rtime,rdt,tem,tem0,thrad,prad,ql
      real :: cpm,cvm
      real*8 afoo,bfoo
      logical :: getdbz

      logical :: doirrad,dosorad
      real :: saltitude,sazimuth,zen
      real, dimension(2) :: x1
      real, dimension(2) :: y1
      real, dimension(rbufsz) :: radbuf
      real, dimension(nkr) :: swtmp,lwtmp
      real, dimension(nkr) :: tem1,tem2,tem3,tem4,tem5,   &
                              tem6,tem7,tem8,tem9,tem10,   &
                              tem11,tem12,tem13,tem14,tem15,   &
                              tem16,tem17
      real, dimension(nkr) :: ptprt,pprt,qv,qc,qr,qi,qs,qh,   &
                              ptbar,pbar,appi,rhostr,zpp,o31

      real :: rad2deg,albedo,albedoz,tema,temb,frac_snowcover
      real :: dtsfc0,dtsfc

      logical :: flag_qi
      integer :: isfflx
      real :: ep1,ep2,rovg
      real :: SVP1,SVP2,SVP3,SVPT0,p1000mb,eomeg,stbolt
      integer :: ifsnow
      real :: dtmin,dfoo,tmp
      real, dimension(ib:ie,jb:je,kb:ke) :: thsave,ppsave
      real :: r1,r2

#ifdef MPI
      ! for mpi comm of cda:
      integer, dimension(8) :: reqs
      real, dimension(3,nj) :: west,newwest,east,neweast
      real, dimension(ni,3) :: south,newsouth,north,newnorth
#endif

!--------------------------------------------------------------------

#ifdef MPI
      nf=0
      nu=0
      nv=0
      nw=0
#endif

      afoo=0.0d0
      bfoo=0.0d0

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Adaptive timestepping:
!   (assumes cflmax has already been calculated)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF( (adapt_dt.eq.1) .and. (myid.eq.0) )THEN
        ! only processor 0 does this:

        cfl_limit  = 1.00    ! maximum CFL allowed  (actually a "target" value)
        max_change = 0.25    ! maximum (percentage) change in timestep

        IF( cflmax.gt.cfl_limit )THEN
          ! decrease timestep timestep:
          dttmp = max( 1.0-max_change , cfl_limit/cflmax )*dt
        ELSE
          ! increase timestep:
          dttmp = min( 1.0+max_change , cfl_limit/cflmax )*dt
        ENDIF

        ! don't allow dt to exceed twice initial timestep
        dttmp = min( dttmp , 2.0*dtl )

        ! ramp-down timestep when approaching output time
      IF( taptim.gt.0.0 )THEN
        tout = sngl( taptim - mtime )
        if( tout.gt.(2.0*dttmp) .and. tout.le.(3.0*dttmp)  )then
          dttmp = 0.33333333*tout
        elseif( tout.gt.dttmp .and. tout.le.(2.0*dttmp)  )then
          dttmp = 0.5*tout
        elseif( tout.le.dttmp )then
          dttmp = tout
        endif
      ENDIF

      IF( rsttim.gt.0.0 )THEN
        ! ramp-down timestep when approaching restart time
        tout = sngl( rsttim - mtime )
        if( tout.gt.(2.0*dttmp) .and. tout.le.(3.0*dttmp)  )then
          dttmp = 0.33333333*tout
        elseif( tout.gt.dttmp .and. tout.le.(2.0*dttmp)  )then
          dttmp = 0.5*tout
        elseif( tout.le.dttmp )then
          dttmp = tout
        endif
      ENDIF

      IF( stattim.gt.0.0 )THEN
        ! ramp-down timestep when approaching stat time
        tout = sngl( stattim - mtime )
        if( tout.gt.(2.0*dttmp) .and. tout.le.(3.0*dttmp)  )then
          dttmp = 0.33333333*tout
        elseif( tout.gt.dttmp .and. tout.le.(2.0*dttmp)  )then
          dttmp = 0.5*tout
        elseif( tout.le.dttmp )then
          dttmp = tout
        endif
      ENDIF

        dt = dttmp

        ! Algorithm to determine number of small steps:
        IF( psolver.eq.2 )THEN
          ! check dx,dy,dz:
          IF( ny.eq.1 )THEN
            ! 2D sims (x-z):
            dtsm = 0.50*min( min_dx , min_dz )/350.0
          ELSEIF( nx.eq.1 )THEN
            ! 2D sims (y-z):
            dtsm = 0.50*min( min_dy , min_dz )/350.0
          ELSE
            ! 3D sims:
            dtsm = 0.50*min( min_dx , min_dy , min_dz )/350.0
          ENDIF
        ELSEIF( psolver.eq.3 )THEN
          ! check dx,dy:
          IF( ny.eq.1 )THEN
            ! 2D sims (x-z):
            dtsm = 0.60*min_dx/350.0
          ELSEIF( nx.eq.1 )THEN
            ! 2D sims (y-z):
            dtsm = 0.60*min_dy/350.0
          ELSE
            ! 3D sims:
            dtsm = 0.60*min( min_dx , min_dy )/350.0
          ENDIF
        ENDIF
        nsound = max( nint( dt/dtsm ) , 4 )
        if( mod(nsound,2).ne.0 ) nsound = nsound + 1
        if( dt/float(nsound).gt.dtsm ) nsound = nsound + 2

        if( nsound.gt.24 )then
          nsound = 24
          dt = 24*dtsm
        endif

        print *,'cflmax,dt,nsound:',cflmax,dt,nsound

        ! end of processor 0 stuff
      ENDIF

      IF( adapt_dt.eq.1 )THEN
        ! all processors:
#ifdef MPI
        call MPI_BCAST(dt    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nsound,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif
        ndt = ndt + 1
        adt = adt + dt
        acfl = acfl + cflmax
        if(timestats.ge.1) time_misc=time_misc+mytime()
        IF( dt.ne.dtlast )THEN
          IF( (imoist.eq.1).and.(ptype.eq.2) )then
            call consat2(dt)
            if(timestats.ge.1) time_microphy=time_microphy+mytime()
          ENDIF
          IF( (imoist.eq.1).and.(ptype.eq.4) )then
            call lfoice_init(dt)
            if(timestats.ge.1) time_microphy=time_microphy+mytime()
          ENDIF
          dtlast = dt
        ENDIF
      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   radiation  ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      IF( radopt.eq.1 )THEN

        ! time at end of timestep:
        rtime=sngl(mtime+dt)
        dorad = .false.
        IF( rtime.ge.sngl(radtim) ) dorad = .true.
        dtrad = max( dtrad , dt )

        IF( dorad )THEN
!$omp parallel do default(shared)  &
!$omp private(i)
          do i=1,ni+1
            x(i)=xf(i)
          enddo
!$omp parallel do default(shared)  &
!$omp private(j)
          do j=1,nj+1
            y(j)=yf(j)
          enddo
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk+3
          do j=1,nj
          do i=1,ni
            zp(i,j,k) =   zf(i,j,k-1)
          enddo
          enddo
          enddo
          i = 1
          j = 1
!$omp parallel do default(shared)  &
!$omp private(k)
          do k=1,nk+3
            z(k)=zf(1,1,k-1)
            za(k)=zh(1,1,min(k-1,ke))
          enddo
          rtime=sngl(mtime+dt)
          if(dowr) write(outfile,*) '  Calculating radiation tendency:'
          if(timestats.ge.1) time_rad=time_rad+mytime()
          call bcs(prs)
          CALL zenangl( ni,nj, x,y, zp(1,1,2),    &
                rad2d(1,1,ncosz), rad2d(1,1,ncosss), radsw,              &
                dum1(1,1,1),dum1(1,1,2),dum1(1,1,3),dum1(1,1,4),        &
                dum2(1,1,1),dum2(1,1,2),dum2(1,1,3),dum2(1,1,4),        &
                saltitude,sazimuth,dx,dy,dt,rtime,                     &
                ctrlat,ctrlon,year,month,day,hour,minute,second,jday )
          if(myid.eq.0)then
            print *,'    solar zenith angle  (degrees) = ',   &
                                   acos(rad2d(ni,nj,ncosz))*degdpi
            print *,'    solar azimuth angle (degrees) = ',sazimuth*degdpi
          endif
!-----------------------------------------------------------------------
!
!  Calculate surface albedo which is dependent on solar zenith angle
!  and soil moisture. Set the albedo for different types of solar
!  flux to be same.
!
!    rsirbm   Solar IR surface albedo for beam radiation
!    rsirdf   Solar IR surface albedo for diffuse radiation
!    rsuvbm   Solar UV surface albedo for beam radiation
!    rsuvdf   Solar UV surface albedo for diffuse radiation
!
!-----------------------------------------------------------------------
!
  rad2deg = 180.0/3.141592654

!$omp parallel do default(shared)  &
!$omp private(i,j,albedo,albedoz,frac_snowcover,tema)
  DO j=1,nj
    DO i=1,ni

      ! let's just use MM5/WRF value, instead:
      albedo = albd(i,j)

      ! arps code for albedo:
      ! (not sure I trust this.....)

!      albedoz = 0.01 * ( EXP( 0.003286         & ! zenith dependent albedo
!          * SQRT( ( ACOS(rad2d(i,j,ncosz))*rad2deg ) ** 3 ) ) - 1.0 )
!
!      IF ( soilmodel == 0 ) THEN             ! soil type not defined
!!!!        stop 12321
!        tema = 0
!      ELSE
!        tema = qsoil(i,j,1)/wsat(soiltyp(i,j))
!      END IF
!
!      frac_snowcover = MIN(snowdpth(i,j)/snowdepth_crit, 1.0)
!
!      IF ( tema > 0.5 ) THEN
!        albedo = albedoz + (1.-frac_snowcover)*0.14                     &
!                         + frac_snowcover*snow_albedo
!      ELSE
!        albedo = albedoz + (1.-frac_snowcover)*(0.31 - 0.34 * tema)     &
!                         + frac_snowcover*snow_albedo
!      END IF
!        albedo = albedoz

      rad2d(i,j,nrsirbm) = albedo
      rad2d(i,j,nrsirdf) = albedo
      rad2d(i,j,nrsuvbm) = albedo
      rad2d(i,j,nrsuvdf) = albedo

    END DO
  END DO
          ! big OpenMP parallelization loop:
!$omp parallel do default(shared)  &
!$omp private(i,j,k,ptprt,pprt,qv,qc,qr,qi,qs,qh,appi,o31,                 &
!$omp tem1,tem2,tem3,tem4,tem5,tem6,tem7,tem8,tem9,tem10,        &
!$omp tem11,tem12,tem13,tem14,tem15,tem16,tem17,radbuf,swtmp,lwtmp,   &
!$omp doirrad,dosorad,z,za,zpp,ptbar,pbar,rhostr,x1,y1)
        do j=1,nj
        do i=1,ni
          swtmp = 0.0
          lwtmp = 0.0
          do k=1,nk+2
            ptprt(k) =  tha(i,j,k-1)
             pprt(k) =  prs(i,j,k-1) - prs0(i,j,k-1)
               qv(k) =   qa(i,j,k-1,nqv)
               qc(k) =   qa(i,j,k-1,nqc)
               qr(k) =   qa(i,j,k-1,nqr)
               qi(k) =   qa(i,j,k-1,nqi)
               qs(k) =   qa(i,j,k-1,nqs)
               qh(k) =   qa(i,j,k-1,nqg)
             appi(k) =  pi0(i,j,k-1) + ppi(i,j,k-1)
              o31(k) =  o30(i,j,k-1)
          enddo
          ptprt(1) = ptprt(2)
           pprt(1) =  pprt(2)
          ptprt(nk+2) = ptprt(nk+1)
           pprt(nk+2) =  pprt(nk+1)
          x1(1) = xf(i)
          x1(2) = xf(i+1)
          y1(1) = yf(j)
          y1(2) = yf(j+1)
          do k=1,nk+3
            z(k)=zf(i,j,k-1)
            za(k)=zh(i,j,min(k-1,ke))
            zpp(k) =   zp(i,j,k)
          enddo
          do k=2,nk+2
            ptbar(k) =  th0(i,j,k-1)
             pbar(k) = prs0(i,j,k-1)
           rhostr(k) = rho0(i,j,k-1)
          enddo
            ptbar(1) = rth0s(i,j)**(-1)
             pbar(1) = prs0s(i,j)
           rhostr(1) = rho0s(i,j)
            doirrad = .true.
            dosorad = .true.
          CALL radtrns(nir,njr,nkr, rbufsz, 0,myid,dx,dy,            &
                 ib,ie,jb,je,kb,ke,xh,yh,prs0s(i,j),                  &
                 ptprt,pprt,qv,qc,qr,qi,qs,qh,                          &
                 ptbar,pbar,appi,o31,rhostr, tsk(i,j), zpp ,                                 &
                 radsw(i,j),rnflx(i,j),radswnet(i,j),radlwin(i,j), rad2d(i,j,ncosss),            &
                 rad2d(i,j,nrsirbm),rad2d(i,j,nrsirdf),rad2d(i,j,nrsuvbm),                       &
                 rad2d(i,j,nrsuvdf), rad2d(i,j,ncosz),sazimuth,                                  &
                 rad2d(i,j,nfdirir),rad2d(i,j,nfdifir),rad2d(i,j,nfdirpar),rad2d(i,j,nfdifpar),  &
                 tem1, tem2, tem3, tem4, tem5,                &
                 tem6, tem7, tem8, tem9, tem10,               &
                 tem11,tem12,tem13,tem14,tem15,tem16,  &
                 radbuf(1), tem17,swtmp,lwtmp,doirrad,dosorad)
          do k=1,nk
            swten(i,j,k) = swtmp(k+1)
            lwten(i,j,k) = lwtmp(k+1)
          enddo
        enddo
        enddo
          radtim=radtim+dtrad
        ENDIF
        if(timestats.ge.1) time_rad=time_rad+mytime()

      ENDIF


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   surface  ccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!-------------------------------------------------------------------
!  prepare some arrays for WRF surface/pbl physics:

      ! between here and call to ysu:
      ! DO NOT CHANGE:  dum1,dum2,dum3,dum4,sten,t11,t23

      IF((oceanmodel.eq.2).or.(ipbl.eq.1).or.(sfcmodel.eq.2))THEN

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=0.5*(ua(i,j,k)+ua(i+1,j,k))
          dum2(i,j,k)=0.5*(va(i,j,k)+va(i,j+1,k))
          dum3(i,j,k)=th0(i,j,k)+tha(i,j,k)
          sten(i,j,k)=pi0(i,j,k)+ppi(i,j,k)
          dum4(i,j,k)=dum3(i,j,k)*sten(i,j,k)
        enddo
        enddo
        enddo

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          t11(i,j,k) = dz*rmf(i,j,k)
          t23(i,j,k) = prs(i,j,k-1)+(prs(i,j,k)-prs(i,j,k-1))   &
                                   *( zf(i,j,k)- zh(i,j,k-1))   &
                                   /( zh(i,j,k)- zh(i,j,k-1))
        enddo
        enddo
        enddo

        ! values at surface, top of model:
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          t23(i,j,1) =  prs(i,j,1)-zh(i,j,1)*( prs(i,j,2)- prs(i,j,1))   &
                                            /(  zh(i,j,2)-  zh(i,j,1))
          t23(i,j,nk+1)= prs(i,j,nk)+(zf(i,j,nk+1)-zh(i,j,nk))       &
                                    *( prs(i,j,nk)- prs(i,j,nk-1))   &
                                    /(  zh(i,j,nk)-  zh(i,j,nk-1))
          psfc(i,j) = t23(i,j,1)
        enddo
        enddo

        ep1 = rv/rd - 1.0
        ep2 = rd/rv
        rovg = rd/g

        ! dum1 = u at scalars
        ! dum2 = v at scalars
        ! dum3 = th
        ! dum4 = t
        ! sten = pi
        ! t11 = dz8w
        ! t12 = qvten
        ! t13 = qcten
        ! t22 = qiten
        ! t23 = p3di
        ! t33 = exch_h
        ! divx = uten
        ! thterm = vten

        isfflx = 1
        SVP1=0.6112
        SVP2=17.67
        SVP3=29.65
        SVPT0=273.15
        p1000mb      = 100000.
        EOMEG=7.2921E-5
        STBOLT=5.67051E-8

        IF(radopt.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            gsw(i,j)=radsw(i,j)
            glw(i,j)=radlwin(i,j)
          enddo
          enddo
        ELSE
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            gsw(i,j)=0.0
            glw(i,j)=0.0
          enddo
          enddo
        ENDIF

        if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()
      ENDIF

      IF( sfcmodel.ge.1 .and. ipbl.eq.0 )THEN

        call gethpbl(zh,th0,tha,qa,hpbl)

      ENDIF

!-------------------------------------------------------------------
! surface schemes:

!---------------------------------------------------------------------------
! original CM1 formulation:

      if( (sfcmodel.eq.1).or.(idrag.eq.1) )then
        call getcecd(cdu,cdv,ce,u0,v0,u1,v1,w1,ua,va,zh,u10,v10,wspd,xland,znt,ust,cd,  &
                     ww31(1,1,1),ww32(1,1,1),we31(1,1,1),we32(1,1,1),                   &
                     ws31(1,1,1),ws32(1,1,1),wn31(1,1,1),wn32(1,1,1),reqs_s)
      endif

    IF(sfcmodel.eq.1)THEN

      ! get surface flux
      if(isfcflx.eq.1)then
        call sfcflux(dt,ruh,xf,rvh,pi0s,ce,zh,pi0,thv0,th0,u0,v0,tsk,thflux,qvflux,mavail, &
                     rho,rf,u1,v1,w1,ua,va,ppi,tha,qa(ibm,jbm,kbm,nqv), &
                     qbudget(8),psfc,u10,v10,wspd,znt)
      endif

      call sfcdiags(tsk,thflux,qvflux,cdu,cdv,ce,u1,v1,w1,           &
                    xland,psfc,qsfc,u10,v10,hfx,qfx,cda,znt,gz1oz0,  &
                    psim,psih,br,zol,mol,hpbl,dsxy,th2,t2,q2,        &
                    zs,zh,pi0s,pi0,th0,ppi,tha,rho,rf,qa,ua,va)

!-------------------------------------------------------------------

    ELSEIF(sfcmodel.eq.2)THEN

      ! surface layer:
      ! (needed by sfcmodel=2 and ipbl=1)
      call SFCLAY(dum1,dum2,dum4,qa(ib,jb,kb,nqv),prs,t11,       &
                   CP,G,ROVCP,RD,XLV,PSFC,CHS,CHS2,CQS2,CPMM,    &
                   ZNT,UST,hpbl,MAVAIL,ZOL,MOL,REGIME,PSIM,PSIH, &
                   XLAND,HFX,QFX,LH,TSK,FLHC,FLQC,QGH,QSFC,RMOL, &
                   U10,V10,TH2,T2,Q2,                            &
                   GZ1OZ0,WSPD,BR,ISFFLX,dsxy,                   &
                   SVP1,SVP2,SVP3,SVPT0,EP1,EP2,                 &
                   KARMAN,EOMEG,STBOLT,                          &
                   P1000mb,                                      &
                   1  ,ni+1 , 1  ,nj+1 , 1  ,nk+1 ,                    &
                   ib ,ie , jb ,je , kb ,ke ,                    &
                   1  ,ni , 1  ,nj , 1  ,nk ,                    &
                   ustm,ck,cka,cd,cda,isftcflx,iz0tlnd           )

      ifsnow = 0
      dtmin = dt/60.0

      ! slab scheme (MM5/WRF):
      call SLAB(dum4,qa(ib,jb,kb,nqv),prs,FLHC,FLQC,                      &
                   PSFC,XLAND,TMN,HFX,QFX,LH,TSK,QSFC,CHKLOWQ,  &
                   GSW,GLW,CAPG,THC,SNOWC,EMISS,MAVAIL,         &
                   DT,ROVCP,XLV,DTMIN,IFSNOW,               &
                   SVP1,SVP2,SVP3,SVPT0,EP2,                    &
                   KARMAN,EOMEG,STBOLT,                         &
                   TSLB,slab_ZS,slab_DZS,num_soil_layers, .true. ,       &
                   P1000mb,                                     &
                     1, ni+1,   1, nj+1,   1, nk+1,             &
                    ib, ie,  jb, je,  kb, ke,                   &
                     1, ni,   1, nj,   1, nk                    )

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()
    ENDIF

!-------------------------------------------------------------------
! simple ocean mixed layer model based Pollard, Rhines and Thompson (1973)
!   (from WRF)

    IF(oceanmodel.eq.2)THEN

        CALL oceanml(tml,t0ml,hml,h0ml,huml,hvml,ust,dum1,dum2, &
                     tmoml,f2d,g,oml_gamma,                     &
                     xland,hfx,lh,tsk,gsw,glw,emiss,            &
                     dt,STBOLT,                                 &
                       1, ni+1,   1, nj+1,   1, nk+1,           &
                      ib, ie,  jb, je,  kb, ke,                 &
                       1, ni,   1, nj,   1, nk                  )

      if(timestats.ge.1) time_sfcphys=time_sfcphys+mytime()
    ENDIF

!-------------------------------------------------------------------
!  PBL scheme:

      IF(ipbl.eq.1)THEN

        divx = 0.0
        thterm = 0.0
        thten = 0.0
        t12 = 0.0
        t13 = 0.0
        t22 = 0.0

        if( iice.eq.1 )then
          flag_qi = .true.
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nqv) = qa(i,j,k,nqi)
          enddo
          enddo
          enddo
        else
          flag_qi = .false.
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,nqv) = 0.0
          enddo
          enddo
          enddo
        endif

        IF(output_km.eq.1.or.output_kh.eq.1)THEN
          ! ppten = exch_m
          t33=0.0
          ppten=0.0
        ENDIF

        ! PBL:
        call ysu(dum1,dum2,dum3,dum4,qa(ib,jb,kb,nqv),         &
                  qa(ib,jb,kb,nqc),qten(ib,jb,kb,nqv),prs,t23,sten,  &
                  divx,thterm,thten,                           &
                  t12,t13,t22,flag_qi,                         &
                  cp,g,rovcp,rd,rovg,ep1,ep2,karman,xlv,rv,    &
                  t11 ,psfc,                                   &
!!!                  znu,znw,mut,p_top,                        &
                  znt,ust,hpbl,psim,psih,                      &
                  xland,hfx,qfx,gz1oz0,wspd,br,                &
                  dt,kpbl2d,                                   &
                  t33,ppten,                                   &
                  u10,v10,                                     &
                  1  ,ni+1 , 1  ,nj+1 , 1  ,nk+1 ,             &
                  ib ,ie , jb ,je , kb ,ke ,                   &
                  1  ,ni , 1  ,nj , 1  ,nk ,                   &
                  regime                                       )
        if(timestats.ge.1) time_pbl=time_pbl+mytime()

        call bcs(divx)
#ifdef MPI
        call comm_1s_start(divx,pw1,pw2,pe1,pe2,   &
                                ps1,ps2,pn1,pn2,reqs_s)
#endif
        call bcs(thterm)
#ifdef MPI
        call comm_1s_start(thterm,vw1,vw2,ve1,ve2,   &
                                  vs1,vs2,vn1,vn2,reqs_p)
#endif
        IF(output_km.eq.1.or.output_kh.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            khv(i,j,k) = t33(i,j,k)
            kmv(i,j,k) = ppten(i,j,k)
          enddo
          enddo
          enddo
        ENDIF
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          thpten(i,j,k) = thten(i,j,k)
          qvpten(i,j,k) =   t12(i,j,k)
          qcpten(i,j,k) =   t13(i,j,k)
          qipten(i,j,k) =   t22(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
#ifdef MPI
        call comm_1s_end(divx,pw1,pw2,pe1,pe2,   &
                              ps1,ps2,pn1,pn2,reqs_s)
#endif
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=0,ni+1
           upten(i,j,k) =  divx(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
#ifdef MPI
        call comm_1s_end(thterm,vw1,vw2,ve1,ve2,   &
                                vs1,vs2,vn1,vn2,reqs_p)
#endif
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=0,nj+1
        do i=1,ni
           vpten(i,j,k) =thterm(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()

      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   subgrid turbulence schemes  cccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Misc prep:
!  Also, set surface stresses:

      IF( sfcmodel.ge.2 )THEN
        ! put WRF parameters into CM1 arrays:

        IF( dosfcflx )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj
          do i=1,ni
            thflux(i,j) = hfx(i,j)/(cp*rho(i,j,1))
            qvflux(i,j) = qfx(i,j)/rho(i,j,1)
          enddo
          enddo
        ENDIF
        IF( dodrag )THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=0,nj+1
          do i=0,ni+1
            u1(i,j) = 0.5*(ua(i,j,1)+ua(i+1,j,1))
            v1(i,j) = 0.5*(va(i,j,1)+va(i,j+1,1))
            w1(i,j) = sqrt( u1(i,j)**2 + v1(i,j)**2 )
            ce(i,j) = cka(i,j)
          enddo
          enddo
          call bc2d(cda)
          call bc2d(ust)
#ifdef MPI
          call comm_2d_start(cda,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
          call comm_2dew_end(cda,west,newwest,east,neweast,reqs)
          call comm_2dns_end(cda,south,newsouth,north,newnorth,reqs)
          call comm_2d_start(ust,west,newwest,east,neweast,   &
                                 south,newsouth,north,newnorth,reqs)
          call comm_2dew_end(ust,west,newwest,east,neweast,reqs)
          call comm_2dns_end(ust,south,newsouth,north,newnorth,reqs)
#endif
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,nj+1
          do i=1,ni+1
            cdu(i,j) = 0.5*(cda(i-1,j)+cda(i,j))
            cdv(i,j) = 0.5*(cda(i,j-1)+cda(i,j))
            t13(i,j,1) = ((0.5*(ust(i-1,j)+ust(i,j)))**2)*ua(i,j,1)/max(0.5*(w1(i-1,j)+w1(i,j)),0.1)
            t23(i,j,1) = ((0.5*(ust(i,j-1)+ust(i,j)))**2)*va(i,j,1)/max(0.5*(w1(i,j-1)+w1(i,j)),0.1)
          enddo
          enddo
        ENDIF

      ENDIF

        ! get surface drag
        if(idrag.eq.1)then
          call sfcdrag(cdu,cdv,u0,v0,u1,v1,t13,t23,ua,va,u10,v10,wspd)
        endif

!--------------------------------------------------------------------
!                 For turbulence section only:
!  dum1 = squared Brunt-Vaisala frequency (N_m^2) (nm)
!  dum2 = Vertical deformation terms (S_v^2) (defsq)
!  dum3 = Horizontal deformation terms (S_h^2) (defh)
!
!  Arrays available for temporary storage:
!  dum4,divx,ppten,sten,thterm

      IF(iturb.ge.1)THEN

        ! squared Brunt-Vaisala frequency:
        call calcnm(c1,c2,mf,pi0,thv0,th0,cloudvar,dum1,dum2,dum3,dum4,divx,   &
                    prs,ppi,tha,qa)

        ! deformation:
        call calcdef(dodrag,rds,sigma,rdsf,sigmaf,zs,gz,rgz,gzu,rgzu,gzv,rgzv, &
                     xh,rxh,uh,xf,rxf,uf,vh,vf,mh,c1,c2,mf,dum2,dum3,       &
                     divx,ppten,ua,va,wa,t11,t12,t13,t22,t23,t33,gx,gy)

      ENDIF

!--------------------------------------------------------------------
!  Note:  store dissten on w (full) levels, initially:

      IF(idiss.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dissten(i,j,k)=0.0
        enddo
        enddo
        enddo
      ENDIF

!--------------------------------------------------------------------
!  iturb=1:  tke scheme  (for large eddy simulation)
!    Reference:  Deardorff, 1980, Bound Layer Meteor, p. 495
!                see also Stevens, Moeng, Sullivan, 1999, JAS, p. 3963

      IF(iturb.eq.1)THEN
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          tketen(i,j,k)=0.0
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()

        call turbtke(dt,dodrag,dosfcflx,ruh,rvh,rmh,zf,mf,rmf,th0,thflux,qvflux,rth0s,   &
                     dum1,dum2,dum3,dum4,divx,ppten,sten,thterm,        &
                     kmh,kmv,khh,khv,tkea,tketen,t13,t23,ua,va,dissten, &
                     nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                   &
                     kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),       &
                     ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),       &
                     kw1(1,1,2),kw2(1,1,2),ke1(1,1,2),ke2(1,1,2),       &
                     ks1(1,1,2),ks2(1,1,2),kn1(1,1,2),kn2(1,1,2),       &
                     kw1(1,1,3),kw2(1,1,3),ke1(1,1,3),ke2(1,1,3),       &
                     ks1(1,1,3),ks2(1,1,3),kn1(1,1,3),kn2(1,1,3),       &
                     kw1(1,1,4),kw2(1,1,4),ke1(1,1,4),ke2(1,1,4),       &
                     ks1(1,1,4),ks2(1,1,4),kn1(1,1,4),kn2(1,1,4))

!-------------------------------------------------
!  iturb=2:  Smagorinsky scheme  (for large eddy simulation)
!    Reference:  see, e.g., Stevens, Moeng, Sullivan, 1999, JAS, p. 3963

      ELSEIF(iturb.eq.2)THEN

        call turbsmag(dt,dodrag,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s,  &
                      dum1,dum2,dum3,dum4,divx,sten,                   &
                      kmh,kmv,khh,khv,t13,t23,ua,va,dissten,           &
                      nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                      kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),     &
                      ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),     &
                      kw1(1,1,2),kw2(1,1,2),ke1(1,1,2),ke2(1,1,2),     &
                      ks1(1,1,2),ks2(1,1,2),kn1(1,1,2),kn2(1,1,2))

!-------------------------------------------------
!  iturb=3:  parameterized turbulence  (no explicit turbulence)
!    Reference:  Rotunno and Emanuel, 1987, JAS, p. 542
!                Bryan and Rotunno, 2009, MWR, p. 1770

      ELSEIF(iturb.eq.3)THEN

        IF( l_inf.le.tsmall .and. ipbl.eq.1 )THEN
          ! save kmv,khv:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            sten(i,j,k)=kmv(i,j,k)
            ppten(i,j,k)=khv(i,j,k)
          enddo
          enddo
          enddo
        ENDIF

        call turbparam(nstep,zf,dt,dodrag,dosfcflx,ruh,rvh,rmh,mf,rmf,th0,thflux,qvflux,rth0s, &
                      dum1,dum2,dum3,dum4,kmh,kmv,khh,khv,t13,t23,ua,va,dissten,znt, &
                      nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                      kw1(1,1,1),kw2(1,1,1),ke1(1,1,1),ke2(1,1,1),     &
                      ks1(1,1,1),ks2(1,1,1),kn1(1,1,1),kn2(1,1,1),     &
                      kw1(1,1,2),kw2(1,1,2),ke1(1,1,2),ke2(1,1,2),     &
                      ks1(1,1,2),ks2(1,1,2),kn1(1,1,2),kn2(1,1,2))

        IF( l_inf.le.tsmall .and. ipbl.eq.1 )THEN
          ! restore kmv,khv:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            kmv(i,j,k)=sten(i,j,k)
            khv(i,j,k)=ppten(i,j,k)
          enddo
          enddo
          enddo
!!!          ! diagnostic:  effective khv
!!!          !   dum1 is theta flux:
!!!          do j=1,nj
!!!          do i=1,ni
!!!            dum1(i,j,nk+1) = 0.0
!!!            do k=nk,1,-1
!!!              dum1(i,j,k) = dum1(i,j,k+1)+thpten(i,j,k)*rho0(i,j,k)*dz*rmh(i,j,k)
!!!            enddo
!!!            do k=1,nk+1
!!!              khv(i,j,k) = -dum1(i,j,k)/((th0(i,j,k)+tha(i,j,k))-(th0(i,j,k-1)+tha(i,j,k-1)))*rdz*mf(i,j,k)*rf0(i,j,k)
!!!            enddo
!!!          enddo
!!!          enddo
        ENDIF

!-------------------------------------------------

      ELSEIF(iturb.ne.0)THEN

        print *,'  unknown turbulence setting ... '
        call stopcm1

      ENDIF

!--------------------------------------------------------------------
!  Note:  average dissten to scalar (half) levels:

      IF(idiss.eq.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do k=1,nk
        do i=1,ni
          dissten(i,j,k)=0.5*(dissten(i,j,k)+dissten(i,j,k+1))
        enddo
        enddo
        enddo
      ENDIF

!-------------------------------------------------
!  check for columns that need vertically implicit diffusion:

      IF(iturb.ge.1)THEN

!!!        tem0 = 0.125*dz*dz/dt
        tem0 = 0.200*dz*dz/dt

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=0,nj+1
        do i=0,ni+1
          doimpl(i,j) = .false.
          k = 2
          do while( ( .not. doimpl(i,j) ) .and. (k.le.nk) )
            if( khv(i,j,k) .gt. tem0*rmf(i,j,k)*rmf(i,j,k) )then
              doimpl(i,j) = .true.
            endif
            k = k + 1
          enddo
        enddo
        enddo
        if(timestats.ge.1) time_turb=time_turb+mytime()

      ENDIF

!-------------------------------------------------
!  some more calculations for TKE scheme:

      IF(iturb.eq.1)THEN

        call turbt(dt,xh,rxh,uh,xf,uf,vh,vf,mh,mf,rho,rr,rf,  &
                   rds,sigma,gz,rgz,gzu,rgzu,gzv,rgzv,        &
                   dum1,dum2,dum3,dum4,divx,sten,tkea,tketen,kmh,kmv,doimpl)

        if(idiff.eq.1)then
          if(difforder.eq.2)then
            call diff2w(1,0,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,thterm,tkea,tketen)
          elseif(difforder.eq.6)then
            ! for diff6, use '3d' array
            call diff6w(dt,dum1,dum2,dum3,tke3d,tketen)
          endif
        endif

      ENDIF

!-------------------------------------------------
!  Get turbulent stresses:

      IF(iturb.ge.1)THEN

        call gettau(dodrag,xf,rxf,rho,rf,kmh,kmv,t11,t12,t13,t22,t23,t33,ua)

      ENDIF

!--------------------------------------------------------------------
!  Dissipative heating term for pbl scheme:

      IF(ipbl.eq.1.and.idiss.eq.1)THEN
        ! Dissipative heating from ysu scheme:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          ! assume t13,t23 are zero at top of domain:
          t13(i,j,nk+1) = 0.0
          t23(i,j,nk+1) = 0.0
          do k=nk,1,-1
            t13(i,j,k) = t13(i,j,k+1)-upten(i,j,k)*rho0(i,j,k)*dz*rmh(i,j,k)
            t23(i,j,k) = t23(i,j,k+1)-vpten(i,j,k)*rho0(i,j,k)*dz*rmh(i,j,k)
          enddo
          do k=2,nk
            dum2(i,j,k)=0.5*((ua(i,j,k  )+ua(i+1,j,k  ))  &
                            -(ua(i,j,k-1)+ua(i+1,j,k-1)))*rdz*mf(i,j,k)
            dum3(i,j,k)=0.5*((va(i,j,k  )+va(i,j+1,k  ))  &
                            -(va(i,j,k-1)+va(i,j+1,k-1)))*rdz*mf(i,j,k)
          enddo
          dum2(i,j,1)=2.0*0.5*(ua(i,j,1)+ua(i+1,j,1))*rdz*mf(i,j,1)
          dum3(i,j,1)=2.0*0.5*(va(i,j,1)+va(i,j+1,1))*rdz*mf(i,j,1)
          dum2(i,j,nk+1)=0.0
          dum3(i,j,nk+1)=0.0
          do k=1,nk
            dissten(i,j,k)=dissten(i,j,k)+rr0(i,j,k)*0.5*(  &
                             ( t13(i,j,k  )*dum2(i,j,k  )   &
                              +t13(i,j,k+1)*dum2(i,j,k+1) ) &
                            +( t23(i,j,k  )*dum3(i,j,k  )   &
                              +t23(i,j,k+1)*dum3(i,j,k+1) ) )
          enddo
        enddo
        enddo
        t13 = 0.0
        t23 = 0.0
      ENDIF

!--------------------------------------------------------------------


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   Pre-RK calculations   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


!--------------------------------------------------------------------
!  radbc
 
      if(irbc.eq.1)then

        if(ibw.eq.1 .or. ibe.eq.1) call radbcew(radbcw,radbce,ua)
 
        if(ibs.eq.1 .or. ibn.eq.1) call radbcns(radbcs,radbcn,va)

      endif

!--------------------------------------------------------------------
!  U-equation
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni+1
!!!        uten1(i,j,k)=0.
        uten1(i,j,k)=-rdalpha*0.5*(tauh(i-1,j,k)+tauh(i,j,k))*(ua(i,j,k)-u0(i,j,k))
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if(idiff.ge.1)then
        if(difforder.eq.2)then
          call diff2u(1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,dissten,ua,uten1)
        elseif(difforder.eq.6)then
          ! for diff6, use '3d' array
          call diff6u(dt,u0,dum1,dum2,dum3,u3d,uten1)
        endif
      endif

      if(dns.eq.1)then
        call diff2u(2,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,dissten,ua,uten1)
      endif
 
      if(iturb.ge.1)then
        call turbu(dt,dodrag,xh,ruh,xf,rxf,uf,vh,mh,mf,rmf,rho,rf,  &
                   zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gxu,     &
                   dum1,dum2,dum3,dum4,divx,ua,uten1,wa,t11,t12,t13,t22,kmv,doimpl)
      endif

      if(ipbl.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          uten1(i,j,k) = uten1(i,j,k) + 0.5*( upten(i-1,j,k)+ upten(i,j,k))
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
      endif

!--------------------------------------------------------------------
!  V-equation
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj+1
      do i=1,ni
!!!        vten1(i,j,k)=0.
        vten1(i,j,k)=-rdalpha*0.5*(tauh(i,j-1,k)+tauh(i,j,k))*(va(i,j,k)-v0(i,j,k))
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if(idiff.ge.1)then
        if(difforder.eq.2)then
          call diff2v(1,xh,uh,rxf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,dissten,va,vten1)
        elseif(difforder.eq.6)then
          ! for diff6, use '3d' array
          call diff6v(dt,v0,dum1,dum2,dum3,v3d,vten1)
        endif
      endif

      if(dns.eq.1)then
        call diff2v(2,xh,uh,rxf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,dissten,va,vten1)
      endif
 
      if(iturb.ge.1)then
        call turbv(dt,dodrag,xh,rxh,uh,xf,rvh,vf,mh,mf,rho,rf,   &
                   zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gyv,  &
                   dum1,dum2,dum3,dum4,divx,va,vten1,wa,t12,t22,t23,kmv,doimpl)
      endif

      if(ipbl.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vten1(i,j,k) = vten1(i,j,k) + 0.5*( vpten(i,j-1,k)+ vpten(i,j,k))
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
      endif
 
!--------------------------------------------------------------------
!  W-equation
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk+1
      do j=1,nj
      do i=1,ni
!!!        wten1(i,j,k)=0.0
        wten1(i,j,k)=-rdalpha*tauf(i,j,k)*wa(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if(idiff.ge.1)then
        if(difforder.eq.2)then
          call diff2w(1,1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,dissten,wa,wten1)
        elseif(difforder.eq.6)then
          ! for diff6, use '3d' array
          call diff6w(dt,dum1,dum2,dum3,w3d,wten1)
        endif
      endif

      if(dns.eq.1)then
        call diff2w(2,1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,dum4,dissten,wa,wten1)
      endif
 
      if(iturb.ge.1)then
        call turbw(dt,xh,rxh,uh,xf,vh,mh,mf,rho,rf,gz,rgzu,rgzv,rds,sigma,   &
                   dum1,dum2,dum3,dum4,divx,wa,wten1,t13,t23,t33,kmv,doimpl)
      endif

!--------------------------------------------------------------------
!  THETA-equation
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
!!!        thten1(i,j,k)=0.0
        thten1(i,j,k)=-rdalpha*taus(i,j,k)*tha(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_rdamp=time_rdamp+mytime()

      if(idiff.eq.1)then
        if(difforder.eq.2)then
          call diff2s(1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,tha,thten1)
        elseif(difforder.eq.6)then
          ! for diff6, use '3d' array
          call diff6s(dt,ql0,dum1,dum2,dum3,th3d,thten1)
        endif
      endif

!----- cvm (if needed) -----!

      IF( neweqts.ge.1 .and. (idiss.eq.1.or.rterm.eq.1) )THEN
        ! store cvm in dum1:
        ! store ql  in dum2:
        ! store qi  in dum3:
!$omp parallel do default(shared)  &
!$omp private(i,j,k,n)
        DO k=1,nk
          do j=1,nj
          do i=1,ni
            dum2(i,j,k)=0.0
            dum3(i,j,k)=0.0
          enddo
          enddo
          do n=nql1,nql2
            do j=1,nj
            do i=1,ni
              dum2(i,j,k)=dum2(i,j,k)+qa(i,j,k,n)
            enddo
            enddo
          enddo
          IF(iice.eq.1)THEN
            do n=nqs1,nqs2
              do j=1,nj
              do i=1,ni
                dum3(i,j,k)=dum3(i,j,k)+qa(i,j,k,n)
              enddo
              enddo
            enddo
          ENDIF
          do j=1,nj
          do i=1,ni
            dum1(i,j,k)=cv+cvv*qa(i,j,k,nqv)+cpl*dum2(i,j,k)+cpi*dum3(i,j,k)
          enddo
          enddo
        ENDDO
      ENDIF

!----- store appropriate rho for budget calculations in dum2 -----!

      IF(axisymm.eq.1)THEN
       ! for axisymmetric grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
        enddo
        enddo
        enddo
      ELSE
       ! for Cartesian grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum2(i,j,k) = rho(i,j,k)
        enddo
        enddo
        enddo
      ENDIF

!---- Dissipative heating term:

      IF(idiss.eq.1)THEN
        ! use dissten array to store epsilon
        if(imoist.eq.1.and.neweqts.ge.1)then
          ! moist, new equations:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            thten1(i,j,k)=thten1(i,j,k)   &
                        +dissten(i,j,k)/( cpdcv*dum1(i,j,k)*(pi0(i,j,k)+ppi(i,j,k)) )
          enddo
          enddo
          enddo
        else
          ! traditional equations:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            thten1(i,j,k)=thten1(i,j,k)   &
                        +dissten(i,j,k)/( cp*(pi0(i,j,k)+ppi(i,j,k)) )
          enddo
          enddo
          enddo
        endif
      ENDIF

!---- Rotunno-Emanuel "radiation" term
!---- (currently capped at 2 K/day ... see RE87 p 546)

      IF(rterm.eq.1)THEN
        tem0 = dt*dx*dy*dz
!$omp parallel do default(shared)  &
!$omp private(i,j,k,thrad,prad)
        do k=1,nk
        bud(k)=0.0d0
        do j=1,nj
        do i=1,ni
          ! NOTE:  thrad is a POTENTIAL TEMPERATURE tendency
          thrad = -tha(i,j,k)/(12.0*3600.0)
          if( tha(i,j,k).gt. 1.0 ) thrad = -1.0/(12.0*3600.0)
          if( tha(i,j,k).lt.-1.0 ) thrad =  1.0/(12.0*3600.0)
          thten1(i,j,k)=thten1(i,j,k)+thrad
          ! associated pressure tendency:
          prad = (pi0(i,j,k)+ppi(i,j,k))*rddcv*thrad/(th0(i,j,k)+tha(i,j,k))
          ! budget:
          bud(k) = bud(k) + dum1(i,j,k)*dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*( &
                            thrad*(pi0(i,j,k)+ppi(i,j,k))    &
                           + prad*(th0(i,j,k)+tha(i,j,k)) )
        enddo
        enddo
        enddo
        do k=1,nk
          qbudget(10) = qbudget(10) + tem0*bud(k)
        enddo
      ENDIF
      if(timestats.ge.1) time_misc=time_misc+mytime()

      if( (iturb.ge.1).or.(dns.eq.1) )then
        ! use thten to store total potential temperature:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
          thten(i,j,k)=th0(i,j,k)+tha(i,j,k)
        enddo
        enddo
        enddo
      endif

      if(dns.eq.1)then
        call diff2s(2,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,thten,thten1)
      endif

      IF( radopt.eq.1 )THEN
        ! tendency from radiation scheme:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          thten1(i,j,k)=thten1(i,j,k)+(swten(i,j,k)+lwten(i,j,k))
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_rad=time_rad+mytime()
      ENDIF

      IF( ipbl.eq.1 )THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          thten1(i,j,k) = thten1(i,j,k) + thpten(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_pbl=time_pbl+mytime()
      ENDIF

      if(iturb.ge.1)then
        call turbs(1,dt,dosfcflx,xh,rxh,uh,xf,uf,vh,vf,thflux,   &
                   rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                   dum1,dum2,dum3,dum4,divx,sten,rho,rr,rf,thten,thten1,khh,khv,doimpl)
      endif

!-------------------------------------------------------------------
!  contribution to pressure tendency from potential temperature:

      IF(neweqts.ge.1.or.imoist.eq.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ppten(i,j,k) = thten1(i,j,k)*rddcv   &
                        *(pi0(i,j,k)+ppi(i,j,k))/(th0(i,j,k)+tha(i,j,k))
        enddo
        enddo
        enddo
      ELSE
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          ppten(i,j,k) = 0.0
        enddo
        enddo
        enddo
      ENDIF
      if(timestats.ge.1) time_misc=time_misc+mytime()

!-------------------------------------------------------------------
!  budget calculations:

      if(dosfcflx.and.imoist.eq.1)then
        tem0 = dt*dx*dy*dz
!$omp parallel do default(shared)  &
!$omp private(i,j,k,delpi,delth,delqv,delt,n)
        do j=1,nj
        bud2(j) = 0.0d0
        do i=1,ni
          k = 1
          delth = rf0(i,j,1)*rr0(i,j,1)*rdz*mh(i,j,1)*thflux(i,j)
          delqv = rf0(i,j,1)*rr0(i,j,1)*rdz*mh(i,j,1)*qvflux(i,j)
          delpi = rddcv*(pi0(i,j,1)+ppi(i,j,1))*(           &
                                delqv/(eps+qa(i,j,1,nqv))   &
                               +delth/(th0(i,j,1)+tha(i,j,1))  )
          delt = (pi0(i,j,k)+ppi(i,j,k))*delth   &
                +(th0(i,j,k)+tha(i,j,k))*delpi
          bud2(j) = bud2(j) + dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*(        &
                  cv*delt                                                   &
                + cvv*qa(i,j,k,nqv)*delt                                    &
                + cvv*(pi0(i,j,k)+ppi(i,j,k))*(th0(i,j,k)+tha(i,j,k))*delqv &
                + g*zh(i,j,k)*delqv   )
          do n=nql1,nql2
            bud2(j) = bud2(j) + dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*cpl*qa(i,j,k,n)*delt
          enddo
          if(iice.eq.1)then
            do n=nqs1,nqs2
              bud2(j) = bud2(j) + dum2(i,j,k)*ruh(i)*rvh(j)*rmh(i,j,k)*cpi*qa(i,j,k,n)*delt
            enddo
          endif
        enddo
        enddo
        do j=1,nj
          qbudget(9) = qbudget(9) + tem0*bud2(j)
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()
      endif

!-------------------------------------------------------------------
!  Passive Tracers

      if(iptra.eq.1)then
        do n=1,npt
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ptten(i,j,k,n)=0.0
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()
          if(idiff.eq.1)then
            if(difforder.eq.2)then
              call diff2s(1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,pta(ib,jb,kb,n),ptten(ib,jb,kb,n))
            elseif(difforder.eq.6)then
              ! for diff6, use '3d' array
              call diff6s(dt,ql0,dum1,dum2,dum3,pt3d(ib,jb,kb,n),ptten(ib,jb,kb,n))
            endif
          endif
          if(iturb.ge.1)then
            call turbs(0,dt,dosfcflx,xh,rxh,uh,xf,uf,vh,vf,qvflux,   &
                       rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                       dum1,dum2,dum3,dum4,divx,sten,rho,rr,rf,pta(ib,jb,kb,n),ptten(ib,jb,kb,n),khh,khv,doimpl)
          endif
        enddo
      endif

!-------------------------------------------------------------------
!  Moisture

      if(imoist.eq.1)then
        DO n=1,numq
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            qten(i,j,k,n)=0.0
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()
!---------------------------
          ! qv:
          if(n.eq.nqv)then
            if(idiff.eq.1)then
              if(difforder.eq.2)then
                call diff2s(1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,qa(ib,jb,kb,n),qten(ib,jb,kb,n))
              elseif(difforder.eq.6)then
                ! for diff6, use '3d' array
                call diff6s(dt,qv0,dum1,dum2,dum3,q3d(ib,jb,kb,n),qten(ib,jb,kb,n))
              endif
            endif
            if(iturb.ge.1)then
              call turbs(1,dt,dosfcflx,xh,rxh,uh,xf,uf,vh,vf,qvflux,   &
                         rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                         dum1,dum2,dum3,dum4,divx,sten,rho,rr,rf,qa(ib,jb,kb,n),qten(ib,jb,kb,n),khh,khv,doimpl)
            endif
            if(ipbl.eq.1)then
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
              do k=1,nk
              do j=1,nj
              do i=1,ni
                qten(i,j,k,nqv) = qten(i,j,k,nqv) + qvpten(i,j,k)
              enddo
              enddo
              enddo
              if(timestats.ge.1) time_pbl=time_pbl+mytime()
            endif
            IF(neweqts.ge.1)THEN
              ! contribution to pressure tendency from water vapor:
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
              do k=1,nk
              do j=1,nj
              do i=1,ni
                ppten(i,j,k)=ppten(i,j,k)+qten(i,j,k,n)   &
                            *rddcv*(pi0(i,j,k)+ppi(i,j,k))/(eps+qa(i,j,k,n))
              enddo
              enddo
              enddo
              if(timestats.ge.1) time_misc=time_misc+mytime()
            ENDIF
!---------------------------
          ! not qv:
          else
            if(idiff.eq.1)then
              if(difforder.eq.2)then
                call diff2s(1,rxh,uh,xf,uf,vh,vf,mh,mf,dum1,dum2,dum3,qa(ib,jb,kb,n),qten(ib,jb,kb,n))
              elseif(difforder.eq.6)then
                ! for diff6, use '3d' array
                call diff6s(dt,ql0,dum1,dum2,dum3,q3d(ib,jb,kb,n),qten(ib,jb,kb,n))
              endif
            endif
            if(iturb.ge.1)then
              call turbs(0,dt,dosfcflx,xh,rxh,uh,xf,uf,vh,vf,qvflux,   &
                         rds,sigma,rdsf,sigmaf,mh,mf,gz,rgz,gzu,rgzu,gzv,rgzv,gx,gxu,gy,gyv, &
                         dum1,dum2,dum3,dum4,divx,sten,rho,rr,rf,qa(ib,jb,kb,n),qten(ib,jb,kb,n),khh,khv,doimpl)
            endif
          endif
!---------------------------
        ENDDO
        IF(ipbl.eq.1)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if(nqc.ne.0)   &
            qten(i,j,k,nqc) = qten(i,j,k,nqc) + qcpten(i,j,k)
            if(nqi.ne.0)   &
            qten(i,j,k,nqi) = qten(i,j,k,nqi) + qipten(i,j,k)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_pbl=time_pbl+mytime()
        ENDIF
      endif

!-------------------------------------------------------------------

        IF(imoist.eq.1.and.neweqts.ge.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
            do j=1,nj
            do i=1,ni
              ppten(i,j,k)     = ppten(i,j,k)    +qpten(i,j,k)
              thten1(i,j,k)    = thten1(i,j,k)   +qtten(i,j,k)
               qten(i,j,k,nqv) =  qten(i,j,k,nqv)+qvten(i,j,k)
               qten(i,j,k,nqc) =  qten(i,j,k,nqc)+qcten(i,j,k)
            enddo
            enddo
            IF(ptype.eq.2)THEN
              do j=1,nj
              do i=1,ni
                qten(i,j,k,nqi) =  qten(i,j,k,nqi)+qiten(i,j,k)
              enddo
              enddo
            ENDIF
          enddo
        ENDIF

        IF(imoist.eq.0.or.neweqts.eq.0)THEN
          ! use simple (non-conserving) equations:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            thterm(i,j,k)=0.0
            t22(i,j,k)=rddcv
          enddo
          enddo
          enddo
        ENDIF
 
!-------------------------------------------------------------------
!  Parcel update

      if(iprcl.eq.1)then
        call  parcel_driver(dt,xh,uh,ruh,yh,vh,rvh,zh,mh,rmh,mf,    &
                            ua,va,wa,pdata,ploc,packet(1,1),reqs_p, &
                            pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,        &
                            nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
      endif


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   Begin RK section   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ! time at end of full timestep:
      rtime=sngl(mtime+dt)

!--------------------------------------------------------------------
! RK3 begin

      DO NRK=1,3

        dttmp=dt/float(4-nrk)

!--------------------------------------------------------------------
        IF(nrk.ge.2)THEN
#ifdef MPI
          call comm_3u_end(u3d,uw31,uw32,ue31,ue32,   &
                               us31,us32,un31,un32,reqs_u)
          call comm_3v_end(v3d,vw31,vw32,ve31,ve32,   &
                               vs31,vs32,vn31,vn32,reqs_v)
          call comm_3w_end(w3d,ww31,ww32,we31,we32,   &
                               ws31,ws32,wn31,wn32,reqs_w)
#endif
          if(terrain_flag)then
            call bcwsfc(zh,zf,sigma,sigmaf,dzdx,dzdy,dum1,dum2,u3d,v3d,w3d)
            call bc2d(w3d(ib,jb,1))
          endif
        ENDIF
!--------------------------------------------------------------------

    IF(.not.terrain_flag)THEN
      ! without terrain:

!$omp parallel do default(shared)  &
!$omp private(i,j,k,r1,r2)
      DO j=0,nj+1
        do k=1,nk
        do i=0,ni+2
          rru(i,j,k)=rho0(1,1,k)*u3d(i,j,k)
        enddo
        enddo
        do k=1,nk
        do i=0,ni+1
          rrv(i,j,k)=rho0(1,1,k)*v3d(i,j,k)
        enddo
        enddo
        IF(j.eq.(nj+1))THEN
          do k=1,nk
          do i=0,ni+1
            rrv(i,j+1,k)=rho0(1,1,k)*v3d(i,j+1,k)
          enddo
          enddo
        ENDIF
        do k=2,nk
        do i=0,ni+1
          rrw(i,j,k)=rf0(1,1,k)*w3d(i,j,k)
        enddo
        enddo
      ENDDO

    ELSE
      ! with terrain:

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      DO j=0,nj+1
        do k=1,nk
        do i=0,ni+2
          rru(i,j,k)=0.5*(rho0(i-1,j,k)+rho0(i,j,k))*u3d(i,j,k)*rgzu(i,j)
        enddo
        enddo
        do k=1,nk
        do i=0,ni+1
          rrv(i,j,k)=0.5*(rho0(i,j-1,k)+rho0(i,j,k))*v3d(i,j,k)*rgzv(i,j)
        enddo
        enddo
        IF(j.eq.(nj+1))THEN
          do k=1,nk
          do i=0,ni+1
            rrv(i,j+1,k)=0.5*(rho0(i,j,k)+rho0(i,j+1,k))*v3d(i,j+1,k)*rgzv(i,j+1)
          enddo
          enddo
        ENDIF
      ENDDO

!$omp parallel do default(shared)  &
!$omp private(i,j,k,r1,r2)
      DO j=0,nj+1
        do k=2,nk
        r2 = (sigmaf(k)-sigma(k-1))*rds(k)
        r1 = 1.0-r2
        do i=0,ni+1
          rrw(i,j,k)=rf0(i,j,k)*w3d(i,j,k)                                  &
                    +0.5*( ( r2*(rru(i,j,k  )+rru(i+1,j,k  ))               &
                            +r1*(rru(i,j,k-1)+rru(i+1,j,k-1)) )*dzdx(i,j)   &
                          +( r2*(rrv(i,j,k  )+rrv(i,j+1,k  ))               &
                            +r1*(rrv(i,j,k-1)+rrv(i,j+1,k-1)) )*dzdy(i,j)   &
                         )*(sigmaf(k)-zt)*gz(i,j)*rzt
        enddo
        enddo
        do i=0,ni+1
          rrw(i,j,   1) = 0.0
          rrw(i,j,nk+1) = 0.0
        enddo
      ENDDO

    ENDIF
    if(timestats.ge.1) time_advs=time_advs+mytime()

        IF(terrain_flag)THEN
          call bcw(rrw,0)
#ifdef MPI
          call comm_1w_start(rrw,ww1,ww2,we1,we2,   &
                                 ws1,ws2,wn1,wn2,reqs_w)
          call comm_1w_end(rrw,ww1,ww2,we1,we2,   &
                               ws1,ws2,wn1,wn2,reqs_w)
#endif
        ENDIF

      IF(.not.terrain_flag)THEN

!------------
        IF(axisymm.eq.0)THEN

          ! Cartesian without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            divx(i,j,k)=(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)        &
                       +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)        &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)
          enddo
          enddo
          enddo

        ELSE

          ! axisymmetric:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            divx(i,j,k)=(xf(i+1)*rru(i+1,j,k)-xf(i)*rru(i,j,k))*rdx*uh(i)*rxh(i)   &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdz*mh(1,1,k)
          enddo
          enddo
          enddo

        ENDIF
!------------

      ELSE

          ! Cartesian with terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            divx(i,j,k)=(rru(i+1,j,k)-rru(i,j,k))*rdx*uh(i)        &
                       +(rrv(i,j+1,k)-rrv(i,j,k))*rdy*vh(j)        &
                       +(rrw(i,j,k+1)-rrw(i,j,k))*rdsf(k)
          enddo
          enddo
          enddo

      ENDIF
      if(timestats.ge.1) time_divx=time_divx+mytime()

!--------------------------------------------------------------------
!  finish comms for q/theta:
#ifdef MPI
        IF(nrk.ge.2)THEN
          if(imoist.eq.1)then
            call comm_3q_end(q3d,qw31,qw32,qe31,qe32,   &
                                 qs31,qs32,qn31,qn32,reqs_q(1,1))
          endif
          call comm_3r_end(th3d,pp3d,rw31,rw32,re31,re32,   &
                                     rs31,rs32,rn31,rn32,reqs_p)
        ENDIF
#endif
!--------------------------------------------------------------------
!  Calculate misc. variables
!
!    These arrays store variables that are used later in the
!    SOUND subroutine.  Do not modify t11 or t22 until after sound!
!    dum1 = vapor
!    dum2 = all liquid
!    dum3 = all solid

        IF(imoist.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,cpm,cvm)
          do k=1,nk

            do j=0,nj+1
            do i=0,ni+1
              dum2(i,j,k)=0.0
              dum3(i,j,k)=0.0
            enddo
            enddo
            do n=nql1,nql2
              do j=0,nj+1
              do i=0,ni+1
                dum2(i,j,k)=dum2(i,j,k)+q3d(i,j,k,n)
              enddo
              enddo
            enddo
            IF(iice.eq.1)THEN
              do n=nqs1,nqs2
                do j=0,nj+1
                do i=0,ni+1
                  dum3(i,j,k)=dum3(i,j,k)+q3d(i,j,k,n)
                enddo
                enddo
              enddo
            ENDIF
            ! save qv,ql,qi for buoyancy calculation:
            do j=0,nj+1
            do i=0,ni+1
              t12(i,j,k)=max(q3d(i,j,k,nqv),0.0)
              t13(i,j,k)=max(0.0,dum2(i,j,k))+max(0.0,dum3(i,j,k))
              t11(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*(1.0+reps*t12(i,j,k))     &
                         /(1.0+t12(i,j,k)+t13(i,j,k))
      ! terms in theta and pi equations for proper mass/energy conservation
      ! Reference:  Bryan and Fritsch (2002, MWR)
              dum4(i,j,k)=cpl*max(0.0,dum2(i,j,k))+cpi*max(0.0,dum3(i,j,k))
              cpm=cp+cpv*t12(i,j,k)+dum4(i,j,k)
              cvm=cv+cvv*t12(i,j,k)+dum4(i,j,k)
              thterm(i,j,k)=( rd+rv*q3d(i,j,k,nqv)-rovcp*cpm )/cvm
              t22(i,j,k)=rovcp*cpm/cvm
            enddo
            enddo

          enddo


        ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=0,nj+1
          do i=0,ni+1
            t11(i,j,k)=th0(i,j,k)+th3d(i,j,k)
          enddo
          enddo
          enddo

        ENDIF

        if(timestats.ge.1) time_buoyan=time_buoyan+mytime()


!--------------------------------------------------------------------
! Moisture

  IF(imoist.eq.1)THEN

    DO n=1,numq

      ! t33 = dummy

      bflag=0
      if(stat_qsrc.eq.1 .and. nrk.eq.3) bflag=1

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=qten(i,j,k,n)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_misc=time_misc+mytime()


      if(nrk.eq.3)then
        pdef = 1
      else
        pdef = 0
      endif

      weps = 0.01*epsilon
      IF( idm.eq.1 .and. n.ge.nnc1 ) weps = 1.0e5*epsilon
      call advs(nrk,1,bflag,bsq(n),xh,rxh,uh,ruh,xf,vh,rvh,gz,rgz,mh,rmh,    &
                 rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,             &
                 rru,rrv,rrw,qa(ib,jb,kb,n),q3d(ib,jb,kb,n),sten,pdef,dttmp,weps, &
                 flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,rdsf)

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        q3d(i,j,k,n)=qa(i,j,k,n)+dttmp*sten(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

      IF(nrk.lt.3)THEN
        call bcs(q3d(ib,jb,kb,n))
      ENDIF

    ENDDO   ! enddo for n loop

#ifdef MPI
    ! dont communicate on last rk step
    IF(nrk.le.2)THEN
      call comm_3q_start(q3d,qw31,qw32,qe31,qe32,   &
                             qs31,qs32,qn31,qn32,reqs_q(1,1))
    endif
#endif

  ENDIF    ! endif for imoist=1

!--------------------------------------------------------------------
!  Passive Tracers

    if(iptra.eq.1)then
    DO n=1,npt

      ! t33 = dummy

      bflag=0
      if(stat_qsrc.eq.1 .and. nrk.eq.3) bflag=1

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        sten(i,j,k)=ptten(i,j,k,n)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_misc=time_misc+mytime()


      if(nrk.eq.3)then
        pdef = 1
      else
        pdef = 0
      endif

#ifdef MPI
          IF(nrk.ge.2)THEN
            call comm_3s_end(pt3d(ib,jb,kb,n),                           &
                  tw1(1,1,1,n),tw2(1,1,1,n),te1(1,1,1,n),te2(1,1,1,n),   &
                  ts1(1,1,1,n),ts2(1,1,1,n),tn1(1,1,1,n),tn2(1,1,1,n),   &
                  reqs_t(1,n))
          ENDIF
#endif

      weps = 1.0*epsilon
      call advs(nrk,1,bflag,bfoo,xh,rxh,uh,ruh,xf,vh,rvh,gz,rgz,mh,rmh,        &
                 rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,               &
                 rru,rrv,rrw,pta(ib,jb,kb,n),pt3d(ib,jb,kb,n),sten,pdef,dttmp,weps, &
                 flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,rdsf)

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        pt3d(i,j,k,n)=pta(i,j,k,n)+dttmp*sten(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()

        if(nrk.eq.3) call pdefq(0.0,afoo,ruh,rvh,rmh,rho,pt3d(ib,jb,kb,n))

        call bcs(pt3d(ib,jb,kb,n))
#ifdef MPI
        call comm_3s_start(pt3d(ib,jb,kb,n)   &
                     ,tw1(1,1,1,n),tw2(1,1,1,n),te1(1,1,1,n),te2(1,1,1,n)     &
                     ,ts1(1,1,1,n),ts2(1,1,1,n),tn1(1,1,1,n),tn2(1,1,1,n)     &
                     ,reqs_t(1,n) )
#endif

    ENDDO
    endif

!--------------------------------------------------------------------
!  TKE advection
 
        IF(iturb.eq.1)THEN

          ! use wten for tke tendency, step tke forward:

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            wten(i,j,k)=tketen(i,j,k)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()

#ifdef MPI
        IF(nrk.ge.2)THEN
          call comm_3t_end(tke3d,tkw1,tkw2,tke1,tke2,   &
                                 tks1,tks2,tkn1,tkn2,reqs_tk)
        ENDIF
#endif

            call advw(xh,rxh,uh,xf,vh,gz,rgz,mf,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,  &
                       rru,rrv,rrw,tke3d,wten,rds)

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
        do k=1,nk+1
          do j=1,nj
          do i=1,ni
            tke3d(i,j,k)=tkea(i,j,k)+dttmp*wten(i,j,k)
          enddo
          enddo
          if(nrk.eq.3)then
            do j=1,nj
            do i=1,ni
              if(tke3d(i,j,k).lt.1.0e-6) tke3d(i,j,k)=0.0
            enddo
            enddo
          endif
        enddo
        if(timestats.ge.1) time_integ=time_integ+mytime()


          call bct(tke3d)
#ifdef MPI
          call comm_3t_start(tke3d,tkw1,tkw2,tke1,tke2,   &
                                   tks1,tks2,tkn1,tkn2,reqs_tk)
#endif

        ENDIF


!--------------------------------------------------------------------
!  THETA-equation

      ! t23  = theta used for advection
      ! t33  = dummy

!$omp parallel do default(shared)   &
!$omp private(i,j,k)
    do k=1,nk
      do j=1,nj
      do i=1,ni
        thten(i,j,k)=thten1(i,j,k)
      enddo
      enddo
      do j=jb,je
      do i=ib,ie
        ! new approach (cm1r16)
        t23(i,j,k)=(th0(i,j,k)-th00(i,j,k))+th3d(i,j,k)
      enddo
      enddo
    enddo
      if(timestats.ge.1) time_misc=time_misc+mytime()


        weps = 10.0*epsilon
        call advs(nrk,1,0,bfoo,xh,rxh,uh,ruh,xf,vh,rvh,gz,rgz,mh,rmh, &
                   rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,   &
                   rru,rrv,rrw,tha,t23,thten,0,dttmp,weps,           &
                   flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,rdsf)


!--------------------------------------------------------------------
!  Pressure equation

      ! t22  = ppterm
      ! t23  = theta used for advection
      ! t33  = dummy

      IF(psolver.le.3)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
        do j=1,nj
        do i=1,ni
          sten(i,j,k)=ppten(i,j,k)
        enddo
        enddo
        do j=jb,je
        do i=ib,ie
          ! new approach (cm1r16)
          t23(i,j,k)=(pi0(i,j,k)-pi00(i,j,k))+pp3d(i,j,k)
        enddo
        enddo
      enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()

          weps = epsilon
          call advs(nrk,0,0,bfoo,xh,rxh,uh,ruh,xf,vh,rvh,gz,rgz,mh,rmh, &
                     rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,t33,   &
                     rru,rrv,rrw,ppi,t23,sten,0,dttmp,weps,            &
                     flag,sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,rdsf)

      ENDIF

!--------------------------------------------------------------------
!  U-equation

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          uten(i,j,k)=uten1(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()
 
        if(icor.eq.1)then

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
            if(pertcor.eq.1)then
              do j=1,nj+1
              do i=0,ni+1
                dum1(i,j,k)=v3d(i,j,k)-v0(i,j,k)
              enddo
              enddo
            else
              do j=1,nj+1
              do i=0,ni+1
                dum1(i,j,k)=v3d(i,j,k)
              enddo
              enddo
            endif
            do j=1,nj
            do i=1,ni+1
              uten(i,j,k)=uten(i,j,k)+fcor*             &
               0.25*( (dum1(i  ,j,k)+dum1(i  ,j+1,k))   &
                     +(dum1(i-1,j,k)+dum1(i-1,j+1,k)) )
            enddo
            enddo
          enddo
          if(timestats.ge.1) time_cor=time_cor+mytime()

        endif

        if(axisymm.eq.1)then

          ! for axisymmetric grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=2,ni+1
            uten(i,j,k)=uten(i,j,k)+0.5*(   &
                 ( v3d(i-1,j,k)**2)*rxh(i-1)+(v3d(i,j,k)**2)*rxh(i) )
          enddo
          enddo
          enddo

        endif

          call advu(xf,rxf,uf,vh,gz,rgz,gzu,mh,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,  &
                     rru,u3d,uten,rrv,rrw,rdsf)

!--------------------------------------------------------------------
!  V-equation
 
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          vten(i,j,k)=vten1(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()
 
        if(icor.eq.1)then

          IF(axisymm.eq.0)THEN
            ! for Cartesian grid:

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
            if(pertcor.eq.1)then
              do j=0,nj+1
              do i=1,ni+1
                dum1(i,j,k)=u3d(i,j,k)-u0(i,j,k)
              enddo
              enddo
            else
              do j=0,nj+1
              do i=1,ni+1
                dum1(i,j,k)=u3d(i,j,k)
              enddo
              enddo
            endif
            do j=1,nj+1
            do i=1,ni
              vten(i,j,k)=vten(i,j,k)-fcor*             &
               0.25*( (dum1(i,j  ,k)+dum1(i+1,j  ,k))   &
                     +(dum1(i,j-1,k)+dum1(i+1,j-1,k)) )
            enddo
            enddo
          enddo

          ELSE
            ! for axisymmetric grid:

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
            do j=1,nj
            do i=1,ni
              vten(i,j,k)=vten(i,j,k)-fcor*0.5*(xf(i)*u3d(i,j,k)+xf(i+1)*u3d(i+1,j,k))*rxh(i)
            enddo
            enddo
          enddo

          ENDIF
          if(timestats.ge.1) time_cor=time_cor+mytime()

        endif

        if(axisymm.eq.1)then
          ! for axisymmetric grid:

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            vten(i,j,k)=vten(i,j,k)-(v3d(i,j,k)*rxh(i))*0.5*(xf(i)*u3d(i,j,k)+xf(i+1)*u3d(i+1,j,k))*rxh(i)
          enddo
          enddo
          enddo

        endif

          call advv(xh,rxh,uh,xf,vf,gz,rgz,gzv,mh,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,  &
                     rru,rrv,v3d,vten,rrw,rdsf)

!--------------------------------------------------------------------
!  W-equation

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          wten(i,j,k)=wten1(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_misc=time_misc+mytime()
 
        if( imoist.eq.1 )then
          ! buoyancy (moisture terms):
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do j=1,nj
            do k=1,nk
            do i=1,ni
              dum2(i,j,k) = repsm1*(t12(i,j,k)-qv0(i,j,k)) - (t13(i,j,k)-qc0(i,j,k))
            enddo
            enddo
            do k=2,nk
            do i=1,ni
              wten(i,j,k)=wten(i,j,k)+g*(c1(i,j,k)*dum2(i,j,k-1)+c2(i,j,k)*dum2(i,j,k))
            enddo
            enddo
          enddo
          if(timestats.ge.1) time_buoyan=time_buoyan+mytime()
        endif

        if(psolver.eq.1.or.psolver.eq.4.or.psolver.eq.5)then
          ! buoyancy for non-time-split solvers
          ! (i.e.,truly-compressible/anelastic/incompressible solvers:)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do j=1,nj
            do k=1,nk
            do i=1,ni
              dum2(i,j,k) = th3d(i,j,k)/th0(i,j,k)
            enddo
            enddo
            do k=2,nk
            do i=1,ni
              wten(i,j,k)=wten(i,j,k)+g*(c1(i,j,k)*dum2(i,j,k-1)+c2(i,j,k)*dum2(i,j,k))
            enddo
            enddo
          enddo
          if(timestats.ge.1) time_buoyan=time_buoyan+mytime()
        endif

          call advw(xh,rxh,uh,xf,vh,gz,rgz,mf,rho0,rr0,rf0,rrf0,dum1,dum2,dum3,dum4,divx,  &
                     rru,rrv,rrw,w3d,wten,rds)


!--------------------------------------------------------------------
!  Update v for axisymmetric model simulations:

        IF(axisymm.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            v3d(i,j,k)=va(i,j,k)+dttmp*vten(i,j,k)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_misc=time_misc+mytime()

          call bcv(v3d)

        ENDIF

!--------------------------------------------------------------------
!  call sound

        IF(psolver.eq.1)THEN

          call soundns(xh,rxh,uh,xf,uf,yh,vh,yf,vf,zh,mh,c1,c2,mf, &
                       pi0,thv0,rr0,rf0,                           &
                       radbcw,radbce,radbcs,radbcn,                &
                       divx,u0,ua,u3d,uten,v0,va,v3d,vten,wa,w3d,wten,   &
                       ppi,pp3d,sten,t11,t22,dttmp,nrk,rtime,      &
                       th0,th00,tha,th3d,thten,thterm,             &
                       reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,         &
                       uw31,uw32,ue31,ue32,us31,us32,un31,un32,    &
                       vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,    &
                       ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,    &
                       sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,    &
                       rw31,rw32,re31,re32,rs31,rs32,rn31,rn32,    &
                       pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)

        ELSEIF(psolver.eq.2)THEN

          call sounde(dt,th00s,thlr,xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf, &
                     rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,zf,mf,rmf,    &
                     pi0,th00,thv0,rho0,rr0,rf0,th0,zs,                  &
                     gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,   &
                     radbcw,radbce,radbcs,radbcn,                    &
                     dum1,dum2,dum3,dum4,divx,t12,t13,               &
                     u0,ua,u3d,uten,v0,va,v3d,vten,wa,w3d,wten,      &
                     ppi,pp3d,sten,tha,th3d,thten,thterm,            &
                     t11,t22,nrk,rtime,                              &
                     reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,             &
                     uw31,uw32,ue31,ue32,us31,us32,un31,un32,        &
                     vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,        &
                     ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,        &
                     sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,        &
                     rw31,rw32,re31,re32,rs31,rs32,rn31,rn32,        &
                     pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)

        ELSEIF(psolver.eq.3)THEN

          ! rho,prs,divx are used as dummy arrays by sound

          call sound(dt,th00s,thlr,xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,yf,vf, &
                     rds,sigma,rdsf,sigmaf,zh,mh,rmh,c1,c2,zf,mf,rmf,   &
                     pi0,th00,thv0,rho0,rr0,rf0,th0,zs,                 &
                     gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,  &
                     radbcw,radbce,radbcs,radbcn,                    &
                     dum1,dum2,dum3,dum4,t12,t13,t23,t33,            &
                     u0,ua,u3d,uten,v0,va,v3d,vten,wa,w3d,wten,      &
                     ppi,pp3d,sten,tha,th3d,thten,thterm,            &
                     t11,t22,rho,prs,divx,nrk,rtime,                 &
                     reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,             &
                     uw31,uw32,ue31,ue32,us31,us32,un31,un32,        &
                     vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,        &
                     ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,        &
                     sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,        &
                     rw31,rw32,re31,re32,rs31,rs32,rn31,rn32,        &
                     pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2)

        ELSEIF(psolver.eq.4.or.psolver.eq.5)THEN
          ! anelastic/incompressible solver:

          call anelp(xh,uh,xf,uf,yh,vh,yf,vf,                     &
                     zh,mh,rmh,mf,rmf,pi0,thv0,rho0,prs0,rf0,     &
                     radbcw,radbce,radbcs,radbcn,dum1,divx,       &
                     u0,ua,u3d,uten,v0,va,v3d,vten,wa,w3d,wten,   &
                     ppi,pp3d,tha,th3d,thten,t11,cfb,cfa,cfc,     &
                     ad1,ad2,pdt,deft,rhs,trans,dttmp,nrk,rtime)

        ENDIF

!--------------------------------------------------------------------
!  radbc

        if(irbc.eq.4)then

          if(ibw.eq.1 .or. ibe.eq.1)then
            call radbcew4(ruf,radbcw,radbce,ua,u3d,dttmp)
          endif

          if(ibs.eq.1 .or. ibn.eq.1)then
            call radbcns4(rvf,radbcs,radbcn,va,v3d,dttmp)
          endif

        endif

!--------------------------------------------------------------------
! RK loop end

      ENDDO


!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   End of RK section   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(imoist.eq.1.and.neweqts.ge.1)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
          do j=1,nj
          do i=1,ni
            pp3d(i,j,k)     = pp3d(i,j,k)    -dt*qpten(i,j,k)
            th3d(i,j,k)     = th3d(i,j,k)    -dt*qtten(i,j,k)
             q3d(i,j,k,nqv) =  q3d(i,j,k,nqv)-dt*qvten(i,j,k)
             q3d(i,j,k,nqc) =  q3d(i,j,k,nqc)-dt*qcten(i,j,k)
          enddo
          enddo
          IF(ptype.eq.2)THEN
            do j=1,nj
            do i=1,ni
              q3d(i,j,k,nqi) =  q3d(i,j,k,nqi)-dt*qiten(i,j,k)
            enddo
            enddo
          ENDIF
        enddo
      ENDIF

!--------------------------------------------------------------------
!  Get pressure
!  Get density

    IF(psolver.eq.4.or.psolver.eq.5)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=1,nk
      do j=1,nj
      do i=1,ni
        prs(i,j,k)=prs0(i,j,k)
        rho(i,j,k)=rho0(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_misc=time_misc+mytime()

    ELSE

      IF(imoist.eq.1)THEN

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
          rho(i,j,k)=prs(i,j,k)                         &
             /( rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))     &
                  *(1.0+max(0.0,q3d(i,j,k,nqv))*reps) )
        enddo
        enddo
        enddo

      ELSE

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
          rho(i,j,k)=prs(i,j,k)   &
             /(rd*(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k)))
        enddo
        enddo
        enddo

      ENDIF

    ENDIF

      if(timestats.ge.1) time_prsrho=time_prsrho+mytime()

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   BEGIN microphysics   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      IF(imoist.eq.1)THEN

        getdbz = .false.
        IF(output_dbz.eq.1)THEN
          rtime=sngl(mtime+dt)
          if( (rtime.ge.sngl(taptim)).or.stopit )then
            getdbz = .true.
          endif
          if(getdbz)then
            if(dowr) write(outfile,*) '  Getting dbz ... '
          endif
        ENDIF

        ! sten = dbz
        ! dum1 = t
        ! dum2 = cvm
        ! dum3 = appropriate rho for budget calculations

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n)
        DO k=1,nk

          do j=1,nj
          do i=1,ni
            dum1(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*(pi0(i,j,k)+pp3d(i,j,k))
          enddo
          enddo

          if(getdbz)then
            do j=1,nj
            do i=1,ni
              sten(i,j,k)=0.0
            enddo
            enddo
          endif

          IF(efall.eq.1)THEN
            if(neweqts.eq.0)then
              do j=1,nj
              do i=1,ni
                dum2(i,j,k)=0.0
                dum3(i,j,k)=0.0
                dum4(i,j,k)=0.0
              enddo
              enddo
            else
              do j=1,nj
              do i=1,ni
                dum2(i,j,k)=q3d(i,j,k,nqv)
                dum3(i,j,k)=0.0
                dum4(i,j,k)=0.0
              enddo
              enddo
              do n=nql1,nql2
                do j=1,nj
                do i=1,ni
                  dum3(i,j,k)=dum3(i,j,k)+q3d(i,j,k,n)
                enddo
                enddo
              enddo
              IF(iice.eq.1)THEN
                do n=nqs1,nqs2
                do j=1,nj
                do i=1,ni
                  dum4(i,j,k)=dum4(i,j,k)+q3d(i,j,k,n)
                enddo
                enddo
                enddo
              ENDIF
            endif
            do j=1,nj
            do i=1,ni
              dum2(i,j,k)=cv+cvv*dum2(i,j,k)+cpl*dum3(i,j,k)+cpi*dum4(i,j,k)
            enddo
            enddo
          ENDIF

          IF(axisymm.eq.0)THEN
            ! for Cartesian grid:
            do j=1,nj
            do i=1,ni
              dum3(i,j,k)=rho(i,j,k)
            enddo
            enddo
          ELSE
            ! for axisymmetric grid:
            do j=1,nj
            do i=1,ni
              dum3(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
            enddo
            enddo
          ENDIF

        ENDDO


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  NOTES:
!           sten       is used for     dbz
!
!           dum1   is   T
!           dum2   is   cvm
!           dum3   is   rho for budget calculations
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Kessler scheme   cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        IF(ptype.eq.1)THEN
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq(1.0e-14,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq(1.0e-14,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call k_fallout(rho,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call geterain(dt,cpl,lv1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          if(efall.ge.1)then
            call getefall(dt,cpl,ruh,rvh,mf,pi0,th0,dum1,dum2,dum3,   &
                          pp3d,th3d,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call kessler(dt,qbudget(3),qbudget(4),qbudget(5),ruh,rvh,rmh,pi0,th0,dum1,   &
                       rho,dum3,pp3d,th3d,prs,                            &
                       q3d(ib,jb,kb,nqv),q3d(ib,jb,kb,2),q3d(ib,jb,kb,3))
          call satadj(4,dt,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                      rho,dum3,pp3d,prs,th3d,q3d,qpten,qtten,qvten,qcten)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Goddard LFO scheme   cccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ELSEIF(ptype.eq.2)THEN
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq(1.0e-14,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq(1.0e-14,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq(1.0e-14,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq(1.0e-14,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq(1.0e-14,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
          call goddard(dt,qbudget(3),qbudget(4),qbudget(5),ruh,rvh,rmh,pi0,th0,             &
                       rho,dum3,prs,pp3d,th3d,                            &
     q3d(ib,jb,kb,1), q3d(ib,jb,kb,2),q3d(ib,jb,kb,3),qten(ib,jb,kb,3),   &
     q3d(ib,jb,kb,4),qten(ib,jb,kb,4),q3d(ib,jb,kb,5),qten(ib,jb,kb,5),   &
     q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          call satadj_ice(4,dt,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,     &
                          rho,dum3,pp3d,prs,th3d,                     &
              q3d(ib,jb,kb,1),q3d(ib,jb,kb,2),q3d(ib,jb,kb,3),   &
              q3d(ib,jb,kb,4),q3d(ib,jb,kb,5),q3d(ib,jb,kb,6),   &
              qpten,qtten,qvten,qcten,qiten)
          call geterain(dt,cpl,lv1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call geterain(dt,cpi,ls1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,4),qten(ib,jb,kb,4))
          call geterain(dt,cpi,ls1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,5),qten(ib,jb,kb,5))
          call geterain(dt,cpi,ls1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          if(efall.ge.1)then
            call getefall(dt,cpl,ruh,rvh,mf,pi0,th0,dum1,dum2,dum3,   &
                          pp3d,th3d,q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
            call getefall(dt,cpi,ruh,rvh,mf,pi0,th0,dum1,dum2,dum3,   &
                          pp3d,th3d,q3d(ib,jb,kb,4),qten(ib,jb,kb,4))
            call getefall(dt,cpi,ruh,rvh,mf,pi0,th0,dum1,dum2,dum3,   &
                          pp3d,th3d,q3d(ib,jb,kb,5),qten(ib,jb,kb,5))
            call getefall(dt,cpi,ruh,rvh,mf,pi0,th0,dum1,dum2,dum3,   &
                          pp3d,th3d,q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,3),qten(ib,jb,kb,3))
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,4),qten(ib,jb,kb,4))
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,5),qten(ib,jb,kb,5))
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,6),qten(ib,jb,kb,6))
          if(getdbz) call calcdbz(rho,q3d(ib,jb,kb,3),q3d(ib,jb,kb,5),q3d(ib,jb,kb,6),sten)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Thompson scheme   ccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ELSEIF(ptype.eq.3)THEN
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq(1.0e-12,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq(1.0e-12,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq(1.0e-12,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq(1.0e-12,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq(1.0e-12,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
!!!          call pdefq(    1.0,asq(7),ruh,rvh,rmh,rho,q3d(ib,jb,kb,7))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! dum1 = pi
            ! dum2 = dz
            ! dum4 = T
            dum1(i,j,k)=pi0(i,j,k)+pp3d(i,j,k)
            dum2(i,j,k)=dz*rmh(i,j,k)
            dum4(i,j,k)=(th0(i,j,k)+th3d(i,j,k))*dum1(i,j,k)
            ! store old T in thten array:
            thten(i,j,k)=dum4(i,j,k)
          enddo
          enddo
          enddo
          call mp_gt_driver(q3d(ib,jb,kb,1),q3d(ib,jb,kb,2),q3d(ib,jb,kb,3), &
                            q3d(ib,jb,kb,4),q3d(ib,jb,kb,5),q3d(ib,jb,kb,6), &
                            q3d(ib,jb,kb,7),q3d(ib,jb,kb,8),                 &
                            th0,dum4,dum1,prs,dum2,dt,rain,                 &
                            qbudget(5),qbudget(6),                           &
                            ruh,rvh,rmh,dum3,sten,getdbz)
        IF(neweqts.ge.1)THEN
          ! for mass conservation:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( abs(dum4(i,j,k)-thten(i,j,k)).ge.tsmall )then
!!!              th3d(i,j,k)=th3d(i,j,k)+(dum4(i,j,k)-thten(i,j,k))
!!!              pp3d(i,j,k)=((rho(i,j,k)*(rd+rv*q3d(i,j,k,nqv))   &
!!!                                      *(th0(i,j,k)+th3d(i,j,k))*rp00)**rddcv)-pi0(i,j,k)
!!!              prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
              prs(i,j,k)=rho(i,j,k)*rd*dum4(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps)
              pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
              th3d(i,j,k)=dum4(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
            endif
          enddo
          enddo
          enddo
        ELSE
          ! traditional thermodynamics:  p,pi remain unchanged
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( abs(dum4(i,j,k)-thten(i,j,k)).ge.tsmall )then
              th3d(i,j,k)=th3d(i,j,k)+(dum4(i,j,k)-thten(i,j,k))/dum1(i,j,k)
              rho(i,j,k)=prs(i,j,k)/(rd*dum4(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps))
            endif
          enddo
          enddo
          enddo
        ENDIF
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
          call satadj(4,dt,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                      rho,dum3,pp3d,prs,th3d,q3d,qpten,qtten,qvten,qcten)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   GSR LFO scheme   cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.4)THEN
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq(1.0e-14,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq(1.0e-14,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq(1.0e-14,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq(1.0e-14,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq(1.0e-14,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
          call lfo_ice_drive(dt, mf, pi0, prs0, pp3d, prs, th0, th3d,    &
                             qv0, rho0, q3d, qten, dum1)
          do n=2,numq
            call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                         q3d(ib,jb,kb,n),qten(ib,jb,kb,n))
          enddo

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   Morrison scheme   cccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.5)THEN
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq(1.0e-12,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
          call pdefq(1.0e-12,asq(3),ruh,rvh,rmh,rho,q3d(ib,jb,kb,3))
          call pdefq(1.0e-12,asq(4),ruh,rvh,rmh,rho,q3d(ib,jb,kb,4))
          call pdefq(1.0e-12,asq(5),ruh,rvh,rmh,rho,q3d(ib,jb,kb,5))
          call pdefq(1.0e-12,asq(6),ruh,rvh,rmh,rho,q3d(ib,jb,kb,6))
!!!          call pdefq(    1.0,asq(7),ruh,rvh,rmh,rho,q3d(ib,jb,kb,7))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ! dum1 = T  (this should have been calculated already)
            ! dum4 = dz
            dum4(i,j,k)=dz*rmh(i,j,k)
            ! store old T in thten array:
            thten(i,j,k)=dum1(i,j,k)
          enddo
          enddo
          enddo
          IF(numq.eq.11)THEN
            ! ppten stores ncc:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              ppten(i,j,k) = q3d(i,j,k,11)
            enddo
            enddo
            enddo
          ENDIF
          call MP_GRAUPEL(nstep,dum1,                                 &
                          q3d(ib,jb,kb, 1),q3d(ib,jb,kb, 2),q3d(ib,jb,kb, 3), &
                          q3d(ib,jb,kb, 4),q3d(ib,jb,kb, 5),q3d(ib,jb,kb, 6), &
                          q3d(ib,jb,kb, 7),q3d(ib,jb,kb, 8),q3d(ib,jb,kb, 9), &
                          q3d(ib,jb,kb,10),ppten,                             &
                               prs,rho,dt,dum4,w3d,rain,                      &
                          qbudget(1),qbudget(2),qbudget(5),qbudget(6),        &
                          ruh,rvh,rmh,dum3,sten,getdbz)
          IF(numq.eq.11)THEN
            ! ppten stores ncc:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
            do k=1,nk
            do j=1,nj
            do i=1,ni
              q3d(i,j,k,11) = ppten(i,j,k)
            enddo
            enddo
            enddo
          ENDIF
        IF(neweqts.ge.1)THEN
          ! for mass conservation:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall )then
!!!              th3d(i,j,k)=th3d(i,j,k)+(dum1(i,j,k)-thten(i,j,k))
!!!              pp3d(i,j,k)=((rho(i,j,k)*(rd+rv*q3d(i,j,k,nqv))   &
!!!                                      *(th0(i,j,k)+th3d(i,j,k))*rp00)**rddcv)-pi0(i,j,k)
!!!              prs(i,j,k)=p00*((pi0(i,j,k)+pp3d(i,j,k))**cpdrd)
              prs(i,j,k)=rho(i,j,k)*rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps)
              pp3d(i,j,k)=(prs(i,j,k)*rp00)**rovcp - pi0(i,j,k)
              th3d(i,j,k)=dum1(i,j,k)/(pi0(i,j,k)+pp3d(i,j,k)) - th0(i,j,k)
            endif
          enddo
          enddo
          enddo
        ELSE
          ! traditional thermodynamics:  p,pi remain unchanged
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if( abs(dum1(i,j,k)-thten(i,j,k)).ge.tsmall )then
              th3d(i,j,k)=th3d(i,j,k)+(dum1(i,j,k)-thten(i,j,k))/(pi0(i,j,k)+pp3d(i,j,k))
              rho(i,j,k)=prs(i,j,k)/(rd*dum1(i,j,k)*(1.0+q3d(i,j,k,nqv)*reps))
            endif
          enddo
          enddo
          enddo
        ENDIF
          if(timestats.ge.1) time_microphy=time_microphy+mytime()
          call satadj(4,dt,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                      rho,dum3,pp3d,prs,th3d,q3d,qpten,qtten,qvten,qcten)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccc   RE87 scheme   ccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        ELSEIF(ptype.eq.6)THEN
          call pdefq(    0.0,asq(1),ruh,rvh,rmh,rho,q3d(ib,jb,kb,1))
          call pdefq(1.0e-14,asq(2),ruh,rvh,rmh,rho,q3d(ib,jb,kb,2))
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            if(q3d(i,j,k,2).gt.0.001)then
              qten(i,j,k,2) = v_t
            else
              qten(i,j,k,2) = 0.0
            endif
          enddo
          enddo
          enddo
          call geterain(dt,cpl,lv1,qbudget(7),ruh,rvh,dum1,dum3,q3d(ib,jb,kb,2),qten(ib,jb,kb,2))
          if(efall.ge.1)then
            call getefall(dt,cpl,ruh,rvh,mf,pi0,th0,dum1,dum2,dum3,   &
                          pp3d,th3d,q3d(ib,jb,kb,2),qten(ib,jb,kb,2))
          endif
          call fallout(dt,qbudget(6),ruh,rvh,zh,mh,mf,rain,dum3,rho,   &
                       q3d(ib,jb,kb,2),qten(ib,jb,kb,2))
          call satadj(4,dt,qbudget(1),qbudget(2),ruh,rvh,rmh,pi0,th0,   &
                      rho,dum3,pp3d,prs,th3d,q3d,qpten,qtten,qvten,qcten)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Milbrandt & Yao scheme
!
!        ELSEIF(ptype.eq.7)THEN
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  Ziegler/Mansell two-moment scheme
!
!        ELSEIF(ptype.ge.26)THEN
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  insert new microphysics schemes here
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!        ELSEIF(ptype.eq.8)THEN
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! otherwise, stop for undefined ptype
        ELSE
          print *,'  Undefined ptype!'
          call stopcm1
        ENDIF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CC   END microphysics   CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!Begin:  message passing

      IF(imoist.eq.1)THEN
          call bcs(th3d)
          call bcs(pp3d)
#ifdef MPI
          call comm_3r_start(th3d,pp3d,rw31,rw32,re31,re32,   &
                                       rs31,rs32,rn31,rn32,reqs_p)
#endif
          DO n=1,numq
            call bcs(q3d(ib,jb,kb,n))
          ENDDO
#ifdef MPI
          call comm_3q_start(q3d,qw31,qw32,qe31,qe32,   &
                                 qs31,qs32,qn31,qn32,reqs_q(1,1))
#endif

      ENDIF

      if(iturb.ge.1)then
        call bcs(rho)
#ifdef MPI
        call comm_1s_start(rho,pw1,pw2,pe1,pe2,   &
                               ps1,ps2,pn1,pn2,reqs_s)
#endif
      endif

!Done:  message passing
!-----------------------------------------------------------------
!  Equate the two arrays

      if(iturb.eq.1)then
#ifdef MPI
        call comm_3t_end(tke3d,tkw1,tkw2,tke1,tke2,   &
                               tks1,tks2,tkn1,tkn2,reqs_tk)
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=0,nk+2
        do j=0,nj+1
        do i=0,ni+1
          tkea(i,j,k)=tke3d(i,j,k)
        enddo
        enddo
        enddo
        if(timestats.ge.1) time_integ=time_integ+mytime()
      endif

      if(iptra.eq.1)then
        do n=1,npt
#ifdef MPI
          call comm_3s_end(pt3d(ib,jb,kb,n),                           &
                tw1(1,1,1,n),tw2(1,1,1,n),te1(1,1,1,n),te2(1,1,1,n),   &
                ts1(1,1,1,n),ts2(1,1,1,n),tn1(1,1,1,n),tn2(1,1,1,n),   &
                reqs_t(1,n))
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=0,nk+1
          do j=0,nj+1
          do i=0,ni+1
            pta(i,j,k,n)=pt3d(i,j,k,n)
          enddo
          enddo
          enddo
          if(timestats.ge.1) time_integ=time_integ+mytime()
        enddo
      endif

#ifdef MPI
      call comm_3u_end(u3d,uw31,uw32,ue31,ue32,   &
                           us31,us32,un31,un32,reqs_u)
      call comm_3v_end(v3d,vw31,vw32,ve31,ve32,   &
                           vs31,vs32,vn31,vn32,reqs_v)
      call comm_3w_end(w3d,ww31,ww32,we31,we32,   &
                           ws31,ws32,wn31,wn32,reqs_w)
#endif

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+2
          ua(i,j,k)=u3d(i,j,k)
        enddo
        enddo
        do j=0,nj+2
        do i=0,ni+1
          va(i,j,k)=v3d(i,j,k)
        enddo
        enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()
 
      if(terrain_flag)then
        call bcwsfc(zh,zf,sigma,sigmaf,dzdx,dzdy,dum1,dum2,u3d,v3d,w3d)
        call bc2d(w3d(ib,jb,1))
      endif

!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+2
      do j=0,nj+1
      do i=0,ni+1
        wa(i,j,k)=w3d(i,j,k)
      enddo
      enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()


#ifdef MPI
      call comm_3r_end(th3d,pp3d,rw31,rw32,re31,re32,   &
                                 rs31,rs32,rn31,rn32,reqs_p)
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
      do k=0,nk+1
        do j=0,nj+1
        do i=0,ni+1
          ppi(i,j,k)=pp3d(i,j,k)
        enddo
        enddo
        do j=0,nj+1
        do i=0,ni+1
          tha(i,j,k)=th3d(i,j,k)
        enddo
        enddo
      enddo
      if(timestats.ge.1) time_integ=time_integ+mytime()



      if(imoist.eq.1)then
#ifdef MPI
        call comm_3q_end(q3d,qw31,qw32,qe31,qe32,   &
                             qs31,qs32,qn31,qn32,reqs_q(1,1))
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k,n)
        do k=0,nk+1
          do n=1,numq
          do j=0,nj+1
          do i=0,ni+1
            qa(i,j,k,n)=q3d(i,j,k,n)
          enddo
          enddo
          enddo
        enddo
        if(timestats.ge.1) time_integ=time_integ+mytime()
      endif


      if(iturb.ge.1)then
#ifdef MPI
        call comm_1s_end(rho,pw1,pw2,pe1,pe2,   &
                             ps1,ps2,pn1,pn2,reqs_s)
        call bcs2(rho)
        call getcorner(rho,nw1(1),nw2(1),ne1(1),ne2(1),sw1(1),sw2(1),se1(1),se2(1))
#endif
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=0,nj+1
          k = 1
          do i=0,ni+1
            rr(i,j,k) = 1.0/rho(i,j,k)
            ! 2nd-order extrapolation:
            rf(i,j,1) = (1.75*rho(i,j,1)-rho(i,j,2)+0.25*rho(i,j,3))
          enddo
          do k=2,nk
          do i=0,ni+1
            rr(i,j,k) = 1.0/rho(i,j,k)
            rf(i,j,k) = (c1(i,j,k)*rho(i,j,k-1)+c2(i,j,k)*rho(i,j,k))
          enddo
          enddo
          do i=0,ni+1
            ! 2nd-order extrapolation:
            rf(i,j,nk+1) = (1.75*rho(i,j,nk)-rho(i,j,nk-1)+0.25*rho(i,j,nk-2))
          enddo
        enddo
        if(timestats.ge.1) time_turb=time_turb+mytime()
      endif



!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!cc   All done   cccccccccccccccccccccccccccccccccccccccccccccccccccc
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!!!#ifdef MPI
!!!      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
!!!      if(timestats.ge.1) time_mpb=time_mpb+mytime()
!!!#endif

!  Calculate surface "swaths."  Move surface (if necessary). 

      if(imove.eq.1.and.imoist.eq.1)then
        weps = 10.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,rain(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

    IF( output_sws.eq.1 )THEN

!--------------------------------------------------------------------
! Maximum horizontal wind speed at lowest model level: 
! (include domain movement in calculation)

!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = sqrt( (umove+0.5*(ua(i,j,1)+ua(i+1,j,1)))**2    &
                   +(vmove+0.5*(va(i,j,1)+va(i,j+1,1)))**2 ) 
        do n=1,nrain
          sws(i,j,n)=max(sws(i,j,n),tem)
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        weps = 10.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,sws(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3),  &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

!--------------------------------------------------------------------
!  Maximum vertical vorticity at lowest model level:

  IF(axisymm.eq.0)THEN
    IF(.not.terrain_flag)THEN
      ! Cartesian grid, without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj+1
      do i=1,ni+1
        tem = (va(i,j,1)-va(i-1,j,1))*rdx*uf(i)   &
             -(ua(i,j,1)-ua(i,j-1,1))*rdy*vf(j)
        do n=1,nrain
          svs(i,j,n)=max(svs(i,j,n),tem)
        enddo
      enddo
      enddo
    ELSE
      ! Cartesian grid, with terrain:
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
        ! interior:
        do k=2,2
        r2 = (sigmaf(k)-sigma(k-1))*rds(k)
        r1 = 1.0-r2
        do i=0,ni+2
          dum1(i,j,k) = r1*ua(i,j,k-1)+r2*ua(i,j,k)
          dum2(i,j,k) = r1*va(i,j,k-1)+r2*va(i,j,k)
        enddo
        enddo
      enddo
      k = 1
!$omp parallel do default(shared)  &
!$omp private(i,j,n,r1,tem)
      do j=1,nj+1
      do i=1,ni+1
        r1 = zt/(zt-0.25*((zs(i-1,j-1)+zs(i,j))+(zs(i-1,j)+zs(i,j-1))))
        tem = ( r1*(va(i,j,k)*rgzv(i,j)-va(i-1,j,k)*rgzv(i-1,j))*rdx*uf(i)  &
               +0.5*( (zt-sigmaf(k+1))*(dum2(i-1,j,k+1)+dum2(i,j,k+1))      &
                     -(zt-sigmaf(k  ))*(dum2(i-1,j,k  )+dum2(i,j,k  ))      &
                    )*rdsf(k)*r1*(rgzv(i,j)-rgzv(i-1,j))*rdx*uf(i) )        &
             -( r1*(ua(i,j,k)*rgzu(i,j)-ua(i,j-1,k)*rgzu(i,j-1))*rdy*vf(j)  &
               +0.5*( (zt-sigmaf(k+1))*(dum1(i,j-1,k+1)+dum1(i,j,k+1))      &
                     -(zt-sigmaf(k  ))*(dum1(i,j-1,k  )+dum1(i,j,k  ))      &
                    )*rdsf(k)*r1*(rgzu(i,j)-rgzu(i,j-1))*rdy*vf(j) )
        do n=1,nrain
          svs(i,j,n)=max(svs(i,j,n),tem)
        enddo
      enddo
      enddo
    ENDIF
  ELSE
      ! Axisymmetric grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj+1
      do i=1,ni+1
        tem = (va(i,j,1)*xh(i)-va(i-1,j,1)*xh(i-1))*rdx*uf(i)*rxf(i)
        do n=1,nrain
          svs(i,j,n)=max(svs(i,j,n),tem)
        enddo
      enddo
      enddo
  ENDIF

      if(imove.eq.1)then
        weps = 0.01*epsilon
        call movesfc(-1000.0,dt,weps,uh,vh,svs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

!--------------------------------------------------------------------
!  Minimum pressure perturbation at lowest model level:

!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = prs(i,j,1)-prs0(i,j,1)
        do n=1,nrain
          sps(i,j,n)=min(sps(i,j,n),tem)
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        weps = 1000.0*epsilon
        call movesfc(-200000.0,dt,weps,uh,vh,sps(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

!--------------------------------------------------------------------
!  Maximum rainwater mixing ratio (qr) at lowest model level:

    IF(imoist.eq.1.and.nqr.ne.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = qa(i,j,1,nqr)
        do n=1,nrain
          srs(i,j,n)=max(srs(i,j,n),tem)
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        weps = 0.01*epsilon
        call movesfc(0.0,dt,weps,uh,vh,srs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3),  &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif
    ENDIF

!--------------------------------------------------------------------
!  Maximum graupel/hail mixing ratio (qg) at lowest model level:

    IF(imoist.eq.1.and.nqg.ne.0)THEN
!$omp parallel do default(shared)  &
!$omp private(i,j,n,tem)
      do j=1,nj
      do i=1,ni
        tem = qa(i,j,1,nqg)
        do n=1,nrain
          sgs(i,j,n)=max(sgs(i,j,n),tem)
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        weps = 0.01*epsilon
        call movesfc(0.0,dt,weps,uh,vh,sgs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3),  &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif
    ENDIF

!--------------------------------------------------------------------

      ! get height AGL:
      if( terrain_flag )then
        ! get height AGL:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dum3(i,j,k) = zh(i,j,k)-zs(i,j)
          wten(i,j,k) = zf(i,j,k)-zs(i,j)
        enddo
        enddo
        enddo
      else
        ! without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk+1
        do j=1,nj
        do i=1,ni
          dum3(i,j,k) = zh(i,j,k)
          wten(i,j,k) = zf(i,j,k)
        enddo
        enddo
        enddo
      endif

!--------------------------------------------------------------------
!  Maximum updraft velocity (w) at 5 km AGL:

!$omp parallel do default(shared)  &
!$omp private(i,j,k,n,tem)
      do j=1,nj
      do i=1,ni
        k = 2
        ! wten is height AGL:
        do while( wten(i,j,k).lt.5000.0 .and. k.lt.nk )
          k = k + 1
        enddo
        tem = w3d(i,j,k)
        do n=1,nrain
          sus(i,j,n)=max(sus(i,j,n),tem)
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        weps = 10.0*epsilon
        call movesfc(-1000.0,dt,weps,uh,vh,sus(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3), &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

!--------------------------------------------------------------------
!  Maximum integrated updraft helicity:

      ! dum3 is zh (agl), wten is zf (agl)
      call calcuh(uf,vf,dum3,wten,ua,va,wa,dum1(ib,jb,1),dum2,t11,t12, &
                  zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf)
!$omp parallel do default(shared)  &
!$omp private(i,j,n)
      do j=1,nj
      do i=1,ni
        do n=1,nrain
          shs(i,j,n)=max(shs(i,j,n),dum1(i,j,1))
        enddo
      enddo
      enddo

      if(imove.eq.1)then
        weps = 100.0*epsilon
        call movesfc(0.0,dt,weps,uh,vh,shs(ib,jb,2),dum1(ib,jb,1),dum1(ib,jb,2),dum1(ib,jb,3),  &
                     reqs_s,sw31(1,1,1),sw32(1,1,1),se31(1,1,1),se32(1,1,1),               &
                            ss31(1,1,1),ss32(1,1,1),sn31(1,1,1),sn32(1,1,1))
      endif

    ENDIF

!  Done with "swaths"
!--------------------------------------------------------------------
!  Step time forward, Get statistics

      mtime = mtime + dt

      if( convinit.eq.1 )then
        if( mtime.gt.convtime ) convinit = 0
      endif

      rtime=sngl(mtime)
      if( rtime.ge.sngl(stattim) .or. statfrq.le.0.0 )then
        IF(axisymm.eq.0)THEN
          ! for Cartesian grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ppten(i,j,k)=rho(i,j,k)
          enddo
          enddo
          enddo
        ELSE
          ! for axisymmetric grid:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            ppten(i,j,k) = rho(i,j,k)*pi*(xf(i+1)**2-xf(i)**2)/(dx*dy)
          enddo
          enddo
          enddo
        ENDIF
        call statpack(nrec,ndt,dt,rtime,adt,acfl,cloudvar,qname,budname,qbudget,asq,bsq, &
                      xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,vf,zh,mh,rmh,mf,    &
                      zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf,               &
                      rstat,pi0,rho0,thv0,th0,qv0,u0,v0,                &
                      dum1,dum2,dum3,dum4,divx,ppten,prs,               &
                      ua,va,wa,ppi,tha,qa,qten,kmh,kmv,khh,khv,tkea,pta,u10,v10)
        stattim=stattim+statfrq
      else
        if( adapt_dt.eq.1 ) call calccfl(1,rstat,dt,acfl,uf,vf,mf,ua,va,wa,0)
      endif

!--------------------------------------------------------------------
!  Writeout and stuff

      rtime=sngl(mtime)
      if(myid.eq.0)then
        if(timeformat.eq.1)then
          write(6,110) nstep,rtime,' sec '
        elseif(timeformat.eq.2)then
          write(6,110) nstep,rtime/60.0,' min '
        elseif(timeformat.eq.3)then
          write(6,110) nstep,rtime/3600.0,' hour'
        elseif(timeformat.eq.4)then
          write(6,110) nstep,rtime/86400.0,' day '
        else
          write(6,110) nstep,rtime,' sec'
        endif
110     format(2x,i12,4x,f18.6,a5)
      endif
      if(timestats.ge.1) time_misc=time_misc+mytime()

      if( (rtime.ge.sngl(taptim)).or.stopit )then
        nwrite=nwrite+1
      IF(output_format.eq.1.or.output_format.eq.2)THEN
        nn = 1
        if(terrain_flag .and. output_interp.eq.1) nn = 2
        if(output_format.eq.2) nn = 1
        DO n=1,nn
          if(n.eq.1)then
            fnum = 51
          else
            fnum = 71
          endif
          call writeout(fnum,nwrite,qname,xh,xf,uf,yh,yf,vf,xfref,yfref,                       &
                        rds,sigma,rdsf,sigmaf,zh,zf,mf,pi0,prs0,rho0,th0,qv0,u0,v0,            &
                        zs,rgzu,rgzv,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,           &
                        cdu,cdv,ce,dum1,dum2,dum3,dum4,                                        &
                        t11,t12,rho,prs,sten,ua,uten,va,vten,wa,wten,ppi,tha,                  &
                        dissten,thpten,qvpten,qcpten,qipten,upten,vpten,                       &
                        lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,  &
                        qa,kmh,kmv,khh,khv,tkea,swten,lwten,radsw,rnflx,radswnet,radlwin,pta,  &
                        num_soil_layers,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br,dat1,dat2,dat3)
        ENDDO
#ifdef HDFOUT
      ELSEIF(output_format.ge.3)THEN
        call writeout_mult_hdf5(rtime,qname,rds,sigma,rdsf,sigmaf,xh,xf,uf,yh,yf,vf,mh,zh,mf,zf, &
                      pi0,prs0,rho0,th0,qv0,u0,v0,                     &
                      zs,rgzu,rgzv,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,cdu,cdv,ce,dum1,dum2,dum3,dum4,  &
                      t11,t12,rho,prs,sten,ua,uten,va,vten,wa,wten,ppi,tha,    &
                      qa,kmh,kmv,khh,khv,tkea,pta,num_soil_layers,   &
                      lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,   &
                      radsw,rnflx,radswnet,radlwin,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br,   &
                      dissten,thpten,qvpten,qcpten,qipten,upten,vpten,swten,lwten)
! See above note about sten being used for dbz calculation
#endif
      ENDIF
        taptim=taptim+tapfrq
        if(timestats.ge.1) time_write=time_write+mytime()
      endif

!-------------------------------------------------------------------
!  Write parcel data:

      if(iprcl.eq.1)then
      IF( rtime.ge.prcltim .or. prclfrq.lt.0.0 )THEN
        call parcel_interp(dt,xh,uh,ruh,yh,vh,rvh,zh,mh,rmh,mf,        &
                           pi0,thv0,th0,dum1,dum2,dum3,dum4,divx,prs,  &
                           ua,va,wa,ppi,thten,tha,qa,khv,pdata,        &
                           ploc,packet,reqs_p,                         &
                           pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,            &
                           nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        call parcel_write(prec,rtime,pdata)
        prcltim = prcltim + prclfrq
        if(timestats.ge.1) time_parcels=time_parcels+mytime()
      ENDIF
      endif

!-------------------------------------------------------------------

      if(rtime.ge.rsttim .and. rstfrq.gt.0)then
       if(irst /= 0)then 
          call end_timestamp()
       endif
       print *, '#########################################',irst
       if(myid.eq.0)then
         print *, 'reading in from notify file'
         open(456,file='notify/cm1.notify',status='old')
         read(456,101)notify
         if(notify==1)then
           print *,'preparing to terminate...'
           close(456)
           call EXIT(1)
         endif
         close(456)
       endif
       101 format(i6)

        nrst=nrst+1
        call chkpt_all(myid)
        call write_restart(nstep,nrec,prec,nwrite,nrst,nrad2d,num_soil_layers, &
                               dt,mtime,radtim,qbudget,asq,bsq,                &
                               rain,sws,svs,sps,srs,sgs,sus,shs,tsk,radbcw,radbce,radbcs,radbcn,     &
                               rho,prs,ua,va,wa,ppi,tha,                           &
                               qpten,qtten,qvten,qcten,qiten,qa,tkea,swten,lwten,  &
                               radsw,rnflx,radswnet,radlwin,rad2d,                 &
                               lu_index,kpbl2d,psfc,u10,v10,hfx,qfx,xland,znt,ust, &
                               hpbl,wspd,psim,psih,gz1oz0,br,                      &
                               CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                      &
                               MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                   &
                               CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                               f2d,gsw,glw,chklowq,capg,snowc,tslb,                &
                               tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,              &
                               pta,pdata,rtime)
        rsttim=rsttim+rstfrq
        if(timestats.ge.1) time_write=time_write+mytime()
      endif

!-------------------------------------------------------------------

      if(stopit)then
        if(myid.eq.0)then
          print *
          print *,' Courant number has exceeded 1.5 '
          print *
          print *,' Stopping model .... '
          print *
        endif
        call stopcm1
      endif

!--------------------------------------------------------------------

      return
      end


