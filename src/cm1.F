

      program cm1
      use Allocator
      implicit none

!-----------------------------------------------------------------------
!  CM1 Numerical Model, Release 16  (cm1r16)
!  6 February 2012
!  http://www.mmm.ucar.edu/people/bryan/cm1/
!-----------------------------------------------------------------------
!
!  Please see documentation at the top of the "solve.F" file.
!
!  See also documentation at the cm1 website, such as:
!
!    "The governing equations for CM1"
!        http://www.mmm.ucar.edu/people/bryan/cm1/cm1_equations.pdf
!
!-----------------------------------------------------------------------

      include 'input.incl'
      include 'radcst.incl'
      include 'constants.incl'
      include 'timestat.incl'
#ifdef MPI
      include 'mpif.h'
#endif
      CHARACTER(LEN=10) varname
      integer  :: cmtsize
      integer :: notify

      integer :: nstep

      integer :: nrec,prec,nwrite,nrst

      !integer, pointer :: nrec
      !integer, pointer :: prec
      !integer, pointer :: nwrite
      !integer, pointer :: nrst

      integer :: num_soil_layers,rbufsz,ndt

      !integer, pointer :: num_soil_layers

      real :: dt,dtlast,th00s,thlr

      !real, pointer :: dt

      real*8 :: mtime,stattim,taptim,rsttim,radtim,adt,acfl

      !real*8, pointer :: mtime
      !real*8, pointer :: radtim

      logical :: dodrag,dosfcflx
      logical, dimension(maxq) :: cloudvar,rhovar
      character*15 :: tdef
      character*3, dimension(maxq) :: qname
      character*6, dimension(maxq) :: budname
      real*8, dimension(:), allocatable :: bud,bud2
      !real*8, dimension(:), allocatable :: qbudget
      !real*8, dimension(:), allocatable :: asq,bsq

      real*8, pointer :: qbudget(:)
      real*8, pointer :: asq(:)
      real*8, pointer :: bsq(:)

      real, dimension(:), allocatable :: xh,rxh,uh,ruh
      real, dimension(:), allocatable :: xf,rxf,uf,ruf
      real, dimension(:), allocatable :: yh,vh,rvh
      real, dimension(:), allocatable :: yf,vf,rvf
      real, dimension(:), allocatable :: xfref,yfref
      real, dimension(:), allocatable :: rds,sigma,rdsf,sigmaf
      real, dimension(:,:,:), allocatable :: tauh,taus,zh,mh,rmh,c1,c2
      real, dimension(:,:,:), allocatable :: tauf,zf,mf,rmf
      real, dimension(:), allocatable :: rstat
      real, dimension(:,:), allocatable :: rho0s,pi0s,prs0s,rth0s
      real, dimension(:,:,:), allocatable :: pi0,rho0,prs0,thv0,th0,th00,pi00,qv0
      real, dimension(:,:,:), allocatable :: ql0,rr0,rf0,rrf0,u0,v0
      real, dimension(:,:,:), allocatable :: t0,rh0,qc0
      real, dimension(:,:), allocatable :: zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, dimension(:,:,:), allocatable :: gx,gxu,gy,gyv
      !real, dimension(:,:,:), allocatable :: rain,sws,svs,sps,srs,sgs,sus,shs

      real, pointer :: rain(:,:,:)
      real, pointer :: sws(:,:,:)
      real, pointer :: svs(:,:,:)
      real, pointer :: sps(:,:,:)
      real, pointer :: srs(:,:,:)
      real, pointer :: sgs(:,:,:)
      real, pointer :: sus(:,:,:)
      real, pointer :: shs(:,:,:)
    

      logical, dimension(:,:), allocatable :: doimpl
      real, dimension(:,:), allocatable :: thflux,qvflux,cdu,cdv,ce,u1,v1,w1

      real, pointer :: tsk(:,:)

      !real, dimension(:,:), allocatable :: radbcw,radbce
      !real, dimension(:,:), allocatable :: radbcs,radbcn

      real, pointer :: radbcw(:,:)
      real, pointer :: radbce(:,:)
      real, pointer :: radbcs(:,:)
      real, pointer :: radbcn(:,:)

      real, dimension(:,:,:), allocatable :: dum1,dum2,dum3,dum4
      real, dimension(:,:,:), allocatable :: divx,rr,rf

      real,pointer :: rho(:,:,:)
      real,pointer :: prs(:,:,:)

      real, dimension(:,:,:), allocatable :: t11,t12,t13,t22,t23,t33
      real, dimension(:,:,:), allocatable :: rru,u3d,uten,uten1
      real, dimension(:,:,:), allocatable :: rrv,v3d,vten,vten1
      real, dimension(:,:,:), allocatable :: rrw,w3d,wten,wten1
      real, dimension(:,:,:), allocatable :: pp3d,ppten,sten
      real, dimension(:,:,:), allocatable :: th3d,thten,thten1,thterm

      real,pointer :: ua(:,:,:)
      real,pointer :: va(:,:,:)
      real,pointer :: wa(:,:,:)
      real,pointer :: ppi(:,:,:)
      real,pointer :: tha(:,:,:)
      

      !real, dimension(:,:,:), allocatable :: qpten,qtten,qvten,qcten,qiten

      real, pointer :: qpten(:,:,:)
      real, pointer :: qtten(:,:,:)
      real, pointer :: qvten(:,:,:)
      real, pointer :: qcten(:,:,:)
      real, pointer :: qiten(:,:,:)

      real, dimension(:,:,:,:), allocatable :: q3d,qten

      !real, dimension(:,:,:,:), allocatable :: qa
      real, pointer :: qa(:,:,:,:)

      real, dimension(:,:,:,:), allocatable :: zvdarray
      real, dimension(:,:,:), allocatable :: kmh,kmv,khh,khv
      real, dimension(:,:,:), allocatable :: tke3d,tketen
    
      real, pointer :: tkea(:,:,:)

      real, dimension(:,:,:), allocatable :: dissten
      real, dimension(:,:,:), allocatable :: thpten,qvpten,qcpten,qipten,upten,vpten
      real, dimension(:,:,:), allocatable :: o30

      real, pointer :: swten(:,:,:)
      real, pointer :: lwten(:,:,:)

      !real, dimension(:,:), allocatable :: radsw,rnflx,radswnet,radlwin
      !real, dimension(:,:,:), allocatable :: rad2d

      real, pointer :: radsw(:,:)
      real, pointer :: rnflx(:,:)
      real, pointer :: radswnet(:,:)
      real, pointer :: radlwin(:,:)

      real, pointer :: rad2d(:,:,:)

      real, dimension(:), allocatable :: x,y,z,za
      real, dimension(:,:,:), allocatable :: zp
      !integer, dimension(:,:), allocatable :: lu_index,kpbl2d

      integer, pointer :: lu_index(:,:)
      integer, pointer :: kpbl2d(:,:)

      !real, dimension(:,:), allocatable :: dsxy
      real, pointer :: dsxy(:,:)


      real, pointer :: psfc(:,:)
      real, pointer :: u10(:,:)
      real, pointer :: v10(:,:)
      real, pointer :: hfx(:,:)
      real, pointer :: qfx(:,:)
      real, pointer :: xland(:,:)
      real, pointer :: znt(:,:)
      real, pointer :: ust(:,:)
      real, pointer :: hpbl(:,:)
      real, pointer :: wspd(:,:)
      real, pointer :: psim(:,:)
      real, pointer :: psih(:,:)
      real, pointer :: gz1oz0(:,:)
      real, pointer :: br(:,:)
      real, pointer :: CHS(:,:)
      real, pointer :: CHS2(:,:)
      real, pointer :: CQS2(:,:)
      real, pointer :: CPMM(:,:)
      real, pointer :: ZOL(:,:)
      real, pointer :: MAVAIL(:,:)
      real, pointer :: MOL(:,:)
      real, pointer :: RMOL(:,:)
      real, pointer :: REGIME(:,:)
      real, pointer :: LH(:,:)
      real, pointer :: FLHC(:,:)
      real, pointer :: FLQC(:,:)
      real, pointer :: QGH(:,:)
      real, pointer :: CK(:,:)
      real, pointer :: CKA(:,:)
      real, pointer :: CD(:,:)
      real, pointer :: CDA(:,:)
      real, pointer :: USTM(:,:)
      real, pointer :: QSFC(:,:)
      real, pointer :: T2(:,:)
      real, pointer :: Q2(:,:)
      real, pointer :: TH2(:,:)
      real, pointer :: EMISS(:,:)
      real, pointer :: THC(:,:)
      real, pointer :: ALBD(:,:)
      real, pointer :: f2d(:,:)
      real, pointer :: gsw(:,:)
      real, pointer :: glw(:,:)
      real, pointer :: chklowq(:,:)
      real, pointer :: capg(:,:)
      real, pointer :: snowc(:,:)


      real, dimension(:), allocatable :: slab_zs,slab_dzs
      !real, dimension(:,:,:), allocatable :: tslb

      real, pointer :: tslb(:,:,:)

      !real, dimension(:,:), allocatable :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml

      real, pointer :: tmn(:,:)
      real, pointer :: tml(:,:)
      real, pointer :: t0ml(:,:)
      real, pointer :: hml(:,:)
      real, pointer :: h0ml(:,:)
      real, pointer :: huml(:,:)
      real, pointer :: hvml(:,:)
      real, pointer :: tmoml(:,:)

      real, dimension(:,:,:,:),  allocatable :: pt3d,ptten

      !real, dimension(:,:,:,:),  allocatable :: pta
      real, pointer :: pta(:,:,:,:)

      !real, dimension(:,:),      allocatable :: pdata

      real, pointer :: pdata(:,:)

      real, dimension(:,:,:),    allocatable :: cfb
      real, dimension(:),        allocatable :: cfa,cfc,d1,d2
      complex, dimension(:,:,:), allocatable :: pdt,deft
      complex, dimension(:,:),   allocatable :: rhs,trans
      logical, dimension(:,:,:), allocatable :: flag

!--- arrays for MPI ---
      integer, dimension(:), allocatable :: reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_tk
      integer, dimension(:,:),  allocatable :: reqs_q,reqs_t
      real, dimension(:), allocatable :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, dimension(:,:), allocatable :: ww1,ww2,we1,we2
      real, dimension(:,:), allocatable :: ws1,ws2,wn1,wn2
      real, dimension(:,:), allocatable :: pw1,pw2,pe1,pe2
      real, dimension(:,:), allocatable :: ps1,ps2,pn1,pn2
      real, dimension(:,:), allocatable :: vw1,vw2,ve1,ve2
      real, dimension(:,:), allocatable :: vs1,vs2,vn1,vn2
      real, dimension(:,:,:), allocatable :: uw31,uw32,ue31,ue32
      real, dimension(:,:,:), allocatable :: us31,us32,un31,un32
      real, dimension(:,:,:), allocatable :: vw31,vw32,ve31,ve32
      real, dimension(:,:,:), allocatable :: vs31,vs32,vn31,vn32
      real, dimension(:,:,:), allocatable :: ww31,ww32,we31,we32
      real, dimension(:,:,:), allocatable :: ws31,ws32,wn31,wn32
      real, dimension(:,:,:), allocatable :: sw31,sw32,se31,se32
      real, dimension(:,:,:), allocatable :: ss31,ss32,sn31,sn32
      real, dimension(:,:,:,:), allocatable :: rw31,rw32,re31,re32
      real, dimension(:,:,:,:), allocatable :: rs31,rs32,rn31,rn32
      real, dimension(:,:,:,:), allocatable :: qw31,qw32,qe31,qe32
      real, dimension(:,:,:,:), allocatable :: qs31,qs32,qn31,qn32
      real, dimension(:,:,:), allocatable :: tkw1,tkw2,tke1,tke2
      real, dimension(:,:,:), allocatable :: tks1,tks2,tkn1,tkn2
      real, dimension(:,:,:), allocatable :: kw1,kw2,ke1,ke2
      real, dimension(:,:,:), allocatable :: ks1,ks2,kn1,kn2
      real, dimension(:,:,:,:), allocatable :: tw1,tw2,te1,te2
      real, dimension(:,:,:,:), allocatable :: ts1,ts2,tn1,tn2
      real, dimension(:,:), allocatable :: ploc,packet
      real, dimension(:,:), allocatable :: dat1,dat2
      real, dimension(:,:,:), allocatable :: dat3

!-----

      integer count,rate,maxr
      real xtime,time_solve,rtime

      !real, pointer :: rtime

      real steptime1,steptime2
      integer :: i,j,k,n,nn,fnum
      real :: sum
#ifdef MPI
      integer rc
      real mp_total
      real*8 tstart,tend
#endif


      namelist /param0/ nx,ny,nz,nodex,nodey,timeformat,timestats,terrain_flag,procfiles

!----------------------------------------------------------------------

      nstep = 0
      mtime = 0.0d0
      nrec=1
      prec=1
      nwrite=1
      nrst=0
      outfile=6
      stopit = .false.
      smeps = 1.0e-30
      tsmall = 0.0001
#ifdef DP
      smeps = 1.0e-60
      tsmall = 0.000000000001
#endif

!----------------------------------------------------------------------
!  Initialize MPI

      myid=0
      numprocs=1

#ifdef MPI
      call MPI_INIT( ierr )
      call MPI_COMM_RANK( MPI_COMM_WORLD, myid, ierr )
      call MPI_COMM_SIZE( MPI_COMM_WORLD, numprocs, ierr )
#endif

      call px_init(myid)

!----------------------------------------------------------------------
!  Get domain dimensions, allocate some arrays, then call PARAM

      open(unit=20,file='namelist.input',form='formatted',status='old',    &
           access='sequential')
      read(20,nml=param0)
      close(unit=20)

      IF( procfiles )THEN
#ifdef MPI
        outfile=10
#endif
        dowr = .true.
      ELSE
        dowr = .false.
      ENDIF

      IF( myid.eq.0 ) dowr = .true.

      ni = nx / nodex
      nj = ny / nodey
      nk = nz
      nkp1 = nk+1

      ! (The following are needed by ZVD, but are also included for future 
      !  development, e.g., possible distributed-memory decomposition in 
      !  z direction)
      !
      ! number of 'ghost' points in the horizontal directions:
      ngxy  = 3
      ! number of 'ghost' points in the vertical direction:
      ngz   = 1

!---------------------------------------------------------------------
!      For ZVD:
!      ngz   = 3
!      IF( ngz.eq.3 )THEN
!        kb =  1 - ngz
!        ke = nk + ngz
!      ENDIF
!---------------------------------------------------------------------

      ib =  1 - ngxy
      ie = ni + ngxy
      jb =  1 - ngxy
      je = nj + ngxy
      kb =  1 - ngz
      ke = nk + ngz

      allocate(    xh(ib:ie) )
      allocate(   rxh(ib:ie) )
      allocate(    uh(ib:ie) )
      allocate(   ruh(ib:ie) )
      allocate(    xf(ib:ie+1) )
      allocate(   rxf(ib:ie+1) )
      allocate(    uf(ib:ie+1) )
      allocate(   ruf(ib:ie+1) )
      allocate(    yh(jb:je) )
      allocate(    vh(jb:je) )
      allocate(   rvh(jb:je) )
      allocate(    yf(jb:je+1) )
      allocate(    vf(jb:je+1) )
      allocate(   rvf(jb:je+1) )
      allocate( xfref(-2:nx+4) )
      allocate( yfref(-2:ny+4) )
      allocate(   rds(kb:ke) )
      allocate( sigma(kb:ke) )
      allocate(   rdsf(kb:ke+1) )
      allocate( sigmaf(kb:ke+1) )
      allocate(  tauh(ib:ie,jb:je,kb:ke) )
      allocate(  taus(ib:ie,jb:je,kb:ke) )
      allocate(    zh(ib:ie,jb:je,kb:ke) )
      allocate(    mh(ib:ie,jb:je,kb:ke) )
      allocate(   rmh(ib:ie,jb:je,kb:ke) )
      allocate(    c1(ib:ie,jb:je,kb:ke) )
      allocate(    c2(ib:ie,jb:je,kb:ke) )
      allocate(  tauf(ib:ie,jb:je,kb:ke+1) )
      allocate(    mf(ib:ie,jb:je,kb:ke+1) )
      allocate(   rmf(ib:ie,jb:je,kb:ke+1) )

      if(terrain_flag)then
        itb=ib
        ite=ie
        jtb=jb
        jte=je
        ktb=kb
        kte=ke
      else
        itb=1
        ite=1
        jtb=1
        jte=1
        ktb=1
        kte=1
      endif

      allocate(   zs(itb:ite,jtb:jte) )
      allocate(   gz(itb:ite,jtb:jte) )
      allocate(  rgz(itb:ite,jtb:jte) )
      allocate(  gzu(itb:ite,jtb:jte) )
      allocate( rgzu(itb:ite,jtb:jte) )
      allocate(  gzv(itb:ite,jtb:jte) )
      allocate( rgzv(itb:ite,jtb:jte) )
      allocate( dzdx(itb:ite,jtb:jte) )
      allocate( dzdy(itb:ite,jtb:jte) )
      allocate(   gx(itb:ite,jtb:jte,ktb:kte) )
      allocate(  gxu(itb:ite,jtb:jte,ktb:kte) )
      allocate(   gy(itb:ite,jtb:jte,ktb:kte) )
      allocate(  gyv(itb:ite,jtb:jte,ktb:kte) )
      allocate(   zf(ib:ie,jb:je,kb:ke+1) )

!------
! allocate the MPI arrays

#ifdef MPI
      imp = max(1,ni)
      jmp = max(1,nj)
      kmp = max(2,nk)
      kmt = max(2,nk+1)
      rmp = 8
      cmp = 3
#else
      imp = 1
      jmp = 1
      kmp = 2
      kmt = 2
      rmp = 1
      cmp = 1
#endif

      allocate( reqs_u(rmp) )
      allocate( reqs_v(rmp) )
      allocate( reqs_w(rmp) )
      allocate( reqs_s(rmp) )
      allocate( reqs_p(rmp) )
      allocate( reqs_tk(rmp) )

      allocate( nw1(kmt) )
      allocate( nw2(kmt) )
      allocate( ne1(kmt) )
      allocate( ne2(kmt) )
      allocate( sw1(kmt) )
      allocate( sw2(kmt) )
      allocate( se1(kmt) )
      allocate( se2(kmt) )

      allocate( ww1(jmp,kmp-1) )
      allocate( ww2(jmp,kmp-1) )
      allocate( we1(jmp,kmp-1) )
      allocate( we2(jmp,kmp-1) )
      allocate( ws1(imp,kmp-1) )
      allocate( ws2(imp,kmp-1) )
      allocate( wn1(imp,kmp-1) )
      allocate( wn2(imp,kmp-1) )

      allocate( pw1(jmp,kmp) )
      allocate( pw2(jmp,kmp) )
      allocate( pe1(jmp,kmp) )
      allocate( pe2(jmp,kmp) )
      allocate( ps1(imp,kmp) )
      allocate( ps2(imp,kmp) )
      allocate( pn1(imp,kmp) )
      allocate( pn2(imp,kmp) )

      allocate( vw1(jmp,kmp) )
      allocate( vw2(jmp,kmp) )
      allocate( ve1(jmp,kmp) )
      allocate( ve2(jmp,kmp) )
      allocate( vs1(imp,kmp) )
      allocate( vs2(imp,kmp) )
      allocate( vn1(imp,kmp) )
      allocate( vn2(imp,kmp) )

      allocate( uw31(cmp,jmp,kmp) )
      allocate( uw32(cmp,jmp,kmp) )
      allocate( ue31(cmp,jmp,kmp) )
      allocate( ue32(cmp,jmp,kmp) )
      allocate( us31(imp+1,cmp,kmp) )
      allocate( us32(imp+1,cmp,kmp) )
      allocate( un31(imp+1,cmp,kmp) )
      allocate( un32(imp+1,cmp,kmp) )

      allocate( vw31(cmp,jmp+1,kmp) )
      allocate( vw32(cmp,jmp+1,kmp) )
      allocate( ve31(cmp,jmp+1,kmp) )
      allocate( ve32(cmp,jmp+1,kmp) )
      allocate( vs31(imp,cmp,kmp) )
      allocate( vs32(imp,cmp,kmp) )
      allocate( vn31(imp,cmp,kmp) )
      allocate( vn32(imp,cmp,kmp) )

      allocate( ww31(cmp,jmp,kmp-1) )
      allocate( ww32(cmp,jmp,kmp-1) )
      allocate( we31(cmp,jmp,kmp-1) )
      allocate( we32(cmp,jmp,kmp-1) )
      allocate( ws31(imp,cmp,kmp-1) )
      allocate( ws32(imp,cmp,kmp-1) )
      allocate( wn31(imp,cmp,kmp-1) )
      allocate( wn32(imp,cmp,kmp-1) )

      allocate( sw31(cmp,jmp,kmp) )
      allocate( sw32(cmp,jmp,kmp) )
      allocate( se31(cmp,jmp,kmp) )
      allocate( se32(cmp,jmp,kmp) )
      allocate( ss31(imp,cmp,kmp) )
      allocate( ss32(imp,cmp,kmp) )
      allocate( sn31(imp,cmp,kmp) )
      allocate( sn32(imp,cmp,kmp) )

      allocate( rw31(cmp,jmp,kmp,2) )
      allocate( rw32(cmp,jmp,kmp,2) )
      allocate( re31(cmp,jmp,kmp,2) )
      allocate( re32(cmp,jmp,kmp,2) )
      allocate( rs31(imp,cmp,kmp,2) )
      allocate( rs32(imp,cmp,kmp,2) )
      allocate( rn31(imp,cmp,kmp,2) )
      allocate( rn32(imp,cmp,kmp,2) )

      allocate( tkw1(cmp,jmp,kmt) )
      allocate( tkw2(cmp,jmp,kmt) )
      allocate( tke1(cmp,jmp,kmt) )
      allocate( tke2(cmp,jmp,kmt) )
      allocate( tks1(imp,cmp,kmt) )
      allocate( tks2(imp,cmp,kmt) )
      allocate( tkn1(imp,cmp,kmt) )
      allocate( tkn2(imp,cmp,kmt) )

      allocate( kw1(jmp,kmt,4) )
      allocate( kw2(jmp,kmt,4) )
      allocate( ke1(jmp,kmt,4) )
      allocate( ke2(jmp,kmt,4) )
      allocate( ks1(imp,kmt,4) )
      allocate( ks2(imp,kmt,4) )
      allocate( kn1(imp,kmt,4) )
      allocate( kn2(imp,kmt,4) )

      call param(dt,dtlast,stattim,taptim,rsttim,radtim,          &
                 cloudvar,rhovar,qname,budname,                   &
                 xh,rxh,uh,ruh,xf,rxf,uf,ruf,yh,vh,rvh,yf,vf,rvf, &
                 xfref,yfref,                                     &
                 rds,sigma,rdsf,sigmaf,tauh,taus,zh,mh,rmh,c1,c2,tauf,zf,mf,rmf, &
                 zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,       &
                 reqs_u,reqs_v,reqs_s,reqs_p,                     &
                 nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                 sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,         &
                 uw31,uw32,ue31,ue32,us31,us32,un31,un32,         &
                 vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,         &
                 ww31,ww32,we31,we32,ws31,ws32,wn31,wn32)

      allocate( reqs_q(rmp,numq) )
      allocate( reqs_t(rmp,npt) )

      allocate( qw31(cmp,jmp,kmp,numq) )
      allocate( qw32(cmp,jmp,kmp,numq) )
      allocate( qe31(cmp,jmp,kmp,numq) )
      allocate( qe32(cmp,jmp,kmp,numq) )
      allocate( qs31(imp,cmp,kmp,numq) )
      allocate( qs32(imp,cmp,kmp,numq) )
      allocate( qn31(imp,cmp,kmp,numq) )
      allocate( qn32(imp,cmp,kmp,numq) )

      allocate( tw1(cmp,jmp,kmp,npt) )
      allocate( tw2(cmp,jmp,kmp,npt) )
      allocate( te1(cmp,jmp,kmp,npt) )
      allocate( te2(cmp,jmp,kmp,npt) )
      allocate( ts1(imp,cmp,kmp,npt) )
      allocate( ts2(imp,cmp,kmp,npt) )
      allocate( tn1(imp,cmp,kmp,npt) )
      allocate( tn2(imp,cmp,kmp,npt) )

!----------------------------------------------------------------------
!  allocate the base state arrays, then call BASE

      allocate( rstat(stat_out) )
      allocate( rho0s(ib:ie,jb:je) )
      allocate(  pi0s(ib:ie,jb:je) )
      allocate( prs0s(ib:ie,jb:je) )
      allocate( rth0s(ib:ie,jb:je) )
      allocate(  pi0(ib:ie,jb:je,kb:ke) )
      allocate( rho0(ib:ie,jb:je,kb:ke) )
      allocate( prs0(ib:ie,jb:je,kb:ke) )
      allocate( thv0(ib:ie,jb:je,kb:ke) )
      allocate(  th0(ib:ie,jb:je,kb:ke) )
      allocate( th00(ib:ie,jb:je,kb:ke) )
      allocate( pi00(ib:ie,jb:je,kb:ke) )
      allocate(  qv0(ib:ie,jb:je,kb:ke) )
      allocate(  ql0(ib:ie,jb:je,kb:ke) )
      allocate(  rr0(ib:ie,jb:je,kb:ke) )
      allocate(  rf0(ib:ie,jb:je,kb:ke) )
      allocate( rrf0(ib:ie,jb:je,kb:ke) )
      allocate(   u0(ib:ie+1,jb:je,kb:ke) )
      allocate(   v0(ib:ie,jb:je+1,kb:ke) )

      allocate(   t0(ib:ie,jb:je,kb:ke) )
      allocate(  rh0(ib:ie,jb:je,kb:ke) )
      allocate(  qc0(ib:ie,jb:je,kb:ke) )

      call base(th00s,thlr,zh,mh,c1,c2,zf,mf,rho0s,pi0s,prs0s,rth0s,  &
                pi0,prs0,rho0,thv0,th0,th00,pi00,t0,qv0,u0,v0,        &
                rh0,qc0,ql0,rr0,rf0,rrf0,                             &
                reqs_u,reqs_v,reqs_s,nw1,nw2,ne1,ne2,sw1,sw2,se1,se2, &
                uw31,uw32,ue31,ue32,us31,us32,un31,un32,              &
                vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,              &
                sw31,sw32,se31,se32,ss31,ss32,sn31,sn32)

!----------------------------------------------------------------------
!  Now, allocate the mother lode, then call INIT3D
       
      !allocate(   rain(ib:ie,jb:je,nrain) )
      varname = "rain"
      cmtsize = (ie-ib+1)*(je-jb+1)*nrain     
      call alloc_3d_real(rain,ie-ib+1,je-jb+1,nrain,varname,myid,cmtsize)     

      !allocate(    sws(ib:ie,jb:je,nrain) )
      varname = "sws"
      cmtsize = (ie-ib+1)*(je-jb+1)*nrain     
      call alloc_3d_real(sws,ie-ib+1,je-jb+1,nrain,varname,myid,cmtsize)     

      !allocate(    svs(ib:ie,jb:je,nrain) )
      varname = "svs"
      cmtsize = (ie-ib+1)*(je-jb+1)*nrain     
      call alloc_3d_real(svs,ie-ib+1,je-jb+1,nrain,varname,myid,cmtsize)     

      !allocate(    sps(ib:ie,jb:je,nrain) )
      varname = "sps"
      cmtsize = (ie-ib+1)*(je-jb+1)*nrain     
      call alloc_3d_real(sps,ie-ib+1,je-jb+1,nrain,varname,myid,cmtsize)     

      !allocate(    srs(ib:ie,jb:je,nrain) )
      varname = "srs"
      cmtsize = (ie-ib+1)*(je-jb+1)*nrain     
      call alloc_3d_real(srs,ie-ib+1,je-jb+1,nrain,varname,myid,cmtsize)     

      !allocate(    sgs(ib:ie,jb:je,nrain) )
      varname = "sgs"
      cmtsize = (ie-ib+1)*(je-jb+1)*nrain     
      call alloc_3d_real(sgs,ie-ib+1,je-jb+1,nrain,varname,myid,cmtsize)     

      !allocate(    sus(ib:ie,jb:je,nrain) )
      varname = "sus"
      cmtsize = (ie-ib+1)*(je-jb+1)*nrain     
      call alloc_3d_real(sus,ie-ib+1,je-jb+1,nrain,varname,myid,cmtsize)     

      allocate(    shs(ib:ie,jb:je,nrain) )
      varname = "shs"
      cmtsize = (ie-ib+1)*(je-jb+1)*nrain     
      call alloc_3d_real(shs,ie-ib+1,je-jb+1,nrain,varname,myid,cmtsize)     
      
      !allocate(    tsk(ib:ie,jb:je) )
      varname = "tsk"
      cmtsize = (ie-ib+1)*(je-jb+1) 
      call alloc_2d_real(tsk,(ie-ib+1),(je-jb+1),varname,myid,cmtsize)

      !allocate( radbcw(jb:je,kb:ke) )
      varname = "radbcw"
      cmtsize = (je-jb+1)*(ke-kb+1) 
      call alloc_2d_real(radbcw,(je-jb+1),(ke-kb+1),varname,myid,cmtsize)

      !allocate( radbce(jb:je,kb:ke) )
      varname = "radbce"
      cmtsize = (je-jb+1)*(ke-kb+1) 
      call alloc_2d_real(radbce,(je-jb+1),(ke-kb+1),varname,myid,cmtsize)

      !allocate( radbcs(ib:ie,kb:ke) )
      varname = "radbcs"
      cmtsize = (ie-ib+1)*(ke-kb+1) 
      call alloc_2d_real(radbcs,(ie-ib+1),(ke-kb+1),varname,myid,cmtsize)

      !allocate( radbcn(ib:ie,kb:ke) )
      varname = "radbcn"
      cmtsize = (ie-ib+1)*(ke-kb+1) 
      call alloc_2d_real(radbcn,(ie-ib+1),(ke-kb+1),varname,myid,cmtsize)

      !allocate(  rho(ib:ie,jb:je,kb:ke) )
      varname = "rho"
      cmtsize = (ie-ib+1)*(je-jb+1)*(ke-kb+1)     
      call alloc_3d_real(rho,ie-ib+1,je-jb+1,ke-kb+1,varname,myid,cmtsize)     


      !allocate(  prs(ib:ie,jb:je,kb:ke) )
      varname = "prs"
      cmtsize = (ie-ib+1)*(je-jb+1)*(ke-kb+1)     
      call alloc_3d_real(prs,ie-ib+1,je-jb+1,ke-kb+1,varname,myid,cmtsize)     

      !allocate(    ua(ib:ie+1,jb:je,kb:ke) )
      varname = "ua"
      cmtsize = (ie-ib+2)*(je-jb+1)*(ke-kb+1)     
      call alloc_3d_real(ua,ie-ib+2,je-jb+1,ke-kb+1,varname,myid,cmtsize)     


      !allocate(    va(ib:ie,jb:je+1,kb:ke) )
      varname = "va"
      cmtsize = (ie-ib+1)*(je-jb+2)*(ke-kb+1)     
      call alloc_3d_real(va,ie-ib+1,je-jb+2,ke-kb+1,varname,myid,cmtsize)     


      !allocate(    wa(ib:ie,jb:je,kb:ke+1) )
      varname = "wa"
      cmtsize = (ie-ib+1)*(je-jb+1)*(ke-kb+2)     
      call alloc_3d_real(wa,ie-ib+1,je-jb+1,ke-kb+2,varname,myid,cmtsize)     


      !allocate(   ppi(ib:ie,jb:je,kb:ke) )
      varname = "ppi"
      cmtsize = (ie-ib+1)*(je-jb+1)*(ke-kb+1)     
      call alloc_3d_real(ppi,ie-ib+1,je-jb+1,ke-kb+1,varname,myid,cmtsize)     


      !allocate(   tha(ib:ie,jb:je,kb:ke) )
      varname = "tha"
      cmtsize = (ie-ib+1)*(je-jb+1)*(ke-kb+1)     
      call alloc_3d_real(tha,ie-ib+1,je-jb+1,ke-kb+1,varname,myid,cmtsize)     

      cmtsize = nbudget
      varname = "qbudget"
      call alloc_1d_real8(qbudget,nbudget,varname,myid,cmtsize)

      !allocate(    asq(numq) )
      cmtsize = numq
      varname = "asq"
      call alloc_1d_real8(asq,numq,varname,myid,cmtsize)

      !allocate(    bsq(numq) )
      cmtsize = numq
      varname = "bsq"
      call alloc_1d_real8(bsq,numq,varname,myid,cmtsize)

      !allocate(  qpten(ibm:iem,jbm:jem,kbm:kem) )
      varname = "qpten"
      cmtsize = (iem-ibm+1)*(jem-jbm+1)*(kem-kbm+1)     
      call alloc_3d_real(qpten,iem-ibm+1,jem-jbm+1,kem-kbm+1,varname,myid,cmtsize)     

      !allocate(  qtten(ibm:iem,jbm:jem,kbm:kem) )
      varname = "qtten"
      cmtsize = (iem-ibm+1)*(jem-jbm+1)*(kem-kbm+1)     
      call alloc_3d_real(qtten,iem-ibm+1,jem-jbm+1,kem-kbm+1,varname,myid,cmtsize)     

      !allocate(  qvten(ibm:iem,jbm:jem,kbm:kem) )
      varname = "qvten"
      cmtsize = (iem-ibm+1)*(jem-jbm+1)*(kem-kbm+1)     
      call alloc_3d_real(qvten,iem-ibm+1,jem-jbm+1,kem-kbm+1,varname,myid,cmtsize)     


      !allocate(  qcten(ibm:iem,jbm:jem,kbm:kem) )
      varname = "qcten"
      cmtsize = (iem-ibm+1)*(jem-jbm+1)*(kem-kbm+1)     
      call alloc_3d_real(qcten,iem-ibm+1,jem-jbm+1,kem-kbm+1,varname,myid,cmtsize)     

      !allocate(  qiten(ibm:iem,jbm:jem,kbm:kem) )
      varname = "qiten"
      cmtsize = (iem-ibm+1)*(jem-jbm+1)*(kem-kbm+1)     
      call alloc_3d_real(qiten,iem-ibm+1,jem-jbm+1,kem-kbm+1,varname,myid,cmtsize)     

      !allocate(     qa(ibm:iem,jbm:jem,kbm:kem,numq) )
      !print *, ibm,iem,jbm,jem,kbm,kem,numq
      !print *,"##################", qa(1,29,7,1)
      varname = "qa"
      cmtsize = (iem-ibm+1)*(jem-jbm+1)*(kem-kbm+1)*numq     
      call alloc_4d_real(qa,iem-ibm+1,jem-jbm+1,kem-kbm+1,numq,varname,myid,cmtsize)     


      !allocate(   tkea(ibt:iet,jbt:jet,kbt:ket) )
      varname = "tkea"
      cmtsize = (iet-ibt+1)*(jet-jbt+1)*(ket-kbt+1)     
      call alloc_3d_real(tkea,iet-ibt+1,jet-jbt+1,ket-kbt+1,varname,myid,cmtsize)     

      !allocate( swten(ibr:ier,jbr:jer,kbr:ker) )
      varname = "swten"
      cmtsize = (ier-ibr+1)*(jer-jbr+1)*(ker-kbr+1)     
      call alloc_3d_real(swten,ier-ibr+1,jer-jbr+1,ker-kbr+1,varname,myid,cmtsize)     

      !allocate( lwten(ibr:ier,jbr:jer,kbr:ker) )
      varname = "lwten"
      cmtsize = (ier-ibr+1)*(jer-jbr+1)*(ker-kbr+1)     
      call alloc_3d_real(lwten,ier-ibr+1,jer-jbr+1,ker-kbr+1,varname,myid,cmtsize)     


      allocate( doimpl(ib:ie,jb:je) )
      allocate( thflux(ib:ie,jb:je) )
      allocate( qvflux(ib:ie,jb:je) )
      allocate(    cdu(ib:ie,jb:je) )
      allocate(    cdv(ib:ie,jb:je) )
      allocate(     ce(ib:ie,jb:je) )
      allocate(     u1(ib:ie,jb:je) )
      allocate(     v1(ib:ie,jb:je) )
      allocate(     w1(ib:ie,jb:je) )

      allocate( dum1(ib:ie,jb:je,kb:ke) )
      allocate( dum2(ib:ie,jb:je,kb:ke) )
      allocate( dum3(ib:ie,jb:je,kb:ke) )
      allocate( dum4(ib:ie,jb:je,kb:ke) )
      allocate( divx(ib:ie,jb:je,kb:ke) )




      allocate(   rr(ib:ie,jb:je,kb:ke) )
      allocate(   rf(ib:ie,jb:je,kb:ke) )


      allocate(  t11(ib:ie,jb:je,kb:ke) )
      allocate(  t12(ib:ie,jb:je,kb:ke) )
      allocate(  t13(ib:ie,jb:je,kb:ke) )
      allocate(  t22(ib:ie,jb:je,kb:ke) )
      allocate(  t23(ib:ie,jb:je,kb:ke) )
      allocate(  t33(ib:ie,jb:je,kb:ke) )

      allocate(   rru(ib:ie+1,jb:je,kb:ke) )



      allocate(   u3d(ib:ie+1,jb:je,kb:ke) )
      allocate(  uten(ib:ie+1,jb:je,kb:ke) )
      allocate( uten1(ib:ie+1,jb:je,kb:ke) )
      allocate(   rrv(ib:ie,jb:je+1,kb:ke) )


      allocate(   v3d(ib:ie,jb:je+1,kb:ke) )
      allocate(  vten(ib:ie,jb:je+1,kb:ke) )
      allocate( vten1(ib:ie,jb:je+1,kb:ke) )
      allocate(   rrw(ib:ie,jb:je,kb:ke+1) )


      allocate(   w3d(ib:ie,jb:je,kb:ke+1) )
      allocate(  wten(ib:ie,jb:je,kb:ke+1) )
      allocate( wten1(ib:ie,jb:je,kb:ke+1) )


      allocate(  pp3d(ib:ie,jb:je,kb:ke) )
      allocate( ppten(ib:ie,jb:je,kb:ke) )
      allocate(  sten(ib:ie,jb:je,kb:ke) )


      allocate(  th3d(ib:ie,jb:je,kb:ke) )
      allocate( thten(ib:ie,jb:je,kb:ke) )
      allocate(thten1(ib:ie,jb:je,kb:ke) )
      allocate(thterm(ib:ie,jb:je,kb:ke) )

      allocate(   bud(nk) )
      allocate(  bud2(nj) )

      !allocate( qbudget(nbudget) )
      !print *, "************", sizeof(mtime)
      



      allocate(    q3d(ibm:iem,jbm:jem,kbm:kem,numq) )
      allocate(   qten(ibm:iem,jbm:jem,kbm:kem,numq) )
      allocate( zvdarray(ibzvd:iezvd,jbzvd:jezvd,kbzvd:kezvd,nqzvd) )
      allocate(    kmh(ibc:iec,jbc:jec,kbc:kec) )
      allocate(    kmv(ibc:iec,jbc:jec,kbc:kec) )
      allocate(    khh(ibc:iec,jbc:jec,kbc:kec) )
      allocate(    khv(ibc:iec,jbc:jec,kbc:kec) )


      allocate(  tke3d(ibt:iet,jbt:jet,kbt:ket) )
      allocate( tketen(ibt:iet,jbt:jet,kbt:ket) )

      allocate( dissten(ib:ie,jb:je,kb:ke) )

      allocate( thpten(ibb:ieb,jbb:jeb,kbb:keb) )
      allocate( qvpten(ibb:ieb,jbb:jeb,kbb:keb) )
      allocate( qcpten(ibb:ieb,jbb:jeb,kbb:keb) )
      allocate( qipten(ibb:ieb,jbb:jeb,kbb:keb) )
      allocate(  upten(ibb:ieb,jbb:jeb,kbb:keb) )
      allocate(  vpten(ibb:ieb,jbb:jeb,kbb:keb) )



      allocate(   o30(ibr:ier,jbr:jer,kbr:ker) )

      nir = 1
      njr = 1
      nkr = nk+3

      IF( radopt .eq. 1 )THEN
        rbufsz = n2d_radiat*nir*njr + n3d_radiat*nir*njr*nkr
      ELSE
        rbufsz = 1
      ENDIF


      !allocate(    rad2d(ni,nj,nrad2d) )
      varname = "rad2d"
      cmtsize = ni*nj*nrad2d     
      call alloc_3d_real(rad2d,ni,nj,nrad2d,varname,myid,cmtsize)     

      !allocate(    radsw(ni,nj) )
      varname = "radsw"
      cmtsize = ni*nj     
      call alloc_2d_real(radsw,ni,nj,varname,myid,cmtsize)     

      !allocate(    rnflx(ni,nj) )
      varname = "rnflx"
      cmtsize = ni*nj     
      call alloc_2d_real(rnflx,ni,nj,varname,myid,cmtsize)     

      !allocate( radswnet(ni,nj) )
      varname = "radswnet"
      cmtsize = ni*nj     
      call alloc_2d_real(radswnet,ni,nj,varname,myid,cmtsize)     

      !allocate(  radlwin(ni,nj) )
      varname = "radlwin"
      cmtsize = ni*nj     
      call alloc_2d_real(radlwin,ni,nj,varname,myid,cmtsize)     


      rad2d = 0.0
      radsw = 0.0
      rnflx = 0.0
      radswnet = 0.0
      radlwin = 0.0

      if(dowr) write(outfile,*) '  rbufsz,nrad2d = ',rbufsz,nrad2d

      allocate( x(ni+1) )
      allocate( y(nj+1) )
      allocate( z(nk+3) )
      allocate( za(nk+3) )
      allocate( zp(ni,nj,nk+3) )


      !allocate( lu_index(ibl:iel,jbl:jel) )
      varname = "lu_index"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_integer(lu_index,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(   kpbl2d(ibl:iel,jbl:jel) )
      varname = "kpbl2d"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_integer(kpbl2d,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     


      !allocate(     psfc(ibl:iel,jbl:jel) )
      varname = "psfc"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(psfc,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      u10(ibl:iel,jbl:jel) )
      varname = "u10"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(u10,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      v10(ibl:iel,jbl:jel) )
      varname = "v10"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(v10,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      hfx(ibl:iel,jbl:jel) )
      varname = "hfx"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(hfx,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      qfx(ibl:iel,jbl:jel) )
      varname = "qfx"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(qfx,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(    xland(ibl:iel,jbl:jel) )
      varname = "xland"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(xland,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      znt(ibl:iel,jbl:jel) )
      varname = "znt"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(znt,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      ust(ibl:iel,jbl:jel) )
      varname = "ust"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(ust,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     hpbl(ibl:iel,jbl:jel) )
      varname = "hpbl"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(hpbl,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     wspd(ibl:iel,jbl:jel) )
      varname = "wspd"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(wspd,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     psim(ibl:iel,jbl:jel) )
      varname = "psim"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(psim,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     psih(ibl:iel,jbl:jel) )
      varname = "psih"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(psih,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(   gz1oz0(ibl:iel,jbl:jel) )
      varname = "gz1oz0"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(gz1oz0,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(       br(ibl:iel,jbl:jel) )
      varname = "br"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(br,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      chs(ibl:iel,jbl:jel) )
      varname = "chs"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(chs,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     chs2(ibl:iel,jbl:jel) )
      varname = "chs2"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(chs2,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     cqs2(ibl:iel,jbl:jel) )
      varname = "cqs2"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(cqs2,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     cpmm(ibl:iel,jbl:jel) )
      varname = "cpmm"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(cpmm,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      zol(ibl:iel,jbl:jel) )
      varname = "zol"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(zol,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(   mavail(ibl:iel,jbl:jel) )
      varname = "mavail"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(mavail,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      mol(ibl:iel,jbl:jel) )
      varname = "mol"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(mol,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     rmol(ibl:iel,jbl:jel) )
      varname = "rmol"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(rmol,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(   regime(ibl:iel,jbl:jel) )
      varname = "regime"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(regime,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(       lh(ibl:iel,jbl:jel) )
      varname = "lh"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(lh,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     flhc(ibl:iel,jbl:jel) )
      varname = "flhc"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(flhc,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     flqc(ibl:iel,jbl:jel) )
      varname = "flqc"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(flqc,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      qgh(ibl:iel,jbl:jel) )
      varname = "qgh"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(qgh,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(       ck(ibl:iel,jbl:jel) )
      varname = "ck"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(ck,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      cka(ibl:iel,jbl:jel) )
      varname = "cka"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(cka,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(       cd(ibl:iel,jbl:jel) )
      varname = "cd"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(cd,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      cda(ibl:iel,jbl:jel) )
      varname = "cda"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(cda,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     ustm(ibl:iel,jbl:jel) )
      varname = "ustm"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(ustm,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     qsfc(ibl:iel,jbl:jel) )
      varname = "qsfc"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(qsfc,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(       t2(ibl:iel,jbl:jel) )
      varname = "t2"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(t2,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(       q2(ibl:iel,jbl:jel) )
      varname = "q2"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(q2,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      th2(ibl:iel,jbl:jel) )
      varname = "th2"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(th2,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(    emiss(ibl:iel,jbl:jel) )
      varname = "emiss"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(emiss,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      thc(ibl:iel,jbl:jel) )
      varname = "thc"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(thc,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     albd(ibl:iel,jbl:jel) )
      varname = "albd"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(albd,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      f2d(ibl:iel,jbl:jel) )
      varname = "f2d"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(f2d,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      gsw(ibl:iel,jbl:jel) )
      varname = "gsw"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(gsw,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(      glw(ibl:iel,jbl:jel) )
      varname = "glw"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(glw,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(  chklowq(ibl:iel,jbl:jel) )
      varname = "chklowq"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(chklowq,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     capg(ibl:iel,jbl:jel) )
      varname = "capg"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(capg,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(    snowc(ibl:iel,jbl:jel) )
      varname = "snowc"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(snowc,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(     dsxy(ibl:iel,jbl:jel) )
      varname = "dsxy"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(dsxy,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     


      ! start with very small, but non-zero, numbers:
      znt = 1.0e-6
      ust = 1.0e-6

      ! start assuming neutral sfclayer:
      mol = 0.0
      zol = 0.0

      num_soil_layers = 5
      allocate(  slab_zs(num_soil_layers) )
      allocate( slab_dzs(num_soil_layers) )
      allocate(  tslb(ibl:iel,jbl:jel,num_soil_layers) )
      allocate(   tmn(ibl:iel,jbl:jel) )



      ! arrays for oml model:
      !allocate(   tml(ibl:iel,jbl:jel) )
      varname = "tml"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(tml,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(  t0ml(ibl:iel,jbl:jel) )
      varname = "t0ml"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(t0ml,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(   hml(ibl:iel,jbl:jel) )
      varname = "hml"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(hml,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(  h0ml(ibl:iel,jbl:jel) )
      varname = "dsxy"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(dsxy,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(  huml(ibl:iel,jbl:jel) )
      varname = "huml"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(huml,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(  hvml(ibl:iel,jbl:jel) )
      varname = "hvml"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(hvml,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate( tmoml(ibl:iel,jbl:jel) )
      varname = "tmoml"
      cmtsize = (iel-ibl+1)*(jel-jbl+1)     
      call alloc_2d_real(tmoml,iel-ibl+1,jel-jbl+1,varname,myid,cmtsize)     

      !allocate(    pta(ibp:iep,jbp:jep,kbp:kep,npt) )
      varname = "pta"
      cmtsize = (iep-ibp+1)*(jep-jbp+1)*(kep-kbp+1)*npt     
      call alloc_4d_real(pta,iel-ibl+1,jel-jbl+1,kep-kbp+1,npt,varname,myid,cmtsize)     

      allocate(   pt3d(ibp:iep,jbp:jep,kbp:kep,npt) )
      allocate(  ptten(ibp:iep,jbp:jep,kbp:kep,npt) )

      !allocate(  pdata(npvals,nparcels) )
      varname = "pdata"
      cmtsize = (npvals)*(nparcels)     
      call alloc_2d_real(pdata,npvals,nparcels,varname,myid,cmtsize)     


      allocate(    cfb(ipb:ipe,jpb:jpe,kpb:kpe) )
      allocate(    cfa(kpb:kpe) )
      allocate(    cfc(kpb:kpe) )
      allocate(     d1(kpb:kpe) )
      allocate(     d2(kpb:kpe) )
      allocate(    pdt(ipb:ipe,jpb:jpe,kpb:kpe) )
      allocate(   deft(ipb:ipe,jpb:jpe,kpb:kpe) )
      allocate(    rhs(ipb:ipe,jpb:jpe) )
      allocate(  trans(ipb:ipe,jpb:jpe) )

      call init3d(num_soil_layers,qbudget,asq,bsq,                  &
                  xh,rxh,uh,ruh,xf,rxf,uf,ruf,yh,vh,rvh,yf,vf,rvf,  &
                  xfref,yfref,sigma,c1,c2,gz,                       &
                  zh,mh,rmh,zf,mf,rmf,rho0s,pi0s,prs0s,             &
                  pi0,prs0,rho0,thv0,th0,t0,qv0,                    &
                  u0,v0,rh0,qc0,ql0,rr0,rf0,rrf0,                   &
                  rain,sws,svs,sps,srs,sgs,sus,shs,                 &
                  thflux,qvflux,cdu,cdv,ce,                         &
                  radbcw,radbce,radbcs,radbcn,                      &
                  dum1,dum2,dum3,dum4,divx,rho,prs,                 &
                  t11,t12,t13,t22,t23,t33,                          &
                  rru,ua,u3d,uten,uten1,rrv,va,v3d,vten,vten1,      &
                  rrw,wa,w3d,wten,wten1,ppi,pp3d,ppten,sten,        &
                  tha,th3d,thten,thten1,thterm,                     &
                  qpten,qtten,qvten,qcten,qiten,qa,q3d,qten,        &
                  kmh,kmv,khh,khv,tkea,tke3d,tketen,dissten,        &
                  pta,pt3d,ptten,                                   &
                  pdata,cfb,cfa,cfc,d1,d2,pdt,deft,rhs,trans)

!----------------------------------------------------------------------

      deallocate( t0 )
      deallocate( rh0 )
!!!      deallocate( qc0 )
      if(ibalance.eq.2 .and.  psolver.ne.4.and.psolver.ne.5 )then
        deallocate( cfb )
        deallocate( cfa )
        deallocate( cfc )
        deallocate( d1 )
        deallocate( d2 )
        deallocate( pdt )
        deallocate( deft )
        deallocate( rhs )
        deallocate( trans )
      endif

      allocate( flag(ib:ie,jb:je,kb:ke) )

      allocate(          ploc(3,nparcels) )
      allocate( packet(npvals+1,nparcels) )

      allocate( dat1(ni+1,nj+1) )
      allocate( dat2(nx+1,ny+1) )
      allocate( dat3(ni+1,nj+1,numprocs) )

!----------------------------------------------------------------------

      call setup_output(tdef,qname,budname,xh,xf,yh,yf,xfref,yfref,sigma,sigmaf,zh,zf)

      call init_physics(prs0,rf0,cdu,cdv,ce,dum1,dum2,dum3,u0,ua,v0,va,o30,   &
                             lu_index,xland,emiss,thc,albd,znt,mavail,f2d,tsk,u1,v1,w1, &
                             zh,u10,v10,wspd,cd)

      call init_surface(num_soil_layers,   &
                        dodrag,dosfcflx,xh,ruh,xf,yh,rvh,yf,   &
                        lu_index,xland,tsk,slab_zs,slab_dzs,tslb, &
                        emiss,thc,albd,znt,mavail,dsxy,prs0s,prs0,   &
                        tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml)

      if( (sfcmodel.eq.1).or.(idrag.eq.1) )then
        call getcecd(cdu,cdv,ce,u0,v0,u1,v1,w1,ua,va,zh,u10,v10,wspd,xland,znt,ust,cd,  &
                     ww31(1,1,1),ww32(1,1,1),we31(1,1,1),we32(1,1,1),                   &
                     ws31(1,1,1),ws32(1,1,1),wn31(1,1,1),wn32(1,1,1),reqs_s)
      endif
      print *, "End of application init. starting restart calc"

      if(irst.eq.1)then
        call read_restart(nstep,nrec,prec,nwrite,nrst,nrad2d,num_soil_layers, &
                              stattim,taptim,rsttim,  &
                              dt,mtime,radtim,qbudget,asq,bsq,              &
                              rain,sws,svs,sps,srs,sgs,sus,shs,tsk,radbcw,radbce,radbcs,radbcn,  &
                              rho,prs,ua,va,wa,ppi,tha,                            &
                              qpten,qtten,qvten,qcten,qiten,qa,tkea,swten,lwten,   &
                              radsw,rnflx,radswnet,radlwin,rad2d,   &
                               lu_index,kpbl2d,psfc,u10,v10,hfx,qfx,xland,znt,ust, &
                               hpbl,wspd,psim,psih,gz1oz0,br,                      &
                               CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                      &
                               MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                   &
                               CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                               f2d,gsw,glw,chklowq,capg,snowc,tslb,                &
                               tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,              &
                              pta,pdata,rtime)
        dtlast = 0.0
      endif

      call getset(sigma,sigmaf,dzdx,dzdy,pi0,th0,rho0,prs0,           &
                  zh,c1,c2,zf,rr,rf,rho,prs,dum1,dum2,                &
                  ua,u3d,va,v3d,wa,w3d,ppi,pp3d,                      &
                  tha,th3d,qa,q3d,tkea,tke3d,pta,pt3d,                &
                  reqs_u,reqs_v,reqs_w,reqs_s,reqs_p,reqs_tk,         &
                  nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                    &
                  pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,                    &
                  uw31,uw32,ue31,ue32,us31,us32,un31,un32,            &
                  vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,            &
                  ww31,ww32,we31,we32,ws31,ws32,wn31,wn32,            &
                  sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,            &
                  tkw1,tkw2,tke1,tke2,tks1,tks2,tkn1,tkn2)

!----------------------------------------------------------------------
!  All done with initialization.  A few more odds and ends ....

      if( adapt_dt.eq.1 )then
        call calccfl(1,rstat,dt,acfl,uf,vf,mf,ua,va,wa,0)
        ndt = 1
        adt = dt
        acfl = cflmax
      endif

      if(irst.ne.1)then
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  initial conditions:'
        if(dowr) write(outfile,*)
      endif

      IF(axisymm.eq.0)THEN
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
      rtime=sngl(mtime)
      call statpack(nrec,ndt,dt,rtime,adt,acfl,cloudvar,qname,budname,qbudget,asq,bsq, &
                    xh,rxh,uh,ruh,xf,uf,yh,vh,rvh,vf,zh,mh,rmh,mf,     &
                    zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf,                &
                    rstat,pi0,rho0,thv0,th0,qv0,u0,v0,                 &
                    dum1,dum2,dum3,dum4,divx,ppten,prs,                &
                    ua,va,wa,ppi,tha,qa,qten,kmh,kmv,khh,khv,tkea,pta,u10,v10)

    if(irst.ne.1)then
!-------------------------------------------------------------------
      ! Initial conditions:  (do not write if this is a restart)
      IF(output_format.eq.1.or.output_format.eq.2)THEN
        sten = 0.0
        nn = 1
        if(terrain_flag .and. output_interp.eq.1) nn = 2
        if(output_format.eq.2) nn = 1
        DO n=1,nn
          if(n.eq.1)then
            fnum = 51
          else
            fnum = 71
          endif
!          call writeout(fnum,1,qname,xh,xf,uf,yh,yf,vf,xfref,yfref,                            &
!                        rds,sigma,rdsf,sigmaf,zh,zf,mf,pi0,prs0,rho0,th0,qv0,u0,v0,            &
!                        zs,rgzu,rgzv,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,           &
!                        cdu,cdv,ce,dum1,dum2,dum3,dum4,                                        &
!                        t11,t12,rho,prs,sten,ua,uten,va,vten,wa,wten,ppi,tha,                  &
!                        dissten,thpten,qvpten,qcpten,qipten,upten,vpten,                       &
!                        lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,  &
!                        qa,kmh,kmv,khh,khv,tkea,swten,lwten,radsw,rnflx,radswnet,radlwin,pta,  &
!                        num_soil_layers,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br,dat1,dat2,dat3)
        ENDDO
#ifdef HDFOUT
      ELSEIF(output_format.ge.3)THEN
       sten = 0.0
        call writeout_mult_hdf5(0.0,qname,rds,sigma,rdsf,sigmaf,xh,xf,uf,yh,yf,vf,mh,zh,mf,zf, &
                      pi0,prs0,rho0,th0,qv0,u0,v0,                     &
                      zs,rgzu,rgzv,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,cdu,cdv,ce,dum1,dum2,dum3,dum4,  &
                      t11,t12,rho,prs,sten,ua,uten,va,vten,wa,wten,ppi,tha,    &
                      qa,kmh,kmv,khh,khv,tkea,pta,num_soil_layers,   &
                      lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,   &
                      radsw,rnflx,radswnet,radlwin,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br,   &
                      dissten,thpten,qvpten,qcpten,qipten,upten,vpten,swten,lwten)
#endif
      ENDIF
!-------------------------------------------------------------------
!  Write parcel data:
      if(iprcl.eq.1)then
        call parcel_interp(dt,xh,uh,ruh,yh,vh,rvh,zh,mh,rmh,mf,        &
                           pi0,thv0,th0,dum1,dum2,dum3,dum4,divx,prs,  &
                           ua,va,wa,ppi,thten,tha,qa,khv,pdata,        &
                           ploc,packet,reqs_p,                         &
                           pw1,pw2,pe1,pe2,ps1,ps2,pn1,pn2,            &
                           nw1,nw2,ne1,ne2,sw1,sw2,se1,se2)
        call parcel_write(prec,rtime,pdata)
      endif
!-------------------------------------------------------------------
    endif  ! endif for irst.ne.1

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

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '-------------Done with Preprocessors-----------'
      if(dowr) write(outfile,*)

      if(iconly.eq.1)then
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  User has requested initial conditions only'
        if(dowr) write(outfile,*) '     (iconly = 1)'
        if(dowr) write(outfile,*) '  ... stopping ... '
        if(dowr) write(outfile,*)
#ifdef MPI
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
        stop 55555
      endif

!----------------------------------------------------------------------

      time_sound=0.
      time_poiss=0.
      time_advs=0.
      time_advu=0.
      time_advv=0.
      time_advw=0.
      time_buoyan=0.
      time_turb=0.
      time_diffu=0.
      time_microphy=0.
      time_stat=0.
      time_bc=0.
      time_misc=0.
      time_integ=0.
      time_rdamp=0.
      time_divx=0.
      time_write=0.
      time_tmix=0.
      time_cor=0.
      time_fall=0.
      time_satadj=0.
      time_sfcphys=0.
      time_parcels=0.0
      time_rad=0.
      time_pbl=0.
      time_swath=0.
      time_pdef=0.
      time_prsrho=0.
#ifdef MPI
      time_mpu1=0.
      time_mpv1=0.
      time_mpw1=0.
      time_mpp1=0.
      time_mpu2=0.
      time_mpv2=0.
      time_mpw2=0.
      time_mpp2=0.
      time_mps1=0.
      time_mpq1=0.
      time_mptk1=0.
      time_mptk2=0.
      time_mps2=0.
      time_mpq2=0.
      time_mpb=0.

      call MPI_BARRIER (MPI_COMM_WORLD,ierr)

      if(myid.eq.0)then
        tstart=mpi_wtime()
      endif
#endif

      ! This initializes timer
      if(timestats.ge.1)then
        call system_clock(count,rate,maxr)
        clock_rate=1.0/rate
        xtime=mytime()
      endif
      
      print *, "End of restart calc. Starting iterations" 

      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
     if(irst /= 0) then 
      if(myid.eq.0)then
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
     endif

!----------------------------------------------------------------------
!  Time loop
!  Iterations

      if(timestats.ge.1)then
        steptime1 = 0.0
        steptime2 = 0.0
      endif

      do while( mtime.lt.timax )
        nstep = nstep + 1
        call solve(nstep,nrec,prec,nwrite,nrst,rbufsz,num_soil_layers,ndt,     &
                   dt,dtlast,th00s,thlr,mtime,stattim,taptim,rsttim,radtim,adt,acfl,  &
                   dodrag,dosfcflx,cloudvar,rhovar,qname,budname,bud,bud2,qbudget,asq,bsq, &
                   xh,rxh,uh,ruh,xf,rxf,uf,ruf,yh,vh,rvh,yf,vf,rvf,   &
                   xfref,yfref,rds,sigma,rdsf,sigmaf,tauh,taus,zh,mh,rmh,c1,c2,tauf,zf,mf,rmf,   &
                   rstat,rho0s,pi0s,prs0s,rth0s,pi0,rho0,prs0,thv0,th0,th00,pi00,qv0,qc0,  &
                   ql0,rr0,rf0,rrf0,                                  &
                   zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv, &
                   rain,sws,svs,sps,srs,sgs,sus,shs,                  &
                   doimpl,tsk,thflux,qvflux,cdu,cdv,ce,u1,v1,w1,      &
                   radbcw,radbce,radbcs,radbcn,                       &
                   dum1,dum2,dum3,dum4,divx,rho,rr,rf,prs,            &
                   t11,t12,t13,t22,t23,t33,                           &
                   u0,rru,ua,u3d,uten,uten1,                          &
                   v0,rrv,va,v3d,vten,vten1,                          &
                   rrw,wa,w3d,wten,wten1,ppi,pp3d,ppten,sten,         &
                   tha,th3d,thten,thten1,thterm,                      &
                   qpten,qtten,qvten,qcten,qiten,qa,q3d,qten,zvdarray, &
                   kmh,kmv,khh,khv,tkea,tke3d,tketen,                 &
                   dissten,thpten,qvpten,qcpten,qipten,upten,vpten,   &
                   swten,lwten,o30,radsw,rnflx,radswnet,radlwin,rad2d, &
                   x,y,z,za,zp,                                       &
                   lu_index,kpbl2d,psfc,u10,v10,hfx,qfx,xland,znt,ust,   &
                   hpbl,wspd,psim,psih,gz1oz0,br,                     &
                   CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                     &
                   MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                  &
                   CK,CKA,CD,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,  &
                   f2d,gsw,glw,chklowq,capg,snowc,dsxy,               &
                   slab_zs,slab_dzs,tslb,tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,        &
                   pta,pt3d,ptten,                                    &
                   pdata,cfb,cfa,cfc,d1,d2,pdt,deft,rhs,trans,flag,   &
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
        if(timestats.eq.2)then
          steptime2=time_sound+time_poiss+time_buoyan+time_turb+            &
                    time_diffu+time_microphy+time_stat+                     &
                    time_bc+time_misc+time_integ+time_rdamp+time_divx+      &
                    time_write+time_tmix+time_cor+time_fall+                &
                    time_satadj+time_sfcphys+time_parcels+                  &
                    time_rad+time_pbl+time_swath+time_pdef+time_prsrho+     &
#ifdef MPI
                    time_mpu1+time_mpv1+time_mpw1+time_mpp1+                &
                    time_mpu2+time_mpv2+time_mpw2+time_mpp2+                &
                    time_mps1+time_mpq1+time_mptk1+                         &
                    time_mps2+time_mpq2+time_mptk2+time_mpb+                &
#endif
                    time_advs+time_advu+time_advv+time_advw
          write(6,157) nstep,steptime2-steptime1
157       format('    timing for time step ',i12,':',f12.4,' s')
          steptime1 = steptime2
        endif
      enddo

!----------------------------------------------------------------------
!  write new stats descriptor file, if necessary:

      IF( output_format.eq.1 .and. myid.eq.0 )THEN
        IF( adapt_dt.eq.1 .and. statfrq.lt.0.0 )THEN
          print *,'  re-writing GrADS stats descriptor file .... '
          call write_statsctl(tdef,qname,budname,nstep+1)
        ENDIF
      ENDIF

!----------------------------------------------------------------------

#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      if(timestats.ge.1) time_mpb=time_mpb+mytime()

      if(myid.eq.0.and.procfiles)then
        tend=mpi_wtime()
        print *
        print *,'Total time (s): ',tend-tstart
        print *
      endif

#endif
!----------------------------------------------------------------------

    IF(timestats.ge.1)THEN

#ifdef MPI
      ! for MPI runs without procfiles, average the timestat terms:
      IF(.not.procfiles)THEN
        sum = 0.0
        call MPI_REDUCE(time_sound   ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_sound = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_poiss   ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_poiss = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_advs    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_advs = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_advu    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_advu = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_advv    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_advv = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_advw    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_advw = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_divx    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_divx = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_buoyan  ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_buoyan = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_turb    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_turb = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_sfcphys ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_sfcphys = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_tmix    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_tmix = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_cor     ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_cor = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_diffu   ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_diffu = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_rdamp   ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_rdamp = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_microphy,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_microphy = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_satadj  ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_satadj = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_fall    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_fall = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_rad     ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_rad = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_pbl     ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_pbl = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_stat    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_stat = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_bc      ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_bc = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_integ   ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_integ = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_write   ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_write = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_misc    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_misc = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_swath   ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_swath = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_pdef    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_pdef = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_prsrho  ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_prsrho = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_parcels ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_parcels = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpu1    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpu1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpv1    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpv1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpw1    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpw1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpp1    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpp1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpu2    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpu2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpv2    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpv2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpw2    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpw2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpp2    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpp2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mps1    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mps1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpq1    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpq1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mps2    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mps2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpq2    ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpq2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mptk1   ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mptk1 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mptk2   ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mptk2 = sum/float(numprocs)
        sum = 0.0
        call MPI_REDUCE(time_mpb     ,sum,1,MPI_REAL,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        time_mpb = sum/float(numprocs)
      ENDIF
#endif

      time_solve=time_sound+time_poiss+time_buoyan+time_turb+             &
                  time_diffu+time_microphy+time_stat+                     &
                  time_bc+time_misc+time_integ+time_rdamp+time_divx+      &
                  time_write+time_tmix+time_cor+time_fall+                &
                  time_satadj+time_sfcphys+time_parcels+                  &
                  time_rad+time_pbl+time_swath+time_pdef+time_prsrho+     &
#ifdef MPI
                  time_mpu1+time_mpv1+time_mpw1+time_mpp1+                &
                  time_mpu2+time_mpv2+time_mpw2+time_mpp2+                &
                  time_mps1+time_mpq1+time_mptk1+                         &
                  time_mps2+time_mpq2+time_mptk2+time_mpb+                &
#endif
                  time_advs+time_advu+time_advv+time_advw

#ifdef MPI
      mp_total=time_mpu1+time_mpv1+time_mpw1+time_mpp1+                   &
               time_mpu2+time_mpv2+time_mpw2+time_mpp2+                   &
               time_mps1+time_mpq1+time_mptk1+                            &
               time_mps2+time_mpq2+time_mptk2+time_mpb
#endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'Total time: ',time_solve
      if(dowr) write(outfile,*)
      time_solve=0.01*time_solve
      if(time_solve.lt.0.0001) time_solve=1.

    IF(dowr)THEN
      write(outfile,100) 'sound   ',time_sound,time_sound/time_solve
      write(outfile,100) 'poiss   ',time_poiss,time_poiss/time_solve
      write(outfile,100) 'advs    ',time_advs,time_advs/time_solve
      write(outfile,100) 'advu    ',time_advu,time_advu/time_solve
      write(outfile,100) 'advv    ',time_advv,time_advv/time_solve
      write(outfile,100) 'advw    ',time_advw,time_advw/time_solve
      write(outfile,100) 'divx    ',time_divx,time_divx/time_solve
      write(outfile,100) 'buoyan  ',time_buoyan,time_buoyan/time_solve
      write(outfile,100) 'turb    ',time_turb,time_turb/time_solve
      write(outfile,100) 'sfcphys ',time_sfcphys,time_sfcphys/time_solve
      write(outfile,100) 'tmix    ',time_tmix,time_tmix/time_solve
      write(outfile,100) 'cor     ',time_cor,time_cor/time_solve
      write(outfile,100) 'diffu   ',time_diffu,time_diffu/time_solve
      write(outfile,100) 'rdamp   ',time_rdamp,time_rdamp/time_solve
      write(outfile,100) 'microphy',time_microphy,time_microphy/time_solve
      write(outfile,100) 'satadj  ',time_satadj,time_satadj/time_solve
      write(outfile,100) 'fallout ',time_fall,time_fall/time_solve
      write(outfile,100) 'radiatio',time_rad,time_rad/time_solve
      write(outfile,100) 'pbl     ',time_pbl,time_pbl/time_solve
      write(outfile,100) 'stat    ',time_stat,time_stat/time_solve
      write(outfile,100) 'bc      ',time_bc,time_bc/time_solve
      write(outfile,100) 'integ   ',time_integ,time_integ/time_solve
      write(outfile,100) 'write   ',time_write,time_write/time_solve
      write(outfile,100) 'misc    ',time_misc,time_misc/time_solve
      write(outfile,100) 'swaths  ',time_swath,time_swath/time_solve
      write(outfile,100) 'pdef    ',time_pdef,time_pdef/time_solve
      write(outfile,100) 'prsrho  ',time_prsrho,time_prsrho/time_solve
      write(outfile,100) 'parcels ',time_parcels,time_parcels/time_solve
#ifdef MPI
      write(outfile,100) 'mp_total',mp_total,mp_total/time_solve
      write(outfile,*)
      write(outfile,100) 'mpu1    ',time_mpu1,time_mpu1/time_solve
      write(outfile,100) 'mpv1    ',time_mpv1,time_mpv1/time_solve
      write(outfile,100) 'mpw1    ',time_mpw1,time_mpw1/time_solve
      write(outfile,100) 'mpp1    ',time_mpp1,time_mpp1/time_solve
      write(outfile,100) 'mpu2    ',time_mpu2,time_mpu2/time_solve
      write(outfile,100) 'mpv2    ',time_mpv2,time_mpv2/time_solve
      write(outfile,100) 'mpw2    ',time_mpw2,time_mpw2/time_solve
      write(outfile,100) 'mpp2    ',time_mpp2,time_mpp2/time_solve
      write(outfile,100) 'mps1    ',time_mps1,time_mps1/time_solve
      write(outfile,100) 'mpq1    ',time_mpq1,time_mpq1/time_solve
      write(outfile,100) 'mps2    ',time_mps2,time_mps2/time_solve
      write(outfile,100) 'mpq2    ',time_mpq2,time_mpq2/time_solve
      write(outfile,100) 'mptk1   ',time_mptk1,time_mptk1/time_solve
      write(outfile,100) 'mptk2   ',time_mptk2,time_mptk2/time_solve
      write(outfile,100) 'mpb     ',time_mpb,time_mpb/time_solve
#endif
      write(outfile,*)
    ENDIF

100   format(3x,a8,' :  ',f10.2,2x,f6.2,'%')

    ENDIF

!  End time loop
!----------------------------------------------------------------------

      close(unit=51)
      close(unit=52)
      close(unit=53)
      close(unit=54)
      close(unit=60)

!----------------------------------------------------------------------

#ifdef MPI
      call MPI_FINALIZE(rc)
      print *,'Program terminated normally:  myid=',myid
#else
      print *,'Program terminated normally'
#endif

      stop
      end


