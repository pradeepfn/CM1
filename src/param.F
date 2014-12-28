

      subroutine param(dt,dtlast,stattim,taptim,rsttim,radtim,          &
                       cloudvar,rhovar,qname,budname,                   &
                       xh,rxh,uh,ruh,xf,rxf,uf,ruf,yh,vh,rvh,yf,vf,rvf, &
                       xfref,yfref,                                     &
                       rds,sigma,rdsf,sigmaf,tauh,taus,zh,mh,rmh,cc1,cc2,tauf,zf,mf,rmf, &
                       zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy,gx,gxu,gy,gyv,         &
                       reqs_u,reqs_v,reqs_s,reqs_p,                     &
                       nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,                 &
                       sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,         &
                       uw31,uw32,ue31,ue32,us31,us32,un31,un32,         &
                       vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,         &
                       ww31,ww32,we31,we32,ws31,ws32,wn31,wn32)
      use module_mp_thompson
      use module_mp_graupel
      implicit none

      include 'input.incl'
      include 'constants.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      real :: dt,dtlast
      real*8 :: stattim,taptim,rsttim,radtim
      logical, dimension(maxq) :: cloudvar,rhovar
      character*3, dimension(maxq) :: qname
      character*6, dimension(maxq) :: budname
      real, dimension(ib:ie) :: xh,rxh,uh,ruh
      real, dimension(ib:ie+1) :: xf,rxf,uf,ruf
      real, dimension(jb:je) :: yh,vh,rvh
      real, dimension(jb:je+1) :: yf,vf,rvf
      real, dimension(-2:nx+4) :: xfref
      real, dimension(-2:ny+4) :: yfref
      real, dimension(kb:ke) :: rds,sigma
      real, dimension(kb:ke+1) :: rdsf,sigmaf
      real, dimension(ib:ie,jb:je,kb:ke) :: tauh,taus,zh,mh,rmh,cc1,cc2
      real, dimension(ib:ie,jb:je,kb:ke+1) :: tauf,zf,mf,rmf
      real, dimension(itb:ite,jtb:jte) :: zs,gz,rgz,gzu,rgzu,gzv,rgzv,dzdx,dzdy
      real, dimension(itb:ite,jtb:jte,ktb:kte) :: gx,gxu,gy,gyv
      integer, intent(inout), dimension(rmp) :: reqs_u,reqs_v,reqs_s,reqs_p
      real, intent(inout), dimension(kmt) :: nw1,nw2,ne1,ne2,sw1,sw2,se1,se2
      real, intent(inout), dimension(cmp,jmp,kmp)   :: sw31,sw32,se31,se32
      real, intent(inout), dimension(imp,cmp,kmp)   :: ss31,ss32,sn31,sn32
      real, intent(inout), dimension(cmp,jmp,kmp)   :: uw31,uw32,ue31,ue32
      real, intent(inout), dimension(imp+1,cmp,kmp) :: us31,us32,un31,un32
      real, intent(inout), dimension(cmp,jmp+1,kmp) :: vw31,vw32,ve31,ve32
      real, intent(inout), dimension(imp,cmp,kmp)   :: vs31,vs32,vn31,vn32
      real, intent(inout), dimension(cmp,jmp,kmp-1) :: ww31,ww32,we31,we32
      real, intent(inout), dimension(imp,cmp,kmp-1) :: ws31,ws32,wn31,wn32

!-----------------------------------------------------------------------

      integer i,j,k,n,kst,ni1,ni2,ni3,nj1,nj2,nj3,nk1,nk2,nk3
      integer ival,jval
      integer iterrain
      integer :: inum
      real :: var
#ifdef MPI
      character*50 fname
      integer nabor
#endif
      real c1,c2,nominal_dx,nominal_dy,nominal_dz,z1,z2,z3,mult
      real x1,x2,y1,y2

      namelist /param1/ dx,dy,dz,dtl,timax,tapfrq,rstfrq,statfrq,prclfrq
      namelist /param2/                                                 &
          adapt_dt,irst,rstnum,iconly,hadvorder,vadvorder,pdscheme,     &
          advweno,idiff,vdiff,mdiff,difforder,imoist,iturb,             &
          tconfig,bcturbu,bcturbs,dns,                                  &
          irdamp,hrdamp,psolver,nsound,ptype,ihail,iautoc,              &
          icor,pertcor,neweqts,idiss,efall,rterm,                       &
          wbc,ebc,sbc,nbc,irbc,roflux,isnd,iwnd,itern,iinit,irandp,     &
          ibalance,iorigin,axisymm,imove,iptra,npt,iprcl,nparcels
      namelist /param3/ kdiff2,kdiff6,fcor,kdiv,alph,rdalpha,zd,xhd,    &
                        umove,vmove,v_t,l_h,l_inf,ndcnst
      namelist /param4/ stretch_x,dx_inner,dx_outer,nos_x_len,tot_x_len
      namelist /param5/ stretch_y,dy_inner,dy_outer,nos_y_len,tot_y_len
      namelist /param6/ stretch_z,ztop,str_bot,str_top,dz_bot,dz_top
      namelist /param7/ bc_wind,bc_temp,ptc_top,ptc_bot,viscosity,pr_num
      namelist /param8/ var1,var2,var3,var4,var5,var6,var7,var8,var9,var10
      namelist /param9/                                                       &
              output_path,output_basename,output_format,output_filetype,      &
              output_interp,output_rain,output_sws,output_coldpool,           &
              output_sfcflx,output_sfcparams,output_sfcdiags,output_zs,       &
              output_zh,output_basestate,                                     &
              output_th,output_thpert,output_prs,output_prspert,              &
              output_pi,output_pipert,output_rho,output_rhopert,output_tke,   &
              output_km,output_kh,                                            &
              output_qv,output_qvpert,output_q,output_dbz,                    &
              output_u,output_upert,output_uinterp,                           &
              output_v,output_vpert,output_vinterp,output_w,output_winterp,   &
              output_vort,output_uh,output_pblten,output_dissten,             &
              output_radten
      namelist /param10/                                                      &
              stat_w,stat_u,stat_v,stat_rmw,stat_pipert,stat_prspert,         &
              stat_thpert,stat_q,                                             &
              stat_tke,stat_km,stat_kh,stat_div,stat_rh,stat_rhi,stat_the,    &
              stat_cloud,stat_sfcprs,stat_wsp,stat_cfl,stat_vort,             &
              stat_tmass,stat_tmois,stat_qmass,stat_tenerg,stat_mo,stat_tmf,  &
              stat_pcn,stat_qsrc
      namelist /param11/                                                      &
              radopt,dtrad,ctrlat,ctrlon,year,month,day,hour,minute,second
      namelist /param12/                                                      &
              idrag,isfcflx,sfcmodel,oceanmodel,ipbl,initsfc,                 &
              tsk0,tmn0,xland0,lu0,season,cecd,pertflx,cnstce,cnstcd,         &
              isftcflx,iz0tlnd,oml_hml0,oml_gamma

!--------------------------------------------------------------

#ifdef MPI
      IF(procfiles)THEN
        fname='procXXXXXX.print.out'
        write(fname(5:10),100) myid
100     format(i6.6)
        open(unit=10,file=fname,status='unknown')
      ENDIF
#endif

      if(dowr) write(outfile,*) 'Inside PARAM'

#ifndef MPI
!--------------------------------------------------------------

      if(nodex.ne.1 .or. nodey.ne.1)then
        print *
        print *,'  For non-MPI runs, nodex and nodey must be = 1 !'
        print *
        call stopcm1
      endif

#endif
!--------------------------------------------------------------

#ifdef MPI
      if(myid.eq.0)then
#endif

      open(unit=20,file='namelist.input',form='formatted',status='old',    &
           access='sequential')
      read(20,nml=param1)
      read(20,nml=param2)
      read(20,nml=param3)
      read(20,nml=param11)
      read(20,nml=param12)
      read(20,nml=param4)
      read(20,nml=param5)
      read(20,nml=param6)
      read(20,nml=param7)
      read(20,nml=param8)
      read(20,nml=param9)
      read(20,nml=param10)
!!!      IF ( ptype .ge. 26 ) THEN
!!!         read(20,nml=micro_params)
!!!      ENDIF
      close(unit=20)

!-----------------------------------------------------------------------
!  Some dummy checks:

      if(imoist.ne.1) neweqts=0
      if(imoist.ne.1) efall=0

      IF( psolver.eq.1 .and. adapt_dt.eq.1 )THEN
        print *
        print *,'  psolver  = ',psolver
        print *,'  adapt_dt = ',adapt_dt
        print *
        print *,'  Cannot use adapt_dt with psolver=1 '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(dns.gt.1.or.dns.lt.0)THEN
        print *
        print *,'  dns   = ',dns
        print *
        print *,'  dns must be either 0 or 1'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(cecd.lt.1.or.cecd.gt.3)THEN
        print *
        print *,'  cecd  = ',cecd
        print *
        print *,'  cecd must be 1,2,3'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(iturb.gt.3.or.iturb.lt.0)THEN
        print *
        print *,'  iturb   = ',iturb
        print *
        print *,'  iturb must be either 0, 1, 2, or 3'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(iturb.ge.1 .and. dns.ge.1)THEN
        print *
        print *,'  iturb = ',iturb
        print *,'  dns   = ',dns
        print *
        print *,'  For dns = 1, iturb must be 0'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(bcturbu.lt.1.or.bcturbu.gt.3)THEN
        print *
        print *,'  bcturbu = ',bcturbu
        print *
        print *,'  bcturbu must be 1, 2, or 3'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(bcturbs.lt.1.or.bcturbs.gt.3)THEN
        print *
        print *,'  bcturbs = ',bcturbs
        print *
        print *,'  bcturbs must be 1, 2, or 3'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(dns.ge.1 .and. imoist.ge.1)THEN
        print *
        print *,'  imoist = ',imoist
        print *,'  dns    = ',dns
        print *
        print *,'  For dns = 1, imoist must be 0'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(imoist.eq.1 .and. (isnd.eq.1.or.isnd.eq.2   &
                        .or.isnd.eq.3.or.isnd.eq.8.) )THEN
        print *
        print *,'  imoist = ',imoist
        print *,'  isnd   = ',isnd
        print *
        print *,'  For this value of isnd, imoist must be 0'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(dns.eq.1 .and. (bc_wind.le.0 .or. bc_wind.ge.3))THEN
        print *
        print *,'  dns     = ',dns
        print *,'  bc_wind = ',bc_wind
        print *
        print *,'  for dns = 1, bc_wind must be either 1 or 2'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(dns.eq.1 .and. (bc_temp.le.0 .or. bc_temp.ge.3))THEN
        print *
        print *,'  dns     = ',dns
        print *,'  bc_temp = ',bc_temp
        print *
        print *,'  for dns = 1, bc_temp must be either 1 or 2'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(ihail.lt.0.or.ihail.gt.1)THEN
        print *
        print *,'  ihail   = ',ihail
        print *
        print *,'  ihail must be 0 or 1'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(imoist.eq.1.and.output_dbz.eq.1.and.ptype.ne.2.and.ptype.ne.3.and.ptype.ne.5)then
        print *
        print *,'  ptype      = ',ptype
        print *,'  output_dbz = ',output_dbz
        print *
        print *,'  output_dbz is only available for ptype=2,3,5'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(imoist.eq.1 .and. neweqts.ge.1 .and. ptype.eq.4)THEN
        print *
        print *,'  neweqts = ',neweqts
        print *,'  ptype   = ',ptype
        print *
        print *,'  neweqts >= 1 is not available for ptype = 4'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(imoist.eq.1 .and. efall.eq.1)THEN
      IF(ptype.ne.1.and.ptype.ne.2.and.ptype.ne.6)THEN
        print *
        print *,'  efall   = ',efall
        print *,'  ptype   = ',ptype
        print *
        print *,'  efall = 1 is only supported with ptype = 1,2,6'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      ENDIF
      IF((imoist.eq.1).and.(ptype.eq.4).and.terrain_flag)THEN
        print *
        print *,'  ptype   = ',ptype
        print *,'  terrain_flag = ',terrain_flag
        print *
        print *,'  ptype = 4 does not work with terrain '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(imoist.eq.1 .and. ndcnst.le.1.0e-6)THEN
        print *
        print *,'  imoist  = ',imoist
        print *,'  ndcnst  = ',ndcnst
        print *
        print *,'  ndcnst is too small.  Please enter a larger value '
        print *,'  in the param3 section of namelist.input '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(terrain_flag .and. (psolver.eq.4.or.psolver.eq.5) )THEN
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *,'  psolver      = ',psolver
        print *
        print *,'  for psolver = 4 or 5, terrain_flag must be .false.'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( (psolver.eq.4.or.psolver.eq.5) .and.    &
          (wbc.eq.2.or.ebc.eq.2.or.sbc.eq.2.or.nbc.eq.2) )THEN
        print *
        print *,'  psolver = ',psolver
        print *
        print *,'  cannot use open boundary conditions for psolver = 4 and 5 (at the moment)'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(terrain_flag .and. ibalance.eq.2)THEN
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *,'  ibalance     = ',ibalance
        print *
        print *,'  for ibalance.eq.2, terrain_flag must be .false.'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(terrain_flag .and. psolver.eq.1)THEN
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *,'  psolver      = ',psolver
        print *
        print *,'  for psolver.eq.1, terrain_flag must be .false.'
        print *,'  (dunno why.  ask George.)'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(iinit.eq.6)THEN
        print *
        print *,'  iinit        = ',iinit
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( (output_format.le.0) .or. (output_format.ge.6) )THEN
        print *
        print *,'  output_format = ',output_format
        print *
        print *,'  only output_format = 1,2,3,4,5 are currently supported'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( (axisymm.eq.1) .and. (iorigin.ne.1) )THEN
        print *
        print *,'  iorigin = ',iorigin
        print *
        print *,'  axisymm=1 requires iorigin=1'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( (axisymm.eq.1) .and. (imove.ne.0) )THEN
        print *
        print *,'  imove = ',imove
        print *
        print *,'  axisymm=1 requires imove=0'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( terrain_flag .and. (imove.ne.0) )THEN
        print *
        print *,'  imove = ',imove
        print *
        print *,'  imove must be 0 when using terrain '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( (axisymm.eq.1) .and. terrain_flag )THEN
        print *
        print *,'  terrain_flag = ',terrain_flag
        print *
        print *,'  axisymm=1 cannot be used with terrain '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( icor.eq.0 ) fcor = 0.0
      IF( (axisymm.eq.1) .and. (wbc.ne.3) )THEN
        print *
        print *,'  wbc = ',wbc
        print *
        print *,'  axisymm=1 requires wbc=3 '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( (axisymm.eq.1) .and. ( (sbc.ne.1).or.(nbc.ne.1) ) )THEN
        print *
        print *,'  sbc = ',sbc
        print *,'  nbc = ',nbc
        print *
        print *,'  axisymm=1 requires sbc=nbc=1 '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( (axisymm.eq.1).and.(ny.gt.1) )THEN
        print *
        print *,'  ny = ',ny
        print *
        print *,'  axisymm=1 requires ny=1'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( (axisymm.eq.1.and.iturb.ge.1).and.iturb.ne.3 )THEN
        print *
        print *,'  iturb    = ',iturb
        print *
        print *,'  axisymm=1 is only available with iturb=3'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( axisymm.eq.1.and.(psolver.lt.2.or.psolver.gt.3) )THEN
        print *
        print *,'  psolver    = ',psolver
        print *
        print *,'  axisymm=1 is only available with psolver=2 or 3'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( (iturb.eq.3) .and. (tconfig.ne.2) )THEN
        print *
        print *,'  iturb    = ',iturb
        print *,'  tconfig  = ',tconfig
        print *
        print *,'  iturb=3 requires tconfig=2'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( (idrag.eq.1).and.(bcturbu.ne.3) )THEN
        print *
        print *,'  idrag    = ',idrag
        print *,'  bcturbu  = ',bcturbu
        print *
        print *,'  idrag=1 requires bcturbu=3'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( (idrag.eq.1).or.(isfcflx.eq.1) )THEN
      IF( iturb.eq.0 .and. ipbl.eq.0 )THEN
        print *
        print *,'  idrag    = ',idrag
        print *,'  isfcflx  = ',isfcflx
        print *
        print *,'  these options require the use of a subgrid turbulence scheme'
        print *
        print *,'  iturb    = ',iturb
        print *,'  ipbl     = ',ipbl
        print *
        print *,'  Use iturb = 1,2,3 or ipbl = 1'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      ENDIF
      IF( ipbl.ge.1 .and. (iturb.eq.1.or.iturb.eq.2) )THEN
        print *
        print *,'  ipbl  = ',ipbl
        print *,'  iturb = ',iturb
        print *
        print *,'  cannot use PBL scheme and LES subgrid turbulence scheme at same time '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( ipbl.ge.1 .and. iturb.eq.3 .and. abs(l_inf).gt.1.0e-6 )THEN
        print *
        print *,'  ipbl  = ',ipbl
        print *,'  iturb = ',iturb
        print *,'  l_inf = ',l_inf
        print *
        print *,'  the PBL scheme requires l_inf = 0 '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
    IF( isfcflx.ne.0 )THEN
      IF( sfcmodel.lt.1 .or. sfcmodel.gt.2 )THEN
        print *
        print *,'  sfcmodel   = ',sfcmodel
        print *
        print *,'  sfcmodel must be 1,2 '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
    ELSE
      IF( sfcmodel.ne.0 )THEN
        print *
        print *,'  isfcflx    = ',isfcflx
        print *,'  sfcmodel   = ',sfcmodel
        print *
        print *,'  sfcmodel must be 0 for isfcflx = 0 '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
    ENDIF
      IF( sfcmodel.eq.2.and.imove.ne.0 )THEN
        print *
        print *,'  sfcmodel = ',sfcmodel
        print *,'  imove     = ',imove
        print *
        print *,'  domain translation is now allowed with sfcmodel = 2 '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( sfcmodel.eq.2.and.(season.le.0.or.season.ge.3) )THEN
        print *
        print *,'  sfcmodel = ',sfcmodel
        print *,'  season   = ',season
        print *
        print *,'  season must have a value of 1 or 2 '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( pertflx.eq.1 .and. sfcmodel.ge.2 )THEN
        print *
        print *,'  pertflx  = ',pertflx
        print *,'  sfcmodel = ',sfcmodel
        print *
        print *,'  pertflx can only be used with sfcmodel = 1  '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( sfcmodel.eq.1 .and. oceanmodel.ne.1 )THEN
        print *
        print *,'  sfcmodel   = ',sfcmodel
        print *,'  oceanmodel = ',oceanmodel
        print *
        print *,'  sfcmodel = 1 requires oceanmodel = 1 '
        print *,'  (oceanmodel = 2 requires sfcmodel = 2 ) '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( radopt.lt.0 .or. radopt.gt.1 )THEN
        print *
        print *,'  radopt   = ',radopt
        print *
        print *,'  radopt must be 0 or 1 '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( radopt.eq.1 .and. imoist.eq.0 )THEN
        print *
        print *,'  radopt   = ',radopt
        print *,'  imoist   = ',imoist
        print *
        print *,'  radopt=1 requires imoist=1 (for now) '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( ipbl.eq.1 .and. imoist.eq.0 )THEN
        print *
        print *,'  ipbl     = ',ipbl
        print *,'  imoist   = ',imoist
        print *
        print *,'  ipbl=1 requires imoist=1 (for now) '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( radopt.eq.1 .and. rterm.eq.1 )THEN
        print *
        print *,'  radopt   = ',radopt
        print *,'  rterm    = ',rterm
        print *
        print *,'  cannot use radopt and rterm at the same time '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( radopt.eq.1 .and. (ptype.eq.1.or.ptype.eq.6)  )THEN
        print *
        print *,'  radopt   = ',radopt
        print *,'  ptype    = ',ptype
        print *
        print *,'  radopt=1 requires an ice microphysics scheme (for now) '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( radopt.eq.1 .and. sfcmodel.eq.0 )THEN
        print *
        print *,'  radopt   = ',radopt
        print *,'  sfcmodel = ',sfcmodel
        print *
        print *,'  radopt=1 requires a surface model '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF( sfcmodel.eq.2 .and. imoist.eq.0 )THEN
        print *
        print *,'  sfcmodel = ',sfcmodel
        print *,'  imoist   = ',imoist
        print *
        print *,'  sfcmodel=2 requires imoist=1 (for now) '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
#ifndef NETCDF
      IF(output_format.eq.2)THEN
        print *
        print *,'  output_format = ',output_format
        print *
        print *,'  You have requested netcdf output, but you have not'
        print *,'  compiled the code with netcdf capability.  Modify the'
        print *,'  Makefile, clean, and recompile'
        print *
        call stopcm1
      ENDIF
#endif
#ifdef NETCDF
      IF( terrain_flag .and. output_interp.ne.0 .and. output_format.eq.2 )THEN
        print *
        print *,'  output_interp = ',output_interp
        print *
        print *,'  output_interp=1 is not currently available for netcdf output'
        print *
        call stopcm1
      ENDIF
#endif
#ifndef HDFOUT
      IF(output_format.eq.3.or.output_format.eq.4.or.output_format.eq.5)THEN
        print *
        print *,'  output_format = ',output_format
        print *
        print *,'  You have requested hdf output, but you have not'
        print *,'  compiled the code with hdf capability.  Modify the'
        print *,'  Makefile, clean, and recompile'
        print *
        call stopcm1
      ENDIF
#endif
#ifdef MPI
      IF(psolver.eq.4.or.psolver.eq.5)THEN
        print *
        print *,'  psolver = ',psolver
        print *
        print *,'  psolver = 4 and 5 are not supported in MPI mode'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
      IF(axisymm.eq.1)THEN
        print *
        print *,'  axisymm = ',axisymm
        print *
        print *,'  axisymm is not supported in MPI mode'
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF
#endif

!-----------------------------------------------------------------------

#ifdef MPI
      endif

      call MPI_BCAST(dx    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dy    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dz    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dtl   ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(timax ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tapfrq,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rstfrq,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(statfrq,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prclfrq,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(adapt_dt ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(irst     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rstnum   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iconly   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(hadvorder,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vadvorder,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pdscheme ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(advweno  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(idiff    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vdiff    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mdiff    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(difforder,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(imoist   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iturb    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tconfig  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(bcturbu  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(bcturbs  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dns      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(irdamp   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(hrdamp   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(psolver  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nsound   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ptype    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ihail    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iautoc   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(icor     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pertcor  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(neweqts  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(idiss    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(efall    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rterm    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(wbc      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ebc      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sbc      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nbc      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(irbc     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(roflux   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(isnd     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iwnd     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(itern    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iinit    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(irandp   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ibalance ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iorigin  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(axisymm  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(imove    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iptra    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(npt      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iprcl    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nparcels ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(kdiff2 ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kdiff6 ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fcor   ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kdiv   ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(alph   ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rdalpha,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(zd     ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xhd    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(umove  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(vmove  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(v_t    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(l_h    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(l_inf  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ndcnst ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(stretch_x,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dx_inner ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dx_outer ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nos_x_len,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tot_x_len,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(stretch_y,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dy_inner ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dy_outer ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nos_y_len,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tot_y_len,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(stretch_z,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ztop     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(str_bot  ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(str_top  ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dz_bot   ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dz_top   ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(bc_wind  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(bc_temp  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ptc_top  ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ptc_bot  ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(viscosity,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pr_num   ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(var1     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var2     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var3     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var4     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var5     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var6     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var7     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var8     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var9     ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(var10    ,1,MPI_REAL   ,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(output_path    ,70,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_basename,70,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_format  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_filetype,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_interp ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_rain   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_sws    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_coldpool,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_sfcflx ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_sfcparams,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_sfcdiags,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_zs     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_zh     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_basestate,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_th     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_thpert ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_prs    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_prspert,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_pi     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_pipert ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_rho    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_rhopert,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_tke    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_km     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_kh     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_qv     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_qvpert ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_q      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_dbz    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_u      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_upert  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_uinterp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_v      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_vpert  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_vinterp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_w      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_winterp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_vort   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_uh     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_pblten ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_dissten,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(output_radten ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(stat_w      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_u      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_v      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_rmw    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_pipert ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_prspert,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_thpert ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_q      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_tke    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_km     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_kh     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_div    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_rh     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_rhi    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_the    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_cloud  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_sfcprs ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_wsp    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_cfl    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_vort   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_tmass  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_tmois  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_qmass  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_tenerg ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_mo     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_tmf    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_pcn    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stat_qsrc   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(radopt   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dtrad  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ctrlat ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ctrlon ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(year   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(month  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(day    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(hour   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(minute ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(second ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(idrag     ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(isfcflx   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sfcmodel  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(oceanmodel,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ipbl      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(initsfc   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tsk0      ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tmn0      ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xland0    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lu0       ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(season    ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cecd      ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pertflx   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cnstce    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cnstcd    ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(isftcflx  ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iz0tlnd   ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(oml_hml0  ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(oml_gamma ,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)

!--------------------------------------------------------------
!  Check domain size

      if(myid.eq.0)then

        if((nodex*nodey).ne.numprocs)then
          print *,'  WARNING!!! '
          print *,'  nodes does not equal numprocs!'
          print *,'  nodex,nodey,nodes=',nodex,nodey,nodex*nodey
          print *,'  numprocs=',numprocs
          call stopcm1
        endif

        if(mod(nx,nodex).ne.0)then
          print *,'  nx does not divide exactly by nodex! '
          print *,'  nx,nodex,mod(nx,nodex)=',nx,nodex,mod(nx,nodex)
          call stopcm1
        endif

        if(mod(ny,nodey).ne.0)then
          print *,'  ny does not divide exactly by nodey! '
          print *,'  ny,nodey,mod(ny,nodey)=',ny,nodey,mod(ny,nodey)
          call stopcm1
        endif

        print *
        print *,'  Everything is cool!'
        print *,'  ni,nj=',ni,nj
        print *

      endif

#endif
!--------------------------------------------------------------

      if(ebc.eq.1 .and. wbc.ne.1)then
        print *,"Can not have periodic b.c.'s on one side only!"
        call stopcm1
      endif
      if(wbc.eq.1 .and. ebc.ne.1)then
        print *,"Can not have periodic b.c.'s on one side only!"
        call stopcm1
      endif
      if(nbc.eq.1 .and. sbc.ne.1)then
        print *,"Can not have periodic b.c.'s on one side only!"
        call stopcm1
      endif
      if(sbc.eq.1 .and. nbc.ne.1)then
        print *,"Can not have periodic b.c.'s on one side only!"
        call stopcm1
      endif

!--------------------------------------------------------------
!  Some basic checks:

      iptra    = max(0,min(1,iptra))
      if(iptra.eq.1)then
        npt      = max(1,npt)
      else
        npt      = 1
      endif
      nparcels = max(1,nparcels)

      if(iprcl.eq.1)then

        npvals = 14

      else

        npvals = 1
        nparcels = 1

      endif

      if(stretch_z.ne.1) ztop = dz*float(nk)

!!!      IF ( ptype .ge. 26 ) THEN
!!!      open(unit=20,file='namelist.input',form='formatted',status='old',    &
!!!           access='sequential')
!!!         read(20,nml=micro_params)
!!!      close(unit=20)
!!!      ENDIF

      IF( advweno.eq.2 )THEN
        hadvorder = 5
        vadvorder = 5
      ENDIF
      IF( advweno.lt.0 .or. advweno.gt.2 )THEN
        print *
        print *,'  advweno = ',advweno
        print *
        print *,'  unrecognized value for advweno '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF

!--------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'dx        =',dx
      if(dowr) write(outfile,*) 'dy        =',dy
      if(dowr) write(outfile,*) 'dz        =',dz
      if(dowr) write(outfile,*) 'dtl       =',dtl
      if(dowr) write(outfile,*) 'timax     =',timax
      if(dowr) write(outfile,*) 'tapfrq    =',tapfrq
      if(dowr) write(outfile,*) 'rstfrq    =',rstfrq
      if(dowr) write(outfile,*) 'statfrq   =',statfrq
      if(dowr) write(outfile,*) 'prclfrq   =',prclfrq
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'adapt_dt  =',adapt_dt
      if(dowr) write(outfile,*) 'irst      =',irst
      if(dowr) write(outfile,*) 'rstnum    =',rstnum
      if(dowr) write(outfile,*) 'iconly    =',iconly
      if(dowr) write(outfile,*) 'hadvorder =',hadvorder
      if(dowr) write(outfile,*) 'vadvorder =',vadvorder
      if(dowr) write(outfile,*) 'pdscheme  =',pdscheme
      if(dowr) write(outfile,*) 'advweno   =',advweno
      if(dowr) write(outfile,*) 'idiff     =',idiff
      if(dowr) write(outfile,*) 'vdiff     =',vdiff
      if(dowr) write(outfile,*) 'mdiff     =',mdiff
      if(dowr) write(outfile,*) 'difforder =',difforder
      if(dowr) write(outfile,*) 'imoist    =',imoist
      if(dowr) write(outfile,*) 'iturb     =',iturb
      if(dowr) write(outfile,*) 'tconfig   =',tconfig
      if(dowr) write(outfile,*) 'bcturbu   =',bcturbu
      if(dowr) write(outfile,*) 'bcturbs   =',bcturbs
      if(dowr) write(outfile,*) 'dns       =',dns
      if(dowr) write(outfile,*) 'irdamp    =',irdamp
      if(dowr) write(outfile,*) 'hrdamp    =',hrdamp
      if(dowr) write(outfile,*) 'psolver   =',psolver
      if(dowr) write(outfile,*) 'nsound    =',nsound
      if(dowr) write(outfile,*) 'ptype     =',ptype
      if(dowr) write(outfile,*) 'ihail     =',ihail
      if(dowr) write(outfile,*) 'iautoc    =',iautoc
      if(dowr) write(outfile,*) 'icor      =',icor
      if(dowr) write(outfile,*) 'pertcor   =',pertcor
      if(dowr) write(outfile,*) 'neweqts   =',neweqts
      if(dowr) write(outfile,*) 'idiss     =',idiss
      if(dowr) write(outfile,*) 'efall     =',efall
      if(dowr) write(outfile,*) 'rterm     =',rterm
      if(dowr) write(outfile,*) 'wbc       =',wbc
      if(dowr) write(outfile,*) 'ebc       =',ebc
      if(dowr) write(outfile,*) 'sbc       =',sbc
      if(dowr) write(outfile,*) 'nbc       =',nbc
      if(dowr) write(outfile,*) 'irbc      =',irbc
      if(dowr) write(outfile,*) 'roflux    =',roflux
      if(dowr) write(outfile,*) 'isnd      =',isnd
      if(dowr) write(outfile,*) 'iwnd      =',iwnd
      if(dowr) write(outfile,*) 'itern     =',itern
      if(dowr) write(outfile,*) 'iinit     =',iinit
      if(dowr) write(outfile,*) 'irandp    =',irandp
      if(dowr) write(outfile,*) 'ibalance  =',ibalance
      if(dowr) write(outfile,*) 'iorigin   =',iorigin
      if(dowr) write(outfile,*) 'axisymm   =',axisymm
      if(dowr) write(outfile,*) 'imove     =',imove
      if(dowr) write(outfile,*) 'iptra     =',iptra
      if(dowr) write(outfile,*) 'npt       =',npt
      if(dowr) write(outfile,*) 'iprcl     =',iprcl
      if(dowr) write(outfile,*) 'nparcels  =',nparcels
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'kdiff2    =',kdiff2
      if(dowr) write(outfile,*) 'kdiff6    =',kdiff6
      if(dowr) write(outfile,*) 'fcor      =',fcor
      if(dowr) write(outfile,*) 'kdiv      =',kdiv
      if(dowr) write(outfile,*) 'alph      =',alph
      if(dowr) write(outfile,*) 'rdalpha   =',rdalpha
      if(dowr) write(outfile,*) 'zd        =',zd
      if(dowr) write(outfile,*) 'xhd       =',xhd
      if(dowr) write(outfile,*) 'umove     =',umove
      if(dowr) write(outfile,*) 'vmove     =',vmove
      if(dowr) write(outfile,*) 'v_t       =',v_t
      if(dowr) write(outfile,*) 'l_h       =',l_h
      if(dowr) write(outfile,*) 'l_inf     =',l_inf
      if(dowr) write(outfile,*) 'ndcnst    =',ndcnst
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'radopt    =',radopt
      if(dowr) write(outfile,*) 'dtrad     =',dtrad
      if(dowr) write(outfile,*) 'ctrlat    =',ctrlat
      if(dowr) write(outfile,*) 'ctrlon    =',ctrlon
      if(dowr) write(outfile,*) 'year      =',year
      if(dowr) write(outfile,*) 'month     =',month
      if(dowr) write(outfile,*) 'day       =',day
      if(dowr) write(outfile,*) 'hour      =',hour
      if(dowr) write(outfile,*) 'minute    =',minute
      if(dowr) write(outfile,*) 'second    =',second
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'idrag     =',idrag
      if(dowr) write(outfile,*) 'isfcflx   =',isfcflx
      if(dowr) write(outfile,*) 'sfcmodel  =',sfcmodel
      if(dowr) write(outfile,*) 'oceanmodel=',oceanmodel
      if(dowr) write(outfile,*) 'ipbl      =',ipbl
      if(dowr) write(outfile,*) 'initsfc   =',initsfc
      if(dowr) write(outfile,*) 'tsk0      =',tsk0
      if(dowr) write(outfile,*) 'tmn0      =',tmn0
      if(dowr) write(outfile,*) 'xland0    =',xland0
      if(dowr) write(outfile,*) 'lu0       =',lu0
      if(dowr) write(outfile,*) 'season    =',season
      if(dowr) write(outfile,*) 'cecd      =',cecd
      if(dowr) write(outfile,*) 'pertflx   =',pertflx
      if(dowr) write(outfile,*) 'cnstce    =',cnstce
      if(dowr) write(outfile,*) 'cnstcd    =',cnstcd
      if(dowr) write(outfile,*) 'isftcflx  =',isftcflx
      if(dowr) write(outfile,*) 'iz0tlnd   =',iz0tlnd
      if(dowr) write(outfile,*) 'oml_hml0  =',oml_hml0
      if(dowr) write(outfile,*) 'oml_gamma =',oml_gamma
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'stretch_x =',stretch_x
      if(dowr) write(outfile,*) 'dx_inner  =',dx_inner
      if(dowr) write(outfile,*) 'dx_outer  =',dx_outer
      if(dowr) write(outfile,*) 'nos_x_len =',nos_x_len
      if(dowr) write(outfile,*) 'tot_x_len =',tot_x_len
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'stretch_y =',stretch_y
      if(dowr) write(outfile,*) 'dy_inner  =',dy_inner
      if(dowr) write(outfile,*) 'dy_outer  =',dy_outer
      if(dowr) write(outfile,*) 'nos_y_len =',nos_y_len
      if(dowr) write(outfile,*) 'tot_y_len =',tot_y_len
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'stretch_z =',stretch_z
      if(dowr) write(outfile,*) 'ztop      =',ztop
      if(dowr) write(outfile,*) 'str_bot   =',str_bot
      if(dowr) write(outfile,*) 'str_top   =',str_top
      if(dowr) write(outfile,*) 'dz_bot    =',dz_bot
      if(dowr) write(outfile,*) 'dz_top    =',dz_top
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'bc_wind   =',bc_wind
      if(dowr) write(outfile,*) 'bc_temp   =',bc_temp
      if(dowr) write(outfile,*) 'ptc_top   =',ptc_top
      if(dowr) write(outfile,*) 'ptc_bot   =',ptc_bot
      if(dowr) write(outfile,*) 'viscosity =',viscosity
      if(dowr) write(outfile,*) 'pr_num    =',pr_num
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'var1      =',var1
      if(dowr) write(outfile,*) 'var2      =',var2
      if(dowr) write(outfile,*) 'var3      =',var3
      if(dowr) write(outfile,*) 'var4      =',var4
      if(dowr) write(outfile,*) 'var5      =',var5
      if(dowr) write(outfile,*) 'var6      =',var6
      if(dowr) write(outfile,*) 'var7      =',var7
      if(dowr) write(outfile,*) 'var8      =',var8
      if(dowr) write(outfile,*) 'var9      =',var9
      if(dowr) write(outfile,*) 'var10     =',var10
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  Configuration for simulations with moisture
!

      !--- begin: define defaults (please do not change) ---------
      iice     = 0
      idm      = 0
      numq     = 1
      nqv      = 1
      nql1     = 1
      nql2     = 1
      nqs1     = 1
      nqs2     = 1
      nnc1     = 1
      nnc2     = 1
      nbudget  = 10
      budrain  = 1
      cloudvar = .false.
      rhovar   = .false.
      !--- end: define defaults ----------------------------------

      IF(imoist.eq.1)THEN

!-----------------------------------------------------------------------
!-------   BEGIN:  modify stuff below here -----------------------------
!-----------------------------------------------------------------------

!-----------------------------------------------------------------------
        IF(ptype.eq.1)THEN        ! Kessler scheme

          numq = 3    ! there are 3 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... using Kessler microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
        ELSEIF((ptype.eq.2).or.(ptype.eq.4))THEN    ! Goddard-LFO or 
                                                    ! GSR-LFO scheme

          iice = 1    ! this means that ptype=2,4 are ice schemes

          numq = 6    ! there are 6 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.
          cloudvar(4) = .true.
          cloudvar(5) = .false.
          cloudvar(6) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '
          qname(4) = 'qi '
          qname(5) = 'qs '
          qname(6) = 'qg '

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the Goddard or GSR LFO scheme -----

          if(ptype.eq.2)THEN

            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) 'Calling CONSAT'
            if(dowr) write(outfile,*)

            call consat
            call consat2(dtl)

            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) '  ... using Goddard LFO microphysics scheme ... '
            if(dowr) write(outfile,*) '         numq   = ',numq
            if(dowr) write(outfile,*) '         ihail  = ',ihail
            if(dowr) write(outfile,*)

          endif

          if(ptype.eq.4)then

            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) 'Calling lfoice_init'
            if(dowr) write(outfile,*)

            call lfoice_init(dtl)

            if(dowr) write(outfile,*)
            if(dowr) write(outfile,*) '  ... using GSR LFO microphysics scheme ... '
            if(dowr) write(outfile,*) '         numq   = ',numq
            if(dowr) write(outfile,*)

          endif

!-----------------------------------------------------------------------
        ELSEIF(ptype.eq.3)THEN    ! Thompson scheme

          iice = 1    ! this means that ptype=3 is an ice scheme
          idm  = 1    ! this means that ptype=3 has at least one double moment

          numq = 8    ! there are 8 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array
          nnc1 = 7    ! the first number concentration var is the seventh array
          nnc2 = 8    ! the last number concentration var is the eighth array

          cloudvar(1) = .false.
          cloudvar(2) = .true.
          cloudvar(3) = .false.
          cloudvar(4) = .true.
          cloudvar(5) = .false.
          cloudvar(6) = .false.
          cloudvar(7) = .false.
          cloudvar(8) = .false.

          qname(1) = 'qv '
          qname(2) = 'qc '
          qname(3) = 'qr '
          qname(4) = 'qi '
          qname(5) = 'qs '
          qname(6) = 'qg '
          qname(7) = 'nci'
          qname(8) = 'ncr'

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the Thompson scheme -----

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Calling thompson_init'
          if(dowr) write(outfile,*) '(this can take several minutes ... please be patient)'

          call thompson_init

          if(dowr) write(outfile,*) 'Done with thompson_init'
          if(dowr) write(outfile,*)

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... using Thompson microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------

        ELSEIF(ptype.eq.5)THEN    ! Morrison scheme

          !----- initialize the Morrison scheme -----

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Calling GRAUPEL_INIT'

          call graupel_init(ihail,inum,ndcnst)

          if(dowr) write(outfile,*) 'Returned from GRAUPEL_INIT'
          if(dowr) write(outfile,*)

          !------------------------------------------

          iice = 1    ! this means that ptype=5 is an ice scheme
          idm  = 1    ! this means that ptype=5 has at least one double moment

        if(inum.eq.1)then
          ! constant cloud-drop concentration

          numq = 10   ! there are 10 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array
          nnc1 = 7    ! the first number concentration var is the seventh array
          nnc2 = 10   ! the last number concentration var is the tenth array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '
          qname( 7) = 'nci'
          qname( 8) = 'ncs'
          qname( 9) = 'ncr'
          qname(10) = 'ncg'

        elseif(inum.eq.0)then
          ! cloud-droplet concentration is a predicted variable

          numq = 11   ! there are 11 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array
          nnc1 = 7    ! the first number concentration var is the seventh array
          nnc2 = 11   ! the last number concentration var is the eleventh array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.
          cloudvar(11) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '
          qname( 7) = 'nci'
          qname( 8) = 'ncs'
          qname( 9) = 'ncr'
          qname(10) = 'ncg'
          qname(11) = 'ncc'

        else

          print *,'  unrecognized value for inum '
          call stopcm1

        endif

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... using Morrison microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*) '         ihail  = ',ihail
        if(inum.eq.1)then
          if(dowr) write(outfile,*) '         assuming constant cloud droplet concentration' 
          if(dowr) write(outfile,*) '         ndcnst = ',ndcnst
        else
          if(dowr) write(outfile,*) '         predicting cloud droplet concentration'
        endif
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------

        ELSEIF(ptype.eq.6)THEN        ! Rotunno-Emanuel scheme

          numq = 2    ! there are 2 q variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 2    ! the last liquid variable is the second array

          cloudvar(1) = .false.
          cloudvar(2) = .true.

          qname(1) = 'qv '
          qname(2) = 'ql '

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... using Rotunno-Emanuel microphysics scheme ... '
          if(dowr) write(outfile,*) '         numq   = ',numq
          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------

!!!        ELSEIF(ptype.eq.7)THEN    ! Milbrandt & Yao double-moment scheme

!-----------------------------------------------------------------------

        ELSEIF ( ptype .eq. 26 ) THEN    ! ZVD scheme (without hail)

          iice = 1    ! this means that ptype=26 is an ice scheme
          idm  = 1    ! this means that ptype=26 has at least one double moment

          numq = 14   ! number of variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array
          nnc1 = 7    ! the first number concentration var is the seventh array
          nnc2 = 12   ! the last number concentration var is the eleventh array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.
          cloudvar(11) = .false.
          cloudvar(12) = .false.
          cloudvar(13) = .false.
          cloudvar(14) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '
          qname( 7) = 'ccn' ! CCN concentration
          qname( 8) = 'ccw' ! droplet conc
          qname( 9) = 'crw' ! rain conc
          qname(10) = 'cci' ! ice crystal conc
          qname(11) = 'csw' ! snow conc
          qname(12) = 'chw' ! graupel conc
          qname(13) = 'ss ' ! max supersaturation
          qname(14) = 'vhw' ! graupel volume

          rhovar( 1) = .false.
          rhovar( 2) = .false.
          rhovar( 3) = .false.
          rhovar( 4) = .false.
          rhovar( 5) = .false.
          rhovar( 6) = .false.
          rhovar( 7) = .true.
          rhovar( 8) = .true.
          rhovar( 9) = .true.
          rhovar(10) = .true.
          rhovar(11) = .true.
          rhovar(12) = .true.
          rhovar(13) = .false.
          rhovar(14) = .true.

!          ipconc = 5
!          lr = 4
!          li = 5
!          ls = 6
!          lh = 7
!          lg = lh
!          lhab = lh
!          lhl = 0
!          lqe  = lhab
!
!          lccn = 8
!          lnc  = 9
!          lnr  = 10
!          lni  = 11
!          lns  = 12
!          lnh  = 13
!          lnhl = 0
!          lss  = 14
!          lvh  = 15
!
!          lsch = 0
!          lschab = 0
!          lscw = 0
!          lscb = lscw
!          lscni = 0
!          lscpi = 0
!          lsce = lscni
!          lsceq= lschab
!
!          lsw  = 0
!          lhw  = 0
!          lhlw = 0

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the ZVD scheme -----

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Calling index_module_init'
          if(dowr) write(outfile,*)

!          call INDEX_MODULE_INIT(ptype)

          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) 'Returned from index_module_init'
          if(dowr) write(outfile,*)

        ELSEIF ( ptype .eq. 27 ) THEN    ! ZVDHV scheme (with hail)

          iice = 1    ! this means that ptype=27 is an ice scheme
          idm  = 1    ! this means that ptype=27 has at least one double moment

          numq = 17   ! number of variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 7    ! the last solid variable is the sixth array
          nnc1 = 8    ! the first number concentration var is the seventh array
          nnc2 = 14   ! the last number concentration var is the eleventh array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.
          cloudvar( 7) = .false.
          cloudvar( 8) = .false.
          cloudvar( 9) = .false.
          cloudvar(10) = .false.
          cloudvar(11) = .false.
          cloudvar(12) = .false.
          cloudvar(13) = .false.
          cloudvar(14) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '
          qname( 7) = 'qhl'
          qname( 8) = 'ccn' ! CCN concentration
          qname( 9) = 'ccw' ! droplet conc
          qname(10) = 'crw' ! rain conc
          qname(11) = 'cci' ! ice crystal conc
          qname(12) = 'csw' ! snow conc
          qname(13) = 'chw' ! graupel conc
          qname(14) = 'chl' ! hail conc
          qname(15) = 'ss ' ! max supersaturation
          qname(16) = 'vhw' ! graupel volume
          qname(17) = 'vhl' ! hail volume

          rhovar( 1) = .false.
          rhovar( 2) = .false.
          rhovar( 3) = .false.
          rhovar( 4) = .false.
          rhovar( 5) = .false.
          rhovar( 6) = .false.
          rhovar( 7) = .false.
          rhovar( 8) = .true.
          rhovar( 9) = .true.
          rhovar(10) = .true.
          rhovar(11) = .true.
          rhovar(12) = .true.
          rhovar(13) = .true.
          rhovar(14) = .true.
          rhovar(15) = .false.
          rhovar(16) = .true.
          rhovar(17) = .true.

!          ipconc = 5
!          lr = 4
!          li = 5
!          ls = 6
!          lh = 7
!          lhl = 8
!          lg = lh
!          lhab = lhl
!          lqe  = lhab
!
!          lccn = 9
!          lnc  = 10
!          lnr  = 11
!          lni  = 12
!          lns  = 13
!          lnh  = 14
!          lnhl = 15
!          lss  = 16
!          lvh  = 17
!          lvhl = 18
!
!          lsch = 0
!          lschab = 0
!          lscw = 0
!          lscb = lscw
!          lscni = 0
!          lscpi = 0
!          lsce = lscni
!          lsceq= lschab
!
!          lsw  = 0
!          lhw  = 0
!          lhlw = 0

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the ZVD scheme -----

!          if(dowr) write(outfile,*)
!          if(dowr) write(outfile,*) 'Calling graupel_init'
!          if(dowr) write(outfile,*)

!          call INDEX_MODULE_INIT(ptype)

!          if(dowr) write(outfile,*)
!          if(dowr) write(outfile,*) 'Returned from graupel_init'
!          if(dowr) write(outfile,*)

        ELSEIF ( ptype .eq. 28 ) THEN    ! single moment ZIEG scheme (without hail)

          iice = 1    ! this means that ptype=28 is an ice scheme
          idm  = 0    ! this means that ptype=28 has at least one double moment

          numq = 6   ! number of variables

          nqv  = 1    ! qv is the first array
          nql1 = 2    ! the first liquid variable is the second array
          nql2 = 3    ! the last liquid variable is the third array
          nqs1 = 4    ! the first solid variable is the fourth array
          nqs2 = 6    ! the last solid variable is the sixth array

          cloudvar( 1) = .false.
          cloudvar( 2) = .true.
          cloudvar( 3) = .false.
          cloudvar( 4) = .true.
          cloudvar( 5) = .false.
          cloudvar( 6) = .false.

          qname( 1) = 'qv '
          qname( 2) = 'qc '
          qname( 3) = 'qr '
          qname( 4) = 'qi '
          qname( 5) = 'qs '
          qname( 6) = 'qg '

          rhovar( 1) = .false.
          rhovar( 2) = .false.
          rhovar( 3) = .false.
          rhovar( 4) = .false.
          rhovar( 5) = .false.
          rhovar( 6) = .false.

!          ipconc = 0
!          lr = 4
!          li = 5
!          ls = 6
!          lh = 7
!          lg = lh
!          lhab = lh
!          lhl = 0
!          lqe  = lhab
!
!          lccn = 0
!          lnc  = 0
!          lnr  = 0
!          lni  = 0
!          lns  = 0
!          lnh  = 0
!          lnhl = 0
!          lss  = 0
!          lvh  = 0
!
!          lsch = 0
!          lschab = 0
!          lscw = 0
!          lscb = lscw
!          lscni = 0
!          lscpi = 0
!          lsce = lscni
!          lsceq= lschab
!
!          lsw  = 0
!          lhw  = 0
!          lhlw = 0

          !----- budget stuff below here -----

          nbudget = 10

          budname(1) = 'tcond '
          budname(2) = 'tevac '
          budname(3) = 'tauto '
          budname(4) = 'taccr '
          budname(5) = 'tevar '
          budname(6) = 'train '
          budname(7) = 'erain '
          budname(8) = 'qsfc  '
          budname(9) = 'esfc  '
          budname(10) = 'erad  '

          budrain = 6

          !----- initialize the ZVD scheme -----

!          if(dowr) write(outfile,*)
!          if(dowr) write(outfile,*) 'Calling graupel_init'
!          if(dowr) write(outfile,*)

!          call INDEX_MODULE_INIT(ptype)

!          if(dowr) write(outfile,*)
!          if(dowr) write(outfile,*) 'Returned from graupel_init'
!          if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!  insert new ptype here

!!!        ELSEIF(ptype.eq.8)THEN    ! new microphysics scheme

!-----------------------------------------------------------------------

        ELSE

          IF(myid.eq.0)THEN
            print *
            print *,'  ptype = ',ptype
            print *
            print *,'  Unrecognized value for ptype '
            print *
            print *,'  ... stopping cm1 ... '
            print *
          ENDIF

          call stopcm1

        ENDIF    ! endif for ptype

      ENDIF    ! endif for imoist=1

!-----------------------------------------------------------------------
!-------   END:  modify stuff above here -------------------------------
!-----------------------------------------------------------------------

      IF( radopt.eq.1 .and. iice.ne.1 )THEN
        print *
        print *,'  radopt   = ',radopt
        print *,'  iice     = ',iice
        print *
        print *,'  radopt=1 requires an ice microphysics scheme '
        print *
        print *,'   stopping model .... '
        print *
        call stopcm1
      ENDIF

!-----------------------------------------------------------------------

      nqc = 0
      nqr = 0
      nqi = 0
      nqs = 0
      nqg = 0

      do n=1,numq
        if( qname(n).eq.'qc ' .or. qname(n).eq.'ql ' ) nqc = n
        if( qname(n).eq.'qr ' ) nqr = n
        if( qname(n).eq.'qi ' ) nqi = n
        if( qname(n).eq.'qs ' ) nqs = n
        if( qname(n).eq.'qg ' ) nqg = n
      enddo

      if(numq .gt. maxq)then
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  WARNING!   numq > maxq'
        if(dowr) write(outfile,*) '  You need to increase maxq in input.incl and recompile'
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Stopping model ....'
        if(dowr) write(outfile,*)
        call stopcm1
      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'iice      =',iice
      if(dowr) write(outfile,*) 'idm       =',idm
      if(dowr) write(outfile,*) 'numq      =',numq
      if(dowr) write(outfile,*) 'nqv       =',nqv
      if(dowr) write(outfile,*) 'nqc       =',nqc
      if(dowr) write(outfile,*) 'nqr       =',nqr
      if(dowr) write(outfile,*) 'nqi       =',nqi
      if(dowr) write(outfile,*) 'nqs       =',nqs
      if(dowr) write(outfile,*) 'nqg       =',nqg
      if(dowr) write(outfile,*) 'nql1      =',nql1
      if(dowr) write(outfile,*) 'nql2      =',nql2
      if(dowr) write(outfile,*) 'nqs1      =',nqs1
      if(dowr) write(outfile,*) 'nqs2      =',nqs2
      if(dowr) write(outfile,*) 'nnc1      =',nnc1
      if(dowr) write(outfile,*) 'nnc2      =',nnc2
      if(dowr) write(outfile,*)

!--------------------------------------------------------------

      iterrain = 0
      if(terrain_flag) iterrain = 1

      output_interp  = max(0,min(1,output_interp))*iterrain
      output_rain    = max(0,min(1,output_rain))
      output_sws     = max(0,min(1,output_sws))
      output_coldpool= max(0,min(1,output_coldpool))
      output_sfcflx  = max(0,min(1,output_sfcflx))
      output_sfcparams = max(0,min(1,output_sfcparams))
      output_sfcdiags = max(0,min(1,output_sfcdiags))
      output_zs      = max(0,min(1,output_zs))*iterrain
      output_zh      = max(0,min(1,output_zh))
      output_basestate = max(0,min(1,output_basestate))
      output_th      = max(0,min(1,output_th))
      output_thpert  = max(0,min(1,output_thpert))
      output_prs     = max(0,min(1,output_prs))
      output_prspert = max(0,min(1,output_prspert))
      output_pi      = max(0,min(1,output_pi))
      output_pipert  = max(0,min(1,output_pipert))
      output_rho     = max(0,min(1,output_rho))
      output_rhopert = max(0,min(1,output_rhopert))
      output_tke     = max(0,min(1,output_tke))
      output_km      = max(0,min(1,output_km))
      output_kh      = max(0,min(1,output_kh))
      output_qv      = max(0,min(1,output_qv))
      output_qvpert  = max(0,min(1,output_qvpert))
      output_q       = max(0,min(1,output_q))
      output_dbz     = max(0,min(1,output_dbz))
      output_u       = max(0,min(1,output_u))
      output_upert   = max(0,min(1,output_upert))
      output_uinterp = max(0,min(1,output_uinterp))
      output_v       = max(0,min(1,output_v))
      output_vpert   = max(0,min(1,output_vpert))
      output_vinterp = max(0,min(1,output_vinterp))
      output_w       = max(0,min(1,output_w))
      output_winterp = max(0,min(1,output_winterp))
      output_vort    = max(0,min(1,output_vort))
      output_uh      = max(0,min(1,output_uh))
      output_pblten  = max(0,min(1,output_pblten))
      output_dissten = max(0,min(1,output_dissten))
      output_radten  = max(0,min(1,output_radten))


      nrain = 1
!!!      if(imove.eq.1.and.imoist.eq.1) nrain = 2
      if(imove.eq.1)                 nrain = 2

      if(dowr) write(outfile,*) 'nrain     =',nrain
      if(dowr) write(outfile,*)

      if(imoist.eq.0)then
        output_rain=0
        output_qv=0
        output_qvpert=0
        output_q=0
        output_dbz=0
      endif
      if( (iturb.eq.0.or.dns.eq.1) )then
        output_tke=0
      endif
      if( (iturb.eq.0.or.dns.eq.1).and.(ipbl.eq.0) )then
        output_km=0
        output_kh=0
      endif
      if(iturb.eq.2.or.iturb.eq.3)then
        output_tke=0
      endif
      if(idiss.ne.1)then
        output_dissten=0
      endif
      if(ipbl.ne.1)then
        output_pblten=0
      endif
      if(radopt.eq.0)then
        output_radten=0
      endif
      if( idrag.eq.0 .and. isfcflx.eq.0 )then
        output_sfcflx = 0
        output_sfcparams = 0
        output_sfcdiags = 0
      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'output_path      = ',output_path
      if(dowr) write(outfile,*) 'output_basename  = ',output_basename
      if(dowr) write(outfile,*) 'output_format    =',output_format
      if(dowr) write(outfile,*) 'output_filetype  =',output_filetype
      if(dowr) write(outfile,*) 'output_interp    =',output_interp
      if(dowr) write(outfile,*) 'output_rain      =',output_rain
      if(dowr) write(outfile,*) 'output_sws       =',output_sws
      if(dowr) write(outfile,*) 'output_coldpool  =',output_coldpool
      if(dowr) write(outfile,*) 'output_sfcflx    =',output_sfcflx
      if(dowr) write(outfile,*) 'output_sfcparams =',output_sfcparams
      if(dowr) write(outfile,*) 'output_sfcdiags  =',output_sfcdiags
      if(dowr) write(outfile,*) 'output_zs        =',output_zs
      if(dowr) write(outfile,*) 'output_zh        =',output_zh
      if(dowr) write(outfile,*) 'output_basestate =',output_basestate
      if(dowr) write(outfile,*) 'output_th        =',output_th
      if(dowr) write(outfile,*) 'output_thpert    =',output_thpert
      if(dowr) write(outfile,*) 'output_prs       =',output_prs
      if(dowr) write(outfile,*) 'output_prspert   =',output_prspert
      if(dowr) write(outfile,*) 'output_pi        =',output_pi
      if(dowr) write(outfile,*) 'output_pipert    =',output_pipert
      if(dowr) write(outfile,*) 'output_rho       =',output_rho
      if(dowr) write(outfile,*) 'output_rhopert   =',output_rhopert
      if(dowr) write(outfile,*) 'output_tke       =',output_tke
      if(dowr) write(outfile,*) 'output_km        =',output_km
      if(dowr) write(outfile,*) 'output_kh        =',output_kh
      if(dowr) write(outfile,*) 'output_qv        =',output_qv
      if(dowr) write(outfile,*) 'output_qvpert    =',output_qvpert
      if(dowr) write(outfile,*) 'output_q         =',output_q
      if(dowr) write(outfile,*) 'output_dbz       =',output_dbz
      if(dowr) write(outfile,*) 'output_u         =',output_u
      if(dowr) write(outfile,*) 'output_upert     =',output_upert
      if(dowr) write(outfile,*) 'output_uinterp   =',output_uinterp
      if(dowr) write(outfile,*) 'output_v         =',output_v
      if(dowr) write(outfile,*) 'output_vpert     =',output_vpert
      if(dowr) write(outfile,*) 'output_vinterp   =',output_vinterp
      if(dowr) write(outfile,*) 'output_w         =',output_w
      if(dowr) write(outfile,*) 'output_winterp   =',output_winterp
      if(dowr) write(outfile,*) 'output_vort      =',output_vort
      if(dowr) write(outfile,*) 'output_uh        =',output_uh
      if(dowr) write(outfile,*) 'output_pblten    =',output_pblten
      if(dowr) write(outfile,*) 'output_dissten   =',output_dissten
      if(dowr) write(outfile,*) 'output_radten    =',output_radten
      if(dowr) write(outfile,*)

!--------------------------------------------------------------


      stat_w       = max(0,min(1,stat_w))
      stat_u       = max(0,min(1,stat_u))
      stat_v       = max(0,min(1,stat_v))
      stat_rmw     = max(0,min(1,stat_rmw))
      IF(axisymm.ne.1) stat_rmw = 0
      stat_pipert  = max(0,min(1,stat_pipert))
      stat_prspert = max(0,min(1,stat_prspert))
      stat_thpert  = max(0,min(1,stat_thpert))
      stat_q       = max(0,min(1,stat_q))
      stat_tke     = max(0,min(1,stat_tke))
      stat_km      = max(0,min(1,stat_km))
      stat_kh      = max(0,min(1,stat_kh))
      stat_div     = max(0,min(1,stat_div))
      stat_rh      = max(0,min(1,stat_rh))
      stat_rhi     = max(0,min(1,stat_rhi))
      stat_the     = max(0,min(1,stat_the))
      stat_cloud   = max(0,min(1,stat_cloud))
      stat_sfcprs  = max(0,min(1,stat_sfcprs))
      stat_wsp     = max(0,min(1,stat_wsp))
      stat_cfl     = max(0,min(1,stat_cfl))
      stat_vort    = max(0,min(1,stat_vort))
      stat_tmass   = max(0,min(1,stat_tmass))
      stat_tmois   = max(0,min(1,stat_tmois))
      stat_qmass   = max(0,min(1,stat_qmass))
      stat_tenerg  = max(0,min(1,stat_tenerg))
      stat_mo      = max(0,min(1,stat_mo))
      stat_tmf     = max(0,min(1,stat_tmf))
      stat_pcn     = max(0,min(1,stat_pcn))
      stat_qsrc    = max(0,min(1,stat_qsrc))


      if(imoist.eq.0)then
        stat_q=0
        stat_rh=0
        stat_rhi=0
        stat_the=0
        stat_cloud=0
        stat_tmois=0
        stat_qmass=0
        stat_pcn=0
        stat_qsrc=0
      endif
      if(iice.eq.0)then
        stat_rhi=0
      endif
      if(iturb.eq.0.or.dns.eq.1)then
        stat_tke=0
        stat_km=0
        stat_kh=0
      endif 
      if(iturb.eq.2.or.iturb.eq.3)then
        stat_tke=0
      endif


      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'stat_w       = ',stat_w
      if(dowr) write(outfile,*) 'stat_u       = ',stat_u
      if(dowr) write(outfile,*) 'stat_v       = ',stat_v
      if(dowr) write(outfile,*) 'stat_rmw     = ',stat_rmw
      if(dowr) write(outfile,*) 'stat_pipert  = ',stat_pipert
      if(dowr) write(outfile,*) 'stat_prspert = ',stat_prspert
      if(dowr) write(outfile,*) 'stat_thpert  = ',stat_thpert
      if(dowr) write(outfile,*) 'stat_q       = ',stat_q
      if(dowr) write(outfile,*) 'stat_tke     = ',stat_tke
      if(dowr) write(outfile,*) 'stat_km      = ',stat_km
      if(dowr) write(outfile,*) 'stat_kh      = ',stat_kh
      if(dowr) write(outfile,*) 'stat_div     = ',stat_div
      if(dowr) write(outfile,*) 'stat_rh      = ',stat_rh
      if(dowr) write(outfile,*) 'stat_rhi     = ',stat_rhi
      if(dowr) write(outfile,*) 'stat_the     = ',stat_the
      if(dowr) write(outfile,*) 'stat_cloud   = ',stat_cloud
      if(dowr) write(outfile,*) 'stat_sfcprs  = ',stat_sfcprs
      if(dowr) write(outfile,*) 'stat_wsp     = ',stat_wsp
      if(dowr) write(outfile,*) 'stat_cfl     = ',stat_cfl
      if(dowr) write(outfile,*) 'stat_vort    = ',stat_vort
      if(dowr) write(outfile,*) 'stat_tmass   = ',stat_tmass
      if(dowr) write(outfile,*) 'stat_tmois   = ',stat_tmois
      if(dowr) write(outfile,*) 'stat_qmass   = ',stat_qmass
      if(dowr) write(outfile,*) 'stat_tenerg  = ',stat_tenerg
      if(dowr) write(outfile,*) 'stat_mo      = ',stat_mo
      if(dowr) write(outfile,*) 'stat_tmf     = ',stat_tmf
      if(dowr) write(outfile,*) 'stat_pcn     = ',stat_pcn
      if(dowr) write(outfile,*) 'stat_qsrc    = ',stat_qsrc
      if(dowr) write(outfile,*)

      stat_out=2*(stat_w+stat_pipert+stat_prspert+numq*stat_q+              &
              stat_tke+2*stat_km+2*stat_kh+stat_div+stat_rh+stat_rhi+       &
              stat_cloud+stat_sfcprs+2*stat_wsp)  +                         &
              4*(stat_thpert+stat_u+stat_v)  + 1*stat_rmw +                 &
              3*stat_cfl  +  6*stat_vort  +  stat_tmass  +  stat_tmois  +   &
              (1+(1+nql2-nql1)+iice*(1+nqs2-nqs1))*stat_qmass +             &
              5*stat_tenerg  +  3*stat_mo  +                                &
              nbudget*stat_pcn  + numq*2*stat_qsrc +                        &
              4*stat_the  +  2*stat_tmf + 2*iptra*npt
      IF( adapt_dt.eq.1 ) stat_out = stat_out + 1
      IF( stat_wsp.eq.1 .and. idrag.eq.1 ) stat_out = stat_out + 2
      IF( stat_cfl.ge.1 .and. iturb.lt.1 ) stat_out = stat_out - 2
      stat_out = max(1,stat_out)

      if(dowr) write(outfile,*) 'stat_out = ',stat_out
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  Define dimensions for allocatable arrays

      if(imoist.eq.1)then
        ibm=ib
        iem=ie
        jbm=jb
        jem=je
        kbm=kb
        kem=ke
        if(ptype.ge.26)then
          ibzvd=ib
          iezvd=ie
          jbzvd=jb
          jezvd=je
          kbzvd=kb
          kezvd=ke
          nqzvd = numq + 1
        else
          ibzvd=1
          iezvd=1
          jbzvd=1
          jezvd=1
          kbzvd=1
          kezvd=1
          nqzvd=1
        endif
      else
        ibm=1
        iem=1
        jbm=1
        jem=1
        kbm=1
        kem=1
        ibzvd=1
        iezvd=1
        jbzvd=1
        jezvd=1
        kbzvd=1
        kezvd=1
        nqzvd=1
      endif

      if(iice.eq.1)then
        ibi=ib
        iei=ie
        jbi=jb
        jei=je
        kbi=kb
        kei=ke
      else
        ibi=1
        iei=1
        jbi=1
        jei=1
        kbi=1
        kei=1
      endif

      if(radopt.eq.1)then
        ibr=ib
        ier=ie
        jbr=jb
        jer=je
        kbr=kb
        ker=ke
      else
        ibr=1
        ier=1
        jbr=1
        jer=1
        kbr=1
        ker=1
      endif

      if(ipbl.eq.1)then
        ibb=ib
        ieb=ie
        jbb=jb
        jeb=je
        kbb=kb
        keb=ke
      else
        ibb=1
        ieb=1
        jbb=1
        jeb=1
        kbb=1
        keb=1
      endif

      if( (sfcmodel.ge.1) .or. (oceanmodel.eq.2) .or. (ipbl.eq.1) .or. (radopt.eq.1) .or. (idrag.eq.1) )then
        ibl=ib
        iel=ie
        jbl=jb
        jel=je
      else
        ibl=1
        iel=1
        jbl=1
        jel=1
      endif

      if((iturb.ge.1).or.(ipbl.eq.1))then
        ibc=ib
        iec=ie
        jbc=jb
        jec=je
        kbc=kb
        kec=ke+1
      else
        ibc=1
        iec=1
        jbc=1
        jec=1
        kbc=1
        kec=1
      endif

      if(iturb.eq.1)then
        ibt=ib
        iet=ie
        jbt=jb
        jet=je
        kbt=kb
        ket=ke+1
      else
        ibt=1
        iet=1
        jbt=1
        jet=1
        kbt=1
        ket=1
      endif

      if(iptra.eq.1)then
        ibp=ib
        iep=ie
        jbp=jb
        jep=je
        kbp=kb
        kep=ke
      else
        ibp=1
        iep=1
        jbp=1
        jep=1
        kbp=1
        kep=1
      endif

      if(psolver.eq.4.or.psolver.eq.5.or.ibalance.eq.2)then

        imirror = 0
        jmirror = 0

        ipb=1
        ipe=ni

        jpb=1
        jpe=nj

        if( (wbc.eq.2.or.wbc.eq.3).or.(ebc.eq.2.or.ebc.eq.3) )then

          imirror = 1
          ipe = ni*2

        endif

        if( (sbc.eq.2.or.sbc.eq.3).or.(nbc.eq.2.or.nbc.eq.3) )then

          jmirror = 1
          jpe = nj*2

        endif

        kpb=0
        kpe=nk+1

      else

        ipb=1
        ipe=1
        jpb=1
        jpe=1
        kpb=1
        kpe=1

      endif

!--------------------------------------------------------------

      rdx=1.0/dx
      rdy=1.0/dy
      rdz=1.0/dz
      rdx2=1.0/(2.0*dx)
      rdy2=1.0/(2.0*dy)
      rdz2=1.0/(2.0*dz)
      rdx4=1.0/(4.0*dx)
      rdy4=1.0/(4.0*dy)
      rdz4=1.0/(4.0*dz)

      thec_mb=0.0
      qt_mb=0.0

      stattim=statfrq
      taptim=tapfrq
      rsttim=rstfrq
      prcltim=prclfrq
      radtim=0.0

!--------------------------------------------------------------
!  Get identity

      ibw=0
      ibe=0
      ibs=0
      ibn=0

      patchsws = .false.
      patchsww = .false.
      patchses = .false.
      patchsee = .false.
      patchnwn = .false.
      patchnww = .false.
      patchnen = .false.
      patchnee = .false.

      p2tchsws = .false.
      p2tchsww = .false.
      p2tchses = .false.
      p2tchsee = .false.
      p2tchnwn = .false.
      p2tchnww = .false.
      p2tchnen = .false.
      p2tchnee = .false.

      myi=1
      myj=1

#ifdef MPI
      myj = myid / nodex + 1
      myi = myid - (myj-1)*nodex  + 1

      if(dowr) write(outfile,*) '  myi,myj=',myi,myj
      if(dowr) write(outfile,*)

      mynorth = nabor(myi,   myj+1, nodex, nodey)
      mysouth = nabor(myi,   myj-1, nodex, nodey)
      myeast  = nabor(myi+1, myj,   nodex, nodey)
      mywest  = nabor(myi-1, myj,   nodex, nodey)

      if(dowr) write(outfile,*) '  mywest  =',mywest
      if(dowr) write(outfile,*) '  myeast  =',myeast
      if(dowr) write(outfile,*) '  mysouth =',mysouth
      if(dowr) write(outfile,*) '  mynorth =',mynorth
      if(dowr) write(outfile,*)

      mysw = nabor(myi-1, myj-1,   nodex, nodey)
      mynw = nabor(myi-1, myj+1,   nodex, nodey)
      myne = nabor(myi+1, myj+1,   nodex, nodey)
      myse = nabor(myi+1, myj-1,   nodex, nodey)

      if(dowr) write(outfile,*) '  mysw    =',mysw
      if(dowr) write(outfile,*) '  mynw    =',mynw
      if(dowr) write(outfile,*) '  myne    =',myne
      if(dowr) write(outfile,*) '  myse    =',myse

      cs1we = (nj)*(nk)
      cs1sn = (ni)*(nk)
      ct1we = (nj)*(nk+1)
      ct1sn = (ni)*(nk+1)
      cv1we = (nj+1)*(nk)
      cu1sn = (ni+1)*(nk)
      cw1we = (nj)*(nk-1)
      cw1sn = (ni)*(nk-1)
      cs2we = 2*(nj)*(nk)
      cs2sn = (ni)*2*(nk)
      cs3we = 3*(nj)*(nk)
      cs3sn = (ni)*3*(nk)
      ct3we = 3*(nj)*(nk+1)
      ct3sn = (ni)*3*(nk+1)
      cv3we = 3*(nj+1)*(nk)
      cu3sn = 3*(ni+1)*(nk)
      cw3we = 3*(nj)*(nk-1)
      cw3sn = (ni)*3*(nk-1)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  cs1we   =',cs1we
      if(dowr) write(outfile,*) '  cs1sn   =',cs1sn
      if(dowr) write(outfile,*) '  ct1we   =',ct1we
      if(dowr) write(outfile,*) '  ct1sn   =',ct1sn
      if(dowr) write(outfile,*) '  cv1we   =',cv1we
      if(dowr) write(outfile,*) '  cu1sn   =',cu1sn
      if(dowr) write(outfile,*) '  cw1we   =',cw1we
      if(dowr) write(outfile,*) '  cw1sn   =',cw1sn
      if(dowr) write(outfile,*) '  cs2we   =',cs2we
      if(dowr) write(outfile,*) '  cs2sn   =',cs2sn
      if(dowr) write(outfile,*) '  cs3we   =',cs3we
      if(dowr) write(outfile,*) '  cs3sn   =',cs3sn
      if(dowr) write(outfile,*) '  ct3we   =',ct3we
      if(dowr) write(outfile,*) '  ct3sn   =',ct3sn
      if(dowr) write(outfile,*) '  cv3we   =',cv3we
      if(dowr) write(outfile,*) '  cu3sn   =',cu3sn
      if(dowr) write(outfile,*) '  cw3we   =',cw3we
      if(dowr) write(outfile,*) '  cw3sn   =',cw3sn
      if(dowr) write(outfile,*)
#endif

    IF(iorigin.eq.1)THEN

      do i=ib,ie
        xh(i)=dx*(i+(myi-1)*nx/nodex)-0.5*dx
        if(i+(myi-1)*nx/nodex.lt.1 .and. wbc.ne.1) ibw=1
        if(i+(myi-1)*nx/nodex.gt.nx .and. ebc.ne.1) ibe=1
      enddo

      do i=ib,ie+1
        xf(i)=dx*(i+(myi-1)*nx/nodex-1)
      enddo
      do i=-2,nx+4
        xfref(i)=dx*(i-1)
      enddo

      do j=jb,je
        yh(j)=dy*(j+(myj-1)*ny/nodey)-0.5*dy
        if(j+(myj-1)*ny/nodey.lt.1 .and. sbc.ne.1) ibs=1
        if(j+(myj-1)*ny/nodey.gt.ny .and. nbc.ne.1) ibn=1
      enddo

      do j=jb,je+1
        yf(j)=dy*(j+(myj-1)*ny/nodey-1)
      enddo
      do j=-2,ny+4
        yfref(j)=dy*(j-1)
      enddo

    ELSEIF(iorigin.eq.2)THEN

      do i=ib,ie
        xh(i)=dx*(i+(myi-1)*nx/nodex)-0.5*dx-0.5*dx*nx
        if(i+(myi-1)*nx/nodex.lt.1 .and. wbc.ne.1) ibw=1
        if(i+(myi-1)*nx/nodex.gt.nx .and. ebc.ne.1) ibe=1
      enddo

      do i=ib,ie+1
        xf(i)=dx*(i+(myi-1)*nx/nodex-1)-0.5*dx*nx
      enddo
      do i=-2,nx+4
        xfref(i)=dx*(i-1)-0.5*dx*nx
      enddo

      do j=jb,je
        yh(j)=dy*(j+(myj-1)*ny/nodey)-0.5*dy-0.5*dy*ny
        if(j+(myj-1)*ny/nodey.lt.1 .and. sbc.ne.1) ibs=1
        if(j+(myj-1)*ny/nodey.gt.ny .and. nbc.ne.1) ibn=1
      enddo

      do j=jb,je+1
        yf(j)=dy*(j+(myj-1)*ny/nodey-1)-0.5*dy*ny
      enddo
      do j=-2,ny+4
        yfref(j)=dy*(j-1)-0.5*dy*ny
      enddo

    ELSE

      print *,'  invalid option for iorigin'
      call stopcm1

    ENDIF

#ifndef MPI
      if(wbc.eq.2)then
        ibw=1
      endif

      if(ebc.eq.2)then
        ibe=1
      endif

      if(sbc.eq.2)then
        ibs=1
      endif

      if(nbc.eq.2)then
        ibn=1
      endif
#endif

!--------------------------------------------------------------

      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) 'g     =',g
      if(dowr) write(outfile,*) 'to    =',to
      if(dowr) write(outfile,*) 'rd    =',rd
      if(dowr) write(outfile,*) 'rv    =',rv
      if(dowr) write(outfile,*) 'cp    =',cp
      if(dowr) write(outfile,*) 'cv    =',cv
      if(dowr) write(outfile,*) 'cpv   =',cpv
      if(dowr) write(outfile,*) 'cvv   =',cvv
      if(dowr) write(outfile,*) 'p00   =',p00
      if(dowr) write(outfile,*) 'rp00  =',rp00
      if(dowr) write(outfile,*) 'rcp   =',rcp
      if(dowr) write(outfile,*) 'pi    =',pi

      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) 'cpdcv =',cpdcv
      if(dowr) write(outfile,*) 'rovcp =',rovcp
      if(dowr) write(outfile,*) 'rddcv =',rddcv
      if(dowr) write(outfile,*) 'cvdrd =',cvdrd
      if(dowr) write(outfile,*) 'cpdrd =',cpdrd
      if(dowr) write(outfile,*) 'eps   =',eps
      if(dowr) write(outfile,*) 'reps  =',reps
      if(dowr) write(outfile,*) 'repsm1=',repsm1
      if(dowr) write(outfile,*) 'cpt   =',cpt
      if(dowr) write(outfile,*) 'cvt   =',cvt
      if(dowr) write(outfile,*) 'pnum  =',pnum
      if(dowr) write(outfile,*) 'xlv   =',xlv
      if(dowr) write(outfile,*) 'xls   =',xls
      if(dowr) write(outfile,*) 'lvdcp =',lvdcp
      if(dowr) write(outfile,*) 'condc =',condc
      if(dowr) write(outfile,*) 'cpl   =',cpl
      if(dowr) write(outfile,*) 'cpi   =',cpi
      if(dowr) write(outfile,*) 'lv1   =',lv1
      if(dowr) write(outfile,*) 'lv2   =',lv2
      if(dowr) write(outfile,*) 'ls1   =',ls1
      if(dowr) write(outfile,*) 'ls2   =',ls2
      if(dowr) write(outfile,*) 'karman=',karman

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'timeformat   =',timeformat
      if(dowr) write(outfile,*) 'timestats    =',timestats
      if(dowr) write(outfile,*) 'terrain_flag =',terrain_flag
      if(dowr) write(outfile,*) 'procfiles    =',procfiles


      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'nx    =',nx
      if(dowr) write(outfile,*) 'ny    =',ny
      if(dowr) write(outfile,*) 'nz    =',nz

#ifdef MPI
      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) 'nodex =',nodex
      if(dowr) write(outfile,*) 'nodey =',nodey
#endif

      if(dowr) write(outfile,*)
 
      if(dowr) write(outfile,*) 'ni    =',ni
      if(dowr) write(outfile,*) 'nj    =',nj
      if(dowr) write(outfile,*) 'nk    =',nk
      if(dowr) write(outfile,*) 'nkp1  =',nkp1

      if(dowr) write(outfile,*)
 
      if(dowr) write(outfile,130) 'ib,ibm,ibi,ibc,ibt=',ib,ibm,ibi,ibc,ibt
      if(dowr) write(outfile,130) 'ie,iem,iei,iec,iet=',ie,iem,iei,iec,iet
      if(dowr) write(outfile,130) 'jb,jbm,jbi,jbc,jbt=',jb,jbm,jbi,jbc,jbt
      if(dowr) write(outfile,130) 'je,jem,jei,jec,jet=',je,jem,jei,jec,jet
      if(dowr) write(outfile,130) 'kb,kbm,kbi,kbc,kbt=',kb,kbm,kbi,kbc,kbt
      if(dowr) write(outfile,130) 'ke,kem,kei,kec,ket=',ke,kem,kei,kec,ket

130   format(1x,a19,5(4x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'imirror,jmirror = ',imirror,jmirror
      if(dowr) write(outfile,*)

      if(dowr) write(outfile,131) 'ibp,itb,ipb,ibr,ibb=',ibp,itb,ipb,ibr,ibb
      if(dowr) write(outfile,131) 'iep,ite,ipe,ier,ieb=',iep,ite,ipe,ier,ieb
      if(dowr) write(outfile,131) 'jbp,jtb,jpb,jbr,jbb=',jbp,jtb,jpb,jbr,jbb
      if(dowr) write(outfile,131) 'jep,jte,jpe,jer,jeb=',jep,jte,jpe,jer,jeb
      if(dowr) write(outfile,131) 'kbp,ktb,kpb,kbr,kbb=',kbp,ktb,kpb,kbr,kbb
      if(dowr) write(outfile,131) 'kep,kte,kpe,ker,keb=',kep,kte,kpe,ker,keb

131   format(1x,a20,5(4x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,132) 'ibl               =',ibl
      if(dowr) write(outfile,132) 'iel               =',iel
      if(dowr) write(outfile,132) 'jbl               =',jbl
      if(dowr) write(outfile,132) 'jel               =',jel

132   format(1x,a19,1(4x,i5))

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,133) imp,jmp,kmp,kmt,rmp,cmp

133   format(' imp,jmp,kmp,kmt,rmp,cmp =',6(2x,i5))

!----------

      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) 'rdx    =',rdx
      if(dowr) write(outfile,*) 'rdy    =',rdy
      if(dowr) write(outfile,*) 'rdz    =',rdz
      if(dowr) write(outfile,*) 'rdx2   =',rdx2
      if(dowr) write(outfile,*) 'rdy2   =',rdy2
      if(dowr) write(outfile,*) 'rdz2   =',rdz2
      if(dowr) write(outfile,*) 'rdx4   =',rdx4
      if(dowr) write(outfile,*) 'rdy4   =',rdy4
      if(dowr) write(outfile,*) 'rdz4   =',rdz4
      if(dowr) write(outfile,*) 'govtwo =',govtwo
      if(dowr) write(outfile,*) 'clwsat =',clwsat
      if(dowr) write(outfile,*) 'smeps  =',smeps
      if(dowr) write(outfile,*) 'tsmall =',tsmall
      if(dowr) write(outfile,*) 'cstar  =',cstar
      if(dowr) write(outfile,*) 'csmax  =',csmax
      if(dowr) write(outfile,*) 'epsilon=',epsilon

      if(dowr) write(outfile,*)

!--------------------------------------------------------------

      do i=ib,ie
        uh(i)=1.0
      enddo

      do i=ib,ie+1
        uf(i)=1.0
      enddo

      IF(stretch_x.ge.1)THEN

!!!        ibw=0
!!!        ibe=0

!-----------------------------------------------------------------------
!  Begin hard-wired analytic stretching function

        nominal_dx = 0.5*( dx_inner + dx_outer )

      IF(stretch_x.eq.1)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_x = 1 ... stretching on both west and east sides of domain:'
        if(dowr) write(outfile,*)
        ni1=(tot_x_len-nos_x_len)*0.5/nominal_dx
        ni2=nos_x_len/dx_inner
        ni3=ni1
        if(dowr) write(outfile,*) '  ni1,ni2,ni3 = ',(tot_x_len-nos_x_len)*0.5/nominal_dx,   &
                         nos_x_len/dx_inner,(tot_x_len-nos_x_len)*0.5/nominal_dx
        if(dowr) write(outfile,*) '    (note:  ni1,ni2,ni3 need to be exact integers for this to work correctly)'
      ELSEIF(stretch_x.eq.2)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_x = 2 ... stretching on east side of domain only:'
        if(dowr) write(outfile,*)
        ni1=0
        ni2=nos_x_len/dx_inner
        ni3=(tot_x_len-nos_x_len)/nominal_dx
        if(dowr) write(outfile,*) '  ni1,ni2,ni3 = ',0.0,nos_x_len/dx_inner,(tot_x_len-nos_x_len)/nominal_dx
        if(dowr) write(outfile,*) '    (note:  ni1,ni2,ni3 need to be exact integers for this to work correctly)'
      ELSE
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_x must be either 1 or 2'
        if(dowr) write(outfile,*)
        call stopcm1
      ENDIF

        c2=(nominal_dx-dx_inner)/(nominal_dx*nominal_dx*float(ni3-1))
        c1=(dx_inner/nominal_dx)-c2*nominal_dx

        if(dowr) write(outfile,*) '  nominal_dx  = ',nominal_dx
        if(dowr) write(outfile,*) '  c1,c2       = ',c1,c2
        if(dowr) write(outfile,*)

        ! Test to see if nx is kosher.
      IF(stretch_x.eq.1)THEN
        if(nx.ne.ni1+ni2+ni3)then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  There is a problem with the settings for horizontal grid stretching'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of nx = ',nx
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value needed for these settings ...'
          if(dowr) write(outfile,*) '       dx_inner  = ',dx_inner
          if(dowr) write(outfile,*) '       dx_outer  = ',dx_outer
          if(dowr) write(outfile,*) '       nos_x_len = ',nos_x_len
          if(dowr) write(outfile,*) '       tot_x_len = ',tot_x_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... would be nx = ',(nos_x_len/dx_inner)+(tot_x_len-nos_x_len)/(0.5*(dx_inner+dx_outer))
          if(dowr) write(outfile,*) '  (if this number is an integer) '
          if(dowr) write(outfile,*) '  (and if ni1,ni2,ni3 are all integers) '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
          call stopcm1
        endif
      ELSEIF(stretch_x.eq.2)THEN
        if(nx.ne.ni1+ni2+ni3)then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  There is a problem with the settings for horizontal grid stretching'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of nx = ',nx
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value for these settings ...'
          if(dowr) write(outfile,*) '       dx_inner  = ',dx_inner
          if(dowr) write(outfile,*) '       dx_outer  = ',dx_outer
          if(dowr) write(outfile,*) '       nos_x_len = ',nos_x_len
          if(dowr) write(outfile,*) '       tot_x_len = ',tot_x_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... would be nx = ',(nos_x_len/dx_inner)+(tot_x_len-nos_x_len)/(0.5*(dx_inner+dx_outer))
          if(dowr) write(outfile,*) '  (if this number is an integer) '
          if(dowr) write(outfile,*) '  (and if ni1,ni2,ni3 are all integers) '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
          call stopcm1
        endif
      ENDIF

        mult = 0.0
        if(iorigin.eq.2) mult = 0.5

      IF(stretch_x.eq.1)THEN

        do i=ni1+1,ni1+ni2+1
            xfref(i)=ni1*nominal_dx+(i-ni1-1)*dx_inner - mult*tot_x_len
        enddo
        do i=ni1+ni2+2,ni1+ni2+ni3+4
            xfref(i)=ni1*nominal_dx+(ni1+ni2+1-ni1-1)*dble(dx_inner)   &
                 +(c1+c2*dble(i-1-ni1-ni2)*nominal_dx)   &
                 *dble(i-1-ni1-ni2)*nominal_dx - mult*tot_x_len
        enddo
        do i=-2,ni1
            xfref(i)=ni1*nominal_dx+(ni1+1-ni1-1)*dble(dx_inner)    &
                 -(c1+c2*dble(ni1+1-i)*nominal_dx)   &
                 *dble(ni1+1-i)*nominal_dx - mult*tot_x_len
        enddo

      ELSEIF(stretch_x.eq.2)THEN

        do i=ni1+1,ni1+ni2+1
            xfref(i)=ni1*nominal_dx+(i-ni1-1)*dx_inner - mult*tot_x_len
        enddo
        do i=ni1+ni2+2,ni1+ni2+ni3+3
            xfref(i)=ni1*nominal_dx+(ni1+ni2+1-ni1-1)*dble(dx_inner)   &
                 +(c1+c2*dble(i-1-ni1-ni2)*nominal_dx)   &
                 *dble(i-1-ni1-ni2)*nominal_dx - mult*tot_x_len
        enddo
        do i=-2,ni1
            xfref(i)=ni1*nominal_dx+(ni1+1-ni1-1)*dble(dx_inner)    &
                 -(c1+c2*dble(ni1+1-i)*nominal_dx)   &
                 *dble(ni1+1-i)*nominal_dx - mult*tot_x_len
        enddo

      ENDIF

!!!        if( xf(ib).lt.0.0  .and. wbc.ne.1 ) ibw=1
!!!        if( xf(ie).gt.maxx .and. ebc.ne.1 ) ibe=1

        IF(stretch_x.eq.1)THEN
          xfref( 0)=xfref(1)-1*dx_outer
          xfref(-1)=xfref(1)-2*dx_outer
          xfref(-2)=xfref(1)-3*dx_outer
        ELSEIF(stretch_x.eq.2)THEN
          xfref( 0)=xfref(1)-1*dx_inner
          xfref(-1)=xfref(1)-2*dx_inner
          xfref(-2)=xfref(1)-3*dx_inner
        ENDIF

          xfref(nx+2)=xfref(nx+1)+1*dx_outer
          xfref(nx+3)=xfref(nx+1)+2*dx_outer
          xfref(nx+4)=xfref(nx+1)+3*dx_outer

!  End hard-wired analytic stretching function
!-----------------------------------------------------------------------
!
!  Optional:  to use a different stretching function, or to use 
!  arbitrarily located grid points, simply comment out the 
!  "hard-wired" section above, and then specify values for xfref
!  here.  Do not change anything below here!
!
!  Note:  xfref stores the location of the staggered u points for
!  the entire domain (from x=-2 to x=nx+4) (note: this includes
!  the boundary points that extend 3 gridpoints beyond the
!  computational domain.
!
!-----------------------------------------------------------------------

        do i=ib,ie+1
          xf(i)=xfref(i+(myi-1)*ni)
        enddo

        do i=ib,ie
          xh(i)=0.5*(xf(i+1)+xf(i))
          uh(i)=dx/(xf(i+1)-xf(i))
        enddo

        do i=ib+1,ie
          uf(i)=dx/(xh(i)-xh(i-1))
        enddo

        if(ibw.eq.1)then
          uf( 0)=uf(1)
          uf(-1)=uf(1)
          uf(-2)=uf(1)
        endif

        if(ibe.eq.1)then
          uf(ni+2)=uf(ni+1)
          uf(ni+3)=uf(ni+1)
          uf(ni+4)=uf(ni+1)
        endif

      ENDIF

      do i=ib,ie
        rxh(i)=1.0/(smeps+xh(i))
        ruh(i)=1.0/uh(i)
      enddo

      do i=ib,ie+1
        rxf(i)=1.0/(smeps+xf(i))
        ruf(i)=1.0/uf(i)
      enddo

      minx = xfref(1)
      maxx = xfref(nx+1)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'x:'
      if(dowr) write(outfile,124)
      if(dowr) write(outfile,125)
#ifdef MPI
    ! for MPI runs without procfiles, print the entire domain info:
    IF(.not.procfiles)THEN
      do i=1-ngxy,1-1
        if(dowr) write(outfile,122) i,xfref(i),0.5*(xfref(i)+xfref(i+1)),xfref(i+1)-xfref(i),dx/(0.5*(xfref(1)+xfref(2))-0.5*(xfref(0)+xfref(1))),dx/(xfref(i+1)-xfref(i)),'   x'
      enddo
      do i=1,nx
        if(dowr) write(outfile,122) i,xfref(i),0.5*(xfref(i)+xfref(i+1)),xfref(i+1)-xfref(i),dx/(0.5*(xfref(i)+xfref(i+1))-0.5*(xfref(i-1)+xfref(i))),dx/(xfref(i+1)-xfref(i)),'    '
      enddo
      do i=nx+1,nx+ngxy
        if(dowr) write(outfile,122) i,xfref(i),0.5*(xfref(i)+xfref(i+1)),xfref(i+1)-xfref(i),dx/(0.5*(xfref(nx+1)+xfref(nx+2))-0.5*(xfref(nx)+xfref(nx+1))),dx/(xfref(i+1)-xfref(i)),'   x'
      enddo
      if(dowr) write(outfile,123) nx+1+ngxy,xfref(nx+1+ngxy),dx/(0.5*(xfref(nx+1)+xfref(nx+2))-0.5*(xfref(nx)+xfref(nx+1)))
    ELSE
#endif
      do i=ib,ib+2
        if(dowr) write(outfile,122) i,xf(i),xh(i),xf(i+1)-xf(i),uf(i),uh(i),'   x'
      enddo
      do i=ib+3,ie-3
        if(dowr) write(outfile,122) i,xf(i),xh(i),xf(i+1)-xf(i),uf(i),uh(i),'    '
      enddo
      do i=ie-2,ie
        if(dowr) write(outfile,122) i,xf(i),xh(i),xf(i+1)-xf(i),uf(i),uh(i),'   x'
      enddo
      if(dowr) write(outfile,123) ie+1,xf(ie+1),uf(ie+1)
#ifdef MPI
    ENDIF
#endif
      if(dowr) write(outfile,*)

122   format(3x,i5,3x,f11.2,3x,f11.2,3x,f9.2,3x,f8.4,3x,f8.4,a4)
123   format(3x,i5,3x,f11.2,29x,f8.4)
124   format('      i         xf           xh         dx         uf         uh')
125   format(' ---------------------------------------------------------------')

!--------------------------------------------------------------

      do j=jb,je
        vh(j)=1.0
      enddo

      do j=jb,je+1
        vf(j)=1.0
      enddo

      IF(stretch_y.ge.1)THEN

!!!        ibs=0
!!!        ibn=0

!-----------------------------------------------------------------------
!  Begin hard-wired analytic stretching function

        nominal_dy = 0.5*( dy_inner + dy_outer )

      IF(stretch_y.eq.1)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_y = 1 ... stretching on both south and north sides of domain:'
        if(dowr) write(outfile,*)
        nj1=(tot_y_len-nos_y_len)*0.5/nominal_dy
        nj2=nos_y_len/dy_inner
        nj3=nj1
        if(dowr) write(outfile,*) '  nj1,nj2,nj3 = ',(tot_y_len-nos_y_len)*0.5/nominal_dy,   &
                         nos_y_len/dy_inner,(tot_y_len-nos_y_len)*0.5/nominal_dy
        if(dowr) write(outfile,*) '    (note:  nj1,nj2,nj3 need to be exact integers for this to work correctly)'
      ELSEIF(stretch_y.eq.2)THEN
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_y = 2 ... stretching on north side of domain only:'
        if(dowr) write(outfile,*)
        nj1=0
        nj2=nos_y_len/dy_inner
        nj3=(tot_y_len-nos_y_len)/nominal_dy
        if(dowr) write(outfile,*) '  nj1,nj2,nj3 = ',0.0,nos_y_len/dy_inner,(tot_y_len-nos_y_len)/nominal_dy
        if(dowr) write(outfile,*) '    (note:  nj1,nj2,nj3 need to be exact integers for this to work correctly)'
      ELSE
        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) ' stretch_y must be either 1 or 2'
        if(dowr) write(outfile,*)
        call stopcm1
      ENDIF

        c2=(nominal_dy-dy_inner)/(nominal_dy*nominal_dy*float(nj3-1))
        c1=(dy_inner/nominal_dy)-c2*nominal_dy

        if(dowr) write(outfile,*) '  nominal_dy  = ',nominal_dy
        if(dowr) write(outfile,*) '  c1,c2       = ',c1,c2
        if(dowr) write(outfile,*)

        ! Test to see if ny is kosher.
      IF(stretch_y.eq.1)THEN
        if(ny.ne.nj1+nj2+nj3)then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  There is a problem with the settings for horizontal grid stretching'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of ny = ',ny
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value needed for these settings ...'
          if(dowr) write(outfile,*) '       dy_inner  = ',dy_inner
          if(dowr) write(outfile,*) '       dy_outer  = ',dy_outer
          if(dowr) write(outfile,*) '       nos_y_len = ',nos_y_len
          if(dowr) write(outfile,*) '       tot_y_len = ',tot_y_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... would be ny = ',(nos_y_len/dy_inner)+(tot_y_len-nos_y_len)/(0.5*(dy_inner+dy_outer))
          if(dowr) write(outfile,*) '  (if this number is an integer) '
          if(dowr) write(outfile,*) '  (and if nj1,nj2,nj3 are all integers) '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
          call stopcm1
        endif
      ELSEIF(stretch_y.eq.2)THEN
        if(ny.ne.nj1+nj2+nj3)then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  There is a problem with the settings for horizontal grid stretching'
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of ny = ',ny
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value for these settings ...'
          if(dowr) write(outfile,*) '       dy_inner  = ',dy_inner
          if(dowr) write(outfile,*) '       dy_outer  = ',dy_outer
          if(dowr) write(outfile,*) '       nos_y_len = ',nos_y_len
          if(dowr) write(outfile,*) '       tot_y_len = ',tot_y_len
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... would be ny = ',(nos_y_len/dy_inner)+(tot_y_len-nos_y_len)/(0.5*(dy_inner+dy_outer))
          if(dowr) write(outfile,*) '  (if this number is an integer) '
          if(dowr) write(outfile,*) '  (and if nj1,nj2,nj3 are all integers) '
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
          call stopcm1
        endif
      ENDIF

        mult = 0.0
        if(iorigin.eq.2) mult = 0.5

      IF(stretch_y.eq.1)THEN

        do j=nj1+1,nj1+nj2+1
            yfref(j)=nj1*nominal_dy+(j-nj1-1)*dy_inner - mult*tot_y_len
        enddo
        do j=nj1+nj2+2,nj1+nj2+nj3+4
            yfref(j)=nj1*nominal_dy+(nj1+nj2+1-nj1-1)*dble(dy_inner)   &
                 +(c1+c2*dble(j-1-nj1-nj2)*nominal_dy)   &
                 *dble(j-1-nj1-nj2)*nominal_dy - mult*tot_y_len
        enddo
        do j=-2,nj1
            yfref(j)=nj1*nominal_dy+(nj1+1-nj1-1)*dble(dy_inner)    &
                 -(c1+c2*dble(nj1+1-j)*nominal_dy)   &
                 *dble(nj1+1-j)*nominal_dy - mult*tot_y_len
        enddo

      ELSEIF(stretch_y.eq.2)THEN

        do j=nj1+1,nj1+nj2+1
            yfref(j)=nj1*nominal_dy+(j-nj1-1)*dy_inner - mult*tot_y_len
        enddo
        do j=nj1+nj2+2,nj1+nj2+nj3+3
            yfref(j)=nj1*nominal_dy+(nj1+nj2+1-nj1-1)*dble(dy_inner)   &
                 +(c1+c2*dble(j-1-nj1-nj2)*nominal_dy)   &
                 *dble(j-1-nj1-nj2)*nominal_dy - mult*tot_y_len
        enddo
        do j=-2,nj1
            yfref(j)=nj1*nominal_dy+(nj1+1-nj1-1)*dble(dy_inner)    &
                 -(c1+c2*dble(nj1+1-j)*nominal_dy)   &
                 *dble(nj1+1-j)*nominal_dy - mult*tot_y_len
        enddo

      ENDIF

!!!        if( yf(jb).lt.0.0  .and. sbc.ne.1 ) ibs=1
!!!        if( yf(je).gt.maxy .and. nbc.ne.1 ) ibn=1

        IF(stretch_y.eq.1)THEN
          yfref( 0)=yfref(1)-1*dy_outer
          yfref(-1)=yfref(1)-2*dy_outer
          yfref(-2)=yfref(1)-3*dy_outer
        ELSEIF(stretch_y.eq.2)THEN
          yfref( 0)=yfref(1)-1*dy_inner
          yfref(-1)=yfref(1)-2*dy_inner
          yfref(-2)=yfref(1)-3*dy_inner
        ENDIF

          yfref(ny+2)=yfref(ny+1)+1*dy_outer
          yfref(ny+3)=yfref(ny+1)+2*dy_outer
          yfref(ny+4)=yfref(ny+1)+3*dy_outer

!  End hard-wired analytic stretching function
!-----------------------------------------------------------------------
!
!  Optional:  to use a different stretching function, or to use 
!  arbitrarily located grid points, simply comment out the 
!  "hard-wired" section above, and then specify values for yfref
!  here.  Do not change anything below here!
!
!  Note:  yfref stores the location of the staggered v points for
!  the entire domain (from y=-2 to y=ny+4) (note: this includes
!  the boundary points that extend 3 gridpoints beyond the
!  computational domain.
!
!-----------------------------------------------------------------------

        do j=jb,je+1
          yf(j)=yfref(j+(myj-1)*nj)
        enddo

        do j=jb,je
          yh(j)=0.5*(yf(j+1)+yf(j))
          vh(j)=dy/(yf(j+1)-yf(j))
        enddo

        do j=jb+1,je
          vf(j)=dy/(yh(j)-yh(j-1))
        enddo

        if(ibs.eq.1)then
          vf( 0)=vf(1)
          vf(-1)=vf(1)
          vf(-2)=vf(1)
        endif

        if(ibn.eq.1)then
          vf(nj+2)=vf(nj+1)
          vf(nj+3)=vf(nj+1)
          vf(nj+4)=vf(nj+1)
        endif

      ENDIF

      do j=jb,je
        rvh(j)=1.0/vh(j)
      enddo

      do j=jb,je+1
        rvf(j)=1.0/vf(j)
      enddo

      miny = yfref(1)
      maxy = yfref(ny+1)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'y:'
      if(dowr) write(outfile,134)
      if(dowr) write(outfile,125)
#ifdef MPI
    ! for MPI runs without procfiles, print the entire domain info:
    IF(.not.procfiles)THEN
      do j=1-ngxy,1-1
        if(dowr) write(outfile,122) j,yfref(j),0.5*(yfref(j)+yfref(j+1)),yfref(j+1)-yfref(j),dy/(0.5*(yfref(1)+yfref(2))-0.5*(yfref(0)+yfref(1))),dy/(yfref(j+1)-yfref(j)),'   x'
      enddo
      do j=1,ny
        if(dowr) write(outfile,122) j,yfref(j),0.5*(yfref(j)+yfref(j+1)),yfref(j+1)-yfref(j),dy/(0.5*(yfref(j)+yfref(j+1))-0.5*(yfref(j-1)+yfref(j))),dy/(yfref(j+1)-yfref(j)),'    '
      enddo
      do j=ny+1,ny+ngxy
        if(dowr) write(outfile,122) j,yfref(j),0.5*(yfref(j)+yfref(j+1)),yfref(j+1)-yfref(j),dy/(0.5*(yfref(ny+1)+yfref(ny+2))-0.5*(yfref(ny)+yfref(ny+1))),dy/(yfref(j+1)-yfref(j)),'   x'
      enddo
      if(dowr) write(outfile,123) ny+1+ngxy,yfref(ny+1+ngxy),dy/(0.5*(yfref(ny+1)+yfref(ny+2))-0.5*(yfref(ny)+yfref(ny+1)))
    ELSE
#endif
      do j=jb,jb+2
        if(dowr) write(outfile,122) j,yf(j),yh(j),yf(j+1)-yf(j),vf(j),vh(j),'   x'
      enddo
      do j=jb+3,je-3
        if(dowr) write(outfile,122) j,yf(j),yh(j),yf(j+1)-yf(j),vf(j),vh(j),'    '
      enddo
      do j=je-2,je
        if(dowr) write(outfile,122) j,yf(j),yh(j),yf(j+1)-yf(j),vf(j),vh(j),'   x'
      enddo
      if(dowr) write(outfile,123) je+1,yf(je+1),vf(je+1)
#ifdef MPI
    ENDIF
#endif
      if(dowr) write(outfile,*)

134   format('      j         yf           yh         dy         vf         vh')

!--------------------------------------------------------------

      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        zf(i,j,k)=dz*(k-1)
        mf(i,j,k)=1.0
      enddo
      enddo
      enddo

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        zh(i,j,k)=0.5*(zf(i,j,k)+zf(i,j,k+1))
        mh(i,j,k)=1.0
      enddo
      enddo
      enddo

      do k=kb,ke+1
        sigmaf(k)=zf(1,1,k)
      enddo
      do k=kb,ke
        sigma(k)=0.5*(sigmaf(k)+sigmaf(k+1))
      enddo

    IF(stretch_z.eq.1)THEN

!-----------------------------------------------------------------------
!  Begin hard-wired analytic stretching function

        nominal_dz = 0.5*(dz_bot+dz_top)

        ! dummy checks:
        if(dowr) write(outfile,*) '  mod bot:',mod(str_bot,dz_bot),1.0e-6*str_bot
        if(mod(str_bot,dz_bot).gt.1.0e-6*str_bot)then
          print *,'  depth of bottom layer does not exactly divide by dz_bot! '
          print *,'  str_bot = ',str_bot
          print *,'  dz_bot  = ',dz_bot
          print *,'  mod     = ',mod(str_bot,dz_bot)
          print *,'  stopping cm1 ... '
          call stopcm1
        endif
        if(dowr) write(outfile,*) '  mod mid:',mod(str_top-str_bot,nominal_dz),1.0e-6*(str_top-str_bot)
        if(mod(str_top-str_bot,nominal_dz).gt.1.0e-6*(str_top-str_bot))then
          print *,'  depth of middle layer does not exactly divide by nominal_dz! '
          print *,'  str_bot         = ',str_bot
          print *,'  str_top         = ',str_top
          print *,'  str_top-str_bot = ',str_top-str_bot
          print *,'  nominal_dz      = ',nominal_dz
          print *,'  mod             = ',mod(str_top-str_bot,nominal_dz)
          print *,'  stopping cm1 ... '
          call stopcm1
        endif
        if(dowr) write(outfile,*) '  mod top:',mod(ztop-str_top,dz_top),1.0e-6*(ztop-str_top)
        if(mod(ztop-str_top,dz_top).gt.1.0e-6*(ztop-str_top))then
          print *,'  depth of top layer does not exactly divide by dz_top! '
          print *,'  str_top      = ',str_top
          print *,'  ztop         = ',ztop
          print *,'  ztop-str_top = ',ztop-str_top
          print *,'  dz_top       = ',dz_top
          print *,'  mod          = ',mod(ztop-str_top,dz_top)
          print *,'  stopping cm1 ... '
          call stopcm1
        endif

        nk1=str_bot/dz_bot
        nk3=(ztop-str_top)/dz_top
        nk2=nk-(nk1+nk3)

        nominal_dz=(str_top-str_bot)/nk2

        c2=(nominal_dz-dz_bot)/(nominal_dz*nominal_dz*float(nk2-1))
        c1=(dz_bot/nominal_dz)-c2*nominal_dz

        ! Test to see if nk is kosher.
        if(dowr) write(outfile,*) float(nk),float(nk1+nk3)+(str_top-str_bot)/(0.5*(dz_bot+dz_top))
        if( abs(float(nk)-(float(nk1+nk3)+(str_top-str_bot)/(0.5*(dz_bot+dz_top)))).gt.1.0e-3 )then
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  User value of nz = ',nz
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  Value needed for these settings:'
          if(dowr) write(outfile,*) '       ztop      = ',ztop
          if(dowr) write(outfile,*) '       str_bot   = ',str_bot
          if(dowr) write(outfile,*) '       str_top   = ',str_top
          if(dowr) write(outfile,*) '       dz_bot    = ',dz_bot
          if(dowr) write(outfile,*) '       dz_top    = ',dz_top
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  would be nz = ',nk1+nk3+(str_top-str_bot)/(0.5*(dz_bot+dz_top))
          if(dowr) write(outfile,*)
          if(dowr) write(outfile,*) '  ... stopping ...  '
          if(dowr) write(outfile,*)
          call stopcm1
        endif

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  nk1,nk2,nk3,ntot=',nk1,nk2,nk3,(nk1+nk2+nk3)
        if(dowr) write(outfile,*) '  nominal_dz =',nominal_dz
        if(dowr) write(outfile,*) '  c1,c2 = ',c1,c2
        if(dowr) write(outfile,*)

      do j=jb,je
      do i=ib,ie

        do k=1,nk1+1
          zf(i,j,k)=(k-1)*dz_bot
        enddo
        do k=(nk1+1),(nk1+nk2+1)
          zf(i,j,k)=zf(i,j,nk1+1)+(c1+c2*float(k-1-nk1)*nominal_dz)   &
                         *float(k-1-nk1)*nominal_dz
        enddo
        do k=(nk1+nk2+2),(nk1+nk2+nk3+1)
          zf(i,j,k)=zf(i,j,k-1)+dz_top
        enddo

      enddo
      enddo

!!!      if(terrain_flag)then

        do k=1,nk1+1
          sigmaf(k)=(k-1)*dz_bot
        enddo
        do k=(nk1+1),(nk1+nk2+1)
          sigmaf(k)=sigmaf(nk1+1)+(c1+c2*float(k-1-nk1)*nominal_dz)   &
                         *float(k-1-nk1)*nominal_dz
        enddo
        do k=(nk1+nk2+2),(nk1+nk2+nk3+1)
          sigmaf(k)=sigmaf(k-1)+dz_top
        enddo

        sigmaf(0)=-sigmaf(2)
        sigmaf(nk+2)=sigmaf(nk+1)+(sigmaf(nk+1)-sigmaf(nk))

!!!      endif

!  End hard-wired analytic stretching function
!-----------------------------------------------------------------------
!
!  Optional:  to use a different stretching function, or to use 
!  arbitrarily located grid points, simply comment out the 
!  "hard-wired" section above, and then specify values for zf
!  here.  Do not change anything below here!
!
!  Note:  zf stores the location of the staggered w points. 
!
!  Note:  if you are using terrain, you need to also specify the nominal 
!  locations of the zf points in the sigmaf array.
!
!-----------------------------------------------------------------------

      do j=jb,je
      do i=ib,ie

        zf(i,j,0)=-zf(i,j,2)
        zf(i,j,nk+2)=zf(i,j,nk+1)+(zf(i,j,nk+1)-zf(i,j,nk))

        do k=0,nk+1
          zh(i,j,k)=0.5*(zf(i,j,k+1)+zf(i,j,k))
          mh(i,j,k)=dz/(zf(i,j,k+1)-zf(i,j,k))
        enddo
        zh(i,j,0)=-zh(i,j,1)
        zh(i,j,nk+1)=zh(i,j,nk)+2.0*(zf(i,j,nk+1)-zh(i,j,nk))

        do k=1,nk+1
          mf(i,j,k)=dz/(zh(i,j,k)-zh(i,j,k-1))
        enddo
        mf(i,j,0)=mf(i,j,1)
        mf(i,j,nk+2)=mf(i,j,nk+1)

      enddo
      enddo

    ENDIF

! end vertical stretching section
!-----------------------------------------------------------------------

      do k=kb,ke
        sigma(k)=0.5*(sigmaf(k)+sigmaf(k+1))
      enddo

      maxz = sigmaf(nk+1)

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        rmh(i,j,k)=1.0/mh(i,j,k)
      enddo
      enddo
      enddo

      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        rmf(i,j,k)=1.0/mf(i,j,k)
      enddo
      enddo
      enddo

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'model heights:'
      if(dowr) write(outfile,104)
104   format('     k       zf         zh         dz         mf         mh')
      if(dowr) write(outfile,105)
105   format(' ---------------------------------------------------------------')
      do k=1,nk
        if(dowr) write(outfile,102) k,zf(1,1,k),zh(1,1,k),zf(1,1,k+1)-zf(1,1,k),mf(1,1,k),mh(1,1,k)
102     format(3x,i4,3x,f8.2,3x,f8.2,3x,f8.2,3x,f8.4,3x,f8.4)
      enddo
      if(dowr) write(outfile,103) nk+1,zf(1,1,nk+1),mf(1,1,nk+1)
103   format(3x,i4,3x,f8.2,25x,f8.4)
      if(dowr) write(outfile,*)

!-----------------------------------------------------------------------

#ifdef MPI
      if( ibw.eq.1 .and. yh(jb).gt.miny ) patchsww = .true.
      if( ibw.eq.1 .and. yh(je).lt.maxy ) patchnww = .true.
      if( ibe.eq.1 .and. yh(jb).gt.miny ) patchsee = .true.
      if( ibe.eq.1 .and. yh(je).lt.maxy ) patchnee = .true.
      if( ibs.eq.1 .and. xh(ib).gt.minx ) patchsws = .true.
      if( ibs.eq.1 .and. xh(ie).lt.maxx ) patchses = .true.
      if( ibn.eq.1 .and. xh(ib).gt.minx ) patchnwn = .true.
      if( ibn.eq.1 .and. xh(ie).lt.maxx ) patchnen = .true.

      if(dowr) write(outfile,*) '  patchsww =',patchsww
      if(dowr) write(outfile,*) '  patchnww =',patchnww
      if(dowr) write(outfile,*) '  patchsee =',patchsee
      if(dowr) write(outfile,*) '  patchnee =',patchnee
      if(dowr) write(outfile,*) '  patchsws =',patchsws
      if(dowr) write(outfile,*) '  patchses =',patchses
      if(dowr) write(outfile,*) '  patchnwn =',patchnwn
      if(dowr) write(outfile,*) '  patchnen =',patchnen
      if(dowr) write(outfile,*)

      if( ibw.eq.1 .and. yh(jb).lt.miny ) p2tchsww = .true.
      if( ibw.eq.1 .and. yh(je).gt.maxy ) p2tchnww = .true.
      if( ibe.eq.1 .and. yh(jb).lt.miny ) p2tchsee = .true.
      if( ibe.eq.1 .and. yh(je).gt.maxy ) p2tchnee = .true.
      if( ibs.eq.1 .and. xh(ib).lt.minx ) p2tchsws = .true.
      if( ibs.eq.1 .and. xh(ie).gt.maxx ) p2tchses = .true.
      if( ibn.eq.1 .and. xh(ib).lt.minx ) p2tchnwn = .true.
      if( ibn.eq.1 .and. xh(ie).gt.maxx ) p2tchnen = .true.

      if(dowr) write(outfile,*) '  p2tchsww =',p2tchsww
      if(dowr) write(outfile,*) '  p2tchnww =',p2tchnww
      if(dowr) write(outfile,*) '  p2tchsee =',p2tchsee
      if(dowr) write(outfile,*) '  p2tchnee =',p2tchnee
      if(dowr) write(outfile,*) '  p2tchsws =',p2tchsws
      if(dowr) write(outfile,*) '  p2tchses =',p2tchses
      if(dowr) write(outfile,*) '  p2tchnwn =',p2tchnwn
      if(dowr) write(outfile,*) '  p2tchnen =',p2tchnen
      if(dowr) write(outfile,*)
#endif

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                  BEGIN TERRAIN !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      rds  = 1.0
      rdsf = 1.0

      do j=jtb,jte
      do i=itb,ite
        zs(i,j)=0.0
        gz(i,j)=1.0
        rgz(i,j)=1.0
        gzu(i,j)=1.0
        rgzu(i,j)=1.0
        gzv(i,j)=1.0
        rgzv(i,j)=1.0
        dzdx(i,j)=0.0
        dzdy(i,j)=0.0
      enddo
      enddo

      do k=ktb,kte
      do j=jtb,jte
      do i=itb,ite
        gx(i,j,k)=0.0
        gxu(i,j,k)=0.0
        gy(i,j,k)=0.0
        gyv(i,j,k)=0.0
      enddo
      enddo
      enddo


      IF(terrain_flag)THEN

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  Terrain included!'
        if(dowr) write(outfile,*)

        zs = 0.0

        ! moved this section of code to init_terrain in cm1r15:
        call init_terrain(xh,uh,xf,uf,yh,vh,yf,vf,rds,sigma,rdsf,sigmaf,  &
                          zh,zf,zs,gz,rgz,gzu,rgzu,gzv,rgzv,         &
                          dzdx,dzdy,gx,gxu,gy,gyv,                   &
                          reqs_u,reqs_v,reqs_s,reqs_p,               &
                          nw1,nw2,ne1,ne2,sw1,sw2,se1,se2,           &
                          sw31,sw32,se31,se32,ss31,ss32,sn31,sn32,   &
                          uw31,uw32,ue31,ue32,us31,us32,un31,un32,   &
                          vw31,vw32,ve31,ve32,vs31,vs32,vn31,vn32,   &
                          ww31(1,1,1),ww32(1,1,1),we31(1,1,1),we32(1,1,1), &
                          ws31(1,1,1),ws32(1,1,1),wn31(1,1,1),wn32(1,1,1))

      ENDIF

!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
!                  END   TERRAIN !
!CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        mh(i,j,k)=dz/(zf(i,j,k+1)-zf(i,j,k))
        rmh(i,j,k)=1.0/mh(i,j,k)
      enddo
      enddo
      enddo

      do k=kb+1,ke
      do j=jb,je
      do i=ib,ie
        mf(i,j,k)=dz/(zh(i,j,k)-zh(i,j,k-1))
      enddo
      enddo
      enddo

      do j=jb,je
      do i=ib,ie
        mf(i,j,0)=mf(i,j,1)
        mf(i,j,nk+2)=mf(i,j,nk+1)
      enddo
      enddo

      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        rmf(i,j,k)=1.0/mf(i,j,k)
      enddo
      enddo
      enddo

!-----------------------------------------------------------------------

      if(dowr) write(outfile,*) '  minx = ',minx
      if(dowr) write(outfile,*) '  maxx = ',maxx
      if(dowr) write(outfile,*) '  miny = ',miny
      if(dowr) write(outfile,*) '  maxy = ',maxy
      if(dowr) write(outfile,*) '  maxz = ',maxz
      if(dowr) write(outfile,*)

      if(dowr) write(outfile,*) '  ibw =',ibw
      if(dowr) write(outfile,*) '  ibe =',ibe
      if(dowr) write(outfile,*) '  ibs =',ibs
      if(dowr) write(outfile,*) '  ibn =',ibn
      if(dowr) write(outfile,*)

!-----------------------------------------------------------------------
!  Get min/max dx,dy,dz on grid
!  (needed for adapt_dt ... but interesting to report, nontheless)

      min_dx = 1.0e20
      min_dy = 1.0e20
      min_dz = 1.0e20

      max_dx = 0.0
      max_dy = 0.0
      max_dz = 0.0

      do i=1,ni
        min_dx = min( min_dx , xf(i+1)-xf(i) )
        max_dx = max( max_dx , xf(i+1)-xf(i) )
      enddo

      do j=1,nj
        min_dy = min( min_dy , yf(j+1)-yf(j) )
        max_dy = max( max_dy , yf(j+1)-yf(j) )
      enddo

      do k=1,nk
      do j=1,nj
      do i=1,ni
        min_dz = min( min_dz , zf(i,j,k+1)-zf(i,j,k) )
        max_dz = max( max_dz , zf(i,j,k+1)-zf(i,j,k) )
      enddo
      enddo
      enddo

#ifdef MPI
      var=0.0
      call MPI_ALLREDUCE(min_dx,var,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
      min_dx=var
      var=0.0
      call MPI_ALLREDUCE(min_dy,var,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
      min_dy=var
      var=0.0
      call MPI_ALLREDUCE(min_dz,var,1,MPI_REAL,MPI_MIN,MPI_COMM_WORLD,ierr)
      min_dz=var
      var=0.0
      call MPI_ALLREDUCE(max_dx,var,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      max_dx=var
      var=0.0
      call MPI_ALLREDUCE(max_dy,var,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      max_dy=var
      var=0.0
      call MPI_ALLREDUCE(max_dz,var,1,MPI_REAL,MPI_MAX,MPI_COMM_WORLD,ierr)
      max_dz=var
#endif

      if(dowr) write(outfile,*) '  min_dx = ',min_dx
      if(dowr) write(outfile,*) '  max_dx = ',max_dx
      if(dowr) write(outfile,*) '  min_dy = ',min_dy
      if(dowr) write(outfile,*) '  max_dy = ',max_dy
      if(dowr) write(outfile,*) '  min_dz = ',min_dz
      if(dowr) write(outfile,*) '  max_dz = ',max_dz
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  new (cm1r16) arrays for vertical interpolation:

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  k,c1,c2:'
      do k=1,nk+1
      do j=jb,je
      do i=ib,ie
        cc2(i,j,k)=(zf(i,j,k)-zh(i,j,k-1))/(zh(i,j,k)-zh(i,j,k-1))
        cc1(i,j,k)=1.0-cc2(i,j,k)
        if(i.eq.1.and.j.eq.1.and.dowr) write(outfile,141) k,cc1(i,j,k),cc2(i,j,k),zh(i,j,k-1),zf(i,j,k),zh(i,j,k)
141     format(3x,i4,2(3x,f6.4),3(3x,f8.2))
      enddo
      enddo
      enddo
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  Specify coefficient for Rayleigh damper in vertical

      do k=kb,ke
      do j=jb,je
      do i=ib,ie
        tauh(i,j,k)=0.0
        taus(i,j,k)=0.0
      enddo
      enddo
      enddo

      do k=kb,ke+1
      do j=jb,je
      do i=ib,ie
        tauf(i,j,k)=0.0
      enddo
      enddo
      enddo

      if( (irdamp.eq.1).and.(zd.lt.maxz) )then

        do j=jb,je
        do i=ib,ie
          do k=1,nk
            if(zh(i,j,k).gt.zd)then
            tauh(i,j,k)=0.5*(1.0-cos(pi*(zh(i,j,k)-zd)/(zf(i,j,nk+1)-zd)))
            taus(i,j,k)=tauh(i,j,k)
            endif
          enddo
          enddo
        enddo
 
        do j=jb,je
        do i=ib,ie
          do k=1,nk+1
            if(zf(i,j,k).gt.zd)then
            tauf(i,j,k)=0.5*(1.0-cos(pi*(zf(i,j,k)-zd)/(zf(i,j,nk+1)-zd)))
            endif
          enddo
          enddo
        enddo

      endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  ------ tauf, tauh -----'
      do k=1,nk
        if(dowr) write(outfile,*) k,tauf(1,1,k),tauh(1,1,k)
      enddo
      if(dowr) write(outfile,*) nk+1,tauf(1,1,nk+1)
      if(dowr) write(outfile,*)

!--------------------------------------------------------------
!  Rayleigh damping near lateral boundaries:

      IF(hrdamp.eq.1)THEN

        do k=1,nk
        do j=jb,je
        do i=ib,ie
          ! skip this section of code for 2d simulations:
          IF(nx.gt.1)THEN
            ! west boundary:
            IF( axisymm.ne.1 )THEN
              x1 = xhd+minx-xh(i)
              if( x1.gt.0.0 )then
                tauh(i,j,k) = max( tauh(i,j,k) , 0.5*(1.0-cos(pi*x1/xhd)) )
                tauf(i,j,k) = max( tauf(i,j,k) , 0.5*(1.0-cos(pi*x1/xhd)) )
              endif
            ENDIF
            ! east boundary:
            x2 = xh(i)-(maxx-xhd)
            if( x2.gt.0.0 )then
              tauh(i,j,k) = max( tauh(i,j,k) , 0.5*(1.0-cos(pi*x2/xhd)) )
              tauf(i,j,k) = max( tauf(i,j,k) , 0.5*(1.0-cos(pi*x2/xhd)) )
            endif
          ENDIF
          ! skip this section of code for 2d simulations:
          IF(ny.gt.1)THEN
            ! south boundary:
            y1 = xhd+miny-yh(j)
            if( y1.gt.0.0 )then
              tauh(i,j,k) = max( tauh(i,j,k) , 0.5*(1.0-cos(pi*y1/xhd)) )
              tauf(i,j,k) = max( tauf(i,j,k) , 0.5*(1.0-cos(pi*y1/xhd)) )
            endif
            ! north boundary:
            y2 = yh(j)-(maxy-xhd)
            if( y2.gt.0.0 )then
              tauh(i,j,k) = max( tauh(i,j,k) , 0.5*(1.0-cos(pi*y2/xhd)) )
              tauf(i,j,k) = max( tauf(i,j,k) , 0.5*(1.0-cos(pi*y2/xhd)) )
            endif
          ENDIF
        enddo
        enddo
        enddo

      ENDIF

!--------------------------------------------------------------
!  vertically implicit turbulent diffusion:

      ! Set vialpha:
      !      0.0 = forward-in-time (unstable if K dt / (dz^2) > 0.5)
      !      0.5 = centered-in-time (Crank-Nicholson) (stable but oscillatory)
      !      1.0 = backward-in-time (stable)
!      vialpha = 1.0

      ! Do not change this:
!      vibeta  = 1.0 - vialpha

!      NOTE:  these are now set in constants.incl files

        if(dowr) write(outfile,*)
        if(dowr) write(outfile,*) '  vialpha,vibeta = ',vialpha,vibeta
        if(dowr) write(outfile,*)

!--------------------------------------------------------------

      dt = dtl
      dtlast = dt

!--------------------------------------------------------------

      if(dowr) write(outfile,*) 'Leaving PARAM'

      return
      end


