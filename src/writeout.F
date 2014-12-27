

      subroutine setup_output(tdef,qname,budname,xh,xf,yh,yf,xfref,yfref,sigma,sigmaf,zh,zf)
      implicit none

      include 'input.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      !------------------------------------------------------------
      ! This subroutine gets things ready for writeouts
      ! (mostly, but not entirely, related to GrADS-format output)
      !------------------------------------------------------------

      character*15, intent(inout) :: tdef
      character*3, intent(in), dimension(maxq) :: qname
      character*6, intent(in), dimension(maxq) :: budname
      real, intent(in), dimension(ib:ie) :: xh
      real, intent(in), dimension(ib:ie+1) :: xf
      real, intent(in), dimension(jb:je) :: yh
      real, intent(in), dimension(jb:je+1) :: yf
      real, intent(in), dimension(-2:nx+4) :: xfref
      real, intent(in), dimension(-2:ny+4) :: yfref
      real, intent(in), dimension(kb:ke) :: sigma
      real, intent(in), dimension(kb:ke+1) :: sigmaf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: zh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: zf

!-----------------------------------------------------------------------

      integer :: i,j,k,n,flag
      character*8 text1
      character*30 text2
      character*50 fname
      character*8,  dimension(:), allocatable :: varname
      character*30, dimension(:), allocatable :: vardesc

!-----------------------------------------------------------------------
! get length of output_path string

    flag=0
    n=0
    do while( flag.eq.0 .and. n.le.70 )
      n=n+1
      if( output_path(n:n).eq.' ' .or. output_path(n:n).eq.'.' ) flag=1
    enddo

    strlen=n-1

!--------------------------------------
! get length of output_basename string

    flag=0
    n=0
    do while( flag.eq.0 .and. n.le.70 )
      n=n+1
      if( output_basename(n:n).eq.' ' .or. output_basename(n:n).eq.'.' ) flag=1
    enddo

    baselen=n-1

!------

    totlen = strlen + baselen

    IF( totlen .gt. (70-22) )THEN
      IF(myid.eq.0)THEN
      print *
      print *,'  baselen = ',baselen
      print *,'  strlen  = ',strlen
      print *,'  totlen  = ',totlen
      print *
      print *,'  totlen is too long ... make either baselen or strlen shorter '
      print *
      print *,'  stopping cm1 .... '
      print *
      ENDIF
#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
#endif
      call stopcm1
    ENDIF

!------

      string = '                                                                      '
    statfile = '                                                                      '
     sstring = '                                                                      '

  if(strlen.gt.0)then
      string(1:strlen) = output_path(1:strlen)
    statfile(1:strlen) = output_path(1:strlen)
  endif

      string(strlen+1:strlen+baselen) = output_basename(1:baselen)
    statfile(strlen+1:strlen+baselen) = output_basename(1:baselen)
     sstring(1:baselen) = output_basename(1:baselen)

    statfile(totlen+1:totlen+22) = '_stats.dat            '

  IF(output_format.eq.1)THEN
    if(dowr) write(outfile,*)
    if(dowr) write(outfile,*) '  writing ctl files ... '
  ENDIF
    if(dowr) write(outfile,*)
    if(dowr) write(outfile,*) '  strlen          = ',strlen
    if(dowr) write(outfile,*) '  baselen         = ',baselen
    if(dowr) write(outfile,*) '  totlen          = ',totlen
  if(strlen.gt.0)then
    if(dowr) write(outfile,*) '  output_path     = ',output_path(1:strlen)
  endif
    if(dowr) write(outfile,*) '  output_basename = ',output_basename(1:baselen)
    if(dowr) write(outfile,*) '  statfile        = ',statfile
    if(dowr) write(outfile,*)

      IF( myid.eq.0 )THEN
        if(output_filetype.ge.2)then
          tdef = '00:00Z03JUL0001'
        else
          tdef = '00:00Z03JUL2000'
        endif
        IF( radopt.ge.1 )THEN
          write(tdef( 1: 2),237) hour
          write(tdef( 4: 5),237) minute
          write(tdef( 7: 8),237) day
        if(output_filetype.ge.2)then
          write(tdef(12:15),238) 1
        else
          write(tdef(12:15),238) year
        endif
237       format(i2.2)
238       format(i4.4)
          IF( month.eq.1 )THEN
            write(tdef(9:11),239) 'JAN'
          ELSEIF( month.eq.2 )THEN
            write(tdef(9:11),239) 'FEB'
          ELSEIF( month.eq.3 )THEN
            write(tdef(9:11),239) 'MAR'
          ELSEIF( month.eq.4 )THEN
            write(tdef(9:11),239) 'APR'
          ELSEIF( month.eq.5 )THEN
            write(tdef(9:11),239) 'MAY'
          ELSEIF( month.eq.6 )THEN
            write(tdef(9:11),239) 'JUN'
          ELSEIF( month.eq.7 )THEN
            write(tdef(9:11),239) 'JUL'
          ELSEIF( month.eq.8 )THEN
            write(tdef(9:11),239) 'AUG'
          ELSEIF( month.eq.9 )THEN
            write(tdef(9:11),239) 'SEP'
          ELSEIF( month.eq.10 )THEN
            write(tdef(9:11),239) 'OCT'
          ELSEIF( month.eq.11 )THEN
            write(tdef(9:11),239) 'NOV'
          ELSEIF( month.eq.12 )THEN
            write(tdef(9:11),239) 'DEC'
          ELSE
            print *
            print *,'  Invalid value for MONTH '
            print *
            print *,'  Stopping CM1 .... '
            print *
            call stopcm1
          ENDIF
239       format(a3)
        ENDIF
      ENDIF

!-----------------------------------------------------------------------
!  GrADS descriptor files
!-----------------------------------------------------------------------

  grads_descriptors: IF( output_format.eq.1 )THEN

      IF(myid.eq.0)THEN

        allocate( varname(1000) )
        allocate( vardesc(1000) )

!----------------------------
! s file:
! accounts for both 2d and 3d variables:

    sout2d = 0
    s_out = 0

    ! all 2d variables MUST be listed first:

    if(output_rain   .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'rn      '
      vardesc(s_out) = 'accumulated rainfall (cm)     '
    endif
    if(output_sws    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'sws     '
      vardesc(s_out) = 'max wind speed lwst lvl (m/s) '
    endif
    if(output_svs    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'svs     '
      vardesc(s_out) = 'max vert vort lwst lvl (s-1)  '
    endif
    if(output_sps    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'sps     '
      vardesc(s_out) = 'min pressure lowest level (Pa)'
    endif
    if(output_srs    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'srs     '
      vardesc(s_out) = 'max sfc rainwater (kg/kg)     '
    endif
    if(output_sgs    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'sgs     '
      vardesc(s_out) = 'max sfc graupel/hail (kg/kg)  '
    endif
    if(output_sus    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'sus     '
      vardesc(s_out) = 'max w at 5 km AGL (m/s)       '
    endif
    if(output_shs    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'shs     '
      vardesc(s_out) = 'max integrated uh (m2/s2)     '
    endif
    if(nrain.eq.2)then
      if(output_rain   .eq.1)then
        s_out = s_out + 1
        varname(s_out) = 'rn2     '
        vardesc(s_out) = 'translated rainfall (cm)      '
      endif
      if(output_sws    .eq.1)then
        s_out = s_out + 1
        varname(s_out) = 'sws2    '
        vardesc(s_out) = 'translated max wind (m/s)     '
      endif
      if(output_svs    .eq.1)then
        s_out = s_out + 1
        varname(s_out) = 'svs2    '
        vardesc(s_out) = 'translated max vorticity (s-1)'
      endif
      if(output_sps    .eq.1)then
        s_out = s_out + 1
        varname(s_out) = 'sps2    '
        vardesc(s_out) = 'translated min pressure (Pa)  '
      endif
      if(output_srs    .eq.1)then
        s_out = s_out + 1
        varname(s_out) = 'srs2    '
        vardesc(s_out) = 'translated max rainwater      '
      endif
      if(output_sgs    .eq.1)then
        s_out = s_out + 1
        varname(s_out) = 'sgs2    '
        vardesc(s_out) = 'translated max graupel/hail   '
      endif
      if(output_sus    .eq.1)then
        s_out = s_out + 1
        varname(s_out) = 'sus2    '
        vardesc(s_out) = 'translated max w at 5 km (m/s)'
      endif
      if(output_shs    .eq.1)then
        s_out = s_out + 1
        varname(s_out) = 'shs2    '
        vardesc(s_out) = 'translated max integrated uh  '
      endif
    endif
    if(output_uh.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'uh      '
      vardesc(s_out) = 'integ. updraft helicity (m2/s2'
    endif
    if(output_coldpool.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'cpc     '
      vardesc(s_out) = 'cold pool intensity C (m/s)   '
      s_out = s_out + 1
      varname(s_out) = 'cph     '
      vardesc(s_out) = 'cold pool depth h (m AGL)     '
    endif
    if(output_sfcflx .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'thflux  '
      vardesc(s_out) = 'sfc theta flux (K m/s)        '
      s_out = s_out + 1
      varname(s_out) = 'qvflux  '
      vardesc(s_out) = 'sfc water vapor flux (g/g m/s)'
      s_out = s_out + 1
      varname(s_out) = 'cd      '
      vardesc(s_out) = 'exchange coeff for momentum   '
      s_out = s_out + 1
      varname(s_out) = 'ch      '
      vardesc(s_out) = 'exchange coeff for sens. heat '
      s_out = s_out + 1
      varname(s_out) = 'cq      '
      vardesc(s_out) = 'exchange coeff for moisture   '
      s_out = s_out + 1
      varname(s_out) = 'tsk     '
      vardesc(s_out) = 'soil/ocean temperature (K)    '
    endif
    if(output_zs     .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'zs      '
      vardesc(s_out) = 'terrain height (m)            '
    endif
    if(output_dbz    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'cref    '
      vardesc(s_out) = 'composite reflectivity (dBZ)  '
    endif
    if(output_sfcparams.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'xland   '
      vardesc(s_out) = 'land/water flag (1=land,2=wtr)'
      s_out = s_out + 1
      varname(s_out) = 'lu      '
      vardesc(s_out) = 'land use index                '
      s_out = s_out + 1
      varname(s_out) = 'mavail  '
      vardesc(s_out) = 'surface moisture availability '
    endif
    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.oceanmodel.eq.2))then
      s_out = s_out + 1
      varname(s_out) = 'tmn     '
      vardesc(s_out) = 'deep-layer soil temperature (K'
      s_out = s_out + 1
      varname(s_out) = 'hfx     '
      vardesc(s_out) = 'heat flux at surface (W/m^2)  '
      s_out = s_out + 1
      varname(s_out) = 'qfx     '
      vardesc(s_out) = 'moisture flux at sfc (kg/m^2/s'
      s_out = s_out + 1
      varname(s_out) = 'gsw     '
      vardesc(s_out) = 'downward SW flux at sfc (W/m2)'
      s_out = s_out + 1
      varname(s_out) = 'glw     '
      vardesc(s_out) = 'downward LW flux at sfc (W/m2)'
    endif
    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2))then
      s_out = s_out + 1
      varname(s_out) = 'tslb1   '
      vardesc(s_out) = 'soil temp, layer 1 (K)        '
      s_out = s_out + 1
      varname(s_out) = 'tslb2   '
      vardesc(s_out) = 'soil temp, layer 2 (K)        '
      s_out = s_out + 1
      varname(s_out) = 'tslb3   '
      vardesc(s_out) = 'soil temp, layer 3 (K)        '
      s_out = s_out + 1
      varname(s_out) = 'tslb4   '
      vardesc(s_out) = 'soil temp, layer 4 (K)        '
      s_out = s_out + 1
      varname(s_out) = 'tslb5   '
      vardesc(s_out) = 'soil temp, layer 5 (K)        '
    endif
    if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then
      s_out = s_out + 1
      varname(s_out) = 'tml     '
      vardesc(s_out) = 'ocean mixed layer temp (K)    '
      s_out = s_out + 1
      varname(s_out) = 'hml     '
      vardesc(s_out) = 'ocean mixed layer depth (m)   '
      s_out = s_out + 1
      varname(s_out) = 'huml    '
      vardesc(s_out) = 'ocean mixed layer u vel. (m/s)'
      s_out = s_out + 1
      varname(s_out) = 'hvml    '
      vardesc(s_out) = 'ocean mixed layer v vel. (m/s)'
    endif
    if(output_radten.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'radsw   '
      vardesc(s_out) = 'solar radiation at surface    '
      s_out = s_out + 1
      varname(s_out) = 'rnflx   '
      vardesc(s_out) = 'net radiation absorbed by sfc '
      s_out = s_out + 1
      varname(s_out) = 'radswnet'
      vardesc(s_out) = 'net solar radiation           '
      s_out = s_out + 1
      varname(s_out) = 'radlwin '
      vardesc(s_out) = 'incoming longwave radiation   '
! MS addition - toa fluxes
      s_out = s_out + 1
      varname(s_out) = 'olr     '
      vardesc(s_out) = 'TOA net outgoing longwave rad.'
      s_out = s_out + 1
      varname(s_out) = 'dsr     '
      vardesc(s_out) = 'TOA net incoming solar rad.   '
    endif
    IF(output_sfcdiags.eq.1)THEN
      s_out = s_out + 1
      varname(s_out) = 'u10     '
      vardesc(s_out) = 'diagnostic 10m u wind (m/s)   '
      s_out = s_out + 1
      varname(s_out) = 'v10     '
      vardesc(s_out) = 'diagnostic 10m v wind (m/s)   '
      s_out = s_out + 1
      varname(s_out) = 't2      '
      vardesc(s_out) = 'diagnostic 2m temperature (K) '
      s_out = s_out + 1
      varname(s_out) = 'q2      '
      vardesc(s_out) = 'diagnostic 2m mixing ratio g/g'
      s_out = s_out + 1
      varname(s_out) = 'znt     '
      vardesc(s_out) = 'roughness length (m)          '
      s_out = s_out + 1
      varname(s_out) = 'ust     '
      vardesc(s_out) = 'u* in similarity theory (m/s) '
      s_out = s_out + 1
      varname(s_out) = 'hpbl    '
    if(ipbl.eq.1)then
      vardesc(s_out) = 'PBL height (m) (from PBL schem'
    else
      vardesc(s_out) = 'rough estimate of PBL hght (m)'
    endif
      s_out = s_out + 1
      varname(s_out) = 'zol     '
      vardesc(s_out) = 'z/L (z over Monin-Obukhov len)'
      s_out = s_out + 1
      varname(s_out) = 'mol     '
      vardesc(s_out) = 'T* (similarity theory) (K)    '
      s_out = s_out + 1
      varname(s_out) = 'br      '
      vardesc(s_out) = 'bulk Richardson No in sfc lay.'
    ENDIF

    ! done with 2d variables

    sout2d = s_out

    ! Now, all 3d variables:

    if(output_zh     .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'zh      '
      vardesc(s_out) = 'height on model levels (m)    '
    endif
    if(output_th     .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'th      '
      vardesc(s_out) = 'potential temp. (K)           '
    endif
    if(output_thpert .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'thpert  '
      vardesc(s_out) = 'potential temp. pert. (K)     '
    endif
    if(output_prs    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'prs     '
      vardesc(s_out) = 'pressure (Pa)                 '
    endif
    if(output_prspert.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'prspert '
      vardesc(s_out) = 'pressure pert. (Pa)           '
    endif
    if(output_pi     .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'pi      '
      vardesc(s_out) = 'nondimensional pressure       '
    endif
    if(output_pipert .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'pipert  '
      vardesc(s_out) = 'nondimensional pressure pert. '
    endif
    if(output_rho    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'rho     '
      vardesc(s_out) = 'density (kg/m^3)              '
    endif
    if(output_rhopert.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'rhopert '
      vardesc(s_out) = 'density pert. (kg/m^3)        '
    endif
    if(iptra         .eq.1)then
      do n=1,npt
        text1='pt      '
        if(n.lt.10)then
          write(text1(3:3),155) n
155       format(i1.1)
        else
          write(text1(3:4),154) n
154       format(i2.2)
        endif
        s_out = s_out + 1
        varname(s_out) = text1
        vardesc(s_out) = 'passive tracer                '
      enddo
    endif
    if(output_qv     .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'qv      '
      vardesc(s_out) = 'water vapor mixing ratio      '
    endif
    if(output_qvpert .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'qvpert  '
      vardesc(s_out) = 'qv pert                       '
    endif
    if(output_q      .eq.1)then
      do n=1,numq
        if(n.ne.nqv)then
          text1='        '
          text2='                              '
          write(text1(1:3),156) qname(n)
          write(text2(1:3),156) qname(n)
156       format(a3)
          s_out = s_out + 1
          varname(s_out) = text1
          vardesc(s_out) = text2
        endif
      enddo
    endif
    if(output_dbz    .eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'dbz     '
      vardesc(s_out) = 'reflectivity (dBZ)            '
    endif
    if(output_buoyancy.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'buoyancy'
      vardesc(s_out) = 'buoyancy (m s^-2)             '
    endif
    if(output_uinterp.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'uinterp '
      vardesc(s_out) = 'u interp. to scalar points    '
    endif
    if(output_vinterp.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'vinterp '
      vardesc(s_out) = 'v interp. to scalar points    '
    endif
    if(output_winterp.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'winterp '
      vardesc(s_out) = 'w interp. to scalar points    '
    endif
    if(output_vort.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'xvort   '
      vardesc(s_out) = 'horiz vorticity (x) (s^-1)    '
      s_out = s_out + 1
      varname(s_out) = 'yvort   '
      vardesc(s_out) = 'horiz vorticity (y) (s^-1)    '
      s_out = s_out + 1
      varname(s_out) = 'zvort   '
      vardesc(s_out) = 'vertical vorticity (s^-1)     '
    endif
    if(output_basestate.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'pi0     '
      vardesc(s_out) = 'base-state nondim. pressure   '
      s_out = s_out + 1
      varname(s_out) = 'th0     '
      vardesc(s_out) = 'base-state potential temp (K) '
      s_out = s_out + 1
      varname(s_out) = 'prs0    '
      vardesc(s_out) = 'base-state pressure (Pa)      '
      s_out = s_out + 1
      varname(s_out) = 'qv0     '
      vardesc(s_out) = 'base-state qv (kg/kg)         '
    endif
    if(output_pblten.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'thpten  '
      vardesc(s_out) = 'pbl tendency:  theta          '
      s_out = s_out + 1
      varname(s_out) = 'qvpten  '
      vardesc(s_out) = 'pbl tendency:  qv             '
      s_out = s_out + 1
      varname(s_out) = 'qcpten  '
      vardesc(s_out) = 'pbl tendency:  qc             '
      s_out = s_out + 1
      varname(s_out) = 'qipten  '
      vardesc(s_out) = 'pbl tendency:  qi             '
      s_out = s_out + 1
      varname(s_out) = 'upten   '
      vardesc(s_out) = 'pbl tendency:  u              '
      s_out = s_out + 1
      varname(s_out) = 'vpten   '
      vardesc(s_out) = 'pbl tendency:  v              '
    endif
    if(output_radten.eq.1)then
      s_out = s_out + 1
      varname(s_out) = 'swten   '
      vardesc(s_out) = 'pot temp tendency, sw rad (K/s'
      s_out = s_out + 1
      varname(s_out) = 'lwten   '
      vardesc(s_out) = 'pot temp tendency, lw rad (K/s'
    endif

    sout3d = s_out - sout2d

!----------------------------
!  ready to write GrADS descriptor file:

  IF(s_out.ge.1)THEN
    string(totlen+1:totlen+22) = '_s.ctl                '
    if(dowr) write(outfile,*) string
    open(unit=50,file=string,status='unknown')

  if(output_filetype.eq.1)then
    sstring(baselen+1:baselen+1+12) = '_s.dat'
  elseif(output_filetype.ge.2)then
    sstring(baselen+1:baselen+1+12) = '_00%y4_s.dat'
  endif
    write(50,201) sstring
!!!    write(50,222)
    if(output_filetype.ge.2) write(50,221)
    write(50,202)
    write(50,203)
    if(stretch_x.ge.1)then
      write(50,214) nx
      do i=1,nx
        write(50,217) 0.001*0.5*(xfref(i)+xfref(i+1))
      enddo
    else
      write(50,204) nx,xh(1)/1000.0,dx/1000.0
    endif
    if(stretch_y.ge.1)then
      write(50,215) ny
      do j=1,ny
        write(50,217) 0.001*0.5*(yfref(j)+yfref(j+1))
      enddo
    else
      write(50,205) ny,yh(1)/1000.0,dy/1000.0
    endif
    if(stretch_z.eq.0)then
      write(50,206) nz,0.5*dz/1000.0,dz/1000.0
    else
      write(50,216) nz
      do k=1,nz
        write(50,217) 0.001*sigma(k)
      enddo
    endif
  if(output_filetype.eq.1)then
    write(50,207) int(1+timax/tapfrq),tdef,max(1,int(tapfrq/60.0))
  elseif(output_filetype.ge.2)then
    write(50,227) int(1+timax/tapfrq),tdef
  endif
    write(50,208) s_out
    ! account for both 2d and 3d output files:
    do n=1,sout2d
      write(50,209) varname(n), 0,vardesc(n)
    enddo
    do n=sout2d+1,s_out
      write(50,209) varname(n),nk,vardesc(n)
    enddo
    write(50,210)
    close(unit=50)
  ENDIF

!-----------------------------------
! i file:  (for interpolated output when using terrain)
!   follows s file very closely:
!   no need to re-define varname,vardesc...

  IF(s_out.ge.1 .and. terrain_flag .and. output_interp.eq.1)THEN
    string(totlen+1:totlen+22) = '_i.ctl                '
    if(dowr) write(outfile,*) string
    open(unit=50,file=string,status='unknown')

  if(output_filetype.eq.1)then
    sstring(baselen+1:baselen+1+12) = '_i.dat'
  elseif(output_filetype.ge.2)then
    sstring(baselen+1:baselen+1+12) = '_00%y4_i.dat'
  endif

    write(50,201) sstring
!!!    write(50,222)
    if(output_filetype.ge.2) write(50,221)
    write(50,202)
    write(50,203)
    if(stretch_x.ge.1)then
      write(50,214) nx
      do i=1,nx
        write(50,217) 0.001*0.5*(xfref(i)+xfref(i+1))
      enddo
    else
      write(50,204) nx,xh(1)/1000.0,dx/1000.0
    endif
    if(stretch_y.ge.1)then
      write(50,215) ny
      do j=1,ny
        write(50,217) 0.001*0.5*(yfref(j)+yfref(j+1))
      enddo
    else
      write(50,205) ny,yh(1)/1000.0,dy/1000.0
    endif
    if(stretch_z.eq.0)then
      write(50,206) nz,0.5*dz/1000.0,dz/1000.0
    else
      write(50,216) nz
      do k=1,nz
        write(50,217) 0.001*sigma(k)
      enddo
    endif
  if(output_filetype.eq.1)then
    write(50,207) int(1+timax/tapfrq),tdef,max(1,int(tapfrq/60.0))
  elseif(output_filetype.ge.2)then
    write(50,227) int(1+timax/tapfrq),tdef
  endif
    write(50,208) s_out
    ! account for both 2d and 3d output files:
    do n=1,sout2d
      write(50,209) varname(n), 0,vardesc(n)
    enddo
    do n=sout2d+1,s_out
      write(50,209) varname(n),nk,vardesc(n)
    enddo
    write(50,210)
    close(unit=50)
  ENDIF

!-----------------------------------
! u file:
! I have assumed that all variables are 3d for this file.

    u_out = 0

    if(output_u    .eq.1)then
      u_out = u_out + 1
      varname(u_out) = 'u       '
      vardesc(u_out) = 'E-W velocity (m/s)            '
    endif
    if(output_upert.eq.1)then
      u_out = u_out + 1
      varname(u_out) = 'upert   '
      vardesc(u_out) = 'u pert. (m/s)                 '
    endif
    if(output_basestate.eq.1)then
      u_out = u_out + 1
      varname(u_out) = 'u0      '
      vardesc(u_out) = 'base-state u (m/s)            '
    endif
    if(output_turbten.eq.1)then
      u_out = u_out + 1
      varname(u_out) = 'ftu     '
      vardesc(u_out) = 'u tendency: turbulence scheme '
    endif
    if(output_impdiften.eq.1)then
      u_out = u_out + 1
      varname(u_out) = 'fdu     '
      vardesc(u_out) = 'u tendency: implicit diffusion'
    endif

  IF(u_out.ge.1)THEN
    string(totlen+1:totlen+22) = '_u.ctl                '
    if(dowr) write(outfile,*) string
    open(unit=50,file=string,status='unknown')

  if(output_filetype.eq.1)then
    sstring(baselen+1:baselen+1+12) = '_u.dat'
  elseif(output_filetype.ge.2)then
    sstring(baselen+1:baselen+1+12) = '_00%y4_u.dat'
  endif

    write(50,201) sstring
    if(output_filetype.ge.2) write(50,221)
    write(50,202)
    write(50,203)
    if(stretch_x.ge.1)then
      write(50,214) nx+1
      do i=1,nx+1
        write(50,217) 0.001*xfref(i)
      enddo
    else
      write(50,204) nx+1,xf(1)/1000.0,dx/1000.0
    endif
    if(stretch_y.ge.1)then
      write(50,215) ny
      do j=1,ny
        write(50,217) 0.001*0.5*(yfref(j)+yfref(j+1))
      enddo
    else
      write(50,205) ny,yh(1)/1000.0,dy/1000.0
    endif
    if(stretch_z.eq.0)then
      write(50,206) nz,0.5*dz/1000.0,dz/1000.0
    else
      write(50,216) nz
      do k=1,nz
        write(50,217) 0.001*sigma(k)
      enddo
    endif
  if(output_filetype.eq.1)then
    write(50,207) int(1+timax/tapfrq),tdef,max(1,int(tapfrq/60.0))
  elseif(output_filetype.ge.2)then
    write(50,227) int(1+timax/tapfrq),tdef
  endif
    write(50,208) u_out
    ! assumes all variables are 3d:
    do n=1,u_out
      write(50,209) varname(n),nk,vardesc(n)
    enddo
    write(50,210)
    close(unit=50)
  ENDIF

!-----------------------------------
! v file:
! I have assumed that all variables are 3d for this file.

    v_out = 0

    if(output_v    .eq.1)then
      v_out = v_out + 1
      varname(v_out) = 'v       '
      vardesc(v_out) = 'N-S velocity (m/s)            '
    endif
    if(output_vpert.eq.1)then
      v_out = v_out + 1
      varname(v_out) = 'vpert   '
      vardesc(v_out) = 'v pert (m/s)                  '
    endif
    if(output_basestate.eq.1)then
      v_out = v_out + 1
      varname(v_out) = 'v0      '
      vardesc(v_out) = 'base-state v (m/s)            '
    endif
    if(output_turbten.eq.1)then
      v_out = v_out + 1
      varname(v_out) = 'ftv     '
      vardesc(v_out) = 'v tendency: turbulence scheme '
    endif
    if(output_impdiften.eq.1)then
      v_out = v_out + 1
      varname(v_out) = 'fdv     '
      vardesc(v_out) = 'v tendency: implicit diffusion'
    endif

  IF(v_out.ge.1)THEN
    string(totlen+1:totlen+22) = '_v.ctl                '
    if(dowr) write(outfile,*) string
    open(unit=50,file=string,status='unknown')

  if(output_filetype.eq.1)then
    sstring(baselen+1:baselen+1+12) = '_v.dat'
  elseif(output_filetype.ge.2)then
    sstring(baselen+1:baselen+1+12) = '_00%y4_v.dat'
  endif

    write(50,201) sstring
    if(output_filetype.ge.2) write(50,221)
    write(50,202)
    write(50,203)
    if(stretch_x.ge.1)then
      write(50,214) nx
      do i=1,nx
        write(50,217) 0.001*0.5*(xfref(i)+xfref(i+1))
      enddo
    else
      write(50,204) nx,xh(1)/1000.0,dx/1000.0
    endif
    if(stretch_y.ge.1)then
      write(50,215) ny+1
      do j=1,ny+1
        write(50,217) 0.001*yfref(j)
      enddo
    else
      write(50,205) ny+1,yf(1)/1000.0,dy/1000.0
    endif
    if(stretch_z.eq.0)then
      write(50,206) nz,0.5*dz/1000.0,dz/1000.0
    else
      write(50,216) nz
      do k=1,nz
        write(50,217) 0.001*sigma(k)
      enddo
    endif
  if(output_filetype.eq.1)then
    write(50,207) int(1+timax/tapfrq),tdef,max(1,int(tapfrq/60.0))
  elseif(output_filetype.ge.2)then
    write(50,227) int(1+timax/tapfrq),tdef
  endif
    write(50,208) v_out
    ! assumes all variables are 3d:
    do n=1,v_out
      write(50,209) varname(n),nk,vardesc(n)
    enddo
    write(50,210)
    close(unit=50)
  ENDIF

!-----------------------------------
! w file:
! I have assumed that all variables are 3d for this file.

    w_out = 0

    if(output_w  .eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'w       '
      vardesc(w_out) = 'vertical velocity (m/s)       '
    endif
    if(output_tke.eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'tke     '
      vardesc(w_out) = 'turb. kinetic energy (m^2/s^2)'
    endif
    if(output_km .eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'kmh     '
      IF( ipbl.eq.1 )THEN
        vardesc(w_out) = 'turb. coef. for mo. (2D Smag.)'
      ELSE
        vardesc(w_out) = 'turb. coef. for mo. (m^2/s)   '
      ENDIF
      w_out = w_out + 1
      varname(w_out) = 'kmv     '
      IF( ipbl.eq.1 )THEN
        vardesc(w_out) = 'turb. coef. for mo. (from YSU)'
      ELSE
        vardesc(w_out) = 'turb. coef. for mo. (m^2/s)   '
      ENDIF
    endif
    if(output_kh .eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'khh     '
      IF( ipbl.eq.1 )THEN
        vardesc(w_out) = 'turb. coef. for scalar (2D Sm)'
      ELSE
        vardesc(w_out) = 'turb. coef. for scalar (m^2/s)'
      ENDIF
      w_out = w_out + 1
      varname(w_out) = 'khv     '
      IF( ipbl.eq.1 )THEN
        vardesc(w_out) = 'turb. coef. for scalar (YSU)  '
      ELSE
        vardesc(w_out) = 'turb. coef. for scalar (m^2/s)'
      ENDIF
    endif
    if(output_dissten.eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'dissten '
      vardesc(w_out) = 'dissipation rate (m^2/s^3)    '
    endif
    if(output_nm.eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'nm      '
      vardesc(w_out) = 'squared Brunt-Vaisala freq    '
    endif
    if(output_def.eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'defv    '
      vardesc(w_out) = 'vertical deformation          '
      w_out = w_out + 1
      varname(w_out) = 'defh    '
      vardesc(w_out) = 'horizontal deformation        '
    endif
    if(output_turbten.eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'ftw     '
      vardesc(w_out) = 'w tendency: turbulence scheme '
    endif
    if(output_impdiften.eq.1)then
      w_out = w_out + 1
      varname(w_out) = 'fdw     '
      vardesc(w_out) = 'w tendency: implicit diffusion'
    endif

  IF(w_out.ge.1)THEN
    string(totlen+1:totlen+22) = '_w.ctl                '
    if(dowr) write(outfile,*) string
    open(unit=50,file=string,status='unknown')

  if(output_filetype.eq.1)then
    sstring(baselen+1:baselen+1+12) = '_w.dat'
  elseif(output_filetype.ge.2)then
    sstring(baselen+1:baselen+1+12) = '_00%y4_w.dat'
  endif

    write(50,201) sstring
    if(output_filetype.ge.2) write(50,221)
    write(50,202)
    write(50,203)
    if(stretch_x.ge.1)then
      write(50,214) nx
      do i=1,nx
        write(50,217) 0.001*0.5*(xfref(i)+xfref(i+1))
      enddo
    else
      write(50,204) nx,xh(1)/1000.0,dx/1000.0
    endif
    if(stretch_y.ge.1)then
      write(50,215) ny
      do j=1,ny
        write(50,217) 0.001*0.5*(yfref(j)+yfref(j+1))
      enddo
    else
      write(50,205) ny,yh(1)/1000.0,dy/1000.0
    endif
    if(stretch_z.eq.0)then
      write(50,206) nz+1,0.0,dz/1000.0
    else
      write(50,216) nz+1
      do k=1,nz+1
        write(50,217) 0.001*sigmaf(k)
      enddo
    endif
  if(output_filetype.eq.1)then
    write(50,207) int(1+timax/tapfrq),tdef,max(1,int(tapfrq/60.0))
  elseif(output_filetype.ge.2)then
    write(50,227) int(1+timax/tapfrq),tdef
  endif
    write(50,208) w_out
    ! assumes all variables are 3d:
    do n=1,w_out
      write(50,209) varname(n),nk+1,vardesc(n)
    enddo
    write(50,210)
    close(unit=50)
  ENDIF

!-----------------------------------

    if(dowr) write(outfile,*)

201   format('dset ^',a70)
202   format('title CM1 output')
221   format('options template')
222   format('byteswapped')
203   format('undef -99999999.')
204   format('xdef ',i6,' linear ',f13.6,1x,f13.6)
214   format('xdef ',i6,' levels ')
205   format('ydef ',i6,' linear ',f13.6,1x,f13.6)
215   format('ydef ',i6,' levels ')
206   format('zdef ',i6,' linear ',f13.6,1x,f13.6)
216   format('zdef ',i6,' levels ')
217   format(2x,f13.6)
207   format('tdef ',i10,' linear ',a15,' ',i5,'MN')
227   format('tdef ',i10,' linear ',a15,' 1YR')
208   format('vars ',i4)
209   format(a8,2x,i6,'  99  ',a30)
210   format('endvars')

211   format(2x,f7.3)

!-----------------------------------------------------------------------

      call write_statsctl(tdef,qname,budname,1+nint(timax/max(statfrq,dtl)))

!-----------------------------------------------------------------------
!  Parcel data file:

      if(iprcl.eq.1.and.myid.eq.0)then

        string(totlen+1:totlen+22) = '_pdata.ctl            '
        if(dowr) write(outfile,*) string
        open(unit=50,file=string,status='unknown')

        sstring(baselen+1:baselen+1+12) = '_pdata.dat  '

        write(50,401) sstring
        write(50,402)
        write(50,403)
        write(50,404) nparcels
        write(50,405)
        write(50,406)
        write(50,407) 1+int(timax/prclfrq),tdef,max(1,int(prclfrq/60.0))
        write(50,408) npvals
                      write(50,409) 'x       ','x (m)                         '
                      write(50,409) 'y       ','y (m)                         '
                      write(50,409) 'z       ','z (m)                         '
                      write(50,409) 'u       ','u (m/s)                       '
                      write(50,409) 'v       ','v (m/s)                       '
                      write(50,409) 'w       ','w (m/s)                       '
                      write(50,409) 'th      ','potential temperature (K)     '
                      write(50,409) 'prs     ','pressure (Pa)                 '
                      write(50,409) 'b       ','buoyancy (m s^-2)             '
        if(iptra.eq.1)then
          do n=1,npt
            text1='pt      '
            if(n.lt.10)then
              write(text1(3:3),155) n
            else
              write(text1(3:4),154) n
            endif
                      write(50,409) text1     ,'passive tracer conc. (kg/kg)  '
          enddo
        endif
        if(prqv.ge.1) write(50,409) 'qv      ','mixing ratio of water vapor   '
        if(prqc.ge.1) write(50,409) 'qc      ','mixing ratio of cloud water   '
        if(prqr.ge.1) write(50,409) 'qr      ','mixing ratio of rain water    '
        if(prqi.ge.1) write(50,409) 'qi      ','mixing ratio of ice crystals  '
        if(prqs.ge.1) write(50,409) 'qs      ','mixing ratio of snow          '
        if(prqg.ge.1) write(50,409) 'qg      ','mixing ratio of graupel       '
        if(prkh.ge.1) write(50,409) 'kh      ','turb. coef. for scalar (m^2/s)'
        write(50,410)

401     format('dset ^',a70)
402     format('undef -9999.')
403     format('title ctl file for pdata.dat')
404     format('xdef ',i10,' linear 1 1')
405     format('ydef          1 linear 1 1')
406     format('zdef          1 linear 1 1')
407     format('tdef ',i10,' linear ',a15,' ',i5,'MN')
408     format('vars ',i6)
409     format(a8,' 1 99 ',a30)
410     format('endvars')

        close(unit=50)

      endif

!-----------------------------------------------------------------------

        deallocate( varname )
        deallocate( vardesc )

      ENDIF     ! endif for myid=0

#ifdef MPI
      call MPI_BCAST(sout2d,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sout3d,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(s_out ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(u_out ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(v_out ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(w_out ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(z_out ,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#endif

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  sout2d = ',sout2d
      if(dowr) write(outfile,*) '  sout3d = ',sout3d
      if(dowr) write(outfile,*) '  s_out  = ',s_out
      if(dowr) write(outfile,*) '  u_out  = ',u_out
      if(dowr) write(outfile,*) '  v_out  = ',v_out
      if(dowr) write(outfile,*) '  w_out  = ',w_out
      if(dowr) write(outfile,*) '  z_out  = ',z_out

  ENDIF grads_descriptors

      if(dowr) write(outfile,*)

!-----------------------------------------------------------------------

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write_statsctl(tdef,qname,budname,numt)
      implicit none
      include 'input.incl'

      !---------------------------------------------------------------
      ! This subroutine creates the GrADS-format stats descriptor file
      !---------------------------------------------------------------

      character*15, intent(in) :: tdef
      character*3, intent(in), dimension(maxq) :: qname
      character*6, intent(in), dimension(maxq) :: budname
      integer, intent(in) :: numt

      integer :: n
      character*8 text1
      character*30 text2
      character*50 fname

!  Subroutine to write GrADS stats descriptor file:
!-----------------------------------------------------------------------
!  write descriptors for stats file:

    string(totlen+1:totlen+22) = '_stats.ctl            '
    if(dowr) write(outfile,*) string
    open(unit=50,file=string,status='unknown')

    sstring(baselen+1:baselen+1+12) = '_stats.dat  '

    write(50,301) sstring
    write(50,302)
    write(50,303)
    write(50,304)
    write(50,305)
    write(50,306)
    write(50,307) numt,tdef,max(1,int(max(statfrq,60.0)/60.0))
    write(50,308) stat_out
    IF( adapt_dt.eq.1 )   write(50,309) 'dt      ','average timestep dt (s)       '
    if(stat_w      .eq.1) write(50,309) 'wmax    ','max vertical velocity (m/s)   '
    if(stat_w      .eq.1) write(50,309) 'wmin    ','min vertical velocity (m/s)   '
    if(stat_u      .eq.1) write(50,309) 'umax    ','max E-W velocity (m/s)        '
    if(stat_u      .eq.1) write(50,309) 'umin    ','min E-W velocity (m/s)        '
    if(stat_u      .eq.1) write(50,309) 'sumax   ','max E-W velocity lwst lvl(m/s)'
    if(stat_u      .eq.1) write(50,309) 'sumin   ','min E-W velocity lwst lvl(m/s)'
    if(stat_v      .eq.1) write(50,309) 'vmax    ','max N-S velocity (m/s)        '
    if(stat_v      .eq.1) write(50,309) 'vmin    ','min N-S velocity (m/s)        '
    if(stat_v      .eq.1) write(50,309) 'svmax   ','max N-S velocity lwst lvl(m/s)'
    if(stat_v      .eq.1) write(50,309) 'svmin   ','min N-S velocity lwst lvl(m/s)'
    if(stat_rmw    .eq.1) write(50,309) 'rmw     ','radius (m) of maximum windspd '
    if(stat_rmw    .eq.1) write(50,309) 'zmw     ','height (m) of maximum windspd '
    if(stat_pipert .eq.1) write(50,309) 'ppimax  ','max pi pert.                  '
    if(stat_pipert .eq.1) write(50,309) 'ppimin  ','min pi pert.                  '
    if(stat_prspert.eq.1) write(50,309) 'ppmax   ','max prs pert.(Pa)             '
    if(stat_prspert.eq.1) write(50,309) 'ppmin   ','min prs pert.(Pa)             '
    if(stat_thpert .eq.1) write(50,309) 'thpmax  ','max potential temp. pert. (K) '
    if(stat_thpert .eq.1) write(50,309) 'thpmin  ','min potential temp. pert. (K) '
    if(stat_thpert .eq.1) write(50,309) 'sthpmax ','max pot temp pert lwst lvl (K)'
    if(stat_thpert .eq.1) write(50,309) 'sthpmin ','min pot temp pert lwst lvl (K)'
    if(stat_q      .eq.1)then
      do n=1,numq
        text1='max     '
        text2='max                           '
        write(text1(4:6),156) qname(n)
        write(text2(5:7),156) qname(n)
        write(50,309) text1,text2
        text1='min     '
        text2='min                           '
        write(text1(4:6),156) qname(n)
        write(text2(5:7),156) qname(n)
        write(50,309) text1,text2
      enddo
    endif
    if(stat_tke    .eq.1) write(50,309) 'tkemax  ','max tke (m^2/s^2)             '
    if(stat_tke    .eq.1) write(50,309) 'tkemin  ','min tke (m^2/s^2)             '
    if(stat_km     .eq.1) write(50,309) 'kmhmax  ','max kmh (m^2/s)               '
    if(stat_km     .eq.1) write(50,309) 'kmhmin  ','min kmh (m^2/s)               '
    if(stat_km     .eq.1) write(50,309) 'kmvmax  ','max kmv (m^2/s)               '
    if(stat_km     .eq.1) write(50,309) 'kmvmin  ','min kmv (m^2/s)               '
    if(stat_kh     .eq.1) write(50,309) 'khhmax  ','max khh (m^2/s)               '
    if(stat_kh     .eq.1) write(50,309) 'khhmin  ','min khh (m^2/s)               '
    if(stat_kh     .eq.1) write(50,309) 'khvmax  ','max khv (m^2/s)               '
    if(stat_kh     .eq.1) write(50,309) 'khvmin  ','min khv (m^2/s)               '
    if(stat_div    .eq.1) write(50,309) 'divmax  ','max 3d divergence             '
    if(stat_div    .eq.1) write(50,309) 'divmin  ','min 3d divergence             '
    if(stat_rh     .eq.1) write(50,309) 'rhmax   ','max relative humidity         '
    if(stat_rh     .eq.1) write(50,309) 'rhmin   ','min relative humidity         '
    if(stat_rhi    .eq.1) write(50,309) 'rhimax  ','max relative humidity wrt ice '
    if(stat_rhi    .eq.1) write(50,309) 'rhimin  ','min relative humidity wrt ice '
    if(iptra       .eq.1)then
      do n=1,npt
        text1='maxpt   '
        text2='max pt                        '
        write(text1(6:6),157) n
        write(text2(7:7),157) n
157     format(i1)
        write(50,309) text1,text2
        text1='minpt   '
        text2='min pt                        '
        write(text1(6:6),157) n
        write(text2(7:7),157) n
        write(50,309) text1,text2
      enddo
    endif
    if(stat_the    .eq.1) write(50,309) 'themax  ','max theta-e below 10 km       '
    if(stat_the    .eq.1) write(50,309) 'themin  ','min theta-e below 10 km       '
    if(stat_the    .eq.1) write(50,309) 'sthemax ','max theta-e at lowest level   '
    if(stat_the    .eq.1) write(50,309) 'sthemin ','min theta-e at lowest level   '
    if(stat_cloud  .eq.1) write(50,309) 'qctop   ','max cloud top height (m)      '
    if(stat_cloud  .eq.1) write(50,309) 'qcbot   ','min cloud base height (m)     '
    if(stat_sfcprs .eq.1) write(50,309) 'sprsmax ','max pressure at lowest lvl (Pa'
    if(stat_sfcprs .eq.1) write(50,309) 'sprsmin ','min pressure at lowest lvl (Pa'
    if(stat_wsp    .eq.1) write(50,309) 'wspmax  ','max wind speed (m/s)          '
    if(stat_wsp    .eq.1) write(50,309) 'wspmin  ','min wind speed (m/s)          '
    if(stat_wsp    .eq.1) write(50,309) 'swspmax ','max wind speed lowst lvl (m/s)'
    if(stat_wsp    .eq.1) write(50,309) 'swspmin ','min wind speed lowst lvl (m/s)'
  IF(bbc.eq.3)THEN
    if(stat_wsp    .eq.1) write(50,309) 'wsp10max','max 10 m wind speed (m/s)     '
    if(stat_wsp    .eq.1) write(50,309) 'wsp10min','min 10 m wind speed (m/s)     '
  ENDIF
  IF( adapt_dt.eq.1 )THEN
    if(stat_cfl    .eq.1) write(50,309) 'cflmax  ','max Courant number (average)  '
  ELSE
    if(stat_cfl    .eq.1) write(50,309) 'cflmax  ','max Courant number            '
  ENDIF
    if(stat_cfl    .eq.1.and.iturb.ge.1) write(50,309) 'kshmax  ','max horiz K stability factor  '
    if(stat_cfl    .eq.1.and.iturb.ge.1) write(50,309) 'ksvmax  ','max vert K stability factor   '
    if(stat_vort   .eq.1) write(50,309) 'vortsfc ','max vert. vort. lwst lvl (1/s)'
    if(stat_vort   .eq.1) write(50,309) 'vort1km ','max vert. vort. at 1 km (1/s) '
    if(stat_vort   .eq.1) write(50,309) 'vort2km ','max vert. vort. at 2 km (1/s) '
    if(stat_vort   .eq.1) write(50,309) 'vort3km ','max vert. vort. at 3 km (1/s) '
    if(stat_vort   .eq.1) write(50,309) 'vort4km ','max vert. vort. at 4 km (1/s) '
    if(stat_vort   .eq.1) write(50,309) 'vort5km ','max vert. vort. at 5 km (1/s) '
    if(stat_tmass  .eq.1) write(50,309) 'tmass   ','total mass of (dry) air       '
    if(stat_tmois  .eq.1) write(50,309) 'tmois   ','total moisture                '
    if(stat_qmass  .eq.1)then
      do n=1,numq
        IF( (n.eq.nqv) .or.                                 &
            (n.ge.nql1.and.n.le.nql2) .or.                  &
            (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )THEN
          text1='mass    '
          text2='total mass of                 '
          write(text1( 5: 7),156) qname(n)
          write(text2(15:17),156) qname(n)
          write(50,309) text1,text2
        ENDIF
      enddo
    endif
    if(stat_tenerg .eq.1) write(50,309) 'ek      ','total kinetic energy          '
    if(stat_tenerg .eq.1) write(50,309) 'ei      ','total internal energy         '
    if(stat_tenerg .eq.1) write(50,309) 'ep      ','total potential energy        '
    if(stat_tenerg .eq.1) write(50,309) 'le      ','total latent energy (sort of) '
    if(stat_tenerg .eq.1) write(50,309) 'et      ','total energy                  '
    if(stat_mo     .eq.1) write(50,309) 'tmu     ','total E-W momentum            '
    if(stat_mo     .eq.1) write(50,309) 'tmv     ','total N-S momentum            '
    if(stat_mo     .eq.1) write(50,309) 'tmw     ','total vertical momentum       '
    if(stat_tmf    .eq.1) write(50,309) 'tmfu    ','total upward mass flux        '
    if(stat_tmf    .eq.1) write(50,309) 'tmfd    ','total downward mass flux      '
    if(stat_pcn    .eq.1)then
      do n=1,nbudget
        text1='        '
        text2='                              '
        write(text1(1:6),158) budname(n)
        write(text2(1:6),158) budname(n)
158     format(a6)
        write(50,309) text1,text2
      enddo
    endif
    if(stat_qsrc   .eq.1)then
      do n=1,numq
        text1='as      '
        text2='artificial source of          '
        write(text1( 3: 5),156) qname(n)
        write(text2(22:24),156) qname(n)
        write(50,309) text1,text2
      enddo
      do n=1,numq
        text1='bs      '
        text2='bndry source/sink of          '
        write(text1( 3: 5),156) qname(n)
        write(text2(22:24),156) qname(n)
        write(50,309) text1,text2
      enddo
    endif
    write(50,310)

156   format(a3)
301   format('dset ^',a70)
302   format('undef -9999.')
303   format('title ctl file for stats.dat')
304   format('xdef 1 linear 1 1')
305   format('ydef 1 linear 1 1')
306   format('zdef 1 linear 1 1')
307   format('tdef ',i10,' linear ',a15,' ',i5,'MN')
308   format('vars ',i6)
309   format(a8,' 1 99 ',a30)
310   format('endvars')

      close(unit=50)

      return
      end subroutine write_statsctl


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine writeout(dt,fnum,nwrite,qname,xh,xf,uf,yh,yf,vf,xfref,yfref,                  &
                        rds,sigma,rdsf,sigmaf,zh,zf,mf,                                        &
                        pi0,prs0,rho0,rr0,rf0,rrf0,th0,qv0,u0,v0,                              &
                        zs,rgzu,rgzv,rain,sws,svs,sps,srs,sgs,sus,shs,thflux,qvflux,           &
                        rxh,arh1,arh2,uh,ruh,rxf,arf1,arf2,vh,rvh,mh,rmf,rr,rf,                &
                        gz,rgz,gzu,gzv,gxu,gyv,dzdx,dzdy,doimpl,c1,c2,                         &
                        cd,ch,cq,dum1,dum2,dum3,dum4,dum5,dum6,t11,t12,t13,t22,t23,t33,        &
                        rho,prs,dbz ,rru,ua,u3d,dumu,rrv,va,v3d,dumv,rrw,wa,w3d,dumw,ppi,tha,  &
                        nm,defv,defh,dissten,thpten,qvpten,qcpten,qipten,upten,vpten,          &
                        lu_index,xland,mavail,tsk,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw,tslb,  &
                        qa,kmh,kmv,khh,khv,tkea,swten,lwten,                                   &
                        radsw,rnflx,radswnet,radlwin,dsr,olr,pta,                              &
                        num_soil_layers,u10,v10,t2,q2,znt,ust,hpbl,zol,mol,br,                 &
                        dat1,dat2,dat3)
#ifdef NETCDF
      use netcdf
#endif
      implicit none

      !----------------------------------------------------------
      ! This subroutine organizes writeouts for GrADS-format and
      ! netcdf-format output.
      !----------------------------------------------------------

      include 'input.incl'
      include 'constants.incl'
#ifdef MPI
      include 'mpif.h'
#endif

      real, intent(inout) :: dt
      integer, intent(in) :: fnum,nwrite
      character*3, dimension(maxq), intent(in) :: qname
      real, dimension(ib:ie), intent(in) :: xh
      real, dimension(ib:ie+1), intent(in) :: xf,uf
      real, dimension(jb:je), intent(in) :: yh
      real, dimension(jb:je+1), intent(in) :: yf,vf
      real, intent(in), dimension(-2:nx+4) :: xfref
      real, intent(in), dimension(-2:ny+4) :: yfref
      real, dimension(kb:ke), intent(in) :: rds,sigma
      real, dimension(kb:ke+1), intent(in) :: rdsf,sigmaf
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: zh
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: zf,mf
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: pi0,prs0,rho0,rr0,rf0,rrf0,th0,qv0
      real, dimension(itb:ite,jtb:jte), intent(in) :: zs,rgzu,rgzv
      real, dimension(ib:ie,jb:je,nrain), intent(in) :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, dimension(ib:ie,jb:je), intent(in) :: thflux,qvflux,cd,ch,cq
      real, intent(in), dimension(ib:ie) :: rxh,arh1,arh2,uh,ruh
      real, intent(in), dimension(ib:ie+1) :: rxf,arf1,arf2
      real, intent(in), dimension(jb:je) :: vh,rvh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: mh
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: rmf
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rr,rf
      real, intent(in), dimension(itb:ite,jtb:jte) :: gz,rgz,gzu,gzv
      real, intent(in), dimension(itb:ite,jtb:jte,ktb:kte) :: gxu,gyv
      real, intent(in), dimension(itb:ite,jtb:jte) :: dzdx,dzdy
      logical, intent(in), dimension(ib:ie,jb:je) :: doimpl
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: c1,c2
      real, dimension(ib:ie,jb:je,kb:ke), intent(inout) :: dum1,dum2,dum3,dum4,dum5,dum6
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: t11,t12,t13,t22,t23,t33
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: rho,prs,dbz
      real, dimension(ib:ie+1,jb:je,kb:ke), intent(in) :: u0,ua,u3d
      real, dimension(ib:ie+1,jb:je,kb:ke), intent(inout) :: rru,dumu
      real, dimension(ib:ie,jb:je+1,kb:ke), intent(in) :: v0,va,v3d
      real, dimension(ib:ie,jb:je+1,kb:ke), intent(inout) :: rrv,dumv
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: wa,w3d
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(inout) :: rrw,dumw
      real, dimension(ib:ie,jb:je,kb:ke), intent(in) :: ppi,tha
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: nm,defv,defh,dissten
      real, dimension(ibb:ieb,jbb:jeb,kbb:keb), intent(in) :: thpten,qvpten,qcpten,qipten,upten,vpten
      integer, dimension(ibl:iel,jbl:jel), intent(in) :: lu_index
      real, dimension(ib:ie,jb:je), intent(in) :: tsk
      real, dimension(ibl:iel,jbl:jel), intent(in) :: xland,mavail,tmn,tml,hml,huml,hvml,hfx,qfx,gsw,glw
      real, dimension(ibl:iel,jbl:jel,num_soil_layers), intent(in) :: tslb
      real, dimension(ibm:iem,jbm:jem,kbm:kem,numq), intent(in) :: qa
      real, dimension(ibc:iec,jbc:jec,kbc:kec), intent(in) :: kmh,kmv,khh,khv
      real, dimension(ibt:iet,jbt:jet,kbt:ket), intent(in) :: tkea
      real, dimension(ibr:ier,jbr:jer,kbr:ker), intent(in) :: swten,lwten
      real, dimension(ni,nj), intent(in) :: radsw,rnflx,radswnet,radlwin,dsr,olr
      real, dimension(ibp:iep,jbp:jep,kbp:kep,npt), intent(in) :: pta
      integer, intent(in) :: num_soil_layers
      real, dimension(ibl:iel,jbl:jel), intent(in) :: u10,v10,t2,q2,hpbl,zol,mol,br
      real, dimension(ib:ie,jb:je), intent(in) :: znt,ust
      real, intent(inout), dimension(ni+1,nj+1) :: dat1
      real, intent(inout), dimension((ni+1)*nodex,(nj+1)*nodey) :: dat2
      real, intent(inout), dimension(ni+1,nj+1,numprocs) :: dat3

      integer i,j,k,n,irec
      integer :: ncid,time_index
      real :: ubar,vbar,wbar,cc1,cc2,r1,r2
      character*8 :: chid
#ifdef MPI
      character*80 sname,uname,vname,wname
      integer status(MPI_STATUS_SIZE)
      logical :: doit
      integer :: reqs
      integer, parameter :: nlim = 1000
#endif

!--------------------------------------------------------------
!  writeout data on scalar-points

  IF(output_format.eq.1)THEN
  ! grads stuff:
  IF(output_filetype.eq.1.and.myid.eq.0)THEN
    ! one output file:
    if(dowr) write(outfile,*)
    if(s_out.ge.1)then
      if(fnum.eq.51)then
        string(totlen+1:totlen+22) = '_s.dat                '
      elseif(fnum.eq.71)then
        string(totlen+1:totlen+22) = '_i.dat                '
      endif
      if(dowr) write(outfile,*) string
      open(unit=fnum,file=string,form='unformatted',access='direct',   &
           recl=(nx*ny*4),status='unknown')
      irec=1+(nwrite-1)*( sout2d + nk*sout3d )
    endif
    if(u_out.ge.1.and.fnum.ne.71)then
      string(totlen+1:totlen+22) = '_u.dat                '
      if(dowr) write(outfile,*) string
      open(unit=52,file=string,form='unformatted',access='direct',   &
           recl=((nx+1)*ny*4),status='unknown')
    endif
    if(v_out.ge.1.and.fnum.ne.71)then
      string(totlen+1:totlen+22) = '_v.dat                '
      if(dowr) write(outfile,*) string
      open(unit=53,file=string,form='unformatted',access='direct',   &
           recl=(nx*(ny+1)*4),status='unknown')
    endif
    if(w_out.ge.1.and.fnum.ne.71)then
      string(totlen+1:totlen+22) = '_w.dat                '
      if(dowr) write(outfile,*) string
      open(unit=54,file=string,form='unformatted',access='direct',   &
           recl=(nx*ny*4),status='unknown')
    endif
  ELSEIF(output_filetype.eq.2.and.myid.eq.0)THEN
    ! one output file per output time:
    if(s_out.ge.1)then
      if(fnum.eq.51)then
        string(totlen+1:totlen+22) = '_XXXXXX_s.dat         '
      elseif(fnum.eq.71)then
        string(totlen+1:totlen+22) = '_XXXXXX_i.dat         '
      endif
      write(string(totlen+2:totlen+7),102) nwrite
102   format(i6.6)
      if(dowr) write(outfile,*) string
      open(unit=fnum,file=string,form='unformatted',access='direct',   &
           recl=(nx*ny*4),status='unknown')
      irec=1
    endif
    if(u_out.ge.1.and.fnum.ne.71)then
      string(totlen+1:totlen+22) = '_XXXXXX_u.dat         '
      write(string(totlen+2:totlen+7),102) nwrite
      if(dowr) write(outfile,*) string
      open(unit=52,file=string,form='unformatted',access='direct',   &
           recl=((nx+1)*ny*4),status='unknown')
    endif
    if(v_out.ge.1.and.fnum.ne.71)then
      string(totlen+1:totlen+22) = '_XXXXXX_v.dat         '
      write(string(totlen+2:totlen+7),102) nwrite
      if(dowr) write(outfile,*) string
      open(unit=53,file=string,form='unformatted',access='direct',   &
           recl=(nx*(ny+1)*4),status='unknown')
    endif
    if(w_out.ge.1.and.fnum.ne.71)then
      string(totlen+1:totlen+22) = '_XXXXXX_w.dat         '
      write(string(totlen+2:totlen+7),102) nwrite
      if(dowr) write(outfile,*) string
      open(unit=54,file=string,form='unformatted',access='direct',   &
           recl=(nx*ny*4),status='unknown')
    endif
  ELSEIF(output_filetype.eq.3)THEN
    ! one output file per output time AND one output file per processor:
    ! (MPI only)
#ifdef MPI
      irec=1

      sname = '                                                                                '
      uname = '                                                                                '
      vname = '                                                                                '
      wname = '                                                                                '

    if(strlen.gt.0)then
      sname(1:strlen) = output_path(1:strlen)
      uname(1:strlen) = output_path(1:strlen)
      vname(1:strlen) = output_path(1:strlen)
      wname(1:strlen) = output_path(1:strlen)
    endif

      sname(strlen+1:strlen+baselen) = output_basename(1:baselen)
      uname(strlen+1:strlen+baselen) = output_basename(1:baselen)
      vname(strlen+1:strlen+baselen) = output_basename(1:baselen)
      wname(strlen+1:strlen+baselen) = output_basename(1:baselen)

      if(fnum.eq.51)then
        sname(totlen+1:totlen+22) = '_XXXXXX_YYYYYY_s.dat  '
      elseif(fnum.eq.71)then
        sname(totlen+1:totlen+22) = '_XXXXXX_YYYYYY_i.dat  '
      endif
      uname(totlen+1:totlen+22) = '_XXXXXX_YYYYYY_u.dat  '
      vname(totlen+1:totlen+22) = '_XXXXXX_YYYYYY_v.dat  '
      wname(totlen+1:totlen+22) = '_XXXXXX_YYYYYY_w.dat  '

      write(sname(totlen+2:totlen+ 7),100) myid
      write(sname(totlen+9:totlen+14),100) nwrite

      write(uname(totlen+2:totlen+ 7),100) myid
      write(uname(totlen+9:totlen+14),100) nwrite

      write(vname(totlen+2:totlen+ 7),100) myid
      write(vname(totlen+9:totlen+14),100) nwrite

      write(wname(totlen+2:totlen+ 7),100) myid
      write(wname(totlen+9:totlen+14),100) nwrite

100   format(i6.6)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  myid,sname=',myid,'   ',sname
      open(unit=fnum,file=sname,                   &
           form='unformatted',access='direct',   &
           recl=(ni*nj*4),status='unknown')

      if(u_out.ge.1.and.fnum.ne.71)then
        if(dowr) write(outfile,*) '  myid,uname=',myid,'   ',uname
        open(unit=52,file=uname,                   &
             form='unformatted',access='direct',   &
             recl=((ni+1)*nj*4),status='unknown')
      endif

      if(v_out.ge.1.and.fnum.ne.71)then
        if(dowr) write(outfile,*) '  myid,vname=',myid,'   ',vname
        open(unit=53,file=vname,                   &
             form='unformatted',access='direct',   &
             recl=(ni*(nj+1)*4),status='unknown')
      endif

      if(w_out.ge.1.and.fnum.ne.71)then
        if(dowr) write(outfile,*) '  myid,wname=',myid,'   ',wname
        open(unit=54,file=wname,                   &
             form='unformatted',access='direct',   &
             recl=(ni*nj*4),status='unknown')
      endif
      !  limit to "nlim" writes at a time:
      IF( numprocs.gt.nlim )THEN
        doit = .false.
        IF( myid.ge.nlim )THEN
          call MPI_IRECV(doit,1,mpi_logical,myid-nlim,999999,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ENDIF
      ENDIF
#else
    print *,'  output_filetype = ',output_filetype
    print *,'  This option is only available for MPI runs '
    print *,'  Stopping cm1 .... '
    call stopcm1
#endif
  ENDIF ! endif for outout_filetype
  ENDIF ! endif for output_format=1
#ifdef NETCDF
  IF(output_format.eq.2)THEN
    ! netcdf stuff:
    if( output_filetype.eq.3 .or. myid.eq.0 )                                        &
    call netcdf_prelim(nwrite,ncid,time_index,qname,xh,xf,yh,yf,                     &
                       xfref,yfref,sigma,sigmaf,zs,zh,zf,                            &
                       dum1(ib,jb,kb),dum2(ib,jb,kb),dum3(ib,jb,kb),dum4(ib,jb,kb),  &
                       dum5(ib,jb,kb),dat2(1,1),dat2(1,2))
  ENDIF
#endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--- 2D vars:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(output_rain.eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,rain(ib,jb,1),'rain    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_sws .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,sws(ib,jb,1),'sws     ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_svs .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,svs(ib,jb,1),'svs     ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_sps .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,sps(ib,jb,1),'sps     ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_srs .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,srs(ib,jb,1),'srs     ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_sgs .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,sgs(ib,jb,1),'sgs     ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_sus .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,sus(ib,jb,1),'sus     ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_shs .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,shs(ib,jb,1),'shs     ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    IF(nrain.eq.2)THEN
      if(output_rain.eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,rain(ib,jb,2),'rain2   ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_sws .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,sws(ib,jb,2),'sws2    ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_svs .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,svs(ib,jb,2),'svs2    ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_sps .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,sps(ib,jb,2),'sps2    ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_srs .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,srs(ib,jb,2),'srs2    ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_sgs .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,sgs(ib,jb,2),'sgs2    ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_sus .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,sus(ib,jb,2),'sus2    ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_shs .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,shs(ib,jb,2),'shs2    ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    ENDIF
      if(output_uh.eq.1)then
        ! get height AGL:
        if( terrain_flag )then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            dum3(i,j,k) = zh(i,j,k)-zs(i,j)
            dumw(i,j,k) = zf(i,j,k)-zs(i,j)
          enddo
          enddo
          enddo
        else
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk+1
          do j=1,nj
          do i=1,ni
            dum3(i,j,k) = zh(i,j,k)
            dumw(i,j,k) = zf(i,j,k)
          enddo
          enddo
          enddo
        endif
        call calcuh(uf,vf,dum3,dumw,ua,va,wa,dum1(ib,jb,1),dum2,dum5,dum6, &
                    zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf)
        call writeo(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'uh      ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_coldpool.eq.1)then
        call calccpch(zh,zf,th0,qv0,dum1(ib,jb,1),dum1(ib,jb,2),tha,qa)
        call writeo(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'cpc     ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,dum1(ib,jb,2),'cph     ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_sfcflx.eq.1)                                        &
        call writeo(ni,nj,1,1,nx,ny,thflux(ib,jb),'thflux  ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_sfcflx.eq.1)                                        &
        call writeo(ni,nj,1,1,nx,ny,qvflux(ib,jb),'qvflux  ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_sfcflx.eq.1)then
        call writeo(ni,nj,1,1,nx,ny,cd(ib,jb),'cd      ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_sfcflx.eq.1)                                        &
        call writeo(ni,nj,1,1,nx,ny,ch(ib,jb),'ch      ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_sfcflx.eq.1)                                        &
        call writeo(ni,nj,1,1,nx,ny,cq(ib,jb),'cq      ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_sfcflx.eq.1)                                        &
        call writeo(ni,nj,1,1,nx,ny,tsk(ib,jb),'tsk     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_zs  .eq.1)                                          &
        call writeo(ni,nj,1,1,nx,ny,zs(ib,jb),'zs      ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_dbz   .eq.1)then
        call calccref(dum1(ib,jb,1),dbz)
        call writeo(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'cref    ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if(output_sfcparams.eq.1)then
        call writeo(ni,nj,1,1,nx,ny,xland(ib,jb),'xland   ',        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1)=float(lu_index(i,j))
        enddo
        enddo
        call writeo(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'lu      ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,mavail(ib,jb),'mavail  ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.oceanmodel.eq.2))then
        call writeo(ni,nj,1,1,nx,ny,tmn(ib,jb),'tmn     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,hfx(ib,jb),'hfx     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,qfx(ib,jb),'qfx     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,gsw(ib,jb),'gsw     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,glw(ib,jb),'glw     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    endif

    if((output_sfcparams.eq.1).and.(sfcmodel.eq.2))then
        call writeo(ni,nj,1,1,nx,ny,tslb(ib,jb,1),'tslb1   ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,tslb(ib,jb,2),'tslb2   ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,tslb(ib,jb,3),'tslb3   ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,tslb(ib,jb,4),'tslb4   ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,tslb(ib,jb,5),'tslb5   ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    endif

      if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then
        call writeo(ni,nj,1,1,nx,ny,tml(ib,jb),'tml     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,hml(ib,jb),'hml     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,huml(ib,jb),'huml    ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,hvml(ib,jb),'hvml    ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if( output_radten.eq.1 )then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radsw(i,j)
        enddo
        enddo
        call writeo(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radsw   ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = rnflx(i,j)
        enddo
        enddo
        call writeo(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'rnflx   ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radswnet(i,j)
        enddo
        enddo
        call writeo(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radswnet',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radlwin(i,j)
        enddo
        enddo
        call writeo(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'radlwin ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = olr(i,j)
        enddo
        enddo
        call writeo(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'olr     ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = dsr(i,j)
        enddo
        enddo
        call writeo(ni,nj,1,1,nx,ny,dum1(ib,jb,1),'dsr     ',       &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      IF(output_sfcdiags.eq.1)THEN
        call writeo(ni,nj,1,1,nx,ny,u10(ib,jb),'u10     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,v10(ib,jb),'v10     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,t2(ib,jb),'t2      ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,q2(ib,jb),'q2      ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,znt(ib,jb),'znt     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,ust(ib,jb),'ust     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,hpbl(ib,jb),'hpbl    ',         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,zol(ib,jb),'zol     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,mol(ib,jb),'mol     ',          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,1,nx,ny,br(ib,jb),'br      ',           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--- 3D vars below here:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      dum1=zh
      if(fnum.eq.71)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=sigma(k)-zs(i,j)
        enddo
        enddo
        enddo
      endif
      if(output_zh  .eq.1)                                          &
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'zh      ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_th  .eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=th0(i,j,k)+tha(i,j,k)
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'th      ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_thpert .eq.1)then
        dum1=tha
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'thpert  ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_prs    .eq.1)then
        dum1=prs
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'prs     ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_prspert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=prs(i,j,k)-p00*(pi0(i,j,k)**cpdrd)
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'prspert ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_pi.eq.1)then  
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=pi0(i,j,k)+ppi(i,j,k)
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'pi      ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_pipert .eq.1)then
        dum1=ppi
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'pipert  ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_rho    .eq.1)then
        dum1=rho
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'rho     ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_rhopert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=rho(i,j,k)-rho0(i,j,k)
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'rhopert ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(iptra.eq.1)then
        chid = 'pt      '
        do n=1,npt
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
          do k=1,nk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k)=pta(i,j,k,n)
          enddo
          enddo
          enddo
          if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
          write(chid(3:4),111) n
111       format(i2.2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),chid,            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        enddo
      endif
    IF(imoist.eq.1)THEN
      if(output_qv    .eq.1)then
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=qa(i,j,k,nqv)
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'qv      ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_qvpert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=qa(i,j,k,nqv)-qv0(i,j,k)
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'qvpert  ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_q.eq.1)then
        chid = '        '
        do n=1,numq
          if(n.ne.nqv)then
            do k=1,nk
            do j=1,nj
            do i=1,ni
              dum1(i,j,k)=qa(i,j,k,n)
            enddo
            enddo
            enddo
            if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
            write(chid(1:3),110) qname(n)
110         format(a3)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),chid,            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
          endif
        enddo
      endif
      if(output_dbz   .eq.1)then
        dum1=dbz
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'dbz     ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
    ENDIF  ! endif for imoist=1

      if(output_buoyancy.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k,n)
        do k=1,nk
          do j=1,nj
          do i=1,ni
            dum1(i,j,k)=g*tha(i,j,k)/th0(i,j,k)
          enddo
          enddo
          IF(imoist.eq.1)THEN
            do j=1,nj
            do i=1,ni
              dum1(i,j,k)=dum1(i,j,k)+g*repsm1*(qa(i,j,k,nqv)-qv0(i,j,k))
            enddo
            enddo
            do n=nql1,nql2
              do j=1,nj
              do i=1,ni
                dum1(i,j,k)=dum1(i,j,k)-g*qa(i,j,k,n)
              enddo
              enddo
            enddo
            IF(iice.eq.1)THEN
            do n=nqs1,nqs2
              do j=1,nj
              do i=1,ni
                dum1(i,j,k)=dum1(i,j,k)-g*qa(i,j,k,n)
              enddo
              enddo
            enddo
            ENDIF
          ENDIF
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'buoyancy',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if(output_uinterp.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=0.5*(ua(i,j,k)+ua(i+1,j,k))
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'uinterp ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_vinterp.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=0.5*(va(i,j,k)+va(i,j+1,k))
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'vinterp ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_winterp.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni
          dum1(i,j,k)=0.5*(wa(i,j,k)+wa(i,j,k+1))
        enddo
        enddo
        enddo
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'winterp ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if(output_vort.eq.1)then
        call calcvort(xh,xf,uf,vf,zh,mf,zf,zs,rgzu,rgzv,rds,sigma,rdsf,sigmaf,  &
                      ua,va,wa,dum2,dum3,dum4,dum1,dum5,dum6)
        dum1=dum2
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'xvort   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        dum1=dum3
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'yvort   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        dum1=dum4
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'zvort   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if(output_basestate.eq.1)then
        dum1=pi0
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'pi0     ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        dum1=th0
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'th0     ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        dum1=prs0
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'prs0    ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        dum1=qv0
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'qv0     ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if(output_pblten.eq.1)then
        dum1=thpten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'thpten  ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        dum1=qvpten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'qvpten  ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        dum1=qcpten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'qcpten  ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        dum1=qipten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'qipten  ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        dum1=upten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'upten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        dum1=vpten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'vpten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if( output_radten.eq.1 )then
        dum1=swten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'swten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        dum1=lwten
        if(fnum.eq.71) call zinterp(sigma,zs,zh,dum1,dum2)
        call writeo(ni,nj,1,nk,nx,ny,dum1(ib,jb,1),'lwten   ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fnum, &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

!--------------------------------------------------------------

    IF(output_impdiften.eq.1)THEN
      IF(.not.terrain_flag)THEN
        ! without terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        DO k=1,nk
          do j=0,nj+1
          do i=0,ni+2
            rru(i,j,k)=rho0(1,1,k)*u3d(i,j,k)
          enddo
          enddo
          do j=0,nj+2
          do i=0,ni+1
            rrv(i,j,k)=rho0(1,1,k)*v3d(i,j,k)
          enddo
          enddo
          IF(k.eq.1)THEN
            do j=0,nj+1
            do i=0,ni+1
              rrw(i,j,   1) = 0.0
              rrw(i,j,nk+1) = 0.0
            enddo
            enddo
          ELSE
            do j=0,nj+1
            do i=0,ni+1
              rrw(i,j,k)=rf0(1,1,k)*w3d(i,j,k)
            enddo
            enddo
          ENDIF
        ENDDO
      ELSE
        ! with terrain:
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        DO k=1,nk
          do j=0,nj+1
          do i=0,ni+2
            rru(i,j,k)=0.5*(rho0(i-1,j,k)+rho0(i,j,k))*u3d(i,j,k)*rgzu(i,j)
          enddo
          enddo
          do j=0,nj+2
          do i=0,ni+1
            rrv(i,j,k)=0.5*(rho0(i,j-1,k)+rho0(i,j,k))*v3d(i,j,k)*rgzv(i,j)
          enddo
          enddo
        ENDDO
!$omp parallel do default(shared)  &
!$omp private(i,j,k,r1,r2)
        DO k=1,nk
          IF(k.eq.1)THEN
            do j=0,nj+1
            do i=0,ni+1  
              rrw(i,j,   1) = 0.0
              rrw(i,j,nk+1) = 0.0
            enddo
            enddo
          ELSE 
            r2 = (sigmaf(k)-sigma(k-1))*rds(k)
            r1 = 1.0-r2
            do j=0,nj+1
            do i=0,ni+1
              rrw(i,j,k)=rf0(i,j,k)*w3d(i,j,k)                                  &
                        +0.5*( ( r2*(rru(i,j,k  )+rru(i+1,j,k  ))               &
                                +r1*(rru(i,j,k-1)+rru(i+1,j,k-1)) )*dzdx(i,j)   &
                              +( r2*(rrv(i,j,k  )+rrv(i,j+1,k  ))               &
                                +r1*(rrv(i,j,k-1)+rrv(i,j+1,k-1)) )*dzdy(i,j)   &
                             )*(sigmaf(k)-zt)*gz(i,j)*rzt
            enddo
            enddo
          ENDIF
        ENDDO
      ENDIF  ! endif for terrain_flag
    ENDIF    ! endif for output_impdiften=1

!--------------------------------------------------------------
!  writeout data on u-points

    IF(fnum.ne.71)THEN
      irec=1+(nwrite-1)*nk*u_out
      if(output_filetype.ge.2) irec=1

      if(output_u    .eq.1)                                         &
        call writeo(ni+1,nj,1,nk,nx+1,ny,ua(ib,jb,1),'u       ',    &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,52,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)


      if(output_upert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj
        do i=1,ni+1
          dumu(i,j,k)=ua(i,j,k)-u0(i,j,k)
        enddo
        enddo
        enddo
        call writeo(ni+1,nj,1,nk,nx+1,ny,dumu(ib,jb,1),'upert   ',  &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,52,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if(output_basestate.eq.1)                                     &
        call writeo(ni+1,nj,1,nk,nx+1,ny,u0(ib,jb,1),'u0      ',    &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,52,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)

      if(output_turbten.eq.1)then
        dumu = 0.0
        call turbu(dt,xh,ruh,xf,rxf,arf1,arf2,uf,vh,mh,mf,rmf,rho,rf,  &
                   zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gxu,     &
                   dum1,dum2,dum3,dum4,dum5,ua,dumu,wa,t11,t12,t13,t22,kmv,doimpl)
        call writeo(ni+1,nj,1,nk,nx+1,ny,dumu(ib,jb,1),'ftu     ',  &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,52,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if(output_impdiften.eq.1)then
        dumu = 0.0
        call impldiffu(uf,vh,arh1,arh2,arf1,arf2,mh,rdsf,gzu,rho0,rr0,dum1,dum2,dum3,rru,rrv,rrw,u3d,dumu)
        call writeo(ni+1,nj,1,nk,nx+1,ny,dumu(ib,jb,1),'fdu     ',  &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,52,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

    ENDIF

!--------------------------------------------------------------
!  writeout data on v-points

    IF(fnum.ne.71)THEN
      irec=1+(nwrite-1)*nk*v_out
      if(output_filetype.ge.2) irec=1

      if(output_v    .eq.1)                                         &
        call writeo(ni,nj+1,1,nk,nx,ny+1,va(ib,jb,1),'v       ',    &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,53,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)

      if(output_vpert.eq.1)then
!$omp parallel do default(shared)  &
!$omp private(i,j,k)
        do k=1,nk
        do j=1,nj+1
        do i=1,ni
          dumv(i,j,k)=va(i,j,k)-v0(i,j,k)
        enddo
        enddo
        enddo
        call writeo(ni,nj+1,1,nk,nx,ny+1,dumv(ib,jb,1),'vpert   ',  &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,53,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if(output_basestate.eq.1)                                     &
        call writeo(ni,nj+1,1,nk,nx,ny+1,v0(ib,jb,1),'v0      ',    &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,53,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)

      if(output_turbten.eq.1)then
        dumv = 0.0
        call turbv(dt,xh,rxh,arh1,arh2,uh,xf,rvh,vf,mh,mf,rho,rr,rf,   &
                   zs,gz,rgz,gzu,gzv,rds,sigma,rdsf,sigmaf,gyv,  &
                   dum1,dum2,dum3,dum4,dum5,va,dumv,wa,t12,t22,t23,kmv,doimpl)
        call writeo(ni,nj+1,1,nk,nx,ny+1,dumv(ib,jb,1),'ftv     ',  &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,53,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if(output_impdiften.eq.1)then
        dumv = 0.0
        call impldiffv(uh,vf,arh1,arh2,mh,rdsf,gzv,rho0,rr0,dum1,dum2,dum3,rru,rrv,rrw,v3d,dumv)
        call writeo(ni,nj+1,1,nk,nx,ny+1,dumv(ib,jb,1),'fdv     ',  &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,53,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

    ENDIF

!--------------------------------------------------------------
!  writeout data on w-points

    IF(fnum.ne.71)THEN
      irec=1+(nwrite-1)*(nk+1)*w_out
      if(output_filetype.ge.2) irec=1

      if(output_w  .eq.1)                                           &
        call writeo(ni,nj,1,nk+1,nx,ny,wa(ib,jb,1),'w       ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,54,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_tke.eq.1.and.iturb.eq.1)                            &
        call writeo(ni,nj,1,nk+1,nx,ny,tkea(ib,jb,1),'tke     ',    &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,54,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_km .eq.1)                                           &
        call writeo(ni,nj,1,nk+1,nx,ny,kmh(ib,jb,1),'kmh     ',     &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,54,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_km .eq.1)                                           &
        call writeo(ni,nj,1,nk+1,nx,ny,kmv(ib,jb,1),'kmv     ',     &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,54,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_kh .eq.1)                                           &
        call writeo(ni,nj,1,nk+1,nx,ny,khh(ib,jb,1),'khh     ',     &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,54,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      if(output_kh .eq.1)                                           &
        call writeo(ni,nj,1,nk+1,nx,ny,khv(ib,jb,1),'khv     ',     &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,54,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)


      if(output_dissten.eq.1)then
        call writeo(ni,nj,1,nk+1,nx,ny,dissten(ib,jb,1),'dissten ', &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,54,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_nm.eq.1)then
        call writeo(ni,nj,1,nk+1,nx,ny,nm(ib,jb,1),'nm      ',      &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,54,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(output_def.eq.1)then
        call writeo(ni,nj,1,nk+1,nx,ny,defv(ib,jb,1),'defv    ',    &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,54,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writeo(ni,nj,1,nk+1,nx,ny,defh(ib,jb,1),'defh    ',    &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,54,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if(output_turbten.eq.1)then
        dumw = 0.0
        call turbw(dt,xh,rxh,arh1,arh2,uh,xf,vh,mh,mf,rho,rf,gz,rgzu,rgzv,rds,sigma,   &
                   dum1,dum2,dum3,dum4,dum5,wa,dumw,t13,t23,t33,t22,kmv,doimpl)
        call writeo(ni,nj,1,nk+1,nx,ny,dumw(ib,jb,1),'ftw     ',    &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,54,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

      if(output_impdiften.eq.1)then
        dumw = 0.0
        call impldiffw(uh,vh,arh1,arh2,mf,gz,rds,rrf0,dum1,dum2,dum3,rru,rrv,rrw,w3d,dumw,c1,c2)
        call writeo(ni,nj,1,nk+1,nx,ny,dumw(ib,jb,1),'fdw     ',    &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,54,   &
                    ncid,time_index,output_format,output_filetype,  &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

    ENDIF

!---------------------------------------------------------------
#ifdef MPI
      !  limit to "nlim" writes at a time:
      IF( numprocs.gt.nlim )THEN
        doit = .true.
        IF( myid+nlim .le. (numprocs-1) )THEN
          call MPI_ISEND(doit,1,mpi_logical,myid+nlim,999999,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ENDIF
      ENDIF
#endif
!--------------------------------------------------------------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) 'Done Writing Data to File: nwrite=',nwrite
      if(dowr) write(outfile,*)

    IF(output_format.eq.1)THEN
      close(unit=fnum)
      if(u_out.ge.1.and.fnum.ne.71)then
        close(unit=52)
      endif
      if(v_out.ge.1.and.fnum.ne.71)then
        close(unit=53)
      endif
      if(w_out.ge.1.and.fnum.ne.71)then
        close(unit=54)
      endif
#ifdef NETCDF
    ELSEIF(output_format.eq.2)THEN
      if( output_filetype.eq.3 .or. myid.eq.0 ) &
      call disp_err( nf90_close(ncid) , .true. )
#endif
    ENDIF

#ifdef MPI
      if(timestats.ge.1)then
        ! this is needed for proper accounting of timing:
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      endif
#endif

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine writeo(numi,numj,numk1,numk2,nxr,nyr,var,aname,             &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,irec,fileunit,  &
                      ncid,time_index,output_format,output_filetype,       &
                      dat1,dat2,dat3,ppnode)
    implicit none

    !-------------------------------------------------------------------
    ! This subroutine collects data (from other processors if this is a
    ! MPI run) and does the actual writing to disk.
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr
    integer, intent(in) :: ppnode
    real, intent(in), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    character*8, intent(in) :: aname
    integer, intent(in) :: ni,nj,ngxy,myid,numprocs,nodex,nodey,fileunit
    integer, intent(inout) :: irec
    integer, intent(in) :: ncid,time_index,output_format,output_filetype
    real, intent(inout), dimension(numi,numj) :: dat1
    real, intent(inout), dimension(numi*nodex,numj*nodey) :: dat2
    real, intent(inout), dimension(numi,numj,0:numprocs-1) :: dat3

    integer :: i,j,k,ii,jj
#ifdef DP
    real*4, dimension(nxr,nyr) :: vout
#endif
#ifdef MPI
    include 'mpif.h'
    integer :: n,nn,nnn,index,fooi,fooj,proc,ierr,reqs
    integer, dimension(numprocs-1) :: reqt
    integer status(MPI_STATUS_SIZE)
    integer :: mynode,nodemaster
#endif

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  this section writes one file:
!    (for MPI, use processor 0 to collect data)

#ifdef MPI
    mynode = int( float(myid)/float(ppnode) )
    nodemaster = mynode * ppnode
#endif

  IF(output_filetype.le.2)THEN

    DO k=numk1,numk2

#ifdef MPI
      !-------------------- MPI section --------------------!
      IF(myid.ne.nodemaster)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat1(i,j)=var(i,j,k)
        enddo
        enddo
        call MPI_ISEND(dat1(1,1),numi*numj,MPI_REAL,nodemaster,k,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,status,ierr)
      ELSE
        ! begin nodemaster section ...
        do proc=myid+1,myid+(ppnode-1)
          call MPI_IRECV(dat3(1,1,proc),numi*numj,MPI_REAL,proc,k,MPI_COMM_WORLD,reqt(proc),ierr)
        enddo
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat3(i,j,myid)=var(i,j,k)
        enddo
        enddo
        nn = 1
        do while( nn.le.(ppnode-1) )
          nn = nn + 1
          call mpi_waitany(ppnode-1,reqt(myid+1),index,status,ierr)
        enddo
        IF(myid.ne.0)THEN
          call MPI_ISEND(dat3(1,1,myid),numi*numj*ppnode,MPI_REAL,0,10000+k,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ELSE
          ! begin processor 0 section ...
          nn = 1
          do proc = ppnode , nodex*nodey-ppnode , ppnode
            call MPI_IRECV(dat3(1,1,proc),numi*numj*ppnode,MPI_REAL,proc,10000+k,MPI_COMM_WORLD,reqt(nn),ierr)
            nn = nn + 1
          enddo
          do nnn=0,ppnode-1
            proc = nnn
            fooj = proc / nodex + 1
            fooi = proc - (fooj-1)*nodex  + 1
            fooi = (fooi-1)*ni
            fooj = (fooj-1)*nj
!$omp parallel do default(shared)  &
!$omp private(i,j)
            do j=1,numj
            do i=1,numi
              dat2(fooi+i,fooj+j) = dat3(i,j,nnn)
            enddo
            enddo
          enddo
          nn = 1
          do while( nn.le.nodex*nodey/ppnode-1 )
            nn = nn + 1
            call mpi_waitany(nodex*nodey/ppnode-1,reqt(1),index,status,ierr)
            do nnn=index*ppnode,(index+1)*ppnode-1
              proc = nnn
              fooj = proc / nodex + 1
              fooi = proc - (fooj-1)*nodex  + 1
              fooi = (fooi-1)*ni
              fooj = (fooj-1)*nj
!$omp parallel do default(shared)  &
!$omp private(i,j)
              do j=1,numj
              do i=1,numi
                dat2(fooi+i,fooj+j) = dat3(i,j,nnn)
              enddo
              enddo
            enddo
          enddo
          ! ... end processor 0 section
        ENDIF
        ! ... end nodemaster section
      ENDIF
#else
      !-------------------- non-MPI section --------------------!
!$omp parallel do default(shared)   &
!$omp private(i,j)
      do j=1,numj
      do i=1,numi
        dat2(i,j)=var(i,j,k)
      enddo
      enddo
#endif

        !-------------------- done MPI/non-MPI sections --------------------!

      IF(myid.eq.0)THEN
        ! only processor 0 writes:
        IF(output_format.eq.1)THEN
          ! grads format:
#ifdef DP
          do j=1,nyr
          do i=1,nxr
            vout(i,j) = dat2(i,j)
          enddo
          enddo
          write(fileunit,rec=irec) ((vout(i,j),i=1,nxr),j=1,nyr)
#else
          write(fileunit,rec=irec) ((dat2(i,j),i=1,nxr),j=1,nyr)
#endif
          irec=irec+1
#ifdef NETCDF
        ELSEIF(output_format.eq.2)THEN
          ! netcdf format:
          IF( numi*nodex .ne. nxr )THEN
            ! u-staggered vars:
            ii = 1
            jj = 1
            do j=1,nyr
            do i=1,nxr
              dat2(ii,jj) = dat2(i,j)
              ii = ii+1
              if( ii.gt.(numi*nodex) )then
                ii = 1
                jj = jj+1
              endif
            enddo
            enddo
          ENDIF
          if(numk1.eq.numk2)then
            call write2d_nc(aname,ncid,time_index,nxr,nyr,dat2(1,1))
          else
            call write3d_nc(aname,k,ncid,time_index,nxr,nyr,dat2(1,1))
          endif
#endif
        ENDIF
      ENDIF

    ENDDO  ! enddo for k-loop

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!  this section wites one output file per MPI process:
!  (for MPI runs only)

  ELSEIF(output_filetype.eq.3)THEN
    IF(output_format.eq.1)THEN
      ! grads format:
      DO k=numk1,numk2
        write(fileunit,rec=irec) ((var(i,j,k),i=1,numi),j=1,numj)
        irec=irec+1
      ENDDO
#ifdef NETCDF
    ELSEIF(output_format.eq.2)THEN
      ! netcdf format:
      DO k=numk1,numk2
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat1(i,j)=var(i,j,k)
        enddo
        enddo
        if(numk1.eq.numk2)then
          call write2d_nc(aname,ncid,time_index,numi,numj,dat1(1,1))
        else
          call write3d_nc(aname,k,ncid,time_index,numi,numj,dat1(1,1))
        endif
      ENDDO
#endif
    ENDIF
  ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      return
      end subroutine writeo


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine write_restart(nstep,nrec,prec,nwrite,nrst,nrad2d,num_soil_layers,  &
                               dt,mtime,stattim,taptim,rsttim,radtim,prcltim,       &
                               qbudget,asq,bsq,                                     &
                               rain,sws,svs,sps,srs,sgs,sus,shs,                    &
                               tsk,znt,ust,cd,ch,cq,u1,v1,w1,thflux,qvflux,         &
                               radbcw,radbce,radbcs,radbcn,                         &
                               rho,prs,ua,va,wa,ppi,tha,qa,tkea,swten,lwten,        &
                               radsw,rnflx,radswnet,radlwin,rad2d,                  &
                               effc,effi,effs,effr,effg,effis,                      &
                               lu_index,kpbl2d,psfc,u10,v10,w10,hfx,qfx,xland,      &
                               hpbl,wspd,psim,psih,gz1oz0,br,                       &
                               CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                       &
                               MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                    &
                               CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,       &
                               f2d,gsw,glw,chklowq,capg,snowc,tslb,                 &
                               tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,               &
                               qpten,qtten,qvten,qcten,pta,pdata,ploc,              &
                               dum1,radwebc,radsnbc,dat1,dat2,dat3)
      implicit none

      include 'input.incl'
      include 'constants.incl'

      !----------------------------------------------------------
      ! This subroutine organizes the writing of restart files
      !----------------------------------------------------------

      integer, intent(in) :: nstep,nrec,prec,nwrite,nrst
      integer, intent(in) :: nrad2d,num_soil_layers
      real, intent(in) :: dt
      double precision, intent(in) :: mtime,stattim,taptim,rsttim,radtim,prcltim
      real*8, intent(in), dimension(nbudget) :: qbudget
      real*8, intent(in), dimension(numq) :: asq,bsq
      real, intent(in), dimension(ib:ie,jb:je,nrain) :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, intent(in), dimension(ib:ie,jb:je) :: tsk,znt,ust,cd,ch,cq,u1,v1,w1,thflux,qvflux
      real, intent(in), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(in), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: rho,prs
      real, intent(in), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(in), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(in), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(in), dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(in), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(in), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,effc,effi,effs,effr,effg,effis
      real, intent(in), dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin
      real, intent(in), dimension(ni,nj,nrad2d) :: rad2d
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(in), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(in), dimension(ibl:iel,jbl:jel) :: psfc,u10,v10,w10,hfx,qfx,xland, &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,   &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      f2d,gsw,glw,chklowq,capg,snowc
      real, intent(in), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(in), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, intent(in), dimension(ibm:iem,jbm:jem,kbm:kem) :: qpten,qtten,qvten,qcten
      real, intent(in), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, intent(in), dimension(npvals,nparcels) :: pdata
      real, intent(inout), dimension(3,nparcels) :: ploc
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1
      real, intent(inout), dimension(jb:je,kb:ke) :: radwebc
      real, intent(inout), dimension(ib:ie,kb:ke) :: radsnbc
      real, intent(inout), dimension(ni+1,nj+1) :: dat1
      real, intent(inout), dimension((ni+1)*nodex,(nj+1)*nodey) :: dat2
      real, intent(inout), dimension(ni+1,nj+1,numprocs) :: dat3

      character*80 fname
      integer :: i,j,n,np,mynode,nodemaster,nodes
      real*8, dimension(nbudget,0:numprocs-1) :: sbudget
      real*8, dimension(numq,0:numprocs-1) :: csq,dsq
#ifdef MPI
      include 'mpif.h'
      integer :: proc,req1,req2,req3,index,count,reqs,reqt
      integer, dimension(0:numprocs-1) :: reqa,reqb,reqc
      integer status(MPI_STATUS_SIZE)
      logical :: doit
      integer, parameter :: nlim = 1000
#endif

!-----------------------------------------------------------------------

#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      mynode = int( float(myid)/float(ppnode) )
      nodemaster = mynode * ppnode
      nodes = nodex * nodey / ppnode
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  myid,mynode,nodemaster,nodes = ',myid,mynode,nodemaster,nodes
      if(dowr) write(outfile,*)
#else
      mynode = 0
      nodemaster = 0
      nodes = 1
#endif

    IF(myid.eq.nodemaster)THEN
      fname = '                                                                                '
    if(strlen.gt.0)then
      fname(1:strlen) = output_path(1:strlen)
    endif
      fname(strlen+1:strlen+baselen) = output_basename(1:baselen)
      fname(totlen+1:totlen+1+22) = '_rst_XXXXXX_YYYYYY.dat'

      write(fname(totlen+ 6:totlen+11),101) mynode
      write(fname(totlen+13:totlen+18),101) nrst
101   format(i6.6)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Writing to restart file!'
      if(dowr) write(outfile,*) '  fname=',fname
      if(dowr) write(outfile,*)

      open(unit=50,file=fname,form='unformatted',status='unknown')
    ENDIF

!---------------------------------------------------------------
! metadata:

    IF(myid.eq.0)THEN
      ! only processor 0 really need to do this:
      write(50) nstep
      write(50) nrec
      write(50) prec
      write(50) nwrite
      write(50) nrst
      write(50) dt
      write(50) mtime
      write(50) stattim
      write(50) taptim
      write(50) rsttim
      write(50) radtim
      write(50) prcltim
    ENDIF

!---------------------------------------------------------------
! budget variables:

#ifdef MPI
    IF(myid.eq.0)THEN
      do proc=1,numprocs-1
        call MPI_IRECV(sbudget(1,proc),nbudget,MPI_DOUBLE_PRECISION,proc,1,MPI_COMM_WORLD,reqa(proc),ierr)
        call MPI_IRECV(csq(1,proc)    ,numq   ,MPI_DOUBLE_PRECISION,proc,2,MPI_COMM_WORLD,reqb(proc),ierr)
        call MPI_IRECV(dsq(1,proc)    ,numq   ,MPI_DOUBLE_PRECISION,proc,3,MPI_COMM_WORLD,reqc(proc),ierr)
      enddo
#endif
      do i=1,nbudget
        sbudget(i,0) = qbudget(i)
      enddo
      do i=1,numq
        csq(i,0) = asq(i)
        dsq(i,0) = bsq(i)
      enddo
#ifdef MPI
      n = 1
      do while( n.le.(numprocs-1) )
        n = n + 1
        call mpi_waitany(numprocs-1,reqa(1),index,status,ierr)
      enddo
#endif
      !--------------------
      !cccccccccccccccccccc
      !--------------------
      write(50) sbudget
      !--------------------
      !cccccccccccccccccccc
      !--------------------
#ifdef MPI
      n = 1
      do while( n.le.(numprocs-1) )
        n = n + 1
        call mpi_waitany(numprocs-1,reqb(1),index,status,ierr)
      enddo
#endif
      !--------------------
      !cccccccccccccccccccc
      !--------------------
      write(50) csq
      !--------------------
      !cccccccccccccccccccc
      !--------------------
#ifdef MPI
      n = 1
      do while( n.le.(numprocs-1) )
        n = n + 1
        call mpi_waitany(numprocs-1,reqc(1),index,status,ierr)
      enddo
#endif
      !--------------------
      !cccccccccccccccccccc
      !--------------------
      write(50) dsq
      !--------------------
      !cccccccccccccccccccc
      !--------------------
#ifdef MPI
    ELSE
      call MPI_ISEND(qbudget(1),nbudget,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,req1,ierr)
      call MPI_ISEND(asq(1)    ,numq   ,MPI_DOUBLE_PRECISION,0,2,MPI_COMM_WORLD,req2,ierr)
      call MPI_ISEND(bsq(1)    ,numq   ,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,req3,ierr)
      call MPI_WAIT(req1,status,ierr)
      call MPI_WAIT(req2,status,ierr)
      call MPI_WAIT(req3,status,ierr)
    ENDIF
#endif

!---------------------------------------------------------------
#ifdef MPI
    !  limit to "nlim" writes at a time:
    IF( myid.eq.nodemaster .and. nodes.ge.nlim )THEN
      doit = .false.
      IF( mynode.ge.nlim )THEN
        call MPI_IRECV(doit,1,mpi_logical,myid-nlim*ppnode,999998,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,status,ierr)
      ENDIF
    ENDIF
#endif
!---------------------------------------------------------------
! standard 2D:

    do n=1,nrain
      call writer(ni,nj,1,1,nx,ny,rain(ib,jb,n),          &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,1,nx,ny,sws(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,1,nx,ny,svs(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,1,nx,ny,sps(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,1,nx,ny,srs(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,1,nx,ny,sgs(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,1,nx,ny,sus(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,1,nx,ny,shs(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    enddo
      call writer(ni,nj,1,1,nx,ny,tsk(ib,jb),             &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)

!---------------------------------------------------------------
! standard 3D:

      call writer(ni,nj,1,nk,nx,ny,rho(ib,jb,1),          &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,nk,nx,ny,prs(ib,jb,1),          &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni+1,nj,1,nk,nx+1,ny,ua(ib,jb,1),       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj+1,1,nk,nx,ny+1,va(ib,jb,1),       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,nk+1,nx,ny,wa(ib,jb,1),         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,nk,nx,ny,ppi(ib,jb,1),          &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,nk,nx,ny,tha(ib,jb,1),          &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    IF(imoist.eq.1)THEN
    do n=1,numq
      call writer(ni,nj,1,nk,nx,ny,qa(ib,jb,1,n),         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    enddo
    ENDIF
    if(imoist.eq.1.and.eqtset.eq.2)then
      call writer(ni,nj,1,nk,nx,ny,qpten(ib,jb,1),        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,nk,nx,ny,qtten(ib,jb,1),        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,nk,nx,ny,qvten(ib,jb,1),        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call writer(ni,nj,1,nk,nx,ny,qcten(ib,jb,1),        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    endif
    if(iturb.eq.1)then
      call writer(ni,nj,1,nk+1,nx,ny,tkea(ib,jb,1),       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    endif

!---------------------------------------------------------------
!  radiation:

      if(radopt.eq.1)then
        call writer(ni,nj,1,nk,nx,ny,lwten(ib,jb,1),        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,nk,nx,ny,swten(ib,jb,1),        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radsw(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = rnflx(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radswnet(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = radlwin(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        do n=1,nrad2d
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = rad2d(i,j,n)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        enddo
      endif
      if( radopt.ge.1 .and. ptype.eq.5 )then
        call writer(ni,nj,1,nk,nx,ny,effc(ib,jb,1),         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,nk,nx,ny,effi(ib,jb,1),         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,nk,nx,ny,effs(ib,jb,1),         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,nk,nx,ny,effr(ib,jb,1),         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,nk,nx,ny,effg(ib,jb,1),         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,nk,nx,ny,effis(ib,jb,1),        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

!---------------------------------------------------------------
!  surface:
!     I don't know how many of these are really needed in restart
!     files, but let's include them all for now ... just to be safe

      if((oceanmodel.eq.2).or.(ipbl.eq.1).or.(sfcmodel.ge.1))then
        !---- (1) ----!
      if(sfcmodel.ge.1)then
        call writer(ni,nj,1,1,nx,ny,ust(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,znt(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,cd(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,ch(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,cq(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,u1(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,v1(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,w1(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,u10(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,v10(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,w10(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,xland(ib,jb),           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,thflux(ib,jb),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,qvflux(ib,jb),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(sfcmodel.eq.2)then
        !---- (2) ----!
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = lu_index(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          dum1(i,j,1) = kpbl2d(i,j)
        enddo
        enddo
        call writer(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,psfc(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,hfx(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,qfx(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,hpbl(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,wspd(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,psim(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,psih(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,gz1oz0(ib,jb),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,br(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,CHS(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,CHS2(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,CQS2(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,CPMM(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,ZOL(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,MAVAIL(ib,jb),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,MOL(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,RMOL(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,REGIME(ib,jb),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,LH(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,tmn(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,FLHC(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,FLQC(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,QGH(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,CK(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,CKA(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,CDA(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,USTM(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,QSFC(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,T2(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,Q2(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,TH2(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,EMISS(ib,jb),           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,THC(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,ALBD(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,gsw(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,glw(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,chklowq(ib,jb),         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,capg(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,snowc(ib,jb),           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        do n=1,num_soil_layers
        call writer(ni,nj,1,1,nx,ny,tslb(ib,jb,n),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        enddo
      endif
      endif
      if(oceanmodel.eq.2)then
        call writer(ni,nj,1,1,nx,ny,tml(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,t0ml(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,hml(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,h0ml(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,huml(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,hvml(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call writer(ni,nj,1,1,nx,ny,tmoml(ib,jb),           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

!---------------------------------------------------------------
!  passive tracers and parcels:

      if(iptra.eq.1)then
        do n=1,npt
          call writer(ni,nj,1,nk,nx,ny,pta(ib,jb,1,n),        &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        enddo
      endif

      if(iprcl.eq.1)then
        ! only write position info:
        if(myid.eq.0)then
        DO np=1,nparcels
          ploc(1,np)=pdata(1,np)
          ploc(2,np)=pdata(2,np)
          ploc(3,np)=pdata(3,np)
        ENDDO
        write(50) ploc
        endif
      endif

!---------------------------------------------------------------
!  open bc:

      if(irbc.eq.4.and.wbc.eq.2)then
#ifndef MPI
        write(50) radbcw
#else
        count = (je-jb+1)*(ke-kb+1)
        IF(myid.ne.nodemaster)THEN
          call MPI_ISEND(radbcw,count,MPI_REAL,nodemaster,31,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ELSE
          write(50) radbcw
          do proc=myid+1,myid+(ppnode-1)
            call MPI_IRECV(radwebc,count,MPI_REAL,proc,31,MPI_COMM_WORLD,reqt,ierr)
            call MPI_WAIT(reqt,status,ierr)
            write(50) radwebc
          enddo
        ENDIF
#endif
      endif
      if(irbc.eq.4.and.ebc.eq.2)then
#ifndef MPI
        write(50) radbce
#else
        count = (je-jb+1)*(ke-kb+1)
        IF(myid.ne.nodemaster)THEN
          call MPI_ISEND(radbce,count,MPI_REAL,nodemaster,32,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ELSE
          write(50) radbce
          do proc=myid+1,myid+(ppnode-1)
            call MPI_IRECV(radwebc,count,MPI_REAL,proc,32,MPI_COMM_WORLD,reqt,ierr)
            call MPI_WAIT(reqt,status,ierr)
            write(50) radwebc
          enddo
        ENDIF
#endif
      endif
      if(irbc.eq.4.and.sbc.eq.2)then
#ifndef MPI
        write(50) radbcs
#else
        count = (ie-ib+1)*(ke-kb+1)
        IF(myid.ne.nodemaster)THEN
          call MPI_ISEND(radbcs,count,MPI_REAL,nodemaster,33,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ELSE
          write(50) radbcs
          do proc=myid+1,myid+(ppnode-1)
            call MPI_IRECV(radsnbc,count,MPI_REAL,proc,33,MPI_COMM_WORLD,reqt,ierr)
            call MPI_WAIT(reqt,status,ierr)
            write(50) radsnbc
          enddo
        ENDIF
#endif
      endif
      if(irbc.eq.4.and.nbc.eq.2)then
#ifndef MPI
        write(50) radbcn
#else
        count = (ie-ib+1)*(ke-kb+1)
        IF(myid.ne.nodemaster)THEN
          call MPI_ISEND(radbcn,count,MPI_REAL,nodemaster,34,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ELSE
          write(50) radbcn
          do proc=myid+1,myid+(ppnode-1)
            call MPI_IRECV(radsnbc,count,MPI_REAL,proc,34,MPI_COMM_WORLD,reqt,ierr)
            call MPI_WAIT(reqt,status,ierr)
            write(50) radsnbc
          enddo
        ENDIF
#endif
      endif

!---------------------------------------------------------------
#ifdef MPI
    !  limit to "nlim" writes at a time:
    IF( myid.eq.nodemaster .and. nodes.ge.nlim )THEN
      doit = .true.
      IF( mynode+nlim .le. (nodes-1) )THEN
        call MPI_ISEND(doit,1,mpi_logical,myid+nlim*ppnode,999998,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,status,ierr)
      ENDIF
    ENDIF
#endif
!---------------------------------------------------------------

      IF(myid.eq.nodemaster) close(unit=50)

#ifdef MPI
      if(timestats.ge.1)then
        ! this is needed for proper accounting of timing:
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      endif
#endif

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      subroutine read_restart(nstep,nrec,prec,nwrite,nrst,nrad2d,num_soil_layers,  &
                              dt,mtime,stattim,taptim,rsttim,radtim,prcltim,       &
                              qbudget,asq,bsq,                                     &
                              rain,sws,svs,sps,srs,sgs,sus,shs,                    &
                              tsk,znt,ust,cd,ch,cq,u1,v1,w1,thflux,qvflux,         &
                              radbcw,radbce,radbcs,radbcn,                         &
                              rho,prs,ua,va,wa,ppi,tha,qa,tkea,swten,lwten,        &
                              radsw,rnflx,radswnet,radlwin,rad2d,                  &
                              effc,effi,effs,effr,effg,effis,                      &
                              lu_index,kpbl2d,psfc,u10,v10,w10,hfx,qfx,xland,      &
                              hpbl,wspd,psim,psih,gz1oz0,br,                       &
                              CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,                       &
                              MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,                    &
                              CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,       &
                              f2d,gsw,glw,chklowq,capg,snowc,tslb,                 &
                              tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml,               &
                              qpten,qtten,qvten,qcten,pta,pdata,ploc,              &
                              dum1,radwebc,radsnbc,dat1,dat2,dat3)
      implicit none

      !----------------------------------------------------------
      ! This subroutine organizes the reading of restart files
      !----------------------------------------------------------

      include 'input.incl'
      include 'constants.incl'

      integer, intent(inout) :: nstep,nrec,prec,nwrite,nrst
      integer, intent(in) :: nrad2d,num_soil_layers
      real, intent(inout) :: dt
      double precision, intent(inout) :: mtime,stattim,taptim,rsttim,radtim,prcltim
      real*8, intent(inout), dimension(nbudget) :: qbudget
      real*8, intent(inout), dimension(numq) :: asq,bsq
      real, intent(inout), dimension(ib:ie,jb:je,nrain) :: rain,sws,svs,sps,srs,sgs,sus,shs
      real, intent(inout), dimension(ib:ie,jb:je) :: tsk,znt,ust,cd,ch,cq,u1,v1,w1,thflux,qvflux
      real, intent(inout), dimension(jb:je,kb:ke) :: radbcw,radbce
      real, intent(inout), dimension(ib:ie,kb:ke) :: radbcs,radbcn
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: rho,prs
      real, intent(inout), dimension(ib:ie+1,jb:je,kb:ke) :: ua
      real, intent(inout), dimension(ib:ie,jb:je+1,kb:ke) :: va
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke+1) :: wa
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: ppi,tha
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem,numq) :: qa
      real, intent(inout), dimension(ibt:iet,jbt:jet,kbt:ket) :: tkea
      real, intent(inout), dimension(ibr:ier,jbr:jer,kbr:ker) :: swten,lwten,effc,effi,effs,effr,effg,effis
      real, intent(inout), dimension(ni,nj) :: radsw,rnflx,radswnet,radlwin
      real, intent(inout), dimension(ni,nj,nrad2d) :: rad2d
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: lu_index
      integer, intent(inout), dimension(ibl:iel,jbl:jel) :: kpbl2d
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: psfc,u10,v10,w10,hfx,qfx,xland, &
                                      hpbl,wspd,psim,psih,gz1oz0,br,          &
                                      CHS,CHS2,CQS2,CPMM,ZOL,MAVAIL,          &
                                      MOL,RMOL,REGIME,LH,FLHC,FLQC,QGH,   &
                                      CK,CKA,CDA,USTM,QSFC,T2,Q2,TH2,EMISS,THC,ALBD,   &
                                      f2d,gsw,glw,chklowq,capg,snowc
      real, intent(inout), dimension(ibl:iel,jbl:jel,num_soil_layers) :: tslb
      real, intent(inout), dimension(ibl:iel,jbl:jel) :: tmn,tml,t0ml,hml,h0ml,huml,hvml,tmoml
      real, intent(inout), dimension(ibm:iem,jbm:jem,kbm:kem) :: qpten,qtten,qvten,qcten
      real, intent(inout), dimension(ibp:iep,jbp:jep,kbp:kep,npt) :: pta
      real, intent(inout), dimension(npvals,nparcels) :: pdata
      real, intent(inout), dimension(3,nparcels) :: ploc
      real, intent(inout), dimension(ib:ie,jb:je,kb:ke) :: dum1
      real, intent(inout), dimension(jb:je,kb:ke) :: radwebc
      real, intent(inout), dimension(ib:ie,kb:ke) :: radsnbc
      real, intent(inout), dimension(ni+1,nj+1) :: dat1
      real, intent(inout), dimension((ni+1)*nodex,(nj+1)*nodey) :: dat2
      real, intent(inout), dimension(ni+1,nj+1,numprocs) :: dat3

      character*80 fname
      integer :: i,j,n,np,mynode,nodemaster,nodes
      real*8, dimension(nbudget,0:numprocs-1) :: sbudget
      real*8, dimension(numq,0:numprocs-1) :: csq,dsq
#ifdef MPI
      include 'mpif.h'
      integer :: proc,req1,req2,req3,index,count,reqs,reqt
      integer, dimension(0:numprocs-1) :: reqa,reqb,reqc
      integer status(MPI_STATUS_SIZE)
      logical :: doit
      integer, parameter :: nlim = 1000
#endif

!-----------------------------------------------------------------------

#ifdef MPI
      call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      mynode = int( float(myid)/float(ppnode) )
      nodemaster = mynode * ppnode
      nodes = nodex * nodey / ppnode
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  myid,mynode,nodemaster,nodes = ',myid,mynode,nodemaster,nodes
      if(dowr) write(outfile,*)
#else
      mynode = 0
      nodemaster = 0
      nodes = 1
#endif

    IF(myid.eq.nodemaster)THEN
      fname = '                                                                                '
    if(strlen.gt.0)then
      fname(1:strlen) = output_path(1:strlen)
    endif
      fname(strlen+1:strlen+baselen) = output_basename(1:baselen)
      fname(totlen+1:totlen+1+22) = '_rst_XXXXXX_YYYYYY.dat'

      write(fname(totlen+ 6:totlen+11),101) mynode
      write(fname(totlen+13:totlen+18),101) rstnum
101   format(i6.6)

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  Reading from restart file!'
      if(dowr) write(outfile,*) '  fname=',fname
      if(dowr) write(outfile,*)

      open(unit=50,file=fname,form='unformatted',status='old')
    ENDIF

!---------------------------------------------------------------
! metadata:

    IF(myid.eq.0)THEN
      ! only processor 0 has these variables:
      read(50) nstep
      read(50) nrec
      read(50) prec
      read(50) nwrite
      read(50) nrst
      read(50) dt
      read(50) mtime
      read(50) stattim
      read(50) taptim
      read(50) rsttim
      read(50) radtim
      read(50) prcltim
    ENDIF

#ifdef MPI
      ! communicate to all other processors:
      call MPI_BCAST(nstep  ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nrec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prec   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nwrite ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nrst   ,1,MPI_INTEGER         ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dt     ,1,MPI_REAL            ,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mtime  ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(stattim,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(taptim ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rsttim ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(radtim ,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prcltim,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#endif

!---------------------------------------------------------------
! budget variables:

#ifdef MPI
    IF(myid.eq.0)THEN
#endif
      !--------------------
      !cccccccccccccccccccc
      !--------------------
      read(50) sbudget
      read(50) csq
      read(50) dsq
      !--------------------
      !cccccccccccccccccccc
      !--------------------
#ifdef MPI
      do proc=1,numprocs-1
        call MPI_ISEND(sbudget(1,proc),nbudget,MPI_DOUBLE_PRECISION,proc,1,MPI_COMM_WORLD,reqa(proc),ierr)
        call MPI_ISEND(csq(1,proc)    ,numq   ,MPI_DOUBLE_PRECISION,proc,2,MPI_COMM_WORLD,reqb(proc),ierr)
        call MPI_ISEND(dsq(1,proc)    ,numq   ,MPI_DOUBLE_PRECISION,proc,3,MPI_COMM_WORLD,reqc(proc),ierr)
      enddo
#endif
      do i=1,nbudget
        qbudget(i) = sbudget(i,0)
      enddo
      do i=1,numq
        asq(i) = csq(i,0)
        bsq(i) = dsq(i,0)
      enddo
#ifdef MPI
      n = 1
      do while( n.le.(numprocs-1) )
        n = n + 1
        call mpi_waitany(numprocs-1,reqa(1),index,status,ierr)
      enddo
      n = 1
      do while( n.le.(numprocs-1) )
        n = n + 1
        call mpi_waitany(numprocs-1,reqb(1),index,status,ierr)
      enddo
      n = 1
      do while( n.le.(numprocs-1) )
        n = n + 1
        call mpi_waitany(numprocs-1,reqc(1),index,status,ierr)
      enddo
    ELSE
      call MPI_IRECV(qbudget(1),nbudget,MPI_DOUBLE_PRECISION,0,1,MPI_COMM_WORLD,req1,ierr)
      call MPI_IRECV(asq(1)    ,numq   ,MPI_DOUBLE_PRECISION,0,2,MPI_COMM_WORLD,req2,ierr)
      call MPI_IRECV(bsq(1)    ,numq   ,MPI_DOUBLE_PRECISION,0,3,MPI_COMM_WORLD,req3,ierr)
      call MPI_WAIT(req1,status,ierr)
      call MPI_WAIT(req2,status,ierr)
      call MPI_WAIT(req3,status,ierr)
    ENDIF
#endif

!---------------------------------------------------------------
#ifdef MPI
    !  limit to "nlim" reads at a time:
    IF( myid.eq.nodemaster .and. nodes.ge.nlim )THEN
      doit = .false.
      IF( mynode.ge.nlim )THEN
        call MPI_IRECV(doit,1,mpi_logical,myid-nlim*ppnode,999997,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,status,ierr)
      ENDIF
    ENDIF
#endif
!---------------------------------------------------------------
! standard 2D:

    do n=1,nrain
      call  readr(ni,nj,1,1,nx,ny,rain(ib,jb,n),          &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,1,nx,ny,sws(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,1,nx,ny,svs(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,1,nx,ny,sps(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,1,nx,ny,srs(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,1,nx,ny,sgs(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,1,nx,ny,sus(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,1,nx,ny,shs(ib,jb,n),           &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    enddo
      call  readr(ni,nj,1,1,nx,ny,tsk(ib,jb),             &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)

!---------------------------------------------------------------
! standard 3D:

      call  readr(ni,nj,1,nk,nx,ny,rho(ib,jb,1),          &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,nk,nx,ny,prs(ib,jb,1),          &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni+1,nj,1,nk,nx+1,ny,ua(ib,jb,1),       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj+1,1,nk,nx,ny+1,va(ib,jb,1),       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,nk+1,nx,ny,wa(ib,jb,1),         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,nk,nx,ny,ppi(ib,jb,1),          &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,nk,nx,ny,tha(ib,jb,1),          &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    IF(imoist.eq.1)THEN
    do n=1,numq
      call  readr(ni,nj,1,nk,nx,ny,qa(ib,jb,1,n),         &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    enddo
    ENDIF
    if(imoist.eq.1.and.eqtset.eq.2)then
      call  readr(ni,nj,1,nk,nx,ny,qpten(ib,jb,1),        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,nk,nx,ny,qtten(ib,jb,1),        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,nk,nx,ny,qvten(ib,jb,1),        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      call  readr(ni,nj,1,nk,nx,ny,qcten(ib,jb,1),        &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    endif
    if(iturb.eq.1)then
      call  readr(ni,nj,1,nk+1,nx,ny,tkea(ib,jb,1),       &
                  ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                  dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
    endif

!---------------------------------------------------------------
!  radiation:

      if(radopt.eq.1)then
        call  readr(ni,nj,1,nk,nx,ny,lwten(ib,jb,1),        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,nk,nx,ny,swten(ib,jb,1),        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          radsw(i,j) = dum1(i,j,1)
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          rnflx(i,j) = dum1(i,j,1)
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          radswnet(i,j) = dum1(i,j,1)
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          radlwin(i,j) = dum1(i,j,1)
        enddo
        enddo
        do n=1,nrad2d
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          rad2d(i,j,n) = dum1(i,j,1)
        enddo
        enddo
        enddo
      endif
      if( radopt.ge.1 .and. ptype.eq.5 )then
        call  readr(ni,nj,1,nk,nx,ny,effc(ib,jb,1),         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,nk,nx,ny,effi(ib,jb,1),         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,nk,nx,ny,effs(ib,jb,1),         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,nk,nx,ny,effr(ib,jb,1),         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,nk,nx,ny,effg(ib,jb,1),         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,nk,nx,ny,effis(ib,jb,1),        &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

!---------------------------------------------------------------
!  surface:
!     I don't know how many of these are really needed in restart
!     files, but let's include them all for now ... just to be safe

      if((oceanmodel.eq.2).or.(ipbl.eq.1).or.(sfcmodel.ge.1))then
        !---- (1) ----!
      if(sfcmodel.ge.1)then
        call  readr(ni,nj,1,1,nx,ny,ust(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,znt(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,cd(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,ch(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,cq(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,u1(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,v1(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,w1(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,u10(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,v10(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,w10(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,xland(ib,jb),           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,thflux(ib,jb),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,qvflux(ib,jb),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif
      if(sfcmodel.eq.2)then
        !---- (2) ----!
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          lu_index(i,j) = nint(dum1(i,j,1))
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,dum1(ib,jb,1),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
!$omp parallel do default(shared)  &
!$omp private(i,j)
        do j=1,nj
        do i=1,ni
          kpbl2d(i,j) = nint(dum1(i,j,1))
        enddo
        enddo
        call  readr(ni,nj,1,1,nx,ny,psfc(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,hfx(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,qfx(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,hpbl(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,wspd(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,psim(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,psih(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,gz1oz0(ib,jb),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,br(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,CHS(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,CHS2(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,CQS2(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,CPMM(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,ZOL(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,MAVAIL(ib,jb),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,MOL(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,RMOL(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,REGIME(ib,jb),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,LH(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,tmn(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,FLHC(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,FLQC(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,QGH(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,CK(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,CKA(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,CDA(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,USTM(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,QSFC(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,T2(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,Q2(ib,jb),              &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,TH2(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,EMISS(ib,jb),           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,THC(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,ALBD(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,gsw(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,glw(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,chklowq(ib,jb),         &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,capg(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,snowc(ib,jb),           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        do n=1,num_soil_layers
        call  readr(ni,nj,1,1,nx,ny,tslb(ib,jb,n),          &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        enddo
      endif
      endif
      if(oceanmodel.eq.2)then
        call  readr(ni,nj,1,1,nx,ny,tml(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,t0ml(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,hml(ib,jb),             &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,h0ml(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,huml(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,hvml(ib,jb),            &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        call  readr(ni,nj,1,1,nx,ny,tmoml(ib,jb),           &
                    ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                    dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
      endif

!---------------------------------------------------------------
!  passive tracers and parcels:

      if(iptra.eq.1)then
        do n=1,npt
          call  readr(ni,nj,1,nk,nx,ny,pta(ib,jb,1,n),        &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                      dat1(1,1),dat2(1,1),dat3(1,1,1),ppnode)
        enddo
      endif

      if(iprcl.eq.1)then
        ! only read position info:
        if(myid.eq.0)then
        read(50) ploc
        endif
#ifdef MPI
        call MPI_BCAST(ploc,3*nparcels,MPI_REAL,0,MPI_COMM_WORLD,ierr)
#endif
        DO np=1,nparcels
          pdata(1,np)=ploc(1,np)
          pdata(2,np)=ploc(2,np)
          pdata(3,np)=ploc(3,np)
        ENDDO
      endif

!---------------------------------------------------------------
!  open bc:

      if(irbc.eq.4.and.wbc.eq.2)then
#ifndef MPI
        read(50) radbcw
#else
        count = (je-jb+1)*(ke-kb+1)
        IF(myid.ne.nodemaster)THEN
          call MPI_IRECV(radbcw,count,MPI_REAL,nodemaster,31,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ELSE
          read(50) radbcw
          do proc=myid+1,myid+(ppnode-1)
            read(50) radwebc
            call MPI_ISEND(radwebc,count,MPI_REAL,proc,31,MPI_COMM_WORLD,reqt,ierr)
            call MPI_WAIT(reqt,status,ierr)
          enddo
        ENDIF
#endif
      endif
      if(irbc.eq.4.and.ebc.eq.2)then
#ifndef MPI
        read(50) radbce
#else
        count = (je-jb+1)*(ke-kb+1)
        IF(myid.ne.nodemaster)THEN
          call MPI_IRECV(radbce,count,MPI_REAL,nodemaster,32,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ELSE
          read(50) radbce
          do proc=myid+1,myid+(ppnode-1)
            read(50) radwebc
            call MPI_ISEND(radwebc,count,MPI_REAL,proc,32,MPI_COMM_WORLD,reqt,ierr)
            call MPI_WAIT(reqt,status,ierr)
          enddo
        ENDIF
#endif
      endif
      if(irbc.eq.4.and.sbc.eq.2)then
#ifndef MPI
        read(50) radbcs
#else
        count = (ie-ib+1)*(ke-kb+1)
        IF(myid.ne.nodemaster)THEN
          call MPI_IRECV(radbcs,count,MPI_REAL,nodemaster,33,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ELSE
          read(50) radbcs
          do proc=myid+1,myid+(ppnode-1)
            read(50) radsnbc
            call MPI_ISEND(radsnbc,count,MPI_REAL,proc,33,MPI_COMM_WORLD,reqt,ierr)
            call MPI_WAIT(reqt,status,ierr)
          enddo
        ENDIF
#endif
      endif
      if(irbc.eq.4.and.nbc.eq.2)then
#ifndef MPI
        read(50) radbcn
#else
        count = (ie-ib+1)*(ke-kb+1)
        IF(myid.ne.nodemaster)THEN
          call MPI_IRECV(radbcn,count,MPI_REAL,nodemaster,34,MPI_COMM_WORLD,reqs,ierr)
          call MPI_WAIT(reqs,status,ierr)
        ELSE
          read(50) radbcn
          do proc=myid+1,myid+(ppnode-1)
            read(50) radsnbc
            call MPI_ISEND(radsnbc,count,MPI_REAL,proc,34,MPI_COMM_WORLD,reqt,ierr)
            call MPI_WAIT(reqt,status,ierr)
          enddo
        ENDIF
#endif
      endif

!---------------------------------------------------------------
#ifdef MPI
    !  limit to "nlim" reads at a time:
    IF( myid.eq.nodemaster .and. nodes.ge.nlim )THEN
      doit = .true.
      IF( mynode+nlim .le. (nodes-1) )THEN
        call MPI_ISEND(doit,1,mpi_logical,myid+nlim*ppnode,999997,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,status,ierr)
      ENDIF
    ENDIF
#endif
!---------------------------------------------------------------

      IF(myid.eq.nodemaster) close(unit=50)

!---------

      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '  From restart file: '
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   mtime   = ',mtime
      if(dowr) write(outfile,*) '   stattim = ',stattim
      if(dowr) write(outfile,*) '   taptim  = ',taptim
      if(dowr) write(outfile,*) '   rsttim  = ',rsttim
      if(dowr) write(outfile,*) '   radtim  = ',radtim
      if(dowr) write(outfile,*) '   prcltim = ',prcltim
      if(dowr) write(outfile,*)
      if(dowr) write(outfile,*) '   nstep   = ',nstep
      if(dowr) write(outfile,*) '   nrec    = ',nrec
      if(dowr) write(outfile,*) '   prec    = ',prec
      if(dowr) write(outfile,*) '   nwrite  = ',nwrite
      if(dowr) write(outfile,*) '   nrst    = ',nrst
      if(dowr) write(outfile,*)

!---------

      if( adapt_dt.eq.0 ) dt = dtl

      ! this is needed for stats files:
      if( output_format.eq.2 )then
        nrec=nrec-1
      else
        nrec=nrec-stat_out
      endif

!---------

#ifdef MPI
      if(timestats.ge.1)then
        ! this is needed for proper accounting of timing:
        call MPI_BARRIER (MPI_COMM_WORLD,ierr)
      endif
#endif

      return
      end


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine writer(numi,numj,numk1,numk2,nxr,nyr,var,      &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                      dat1,dat2,dat3,ppnode)
    implicit none

    !-------------------------------------------------------------------
    ! This subroutine collects data (from other processors if this is a
    ! MPI run) and does the actual writing of restart files.
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr
    integer, intent(in) :: ppnode
    real, intent(in), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    integer, intent(in) :: ni,nj,ngxy,myid,numprocs,nodex,nodey
    real, intent(inout), dimension(numi,numj) :: dat1
    real, intent(inout), dimension(numi*ppnode,numj) :: dat2
    real, intent(inout), dimension(numi,numj,0:numprocs-1) :: dat3

    integer :: i,j,k
    integer :: mynode,nodemaster
#ifdef MPI
    include 'mpif.h'
    integer :: n,nn,nnn,index,fooi,fooj,proc,ierr,reqs
    integer, dimension(numprocs-1) :: reqt
    integer status(MPI_STATUS_SIZE)
#endif

    mynode = int( float(myid)/float(ppnode) )
    nodemaster = mynode * ppnode

#ifndef MPI
    write(50) var
#else
    DO k=numk1,numk2
      IF(myid.ne.nodemaster)THEN
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat1(i,j) = var(i,j,k)
        enddo
        enddo
        call MPI_ISEND(dat1,numi*numj,MPI_REAL,nodemaster,k,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,status,ierr)
      ELSE
        do proc=myid+1,myid+(ppnode-1)
          call MPI_IRECV(dat3(1,1,proc),numi*numj,MPI_REAL,proc,k,MPI_COMM_WORLD,reqt(proc),ierr)
        enddo
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,numj
        do i=1,numi
          dat2(i,j)=var(i,j,k)
        enddo
        enddo
        nn = 1
        do while( nn.le.(ppnode-1) )
          nn = nn + 1
          call mpi_waitany(ppnode-1,reqt(myid+1),index,status,ierr)
          fooi = numi*index
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,numj
          do i=1,numi
            dat2(fooi+i,j)=dat3(i,j,nodemaster+index)
          enddo
          enddo
        enddo
        write(50) dat2
      ENDIF
    ENDDO
#endif

    return
    end subroutine writer


!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc


    subroutine  readr(numi,numj,numk1,numk2,nxr,nyr,var,      &
                      ni,nj,ngxy,myid,numprocs,nodex,nodey,   &
                      dat1,dat2,dat3,ppnode)
    implicit none

    !-------------------------------------------------------------------
    ! This subroutine reads restart files and then passes data 
    ! to other processors if this is a MPI run. 
    !-------------------------------------------------------------------

    integer, intent(in) :: numi,numj,numk1,numk2,nxr,nyr
    integer, intent(in) :: ppnode
    real, intent(inout), dimension(1-ngxy:numi+ngxy,1-ngxy:numj+ngxy,numk1:numk2) :: var
    integer, intent(in) :: ni,nj,ngxy,myid,numprocs,nodex,nodey
    real, intent(inout), dimension(numi,numj) :: dat1
    real, intent(inout), dimension(numi*ppnode,numj) :: dat2
    real, intent(inout), dimension(numi,numj,0:numprocs-1) :: dat3

    integer :: i,j,k
    integer :: mynode,nodemaster
#ifdef MPI
    include 'mpif.h'
    integer :: n,nn,nnn,index,fooi,fooj,proc,ierr,reqs
    integer, dimension(numprocs-1) :: reqt
    integer status(MPI_STATUS_SIZE)
#endif

    mynode = int( float(myid)/float(ppnode) )
    nodemaster = mynode * ppnode

#ifndef MPI
    read(50) var
#else
    DO k=numk1,numk2
      IF(myid.ne.nodemaster)THEN
        call MPI_IRECV(dat1,numi*numj,MPI_REAL,nodemaster,k,MPI_COMM_WORLD,reqs,ierr)
        call MPI_WAIT(reqs,status,ierr)
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,numj
        do i=1,numi
          var(i,j,k) = dat1(i,j)
        enddo
        enddo
      ELSE
        read(50) dat2
        do proc=myid+1,myid+(ppnode-1)
          fooi = numi*(proc-myid)
!$omp parallel do default(shared)   &
!$omp private(i,j)
          do j=1,numj
          do i=1,numi
            dat3(i,j,proc)=dat2(fooi+i,j)
          enddo
          enddo
          call MPI_ISEND(dat3(1,1,proc),numi*numj,MPI_REAL,proc,k,MPI_COMM_WORLD,reqt(proc),ierr)
        enddo
!$omp parallel do default(shared)   &
!$omp private(i,j)
        do j=1,numj
        do i=1,numi
          var(i,j,k)=dat2(i,j)
        enddo
        enddo
        nn = 1
        do while( nn.le.(ppnode-1) )
          nn = nn + 1
          call mpi_waitany(ppnode-1,reqt(myid+1),index,status,ierr)
        enddo
      ENDIF
    ENDDO
#endif

    return
    end subroutine  readr

