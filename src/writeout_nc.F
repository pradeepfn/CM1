
!-------------------------------------------------------------
!
!  This subroutine writes data in NetCDF files.
!
!  Code originally written by Daniel Kirshbaum
!  Code converted to netcdf4 (f90) by George Bryan, May 2013
!  Code last modified by George Bryan, 130910
!
!-------------------------------------------------------------


      subroutine netcdf_prelim(nwrite,ncid,time_index,qname,xh,xf,yh,yf, &
                               xfref,yfref,sigma,sigmaf,zs,zh,zf,        &
                               d2d,ds,du,dv,                             &
                               dumz,dumx,dumy)
#ifdef NETCDF
      use netcdf
#endif
      implicit none
      include 'input.incl'
      include 'constants.incl'

      integer, intent(in) :: nwrite
      integer, intent(inout) :: ncid,time_index
      character*3, dimension(maxq), intent(in) :: qname
      real, dimension(ib:ie),   intent(in) :: xh
      real, dimension(ib:ie+1), intent(in) :: xf
      real, dimension(jb:je),   intent(in) :: yh
      real, dimension(jb:je+1), intent(in) :: yf
      real, intent(in), dimension(-2:nx+4) :: xfref
      real, intent(in), dimension(-2:ny+4) :: yfref
      real, dimension(kb:ke)  , intent(in) :: sigma
      real, dimension(kb:ke+1), intent(in) :: sigmaf
      real, dimension(itb:ite,jtb:jte), intent(in) :: zs
      real, dimension(ib:ie,jb:je,kb:ke),   intent(in) :: zh
      real, dimension(ib:ie,jb:je,kb:ke+1), intent(in) :: zf
      real, dimension(ni,nj) :: d2d
      real, dimension(ni,nj,nk) :: ds
      real, dimension(ni+1,nj,nk) :: du
      real, dimension(ni,nj+1,nk) :: dv
      real, intent(inout), dimension(nz+1) :: dumz
      real, intent(inout), dimension(nx+1) :: dumx
      real, intent(inout), dimension(ny+1) :: dumy

#ifdef NETCDF
      integer i,j,k,n,irec

      ! Users of GrADS might want to set coards to .true.
      logical, parameter :: coards = .false.

      integer :: cdfid    ! ID for the netCDF file to be created
      integer :: status,dimid,varid
      integer :: niid,njid,nkid,nip1id,njp1id,nkp1id,timeid,oneid
      character*8 chid

!-------------------------------------------------------------
! Declare and set integer values for the netCDF dimensions 
!-------------------------------------------------------------

      integer :: ival,jval,ivalp1,jvalp1,kvalp1
      real :: actual_time

      logical :: allinfo

      integer, parameter :: shuffle       = 1
      integer, parameter :: deflate       = 1
      integer, parameter :: deflate_level = 1

      integer, parameter :: cache_size       = 32000000
      integer, parameter :: cache_nelems     = 37
      integer, parameter :: cache_preemption = 100

      integer :: chkx,chky

      chkx = ni
      chky = nj

!--------------------------------------------------------------
! Initializing some things
!--------------------------------------------------------------

      time_index = nwrite
      if(output_filetype.ge.2) time_index = 1

    if(coards)then
      if(tapfrq.lt.60.0)then
        print *
        print *,'  Output frequency cannot be less than 60 s for coards format'
        print *
        call stopcm1
      endif
      actual_time = (nwrite-1)*tapfrq/60.0
    else
      actual_time = (nwrite-1)*tapfrq
    endif

!--------------------------------------------------------------
!  Write data to cdf file
!--------------------------------------------------------------

    IF(output_filetype.eq.1)THEN
      string(totlen+1:totlen+22) = '.nc                   '
    ELSEIF(output_filetype.eq.2)THEN
      string(totlen+1:totlen+22) = '_XXXXXX.nc            '
      write(string(totlen+2:totlen+7),100) nwrite
    ELSEIF(output_filetype.eq.3)THEN
      string(totlen+1:totlen+17) = '_XXXXXX_YYYYYY.nc     '
      write(string(totlen+ 2:totlen+ 7),100) myid
      write(string(totlen+ 9:totlen+14),100) nwrite
    ELSE
      if(dowr) write(outfile,*) '  for netcdf output, output_filetype must be either 1,2, or 3 '
      call stopcm1
    ENDIF

      if(myid.eq.0) print *,'  string = ',string

100   format(i6.6)

!--------------------------------------------------------------
!  Dimensions of data:

    IF( output_filetype.eq.1 .or. output_filetype.eq.2 )THEN
      ival = nx
      jval = ny
    ELSEIF( output_filetype.eq.3 )THEN
      ival = ni
      jval = nj
    ELSE
      print *,'  unrecognized value for output_filetype '
      call stopcm1
    ENDIF

    ivalp1 = ival+1
    jvalp1 = jval+1
    kvalp1 = nk+1

!--------------------------------------------------------------
!  if this is the start of a file, then do this stuff:

    allinfo = .false.
    IF(nwrite.eq.1) allinfo=.true.
    IF(output_filetype.ge.2) allinfo=.true.

    ifallinfo: IF(allinfo)THEN
!!!      print *,'  allinfo ... '

!-----------------------------------------------------------------------
!  BEGIN NEW:

#ifdef NCFPLUS
!--- works with netcdf 4.2, but not 4.0 (grumble)
      call disp_err( nf90_create(string,IOR(nf90_netcdf4, nf90_classic_model),ncid,   &
         cache_size=cache_size,cache_nelems=cache_nelems,cache_preemption=cache_preemption)   &
                     , .true. )
#else
      ! must do this for netcdf 4.0 (it seems) !
      call disp_err( nf90_create(path=string,cmode=0,ncid=ncid) )
#endif

      status = nf90_def_dim(ncid,'ni',ival,niid)
      status = nf90_def_dim(ncid,'nj',jval,njid)
      status = nf90_def_dim(ncid,'nk',nk,nkid)
      status = nf90_def_dim(ncid,'nip1',ivalp1,nip1id)
      status = nf90_def_dim(ncid,'njp1',jvalp1,njp1id)
      status = nf90_def_dim(ncid,'nkp1',kvalp1,nkp1id)
      status = nf90_def_dim(ncid,'time',nf90_unlimited,timeid)
      status = nf90_def_dim(ncid,'one',1,oneid)

    IF(icor.eq.1)THEN
      status = nf90_def_var(ncid,"f_cor",nf90_float,oneid,varid)
      status = nf90_put_att(ncid,varid,"def","Coriolis parameter")
      status = nf90_put_att(ncid,varid,"units","1/s")
    ENDIF

    if(.not.coards)then
      status = nf90_def_var(ncid,"ztop",nf90_float,oneid,varid)
      status = nf90_put_att(ncid,varid,"units","km")
    endif

    IF(coards)THEN

      status = nf90_def_var(ncid,"time",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"def","time since beginning of simulation")
      status = nf90_put_att(ncid,varid,"units","minutes since 2000-07-03 00:00:00")

      status = nf90_def_var(ncid,"ni",nf90_float,niid,varid)
      status = nf90_put_att(ncid,varid,"def","west-east location of scalar grid points")
      status = nf90_put_att(ncid,varid,"units","degree_east")

      status = nf90_def_var(ncid,"nip1",nf90_float,nip1id,varid)
      status = nf90_put_att(ncid,varid,"def","west-east location of staggered u grid points")
      status = nf90_put_att(ncid,varid,"units","degree_east")

      status = nf90_def_var(ncid,"nj",nf90_float,njid,varid)
      status = nf90_put_att(ncid,varid,"def","south-north location of scalar grid points")
      status = nf90_put_att(ncid,varid,"units","degree_north")

      status = nf90_def_var(ncid,"njp1",nf90_float,njp1id,varid)
      status = nf90_put_att(ncid,varid,"def","south-north location of staggered v grid points")
      status = nf90_put_att(ncid,varid,"units","degree_north")

      status = nf90_def_var(ncid,"nk",nf90_float,nkid,varid)
      status = nf90_put_att(ncid,varid,"def","nominal height of scalar grid points")
      status = nf90_put_att(ncid,varid,"units","km")

      status = nf90_def_var(ncid,"nkp1",nf90_float,nkp1id,varid)
      status = nf90_put_att(ncid,varid,"def","nominal height of staggered w grid points")
      status = nf90_put_att(ncid,varid,"units","km")

    ELSE

      status = nf90_def_var(ncid,"time",nf90_float,timeid,varid)
      status = nf90_put_att(ncid,varid,"def","time since beginning of simulation")
      status = nf90_put_att(ncid,varid,"units","seconds since 2000-07-03 00:00:00")

      status = nf90_def_var(ncid,"xh",nf90_float,niid,varid)
      status = nf90_put_att(ncid,varid,"def","west-east location of scalar grid points")
      status = nf90_put_att(ncid,varid,"units","km")

      status = nf90_def_var(ncid,"xf",nf90_float,nip1id,varid)
      status = nf90_put_att(ncid,varid,"def","west-east location of staggered u grid points")
      status = nf90_put_att(ncid,varid,"units","km")

      status = nf90_def_var(ncid,"yh",nf90_float,njid,varid)
      status = nf90_put_att(ncid,varid,"def","south-north location of scalar grid points")
      status = nf90_put_att(ncid,varid,"units","km")

      status = nf90_def_var(ncid,"yf",nf90_float,njp1id,varid)
      status = nf90_put_att(ncid,varid,"def","south-north location of staggered v grid points")
      status = nf90_put_att(ncid,varid,"units","km")

      status = nf90_def_var(ncid,"z",nf90_float,nkid,varid)
      status = nf90_put_att(ncid,varid,"def","nominal height of scalar grid points")
      status = nf90_put_att(ncid,varid,"units","km")

      status = nf90_def_var(ncid,"zf",nf90_float,nkp1id,varid)
      status = nf90_put_att(ncid,varid,"def","nominal height of staggered w grid points")
      status = nf90_put_att(ncid,varid,"units","km")

    ENDIF

!--------------------------------------------------------
!  Just to be sure:

        status = nf90_inq_dimid(ncid,'time',timeid)
        status = nf90_inq_dimid(ncid,'ni',niid)
        status = nf90_inq_dimid(ncid,'nj',njid)
        status = nf90_inq_dimid(ncid,'nk',nkid)
        status = nf90_inq_dimid(ncid,'nip1',nip1id)
        status = nf90_inq_dimid(ncid,'njp1',njp1id)
        status = nf90_inq_dimid(ncid,'nkp1',nkp1id)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--- 2D vars:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if(output_rain.eq.1)then
        status = nf90_def_var(ncid,"rain",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","accumulated surface rainfall")
        status = nf90_put_att(ncid,varid,"units","cm")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_sws.eq.1) then
        status = nf90_def_var(ncid,"sws",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","max windspeed at lowest level")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_svs.eq.1) then
        status = nf90_def_var(ncid,"svs",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","max vert vorticity at lowest level")
        status = nf90_put_att(ncid,varid,"units","s-1")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_sps.eq.1) then
        status = nf90_def_var(ncid,"sps",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","min pressure at lowest level")
        status = nf90_put_att(ncid,varid,"units","Pa")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_srs.eq.1) then
        status = nf90_def_var(ncid,"srs",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","max surface rainwater")
        status = nf90_put_att(ncid,varid,"units","kg/kg")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_sgs.eq.1) then
        status = nf90_def_var(ncid,"sgs",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","max surface graupel/hail")
        status = nf90_put_att(ncid,varid,"units","kg/kg")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_sus.eq.1) then
        status = nf90_def_var(ncid,"sus",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","max w at 5 km AGL")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_shs.eq.1) then
        status = nf90_def_var(ncid,"shs",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","max integrated updraft helicity")
        status = nf90_put_att(ncid,varid,"units","m2/s2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      IF(nrain.eq.2)THEN
        if(output_rain.eq.1)then
          status = nf90_def_var(ncid,"rain2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"def","accumulated surface rainfall, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","cm")
#ifdef NCFPLUS
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
        endif
        if(output_sws.eq.1) then
          status = nf90_def_var(ncid,"sws2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"def","max windspeed at lowest level, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
        endif
        if(output_svs.eq.1) then
          status = nf90_def_var(ncid,"svs2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"def","max vorticity at lowest level, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","s-1")
#ifdef NCFPLUS
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
        endif
        if(output_sps.eq.1) then
          status = nf90_def_var(ncid,"sps2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"def","min pressure at lowest level, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","Pa")
#ifdef NCFPLUS
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
        endif
        if(output_srs.eq.1) then
          status = nf90_def_var(ncid,"srs2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"def","max surface rainwater, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","kg/kg")
#ifdef NCFPLUS
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
        endif
        if(output_sgs.eq.1) then
          status = nf90_def_var(ncid,"sgs2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"def","max surface graupel/hail, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","kg/kg")
#ifdef NCFPLUS
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
        endif
        if(output_sus.eq.1) then
          status = nf90_def_var(ncid,"sus2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"def","max w at 5 km AGL, translated with moving domain")
          status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
        endif
        if(output_shs.eq.1) then
          status = nf90_def_var(ncid,"shs2",nf90_float,(/niid,njid,timeid/),varid)
          status = nf90_put_att(ncid,varid,"def","translated max integrated updraft helicity")
          status = nf90_put_att(ncid,varid,"units","m2/s2")
#ifdef NCFPLUS
          call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
        endif
      ENDIF

      IF(output_uh.eq.1)THEN
        status = nf90_def_var(ncid,"uh",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","integrated updraft helicity")
        status = nf90_put_att(ncid,varid,"units","m2/s2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      ENDIF

      IF (output_coldpool.eq.1) THEN
        status = nf90_def_var(ncid,"cpc",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","cold pool intensity C")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
        status = nf90_def_var(ncid,"cph",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","cold pool depth h")
        status = nf90_put_att(ncid,varid,"units","m AGL")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      ENDIF

      if (output_sfcflx.eq.1) then
        status = nf90_def_var(ncid,"thflux",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","surface potential temperature flux")
        status = nf90_put_att(ncid,varid,"units","K m s^{-1}")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"qvflux",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","surface water vapor flux")
        status = nf90_put_att(ncid,varid,"units","kg kg^{-1} m s^{-1}")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"cd",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","cd")
        status = nf90_put_att(ncid,varid,"units","nondimensional")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"ch",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","ch")
        status = nf90_put_att(ncid,varid,"units","nondimensional")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"cq",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","cq")
        status = nf90_put_att(ncid,varid,"units","nondimensional")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"tsk",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","soil/ocean temperature")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif

      if (output_zs.eq.1) then
        status = nf90_def_var(ncid,"zs",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","terrain height")
        status = nf90_put_att(ncid,varid,"units","m")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
      endif

      IF (output_dbz.eq.1) THEN
        status = nf90_def_var(ncid,"cref",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","composite reflectivity (dBZ)")
        status = nf90_put_att(ncid,varid,"units","dBZ")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      ENDIF

      if(output_sfcparams.eq.1)then
        status = nf90_def_var(ncid,"xland",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","land/water flag (1=land,2=water)")
        status = nf90_put_att(ncid,varid,"units","integer flag")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"lu",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","land use index")
        status = nf90_put_att(ncid,varid,"units","integer flag")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"mavail",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","surface moisture availability")
        status = nf90_put_att(ncid,varid,"units","integer flag")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif

      if((output_sfcparams.eq.1).and.(sfcmodel.eq.2.or.oceanmodel.eq.2))then
        status = nf90_def_var(ncid,"tmn",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","deep-layer soil temperature")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"hfx",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","heat flux at surface")
        status = nf90_put_att(ncid,varid,"units","W/m^2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"qfx",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","surface moisture flux")
        status = nf90_put_att(ncid,varid,"units","W/m^2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"gsw",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","downward SW flux at surface")
        status = nf90_put_att(ncid,varid,"units","W/m^2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"glw",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","downward LW flux at surface")
        status = nf90_put_att(ncid,varid,"units","W/m^2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
      endif

      if((output_sfcparams.eq.1).and.(sfcmodel.eq.2))then
        status = nf90_def_var(ncid,"tslb1",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","soil temp, layer 1")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"tslb2",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","soil temp, layer 2")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"tslb3",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","soil temp, layer 3")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"tslb4",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","soil temp, layer 4")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"tslb5",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","soil temp, layer 5")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
      endif

      if(output_sfcparams.eq.1.and.oceanmodel.eq.2)then
        status = nf90_def_var(ncid,"tml",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","ocean mixed layer temperature")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"hml",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","ocean mixed layer depth")
        status = nf90_put_att(ncid,varid,"units","m")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"huml",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","ocean mixed layer u vel.")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"hvml",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","ocean mixed layer v vel.")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
      endif

      if(output_radten.eq.1)then
        status = nf90_def_var(ncid,"radsw",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","solar radiation at surface")
        status = nf90_put_att(ncid,varid,"units","w/m^2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"rnflx",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","net radiation absorbed by surface")
        status = nf90_put_att(ncid,varid,"units","W/m^2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"radswnet",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","net solar radiation")
        status = nf90_put_att(ncid,varid,"units","W/m^2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"radlwin",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","incoming longwave radiation")
        status = nf90_put_att(ncid,varid,"units","W/m^2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif

        ! MS add
        status = nf90_def_var(ncid,"olr",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","TOA net outgoing longwave radiation")
        status = nf90_put_att(ncid,varid,"units","W/m^2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
        status = nf90_def_var(ncid,"dsr",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","TOA net incoming solar radiation")
        status = nf90_put_att(ncid,varid,"units","W/m^2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
#endif
        ! MS add
      endif

      IF(output_sfcdiags.eq.1)THEN
        status = nf90_def_var(ncid,"u10",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","diagnostic 10 m u wind")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"v10",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","diagnostic 10 m v wind")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"t2",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","diagnostic 2 m temperature")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"q2",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","diagnostic 2 m mixing ratio")
        status = nf90_put_att(ncid,varid,"units","kg/kg")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"znt",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","roughness length")
        status = nf90_put_att(ncid,varid,"units","m")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"ust",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","u* in similarity theory")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"hpbl",nf90_float,(/niid,njid,timeid/),varid)
      if(ipbl.eq.1)then
        status = nf90_put_att(ncid,varid,"def","PBL height (from PBL scheme)")
      else
        status = nf90_put_att(ncid,varid,"def","rough estimate of PBL height")
      endif
        status = nf90_put_att(ncid,varid,"units","m")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"zol",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","z/L (z over Monin-Obukhov length)")
        status = nf90_put_att(ncid,varid,"units","   ")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"mol",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","T* (similarity theory)")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif

        status = nf90_def_var(ncid,"br",nf90_float,(/niid,njid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","bulk Richardson number in surface layer")
        status = nf90_put_att(ncid,varid,"units"," ")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      ENDIF

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!--- 3D vars:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      if (output_zh.eq.1) then
        status = nf90_def_var(ncid,"zh",nf90_float,(/niid,njid,nkid/),varid)
        status = nf90_put_att(ncid,varid,"def","height (above ground) of scalar grid points")
        status = nf90_put_att(ncid,varid,"units","m")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif

      if(output_th.eq.1)then
        status = nf90_def_var(ncid,"th",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","potential temperature")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_thpert.eq.1)then
        status = nf90_def_var(ncid,"thpert",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","perturbation potential temperature")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif

      if(output_prs.eq.1)then
        status = nf90_def_var(ncid,"prs",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","pressure")
        status = nf90_put_att(ncid,varid,"units","Pa")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_prspert.eq.1)then
        status = nf90_def_var(ncid,"prspert",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","perturbation pressure")
        status = nf90_put_att(ncid,varid,"units","Pa")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif

      if(output_pi.eq.1)then
        status = nf90_def_var(ncid,"pi",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","nondimensional pressure")
        status = nf90_put_att(ncid,varid,"units","dimensionless")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if(output_pipert.eq.1)then
        status = nf90_def_var(ncid,"pipert",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","perturbation nondimensional pressure")
        status = nf90_put_att(ncid,varid,"units","dimensionless")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if(output_rho.eq.1)then
        status = nf90_def_var(ncid,"rho",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","density of dry air")
        status = nf90_put_att(ncid,varid,"units","kg/m^3")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if(output_rhopert.eq.1)then
        status = nf90_def_var(ncid,"rhopert",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","perturbation density of dry air")
        status = nf90_put_att(ncid,varid,"units","kg/m^3")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
    IF(iptra.eq.1)THEN
      do n=1,npt
        chid = 'pt      '
        write(chid(3:4),111) n
111     format(i2.2)
        status = nf90_def_var(ncid,chid,nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","mixing ratio of passive tracer")
        status = nf90_put_att(ncid,varid,"units","kg/kg")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      enddo
    ENDIF
    IF(imoist.eq.1)THEN
      if(output_qv.eq.1)then
        status = nf90_def_var(ncid,"qv",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","water vapor mixing ratio")
        status = nf90_put_att(ncid,varid,"units","kg/kg")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_qvpert.eq.1)then
        status = nf90_def_var(ncid,"qvpert",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","perturbation water vapor mixing ratio")
        status = nf90_put_att(ncid,varid,"units","kg/kg")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
        if(output_q.eq.1)then
        do n=1,numq
          if(n.ne.nqv)then
            status = nf90_def_var(ncid,qname(n),nf90_float,(/niid,njid,nkid,timeid/),varid)
            if(idm.eq.1.and.n.ge.nnc1.and.n.le.nnc2)then
              status = nf90_put_att(ncid,varid,"def","number concentration")
              status = nf90_put_att(ncid,varid,"units","kg^{-1}")
            elseif(idm .eq. 1 .and. n.ge.nzl1 .and. n .le. nzl2)then
              status = nf90_put_att(ncid,varid,"def","reflectivity moment")
              status = nf90_put_att(ncid,varid,"units","Z m^{-3}kg^{-1}")
            elseif(idm .eq. 1 .and. n.ge.nvl1 .and. n .le. nvl2)then
              status = nf90_put_att(ncid,varid,"def","particle volume")
              status = nf90_put_att(ncid,varid,"units","m^{3}kg^{-1}")
            else
              status = nf90_put_att(ncid,varid,"def","mixing ratio")
              status = nf90_put_att(ncid,varid,"units","kg/kg")
            endif
#ifdef NCFPLUS
            call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
            call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
          endif
        enddo
      endif
      if(output_dbz.eq.1)then
        status = nf90_def_var(ncid,"dbz",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","reflectivity")
        status = nf90_put_att(ncid,varid,"units","dBZ")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
    ENDIF

      if(output_buoyancy.eq.1)then
        status = nf90_def_var(ncid,"buoyancy",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","buoyancy")
        status = nf90_put_att(ncid,varid,"units","m s^-2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif

      if(output_uinterp.eq.1)then
        status = nf90_def_var(ncid,"uinterp",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","velocity in x-direction, interpolated to scalar points")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_vinterp.eq.1)then
        status = nf90_def_var(ncid,"vinterp",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","velocity in y-direction, interpolated to scalar points")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_winterp.eq.1)then
        status = nf90_def_var(ncid,"winterp",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","velocity in z-direction, interpolated to scalar points")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif

      if(output_vort.eq.1)then
        status = nf90_def_var(ncid,"xvort",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","horizontal vorticity (x)")
        status = nf90_put_att(ncid,varid,"units","s^-1")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"yvort",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","horizontal vorticity (y)")
        status = nf90_put_att(ncid,varid,"units","s^-1")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"zvort",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","vertical vorticity")
        status = nf90_put_att(ncid,varid,"units","s^-1")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif

      if (output_basestate.eq.1) then

        status = nf90_def_var(ncid,"pi0",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","base-state nondimensional pressure")
        status = nf90_put_att(ncid,varid,"units","dimensionless")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"th0",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","base-state potential temperature")
        status = nf90_put_att(ncid,varid,"units","K")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"prs0",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","base-state pressure")
        status = nf90_put_att(ncid,varid,"units","Pa")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"qv0",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","base-state water vapor mixing ratio")
        status = nf90_put_att(ncid,varid,"units","kg/kg")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif

      endif

      if(output_pblten.eq.1)then
        status = nf90_def_var(ncid,"thpten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","pbl tendency: theta")
        status = nf90_put_att(ncid,varid,"units","    ")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"qvpten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","pbl tendency: qv")
        status = nf90_put_att(ncid,varid,"units","    ")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"qcpten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","pbl tendency: qc")
        status = nf90_put_att(ncid,varid,"units","    ")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"qipten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","pbl tendency: qi")
        status = nf90_put_att(ncid,varid,"units","    ")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"upten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","pbl tendency: u")
        status = nf90_put_att(ncid,varid,"units","    ")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"vpten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","pbl tendency: v")
        status = nf90_put_att(ncid,varid,"units","    ")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif

      if(output_radten.eq.1)then
        status = nf90_def_var(ncid,"swten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","pot temp tendency, sw rad")
        status = nf90_put_att(ncid,varid,"units","    ")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif

        status = nf90_def_var(ncid,"lwten",nf90_float,(/niid,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","pot temp tendency, lw rad")
        status = nf90_put_att(ncid,varid,"units","K/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif

      if(output_u.eq.1)then
        status = nf90_def_var(ncid,"u",nf90_float,(/nip1id,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","velocity in x-direction")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_upert.eq.1)then
        status = nf90_def_var(ncid,"upert",nf90_float,(/nip1id,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","perturbation velocity in x-direction")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if (output_basestate.eq.1) then
        status = nf90_def_var(ncid,"u0",nf90_float,(/nip1id,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","base-state x-component of velocity")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if(output_turbten.eq.1)then
        status = nf90_def_var(ncid,"ftu",nf90_float,(/nip1id,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","u tendency: turbulence scheme")
        status = nf90_put_att(ncid,varid,"units","m/s/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if(output_impdiften.eq.1)then
        status = nf90_def_var(ncid,"fdu",nf90_float,(/nip1id,njid,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","u tendency: implicit diffusion")
        status = nf90_put_att(ncid,varid,"units","m/s/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif

      if(output_v.eq.1)then
        status = nf90_def_var(ncid,"v",nf90_float,(/niid,njp1id,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","velocity in y-direction")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      if(output_vpert.eq.1)then
        status = nf90_def_var(ncid,"vpert",nf90_float,(/niid,njp1id,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","perturbation velocity in y-direction")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if (output_basestate.eq.1) then
        status = nf90_def_var(ncid,"v0",nf90_float,(/niid,njp1id,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","base-state y-component of velocity")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if(output_turbten.eq.1)then
        status = nf90_def_var(ncid,"ftv",nf90_float,(/niid,njp1id,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","v tendency: turbulence scheme")
        status = nf90_put_att(ncid,varid,"units","m/s/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if(output_impdiften.eq.1)then
        status = nf90_def_var(ncid,"fdv",nf90_float,(/niid,njp1id,nkid,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","v tendency: implicit diffusion")
        status = nf90_put_att(ncid,varid,"units","m/s/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif

      if(output_w.eq.1)then
        status = nf90_def_var(ncid,"w",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","velocity in z-direction")
        status = nf90_put_att(ncid,varid,"units","m/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      endif
      IF((iturb.eq.1).and.(output_tke.eq.1))THEN
        status = nf90_def_var(ncid,"tke",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","subgrid turbulence kinetic energy")
        status = nf90_put_att(ncid,varid,"units","m^2/s^2")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
      ENDIF
      IF(output_km.eq.1)THEN
        !----
        status = nf90_def_var(ncid,"kmh",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","eddy mixing coefficient for momentum in the horizontal direction")
        status = nf90_put_att(ncid,varid,"units","m^2/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
        !----
        status = nf90_def_var(ncid,"kmv",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","eddy mixing coefficient for momentum in the vertical direction")
        status = nf90_put_att(ncid,varid,"units","m^2/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
        !----
      ENDIF
      IF(output_kh.eq.1)THEN
        !----
        status = nf90_def_var(ncid,"khh",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","eddy mixing coefficient for scalars in the horizontal direction")
        status = nf90_put_att(ncid,varid,"units","m^2/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
        !----
        status = nf90_def_var(ncid,"khv",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","eddy mixing coefficient for scalars in the vertical direction")
        status = nf90_put_att(ncid,varid,"units","m^2/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
!!!        call disp_err( nf90_def_var_deflate(ncid,varid,shuffle,deflate,deflate_level) , .true. )
#endif
        !----
      ENDIF

      if(output_dissten.eq.1)then
        status = nf90_def_var(ncid,"dissten",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","dissipation rate")
        status = nf90_put_att(ncid,varid,"units","m^2/s^3")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if(output_nm.eq.1)then
        status = nf90_def_var(ncid,"nm",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","squared Brunt-Vaisala freq")
        status = nf90_put_att(ncid,varid,"units","s^{-2}")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if(output_def.eq.1)then
        status = nf90_def_var(ncid,"defv",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","vertical deformation")
        status = nf90_put_att(ncid,varid,"units","s^{-2}")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
        status = nf90_def_var(ncid,"defh",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","horizontal deformation")
        status = nf90_put_att(ncid,varid,"units","s^{-2}")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if(output_turbten.eq.1)then
        status = nf90_def_var(ncid,"ftw",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","w tendency: turbulence scheme")
        status = nf90_put_att(ncid,varid,"units","m/s/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif
      if(output_impdiften.eq.1)then
        status = nf90_def_var(ncid,"fdw",nf90_float,(/niid,njid,nkp1id,timeid/),varid)
        status = nf90_put_att(ncid,varid,"def","w tendency: implicit diffusion")
        status = nf90_put_att(ncid,varid,"units","m/s/s")
#ifdef NCFPLUS
        call disp_err( NF90_DEF_VAR_CHUNKING(ncid,varid,NF90_CHUNKED,(/chkx,chky,1,1/)) , .true. )
#endif
      endif

!--------------------------------------------------

      status = nf90_put_att(ncid,NF90_GLOBAL,'Conventions','COARDS')
      status = nf90_put_att(ncid,NF90_GLOBAL,'x_units','km')
      status = nf90_put_att(ncid,NF90_GLOBAL,'x_label','x')
      status = nf90_put_att(ncid,NF90_GLOBAL,'y_units','km')
      status = nf90_put_att(ncid,NF90_GLOBAL,'y_label','y')
      status = nf90_put_att(ncid,NF90_GLOBAL,'z_units','km')
      status = nf90_put_att(ncid,NF90_GLOBAL,'z_label','z')

      status = nf90_enddef(ncid)

! ... end of defs
!--------------------------------------------------
! begin data ... initial time ...

    IF(icor.eq.1)THEN
      status = nf90_inq_varid(ncid,'f_cor',varid)
      status = nf90_put_var(ncid,varid,fcor)
    ENDIF

    if(.not.coards)then
      status = nf90_inq_varid(ncid,'ztop',varid)
      status = nf90_put_var(ncid,varid,0.001*ztop)
    endif

      if(coards)then
        status = nf90_inq_varid(ncid,'ni',varid)
      else
        status = nf90_inq_varid(ncid,'xh',varid)
      endif
      do i=1,nx
        dumx(i) = 0.001*0.5*(xfref(i)+xfref(i+1))
      enddo
      status = nf90_put_var(ncid,varid,dumx,(/1/),(/nx/))

      if(coards)then
        status = nf90_inq_varid(ncid,'nip1',varid)
      else
        status = nf90_inq_varid(ncid,'xf',varid)
      endif
      do i=1,nx+1
        dumx(i) = 0.001*xfref(i)
      enddo
      status = nf90_put_var(ncid,varid,dumx,(/1/),(/nx+1/))

      if(coards)then
        status = nf90_inq_varid(ncid,'nj',varid)
      else
        status = nf90_inq_varid(ncid,'yh',varid)
      endif
      do j=1,ny
        dumy(j) = 0.001*0.5*(yfref(j)+yfref(j+1))
      enddo
      status = nf90_put_var(ncid,varid,dumy,(/1/),(/ny/))

      if(coards)then
        status = nf90_inq_varid(ncid,'njp1',varid)
      else
        status = nf90_inq_varid(ncid,'yf',varid)
      endif
      do j=1,ny+1
        dumy(j) = 0.001*yfref(j)
      enddo
      status = nf90_put_var(ncid,varid,dumy,(/1/),(/ny+1/))

      if(coards)then
        status = nf90_inq_varid(ncid,'nk',varid)
      else
        status = nf90_inq_varid(ncid,'z',varid)
      endif
      if(terrain_flag)then
        do k=1,nk
          dumz(k) = 0.001*sigma(k)
        enddo
        status = nf90_put_var(ncid,varid,dumz,(/1/),(/nk/))
      else
        do k=1,nk
          dumz(k) = 0.001*zh(1,1,k)
        enddo
        status = nf90_put_var(ncid,varid,dumz,(/1/),(/nk/))
      endif

      if(coards)then
        status = nf90_inq_varid(ncid,'nkp1',varid)
      else
        status = nf90_inq_varid(ncid,'zf',varid)
      endif
      if(terrain_flag)then
        do k=1,nk+1
          dumz(k) = 0.001*sigmaf(k)
        enddo
        status = nf90_put_var(ncid,varid,dumz,(/1/),(/nk+1/))
      else
        do k=1,nk+1
          dumz(k) = 0.001*zf(1,1,k)
        enddo
        status = nf90_put_var(ncid,varid,dumz,(/1/),(/nk+1/))
      endif

      ! ... end if info at initial time only

!----------------------------------------------------------

    ELSE

      ! Not initial time ... open file ...

      call disp_err( nf90_open(string,nf90_write,ncid), .true. )

    ENDIF ifallinfo

      status = nf90_inq_varid(ncid,'time',varid)
      call checkstatus(status)
      status = nf90_put_var(ncid,varid,actual_time,(/time_index/))
      call checkstatus(status)

#endif

    return
    end subroutine netcdf_prelim


#ifdef NETCDF
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine checkstatus(status)
      use netcdf
      implicit none

      integer :: status

      if(status.ne.nf90_noerr)then
        print *,'  Error ... '
        print *,nf90_strerror(status)
        call stopcm1
      endif

      return
      end subroutine checkstatus

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine write2d_nc(chid,ncid,time_index,ni,nj,d2d)
      use netcdf
      implicit none

      character*8, intent(in) :: chid
      integer, intent(in) :: ncid,time_index,ni,nj
      real, dimension(ni,nj), intent(in) :: d2d

      integer :: varid,status

!----------------------------------

      status = nf90_inq_varid(ncid,chid,varid)
      if(status.ne.nf90_noerr)then
        print *,'  Error1 in write2d_nc, chid = ',chid
        print *,nf90_strerror(status)
        call stopcm1
      endif

      status = nf90_put_var(ncid,varid,d2d,(/1,1,time_index/),(/ni,nj,1/))
      if(status.ne.nf90_noerr)then
        print *,'  Error2 in write2d_nc, chid = ',chid
        print *,nf90_strerror(status)
        call stopcm1
      endif

!----------------------------------

      return
      end subroutine write2d_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine write3d_nc(chid,k,ncid,time_index,ni,nj,ds)
      use netcdf
      implicit none

      character*8, intent(in) :: chid
      integer, intent(in) :: k,ncid,time_index,ni,nj
      real, dimension(ni,nj), intent(in) :: ds

      integer :: varid,status

!----------------------------------

      status = nf90_inq_varid(ncid,chid,varid)
      if(status.ne.nf90_noerr)then
        print *,'  Error1 in write3d_nc, chid = ',chid
        print *,nf90_strerror(status)
        call stopcm1
      endif

      status = nf90_put_var(ncid,varid,ds,(/1,1,k,time_index/),(/ni,nj,1,1/))
      if(status.ne.nf90_noerr)then
        print *,'  Error2 in write3d_nc, chid = ',chid
        print *,nf90_strerror(status)
        call stopcm1
      endif

!----------------------------------

      return
      end subroutine write3d_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine writestat_nc(nrec,rtime,nstat,rstat,qname,budname)
      use netcdf
      implicit none

      include 'input.incl'

      integer, intent(inout) :: nrec
      real,    intent(in)    :: rtime
      integer, intent(in)    :: nstat
      real, dimension(stat_out), intent(in) :: rstat
      character*3, dimension(maxq), intent(in) :: qname
      character*6, dimension(maxq), intent(in) :: budname

      integer :: n,ncid,status,dimid,varid,time_index
      character*8  :: text1
      character*30 :: text2

    string(totlen+1:totlen+22) = '_stats.nc             '

    print *,string

  IF(nrec.eq.1)THEN
    ! Definitions/descriptions:

#ifdef NCFPLUS
!--- works with netcdf 4.2, but not 4.0 (grumble)
      call disp_err( nf90_create(string,IOR(nf90_netcdf4, nf90_classic_model),ncid) , .true. )
#else
      ! must do this for netcdf 4.0 (it seems) !
      call disp_err( nf90_create(path=string,cmode=0,ncid=ncid) , .true. )
#endif

    status = nf90_def_dim(ncid,"xh",1,dimid)
    status = nf90_def_dim(ncid,"yh",1,dimid)
    status = nf90_def_dim(ncid,"zh",1,dimid)
    status = nf90_def_dim(ncid,"time",nf90_unlimited,dimid)

    status = nf90_def_var(ncid,"xh",nf90_float,(/1/),varid)
    status = nf90_put_att(ncid,varid,"def","west-east location")
    status = nf90_put_att(ncid,varid,"units","degree_east")

    status = nf90_def_var(ncid,"yh",nf90_float,(/2/),varid)
    status = nf90_put_att(ncid,varid,"def","south-north location")
    status = nf90_put_att(ncid,varid,"units","degree_north")

    status = nf90_def_var(ncid,"zh",nf90_float,(/3/),varid)
    status = nf90_put_att(ncid,varid,"def","height")
    status = nf90_put_att(ncid,varid,"units","m")

    status = nf90_def_var(ncid,"time",nf90_float,(/4/),varid)
    status = nf90_put_att(ncid,varid,"def","time")
    status = nf90_put_att(ncid,varid,"units","seconds")

    IF(adapt_dt.eq.1)THEN
      status = nf90_def_var(ncid,"dt",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","average timestep dt           ")
      status = nf90_put_att(ncid,varid,"units","seconds")
    ENDIF
    IF(stat_w.eq.1)THEN
      status = nf90_def_var(ncid,"wmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","maximum vertical velocity     ")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"wmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","minimum vertical velocity     ")
      status = nf90_put_att(ncid,varid,"units","m/s")
    ENDIF
    IF(stat_u.eq.1)THEN
      status = nf90_def_var(ncid,"umax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max E-W velocity              ")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"umin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min E-W velocity              ")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"sumax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max E-W velocity at lowest lvl")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"sumin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min E-W velocity at lowest lvl")
      status = nf90_put_att(ncid,varid,"units","m/s")
    ENDIF
    IF(stat_v.eq.1)THEN
      status = nf90_def_var(ncid,"vmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max N-S velocity              ")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"vmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min N-S velocity              ")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"svmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max N-S velocity at lowest lvl")
      status = nf90_put_att(ncid,varid,"units","m/s")
      status = nf90_def_var(ncid,"svmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min N-S velocity at lowest lvl")
      status = nf90_put_att(ncid,varid,"units","m/s")
    ENDIF
    IF(stat_rmw.eq.1)THEN
      status = nf90_def_var(ncid,"rmw",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","radius of maximum windspeed")
      status = nf90_put_att(ncid,varid,"units","m")
      status = nf90_def_var(ncid,"zmw",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","height (ASL) of maximum windspeed")
      status = nf90_put_att(ncid,varid,"units","m")
    ENDIF
    IF(stat_pipert.eq.1)THEN
      status = nf90_def_var(ncid,"ppimax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max pi pert.                  ")
      status = nf90_put_att(ncid,varid,"units","nondimensional")
      status = nf90_def_var(ncid,"ppimin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min pi pert.                  ")
      status = nf90_put_att(ncid,varid,"units","nondimensional")
    ENDIF
    IF(stat_prspert.eq.1)THEN
      status = nf90_def_var(ncid,"ppmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max pressure pert.            ")
      status = nf90_put_att(ncid,varid,"units","Pa")
      status = nf90_def_var(ncid,"ppmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min pressure pert.            ")
      status = nf90_put_att(ncid,varid,"units","Pa")
    ENDIF
    IF(stat_thpert.eq.1)THEN
      status = nf90_def_var(ncid,"thpmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max potential temp. pert.     ")
      status = nf90_put_att(ncid,varid,"units","K")
      status = nf90_def_var(ncid,"thpmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min potential temp. pert.     ")
      status = nf90_put_att(ncid,varid,"units","K")
      status = nf90_def_var(ncid,"sthpmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max pot temp pert lowest level")
      status = nf90_put_att(ncid,varid,"units","K")
      status = nf90_def_var(ncid,"sthpmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min pot temp pert lowest level")
      status = nf90_put_att(ncid,varid,"units","K")
    ENDIF
    IF(stat_q.eq.1)THEN
      do n=1,numq
        text1='max     '
        text2='max                           '
        write(text1(4:6),156) qname(n)
        write(text2(5:7),156) qname(n)
156     format(a3)
        status = nf90_def_var(ncid,text1,nf90_float,4,varid)
        status = nf90_put_att(ncid,varid,"def",text2)
        status = nf90_put_att(ncid,varid,"units","kg/kg")
        text1='min     '
        text2='min                           '
        write(text1(4:6),156) qname(n)
        write(text2(5:7),156) qname(n)
        status = nf90_def_var(ncid,text1,nf90_float,4,varid)
        status = nf90_put_att(ncid,varid,"def",text2)
        status = nf90_put_att(ncid,varid,"units","kg/kg")
      enddo
    ENDIF
    IF(stat_tke.eq.1)THEN
      status = nf90_def_var(ncid,"tkemax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max tke                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s^2")
      status = nf90_def_var(ncid,"tkemin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min tke                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s^2")
    ENDIF
    IF(stat_km.eq.1)THEN
      status = nf90_def_var(ncid,"kmhmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max kmh                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
      status = nf90_def_var(ncid,"kmhmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min kmh                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
      status = nf90_def_var(ncid,"kmvmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max kmv                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
      status = nf90_def_var(ncid,"kmvmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min kmv                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
    ENDIF
    IF(stat_kh.eq.1)THEN
      status = nf90_def_var(ncid,"khhmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max khh                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
      status = nf90_def_var(ncid,"khhmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min khh                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
      status = nf90_def_var(ncid,"khvmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max khv                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
      status = nf90_def_var(ncid,"khvmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min khv                       ")
      status = nf90_put_att(ncid,varid,"units","m^2/s")
    ENDIF
    IF(stat_div.eq.1)THEN
      status = nf90_def_var(ncid,"divmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max 3d divergence             ")
      status = nf90_def_var(ncid,"divmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min 3d divergence             ")
    ENDIF
    IF(stat_rh.eq.1)THEN
      status = nf90_def_var(ncid,"rhmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max relative humidity         ")
      status = nf90_def_var(ncid,"rhmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min relative humidity         ")
    ENDIF
    IF(stat_rhi.eq.1)THEN
      status = nf90_def_var(ncid,"rhimax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max relative humidity wrt ice ")
      status = nf90_def_var(ncid,"rhimin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min relative humidity wrt ice ")
    ENDIF
    IF(iptra.eq.1)then
      do n=1,npt
        text1='maxpt   '
        text2='max pt                        '
        write(text1(6:6),157) n
        write(text2(7:7),157) n
157     format(i1)
        status = nf90_def_var(ncid,text1,nf90_float,4,varid)
        status = nf90_put_att(ncid,varid,"def",text2)
        text1='minpt   '
        text2='min pt                        '
        write(text1(6:6),157) n
        write(text2(7:7),157) n
        status = nf90_def_var(ncid,text1,nf90_float,4,varid)
        status = nf90_put_att(ncid,varid,"def",text2)
      enddo
    endif
    IF(stat_the.eq.1)THEN
      status = nf90_def_var(ncid,"themax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max theta-e below 10 km       ")
      status = nf90_put_att(ncid,varid,"units","K")
      status = nf90_def_var(ncid,"themin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min theta-e below 10 km       ")
      status = nf90_put_att(ncid,varid,"units","K")
      status = nf90_def_var(ncid,"sthemax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max theta-e at lowest level   ")
      status = nf90_put_att(ncid,varid,"units","K")
      status = nf90_def_var(ncid,"sthemin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min theta-e at lowest level   ")
      status = nf90_put_att(ncid,varid,"units","K")
    ENDIF
    IF(stat_cloud.eq.1)THEN
      status = nf90_def_var(ncid,"qctop",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max cloud top height          ")
      status = nf90_put_att(ncid,varid,"units","m")
      status = nf90_def_var(ncid,"qcbot",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min cloud base height         ")
      status = nf90_put_att(ncid,varid,"units","m")
    ENDIF
    IF(stat_sfcprs.eq.1)THEN
      status = nf90_def_var(ncid,"sprsmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max pressure at lowest level  ")
      status = nf90_put_att(ncid,varid,"units","Pa")
      status = nf90_def_var(ncid,"sprsmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min pressure at lowest level  ")
      status = nf90_put_att(ncid,varid,"units","Pa")
    ENDIF
    IF(stat_wsp.eq.1)THEN
      status = nf90_def_var(ncid,"wspmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max wind speed                ")
      status = nf90_put_att(ncid,varid,"units","m/s")

      status = nf90_def_var(ncid,"wspmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min wind speed                ")
      status = nf90_put_att(ncid,varid,"units","m/s")

      status = nf90_def_var(ncid,"swspmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max wind speed at lowest level")
      status = nf90_put_att(ncid,varid,"units","m/s")

      status = nf90_def_var(ncid,"swspmin",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min wind speed at lowest level")
      status = nf90_put_att(ncid,varid,"units","m/s")

    IF(bbc.eq.3)THEN
      status = nf90_def_var(ncid,"wsp10max",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max 10 m wind speed           ")
      status = nf90_put_att(ncid,varid,"units","m/s")

      status = nf90_def_var(ncid,"wsp10min",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","min 10 m wind speed           ")
      status = nf90_put_att(ncid,varid,"units","m/s")
    ENDIF
    ENDIF
    IF(stat_cfl.eq.1)THEN
    IF(adapt_dt.eq.1)THEN
      status = nf90_def_var(ncid,"cflmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max Courant number (average)  ")
      status = nf90_put_att(ncid,varid,"units","nondimensional")
    ELSE
      status = nf90_def_var(ncid,"cflmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max Courant number            ")
      status = nf90_put_att(ncid,varid,"units","nondimensional")
    ENDIF
      status = nf90_def_var(ncid,"kshmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max horiz K stability factor  ")
      status = nf90_put_att(ncid,varid,"units","nondimensional")
      status = nf90_def_var(ncid,"ksvmax",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max vert K stability factor   ")
      status = nf90_put_att(ncid,varid,"units","nondimensional")
    ENDIF
    IF(stat_vort.eq.1)THEN
      status = nf90_def_var(ncid,"vortsfc",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max vert. vort. at lowest lvl ")
      status = nf90_put_att(ncid,varid,"units","1/s")
      status = nf90_def_var(ncid,"vort1km",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max vert. vort. at z = 1 km   ")
      status = nf90_put_att(ncid,varid,"units","1/s")
      status = nf90_def_var(ncid,"vort2km",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max vert. vort. at z = 2 km   ")
      status = nf90_put_att(ncid,varid,"units","1/s")
      status = nf90_def_var(ncid,"vort3km",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max vert. vort. at z = 3 km   ")
      status = nf90_put_att(ncid,varid,"units","1/s")
      status = nf90_def_var(ncid,"vort4km",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max vert. vort. at z = 4 km   ")
      status = nf90_put_att(ncid,varid,"units","1/s")
      status = nf90_def_var(ncid,"vort5km",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","max vert. vort. at z = 5 km   ")
      status = nf90_put_att(ncid,varid,"units","1/s")
    ENDIF
    IF(stat_tmass.eq.1)THEN
      status = nf90_def_var(ncid,"tmass",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","total mass of (dry) air       ")
      status = nf90_put_att(ncid,varid,"units","kg")
    ENDIF
    IF(stat_tmois.eq.1)THEN
      status = nf90_def_var(ncid,"tmois",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","total moisture                ")
      status = nf90_put_att(ncid,varid,"units","kg")
    ENDIF
    IF(stat_qmass.eq.1)THEN
      do n=1,numq
        IF( (n.eq.nqv) .or.                                 &
            (n.ge.nql1.and.n.le.nql2) .or.                  &
            (n.ge.nqs1.and.n.le.nqs2.and.iice.eq.1) )THEN
          text1='mass    '
          text2='total mass of                 '
          write(text1( 5: 7),156) qname(n)
          write(text2(15:17),156) qname(n)
          status = nf90_def_var(ncid,text1,nf90_float,4,varid)
          status = nf90_put_att(ncid,varid,"def",text2)
          status = nf90_put_att(ncid,varid,"units","kg")
        ENDIF
      enddo
    ENDIF
    IF(stat_tenerg.eq.1)THEN
      status = nf90_def_var(ncid,"ek",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","total kinetic energy          ")
      status = nf90_def_var(ncid,"ei",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","total internal energy         ")
      status = nf90_def_var(ncid,"ep",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","total potential energy        ")
      status = nf90_def_var(ncid,"le",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","total latent energy           ")
      status = nf90_def_var(ncid,"et",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","total energy                  ")
    ENDIF
    IF(stat_mo.eq.1)THEN
      status = nf90_def_var(ncid,"tmu",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","total E-W momentum            ")
      status = nf90_def_var(ncid,"tmv",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","total N-S momentum            ")
      status = nf90_def_var(ncid,"tmw",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","total vertical momentum       ")
    ENDIF
    IF(stat_tmf.eq.1)THEN
      status = nf90_def_var(ncid,"tmfu",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","total upward mass flux        ")
      status = nf90_def_var(ncid,"tmfd",nf90_float,4,varid)
      status = nf90_put_att(ncid,varid,"def","total downward mass flux      ")
    ENDIF
    IF(stat_pcn.eq.1)THEN
      do n=1,nbudget
        text1='        '
        text2='                              '
        write(text1(1:6),158) budname(n)
        write(text2(1:6),158) budname(n)
158     format(a6)
        status = nf90_def_var(ncid,text1,nf90_float,4,varid)
        status = nf90_put_att(ncid,varid,"def",text2)
      enddo
    ENDIF
    IF(stat_qsrc.eq.1)THEN
      do n=1,numq
        text1='as      '
        text2='artificial source of          '
        write(text1( 3: 5),156) qname(n)
        write(text2(22:24),156) qname(n)
        status = nf90_def_var(ncid,text1,nf90_float,4,varid)
        status = nf90_put_att(ncid,varid,"def",text2)
      enddo
      do n=1,numq
        text1='bs      '
        text2='bndry source/sink of          '
        write(text1( 3: 5),156) qname(n)
        write(text2(22:24),156) qname(n)
        status = nf90_def_var(ncid,text1,nf90_float,4,varid)
        status = nf90_put_att(ncid,varid,"def",text2)
      enddo
    ENDIF

    status = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','COARDS')

    status = nf90_enddef(ncid)

    status = nf90_put_var(ncid,1,0.0)
    status = nf90_put_var(ncid,2,0.0)
    status = nf90_put_var(ncid,3,0.0)

  ELSE

    ! open file:

    call disp_err( nf90_open(string,nf90_write,ncid), .true. )

  ENDIF

    ! Write data:

    time_index = nrec

    status = nf90_put_var(ncid,4,rtime,(/time_index/))

    DO n=1,nstat
      varid = 4 + n
      status = nf90_put_var(ncid,varid,rstat(n),(/time_index/))
    ENDDO

    ! close file

    call disp_err( nf90_close(ncid) , .true. )

    nrec = nrec + 1

    ! all done

      return
      end subroutine writestat_nc

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine disp_err( status , stop_on_error )
      use netcdf
      implicit none

      integer, intent(in) :: status
      logical, intent(in) :: stop_on_error

      IF( status.ne.nf90_noerr )THEN
        IF( stop_on_error )THEN
          print *,'  netcdf status returned an error: ', status,' ... stopping program'
          call stopcm1
        ENDIF
      ENDIF

      return
      end subroutine disp_err

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine writepdata_nc(prec,rtime,pdata)
      use netcdf
      implicit none

      include 'input.incl'

      integer :: prec
      real :: rtime
      real, dimension(npvals,nparcels) :: pdata

      integer :: ncid,status,dimid,varid,time_index,n,np
      character*8 text1

!-----------------------------------------------------------------------

  IF(prec.eq.1)THEN
    ! Definitions/descriptions:

#ifdef NCFPLUS
!--- works with netcdf 4.2, but not 4.0 (grumble)
      call disp_err( nf90_create('cm1out_pdata.nc',IOR(nf90_netcdf4, nf90_classic_model),ncid) , .true. )
#else
      ! must do this for netcdf 4.0 (it seems) !
      call disp_err( nf90_create(path='cm1out_pdata.nc',cmode=0,ncid=ncid) , .true. )
#endif

    status = nf90_def_dim(ncid,"xh",nparcels,dimid)
    status = nf90_def_dim(ncid,"yh",1,dimid)
    status = nf90_def_dim(ncid,"zh",1,dimid)
    status = nf90_def_dim(ncid,"time",nf90_unlimited,dimid)

    status = nf90_def_var(ncid,"xh",nf90_float,(/1/),varid)
    status = nf90_put_att(ncid,varid,"def","west-east location ... actually, really parcel info")
    status = nf90_put_att(ncid,varid,"units","degree_east")

    status = nf90_def_var(ncid,"yh",nf90_float,(/2/),varid)
    status = nf90_put_att(ncid,varid,"def","south-north location")
    status = nf90_put_att(ncid,varid,"units","degree_north")

    status = nf90_def_var(ncid,"zh",nf90_float,(/3/),varid)
    status = nf90_put_att(ncid,varid,"def","height")
    status = nf90_put_att(ncid,varid,"units","m")

    status = nf90_def_var(ncid,"time",nf90_float,(/4/),varid)
    status = nf90_put_att(ncid,varid,"def","time")
    status = nf90_put_att(ncid,varid,"units","seconds")

!------------------------

    status = nf90_def_var(ncid,"x",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","x                             ")
    status = nf90_put_att(ncid,varid,"units","m")

    status = nf90_def_var(ncid,"y",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","y                             ")
    status = nf90_put_att(ncid,varid,"units","m")

    status = nf90_def_var(ncid,"z",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","z                             ")
    status = nf90_put_att(ncid,varid,"units","m")

    status = nf90_def_var(ncid,"u",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","u                             ")
    status = nf90_put_att(ncid,varid,"units","m/s")

    status = nf90_def_var(ncid,"v",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","v                             ")
    status = nf90_put_att(ncid,varid,"units","m/s")

    status = nf90_def_var(ncid,"w",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","w                             ")
    status = nf90_put_att(ncid,varid,"units","m/s")

    status = nf90_def_var(ncid,"th",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","potential temperature         ")
    status = nf90_put_att(ncid,varid,"units","K")

    status = nf90_def_var(ncid,"prs",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","pressure                      ")
    status = nf90_put_att(ncid,varid,"units","Pa")

    status = nf90_def_var(ncid,"b",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","buoyancy                      ")
    status = nf90_put_att(ncid,varid,"units","m/s^2")

!-----

    if(iptra.eq.1)then
      do n=1,npt
        text1='pt      '
        if(n.lt.10)then
          write(text1(3:3),155) n
        else
          write(text1(3:4),154) n
        endif
154     format(i2.2)
155     format(i1.1)
        status = nf90_def_var(ncid,text1,nf90_float,(/1,4/),varid)
        status = nf90_put_att(ncid,varid,"def","passive tracer conc.")
        status = nf90_put_att(ncid,varid,"units","kg/kg")
      enddo
    endif

!-----

  if(prqv.ge.1)then
    status = nf90_def_var(ncid,"qv",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","water vapor mixing ratio")
    status = nf90_put_att(ncid,varid,"units","kg/kg")
  endif
  if(prqc.ge.1)then
    status = nf90_def_var(ncid,"qc",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","cloud water mixing ratio")
    status = nf90_put_att(ncid,varid,"units","kg/kg")
  endif
  if(prqr.ge.1)then
    status = nf90_def_var(ncid,"qr",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","rain water mixing ratio")
    status = nf90_put_att(ncid,varid,"units","kg/kg")
  endif
  if(prqi.ge.1)then
    status = nf90_def_var(ncid,"qi",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","mixing ratio of ice crystals")
    status = nf90_put_att(ncid,varid,"units","kg/kg")
  endif
  if(prqs.ge.1)then
    status = nf90_def_var(ncid,"qs",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","mixing ratio of snow")
    status = nf90_put_att(ncid,varid,"units","kg/kg")
  endif
  if(prqg.ge.1)then
    status = nf90_def_var(ncid,"qg",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","mixing ratio of graupel")
    status = nf90_put_att(ncid,varid,"units","kg/kg")
  endif

  if(prkh.ge.1)then
    status = nf90_def_var(ncid,"kh",nf90_float,(/1,4/),varid)
    status = nf90_put_att(ncid,varid,"def","turb. coef. for scalar        ")
    status = nf90_put_att(ncid,varid,"units","m^2/s")
  endif

!------------------------

    status = nf90_put_att(ncid, NF90_GLOBAL, 'Conventions','COARDS')

    status = nf90_enddef(ncid)

  do np=1,nparcels
    status = nf90_put_var(ncid,1,float(np),(/np/))
  enddo
    status = nf90_put_var(ncid,2,0.0)
    status = nf90_put_var(ncid,3,0.0)

!------------------------

  ELSE

    ! open file:

    call disp_err( nf90_open('cm1out_pdata.nc',nf90_write,ncid), .true. )

  ENDIF

      ! Write data:

      time_index = prec

      status = nf90_put_var(ncid,4,rtime,(/time_index/))

      DO np=1,nparcels
      DO n=1,npvals
        varid = 4 + n
        status = nf90_put_var(ncid,varid,pdata(n,np),(/np,time_index/))
      ENDDO
      ENDDO

      ! close file

      call disp_err( nf90_close(ncid) , .true. )

      prec = prec + 1

      ! all done

      return
      end subroutine writepdata_nc
#endif
