!-----------------------------------------------------------------------
!  $Id$
module stats_subs

  implicit none

  private ! Set Default Scope

  public :: stats_init, stats_begin_timestep, stats_end_timestep, & 
    stats_accumulate, stats_finalize

  private :: stats_zero, stats_avg

  contains

  !-----------------------------------------------------------------------
  subroutine stats_init( iunit, fname_prefix, fdir, l_stats_in, &
                         stats_fmt_in, stats_tsamp_in, stats_tout_in, fnamelist, &
                         nnzp, gzt, gzm, nnrad_zt, &
                         grad_zt, nnrad_zm, grad_zm, day, month, year, &
                         rlat, rlon, time_current, delt )
    !
    !     Description: Initializes the statistics saving functionality of
    !     the CLUBB model.
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
      zt,      & ! Variables
      ztscr01, & 
      ztscr02, & 
      ztscr03, & 
      ztscr04, & 
      ztscr05, & 
      ztscr06, & 
      ztscr07, & 
      ztscr08, & 
      ztscr09, & 
      ztscr10, & 
      ztscr11, & 
      ztscr12, & 
      ztscr13, & 
      ztscr14, & 
      ztscr15, & 
      ztscr16, & 
      ztscr17, & 
      ztscr18, & 
      ztscr19, & 
      ztscr20, & 
      ztscr21

    use stats_variables, only: & 
      zm,      & 
      zmscr01, & 
      zmscr02, & 
      zmscr03, & 
      zmscr04, & 
      zmscr05, & 
      zmscr06, & 
      zmscr07, & 
      zmscr08, & 
      zmscr09, & 
      zmscr10, & 
      zmscr11, & 
      zmscr12, & 
      zmscr13, & 
      zmscr14, & 
      zmscr15, &
      zmscr16, &
      zmscr17, &
      rad_zt,  &
      rad_zm,  &
      sfc,     & 
      l_stats, &
      l_output_rad_files, & 
      stats_tsamp,   & 
      stats_tout,    & 
      l_stats_samp,  & 
      l_stats_first, & 
      l_stats_last, & 
      fname_zt, & 
      fname_zm, &
      fname_rad_zt, &
      fname_rad_zm, & 
      fname_sfc, & 
      l_netcdf, & 
      l_grads

    use stats_precision, only: & 
      time_precision   ! Variable(s)

    use output_grads, only: & 
      open_grads  ! Procedure

#ifdef NETCDF
    use output_netcdf, only: & 
      open_netcdf     ! Procedure
#endif

    use stats_zm, only: &
      nvarmax_zm, & ! Constant(s) 
      stats_init_zm ! Procedure(s)

    use stats_zt, only: & 
      nvarmax_zt, & ! Constant(s)
      stats_init_zt ! Procedure(s)

    use stats_rad_zt, only: & 
      nvarmax_rad_zt, & ! Constant(s)
      stats_init_rad_zt ! Procedure(s)

    use stats_rad_zm, only: & 
      nvarmax_rad_zm, & ! Constant(s)
      stats_init_rad_zm ! Procedure(s)       

    use stats_sfc, only: &
      nvarmax_sfc, & ! Constant(s)
      stats_init_sfc ! Procedure(s)

    use error_code, only: &
      clubb_at_least_debug_level ! Function

    use constants, only: &
      fstdout, fstderr, var_length ! Constants

    implicit none

    ! Input Variables

    integer, intent(in) :: iunit  ! File unit for fnamelist

    character(len=*), intent(in) ::  & 
      fname_prefix, & ! Start of the stats filenames
      fdir            ! Directory to output to

    logical, intent(in) :: l_stats_in ! Stats on? T/F

    character(len=*), intent(in) :: &
      stats_fmt_in    ! Format of the stats file output

    real(kind=time_precision), intent(in) ::  & 
      stats_tsamp_in,  & ! Sampling interval   [s]
      stats_tout_in      ! Output interval     [s]

    character(len=*), intent(in) :: &
      fnamelist          ! Filename holding the &statsnl

    integer, intent(in) :: nnzp ! Grid points in the vertical [count]

    real, intent(in), dimension(nnzp) ::  & 
      gzt, gzm  ! Thermodynamic and momentum levels           [m]

    integer, intent(in) :: nnrad_zt ! Grid points in the radiation grid [count]

    real, intent(in), dimension(nnrad_zt) :: grad_zt ! Radiation levels [m]  

    integer, intent(in) :: nnrad_zm ! Grid points in the radiation grid [count]

    real, intent(in), dimension(nnrad_zm) :: grad_zm ! Radiation levels [m]

    integer, intent(in) :: day, month, year  ! Time of year

    real, dimension(1), intent(in) ::  & 
      rlat, rlon   ! Latitude and Longitude             [Degrees N/E]

    real(kind=time_precision), intent(in) ::  & 
      time_current ! Model time                         [s]

    real(kind=time_precision), intent(in) ::  & 
      delt         ! Timestep (dtmain in CLUBB)         [s]


    ! Local Variables

    ! Namelist Variables

    character(len=10) :: stats_fmt  ! File storage convention

    character(len=var_length), dimension(nvarmax_zt) ::  & 
      vars_zt  ! Variables on the thermodynamic levels

    character(len=var_length), dimension(nvarmax_zm) ::  & 
      vars_zm  ! Variables on the momentum levels

    character(len=var_length), dimension(nvarmax_rad_zt) ::  & 
      vars_rad_zt  ! Variables on the radiation levels

    character(len=var_length), dimension(nvarmax_rad_zm) ::  & 
      vars_rad_zm  ! Variables on the radiation levels

    character(len=var_length), dimension(nvarmax_sfc) ::  &
      vars_sfc ! Variables at the model surface

    namelist /statsnl/ & 
      vars_zt, & 
      vars_zm, &
      vars_rad_zt, &
      vars_rad_zm, & 
      vars_sfc

    ! Local Variables

    logical :: l_error

    character(len=200) :: fname

    integer :: i, ntot, read_status

    ! Initialize
    l_error = .false.

    ! Set stats_variables variables with inputs from calling subroutine
    l_stats = l_stats_in

    stats_tsamp = stats_tsamp_in
    stats_tsamp = stats_tsamp_in
    stats_tout  = stats_tout_in
    stats_fmt   = trim( stats_fmt_in )

    if ( .not. l_stats ) then
      l_stats_samp  = .false.
      l_stats_first = .false.
      l_stats_last  = .false.
      return
    end if

    ! Initialize namelist variables

    vars_zt  = ''
    vars_zm  = ''
    vars_rad_zt = ''
    vars_rad_zm = ''
    vars_sfc = ''

    ! Read namelist

    open(unit=iunit, file=fnamelist)
    read(unit=iunit, nml=statsnl, iostat=read_status, end=100)
    if ( read_status > 0 ) then
      write(fstderr,*) "Error reading stats namelist in file ",  &
                       trim( fnamelist )
      write(fstderr,*) "One cause is having more statistical variables ",  &
                       "listed in the namelist for var_zt, var_zm, or ",  &
                       "var_sfc than allowed by nvarmax_zt, nvarmax_zm, ",  &
                       "or nvarmax_sfc, respectively."
      write(fstderr,*) "Maximum variables allowed for var_zt = ", nvarmax_zt
      write(fstderr,*) "Maximum variables allowed for var_zm = ", nvarmax_zm
      write(fstderr,*) "Maximum variables allowed for var_rad_zt = ", nvarmax_rad_zt
      write(fstderr,*) "Maximum variables allowed for var_rad_zm = ", nvarmax_rad_zm
      write(fstderr,*) "Maximum variables allowed for var_sfc = ", nvarmax_sfc
      stop "stats_init:  error reading stats namelist."
    endif
    close(unit=iunit)

    if ( clubb_at_least_debug_level( 1 ) ) then
      write(fstdout,*) "--------------------------------------------------"

      write(fstdout,*) "Statistics"

      write(fstdout,*) "--------------------------------------------------"
      write(fstdout,*) "vars_zt = "
      i = 1
      do while ( vars_zt(i) /= '' )
        write(fstdout,*) vars_zt(i)
        i = i + 1
      end do

      write(fstdout,*) "vars_zm = "
      i = 1
      do while ( vars_zm(i) /= '' )
        write(fstdout,*) vars_zm(i)
        i = i + 1
      end do

      if (l_output_rad_files) then
        write(fstdout,*) "vars_rad_zt = "
        i = 1
        do while ( vars_rad_zt(i) /= '' )
          write(fstdout,*) vars_rad_zt(i)
          i = i + 1
        end do

        write(fstdout,*) "vars_rad_zm = "
        i = 1
        do while ( vars_rad_zm(i) /= '' )
          write(fstdout,*) vars_rad_zm(i)
          i = i + 1
        end do
      end if ! l_output_rad_files

      write(fstdout,*) "vars_sfc = "
      i = 1
      do while ( vars_sfc(i) /= '' )
        write(fstdout,*) vars_sfc(i)
        i = i + 1
      end do

      write(fstdout,*) "--------------------------------------------------"
    end if ! clubb_at_least_debug_level 1

    ! Determine file names for GrADS or NetCDF files
    fname_zt  = trim( fname_prefix )//"_zt"
    fname_zm  = trim( fname_prefix )//"_zm"
    fname_rad_zt  = trim( fname_prefix )//"_rad_zt"
    fname_rad_zm  = trim( fname_prefix )//"_rad_zm"
    fname_sfc = trim( fname_prefix )//"_sfc"

    ! Parse the file type for stats output.  Currently only GrADS and
    ! NetCDF v3 are supported by this code.

    select case( trim( stats_fmt ) )
    case( "GrADS", "grads", "gr" )
      l_netcdf = .false.
      l_grads  = .true.

    case ( "NetCDF", "netcdf", "nc" )
      l_netcdf = .true.
      l_grads  = .false.

    case default
      write(fstderr,*) "Invalid data format "//trim( stats_fmt )
      stop

    end select

    ! Check sampling and output frequencies

    ! The model time step length, delt (which is dtmain), should multiply
    ! evenly into the statistical sampling time step length, stats_tsamp.
    if ( abs( stats_tsamp/delt - floor(stats_tsamp/delt) )  & 
           > 1.e-8 ) then
      l_error = .true.  ! This will cause the run to stop.
      write(fstderr,*) 'Error:  stats_tsamp should be an even multiple of ',  &
                       'delt (which is dtmain).  Check the appropriate ',  &
                       'model.in file.'
      write(fstderr,*) 'stats_tsamp = ', stats_tsamp
      write(fstderr,*) 'delt = ', delt
    endif

    ! The statistical sampling time step length, stats_tsamp, should multiply
    ! evenly into the statistical output time step length, stats_tout.
    if ( abs( stats_tout/stats_tsamp - floor(stats_tout/stats_tsamp) ) & 
           > 1.e-8 ) then
      l_error = .true.  ! This will cause the run to stop.
      write(fstderr,*) 'Error:  stats_tout should be an even multiple of ',  &
                       'stats_tsamp.  Check the appropriate model.in file.'
      write(fstderr,*) 'stats_tout = ', stats_tout
      write(fstderr,*) 'stats_tsamp = ', stats_tsamp
    endif

    ! Initialize zt (mass points)

    i = 1
    do while ( ichar(vars_zt(i)(1:1)) /= 0  & 
               .and. len_trim(vars_zt(i)) /= 0 & 
               .and. i <= nvarmax_zt )
      i = i + 1
    enddo
    ntot = i - 1
    if ( ntot == nvarmax_zt ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "vars_zt than allowed for by nvarmax_zt."
      write(fstderr,*) "Check the number of variables listed for vars_zt ",  &
                       "in the stats namelist, or change nvarmax_zt."
      write(fstderr,*) "nvarmax_zt = ", nvarmax_zt
      stop "stats_init:  number of zt statistical variables exceeds limit"
    endif

    zt%nn = ntot
    zt%kk = nnzp

    allocate( zt%z( zt%kk ) )
    zt%z = gzt

    allocate( zt%x( 1, 1, zt%kk, zt%nn ) )
    allocate( zt%n( 1, 1, zt%kk, zt%nn ) )
    allocate( zt%l_in_update( 1, 1, zt%kk, zt%nn ) )
    call stats_zero( zt%kk, zt%nn, zt%x, zt%n, zt%l_in_update )

    allocate( zt%f%var( zt%nn ) )
    allocate( zt%f%z( zt%kk ) )

    ! Allocate scratch space

    allocate( ztscr01(zt%kk) )
    allocate( ztscr02(zt%kk) )
    allocate( ztscr03(zt%kk) )
    allocate( ztscr04(zt%kk) )
    allocate( ztscr05(zt%kk) )
    allocate( ztscr06(zt%kk) )
    allocate( ztscr07(zt%kk) )
    allocate( ztscr08(zt%kk) )
    allocate( ztscr09(zt%kk) )
    allocate( ztscr10(zt%kk) )
    allocate( ztscr11(zt%kk) )
    allocate( ztscr12(zt%kk) )
    allocate( ztscr13(zt%kk) )
    allocate( ztscr14(zt%kk) )
    allocate( ztscr15(zt%kk) )
    allocate( ztscr16(zt%kk) )
    allocate( ztscr17(zt%kk) )
    allocate( ztscr18(zt%kk) )
    allocate( ztscr19(zt%kk) )
    allocate( ztscr20(zt%kk) )
    allocate( ztscr21(zt%kk) )

    ztscr01 = 0.0
    ztscr02 = 0.0
    ztscr03 = 0.0
    ztscr04 = 0.0
    ztscr05 = 0.0
    ztscr06 = 0.0
    ztscr07 = 0.0
    ztscr08 = 0.0
    ztscr09 = 0.0
    ztscr10 = 0.0
    ztscr11 = 0.0
    ztscr12 = 0.0
    ztscr13 = 0.0
    ztscr14 = 0.0
    ztscr15 = 0.0
    ztscr16 = 0.0
    ztscr17 = 0.0
    ztscr18 = 0.0
    ztscr19 = 0.0
    ztscr20 = 0.0
    ztscr21 = 0.0

    fname = trim( fname_zt )

    if ( l_grads ) then

      ! Open GrADS file
      call open_grads( iunit, fdir, fname,  & 
                       1, zt%kk, zt%z, & 
                       day, month, year, rlat, rlon, & 
                       time_current+stats_tout, stats_tout, & 
                       zt%nn, zt%f )

    else ! Open NetCDF file
#ifdef NETCDF
      call open_netcdf( 1, 1, fdir, fname, 1, zt%kk, zt%z, &  ! In
                        day, month, year, rlat, rlon, &  ! In
                        time_current+stats_tout, stats_tout, zt%nn, &  ! In
                        zt%f ) ! InOut
#else
      stop "netCDF support was not compiled into this build."
#endif

    end if

    ! Default initialization for array indices for zt

    call stats_init_zt( vars_zt, l_error )

    ! Initialize zm (momentum points)

    i = 1
    do while ( ichar(vars_zm(i)(1:1)) /= 0  & 
               .and. len_trim(vars_zm(i)) /= 0 & 
               .and. i <= nvarmax_zm )
      i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_zm ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "vars_zm than allowed for by nvarmax_zm."
      write(fstderr,*) "Check the number of variables listed for vars_zm ",  &
                       "in the stats namelist, or change nvarmax_zm."
      write(fstderr,*) "nvarmax_zm = ", nvarmax_zm
      stop "stats_init:  number of zm statistical variables exceeds limit"
    endif

    zm%nn = ntot
    zm%kk = nnzp

    allocate( zm%z( zm%kk ) )
    zm%z = gzm

    allocate( zm%x( 1, 1, zm%kk, zm%nn ) )
    allocate( zm%n( 1, 1, zm%kk, zm%nn ) )
    allocate( zm%l_in_update( 1, 1, zm%kk, zm%nn ) )

    call stats_zero( zm%kk, zm%nn, zm%x, zm%n, zm%l_in_update )

    allocate( zm%f%var( zm%nn ) )
    allocate( zm%f%z( zm%kk ) )

    ! Allocate scratch space

    allocate( zmscr01(zm%kk) )
    allocate( zmscr02(zm%kk) )
    allocate( zmscr03(zm%kk) )
    allocate( zmscr04(zm%kk) )
    allocate( zmscr05(zm%kk) )
    allocate( zmscr06(zm%kk) )
    allocate( zmscr07(zm%kk) )
    allocate( zmscr08(zm%kk) )
    allocate( zmscr09(zm%kk) )
    allocate( zmscr10(zm%kk) )
    allocate( zmscr11(zm%kk) )
    allocate( zmscr12(zm%kk) )
    allocate( zmscr13(zm%kk) )
    allocate( zmscr14(zm%kk) )
    allocate( zmscr15(zm%kk) )
    allocate( zmscr16(zm%kk) )
    allocate( zmscr17(zm%kk) )

    zmscr01 = 0.0
    zmscr02 = 0.0
    zmscr03 = 0.0
    zmscr04 = 0.0
    zmscr05 = 0.0
    zmscr06 = 0.0
    zmscr07 = 0.0
    zmscr08 = 0.0
    zmscr09 = 0.0
    zmscr10 = 0.0
    zmscr11 = 0.0
    zmscr12 = 0.0
    zmscr13 = 0.0
    zmscr14 = 0.0
    zmscr15 = 0.0
    zmscr16 = 0.0
    zmscr17 = 0.0


    fname = trim( fname_zm )
    if ( l_grads ) then

      ! Open GrADS files
      call open_grads( iunit, fdir, fname,  & 
                       1, zm%kk, zm%z, & 
                       day, month, year, rlat, rlon, & 
                       time_current+stats_tout, stats_tout, & 
                       zm%nn, zm%f )

    else ! Open NetCDF file
#ifdef NETCDF
      call open_netcdf( 1, 1, fdir, fname, 1, zm%kk, zm%z, &  ! In
                        day, month, year, rlat, rlon, &  ! In
                        time_current+stats_tout, stats_tout, zm%nn, &  ! In
                        zm%f ) ! InOut

#else
      stop "netCDF support was not compiled into this build."
#endif
    end if

    call stats_init_zm( vars_zm, l_error )

    ! Initialize rad_zt (radiation points)

    if (l_output_rad_files) then
    
      i = 1
      do while ( ichar(vars_rad_zt(i)(1:1)) /= 0  & 
                 .and. len_trim(vars_rad_zt(i)) /= 0 & 
                 .and. i <= nvarmax_rad_zt )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_rad_zt ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_rad_zt than allowed for by nvarmax_rad_zt."
        write(fstderr,*) "Check the number of variables listed for vars_rad_zt ",  &
                         "in the stats namelist, or change nvarmax_rad_zt."
        write(fstderr,*) "nvarmax_rad_zt = ", nvarmax_rad_zt
        stop "stats_init:  number of rad_zt statistical variables exceeds limit"
      endif

      rad_zt%nn = ntot
      rad_zt%kk = nnrad_zt

      allocate( rad_zt%z( rad_zt%kk ) )
      rad_zt%z = grad_zt

      allocate( rad_zt%x( 1, 1, rad_zt%kk, rad_zt%nn ) )
      allocate( rad_zt%n( 1, 1, rad_zt%kk, rad_zt%nn ) )
      allocate( rad_zt%l_in_update( 1, 1, rad_zt%kk, rad_zt%nn ) )

      call stats_zero( rad_zt%kk, rad_zt%nn, rad_zt%x, rad_zt%n, rad_zt%l_in_update )

      allocate( rad_zt%f%var( rad_zt%nn ) )
      allocate( rad_zt%f%z( rad_zt%kk ) )

      ! Allocate scratch space

      !allocate( radscr01(rad%kk) )
      !allocate( radscr02(rad%kk) )
      !allocate( radscr03(rad%kk) )
      !allocate( radscr04(rad%kk) )
      !allocate( radscr05(rad%kk) )
      !allocate( radscr06(rad%kk) )
      !allocate( radscr07(rad%kk) )
      !allocate( radscr08(rad%kk) )
      !allocate( radscr09(rad%kk) )
      !allocate( radscr10(rad%kk) )
      !allocate( radscr11(rad%kk) )
      !allocate( radscr12(rad%kk) )
      !allocate( radscr13(rad%kk) )
      !allocate( radscr14(rad%kk) )
      !allocate( radscr15(rad%kk) )
      !allocate( radscr16(rad%kk) )
      !allocate( radscr17(rad%kk) )

      !radscr01 = 0.0
      !radscr02 = 0.0
      !radscr03 = 0.0
      !radscr04 = 0.0
      !radscr05 = 0.0
      !radscr06 = 0.0
      !radscr07 = 0.0
      !radscr08 = 0.0
      !radscr09 = 0.0
      !radscr10 = 0.0
      !radscr11 = 0.0
      !radscr12 = 0.0
      !radscr13 = 0.0
      !radscr14 = 0.0
      !radscr15 = 0.0
      !radscr16 = 0.0
      !radscr17 = 0.0


      fname = trim( fname_rad_zt )
      if ( l_grads ) then

        ! Open GrADS files
        call open_grads( iunit, fdir, fname,  & 
                         1, rad_zt%kk, rad_zt%z, & 
                         day, month, year, rlat, rlon, & 
                         time_current+stats_tout, stats_tout, & 
                         rad_zt%nn, rad_zt%f )

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf( 1, 1, fdir, fname,  & 
                          1, rad_zt%kk, rad_zt%z, & 
                          day, month, year, rlat, rlon, & 
                          time_current+stats_tout, stats_tout, & 
                          rad_zt%nn, rad_zt%f )

#else
        stop "netCDF support was not compiled into this build."
#endif
      end if

      call stats_init_rad_zt( vars_rad_zt, l_error )

      ! Initialize rad_zm (radiation points)

      i = 1
      do while ( ichar(vars_rad_zm(i)(1:1)) /= 0  & 
                 .and. len_trim(vars_rad_zm(i)) /= 0 & 
                 .and. i <= nvarmax_rad_zm )
        i = i + 1
      end do
      ntot = i - 1
      if ( ntot == nvarmax_rad_zm ) then
        write(fstderr,*) "There are more statistical variables listed in ",  &
                         "vars_rad_zm than allowed for by nvarmax_rad_zm."
        write(fstderr,*) "Check the number of variables listed for vars_rad_zm ",  &
                         "in the stats namelist, or change nvarmax_rad_zm."
        write(fstderr,*) "nvarmax_rad_zm = ", nvarmax_rad_zm
        stop "stats_init:  number of rad_zm statistical variables exceeds limit"
      endif

      rad_zm%nn = ntot
      rad_zm%kk = nnrad_zm

      allocate( rad_zm%z( rad_zm%kk ) )
      rad_zm%z = grad_zm

      allocate( rad_zm%x( 1, 1, rad_zm%kk, rad_zm%nn ) )
      allocate( rad_zm%n( 1, 1, rad_zm%kk, rad_zm%nn ) )
      allocate( rad_zm%l_in_update( 1, 1, rad_zm%kk, rad_zm%nn ) )

      call stats_zero( rad_zm%kk, rad_zm%nn, rad_zm%x, rad_zm%n, rad_zm%l_in_update )

      allocate( rad_zm%f%var( rad_zm%nn ) )
      allocate( rad_zm%f%z( rad_zm%kk ) )

      ! Allocate scratch space

      !allocate( radscr01(rad%kk) )
      !allocate( radscr02(rad%kk) )
      !allocate( radscr03(rad%kk) )
      !allocate( radscr04(rad%kk) )
      !allocate( radscr05(rad%kk) )
      !allocate( radscr06(rad%kk) )
      !allocate( radscr07(rad%kk) )
      !allocate( radscr08(rad%kk) )
      !allocate( radscr09(rad%kk) )
      !allocate( radscr10(rad%kk) )
      !allocate( radscr11(rad%kk) )
      !allocate( radscr12(rad%kk) )
      !allocate( radscr13(rad%kk) )
      !allocate( radscr14(rad%kk) )
      !allocate( radscr15(rad%kk) )
      !allocate( radscr16(rad%kk) )
      !allocate( radscr17(rad%kk) )

      !radscr01 = 0.0
      !radscr02 = 0.0
      !radscr03 = 0.0
      !radscr04 = 0.0
      !radscr05 = 0.0
      !radscr06 = 0.0
      !radscr07 = 0.0
      !radscr08 = 0.0
      !radscr09 = 0.0
      !radscr10 = 0.0
      !radscr11 = 0.0
      !radscr12 = 0.0
      !radscr13 = 0.0
      !radscr14 = 0.0
      !radscr15 = 0.0
      !radscr16 = 0.0
      !radscr17 = 0.0


      fname = trim( fname_rad_zm )
      if ( l_grads ) then

        ! Open GrADS files
        call open_grads( iunit, fdir, fname,  & 
                         1, rad_zm%kk, rad_zm%z, & 
                         day, month, year, rlat, rlon, & 
                         time_current+stats_tout, stats_tout, & 
                         rad_zm%nn, rad_zm%f )

      else ! Open NetCDF file
#ifdef NETCDF
        call open_netcdf( 1, 1, fdir, fname,  & 
                          1, rad_zm%kk, rad_zm%z, & 
                          day, month, year, rlat, rlon, & 
                          time_current+stats_tout, stats_tout, & 
                          rad_zm%nn, rad_zm%f )

#else
        stop "netCDF support was not compiled into this build."
#endif
      end if

      call stats_init_rad_zm( vars_rad_zm, l_error )
    end if ! l_output_rad_files


    ! Initialize sfc (surface point)

    i = 1
    do while ( ichar(vars_sfc(i)(1:1)) /= 0  & 
               .and. len_trim(vars_sfc(i)) /= 0 & 
               .and. i <= nvarmax_sfc )
      i = i + 1
    end do
    ntot = i - 1
    if ( ntot == nvarmax_sfc ) then
      write(fstderr,*) "There are more statistical variables listed in ",  &
                       "vars_sfc than allowed for by nvarmax_sfc."
      write(fstderr,*) "Check the number of variables listed for vars_sfc ",  &
                       "in the stats namelist, or change nvarmax_sfc."
      write(fstderr,*) "nvarmax_sfc = ", nvarmax_sfc
      stop "stats_init:  number of sfc statistical variables exceeds limit"
    endif

    sfc%nn = ntot
    sfc%kk = 1

    allocate( sfc%z( sfc%kk ) )
    sfc%z = gzm(1)

    allocate( sfc%x( 1, 1, sfc%kk, sfc%nn ) )
    allocate( sfc%n( 1, 1, sfc%kk, sfc%nn ) )
    allocate( sfc%l_in_update( 1, 1, sfc%kk, sfc%nn ) )

    call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n, sfc%l_in_update )

    allocate( sfc%f%var( sfc%nn ) )
    allocate( sfc%f%z( sfc%kk ) )

    fname = trim( fname_sfc )

    if ( l_grads ) then

      ! Open GrADS files
      call open_grads( iunit, fdir, fname,  & 
                       1, sfc%kk, sfc%z, & 
                       day, month, year, rlat, rlon, & 
                         time_current+stats_tout, stats_tout, & 
                         sfc%nn, sfc%f )

    else ! Open NetCDF files
#ifdef NETCDF
      call open_netcdf( 1, 1, fdir, fname, 1, sfc%kk, sfc%z, &  ! In
                        day, month, year, rlat, rlon, &  ! In
                        time_current+stats_tout, stats_tout, sfc%nn, &  ! In
                        sfc%f ) ! InOut

#else
      stop "netCDF support was not compiled into this build."
#endif
    end if

    call stats_init_sfc( vars_sfc, l_error )

    ! Check for errors

    if ( l_error ) then
      write(fstderr,*) 'stats_init:  errors found'
      stop
    endif

    return

    ! If namelist was not found in input file, turn off statistics

    100 continue
    write(fstderr,*) 'Error with statsnl, statistics is turned off'
    l_stats       = .false.
    l_stats_samp  = .false.
    l_stats_first = .false.
    l_stats_last  = .false.

    return
  end subroutine stats_init
  !-----------------------------------------------------------------------
  subroutine stats_zero( kk, nn, x, n, l_in_update )

    !     Description:
    !     Initialize stats to zero
    !-----------------------------------------------------------------------
    use stats_precision, only: & 
        stat_rknd,   & ! Variable(s)
        stat_nknd

    implicit none

    ! Input
    integer, intent(in) :: kk, nn

    ! Output
    real(kind=stat_rknd), dimension(1,1,kk,nn), intent(out)    :: x
    integer(kind=stat_nknd), dimension(1,1,kk,nn), intent(out) :: n
    logical, dimension(1,1,kk,nn), intent(out) :: l_in_update

    ! Zero out arrays

    if ( nn > 0 ) then
      x(:,:,:,:) = 0.0
      n(:,:,:,:) = 0
      l_in_update(:,:,:,:) = .false.
    end if

    return
  end subroutine stats_zero

  !-----------------------------------------------------------------------
  subroutine stats_avg( kk, nn, x, n )

    !     Description:
    !     Compute the average of stats fields
    !-----------------------------------------------------------------------
    use stats_precision, only: & 
        stat_rknd,   & ! Variable(s)
        stat_nknd

    implicit none

    ! Input
    integer, intent(in) :: nn, kk
    integer(kind=stat_nknd), dimension(1,1,kk,nn), intent(in) :: n

    ! Output
    real(kind=stat_rknd), dimension(1,1,kk,nn), intent(inout)  :: x

    ! Internal

    integer k,m

    ! Compute averages

    do m=1,nn
      do k=1,kk

        if ( n(1,1,k,m) > 0 ) then
          x(1,1,k,m) = x(1,1,k,m) / real( n(1,1,k,m) )
        end if

      end do
    end do

    return
  end subroutine stats_avg

  !-----------------------------------------------------------------------
  subroutine stats_begin_timestep( time_elapsed, delt )

    !     Description:
    !     Begin sampling for the current timestep.
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
        l_stats,  & ! Variable(s)
        l_stats_samp, & 
        l_stats_first, & 
        l_stats_last, & 
        stats_tsamp, & 
        stats_tout

    use stats_precision, only: & 
        time_precision ! Variable(s)

    implicit none

    ! Input

    real(kind=time_precision), intent(in) ::  & 
      time_elapsed ! Elapsed model time       [s]

    real(kind=time_precision), intent(in) ::  & 
      delt         ! Model time step          [s]

    if ( .not. l_stats ) return

    ! Set sample this time step flag
    if ( mod( time_elapsed, stats_tsamp ) < 1.e-8 ) then
      l_stats_samp = .true.
    else
      l_stats_samp = .false.
    end if

    ! Set first time step flag

    if ( mod( time_elapsed - delt, stats_tout ) < 1.e-8 ) then
      l_stats_first = .true.
    else
      l_stats_first = .false.
    end if

    ! Set last time step flag

    if ( mod( time_elapsed, stats_tout ) < 1.e-8 ) then
      l_stats_last = .true.
    else
      l_stats_last = .false.
    end if

    return

  end subroutine stats_begin_timestep

  !-----------------------------------------------------------------------
  subroutine stats_end_timestep( )

    !     Description: Called when the stats timestep has ended. This subroutine
    !     is responsible for calling statistics to be written to the output
    !     format.
    !-----------------------------------------------------------------------

    use constants, only: &
        fstderr ! Constant(s)

    use stats_variables, only: & 
        zt,  & ! Variable(s)
        zm, & 
        rad_zt, &
        rad_zm, &
        sfc, & 
        l_stats_last, & 
        stats_tsamp, & 
        stats_tout, &
        l_output_rad_files, & 
        l_grads

    use stats_precision, only: & 
        time_precision ! Variable(s)

    use output_grads, only: & 
        write_grads ! Procedure(s)

    use error_code, only: &
        clubb_at_least_debug_level ! Procedure(s)

#ifdef NETCDF
    use output_netcdf, only: & 
        write_netcdf ! Procedure(s)
#endif

    implicit none

    ! Local Variables

    integer :: i, k

    logical :: l_error

    ! Check if it is time to write to file

    if ( .not. l_stats_last ) return

    ! Initialize
    l_error = .false.

    ! Check number of sampling points for each variable in the zt statistics
    ! at each vertical level.
    do i = 1, zt%nn
      do k = 1, zt%kk

        if ( zt%n(1,1,k,i) /= 0 .and.  &
             zt%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            ! Made error message more descriptive
            ! Joshua Fasching July 2008
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(zt%f%var(i)%name), ' in zt ',  &
                             'at k = ', k,  &
                             '; zt%n(',k,',',i,') = ', zt%n(1,1,k,i)
          endif

        endif

      enddo
    enddo

    ! Check number of sampling points for each variable in the zm statistics
    ! at each vertical level.
    do i = 1, zm%nn
      do k = 1, zm%kk

        if ( zm%n(1,1,k,i) /= 0 .and.  &
             zm%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            ! Made error message more descriptive
            ! Joshua Fasching July 2008
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(zm%f%var(i)%name), ' in zm ',  &
                             'at k = ', k,  &
                             '; zm%n(',k,',',i,') = ', zm%n(1,1,k,i)
          endif

        endif

      enddo
    enddo

    if (l_output_rad_files) then
      ! Check number of sampling points for each variable in the rad statistics
      ! at each vertical level.
      do i = 1, rad_zt%nn
        do k = 1, rad_zt%kk

          if ( rad_zt%n(1,1,k,i) /= 0 .and.  &
               rad_zt%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

            l_error = .true.  ! This will stop the run

            if ( clubb_at_least_debug_level( 1 ) ) then
              ! Made error message more descriptive
              ! Joshua Fasching July 2008
              write(fstderr,*) 'Possible sampling error for variable ',  &
                               trim(rad_zt%f%var(i)%name), ' in rad_zt ',  &
                               'at k = ', k,  &
                               '; rad_zt%n(',k,',',i,') = ', rad_zt%n(1,1,k,i)
            endif

          endif

        enddo
      enddo
    
      ! Check number of sampling points for each variable in the zm statistics
      ! at each vertical level.
      do i = 1, rad_zm%nn
        do k = 1, rad_zm%kk

          if ( rad_zm%n(1,1,k,i) /= 0 .and.  &
               rad_zm%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

            l_error = .true.  ! This will stop the run

            if ( clubb_at_least_debug_level( 1 ) ) then
              ! Made error message more descriptive
              ! Joshua Fasching July 2008
              write(fstderr,*) 'Possible sampling error for variable ',  &
                               trim(rad_zm%f%var(i)%name), ' in rad_zm ',  &
                               'at k = ', k,  &
                               '; rad_zm%n(',k,',',i,') = ', rad_zm%n(1,1,k,i)
            endif

          endif

        enddo
      enddo
    end if ! l_output_rad_files

    ! Check number of sampling points for each variable in the zm statistics
    ! at each vertical level.
    do i = 1, sfc%nn
      do k = 1, sfc%kk

        if ( sfc%n(1,1,k,i) /= 0 .and.  &
             sfc%n(1,1,k,i) /= floor(stats_tout/stats_tsamp) ) then

          l_error = .true.  ! This will stop the run

          if ( clubb_at_least_debug_level( 1 ) ) then
            ! Made error message more descriptive
            ! Joshua Fasching July 2008
            write(fstderr,*) 'Possible sampling error for variable ',  &
                             trim(sfc%f%var(i)%name), ' in sfc ',  &
                             'at k = ', k,  &
                             '; sfc%n(',k,',',i,') = ', sfc%n(1,1,k,i)
          endif

        endif

      enddo
    enddo

    ! Stop the run if errors are found.
    if ( l_error ) then
      write(fstderr,*) 'Possible statistical sampling error'
      write(fstderr,*) 'For details, set debug_level to a value of at ',  &
                       'least 1 in the appropriate model.in file.'
      stop 'stats_end_timestep:  error(s) found'
    endif

    ! Compute averages

    call stats_avg( zt%kk, zt%nn, zt%x, zt%n )
    call stats_avg( zm%kk, zm%nn, zm%x, zm%n )
    if (l_output_rad_files) then
      call stats_avg( rad_zt%kk, rad_zt%nn, rad_zt%x, rad_zt%n )
      call stats_avg( rad_zm%kk, rad_zm%nn, rad_zm%x, rad_zm%n )
    end if
    call stats_avg( sfc%kk, sfc%nn, sfc%x, sfc%n )

    ! Write to file
    if ( l_grads ) then
      call write_grads( zt%f  )
      call write_grads( zm%f  )
      if (l_output_rad_files) then
        call write_grads( rad_zt%f  )
        call write_grads( rad_zm%f  )
      end if
      call write_grads( sfc%f  )
    else ! l_netcdf
#ifdef NETCDF
      call write_netcdf( zt%f  )
      call write_netcdf( zm%f  )
      if (l_output_rad_files) then
        call write_netcdf( rad_zt%f  )
        call write_netcdf( rad_zm%f  )
      end if  
      call write_netcdf( sfc%f  )
#else
      stop "This program was not compiled with netCDF support"
#endif
    endif

    ! Reset sample fields
    call stats_zero( zt%kk, zt%nn, zt%x, zt%n, zt%l_in_update )
    call stats_zero( zm%kk, zm%nn, zm%x, zm%n, zm%l_in_update )
    if (l_output_rad_files) then
      call stats_zero( rad_zt%kk, rad_zt%nn, rad_zt%x, rad_zt%n, rad_zt%l_in_update )
      call stats_zero( rad_zm%kk, rad_zm%nn, rad_zm%x, rad_zm%n, rad_zm%l_in_update )
    end if
    call stats_zero( sfc%kk, sfc%nn, sfc%x, sfc%n, sfc%l_in_update )


    return
  end subroutine stats_end_timestep

  !----------------------------------------------------------------------
  subroutine stats_accumulate & 
                   ( um, vm, upwp, vpwp, up2, vp2, thlm, &
                     rtm, wprtp, wpthlp, wpthvp, wprcp, &
                     wp2, wp3, rtp2, thlp2, rtpthlp, &
                     p_in_Pa, exner, rho, rho_zm, rho_ds_zm, &
                     rho_ds_zt, Kh_zt, wm_zt, wm_zm, &
                     sigma_sqd_w, tau_zm, rcm, cloud_frac, &
                     rcm_in_layer, cloud_cover, &
                     pdf_params, &
                     sclrm, sclrp2, sclrprtp, sclrpthlp, sclrm_forcing, &
                     wpsclrp, edsclrm, edsclrm_forcing )

    ! Description:
    ! Accumulate those stats variables that are preserved in CLUBB from timestep to
    ! timestep, but not those stats that are not, (e.g. budget terms, longwave and
    ! shortwave components, etc. )
    !----------------------------------------------------------------------

    use stats_variables, only: & 
        zt,      & ! Variables
        zm, & 
        sfc, & 
        l_stats_samp, & 
        ithlm, & 
        iT_in_K, & 
        ithvm, & 
        irtm, & 
        ircm, & 
        ium, & 
        ivm, & 
        iwm_zt, & 
        iwm_zm, & 
        iug, & 
        ivg, & 
        icloud_frac, & 
        ircm_in_layer, &
        icloud_cover, &
        ip_in_Pa, & 
        iexner, & 
        irho_ds_zt, &
        iLscale, & 
        iwp3, & 
        iwpthlp2, & 
        iwp2thlp,  & 
        iwprtp2, & 
        iwp2rtp, & 
        iLscale_up, & 
        iLscale_down, & 
        itau_zt, & 
        iKh_zt, & 
        iwp2thvp, & 
        iwp2rcp, & 
        iwprtpthlp, & 
        isigma_sqd_w_zt,          & 
        irho, & 
        irsat, & 
        iAKm, & 
        iLH_AKm, & 
        iradht

    use stats_variables, only: & 
        ia, & 
        iw1, & 
        iw2, & 
        ivarnce_w1, & 
        ivarnce_w2, & 
        ithl1, & 
        ithl2, & 
        ivarnce_thl1, & 
        ivarnce_thl2, & 
        irt1, & 
        irt2, & 
        ivarnce_rt1, & 
        ivarnce_rt2, & 
        irc1, & 
        irc2, & 
        irsl1, & 
        irsl2, & 
        icloud_frac1, & 
        icloud_frac2, & 
        is1, & 
        is2, & 
        istdev_s1, & 
        istdev_s2, & 
        irrtthl

    use stats_variables, only: & 
        iwp2_zt, & 
        ithlp2_zt, & 
        iwpthlp_zt, & 
        iwprtp_zt, & 
        irtp2_zt, & 
        irtpthlp_zt, &
        iup2_zt, &
        ivp2_zt, &
        iupwp_zt, &
        ivpwp_zt, & 
        iwp2, & 
        irtp2, & 
        ithlp2, & 
        irtpthlp, & 
        iwprtp,  & 
        iwpthlp, & 
        iwp4,  & 
        iwpthvp, & 
        irtpthvp, & 
        ithlpthvp, & 
        itau_zm, & 
        iKh_zm, & 
        iwprcp, & 
        ithlprcp, & 
        irtprcp, & 
        ircp2, & 
        iupwp, & 
        ivpwp, & 
        iup2, & 
        ivp2, & 
        irho_zm, & 
        isigma_sqd_w, &
        irho_ds_zm, &
        iem

    use stats_variables, only: & 
        ishear, & 
        iFrad, & 
        icc, & 
        izb, & 
        ilwp, &
        ivwp, &
        iswp, &
        iiwp, &
        ithlm_vert_avg, &
        irtm_vert_avg, &
        ium_vert_avg, &
        ivm_vert_avg, &
        iwp2_vert_avg, &
        iup2_vert_avg, &
        ivp2_vert_avg, &
        irtp2_vert_avg, &
        ithlp2_vert_avg

    use stats_variables, only: & 
        isclrm, & 
        isclrm_f, & 
        iedsclrm, & 
        iedsclrm_f, & 
        isclrprtp, & 
        isclrp2, & 
        isclrpthvp, & 
        isclrpthlp, & 
        isclrprcp, & 
        iwpsclrp, & 
        iwp2sclrp, & 
        iwpsclrp2, & 
        iwpsclrprtp, & 
        iwpsclrpthlp, & 
        iwpedsclrp

    use stats_variables, only: &
      iAKstd, &
      iAKstd_cld, &
      iAKm_rcm, &
      iAKm_rcc

    use stats_variables, only: &
      iLH_rcm_avg

    use grid_class, only: & 
        gr ! Variable

    use variables_diagnostic_module, only: & 
        hydromet, &
        thvm, & ! Variable(s)
        ug, & 
        vg, & 
        Lscale, & 
        wpthlp2, & 
        wp2thlp, & 
        wprtp2, & 
        wp2rtp, & 
        Lscale_up, & 
        Lscale_down, & 
        tau_zt, & 
        wp2thvp, & 
        wp2rcp, & 
        wprtpthlp, & 
        sigma_sqd_w_zt, & 
        rsat, & 
        AKm, & 
        lh_AKm, & 
        radht

    use variables_diagnostic_module, only: & 
        wp2_zt, & 
        thlp2_zt, & 
        wpthlp_zt, & 
        wprtp_zt, & 
        rtp2_zt, & 
        rtpthlp_zt, &
        up2_zt, &
        vp2_zt, &
        upwp_zt, &
        vpwp_zt, & 
        wp4, & 
        rtpthvp, & 
        thlpthvp, & 
        Kh_zm, & 
        thlprcp, & 
        rtprcp, & 
        rcp2, & 
        em, & 
        shear, & 
        Frad, & 
        sclrpthvp, & 
        sclrprcp, & 
        wp2sclrp, & 
        wpsclrp2, & 
        wpsclrprtp, & 
        wpsclrpthlp, & 
        wpedsclrp

    use variables_diagnostic_module, only: & 
      AKstd, & ! Variable(s)
      lh_rcm_avg, &
      AKstd_cld, &
      AKm_rcm, &
      AKm_rcc

    use variables_prognostic_module, only: & 
      pdf_parameter ! Type

    use T_in_K_mod, only: & 
        thlm2T_in_K ! Procedure

    use constants, only: & 
        rc_tol

    use parameters_model, only: & 
        sclr_dim,  & ! Variable(s)
        edsclr_dim

    use stats_type, only: & 
        stat_update_var,  & ! Procedure(s)
        stat_update_var_pt

    use fill_holes, only: &
        vertical_avg

    use interpolation, only: & 
        lin_int ! Procedure

    use array_index, only: & 
        iirsnowm, iiricem

    implicit none

    ! Input Variable
    real, intent(in), dimension(gr%nnzp) :: & 
      um,      & ! u wind                        [m/s]
      vm,      & ! v wind                        [m/s]
      upwp,    & ! vertical u momentum flux      [m^2/s^2]
      vpwp,    & ! vertical v momentum flux      [m^2/s^2]
      up2,     & ! u'^2                          [m^2/s^2]
      vp2,     & ! v'^2                          [m^2/s^2]
      thlm,    & ! liquid potential temperature  [K]
      rtm,     & ! total water mixing ratio      [kg/kg]
      wprtp,   & ! w'rt'                         [(kg/kg) m/s]
      wpthlp,  & ! w'thl'                        [m K /s]
      wpthvp,  & ! w'thv'                        [m K /s]
      wprcp,   & ! w'rc'                         [(kg/kg) m/s]
      wp2,     & ! w'^2                          [m^2/s^2]
      wp3,     & ! w'^3                          [m^3/s^3]
      rtp2,    & ! rt'^2                         [(kg/kg)^2]
      thlp2,   & ! thl'^2                        [K^2]
      rtpthlp    ! rt'thl'                       [kg/kg K]

    real, intent(in), dimension(gr%nnzp) :: & 
      p_in_Pa,      & ! Pressure (Pa) on thermodynamic points    [Pa]
      exner,        & ! Exner function = ( p / p0 ) ** kappa     [-]
      rho,          & ! Density                                  [kg/m^3]
      rho_zm,       & ! Density                                  [kg/m^3]
      rho_ds_zm,    & ! Dry, static density (momentum levels)    [kg/m^3]
      rho_ds_zt,    & ! Dry, static density (thermo. levs.)      [kg/m^3]
      Kh_zt,        & ! Eddy diffusivity                         [m^2/s]
      wm_zt,        & ! w on thermodynamic levels                [m/s]
      wm_zm,        & ! w on momentum levels                     [m/s]
      sigma_sqd_w,  & ! PDF width paramter                       [-]
      tau_zm          ! Dissipation time                         [s]

    real, intent(in), dimension(gr%nnzp) :: & 
      rcm,         & ! Cloud water mixing ratio                 [kg/kg]
      cloud_frac,  & ! Cloud fraction                           [-]
      rcm_in_layer,& ! Cloud water mixing ratio in cloud layer  [kg/kg]
      cloud_cover    ! Cloud cover                              [-]

    type(pdf_parameter), intent(in) :: & 
      pdf_params ! PDF parameters [units vary]

    real, intent(in), dimension(gr%nnzp,sclr_dim) :: & 
      sclrm,           & ! High-order passive scalar          [units vary]
      sclrp2,          & ! High-order passive scalar variance [units^2]
      sclrprtp,        & ! High-order passive scalar covariance [units kg/kg]
      sclrpthlp,       & ! High-order passive scalar covariance [units K]
      sclrm_forcing,   & ! Large-scale forcing of scalar      [units/s]
      wpsclrp            ! w'sclr'                            [units m/s]

    real, intent(in), dimension(gr%nnzp,edsclr_dim) :: & 
      edsclrm,         & ! Eddy-diff passive scalar      [units vary] 
      edsclrm_forcing    ! Large-scale forcing of edscalar  [units vary]

    ! Local Variables

    integer :: i, k

    real :: xtmp

    real, dimension(gr%nnzp) ::  & 
      rsnowm,  & ! Snow mixing ratio                         [kg/kg]
      ricem      ! Prisitine ice water mixing ratio          [kg/kg]


    ! Sample fields

    if ( l_stats_samp ) then

      ! zt variables

      call stat_update_var( ithlm, thlm, zt )
      call stat_update_var( iT_in_K,  & 
                            thlm2T_in_K( thlm, exner, rcm), zt )
      call stat_update_var( ithvm, thvm, zt )
      call stat_update_var( irtm, rtm, zt )
      call stat_update_var( ircm, rcm, zt )
      call stat_update_var( ium, um, zt )
      call stat_update_var( ivm, vm, zt )
      call stat_update_var( iwm_zt, wm_zt, zt )
      call stat_update_var( iwm_zm, wm_zm, zm )
      call stat_update_var( iug, ug, zt )
      call stat_update_var( ivg, vg, zt )
      call stat_update_var( icloud_frac, cloud_frac, zt )
      call stat_update_var( ircm_in_layer, rcm_in_layer, zt )
      call stat_update_var( icloud_cover, cloud_cover, zt )
      call stat_update_var( ip_in_Pa, p_in_Pa, zt )
      call stat_update_var( iexner, exner, zt )
      call stat_update_var( irho_ds_zt, rho_ds_zt, zt )
      call stat_update_var( iLscale, Lscale, zt )
      call stat_update_var( iwp3, wp3, zt )
      call stat_update_var( iwpthlp2, wpthlp2, zt )
      call stat_update_var( iwp2thlp, wp2thlp, zt )
      call stat_update_var( iwprtp2, wprtp2, zt )
      call stat_update_var( iwp2rtp, wp2rtp, zt )
      call stat_update_var( iLscale_up, Lscale_up, zt )
      call stat_update_var( iLscale_down, Lscale_down, zt )
      call stat_update_var( itau_zt, tau_zt, zt )
      call stat_update_var( iKh_zt, Kh_zt, zt )
      call stat_update_var( iwp2thvp, wp2thvp, zt )
      call stat_update_var( iwp2rcp, wp2rcp, zt )
      call stat_update_var( iwprtpthlp, wprtpthlp, zt )
      call stat_update_var( isigma_sqd_w_zt, sigma_sqd_w_zt, zt )
      call stat_update_var( irho, rho, zt )
      call stat_update_var( irsat, rsat, zt )

      call stat_update_var( iAKm, AKm, zt )
      call stat_update_var( iLH_AKm, lh_AKm, zt)
      call stat_update_var( iLH_rcm_avg, lh_rcm_avg, zt )
      call stat_update_var( iAKstd, AKstd, zt )
      call stat_update_var( iAKstd_cld, AKstd_cld, zt )

      call stat_update_var( iAKm_rcm, AKm_rcm, zt)
      call stat_update_var( iAKm_rcc, AKm_rcc, zt )

      call stat_update_var( iradht, radht, zt )
      call stat_update_var( ia, pdf_params%a, zt )
      call stat_update_var( iw1, pdf_params%w1, zt )
      call stat_update_var( iw2, pdf_params%w2, zt )
      call stat_update_var( ivarnce_w1, pdf_params%varnce_w1, zt )
      call stat_update_var( ivarnce_w2, pdf_params%varnce_w2, zt )
      call stat_update_var( ithl1, pdf_params%thl1, zt )
      call stat_update_var( ithl2, pdf_params%thl2, zt )
      call stat_update_var( ivarnce_thl1, pdf_params%varnce_thl1, zt )
      call stat_update_var( ivarnce_thl2, pdf_params%varnce_thl2, zt )
      call stat_update_var( irt1, pdf_params%rt1, zt )
      call stat_update_var( irt2, pdf_params%rt2, zt )
      call stat_update_var( ivarnce_rt1, pdf_params%varnce_rt1, zt )
      call stat_update_var( ivarnce_rt2, pdf_params%varnce_rt2, zt )
      call stat_update_var( irc1, pdf_params%rc1, zt )
      call stat_update_var( irc2, pdf_params%rc2, zt )
      call stat_update_var( irsl1, pdf_params%rsl1, zt )
      call stat_update_var( irsl2, pdf_params%rsl2, zt )
      call stat_update_var( icloud_frac1, pdf_params%cloud_frac1, zt )
      call stat_update_var( icloud_frac2, pdf_params%cloud_frac2, zt )
      call stat_update_var( is1, pdf_params%s1, zt )
      call stat_update_var( is2, pdf_params%s2, zt )
      call stat_update_var( istdev_s1, pdf_params%stdev_s1, zt )
      call stat_update_var( istdev_s2, pdf_params%stdev_s2, zt )
      call stat_update_var( irrtthl, pdf_params%rrtthl, zt )
      call stat_update_var( iwp2_zt, wp2_zt, zt )
      call stat_update_var( ithlp2_zt, thlp2_zt, zt )
      call stat_update_var( iwpthlp_zt, wpthlp_zt, zt )
      call stat_update_var( iwprtp_zt, wprtp_zt, zt )
      call stat_update_var( irtp2_zt, rtp2_zt, zt )
      call stat_update_var( irtpthlp_zt, rtpthlp_zt, zt )
      call stat_update_var( iup2_zt, up2_zt, zt )
      call stat_update_var( ivp2_zt, vp2_zt, zt )
      call stat_update_var( iupwp_zt, upwp_zt, zt )
      call stat_update_var( ivpwp_zt, vpwp_zt, zt )

      if (sclr_dim > 0 ) then
        do i=1, sclr_dim
          call stat_update_var( isclrm(i), sclrm(:,i), zt )
          call stat_update_var( isclrm_f(i), sclrm_forcing(:,i),  zt )
        end do
      end if

      if ( edsclr_dim > 0 ) then
        do i=1, edsclr_dim
          call stat_update_var( iedsclrm(i), edsclrm(:,i), zt )
          call stat_update_var( iedsclrm_f(i), edsclrm_forcing(:,i), zt )
        end do
      end if

      ! zm variables

      call stat_update_var( iwp2, wp2, zm )
      call stat_update_var( irtp2, rtp2, zm )
      call stat_update_var( ithlp2, thlp2, zm )
      call stat_update_var( irtpthlp, rtpthlp, zm )
      call stat_update_var( iwprtp, wprtp, zm )
      call stat_update_var( iwpthlp, wpthlp, zm )
      call stat_update_var( iwp4, wp4, zm )
      call stat_update_var( iwpthvp, wpthvp, zm )
      call stat_update_var( irtpthvp, rtpthvp, zm )
      call stat_update_var( ithlpthvp, thlpthvp, zm )
      call stat_update_var( itau_zm, tau_zm, zm )
      call stat_update_var( iKh_zm, Kh_zm, zm )
      call stat_update_var( iwprcp, wprcp, zm )
      call stat_update_var( ithlprcp, thlprcp, zm )
      call stat_update_var( irtprcp, rtprcp, zm )
      call stat_update_var( ircp2, rcp2, zm )
      call stat_update_var( iupwp, upwp, zm )
      call stat_update_var( ivpwp, vpwp, zm )
      call stat_update_var( ivp2, vp2, zm )
      call stat_update_var( iup2, up2, zm )
      call stat_update_var( irho_zm, rho_zm, zm )
      call stat_update_var( isigma_sqd_w, sigma_sqd_w, zm )
      call stat_update_var( irho_ds_zm, rho_ds_zm, zm )
      call stat_update_var( iem, em, zm )
      call stat_update_var( ishear, shear, zm )
      call stat_update_var( iFrad, Frad, zm )

      if ( sclr_dim > 0 ) then
        do i=1, sclr_dim
          call stat_update_var( isclrp2(i), sclrp2(:,i), zm )
          call stat_update_var( isclrprtp(i), sclrprtp(:,i), zm )
          call stat_update_var( isclrpthvp(i), sclrpthvp(:,i), zm )
          call stat_update_var( isclrpthlp(i), sclrpthlp(:,i), zm )
          call stat_update_var( isclrprcp(i), sclrprcp(:,i), zm )
          call stat_update_var( iwpsclrp(i), wpsclrp(:,i), zm )
          call stat_update_var( iwp2sclrp(i), wp2sclrp(:,i), zm )
          call stat_update_var( iwpsclrp2(i), wpsclrp2(:,i), zm )
          call stat_update_var( iwpsclrprtp(i), wpsclrprtp(:,i), zm )
          call stat_update_var( iwpsclrpthlp(i), wpsclrpthlp(:,i), zm )
        end do
      endif
      if ( edsclr_dim > 0 ) then
        do i=1, edsclr_dim
          call stat_update_var( iwpedsclrp(i), wpedsclrp(:,i), zm )
        end do
      endif


      ! sfc variables

      ! Cloud cover
      call stat_update_var_pt( icc, 1, maxval( cloud_frac(1:gr%nnzp) ), sfc )

      ! Cloud base
      if ( izb > 0 ) then

        k = 1
        do while ( rcm(k) < rc_tol .and. k < gr%nnzp )
          k = k + 1
        enddo

        if ( k > 1 .AND. k < gr%nnzp) then

          ! Use linear interpolation to find the exact height of the
          ! rc_tol kg/kg level.  Brian.
          call stat_update_var_pt( izb, 1, lin_int( rc_tol, rcm(k),  &
                                   rcm(k-1), gr%zt(k), gr%zt(k-1) ), sfc )

        else

          ! Mark the cloud base at -10 m. if it's clear.
          call stat_update_var_pt( izb, 1, -10.0 , sfc )

        endif

      endif

      ! Liquid Water Path
      if ( ilwp > 0 ) then

        xtmp = 0.
        do i = gr%nnzp-1, 1, -1
          xtmp = xtmp + rho(i+1) * rcm(i+1) / gr%dzt(i+1)
        enddo

        call stat_update_var_pt( ilwp, 1, xtmp, sfc )

      endif

      ! Vapor Water Path (Preciptable Water)
      if ( ivwp > 0 ) then

        xtmp = 0.
        do i = gr%nnzp-1, 1, -1
          xtmp = xtmp + rho(i+1) * (rtm(i+1) - rcm(i+1)) / gr%dzt(i+1)
        enddo

        call stat_update_var_pt( ivwp, 1, xtmp, sfc )

      endif

      ! Snow Water Path
      if ( iswp > 0 ) then

        ! Calculate rsnowm
        if ( iirsnowm > 0 ) then
          rsnowm = hydromet(1:gr%nnzp,iirsnowm)
        else
          rsnowm = 0.0
        end if

        xtmp = 0.
        do i = gr%nnzp-1, 1, -1
          xtmp = xtmp + rho(i+1) * rsnowm(i+1) / gr%dzt(i+1)
        enddo

        call stat_update_var_pt( iswp, 1, xtmp, sfc )

      endif

      ! Ice Water Path
      if ( iiwp > 0 ) then

        ! Calculate ricem
        if ( iiricem > 0 ) then
          ricem = hydromet(1:gr%nnzp,iiricem)
        else
          ricem = 0.0
        end if

        xtmp = 0.
        do i = gr%nnzp-1, 1, -1
          xtmp = xtmp + rho(i+1) * ricem(i+1) / gr%dzt(i+1)
        enddo

        call stat_update_var_pt( iiwp, 1, xtmp, sfc )

      endif

      ! Vertical average of thermodynamic level variables.

      ! Find the vertical average of thermodynamic level variables, averaged from
      ! level 2 (the first thermodynamic level above model surface) through
      ! level gr%nnzp (the top of the model).  Use the vertical averaging function
      ! found in fill_holes.F90.

      ! Vertical average of thlm.
      call stat_update_var_pt( ithlm_vert_avg, 1,  &
           vertical_avg( 2, gr%nnzp, "zt", thlm(2:gr%nnzp) ), sfc )

      ! Vertical average of rtm.
      call stat_update_var_pt( irtm_vert_avg, 1,  &
           vertical_avg( 2, gr%nnzp, "zt", rtm(2:gr%nnzp) ), sfc )

      ! Vertical average of um.
      call stat_update_var_pt( ium_vert_avg, 1,  &
           vertical_avg( 2, gr%nnzp, "zt", um(2:gr%nnzp) ), sfc )

      ! Vertical average of vm.
      call stat_update_var_pt( ivm_vert_avg, 1,  &
           vertical_avg( 2, gr%nnzp, "zt", vm(2:gr%nnzp) ), sfc )

      ! Vertical average of momentum level variables.

      ! Find the vertical average of momentum level variables, averaged over the
      ! entire vertical profile (level 1 through level gr%nnzp).  Use the vertical
      ! averaging function found in fill_holes.F90.

      ! Vertical average of wp2.
      call stat_update_var_pt( iwp2_vert_avg, 1,  &
           vertical_avg( 1, gr%nnzp, "zm", wp2(1:gr%nnzp) ), sfc )

      ! Vertical average of up2.
      call stat_update_var_pt( iup2_vert_avg, 1,  &
           vertical_avg( 1, gr%nnzp, "zm", up2(1:gr%nnzp) ), sfc )

      ! Vertical average of vp2.
      call stat_update_var_pt( ivp2_vert_avg, 1,  &
           vertical_avg( 1, gr%nnzp, "zm", vp2(1:gr%nnzp) ), sfc )

      ! Vertical average of rtp2.
      call stat_update_var_pt( irtp2_vert_avg, 1,  &
           vertical_avg( 1, gr%nnzp, "zm", rtp2(1:gr%nnzp) ), sfc )

      ! Vertical average of thlp2.
      call stat_update_var_pt( ithlp2_vert_avg, 1,  &
           vertical_avg( 1, gr%nnzp, "zm", thlp2(1:gr%nnzp) ), sfc )


    endif  ! l_stats_samp


    return
  end subroutine stats_accumulate

  !-----------------------------------------------------------------------
  subroutine stats_finalize( )

    !     Description:
    !     Close NetCDF files and deallocate scratch space and
    !     stats file structures.
    !-----------------------------------------------------------------------

    use stats_variables, only: & 
        zt,  & ! Variable(s)
        zm, &
        rad_zt, &
        rad_zm, & 
        sfc, & 
        l_netcdf, & 
        l_stats, &
        l_output_rad_files

    use stats_variables, only: & 
        ztscr01, &  ! Variable(s)
        ztscr02, & 
        ztscr03, & 
        ztscr04, & 
        ztscr05, & 
        ztscr06, & 
        ztscr07, & 
        ztscr08, & 
        ztscr09, & 
        ztscr10, & 
        ztscr11, & 
        ztscr12, & 
        ztscr13, & 
        ztscr14, & 
        ztscr15, & 
        ztscr16, & 
        ztscr17, & 
        ztscr18, & 
        ztscr19, & 
        ztscr20, & 
        ztscr21

    use stats_variables, only: & 
        zmscr01, &  ! Variable(s)
        zmscr02, & 
        zmscr03, & 
        zmscr04, & 
        zmscr05, & 
        zmscr06, & 
        zmscr07, & 
        zmscr08, & 
        zmscr09, & 
        zmscr10, & 
        zmscr11, & 
        zmscr12, & 
        zmscr13, & 
        zmscr14, & 
        zmscr15, & 
        zmscr16, & 
        zmscr17

    !use stats_variables, only: & 
    !    radscr01, &  ! Variable(s)
    !    radscr02, & 
    !    radscr03, & 
    !    radscr04, & 
    !    radscr05, & 
    !    radscr06, & 
    !    radscr07, & 
    !    radscr08, & 
    !    radscr09, & 
    !    radscr10, & 
    !    radscr11, & 
    !    radscr12, & 
    !    radscr13, & 
    !    radscr14, & 
    !    radscr15, & 
    !    radscr16, & 
    !    radscr17        

    use stats_variables, only: & 
      isclrm, & 
      isclrm_f, & 
      iedsclrm, & 
      iedsclrm_f, & 
      isclrprtp, & 
      isclrp2, & 
      isclrpthvp, & 
      isclrpthlp, & 
      isclrprcp, & 
      iwpsclrp, & 
      iwp2sclrp, & 
      iwpsclrp2, & 
      iwpsclrprtp, & 
      iwpsclrpthlp, & 
      iwpedsclrp

#ifdef NETCDF
    use output_netcdf, only:  & 
        close_netcdf ! Procedure
#endif

    implicit none

    if ( l_stats .and. l_netcdf ) then
#ifdef NETCDF
      call close_netcdf( zt%f )
      call close_netcdf( zm%f )
      call close_netcdf( rad_zt%f )
      call close_netcdf( rad_zm%f )
      call close_netcdf( sfc%f )
#else
      stop "This program was not compiled with netCDF support"
#endif
    end if

    if ( l_stats ) then
      ! De-allocate all zt variables
      deallocate( zt%z )

      deallocate( zt%x )

      deallocate( zt%n )
      deallocate( zt%l_in_update )


      deallocate( zt%f%var )
      deallocate( zt%f%z )
      deallocate( zt%f%rlat )
      deallocate( zt%f%rlon )

      deallocate ( ztscr01 )
      deallocate ( ztscr02 )
      deallocate ( ztscr03 )
      deallocate ( ztscr04 )
      deallocate ( ztscr05 )
      deallocate ( ztscr06 )
      deallocate ( ztscr07 )
      deallocate ( ztscr08 )
      deallocate ( ztscr09 )
      deallocate ( ztscr10 )
      deallocate ( ztscr11 )
      deallocate ( ztscr12 )
      deallocate ( ztscr13 )
      deallocate ( ztscr14 )
      deallocate ( ztscr15 )
      deallocate ( ztscr16 )
      deallocate ( ztscr17 )
      deallocate ( ztscr18 )
      deallocate ( ztscr19 )
      deallocate ( ztscr20 )
      deallocate ( ztscr21 )

      ! De-allocate all zm variables
      deallocate( zm%z )

      deallocate( zm%x )
      deallocate( zm%n )

      deallocate( zm%f%var )
      deallocate( zm%f%z )
      deallocate( zm%l_in_update )

      deallocate ( zmscr01 )
      deallocate ( zmscr02 )
      deallocate ( zmscr03 )
      deallocate ( zmscr04 )
      deallocate ( zmscr05 )
      deallocate ( zmscr06 )
      deallocate ( zmscr07 )
      deallocate ( zmscr08 )
      deallocate ( zmscr09 )
      deallocate ( zmscr10 )
      deallocate ( zmscr11 )
      deallocate ( zmscr12 )
      deallocate ( zmscr13 )
      deallocate ( zmscr14 )
      deallocate ( zmscr15 )
      deallocate ( zmscr16 )
      deallocate ( zmscr17 )

      if (l_output_rad_files) then
        ! De-allocate all rad_zt variables
        deallocate( rad_zt%z )

        deallocate( rad_zt%x )
        deallocate( rad_zt%n )

        deallocate( rad_zt%f%var )
        deallocate( rad_zt%f%z )
        deallocate( rad_zt%l_in_update )

        ! De-allocate all rad_zm variables
        deallocate( rad_zm%z )

        deallocate( rad_zm%x )
        deallocate( rad_zm%n )

        deallocate( rad_zm%f%var )
        deallocate( rad_zm%f%z )
        deallocate( rad_zm%l_in_update )

        !deallocate ( radscr01 )
        !deallocate ( radscr02 )
        !deallocate ( radscr03 )
        !deallocate ( radscr04 )
        !deallocate ( radscr05 )
        !deallocate ( radscr06 )
        !deallocate ( radscr07 )
        !deallocate ( radscr08 )
        !deallocate ( radscr09 )
        !deallocate ( radscr10 )
        !deallocate ( radscr11 )
        !deallocate ( radscr12 )
        !deallocate ( radscr13 )
        !deallocate ( radscr14 )
        !deallocate ( radscr15 )
        !deallocate ( radscr16 )
        !deallocate ( radscr17 )
      end if ! l_output_rad_files

      ! De-allocate all sfc variables
      deallocate( sfc%z )

      deallocate( sfc%x )
      deallocate( sfc%n )
      deallocate( sfc%l_in_update )

      deallocate( sfc%f%var )
      deallocate( sfc%f%z )

      ! De-allocate scalar indices
      deallocate( isclrm )
      deallocate( isclrm_f )
      deallocate( iedsclrm )
      deallocate( iedsclrm_f )
      deallocate( isclrprtp )
      deallocate( isclrp2 )
      deallocate( isclrpthvp )
      deallocate( isclrpthlp )
      deallocate( isclrprcp )
      deallocate( iwpsclrp )
      deallocate( iwp2sclrp )
      deallocate( iwpsclrp2 )
      deallocate( iwpsclrprtp )
      deallocate( iwpsclrpthlp )
      deallocate( iwpedsclrp )

    endif ! l_stats


    return
  end subroutine stats_finalize

end module stats_subs
