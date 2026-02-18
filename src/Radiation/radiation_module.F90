module radiation_module

  use clubb_precision, only: time_precision ! Variable(s)

  use stats_netcdf, only: stats_type ! Type

  implicit none

  public :: advance_clubb_radiation, &
            update_radiation_variables, &
            silhs_radiation_driver

  contains

!-------------------------------------------------------------------------------
  subroutine advance_clubb_radiation( &
               gr, time_current, time_initial, hydromet_dim, &
               day, month, year, &
               lin_int_buffer, &
               extended_atmos_bottom_level, &
               extended_atmos_top_level, &
               extended_atmos_range_size, &
               lat_vals, lon_vals, &
               rho, rho_zm, p_in_Pa, &
               exner, cloud_frac, ice_supersat_frac, &
               thlm, rtm, rcm, hydromet, &
               hm_metadata, stats,         &
               icol, &
               err_info, &
               radht, Frad, Frad_SW_up, Frad_LW_up, &
               Frad_SW_down, Frad_LW_down )
! Description:
!   Compute a radiation tendency.

! References:
!   None
!-------------------------------------------------------------------------------

    use constants_clubb, only: fstderr  !-------------------------------- Constant(s)

    use grid_class, only: grid ! Type

    use numerical_check, only: is_nan_2d, rad_check !-------------------- Procedure(s)

    use parameters_radiation, only: &
      rad_scheme, & !---------------------------------------------------- Variable(s)
      nparam, &
      l_fix_cos_solar_zen, &
      l_sw_radiation, &
      Fs_values, &
      cos_solar_zen_values, &
      cos_solar_zen_times

    use cos_solar_zen_module, only: cos_solar_zen !--------------------- Procedure(s)

    use simple_rad_module, only: &
      simple_rad, simple_rad_bomex, simple_rad_lba, sunray_sw_wrap

    use grid_class, only: zt2zm_api !--------------------------------------- Procedure

    use interpolation, only: binary_search !--- Procdure(s)

    use clubb_api_module, only: &
        clubb_at_least_debug_level_api, & !-------------------------------- Procedure(s)
        lin_interpolate_on_grid_api, &
        clubb_fatal_error                 !-------------------------------- Constant

#ifdef radoffline
    use bugsrad_driver, only: compute_bugsrad_radiation !--------------- Procedure(s)
#endif

    use variables_radiation_module, only: &
      radht_LW, radht_SW, Frad_SW, Frad_LW

    use clubb_precision, only: &
      dp, & !----------------------------------------------------------- double precision
      time_precision, & !----------------------------------------------- Variable(s)
      core_rknd

    use corr_varnce_module, only: &
        hm_metadata_type

    use err_info_type_module, only: &
      err_info_type        ! Type

    implicit none

    ! External
    intrinsic :: trim

    ! Input Variables
    type (grid), intent(in) :: &
      gr

    real(kind=time_precision), intent(in) :: &
      time_current, & ! Current time (UTC)               [s]
      time_initial    ! Start time of model run (UTC)    [s]

    integer, intent(in) ::  & 
      day, month, year,            & ! Start time the of simulation
      lin_int_buffer,              & ! The number of interpolated levels between the computational grid
                                     ! and the extended atmosphere
      extended_atmos_bottom_level, & ! Bottom level of the extended atmosphere
      extended_atmos_top_level,    & ! Top level of the extended atmosphere
      extended_atmos_range_size      ! The number of levels in the extended atmosphere

    real( kind = core_rknd ), intent(in) ::  & 
      lat_vals, & ! Latitude  [Degrees North]
      lon_vals    ! Longitude [Degrees East]

    integer, intent(in) :: &
      hydromet_dim

    real( kind = core_rknd ), dimension(gr%nzt), intent(in) :: &
      rho,              & ! Density on thermo. grid                          [kg/m^3]
      p_in_Pa,          & ! Pressure.                                        [Pa] 
      exner,            & ! Exner function.                                  [-]
      cloud_frac,       & ! Cloud fraction (thermodynamic levels)            [-]
      ice_supersat_frac,& ! Ice cloud fraction (thermodynamic levels)        [-]
      thlm,             & ! Liquid potential temperature                     [K]
      rtm,              & ! Total water mixing ratio, r_t (thermo. levels)   [kg/kg]
      rcm                 ! Cloud water mixing ratio, r_c (thermo. levels)   [kg/kg]

    real( kind = core_rknd ), dimension(gr%nzm), intent(in) :: &
      rho_zm              ! Density on moment. grid                          [kg/m^3]

    real( kind = core_rknd ), dimension(gr%nzt,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species                                         [units vary]

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    ! Input/Output Variables
    type(stats_type), intent(inout) :: &
      stats

    integer, intent(in) :: &
      icol

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    real( kind = core_rknd ), dimension(gr%nzt), intent(out) :: &
      radht ! Radiative heating rate                                         [K/s]

    ! Output Variables
    real( kind = core_rknd ), dimension(gr%nzm), intent(out) :: &
      Frad,         & ! Total radiative flux                   [W/m^2]
      Frad_SW_up,   & ! Short-wave upwelling radiative flux    [W/m^2]
      Frad_LW_up,   & ! Long-wave upwelling radiative flux     [W/m^2]
      Frad_SW_down, & ! Short-wave upwelling radiative flux    [W/m^2]
      Frad_LW_down    ! Long-wave upwelling radiative flux     [W/m^2]

    ! Local Variables
    real( kind = core_rknd ), dimension(gr%nzt) :: &
      rsm,   & ! Snow mixing ratio                             [kg/kg]
      rim       ! Prisitine ice water mixing ratio             [kg/kg]

    real( kind = core_rknd ) :: Fs0, amu0_core_rknd

    real( kind = dp ) :: amu0 ! Cosine of the solar zenith angle [-]

    integer :: i

    ! ---- Begin Code ----

    ! Initialize all outputs to 0.
    Frad = 0._core_rknd   ! The addition is to prevent an Intel compiler warning of an unused
    Frad_SW_up = 0._core_rknd       ! variable.  May be removed if rho is used below.  -meyern
    Frad_LW_up = 0._core_rknd
    Frad_SW_down = 0._core_rknd
    Frad_LW_down = 0._core_rknd

    radht = 0._core_rknd ! Initialize the radiative heating rate to 0.

    ! If l_fix_cos_solar_zen is not set in the model.in, calculate amu0
    ! Otherwise, it was defined in cos_solar_zen_list file
    if ( l_sw_radiation ) then
      if ( l_fix_cos_solar_zen ) then
        if ( nparam > 1 ) then
          ! Find the closest time value greater than or equal to time_current
          i = binary_search( nparam, cos_solar_zen_times(1:nparam), &
                real( time_current, kind = core_rknd ) )
        else
          i = 1
        end if
        if ( i /= -1 ) then
          amu0 = dble( cos_solar_zen_values(i) )
        else
          write(fstderr,*) "Time not found in cos_solar_zen_times"
          error stop "Critical error."
        end if

      else ! Compute using the formula
        amu0 = cos_solar_zen( day, month, year, time_current, lat_vals, lon_vals )

      end if ! l_fix_cos_solar_zen
    else
      amu0 = 0._dp ! This should disable shortwave radiation

    end if ! l_sw_radiation

    select case ( trim( rad_scheme ) )

    case ( "bugsrad" )
      !----------------------------------------------------------------
      ! BUGSrad Radiation
      !----------------------------------------------------------------
#ifdef radoffline /*This directive is needed for BUGSrad to work with CLUBB.*/

      ! Copy snow and ice
      if ( hm_metadata%iirs > 0 ) then
        rsm = hydromet(1:gr%nzt,hm_metadata%iirs)
      else
        rsm = 0.0_core_rknd
      endif

      if ( hm_metadata%iiri > 0 ) then
        rim = hydromet(1:gr%nzt,hm_metadata%iiri)
      else
        rim = 0.0_core_rknd
      end if

      ! NaN checks added to detect possible errors with BUGSrad
      ! Joshua Fasching November 2007

      if ( clubb_at_least_debug_level_api( 0 ) ) then

        if ( is_nan_2d( thlm ) ) then
          write(fstderr,*) "thlm before BUGSrad is NaN"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( rcm ) ) then
          write(fstderr,*) "rcm before BUGSrad is NaN"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( rtm ) ) then
          write(fstderr,*) "rtm before BUGSrad is NaN"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( rsm ) ) then
          write(fstderr,*) "rsm before BUGSrad is NaN"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( rim ) ) then
          write(fstderr,*) "rim before BUGSrad is NaN"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( cloud_frac ) ) then
          write(fstderr,*) "cloud_frac before BUGSrad is NaN"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( p_in_Pa ) ) then
          write(fstderr,*) "p_in_Pa before BUGSrad is NaN"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( exner ) ) then
          write(fstderr,*) "exner before BUGSrad is NaN"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( rho_zm ) ) then
          write(fstderr,*) "rho_zm before BUGSrad is NaN"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        ! Check for impossible negative values
        call rad_check( gr%nzm, gr%nzt, thlm, rcm, rtm, rim, & ! Intent(in)
                        cloud_frac, p_in_Pa, exner, rho_zm,  & ! Intent(in)
                        err_info )                             ! Intent(inout)

      end if  ! clubb_at_least_debug_level_api( 0 )

      call compute_bugsrad_radiation &
           ( gr%zm(1,:), gr%nzm, gr%nzt, lin_int_buffer,       & ! Intent(in)
             extended_atmos_range_size,                        & ! Intent(in)
             extended_atmos_bottom_level,                      & ! Intent(in)
             extended_atmos_top_level,                         & ! Intent(in)
             amu0,                                             & ! Intent(in)
             thlm, rcm, rtm, rsm, rim,                         & ! Intent(in)
             cloud_frac, ice_supersat_frac,                    & ! Intent(in)
             p_in_Pa, zt2zm_api( gr, p_in_Pa ), exner, rho_zm, & ! Intent(in)
             radht, Frad,                                      & ! Intent(out)
             Frad_SW_up, Frad_LW_up,                           & ! Intent(out)
             Frad_SW_down, Frad_LW_down )                        ! Intent(out)

      if ( clubb_at_least_debug_level_api( 0 ) ) then

        if ( is_nan_2d( Frad ) ) then
          write(fstderr,*) "Frad after BUGSrad is NaN"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

        if ( is_nan_2d( radht ) ) then
          write(fstderr,*) "radht after BUGSrad is NaN"
          ! General error -> set all entries to clubb_fatal_error
          err_info%err_code = clubb_fatal_error
        end if

      end if  ! clubb_at_least_debug_level_api( 2 )

#else

      error stop "Cannot call BUGSrad with these compile options."

#endif /*radoffline*/

    case ( "simplified" )
      !----------------------------------------------------------------
      ! Simplified radiation
      !----------------------------------------------------------------

      ! The sunray_sw code cannot handle negative values of cosine
      ! so we check that the value of amu0 is positive here.
      if ( l_sw_radiation .and. amu0 > 0._dp ) then
        amu0_core_rknd = real( amu0, kind = core_rknd )
        if ( nparam > 1 ) then
          call lin_interpolate_on_grid_api( nparam, cos_solar_zen_values(1:nparam), &
                                            Fs_values(1:nparam), amu0_core_rknd, Fs0 )
        else
          Fs0 = Fs_values(1)
        end if
        call sunray_sw_wrap( gr, Fs0, amu0_core_rknd, rho, rcm, & ! In
                             Frad_SW, radht_SW )              ! Out
      else
        radht_SW = 0._core_rknd
        Frad_SW  = 0._core_rknd
      end if

      call simple_rad( gr, rho, rho_zm, rtm, rcm, exner, & ! In
                       stats, icol, err_info,           & ! Inout
                       Frad_LW, radht_LW )                ! Out


      Frad = Frad_SW + Frad_LW
      radht = radht_SW + radht_LW

    case ( "simplified_bomex" )
      !----------------------------------------------------------------
      ! GCSS BOMEX specifiction radiation
      !----------------------------------------------------------------

      call simple_rad_bomex( gr, radht ) ! Out

    case ( "lba"  )
      call simple_rad_lba( gr, time_current, time_initial, & ! In
                           radht )   ! Out

    case ( "none" )
      radht_SW = 0._core_rknd
      Frad_SW  = 0._core_rknd
      radht    = 0._core_rknd

    case default
      write(fstderr,*) "Undefined value for namelist variable rad_scheme: "//trim( rad_scheme )
      error stop "Fatal error encountered in advance_clubb_radiation."

    end select ! Radiation scheme

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "Fatal error in advance_clubb_radiation:"
        return
      end if
    end if
    if ( stats%l_sample ) then
      continue
    end if

    return
  end subroutine advance_clubb_radiation

  !-----------------------------------------------------------------------------
  subroutine update_radiation_variables( ngrdcol, nzm, nzt, radht, Frad_SW_up, Frad_LW_up, &
                                         Frad_SW_down, Frad_LW_down, &
                                         extended_atmos_range_size, lin_int_buffer, &
                                         stats )

    ! Description:
    !   Updates the radiation variables using the stat_var_update() subroutine.
    !
    ! References:
    !   None
    !---------------------------------------------------------------------------

    use variables_radiation_module, only: &
      radht_LW, radht_SW, Frad_SW, Frad_LW, T_in_k, rcil, o3l, & !---------------------- Variables
      rsm_2d, rcm_in_cloud_2d, cloud_frac_2d, ice_supersat_frac_2d, radht_LW_2d, &
      radht_SW_2d, p_in_mb, sp_humidity, Frad_uLW, Frad_dLW, Frad_uSW, Frad_dSW, &
      fdswcl, fuswcl, fdlwcl, fulwcl

    use grid_class, only: &
      flip !------------------------------------------------------------------------- Prodecure(s)

    use clubb_precision, only: &
      core_rknd

    use stats_netcdf, only: &
      stats_type, &
      stats_update

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      ngrdcol, &
      nzm, & ! Model domain / # of momentum vertical levels          [-]
      nzt    ! Model domain / # of thermodynamic vertical levels     [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      radht           ! SW + LW heating rate               [K/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      Frad_SW_up,   & ! SW radiative upwelling flux        [W/m^2]
      Frad_LW_up,   & ! LW radiative upwelling flux        [W/m^2]
      Frad_SW_down, & ! SW radiative downwelling flux      [W/m^2]
      Frad_LW_down    ! LW radiative downwelling flux      [W/m^2]

    integer, intent(in) :: &
      extended_atmos_range_size, &
      lin_int_buffer

    type(stats_type), intent(inout) :: &
      stats

    ! Local Variables

    integer :: rad_zt_dim, rad_zm_dim, i, rad_col ! Dimensions of the radiation grid

    ! ---- Begin Code ----

    if (.not. stats%l_sample) return

    call stats_update( "radht", radht, stats )
    call stats_update( "Frad_SW_up", Frad_SW_up, stats )
    call stats_update( "Frad_LW_up", Frad_LW_up, stats )
    call stats_update( "Frad_SW_down", Frad_SW_down, stats )
    call stats_update( "Frad_LW_down", Frad_LW_down, stats )

    ! This is duplicating these fields to all columns
    do i = 1, ngrdcol
      call stats_update( "radht_LW", radht_LW, stats, i )
      call stats_update( "radht_SW", radht_SW, stats, i )
      call stats_update( "Frad_SW", Frad_SW, stats, i )
      call stats_update( "Frad_LW", Frad_LW, stats, i )
    end do

    if ( stats%l_output_rad_files ) then
      rad_zt_dim = nzt + lin_int_buffer+extended_atmos_range_size
      rad_zm_dim = nzm + lin_int_buffer+extended_atmos_range_size
      do i = 1, ngrdcol
        rad_col = min( i, size( T_in_K, 1 ) )
        call stats_update( "T_in_K_rad", real( T_in_K(rad_col,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "rcil_rad", real( rcil(rad_col,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "o3l_rad", real( o3l(rad_col,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "rsm_rad", real( rsm_2d(rad_col,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "rcm_in_cloud_rad", real( rcm_in_cloud_2d(rad_col,rad_zt_dim:1:-1), &
                                                     kind=core_rknd ), stats, i )
        call stats_update( "cloud_frac_rad", real( cloud_frac_2d(rad_col,rad_zt_dim:1:-1), &
                                                   kind=core_rknd ), stats, i )
        call stats_update( "ice_supersat_frac_rad", real( ice_supersat_frac_2d(rad_col,rad_zt_dim:1:-1), &
                                                          kind=core_rknd ), stats, i )
        call stats_update( "radht_rad", real( radht_SW_2d(rad_col,rad_zt_dim:1:-1) + &
                                              radht_LW_2d(rad_col,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "radht_LW_rad", real( radht_LW_2d(rad_col,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "p_in_mb_rad", real( p_in_mb(rad_col,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "sp_humidity_rad", real( sp_humidity(rad_col,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )

        rad_col = min( i, size( Frad_uSW, 1 ) )
        call stats_update( "Frad_SW_rad", real( Frad_uSW(rad_col,rad_zm_dim:1:-1) - &
                                                Frad_dSW(rad_col,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "Frad_LW_rad", real( Frad_uLW(rad_col,rad_zm_dim:1:-1) - &
                                                Frad_dLW(rad_col,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "Frad_SW_up_rad", real( Frad_uSW(rad_col,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "Frad_LW_up_rad", real( Frad_uLW(rad_col,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "Frad_SW_down_rad", real( Frad_dSW(rad_col,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "Frad_LW_down_rad", real( Frad_dLW(rad_col,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )

        rad_col = min( i, size( fdswcl, 1 ) )
        call stats_update( "fdswcl", real( fdswcl(rad_col,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "fuswcl", real( fuswcl(rad_col,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "fdlwcl", real( fdlwcl(rad_col,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "fulwcl", real( fulwcl(rad_col,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
      end do
    end if

  end subroutine update_radiation_variables

  !-----------------------------------------------------------------------
  subroutine silhs_radiation_driver( &
               gr, nzm, nzt, lh_num_samples, pdf_dim, hydromet_dim, hm_metadata, &
               day, month, year, &
               lin_int_buffer, &
               extended_atmos_bottom_level, &
               extended_atmos_top_level, &
               extended_atmos_range_size, &
               lat_vals, lon_vals, &
               time_current, time_initial, rho, rho_zm, p_in_Pa, exner, &
               cloud_frac, ice_supersat_frac, X_nl_all_levs, &
               lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, &
               lh_sample_point_weights, hydromet, stats, icol,         &
               err_info, &
               radht, Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down )

  ! Description:
  !   Computes radiation over a set of sample points and averages the
  !   results

  ! References
  !   clubb:ticket:663
  !-----------------------------------------------------------------------

    ! Included Modules
    use grid_class, only: &
      grid
    
    use clubb_precision, only: &
      core_rknd

    use clubb_api_module, only: &
      clubb_at_least_debug_level_api, & ! Procedure
      clubb_fatal_error             ! Constant

    use latin_hypercube_driver_module, only: &
      copy_X_nl_into_hydromet_all_pts   !--------------------- Procedure

    use constants_clubb, only: &
      fstderr        !------------------------------------------- Constant

    use corr_varnce_module, only: &
      hm_metadata_type

    use err_info_type_module, only: &
      err_info_type        ! Type

    implicit none

    ! Input Variables
    type(grid), intent(in) :: &
      gr

    integer, intent(in) :: &
      nzm, &                ! Number of momentum vertical levels
      nzt, &                ! Number of thermodynamic vertical levels
      lh_num_samples, &     ! Number of SILHS sample points
      pdf_dim, &            ! Number of lognormal variates
      hydromet_dim          ! Number of hydrometeor species

    integer, intent(in) :: &
      day, month, year, &
      lin_int_buffer, &
      extended_atmos_bottom_level, &
      extended_atmos_top_level, &
      extended_atmos_range_size

    real( kind = core_rknd ), intent(in) ::  & 
      lat_vals, lon_vals

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = time_precision ), intent(in) :: &
      time_current, & ! Current time of simulation               [s]
      time_initial    ! Start time of simulation                 [s]

    real( kind = core_rknd ), dimension(nzt), intent(in) :: &
      rho,               & ! Density on thermo. grid                   [kg/m^3]
      p_in_Pa,           & ! Pressure.                                 [Pa] 
      exner,             & ! Exner function.                           [-]
      cloud_frac,        & ! Cloud fraction (thermodynamic levels)     [-]
      ice_supersat_frac    ! Ice cloud fraction (thermodynamic levels) [-]

    real( kind = core_rknd ), dimension(nzm), intent(in) :: &
      rho_zm               ! Density on moment. grid                   [kg/m^3]

    real( kind = core_rknd ), dimension(lh_num_samples,nzt,pdf_dim), intent(in) :: &
      X_nl_all_levs        ! Normal-lognormal samples                  [units vary]

    real( kind = core_rknd ), dimension(lh_num_samples,nzt), intent(in) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped     ! rc generated from silhs sample points

    real( kind = core_rknd ), dimension(lh_num_samples,nzt), intent(in) :: &
      lh_sample_point_weights ! Weight of each SILHS sample point      [-]

    real( kind = core_rknd ), dimension(nzt,hydromet_dim), intent(in) :: &
      hydromet             ! Hydrometeor mean fields

    type(stats_type), intent(inout) :: &
      stats

    integer, intent(in) :: &
      icol

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    ! Output Variables
    real( kind = core_rknd ), dimension(nzt), intent(out) :: &
      radht           ! Radiative heating rate                         [K/s]

    real( kind = core_rknd ), dimension(nzm), intent(out) :: &
      Frad,         & ! Total radiative flux                           [W/m^2]
      Frad_SW_up,   & ! Short-wave upwelling radiative flux            [W/m^2]
      Frad_LW_up,   & ! Long-wave upwelling radiative flux             [W/m^2]
      Frad_SW_down, & ! Short-wave downwelling radiative flux          [W/m^2]
      Frad_LW_down    ! Long-wave downwelling radiative flux           [W/m^2]

    ! Local Variables
    real( kind = core_rknd ), dimension(lh_num_samples,nzt,hydromet_dim) :: &
      hydromet_all_pts ! SILHS sample of hydrometeors for each column  [units vary]

    real( kind = core_rknd ), dimension(lh_num_samples,nzt) :: &
      Ncn_all_points   ! SILHS sample of Ncn for each column           [#/kg]
                       ! (not used)

    real( kind = core_rknd ), dimension(lh_num_samples,nzt) :: &
      radht_samples              ! radht evaluated at each sample point

    real( kind = core_rknd ), dimension(lh_num_samples,nzm) :: &
      Frad_samples,         &    ! Frad evaluated at each sample point
      Frad_SW_up_samples,   &    ! Frad_SW_up evaluated at each sample point
      Frad_LW_up_samples,   &    ! Frad_LW_up evaluated at each sample point
      Frad_SW_down_samples, &    ! Frad_SW_down evaluated at each sample point
      Frad_LW_down_samples       ! Frad_LW_down evaluated at each sample point

    type(stats_type) :: stats_dummy
    integer :: isample, k ! Looping variates

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    call copy_X_nl_into_hydromet_all_pts( &
           nzt, pdf_dim, lh_num_samples, &                ! Intent(in)
           X_nl_all_levs, &                               ! Intent(in)
           hydromet_dim, hm_metadata, &                   ! Intent(in)
           hydromet, &                                    ! Intent(in)
           hydromet_all_pts, &                            ! Intent(out)
           Ncn_all_points )                               ! Intent(out)

    stats_dummy%enabled = .false.

    do isample=1, lh_num_samples
      ! Call a radiation scheme
      call advance_clubb_radiation( &
             gr, time_current, time_initial, hydromet_dim, &                        ! Intent(in)
             day, month, year, &
             lin_int_buffer, &
             extended_atmos_bottom_level, &
             extended_atmos_top_level, &
             extended_atmos_range_size, &
             lat_vals, lon_vals, &
             rho, rho_zm, p_in_Pa, &                                                ! Intent(in)
             exner, cloud_frac, ice_supersat_frac, lh_thl_clipped(isample,:), &     ! Intent(in)
             lh_rt_clipped(isample,:), lh_rc_clipped(isample,:), &                  ! Intent(in)
             hydromet_all_pts(isample,:,:), &                                       ! Intent(in)
             hm_metadata, stats_dummy, &                                        ! Intent(inout)
             icol, &                                                                ! Intent(in)
             err_info, &                                                            ! Intent(inout)
             radht_samples(isample,:), Frad_samples(isample,:), &                   ! Intent(out)
             Frad_SW_up_samples(isample,:), Frad_LW_up_samples(isample,:), &        ! Intent(out)
             Frad_SW_down_samples(isample,:), Frad_LW_down_samples(isample,:) )     ! Intent(out)
    end do

    ! Average results
    forall ( k = 1:nzt )

      radht(k) = sum( radht_samples(:,k) * lh_sample_point_weights(:,k) ) / &
                  real( lh_num_samples, kind=core_rknd )
    end forall

    forall ( k = 1:nzm )

      Frad(k)  = sum( Frad_samples(:,k) * lh_sample_point_weights(:,k) ) / &
                  real( lh_num_samples, kind=core_rknd )

      Frad_SW_up(k) = sum( Frad_SW_up_samples(:,k) * lh_sample_point_weights(:,k) ) / &
                       real( lh_num_samples, kind=core_rknd )

      Frad_LW_up(k) = sum( Frad_LW_up_samples(:,k) * lh_sample_point_weights(:,k) ) / &
                       real( lh_num_samples, kind=core_rknd )

      Frad_SW_down(k)  = sum( Frad_SW_down_samples(:,k) * lh_sample_point_weights(:,k) ) / &
                          real( lh_num_samples, kind=core_rknd )
                          
      Frad_LW_down(k)  = sum( Frad_LW_down_samples(:,k) * lh_sample_point_weights(:,k) ) / &
                          real( lh_num_samples, kind=core_rknd )

    end forall

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "Fatal error in silhs_radiation_driver:"
        return
      end if
    end if

    return
  end subroutine silhs_radiation_driver

end module radiation_module
