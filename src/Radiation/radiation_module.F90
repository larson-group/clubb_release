module radiation_module

  use clubb_precision, only: time_precision ! Variable(s)

  use stats_netcdf, only: stats_type, stats_update ! Type / procedure

  implicit none

  public :: advance_clubb_radiation, &
            radiation_driver, &
            update_radiation_variables, &
            silhs_radiation_driver

  contains

!-------------------------------------------------------------------------------
  subroutine advance_clubb_radiation( &
               gr, ngrdcol, hydromet_dim, pdf_dim, lh_num_samples, &
               l_rad_itime, dt_main, day, month, year, &
               lat_vals, lon_vals, time_current, time_initial, &
               rho, rho_zm, p_in_Pa, exner, &
               wpthlp_sfc, wprtp_sfc, p_sfc, &
               cloud_frac, ice_supersat_frac, &
               thlm, rtm, rcm, &
               X_nl_all_levs, lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, &
               lh_sample_point_weights, hydromet, hm_metadata, &
               stats, err_info, &
               deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K, &
               radht )
! Description:
!   Advance the active radiation scheme when needed and update
!   radiation statistics every sampling timestep.
!-------------------------------------------------------------------------------

    use grid_class, only: &
      grid

    use clubb_precision, only: &
      time_precision, &
      core_rknd

    use constants_clubb, only: &
      fstderr

    use model_flags, only: &
      l_silhs_rad

    use soil_vegetation, only: &
      l_soil_veg, &
      advance_soil_veg

    use radiation_variables_module, only: &
      Frad, &
      Frad_SW_up, &
      Frad_LW_up, &
      Frad_SW_down, &
      Frad_LW_down

    use clubb_api_module, only: &
      clubb_at_least_debug_level_api, &
      clubb_fatal_error

    use corr_varnce_module, only: &
      hm_metadata_type

    use err_info_type_module, only: &
      err_info_type

    implicit none

    type(grid), intent(in) :: &
      gr

    integer, intent(in) :: &
      ngrdcol, &
      hydromet_dim, &
      pdf_dim, &
      lh_num_samples, &
      day, month, year

    logical, intent(in) :: &
      l_rad_itime

    real( kind = core_rknd ), intent(in) :: &
      dt_main, &
      lat_vals, &
      lon_vals

    real( kind = time_precision ), intent(in) :: &
      time_current, &
      time_initial

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
      rho, &
      p_in_Pa, &
      exner, &
      cloud_frac, &
      ice_supersat_frac, &
      thlm, &
      rtm, &
      rcm

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(in) :: &
      rho_zm

    real( kind = core_rknd ), dimension(ngrdcol), intent(in) :: &
      wpthlp_sfc, &
      wprtp_sfc, &
      p_sfc

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,gr%nzt,pdf_dim), intent(in) :: &
      X_nl_all_levs

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,gr%nzt), intent(in) :: &
      lh_rt_clipped, &
      lh_thl_clipped, &
      lh_rc_clipped, &
      lh_sample_point_weights

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,hydromet_dim), intent(in) :: &
      hydromet

    type(hm_metadata_type), intent(in) :: &
      hm_metadata

    type(stats_type), intent(inout) :: &
      stats

    type(err_info_type), intent(inout) :: &
      err_info

    real( kind = core_rknd ), dimension(ngrdcol), intent(inout) :: &
      deep_soil_T_in_K, &
      sfc_soil_T_in_K, &
      veg_T_in_K

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(out) :: &
      radht

    !------------------------ Begin Code -----------------------

    ! Soil vegetation is not neccesarily just Radiation, but we put it here to avoid 
    ! having to call it from the main driver, since it is only used in cases that
    ! also use radiation - currently only gabls3. We could move it to the main driver 
    ! if we wanted to use it in more cases, but for now this is sufficient.
    if ( l_soil_veg ) then

      ! This modifies the soil and vegetation temperatures at the surface, so 
      ! we need only _sfc arrays or surface slices, e.g (:,1)
      call advance_soil_veg( ngrdcol, dt_main, rho_zm(:,1), &
                             Frad_SW_up(:,1), Frad_SW_down(:,1), Frad_LW_down(:,1), &
                             wpthlp_sfc, wprtp_sfc, p_sfc, stats, &
                             deep_soil_T_in_K, sfc_soil_T_in_K, veg_T_in_K )

    end if
    
    ! Only advance radiation if l_rad_itime is true
    if ( l_rad_itime ) then

      ! Advance a radiation scheme
      ! With this call ordering, snow and ice water mixing ratio will be
      ! updated by the microphysics, but thlm and rtm will not.  This
      ! somewhat inconsistent, but we would need to move the call to
      ! radiation before the call the microphysics to change this.
      ! -dschanen 17 Aug 2009
      if ( l_silhs_rad ) then

        !$acc update host( rho, rho_zm, p_in_Pa, exner, wpthlp_sfc, wprtp_sfc, p_sfc, &
        !$acc              cloud_frac, ice_supersat_frac, thlm, rtm, rcm, &
        !$acc              X_nl_all_levs, lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, &
        !$acc              lh_sample_point_weights )

        call silhs_radiation_driver( &
              gr, ngrdcol, gr%nzm, gr%nzt, lh_num_samples, pdf_dim, hydromet_dim, hm_metadata, & ! In
              day, month, year,                                                       & ! In
              lat_vals, lon_vals, &
              time_current, time_initial, rho, rho_zm, p_in_Pa, exner,               & ! In
              cloud_frac, ice_supersat_frac, X_nl_all_levs,                           & ! In
              lh_rt_clipped, lh_thl_clipped, lh_rc_clipped,                           & ! In
              lh_sample_point_weights, hydromet,                                      & ! In
              stats,                                                                  & ! InOut
              err_info,                                                               & ! InOut
              radht, Frad, Frad_SW_up, Frad_LW_up,                                   & ! Out
              Frad_SW_down, Frad_LW_down )                                            ! Out

        !$acc update device( radht, Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down  )

      else

        call radiation_driver( &
              gr, time_current, time_initial, hydromet_dim,            & ! In
              ngrdcol,                                                  & ! In
              day, month, year,                                        & ! In
              lat_vals, lon_vals, &
              rho, rho_zm, p_in_Pa,                                    & ! In
              exner, cloud_frac, ice_supersat_frac,                    & ! In
              thlm, rtm, rcm, hydromet,                                & ! In
              hm_metadata, stats,                                      & ! InOut
              err_info,                                                & ! InOut
              radht, Frad, Frad_SW_up, Frad_LW_up,                     & ! Out
              Frad_SW_down, Frad_LW_down )                               ! Out

      end if ! l_silhs_rad

      if ( clubb_at_least_debug_level_api( 0 ) ) then
        if ( any(err_info%err_code == clubb_fatal_error) ) then
          write(fstderr, *) err_info%err_header_global
          write(fstderr, *) "Fatal error in radiation, " &
                            // "check your parameter values and timestep"
          return
        end if
      end if

    end if ! l_rad_itime

    ! We update stats here each sample timestep - even if radiation is not advanced
    if ( stats%l_sample ) then

      call update_radiation_variables( ngrdcol, gr%nzm, gr%nzt, radht, Frad, Frad_SW_up, Frad_LW_up, &
                                       Frad_SW_down, Frad_LW_down, &
                                       stats )
    end if

  end subroutine advance_clubb_radiation

!-------------------------------------------------------------------------------
  subroutine radiation_driver( &
               gr, time_current, time_initial, hydromet_dim, &
               ngrdcol, &
               day, month, year, &
               lat_vals, lon_vals, &
               rho, rho_zm, p_in_Pa, &
               exner, cloud_frac, ice_supersat_frac, &
               thlm, rtm, rcm, hydromet, &
               hm_metadata, stats,         &
               err_info, &
               radht, Frad, Frad_SW_up, Frad_LW_up, &
               Frad_SW_down, Frad_LW_down )
! Description:
!   Compute a radiation tendency.

! References:
!   None
!-------------------------------------------------------------------------------

    use constants_clubb, only: fstderr  !-------------------------------- Constant(s)

    use grid_class, only: grid

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
      simple_rad, simple_rad_bomex, simple_rad_lba

    use rad_lwsw_module, only: sunray_sw

    use parameters_radiation, only: &
      eff_drop_radius, & ! Variable(s)
      omega, &
      alvdr, &
      gc

    use grid_class, only: ddzm ! Procedure(s)

    use constants_clubb, only: Cp ! Variable(s)

    use grid_class, only: zt2zm_api !--------------------------------------- Procedure

    use interpolation, only: binary_search !--- Procdure(s)

    use clubb_api_module, only: &
        clubb_at_least_debug_level_api, & !-------------------------------- Procedure(s)
        lin_interpolate_on_grid_api, &
        clubb_fatal_error                 !-------------------------------- Constant

#ifdef radoffline
    use bugsrad_driver, only: compute_bugsrad_radiation !--------------- Procedure(s)
#endif

    use clubb_precision, only: &
      dp, & !----------------------------------------------------------- double precision
      time_precision, & !----------------------------------------------- Variable(s)
      core_rknd

    use corr_varnce_module, only: &
        hm_metadata_type

    use radiation_variables_module, only: &
      radht_LW, & !--------------------------------------- Variable(s)
      radht_SW, &
      Frad_SW, &
      Frad_LW

    use err_info_type_module, only: &
      err_info_type        ! Type

    implicit none

    !------------------------ Input Variables ------------------------
    type (grid), intent(in) :: &
      gr

    real(kind=time_precision), intent(in) :: &
      time_current, & ! Current time (UTC)               [s]
      time_initial    ! Start time of model run (UTC)    [s]

    integer, intent(in) ::  &
      day, month, year               ! Start time the of simulation

    real( kind = core_rknd ), intent(in) ::  & 
      lat_vals, & ! Latitude  [Degrees North]
      lon_vals    ! Longitude [Degrees East]

    integer, intent(in) :: &
      hydromet_dim, &
      ngrdcol

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(in) :: &
      rho,              & ! Density on thermo. grid                          [kg/m^3]
      p_in_Pa,          & ! Pressure.                                        [Pa] 
      exner,            & ! Exner function.                                  [-]
      cloud_frac,       & ! Cloud fraction (thermodynamic levels)            [-]
      ice_supersat_frac,& ! Ice cloud fraction (thermodynamic levels)        [-]
      thlm,             & ! Liquid potential temperature                     [K]
      rtm,              & ! Total water mixing ratio, r_t (thermo. levels)   [kg/kg]
      rcm                 ! Cloud water mixing ratio, r_c (thermo. levels)   [kg/kg]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(in) :: &
      rho_zm              ! Density on moment. grid                          [kg/m^3]

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species                                         [units vary]

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    !------------------------ Input/Output Variables ------------------------
    type(stats_type), intent(inout) :: &
      stats

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt), intent(out) :: &
      radht ! Radiative heating rate                                         [K/s]

    !------------------------ Output Variables ------------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzm), intent(out) :: &
      Frad,         & ! Total radiative flux                   [W/m^2]
      Frad_SW_up,   & ! Short-wave upwelling radiative flux    [W/m^2]
      Frad_LW_up,   & ! Long-wave upwelling radiative flux     [W/m^2]
      Frad_SW_down, & ! Short-wave upwelling radiative flux    [W/m^2]
      Frad_LW_down    ! Long-wave upwelling radiative flux     [W/m^2]

    !------------------------ Local Variables ------------------------
    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) :: &
      rsm,   & ! Snow mixing ratio                             [kg/kg]
      rim       ! Prisitine ice water mixing ratio             [kg/kg]

    real( kind = core_rknd ) :: Fs0, amu0_core_rknd

    real( kind = dp ) :: amu0 ! Cosine of the solar zenith angle [-]

    integer :: i, k, amu0_index

    real( kind = core_rknd ), dimension(ngrdcol,gr%nzt) :: &
      radht_SW_ddzm   ! Vertical derivative of the shortwave heating rate [K/s/m]

    ! Toggle for centered/forward differencing (in sunray_sw interpolations)
    ! To use centered differencing, set the toggle to .true.
    ! To use forward differencing, set the toggle to  .false.
    logical, parameter :: &
      l_center = .true.

    !------------------------ Begin Code ------------------------

    ! Initialize all outputs to 0.
    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzm
      do i = 1, ngrdcol
        Frad(i,k) = 0._core_rknd
        Frad_SW_up(i,k) = 0._core_rknd
        Frad_LW_up(i,k) = 0._core_rknd
        Frad_SW_down(i,k) = 0._core_rknd
        Frad_LW_down(i,k) = 0._core_rknd
      end do
    end do

    !$acc parallel loop gang vector collapse(2) default(present)
    do k = 1, gr%nzt
      do i = 1, ngrdcol
        radht(i,k) = 0._core_rknd ! Initialize the radiative heating rate to 0.
      end do
    end do

    ! If l_fix_cos_solar_zen is not set in the model.in, calculate amu0
    ! Otherwise, it was defined in cos_solar_zen_list file
    if ( l_sw_radiation ) then
      if ( l_fix_cos_solar_zen ) then
        if ( nparam > 1 ) then
          ! Find the closest time value greater than or equal to time_current
          amu0_index = binary_search( nparam, cos_solar_zen_times(1:nparam), &
                real( time_current, kind = core_rknd ) )
        else
          amu0_index = 1
        end if
        if ( amu0_index /= -1 ) then
          amu0 = dble( cos_solar_zen_values(amu0_index) )
        else
          write(fstderr,*) "Time not found in cos_solar_zen_times"
          error stop "Critical error."
        end if

      else ! Compute using the formula
        amu0 = cos_solar_zen( day, month, year, time_current, lat_vals, lon_vals )

      end if
    else
      amu0 = 0._dp ! This should disable shortwave radiation
    end if ! l_sw_radiation

    select case ( trim( rad_scheme ) )

    case ( "bugsrad" )
      !----------------------------------------------------------------
      ! BUGSrad Radiation
      !----------------------------------------------------------------
#ifdef radoffline /*This directive is needed for BUGSrad to work with CLUBB.*/

      !$acc update host( rho, rho_zm, p_in_Pa, exner, &
      !$acc              cloud_frac, ice_supersat_frac, thlm, rtm, rcm,&
      !$acc              Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down, &
      !$acc              radht_LW, radht_SW, Frad_SW, Frad_LW )

      ! Copy snow and ice
      rsm = 0.0_core_rknd
      rim = 0.0_core_rknd

      if ( hm_metadata%iirs > 0 ) then
        rsm = hydromet(:,:, hm_metadata%iirs)
      end if

      if ( hm_metadata%iiri > 0 ) then
        rim = hydromet(:,:, hm_metadata%iiri)
      end if

      call compute_bugsrad_radiation &
           ( gr, ngrdcol, gr%nzm, gr%nzt, amu0,                    & ! Intent(in)
             thlm, rcm, rtm, rsm, rim,                             & ! Intent(in)
             cloud_frac, ice_supersat_frac,                        & ! Intent(in)
             p_in_Pa, exner, rho_zm,                               & ! Intent(in)
             err_info,                                             & ! Intent(inout)
             radht, Frad,                                          & ! Intent(out)
             Frad_SW_up, Frad_LW_up,                               & ! Intent(out)
             Frad_SW_down, Frad_LW_down )                          ! Intent(out)

      !$acc update device( radht, Frad, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down, &
      !$acc                radht_LW, radht_SW, Frad_SW, Frad_LW )
#else

      error stop "Cannot call BUGSrad with these compile options."

#endif /*radoffline*/

    case ( "simplified" )
      
      !----------------------------------------------------------------
      ! Simplified radiation
      !----------------------------------------------------------------
      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzm
        do i = 1, ngrdcol
          Frad_SW(i,k) = 0._core_rknd
          Frad_LW(i,k) = 0._core_rknd
        end do
      end do

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzt
        do i = 1, ngrdcol
          radht_SW(i,k) = 0._core_rknd
          radht_LW(i,k) = 0._core_rknd
        end do
      end do

      if ( l_sw_radiation .and. amu0 > 0._dp ) then

        amu0_core_rknd = real( amu0, kind = core_rknd )
        
        if ( nparam > 1 ) then
          call lin_interpolate_on_grid_api( nparam, cos_solar_zen_values(1:nparam), &
                                            Fs_values(1:nparam), amu0_core_rknd, Fs0 )
        else
          Fs0 = Fs_values(1)
        end if

        call sunray_sw( ngrdcol, gr%nzt, rcm, rho, amu0_core_rknd, &
                        gr%dzt, gr%zm, gr%zt, &
                        eff_drop_radius, real( alvdr, kind = core_rknd ), &
                        gc, Fs0, omega, l_center, &
                        Frad_SW )

        !$acc data create( radht_SW_ddzm )
        radht_SW_ddzm = ddzm( gr%nzm, gr%nzt, ngrdcol, gr, Frad_SW )

        !$acc parallel loop gang vector collapse(2) default(present)
        do k = 1, gr%nzt
          do i = 1, ngrdcol
            radht_SW(i,k) = - radht_SW_ddzm(i,k) / (rho(i,k) * Cp)
          end do
        end do
        !$acc end data
      end if

      call simple_rad( gr, ngrdcol, rho, rho_zm, rtm, rcm, exner,  & ! In
                       stats, err_info,                   & ! Inout
                       Frad_LW, radht_LW )          ! Out

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzm
        do i = 1, ngrdcol
          Frad(i,k)  = Frad_SW(i,k)  + Frad_LW(i,k)
        end do
      end do

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzt
        do i = 1, ngrdcol
          radht(i,k) = radht_SW(i,k) + radht_LW(i,k)
        end do
      end do
      
    case ( "simplified_bomex" )
      !----------------------------------------------------------------
      ! GCSS BOMEX specifiction radiation
      !----------------------------------------------------------------

      call simple_rad_bomex( gr, ngrdcol, radht ) ! Out

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzm
        do i = 1, ngrdcol
          Frad(i,k) = 0._core_rknd
        end do
      end do 

    case ( "lba"  )

      call simple_rad_lba( gr, ngrdcol, time_current, time_initial, & ! In
                           radht )   ! Out

      !$acc update device( radht )

      !$acc parallel loop gang vector collapse(2) default(present)
      do k = 1, gr%nzm
        do i = 1, ngrdcol
          Frad(i,k) = 0._core_rknd
        end do
      end do 

    case default
      write(fstderr,*) "Undefined value for namelist variable rad_scheme: "//trim( rad_scheme )
      error stop "Fatal error encountered in radiation_driver."

    end select ! Radiation scheme

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( any(err_info%err_code == clubb_fatal_error) ) then
        write(fstderr, *) err_info%err_header_global
        write(fstderr,*) "Fatal error in radiation_driver:"
        return
      end if
    end if

    return
  end subroutine radiation_driver

  !-----------------------------------------------------------------------------
  subroutine update_radiation_variables( ngrdcol, nzm, nzt, radht, Frad, Frad_SW_up, Frad_LW_up, &
                                         Frad_SW_down, Frad_LW_down, &
                                         stats )

    ! Description:
    !   Updates the radiation variables using the stat_var_update() subroutine.
    !
    ! References:
    !   None
    !---------------------------------------------------------------------------

    use radiation_variables_module, only: &
      radht_LW, radht_SW, Frad_SW, Frad_LW, &
      extended_atmos_range_size, lin_int_buffer, &
      T_in_K, rcil, o3l, & !---------------------- Variables
      rsm_rad, rcm_in_cloud_rad, cloud_frac_rad, ice_supersat_frac_rad, radht_LW_rad, &
      radht_SW_rad, p_in_mb, sp_humidity, Frad_uLW, Frad_dLW, Frad_uSW, Frad_dSW, &
      fdswcl, fuswcl, fdlwcl, fulwcl

    use clubb_precision, only: &
      core_rknd

    use stats_netcdf, only: &
      stats_update

    use parameters_radiation, only: &
      rad_scheme ! Variable(s)

    implicit none

    ! Input Variables

    integer, intent(in) :: &
      ngrdcol, &
      nzm, & ! Model domain / # of momentum vertical levels          [-]
      nzt    ! Model domain / # of thermodynamic vertical levels     [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      radht           ! SW + LW heating rate               [K/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      Frad

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      Frad_SW_up,   & ! SW radiative upwelling flux        [W/m^2]
      Frad_LW_up,   & ! LW radiative upwelling flux        [W/m^2]
      Frad_SW_down, & ! SW radiative downwelling flux      [W/m^2]
      Frad_LW_down    ! LW radiative downwelling flux      [W/m^2]

    type(stats_type), intent(inout) :: &
      stats

    ! Local Variables

    integer :: rad_zt_dim, rad_zm_dim, i ! Dimensions of the radiation grid

    ! ---- Begin Code ----

    if (.not. stats%l_sample) return

    !$acc update host( Frad, radht, Frad_SW_up, Frad_LW_up, Frad_SW_down, Frad_LW_down )
    call stats_update( "Frad", Frad, stats )
    call stats_update( "radht", radht, stats )
    call stats_update( "Frad_SW_up", Frad_SW_up, stats )
    call stats_update( "Frad_LW_up", Frad_LW_up, stats )
    call stats_update( "Frad_SW_down", Frad_SW_down, stats )
    call stats_update( "Frad_LW_down", Frad_LW_down, stats )

    select case ( trim( rad_scheme ) )
    case ( "simplified", "bugsrad" )
      !$acc update host( radht_LW, radht_SW, Frad_SW, Frad_LW )
      ! Copy all the last column to all columns in stats - this is a legacy bug artifact
      ! from when clubb was single column, just being temporarily preserved.
      do i = 1, ngrdcol
        call stats_update( "radht_LW", radht_LW(ngrdcol,:), stats, i )
        call stats_update( "radht_SW", radht_SW(ngrdcol,:), stats, i )
        call stats_update( "Frad_SW", Frad_SW(ngrdcol,:), stats, i )
        call stats_update( "Frad_LW", Frad_LW(ngrdcol,:), stats, i )
      end do
    end select

    if ( stats%l_output_rad_files .and. trim( rad_scheme ) == "bugsrad" ) then
      rad_zt_dim = nzt + lin_int_buffer+extended_atmos_range_size
      rad_zm_dim = nzm + lin_int_buffer+extended_atmos_range_size

      ! Copy all the last column to all columns in stats - this is a legacy bug artifact
      ! from when clubb was single column, just being temporarily preserved.
      do i = 1, ngrdcol
        call stats_update( "T_in_K_rad", real( T_in_K(ngrdcol,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "rcil_rad", real( rcil(ngrdcol,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "o3l_rad", real( o3l(ngrdcol,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "rsm_rad", real( rsm_rad(ngrdcol,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "rcm_in_cloud_rad", real( rcm_in_cloud_rad(ngrdcol,rad_zt_dim:1:-1), &
                                                     kind=core_rknd ), stats, i )
        call stats_update( "cloud_frac_rad", real( cloud_frac_rad(ngrdcol,rad_zt_dim:1:-1), &
                                                     kind=core_rknd ), stats, i )
        call stats_update( "ice_supersat_frac_rad", real( ice_supersat_frac_rad(ngrdcol,rad_zt_dim:1:-1), &
                                                            kind=core_rknd ), stats, i )
        call stats_update( "radht_rad", real( radht_SW_rad(ngrdcol,rad_zt_dim:1:-1) + &
                                              radht_LW_rad(ngrdcol,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "radht_LW_rad", real( radht_LW_rad(ngrdcol,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "p_in_mb_rad", real( p_in_mb(ngrdcol,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "sp_humidity_rad", real( sp_humidity(ngrdcol,rad_zt_dim:1:-1), kind=core_rknd ), stats, i )

        call stats_update( "Frad_SW_rad", real( Frad_uSW(ngrdcol,rad_zm_dim:1:-1) - &
                                                Frad_dSW(ngrdcol,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "Frad_LW_rad", real( Frad_uLW(ngrdcol,rad_zm_dim:1:-1) - &
                                                Frad_dLW(ngrdcol,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "Frad_SW_up_rad", real( Frad_uSW(ngrdcol,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "Frad_LW_up_rad", real( Frad_uLW(ngrdcol,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "Frad_SW_down_rad", real( Frad_dSW(ngrdcol,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "Frad_LW_down_rad", real( Frad_dLW(ngrdcol,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )

        call stats_update( "fdswcl", real( fdswcl(ngrdcol,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "fuswcl", real( fuswcl(ngrdcol,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "fdlwcl", real( fdlwcl(ngrdcol,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
        call stats_update( "fulwcl", real( fulwcl(ngrdcol,rad_zm_dim:1:-1), kind=core_rknd ), stats, i )
      end do
    end if

  end subroutine update_radiation_variables

  !-----------------------------------------------------------------------
  subroutine silhs_radiation_driver( &
               gr, ngrdcol, nzm, nzt, lh_num_samples, pdf_dim, hydromet_dim, hm_metadata, &
               day, month, year, &
               lat_vals, lon_vals, &
               time_current, time_initial, rho, rho_zm, p_in_Pa, exner, &
               cloud_frac, ice_supersat_frac, X_nl_all_levs, &
               lh_rt_clipped, lh_thl_clipped, lh_rc_clipped, &
               lh_sample_point_weights, hydromet, stats,            &
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
      ngrdcol, &            ! Number of columns
      nzm, &                ! Number of momentum vertical levels
      nzt, &                ! Number of thermodynamic vertical levels
      lh_num_samples, &     ! Number of SILHS sample points
      pdf_dim, &            ! Number of lognormal variates
      hydromet_dim          ! Number of hydrometeor species

    integer, intent(in) :: &
      day, month, year

    real( kind = core_rknd ), intent(in) ::  & 
      lat_vals, lon_vals

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    real( kind = time_precision ), intent(in) :: &
      time_current, & ! Current time of simulation               [s]
      time_initial    ! Start time of simulation                 [s]

    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(in) :: &
      rho,               & ! Density on thermo. grid                   [kg/m^3]
      p_in_Pa,           & ! Pressure.                                 [Pa] 
      exner,             & ! Exner function.                           [-]
      cloud_frac,        & ! Cloud fraction (thermodynamic levels)     [-]
      ice_supersat_frac    ! Ice cloud fraction (thermodynamic levels) [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(in) :: &
      rho_zm               ! Density on moment. grid                   [kg/m^3]

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,nzt,pdf_dim), intent(in) :: &
      X_nl_all_levs        ! Normal-lognormal samples                  [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,nzt), intent(in) :: &
      lh_rt_clipped,  & ! rt generated from silhs sample points
      lh_thl_clipped, & ! thl generated from silhs sample points
      lh_rc_clipped     ! rc generated from silhs sample points

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,nzt), intent(in) :: &
      lh_sample_point_weights ! Weight of each SILHS sample point      [-]

    real( kind = core_rknd ), dimension(ngrdcol,nzt,hydromet_dim), intent(in) :: &
      hydromet             ! Hydrometeor mean fields

    type(stats_type), intent(inout) :: &
      stats

    type(err_info_type), intent(inout) :: &
      err_info        ! err_info struct containing err_code and err_header

    ! Output Variables
    real( kind = core_rknd ), dimension(ngrdcol,nzt), intent(out) :: &
      radht           ! Radiative heating rate                         [K/s]

    real( kind = core_rknd ), dimension(ngrdcol,nzm), intent(out) :: &
      Frad,         & ! Total radiative flux                           [W/m^2]
      Frad_SW_up,   & ! Short-wave upwelling radiative flux            [W/m^2]
      Frad_LW_up,   & ! Long-wave upwelling radiative flux             [W/m^2]
      Frad_SW_down, & ! Short-wave downwelling radiative flux          [W/m^2]
      Frad_LW_down    ! Long-wave downwelling radiative flux           [W/m^2]

    ! Local Variables
    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,nzt,hydromet_dim) :: &
      hydromet_all_pts ! SILHS sample of hydrometeors for each column  [units vary]

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,nzt) :: &
      Ncn_all_points   ! SILHS sample of Ncn for each column           [#/kg]
                       ! (not used)

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,nzt) :: &
      radht_samples              ! radht evaluated at each sample point

    real( kind = core_rknd ), dimension(ngrdcol,lh_num_samples,nzm) :: &
      Frad_samples,         &    ! Frad evaluated at each sample point
      Frad_SW_up_samples,   &    ! Frad_SW_up evaluated at each sample point
      Frad_LW_up_samples,   &    ! Frad_LW_up evaluated at each sample point
      Frad_SW_down_samples, &    ! Frad_SW_down evaluated at each sample point
      Frad_LW_down_samples       ! Frad_LW_down evaluated at each sample point

    type(stats_type) :: stats_dummy
    real( kind = core_rknd ) :: inv_nsample
    integer :: isample, i ! Looping variates

  !-----------------------------------------------------------------------

    !----- Begin Code -----

    inv_nsample = 1.0_core_rknd / real( lh_num_samples, kind = core_rknd )

    do i = 1, ngrdcol
      call copy_X_nl_into_hydromet_all_pts( &
             nzt, pdf_dim, lh_num_samples, &                ! Intent(in)
             X_nl_all_levs(i,:,:,:), &                      ! Intent(in)
             hydromet_dim, hm_metadata, &                   ! Intent(in)
             hydromet(i,:,:), &                            ! Intent(in)
             hydromet_all_pts(i,:,:,:), &                   ! Intent(out)
             Ncn_all_points(i,:,:) )                        ! Intent(out)
    end do

    ! We don't want radiation_driver updating stats, since
    ! we will be averaging over many samples before updating stats.  
    stats_dummy%enabled = .false.

    do isample=1, lh_num_samples
      ! Call a radiation scheme
      call radiation_driver( &
             gr, time_current, time_initial, hydromet_dim, &                        ! Intent(in)
             ngrdcol, &                                                            ! Intent(in)
             day, month, year, &
             lat_vals, lon_vals, &
             rho, rho_zm, p_in_Pa, & ! Intent(in)
             exner, cloud_frac, ice_supersat_frac, & ! Intent(in)
             lh_thl_clipped(:,isample,:), &     ! Intent(in)
             lh_rt_clipped(:,isample,:), lh_rc_clipped(:,isample,:), &                  ! Intent(in)
             hydromet_all_pts(:,isample,:,:), &                                       ! Intent(in)
             hm_metadata, stats_dummy, &                                        ! Intent(inout)
             err_info, &                                                            ! Intent(inout)
             radht_samples(:,isample,:), Frad_samples(:,isample,:), &                   ! Intent(out)
             Frad_SW_up_samples(:,isample,:), Frad_LW_up_samples(:,isample,:), &        ! Intent(out)
             Frad_SW_down_samples(:,isample,:), Frad_LW_down_samples(:,isample,:) )     ! Intent(out)
    end do

    ! Average results
    radht = inv_nsample * sum( radht_samples * lh_sample_point_weights, dim=2 )
    Frad = inv_nsample * sum( Frad_samples * lh_sample_point_weights, dim=2 )
    Frad_SW_up = inv_nsample * sum( Frad_SW_up_samples * lh_sample_point_weights, dim=2 )
    Frad_LW_up = inv_nsample * sum( Frad_LW_up_samples * lh_sample_point_weights, dim=2 )
    Frad_SW_down = inv_nsample * sum( Frad_SW_down_samples * lh_sample_point_weights, dim=2 )
    Frad_LW_down = inv_nsample * sum( Frad_LW_down_samples * lh_sample_point_weights, dim=2 )

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
