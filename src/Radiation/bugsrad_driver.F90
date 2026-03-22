!-------------------------------------------------------------------------------
! $Id$
module bugsrad_driver

  implicit none

  public :: compute_bugsrad_radiation, init_radiation

  private ! Default Scope

  contains

  subroutine compute_bugsrad_radiation &
             ( gr, ngrdcol, nzm, nzt, amu0, &
               thlm, rcm, rtm, rsm, rim, &
               cloud_frac, ice_supersat_frac, &
               p_in_Pa, exner, rho_zm, &
               err_info, &
               radht, Frad, &
               Frad_SW_up, Frad_LW_up, &
               Frad_SW_down, Frad_LW_down )

! Description:
!   Does the necessary operations to interface the CLUBB model with
!   the BUGS radition scheme.
!
! References:
!   Stevens, et al., (2001) _Journal of Atmospheric Science_, Vol 58, p.3391-3409
!
! Contact for information on BUGSrad (other than this routine)
!   Norm Wood <norm@atmos.colostate.edu>
!
! All code external to this based on the BUGSrad source from 2004/7/10, with
! modifications to use LAPACK rather than Numerical Recipes for band-diagonal
! solvers when -Dnooverlap is not added to the CPPFLAGS.
!-------------------------------------------------------------------------------

    use constants_clubb, only: fstderr, grav, Cp, cloud_frac_min, &
                               pascal_per_mb, g_per_kg ! Variable(s)

    use clubb_precision, only: dp, core_rknd ! Variable(s)

    use grid_class, only: &
      grid, &
      zt2zm_api

    use T_in_K_module, only: thlm2T_in_K_api ! Procedure(s)

    use error_code, only: &
      clubb_at_least_debug_level_api, &
      clubb_fatal_error

    use numerical_check, only: &
      is_nan_2d, &
      rad_check

    use extended_atmosphere_module, only: &
      extended_atmos_dim, extended_alt, extended_p_in_mb, & ! Variable(s)
      extended_T_in_K, extended_sp_hmdty, extended_o3l

    use parameters_radiation, only: &
      sol_const, & ! Variable(s)
      alvdr, &
      alvdf, &
      alndr, &
      alndf, &
      slr

    use radiation_variables_module, only: &
      lin_int_buffer, extended_atmos_range_size, &
      extended_atmos_bottom_level, extended_atmos_top_level, &
      radht_LW, radht_SW, Frad_SW, Frad_LW, &
      T_in_K, rcil, o3l, rsm_rad, rcm_in_cloud_rad, cloud_frac_rad, ice_supersat_frac_rad, &
      radht_SW_rad, radht_LW_rad, Frad_uLW, Frad_dLW, Frad_uSW, Frad_dSW, &
      p_in_mb, sp_humidity, &
      fdswcl, fuswcl, fdlwcl, fulwcl

    use err_info_type_module, only: &
      err_info_type

    implicit none

    ! External
    external :: bugs_rad ! Subroutine

    intrinsic :: real ! Function

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      amu0  ! Cosine of the solar zenith angle  [-]

    integer, intent(in) :: &
      ngrdcol, &
      nzm, & ! Vertical extent;  i.e. nz in the grid class
      nzt

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzt) :: &
      thlm,             & ! Liquid potential temp.              [K]
      rcm,              & ! Liquid water mixing ratio           [kg/kg]
      rsm,              & ! Snow water mixing ratio             [kg/kg]
      rim,              & ! Ice water mixing ratio              [kg/kg]
      rtm,              & ! Total water mixing ratio            [kg/kg]
      cloud_frac,       & ! Cloud fraction                      [-]
      ice_supersat_frac,& ! Ice cloud fraction                  [-]
      p_in_Pa,          & ! Pressure on the t grid              [Pa]
      exner               ! Exner function                      [-]

    type(grid), intent(in) :: &
      gr

    real( kind = core_rknd ), intent(in), dimension(ngrdcol,nzm) :: &
      rho_zm              ! Density                             [kg/m^3]

    type(err_info_type), intent(inout) :: &
      err_info

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nzm) :: &
      Frad,         & ! Total radiative flux          [W/m^2]
      Frad_SW_up,   & ! SW radiative upwelling flux   [W/m^2]
      Frad_LW_up,   & ! LW radiative upwelling flux   [W/m^2]
      Frad_SW_down, & ! SW radiative downwelling flux [W/m^2]
      Frad_LW_down    ! LW radiative downwelling flux [W/m^2]

    real( kind = core_rknd ), intent(out), dimension(ngrdcol,nzt) :: &
      radht           ! Total heating rate            [K/s]

    ! Local Variables

    real( kind = core_rknd ), dimension(ngrdcol,nzm) :: &
      p_in_Pam ! Pressure on the m grid [Pa]

    real( kind = core_rknd ), dimension(ngrdcol,nzt) :: &
      rcm_in_cloud  ! Liquid water mixing ratio in cloud  [kg/kg]

    ! Pressure in millibars for layers (calculated as an average of p_in_mb)
    real( kind = dp ), dimension(ngrdcol,(nzm-1)+lin_int_buffer+extended_atmos_range_size+1) :: &
      playerinmb ! [hPa]

    real( kind = dp ), dimension(ngrdcol,(nzm-1)+lin_int_buffer+extended_atmos_range_size) :: &
      diff_pres_lvls  ! Difference in pressure levels       [hPa]

    real( kind = dp ), dimension(ngrdcol) :: &
      ts  ! Surface temperature [K]

    real( kind = dp ) :: z1_fact, z2_fact ! Temp storage

    integer :: i, z, z1, z2  ! Loop indices

    integer :: buffer ! The sum of the two buffers

    !character(len=40) :: time_char

    !-------------------------------------------------------------------------------

    buffer = lin_int_buffer + extended_atmos_range_size
    if ( gr%zm(1,nzm) > extended_alt(extended_atmos_dim) ) then

      write(fstderr,*) "The CLUBB model grid (for zm levels) contains an ",  &
                       "altitude above the top of the extended atmosphere ",  &
                       "profile."
      write(fstderr,*) "Top of CLUBB model zm grid =", gr%zm(1,nzm), "m."
      write(fstderr,*) "Top of extended atmosphere profile =",  &
                       extended_alt(extended_atmos_dim), "m."
      write(fstderr,*) "Reduce the vertical extent of the CLUBB model grid."
      ! CLUBB zm grid exceeds a 50 km altitude
      error stop "compute_bugsrad_radiation: cannot handle this altitude"

    else
      ! Continue
    end if

    T_in_K = 0._dp
    rcil = 0._dp
    o3l = 0._dp
    rsm_rad = 0._dp
    rcm_in_cloud_rad = 0._dp
    cloud_frac_rad = 0._dp
    ice_supersat_frac_rad = 0._dp
    p_in_mb = 0._dp
    sp_humidity = 0._dp
    radht_SW = 0._core_rknd
    radht_LW = 0._core_rknd
    radht_SW_rad = 0._dp
    radht_LW_rad = 0._dp
    Frad_dSW = 0._dp
    Frad_uSW = 0._dp
    Frad_dLW = 0._dp
    Frad_uLW = 0._dp
    fdswcl = 0._dp
    fuswcl = 0._dp
    fdlwcl = 0._dp
    fulwcl = 0._dp
    playerinmb = 0._dp
    diff_pres_lvls = 0._dp
    ts = 0._dp

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( is_nan_2d( thlm ) ) then
        write(fstderr,*) "thlm before BUGSrad is NaN"
        err_info%err_code = clubb_fatal_error
      end if

      if ( is_nan_2d( rcm ) ) then
        write(fstderr,*) "rcm before BUGSrad is NaN"
        err_info%err_code = clubb_fatal_error
      end if

      if ( is_nan_2d( rtm ) ) then
        write(fstderr,*) "rtm before BUGSrad is NaN"
        err_info%err_code = clubb_fatal_error
      end if

      if ( is_nan_2d( rsm ) ) then
        write(fstderr,*) "rsm before BUGSrad is NaN"
        err_info%err_code = clubb_fatal_error
      end if

      if ( is_nan_2d( rim ) ) then
        write(fstderr,*) "rim before BUGSrad is NaN"
        err_info%err_code = clubb_fatal_error
      end if

      if ( is_nan_2d( cloud_frac ) ) then
        write(fstderr,*) "cloud_frac before BUGSrad is NaN"
        err_info%err_code = clubb_fatal_error
      end if

      if ( is_nan_2d( p_in_Pa ) ) then
        write(fstderr,*) "p_in_Pa before BUGSrad is NaN"
        err_info%err_code = clubb_fatal_error
      end if

      if ( is_nan_2d( exner ) ) then
        write(fstderr,*) "exner before BUGSrad is NaN"
        err_info%err_code = clubb_fatal_error
      end if

      if ( is_nan_2d( rho_zm ) ) then
        write(fstderr,*) "rho_zm before BUGSrad is NaN"
        err_info%err_code = clubb_fatal_error
      end if

      call rad_check( ngrdcol, gr%nzm, gr%nzt, thlm, rcm, rtm, rim, &
                      cloud_frac, p_in_Pa, exner, rho_zm, err_info )
    end if

    p_in_Pam = zt2zm_api( gr%nzm, gr%nzt, ngrdcol, gr, p_in_Pa )

    rcm_in_cloud = rcm / max( cloud_frac, cloud_frac_min )
    T_in_K(:,1:nzt) = thlm2T_in_K_api( nzt, ngrdcol, thlm, exner, rcm )
    p_in_mb(:,1:nzt) = p_in_Pa(:,1:nzt) / pascal_per_mb
    o3l(:,1:nzt) = ( 5.4e-5_core_rknd / rho_zm(:,1:nzt) ) / g_per_kg

    where ( rtm(:,1:nzt) < rcm(:,1:nzt) )
      sp_humidity(:,1:nzt) = 0.0_dp
    elsewhere
      sp_humidity(:,1:nzt) = real( rtm(:,1:nzt) - rcm(:,1:nzt), kind=dp ) &
                             / real( 1.0_core_rknd + rtm(:,1:nzt), kind=dp )
    end where

    if ( clubb_at_least_debug_level_api( 1 ) ) then
      do i = 1, ngrdcol
        do z = 1, nzt
          if ( rtm(i,z) < rcm(i,z) ) then
            write(fstderr,*) "rvm < 0 at ", z, " before BUGSrad, specific humidity set to 0"
          end if
        end do
      end do
    end if

    playerinmb(:,1:nzm) = p_in_Pam(:,1:nzm) / pascal_per_mb

    rcil(:,buffer+1:nzt+buffer) = rim(:,nzt:1:-1)
    rsm_rad(:,buffer+1:nzt+buffer) = rsm(:,nzt:1:-1)
    rcm_in_cloud_rad(:,buffer+1:nzt+buffer) = rcm_in_cloud(:,nzt:1:-1)
    cloud_frac_rad(:,buffer+1:nzt+buffer) = cloud_frac(:,nzt:1:-1)
    ice_supersat_frac_rad(:,buffer+1:nzt+buffer) = ice_supersat_frac(:,nzt:1:-1)

    T_in_K(:,buffer+1:nzt+buffer) = T_in_K(:,nzt:1:-1)
    sp_humidity(:,buffer+1:nzt+buffer) = sp_humidity(:,nzt:1:-1)
    p_in_mb(:,buffer+1:nzt+buffer) = p_in_mb(:,nzt:1:-1)
    playerinmb(:,buffer+1:nzm+buffer) = playerinmb(:,nzm:1:-1)
    o3l(:,buffer+1:nzt+buffer) = o3l(:,nzt:1:-1)

    if ( extended_atmos_range_size > 0 ) then
      do i = 1, ngrdcol
        T_in_K(i,1:extended_atmos_range_size) = &
          extended_T_in_K( extended_atmos_top_level:extended_atmos_bottom_level:-1 )

        sp_humidity(i,1:extended_atmos_range_size) = &
          extended_sp_hmdty( extended_atmos_top_level:extended_atmos_bottom_level:-1 )

        o3l(i,1:extended_atmos_range_size) = &
          extended_o3l( extended_atmos_top_level:extended_atmos_bottom_level:-1 )

        p_in_mb(i,1:extended_atmos_range_size) = &
          extended_p_in_mb( extended_atmos_top_level:extended_atmos_bottom_level:-1 )
      end do
    end if

    z1 = buffer + 1
    z2 = extended_atmos_range_size
    do z = buffer, extended_atmos_range_size+1, -1
      z1_fact = real( z2 - z, kind=dp ) / real( z2 - z1,kind=dp )
      z2_fact = real( z - z1,kind=dp ) / real( z2 - z1, kind=dp )

      T_in_K(:,z) = z1_fact * T_in_K(:,z1) + z2_fact * T_in_K(:,z2)
      sp_humidity(:,z) = z1_fact * sp_humidity(:,z1) + z2_fact * sp_humidity(:,z2)
      o3l(:,z) = z1_fact * o3l(:,z1) + z2_fact * o3l(:,z2)
      p_in_mb(:,z) = z1_fact * p_in_mb(:,z1) + z2_fact * p_in_mb(:,z2)
    end do

    playerinmb(:,2:buffer+1) = ( p_in_mb(:,1:buffer) + p_in_mb(:,2:buffer+1) ) / 2._dp

    where ( 2._dp * playerinmb(:,2) - playerinmb(:,3) > 0._dp )
      playerinmb(:,1) = 2._dp * playerinmb(:,2) - playerinmb(:,3)
    elsewhere
      playerinmb(:,1) = 0.5_dp * playerinmb(:,2)
    end where

    diff_pres_lvls(:,1:(nzm-1)+buffer) = playerinmb(:,2:(nzm-1)+buffer+1) - &
                                         playerinmb(:,1:(nzm-1)+buffer)

    ts(:) = T_in_K(:,nzt+buffer)

    do i = 1, ngrdcol
      do z = 2, nzt + buffer
        if ( p_in_mb(i,z-1) > p_in_mb(i,z) ) then
          write(0,*) "Column = ", i
          write(0,*) "Pressure (z)   [mb] = ", p_in_mb(i,z)
          write(0,*) "Pressure (z-1) [mb] = ", p_in_mb(i,z-1)
          write(0,*) "z level = ", z
          error stop "Fatal error: assertion check for pressure in BUGSrad_driver failed"
        end if
      end do
    end do

    radht_SW_rad(:,:) = 0._dp
    radht_LW_rad(:,:) = 0._dp
    Frad_dSW(:,:) = 0._dp
    Frad_uSW(:,:) = 0._dp
    Frad_dLW(:,:) = 0._dp
    Frad_uLW(:,:) = 0._dp
    fdswcl(:,:) = 0._dp
    fuswcl(:,:) = 0._dp
    fdlwcl(:,:) = 0._dp
    fulwcl(:,:) = 0._dp

    do i = 1, ngrdcol
      call bugs_rad( 1, 1, nzt+buffer, playerinmb(i:i,:), p_in_mb(i:i,:), &
                     diff_pres_lvls(i:i,:), T_in_K(i:i,:), sp_humidity(i:i,:), &
                     rcm_in_cloud_rad(i:i,:), rcil(i:i,:), rsm_rad(i:i,:), o3l(i:i,:), &
                     ts(i:i), amu0, slr, alvdf, alndf, alvdr, alndr, sol_const, &
                     real( grav, kind=dp ), real( Cp, kind=dp ), radht_SW_rad(i:i,:), &
                     radht_LW_rad(i:i,:), Frad_dSW(i:i,:), Frad_uSW(i:i,:), &
                     Frad_dLW(i:i,:), Frad_uLW(i:i,:), fdswcl(i:i,:), fuswcl(i:i,:), &
                     fdlwcl(i:i,:), fulwcl(i:i,:), cloud_frac_rad(i:i,:) )
    end do

    radht_SW(:,1:nzt) = real( radht_SW_rad(:,nzt+buffer:buffer+1:-1), kind = core_rknd ) &
                        * ( 1.0_core_rknd / exner(:,1:nzt) )
    radht_LW(:,1:nzt) = real( radht_LW_rad(:,nzt+buffer:buffer+1:-1), kind = core_rknd ) &
                        * ( 1.0_core_rknd / exner(:,1:nzt) )
    radht(:,1:nzt) = radht_SW(:,1:nzt) + radht_LW(:,1:nzt)

    Frad_SW_up(:,:) = Frad_uSW(:,nzm+buffer:buffer+1:-1)
    Frad_LW_up(:,:) = Frad_uLW(:,nzm+buffer:buffer+1:-1)
    Frad_SW_down(:,:) = Frad_dSW(:,nzm+buffer:buffer+1:-1)
    Frad_LW_down(:,:) = Frad_dLW(:,nzm+buffer:buffer+1:-1)
    Frad_SW(:,:) = Frad_SW_up(:,:) - Frad_SW_down(:,:)
    Frad_LW(:,:) = Frad_LW_up(:,:) - Frad_LW_down(:,:)
    Frad(:,:) = Frad_SW(:,:) + Frad_LW(:,:)

    if ( clubb_at_least_debug_level_api( 0 ) ) then
      if ( is_nan_2d( Frad ) ) then
        write(fstderr,*) "Frad after BUGSrad is NaN"
        err_info%err_code = clubb_fatal_error
      end if

      if ( is_nan_2d( radht ) ) then
        write(fstderr,*) "radht after BUGSrad is NaN"
        err_info%err_code = clubb_fatal_error
      end if
    end if

    return
  end subroutine compute_bugsrad_radiation
!-------------------------------------------------------------------------------

!-------------------------------------------------------------------------------
  subroutine init_radiation( iunit, namelist_file, case_info_file, &
                             l_calc_thlp2_rad )
! Description:
!   Setup radiation parameters

! References:
!   None
!-------------------------------------------------------------------------------
    use parameters_radiation, only: &
      rad_scheme, sol_const, alvdr, alvdf, alndr, alndf, &
      kappa, F0, F1, eff_drop_radius, gc, omega, Fs_values, &
      cos_solar_zen_values, cos_solar_zen_times, &
      radiation_top, l_fix_cos_solar_zen, l_sw_radiation, &
      slr, l_rad_above_cloud, l_use_default_std_atmosphere, nparam

    use error_code, only: clubb_at_least_debug_level_api ! Function

    use text_writer, only: &
      write_text   ! Used to write radiation settings to setup.txt file

    use constants_clubb, only: eps

    use clubb_precision, only: &
      dp, & ! double precision
      core_rknd

    implicit none

    ! Constant parameters
    logical, parameter :: &
      l_write_to_file = .true. ! If true, will write case information to a file.

    ! Input Variables
    integer, intent(in) :: iunit ! File unit

    character(len=*), intent(in) :: &
      namelist_file, & ! Filename containing the namelist
      case_info_file   ! Name of simulation info file (plain text)

    ! Input/Output Variables
    logical, intent(inout) :: &
      l_calc_thlp2_rad ! Include the contribution of radiation to thlp2

    ! Local variables
    integer :: k

    namelist /radiation_setting/ &
      rad_scheme, sol_const, alvdr, alvdf, alndr, alndf, &
      kappa, F0, F1, eff_drop_radius, gc, omega, Fs_values, &
      cos_solar_zen_values, cos_solar_zen_times, &
      radiation_top, l_fix_cos_solar_zen, l_sw_radiation, &
      slr, l_rad_above_cloud, l_use_default_std_atmosphere, &
      l_calc_thlp2_rad

    ! ---- Begin Code ----

    ! Set default values, then read in the namelist
    rad_scheme = "none"

    ! BUGSrad parameters
    sol_const =  1367._dp ! W/m^2

    alvdf = 0.1_dp ! Visible diffuse surface albedo       [-]
    alndr = 0.1_dp ! Near-IR direct surface albedo        [-]
    alndf = 0.1_dp ! Near-IR diffuse surface albedo       [-]

    ! 50000m is the top of the U.S. Standard Atmosphere data used
    ! in CLUBB.
    radiation_top = 50000._core_rknd! [m]

    ! Variables used by both schemes
    alvdr = 0.1_dp ! Visible direct surface albedo        [-]

    ! Simplified radiation parameters
    F0    = 100.0_core_rknd  ! Coefficient for cloud top heating (see Stevens) [W/m^2]
    F1    = 20.0_core_rknd   ! Coefficient for cloud base heating (see Stevens)[W/m^2]
    kappa = 119.0_core_rknd  ! A constant (Duynkerke eqn. 5)                   [m^2/kg]
    gc    = 0.86_core_rknd   ! Asymmetry parameter, "g" in Duynkerke           [-]
    omega = 0.9965_core_rknd ! Single-scattering albedo                        [-]

    slr  = 1.0_dp  ! Fraction of daylight

    l_rad_above_cloud = .false. ! For the heaviside step function
    l_sw_radiation = .false. ! Set to true to enable shortwave radiation

    ! Use the 1976 standard atmsophere table to add a buffer above the model
    ! domain for radiation.  Otherwise, the sounding will be used.
    l_use_default_std_atmosphere = .true. 

    ! Parameters for fixing the value of cosine of the solar zenith angle
    l_fix_cos_solar_zen = .false.

    ! The incident of incoming SW insolation at cloud top the
    ! direction of the incoming beam (not the vertical)   [W/m^2]
    Fs_values(:) = 0.0_core_rknd

    cos_solar_zen_values(:) = -999.0_core_rknd ! Cosine of the solar zenith angle [-]
    cos_solar_zen_times(:)  = -999.0_core_rknd ! Simulation times corresponding to above [s]

    eff_drop_radius = 1.e-5_core_rknd ! Effective droplet radius [m]

    ! Read the namelist values in
    open(unit=iunit, file=namelist_file, status='old',action='read')
    read(iunit, nml=radiation_setting)
    close(unit=iunit)

    if ( clubb_at_least_debug_level_api( 1 ) ) then

      ! This will open the cases setup.txt file and append it to include the
      ! parameters in the radiation_setting namelist. This file was created
      ! and written to from clubb_driver previously.
      if ( l_write_to_file ) open(unit=iunit, file=case_info_file, &
          status='old', action='write', position='append')

      ! Write to file all parameters from the namelist microphysics_seting.
      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )
      call write_text( "&radiation_setting", l_write_to_file, iunit)
      call write_text( "--------------------------------------------------", &
        l_write_to_file, iunit )

      call write_text ( "rad_scheme = " // rad_scheme, l_write_to_file, &
        iunit )
      call write_text ( "sol_const = ", real( sol_const, kind = core_rknd ), &
        l_write_to_file, iunit )
      call write_text ( "alvdr = ", real( alvdr, kind = core_rknd ), l_write_to_file, iunit )
      call write_text ( "alvdf = ", real( alvdf, kind = core_rknd ), l_write_to_file, iunit )
      call write_text ( "alndr = ", real( alndr, kind = core_rknd ), l_write_to_file, iunit )
      call write_text ( "alndf = ", real( alndf, kind = core_rknd ), l_write_to_file, iunit )
      call write_text ( "radiation_top = ", radiation_top, l_write_to_file, iunit )
      call write_text ( "F0 = ", F0, l_write_to_file, iunit )
      call write_text ( "F1 = ", F1, l_write_to_file, iunit )
      call write_text ( "kappa = ", kappa, l_write_to_file, iunit )
      call write_text ( "gc = ", gc, l_write_to_file, iunit )
      call write_text ( "omega = ", omega, l_write_to_file, iunit )
      call write_text ( "slr = ", real( slr, kind = core_rknd ), l_write_to_file, iunit )
      call write_text ( "l_rad_above_cloud = ", l_rad_above_cloud, l_write_to_file, iunit )
      call write_text ( "l_sw_radiation = ", l_sw_radiation, l_write_to_file, iunit )
      call write_text ( "l_fix_cos_solar_zen = ", l_fix_cos_solar_zen, l_write_to_file, iunit )
      call write_text ( "l_use_default_std_atmosphere = ", l_use_default_std_atmosphere, &
                        l_write_to_file, iunit )
      call write_text ( "Fs_values = ", Fs_values, l_write_to_file, iunit )
      call write_text ( "cos_solar_zen_values = ", cos_solar_zen_values, l_write_to_file, iunit )
      call write_text ( "cos_solar_zen_times = ", cos_solar_zen_times, l_write_to_file, iunit )
      call write_text ( "eff_drop_radius = ", eff_drop_radius, l_write_to_file, iunit )

      if ( l_write_to_file ) close(unit=iunit);

    end if ! clubb_at_least_debug_level_api( 1 )

    do k = 1, size( cos_solar_zen_values )
      if ( abs(cos_solar_zen_values(k) - (-999._core_rknd)) &
        < abs(cos_solar_zen_values(k) + (-999._core_rknd)) / 2 * eps ) then
        exit
      else
        nparam = k
      end if
    end do

    return
  end subroutine init_radiation

end module bugsrad_driver
