!-------------------------------------------------------------------------------
! $Id$
module bugsrad_driver

  implicit none

  public :: compute_bugsrad_radiation, init_radiation

  private ! Default Scope

  contains

  subroutine compute_bugsrad_radiation &
             ( alt, nzm, nzt, lin_int_buffer,         &
               extended_atmos_range_size,             &
               extended_atmos_bottom_level,           &
               extended_atmos_top_level,              &
               amu0, &
               thlm, rcm, rtm, rsm, rim,& 
               cloud_frac, ice_supersat_frac,       &
               p_in_Pa, p_in_Pam, exner, rho_zm,    &
               radht, Frad,                         &
               Frad_SW_up, Frad_LW_up,              &
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

    use T_in_K_module, only: thlm2T_in_K_api ! Procedure(s)

    use error_code, only: clubb_at_least_debug_level_api ! Procedure(s)

    use grid_class, only: flip  ! Procedure(s)

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

    use variables_radiation_module, only: &
      slen, nlen ! Constats

    use variables_radiation_module, only: &
      radht_LW, radht_SW, Frad_SW, Frad_LW, & ! Variable(s)
      T_in_K, rcil, o3l, rsm_2d, rcm_in_cloud_2d, cloud_frac_2d, ice_supersat_frac_2d, &
      radht_SW_2d, radht_LW_2d, Frad_uLW, Frad_dLW, Frad_uSW, Frad_dSW, &
      p_in_mb, sp_humidity

    implicit none

    ! External
    external :: bugs_rad ! Subroutine

    intrinsic :: real ! Function

    ! Input Variables
    real( kind = dp ), intent(in) :: &
      amu0  ! Cosine of the solar zenith angle  [-]

    integer, intent(in) :: &
      nzm, & ! Vertical extent;  i.e. nz in the grid class
      nzt

    ! Number of levels to interpolate from the bottom of extended_atmos to the top
    ! of the CLUBB profile, hopefully enough to eliminate cooling spikes, etc.
    integer, intent(in) :: lin_int_buffer, extended_atmos_range_size

    real( kind = core_rknd ), intent(in), dimension(nzt) :: &
      thlm,             & ! Liquid potential temp.              [K]
      rcm,              & ! Liquid water mixing ratio           [kg/kg]
      rsm,              & ! Snow water mixing ratio             [kg/kg]
      rim,              & ! Ice water mixing ratio              [kg/kg]
      rtm,              & ! Total water mixing ratio            [kg/kg]
      cloud_frac,       & ! Cloud fraction                      [-]
      ice_supersat_frac,& ! Ice cloud fraction                  [-]
      p_in_Pa,          & ! Pressure on the t grid              [Pa]
      exner               ! Exner function                      [-]

    real( kind = core_rknd ), intent(in), dimension(nzm) :: &
      alt,              & ! Altitudes of the model              [m]
      rho_zm,           & ! Density                             [kg/m^3]
      p_in_Pam            ! Pressure on the m grid              [Pa]

    integer,intent(in) ::&
      extended_atmos_bottom_level, &
      extended_atmos_top_level

    ! Output Variables
    real( kind = core_rknd ), intent(out), dimension(nzm) :: &
      Frad,         & ! Total radiative flux          [W/m^2]
      Frad_SW_up,   & ! SW radiative upwelling flux   [W/m^2]
      Frad_LW_up,   & ! LW radiative upwelling flux   [W/m^2]
      Frad_SW_down, & ! SW radiative downwelling flux [W/m^2]
      Frad_LW_down    ! LW radiative downwelling flux [W/m^2]

    real( kind = core_rknd ), intent(out), dimension(nzt) :: &
      radht           ! Total heating rate            [K/s]

    ! Local Variables

    real( kind = core_rknd ), dimension(nzt) :: &
      rcm_in_cloud  ! Liquid water mixing ratio in cloud  [kg/kg]

!   real( kind = dp ), dimension(nlen,nzt+lin_int_buffer+extended_atmos_range_size) :: &
!     sp_humidity, & ! Specific humidity      [kg/kg]
!     p_in_mb          ! Pressure in millibars  [hPa]

    ! Pressure in millibars for layers (calculated as an average of p_in_mb)
    real( kind = dp ), dimension(nlen,(nzm-1)+lin_int_buffer+extended_atmos_range_size+1) :: &
      playerinmb ! [hPa]

    real( kind = dp ), dimension(nlen,(nzm-1)+lin_int_buffer+extended_atmos_range_size) :: &
      diff_pres_lvls  ! Difference in pressure levels       [hPa]

    real( kind = dp ), dimension(nlen) :: &
      ts  ! Surface temperature [K]

    real( kind = dp ) :: z1_fact, z2_fact, tmp ! Temp storage

    integer :: i, z, z1, z2  ! Loop indices

    integer :: buffer ! The sum of the two buffers

    !character(len=40) :: time_char

    !-------------------------------------------------------------------------------

    buffer = lin_int_buffer + extended_atmos_range_size

    ! Convert to millibars
    p_in_mb(1,1:nzt)  = real( p_in_Pa(1:nzt) / pascal_per_mb, kind=dp ) ! t grid in CLUBB

    playerinmb(1,1:nzm) = real( p_in_Pam / pascal_per_mb, kind=dp ) ! m grid in CLUBB

    ! Determine rcm in cloud
    rcm_in_cloud = rcm / max( cloud_frac, cloud_frac_min )

    ! Convert theta_l to temperature

    T_in_K(1,1:nzt) = real( thlm2T_in_K_api( nzt, thlm(1:nzt), exner(1:nzt), rcm(1:nzt) ),kind=dp )

    ! Derive Specific humidity from rc & rt.
    do z = 1, nzt
      if ( rtm(z) < rcm(z) ) then
        sp_humidity(1,z) = 0.0_dp
        if ( clubb_at_least_debug_level_api( 1 ) ) then
          write(fstderr,*) "rvm < 0 at ", z, " before BUGSrad, specific humidity set to 0"
        endif
      else
        sp_humidity(1,z) &
          = real( rtm(z) - rcm(z),kind=dp ) / real( 1.0_core_rknd+rtm(z),kind=dp )
      end if
    end do

    ! Setup miscellaneous variables

    ! Ozone at < 1 km = 5.4e-5 g/m^3 from U.S. Standard Atmosphere, 1976.
    !   Convert from g to kg.
    o3l(1,1:nzt) = real( ( 5.4e-5_core_rknd / rho_zm(1:(nzm-1)) ) / &
       g_per_kg, kind=dp ) !Known magic number

    ! Convert and transpose as needed
    rcil(1,buffer+1:nzt+buffer)            = real( flip( rim(1:nzt), nzt ), kind = dp )
    rsm_2d(1,buffer+1:nzt+buffer)          = real( flip( rsm(1:nzt), nzt ), kind = dp )
    rcm_in_cloud_2d(1,buffer+1:nzt+buffer) = real( flip( rcm_in_cloud(1:nzt), nzt ), kind = dp )
    cloud_frac_2d(1,buffer+1:nzt+buffer)   = real( flip( cloud_frac(1:nzt), nzt ), kind = dp )
    ice_supersat_frac_2d(1,buffer+1:nzt+buffer) = real( flip( ice_supersat_frac(1:nzt), nzt ), &
                                                        kind = dp )

    T_in_K(1,buffer+1:nzt+buffer) = real( flip( real( T_in_K(1,1:nzt), kind = core_rknd ), &
                                                nzt ), kind = dp )

    sp_humidity(1,buffer+1:nzt+buffer) = real( flip( real( sp_humidity(1,1:nzt), &
                                                           kind = core_rknd ), &
                                                     nzt ), kind = dp )

    p_in_mb(1,(buffer+1):(nzt+buffer)) = real( flip( real( p_in_mb(1,1:nzt), kind = core_rknd ), &
                                                     nzt ), kind = dp )
                                          
    playerinmb(1,(buffer+1):(nzm+buffer)) = real( flip( real( playerinmb(1,1:nzm), &
                                                              kind = core_rknd ), &
                                                        nzm ), kind = dp )

    o3l(1,buffer+1:nzt+buffer) = real( flip( real( o3l(1,1:nzt), kind = core_rknd ), nzt ), &
                                             kind = dp )

    ! Assume these are all zero above the CLUBB profile
    rsm_2d(1,1:buffer)                = 0.0_dp
    rcil(1,1:buffer)                  = 0.0_dp
    rcm_in_cloud_2d(1,1:buffer)       = 0.0_dp
    cloud_frac_2d(1,1:buffer)         = 0.0_dp
    ice_supersat_frac_2d(1,1:buffer)  = 0.0_dp

    if ( alt(nzm) > extended_alt(extended_atmos_dim) ) then

      write(fstderr,*) "The CLUBB model grid (for zm levels) contains an ",  &
                       "altitude above the top of the extended atmosphere ",  &
                       "profile."
      write(fstderr,*) "Top of CLUBB model zm grid =", alt(nzm), "m."
      write(fstderr,*) "Top of extended atmosphere profile =",  &
                       extended_alt(extended_atmos_dim), "m."
      write(fstderr,*) "Reduce the vertical extent of the CLUBB model grid."
      ! CLUBB zm grid exceeds a 50 km altitude
      error stop "compute_bugsrad_radiation: cannot handle this altitude"

    else
      ! Continue
    end if

    ! Add the extended atmospheric profile above the linear interpolation
    T_in_K(1,1:extended_atmos_range_size) = &
               real( flip( extended_T_in_K( extended_atmos_bottom_level: &
                      extended_atmos_top_level ), extended_atmos_range_size ), kind = dp )

    sp_humidity(1,1:extended_atmos_range_size) = &
               real( flip( extended_sp_hmdty( extended_atmos_bottom_level: &
                      extended_atmos_top_level ), extended_atmos_range_size ), kind = dp )

    o3l(1,1:extended_atmos_range_size) = &
               real( flip( extended_o3l( extended_atmos_bottom_level: &
                      extended_atmos_top_level ), extended_atmos_range_size ), kind = dp )

    p_in_mb(1,1:extended_atmos_range_size) = &
               real( flip( extended_p_in_mb( extended_atmos_bottom_level: &
                      extended_atmos_top_level ), extended_atmos_range_size ), kind = dp )

    ! Do a linear interpolation to produce the levels between the extended
    ! atmospheric levels and the CLUBB levels;
    ! These levels should number the lin_int_buffer parameter
    z1 = buffer + 1
    z2 = extended_atmos_range_size
    do z = buffer, extended_atmos_range_size+1, -1
      z1_fact = real( z2 - z, kind=dp ) / real( z2 - z1,kind=dp )
      z2_fact = real( z - z1,kind=dp ) / real( z2 - z1, kind=dp )

      T_in_K(1,z) = z1_fact * T_in_K(1,z1) + z2_fact * T_in_K(1,z2)

      sp_humidity(1,z) = z1_fact * sp_humidity(1,z1) + z2_fact * sp_humidity(1,z2)
      o3l(1,z) = z1_fact * o3l(1,z1) + z2_fact * o3l(1,z2)

      p_in_mb(1,z) = z1_fact * p_in_mb(1,z1) + z2_fact * p_in_mb(1,z2)
    end do

    ! Do a linear interpolation to find playerinmb.  Since this interpolation
    ! occurs at levels above the top of the CLUBB model, the CLUBB zt2zm_api function
    ! or CLUBB weighted averages do not apply.  The variable playerinmb is being
    ! defined on momentum levels above the top of the CLUBB model, which are
    ! being defined here at points half-way inbetween the thermodynamic levels
    ! above the top of the CLUBB model.  Brian Griffin; May 13, 2008.
    playerinmb(1,2:buffer+1) = ( p_in_mb(1,1:buffer) + p_in_mb(1,2:buffer+1) ) / 2._dp

    ! Do a linear extension to find playerinmb at the uppermost standard
    ! atmosphere momentum level.  The grid is evenly-spaced at these points.
    ! Brian Griffin; May 13, 2008.
    tmp = 2._dp * playerinmb(1,2) - playerinmb(1,3)
    if ( tmp > 0._dp ) then
      playerinmb(1,1) = tmp
    else ! Assuming a linear extension didn't work
      playerinmb(1,1) = 0.5_dp * playerinmb(1,2)
    end if

    ! Calculate the difference in pressure layers (including buffer levels)
    do i = 1, (nzm-1)+buffer
      diff_pres_lvls(1,i) = playerinmb(1,i+1) - playerinmb(1,i)
    end do

    ts(1) = T_in_K(1,nzt+buffer)
!  Write a profile for Kurt's driver program for debugging purposes
!  write(time_char ,*) time
!  time_char =adjustl(time_char)
    !open(10, file="profile"//trim(time_char)//"dat")
    !write(10,'(2i4,a10)') nlen, (nz-1)+buffer, "TROPICAL"
    !do i=1, (nz-1)+buffer
    ! write(10,'(i4,9f12.6)') i, p_in_mb(1,i), playerinmb(1,i),T_in_K(1,i), &
    !sp_humidity(1,i), 100000.0*o3l(1,i), rcm_in_cloud_2d(1,i), &
    !rcil(1,i), cloud_frac_2d(1,i), diff_pres_lvls(1,i)
    !end do
    !write(10,'(a4,a12,3f12.6)') "","", playerinmb(1,nz+buffer), ts(1), amu0
    !close(10)

    ! Assertion check for unrealistic pressure
    do i = 2, size( p_in_mb(1,:) )
      if ( p_in_mb(1,i-1) > p_in_mb(1,i) ) then
        write(0,*) "Pressure (i)  [mb] = ", p_in_mb(1,i)
        write(0,*) "Pressure (i-1) [mb] = ", p_in_mb(1,i-1)
        write(0,*) "i level = ", i
        error stop "Fatal error: assertion check for pressure in BUGSrad_driver failed"
      end if
    end do

    call bugs_rad( nlen, slen, nzt+buffer, playerinmb,              &
                   p_in_mb, diff_pres_lvls, T_in_K, sp_humidity,         &
                   rcm_in_cloud_2d, rcil, rsm_2d, o3l,              &
                   ts, amu0, slr, alvdf,                               &
                   alndf, alvdr, alndr, sol_const,                     &
                   real( grav, kind=dp ), real( Cp, kind=dp ), radht_SW_2d, radht_LW_2d, &
                   Frad_dSW, Frad_uSW, Frad_dLW, Frad_uLW,             &
                   cloud_frac_2d )

    ! Michael pointed out that this was a temperature tendency, not a theta_l
    ! tendency.  The 2nd line should fix both.  -dschanen 28 July 2006
    radht_SW(1:nzt) = flip( real( radht_SW_2d(1,buffer+1:nzt+buffer), &
                                  kind = core_rknd ), nzt ) &
                      * ( 1.0_core_rknd / exner(1:nzt) )

    radht_LW(1:nzt) = flip( real( radht_LW_2d(1,buffer+1:nzt+buffer), &
                                  kind = core_rknd ), nzt ) &
                      * ( 1.0_core_rknd / exner(1:nzt) )

    radht = radht_SW + radht_LW

    Frad_SW_up = flip( real( Frad_uSW(1,buffer+1:nzm+buffer), kind = core_rknd ), nzm )

    Frad_LW_up = flip( real( Frad_uLW(1,buffer+1:nzm+buffer), kind = core_rknd ), nzm )

    Frad_SW_down = flip( real( Frad_dSW(1,buffer+1:nzm+buffer), kind = core_rknd ), nzm )

    Frad_LW_down = flip( real( Frad_dLW(1,buffer+1:nzm+buffer), kind = core_rknd ), nzm )

    Frad_SW(1:nzm) = Frad_SW_up - Frad_SW_down

    Frad_LW(1:nzm) = Frad_LW_up - Frad_LW_down

    Frad(1:nzm) = Frad_SW(1:nzm) + Frad_LW(1:nzm)

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
