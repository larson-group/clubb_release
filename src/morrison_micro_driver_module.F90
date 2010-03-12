! $Id$
module morrison_micro_driver_module

  implicit none

  public :: morrison_micro_driver

  private

  contains
!-------------------------------------------------------------------------------
  subroutine morrison_micro_driver &
             ( dt, nnzp, l_stats_samp, l_local_kk, l_latin_hypercube, &
               thlm, p_in_Pa, exner, rho, pdf_params, &
               wm, w_std_dev, dzq, rcm, s_mellor, rvm, hydromet, hydromet_mc, &
               hydromet_vel, rcm_mc, rvm_mc, thlm_mc )
! Description:
!   Wrapper for the Morrison microphysics
! References:
!   None
!-------------------------------------------------------------------------------

    use variables_prognostic_module, only: pdf_parameter

    use parameters_model, only: hydromet_dim

    ! The version of the Morrison 2005 microphysics that is in SAM.
    use module_MP_graupel, only: &
      M2005MICRO_GRAUPEL  ! Procedure

    use stats_variables, only: &
      zt,  & ! Variables
      sfc

    use stats_variables, only: & 
      irsnowm_sd, & ! Variables
      iricem_sd, & 
      irrainm_sd, & 
      irsnowm_sd, &
      irgraupelm_sd

    use stats_variables, only: & 
      ieff_rad_cloud, & ! Variables
      ieff_rad_ice, &
      ieff_rad_snow, &
      ieff_rad_rain, &
      ieff_rad_graupel

    use stats_variables, only: & 
      imorr_rain_rate, & ! Variables
      imorr_snow_rate

    use stats_type, only:  & 
        stat_update_var, stat_update_var_pt  ! Procedure(s)

    use T_in_K_module, only: &
      T_in_K2thlm, & ! Procedure(s)
      thlm2T_in_K

    use array_index, only:  & 
      iirrainm, iirsnowm, iiricem, iirgraupelm, &
      iiNrm, iiNsnowm, iiNim, iiNgraupelm, iiNcm

    use constants, only: &
      sec_per_day, &
      zero_threshold

    implicit none

    ! External
    intrinsic :: max

    ! Input Variables
    real, intent(in) :: dt ! Model timestep        [s]

    integer, intent(in) :: nnzp ! Points in the Vertical        [-]

    logical, intent(in) :: &
      l_stats_samp,     & ! Whether to accumulate statistics [T/F]
      l_local_kk,       & ! Whether we're using the local formulas
      l_latin_hypercube   ! Whether we're using latin hypercube sampling

    real, dimension(nnzp), intent(in) :: &
      thlm,       & ! Liquid potential temperature       [K]
      p_in_Pa,    & ! Pressure                           [Pa]
      exner,      & ! Exner function                     [-]
      rho           ! Density on thermo. grid            [kg/m^3]

    type(pdf_parameter), intent(in) :: &
      pdf_params ! PDF parameters

    real, dimension(nnzp), intent(in) :: &
      wm, &        ! Mean w                     [m/s]
      w_std_dev, & ! Standard deviation of w    [m/s]
      dzq          ! Difference in heights      [m]

    real, dimension(nnzp), intent(in) :: &
      rcm,      & ! Liquid water mixing ratio        [kg/kg]
      s_mellor, & ! The variable 's' from Mellor     [kg/kg]
      rvm         ! Vapor water mixing ratio         [kg/kg]

    real, dimension(nnzp,hydromet_dim), target, intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    ! Input/Output Variables
    real, dimension(nnzp,hydromet_dim), target, intent(inout) :: &
      hydromet_mc,   & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel     ! Hydrometeor sedimentation velocity [m/s]

    ! Output Variables
    real, dimension(nnzp), intent(out) :: &
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

    ! Local Variables
    real, dimension(nnzp) :: & 
      effc, effi, effg, effs, effr ! Effective droplet radii [μ]

    real, dimension(nnzp) :: & 
      T_in_K,    & ! Temperature                        [K]
      T_in_K_mc, & ! Temperature tendency               [K/s]
      rcm_tmp,   & ! Temporary array for cloud water mixing ratio  [kg/kg]
      rvm_tmp,   & ! Temporary array for vapor water mixing ratio  [kg/kg]
      rcm_sten     ! Cloud dropet sedimentation tendency           [kg/kg/s]

    real, dimension(nnzp) :: & 
      cloud_frac  ! Cloud fraction   [-]

    real, dimension(nnzp,hydromet_dim) :: & 
      hydromet_sten, & ! Hydrometeor sedimentation tendency [(units vary)/s]
      hydromet_tmp     ! Temporary variable

    real, pointer, dimension(:,:) :: &
      dummy

    real :: Morr_snow_rate, Morr_rain_rate

    integer :: i

    ! ---- Begin Code ----
    ! Some dummy assignments to make compiler warnings go away...
    if ( .false. ) then
      dummy => hydromet
      dummy => hydromet_mc
      dummy => hydromet_vel
      cloud_frac = dummy(:,1)
      cloud_frac = s_mellor
      if ( l_local_kk ) stop
    end if

    ! Determine temperature
    T_in_K = thlm2T_in_K( thlm, exner, rcm )

    if ( .not. l_latin_hypercube ) then
      cloud_frac(1:nnzp) = max( zero_threshold, &
                        pdf_params%mixt_frac * pdf_params%cloud_frac1 &
                        + (1.-pdf_params%mixt_frac) * pdf_params%cloud_frac2 )
    else
      cloud_frac(1:nnzp) = 0.0
    end if

    rcm_tmp = rcm
    rvm_tmp = rvm

    do i = 1, hydromet_dim, 1
      hydromet_tmp(1:nnzp,i) = hydromet(1:nnzp,i)
    end do

    ! Initialize tendencies to zero
    T_in_K_mc(1:nnzp) = 0.0
    rcm_mc(1:nnzp) = 0.0
    rvm_mc(1:nnzp) = 0.0
    hydromet_mc(1:nnzp,:) = 0.0

    ! Call the Morrison microphysics
    call M2005MICRO_GRAUPEL &
         ( rcm_mc, hydromet_mc(:,iiricem), hydromet_mc(:,iirsnowm), &
           hydromet_mc(:,iirrainm), hydromet_mc(:,iiNcm), &
           hydromet_mc(:,iiNim), hydromet_mc(:,iiNsnowm), &
           hydromet_mc(:,iiNrm), rcm_tmp, hydromet_tmp(:,iiricem), &
           hydromet_tmp(:,iirsnowm), hydromet_tmp(:,iirrainm), hydromet_tmp(:,iiNcm), &
           hydromet_tmp(:,iiNim), hydromet_tmp(:,iiNsnowm), hydromet_tmp(:,iiNrm), &
           T_in_K_mc, rvm_mc, T_in_K, rvm_tmp, P_in_pa, rho, dzq, wm, w_std_dev, &
           Morr_rain_rate, Morr_snow_rate, effc, effi, effs, effr, dt, &
           1,1, 1,1, 1,nnzp, 1,1, 1,1, 1,nnzp, &
           hydromet_mc(:,iirgraupelm), hydromet_mc(:,iiNgraupelm), &
           hydromet_tmp(:,iirgraupelm), hydromet_tmp(:,iiNgraupelm), effg, &
           hydromet_sten(:,iirgraupelm), hydromet_sten(:,iirrainm), &
           hydromet_sten(:,iiricem), hydromet_sten(:,iirsnowm), &
           rcm_sten, cloud_frac )

    ! Update hydrometeor tendencies
    ! This done because the hydromet_mc arrays that are produced by
    ! M2005MICRO_GRAUPEL don't include the clipping term.
    do i = 1, hydromet_dim, 1
      hydromet_mc(:,i) = ( hydromet_tmp(:,i) - hydromet(:,i)  ) / real( dt )
    end do

    ! Update thetal based on absolute temperature
    thlm_mc = ( T_in_K2thlm( T_in_K, exner, rcm_tmp ) - thlm ) / real( dt )

    ! Sedimentation is handled within the Morrison microphysics
    hydromet_vel(:,:) = 0.0

    if ( l_stats_samp ) then

      ! -------- Sedimentation tendency from Morrison microphysics --------

      ! --- Mixing ratios ---

      call stat_update_var( irgraupelm_sd, hydromet_sten(:,iirgraupelm), zt )

      call stat_update_var( irrainm_sd, hydromet_sten(:,iirrainm), zt )

      call stat_update_var( irsnowm_sd, hydromet_sten(:,iirsnowm), zt )

      call stat_update_var( iricem_sd, hydromet_sten(:,iiricem), zt )

      ! --- Number concentrations ---
      ! No budgets for sedimentation are output

      ! Effective radii of hydrometeor species
      call stat_update_var( ieff_rad_cloud, effc(:), zt )
      call stat_update_var( ieff_rad_ice, effi(:), zt )
      call stat_update_var( ieff_rad_snow, effs(:), zt )
      call stat_update_var( ieff_rad_rain, effr(:), zt )
      call stat_update_var( ieff_rad_graupel, effg(:), zt )

      ! Snow and Rain rates at the bottom of the domain, in mm/day
      call stat_update_var_pt( imorr_rain_rate, 1, &
        Morr_rain_rate * real( sec_per_day ) / dt, sfc )

      call stat_update_var_pt( imorr_snow_rate, 1, &
        Morr_snow_rate * real(  sec_per_day ) / dt, sfc )
    end if ! l_stats_samp

    return
  end subroutine morrison_micro_driver

end module morrison_micro_driver_module
