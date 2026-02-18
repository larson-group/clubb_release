!------------------------------------------------------------------------
! $Id$
!===============================================================================
! Interface specification for any microphysics scheme that we want to use latin
! hypercube sampling for.  New microphysics drivers will need to have the same
! number and type of arguments as this does.
interface
  subroutine microphys_sub( gr, dt, nnzp, &
                            hydromet_dim, hm_metadata, &
                            l_latin_hypercube, &
                            thlm, wm, p_in_Pa, &
                            exner, rho, cloud_frac, w_std_dev, &
                            dzq, rcm, Nc, chi, rvm, hydromet, &
                            saturation_formula, &
                            stats, icol,         &
                            hydromet_mc, hydromet_vel, Ncm_mc, &
                            rcm_mc, rvm_mc, thlm_mc, &
                            microphys_stats_zt, microphys_stats_sfc )

    use clubb_precision, only: core_rknd

    use grid_class, only: grid ! Type

    use stats_netcdf, only: &
      stats_type

    use corr_varnce_module, only: &
        hm_metadata_type

    use microphys_stats_vars_module, only: &
      microphys_stats_vars_type   ! Type

    ! Input Variables

    type(grid), intent(in) :: gr

    real( kind = core_rknd ), intent(in) :: dt ! Model timestep        [s]

    integer, intent(in) :: &
      nnzp, &       ! Points in the Vertical        [-]
      hydromet_dim

    type (hm_metadata_type), intent(in) :: &
      hm_metadata

    logical, intent(in) :: &
      l_latin_hypercube   ! Whether we're using latin hypercube sampling

    real( kind = core_rknd ), dimension(nnzp), intent(in) :: &
      thlm,       & ! Liquid potential temperature       [K]
      wm,         & ! Mean vertical velocity             [m/s]
      p_in_Pa,    & ! Pressure                           [Pa]
      exner,      & ! Exner function                     [-]
      rho,        & ! Density on thermodynamic levels    [kg/m^3]
      cloud_frac    ! Cloud fraction                     [-]

    real( kind = core_rknd ), dimension(nnzp), intent(in) :: &
      w_std_dev, & ! Standard deviation of w    [m/s]
      dzq          ! Difference in height       [m]

    real( kind = core_rknd ), dimension(nnzp), intent(in) :: &
      rcm,           & ! Cloud water mixing ratio            [kg/kg]
      Nc,            & ! Cloud droplet number conc. in cloud [#/kg]
      chi,           & ! The variable 's' from Mellor        [kg/kg]
      rvm              ! Vapor water mixing ratio            [kg/kg]

    real( kind = core_rknd ), dimension(nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    integer, intent(in) :: &
      saturation_formula ! Integer that stores the saturation formula to be used

    type(stats_type), intent(inout) :: &
      stats

    integer, intent(in) :: &
      icol

    ! Output Variables

    real( kind = core_rknd ), dimension(nnzp,hydromet_dim), intent(out) :: &
      hydromet_mc, & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel   ! Hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(nnzp), intent(out) :: &
      Ncm_mc, & ! Time tendency of cloud droplet concentration  [num/kg/s]
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

    type(microphys_stats_vars_type), intent(out) :: &
      microphys_stats_zt, &   ! Output variables for statistical sampling
      microphys_stats_sfc

  end subroutine microphys_sub
end interface
