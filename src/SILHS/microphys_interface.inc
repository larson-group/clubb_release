! Interface specification for any microphysics scheme that we want to use latin
! hypercube sampling for.  New microphysics drivers will need to have the same
! number and type of arguments as this does.
interface
  subroutine microphys_sub( dt, nnzp, l_stats_sample, &
                            l_latin_hypercube, thlm, wm, p_in_Pa, &
                            exner, rho, cloud_frac, w_std_dev, &
                            dzq, rcm, Nc, s_mellor, rvm, hydromet, lh_stat_sample_weight, &
                            hydromet_mc, hydromet_vel, Ncm_mc, &
                            rcm_mc, rvm_mc, thlm_mc, &
                            rrainm_auto, rrainm_accr, rrainm_evap, &
                            Nrm_auto, Nrm_evap, microphys_stats_zt, microphys_stats_sfc )

    use pdf_parameter_module, only: pdf_parameter

    use parameters_model, only: hydromet_dim

    use clubb_precision, only: core_rknd, time_precision

    use microphys_stats_vars_module, only: &
      microphys_stats_vars_type   ! Type

    ! Input Variables

    real( kind = time_precision ), intent(in) :: dt ! Model timestep        [s]

    integer, intent(in) :: nnzp ! Points in the Vertical        [-]

    logical, intent(in) :: &
      l_stats_sample,   & ! Whether to accumulate statistics [T/F]
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
      Nc,           & ! Cloud droplet number conc. in cloud [#/kg]
      s_mellor,      & ! The variable 's' from Mellor        [kg/kg]
      rvm              ! Vapor water mixing ratio            [kg/kg]

    real( kind = core_rknd ), dimension(nnzp,hydromet_dim), intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), intent(in) :: &
      lh_stat_sample_weight ! The weight of this sample, if using SILHS.

    ! Output Variables

    real( kind = core_rknd ), dimension(nnzp,hydromet_dim), intent(out) :: &
      hydromet_mc, & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel   ! Hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(nnzp), intent(out) :: &
      Ncm_mc, & ! Time tendency of cloud droplet concentration  [num/kg/s]
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(nnzp), intent(out) :: &
      rrainm_auto, & ! Autoconversion budget for rain water mixing ratio [kg/kg/s]
      rrainm_accr, & ! Accretion budget for rain water mixing ratio      [kg/kg/s]
      rrainm_evap, & ! Evaporation budget for rain water mixing ratio    [kg/kg/s]
      Nrm_auto,    & ! Change in Nrm due to autoconversion               [num/kg/s]
      Nrm_evap       ! Change in Nrm due to evaporation                  [num/kg/s]

    type(microphys_stats_vars_type), intent(out) :: &
      microphys_stats_zt, &   ! Output variables for statistical sampling
      microphys_stats_sfc

  end subroutine microphys_sub
end interface