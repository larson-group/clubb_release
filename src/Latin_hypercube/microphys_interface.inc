! Interface specification for any microphysics scheme that we want to use latin
! hypercube sampling for.  New microphysics drivers will need to have the same
! number and type of arguments as this does.
interface
  subroutine microphys_sub( dt, nnzp, l_stats_sample, &
                            l_latin_hypercube, thlm, wm, p_in_Pa, &
                            exner, rho, cloud_frac, pdf_params, w_std_dev, &
                            dzq, rcm, Nc, s_mellor, rvm, hydromet, &
                            hydromet_mc, hydromet_vel, &
                            rcm_mc, rvm_mc, thlm_mc, &
                            rtp2_mc, thlp2_mc, &
                            wprtp_mc, wpthlp_mc, &
                            rtpthlp_mc, & 
                            rrainm_auto, rrainm_accr, rrainm_evap, &
                            Nrm_auto, Nrm_evap )

    use pdf_parameter_module, only: pdf_parameter

    use parameters_model, only: hydromet_dim

    use clubb_precision, only: core_rknd, time_precision

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

    type(pdf_parameter), target, dimension(nnzp), intent(in) :: &
      pdf_params ! PDF parameters

    real( kind = core_rknd ), dimension(nnzp), intent(in) :: &
      w_std_dev, & ! Standard deviation of w    [m/s]
      dzq          ! Difference in height       [m]

    real( kind = core_rknd ), dimension(nnzp), intent(in) :: &
      rcm,           & ! Cloud water mixing ratio            [kg/kg]
      Nc,           & ! Cloud droplet number conc. in cloud [#/kg]
      s_mellor,      & ! The variable 's' from Mellor        [kg/kg]
      rvm              ! Vapor water mixing ratio            [kg/kg]

    real( kind = core_rknd ), dimension(nnzp,hydromet_dim), target, intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nnzp,hydromet_dim), target, intent(out) :: &
      hydromet_mc, & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel   ! Hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(nnzp), intent(out) :: &
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(nnzp), intent(out) :: &
      rtp2_mc,  &  ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc, &  ! Microphysics tendency for <thl'^2>  [K^2/s]
      wprtp_mc, &  ! Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
      wpthlp_mc, & ! Microphysics tendency for <w'thl'>  [m*K/s^2]
      rtpthlp_mc   ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]

    real( kind = core_rknd ), dimension(nnzp), intent(out) :: &
      rrainm_auto, & ! Autoconversion budget for rain water mixing ratio [kg/kg/s]
      rrainm_accr, & ! Accretion budget for rain water mixing ratio      [kg/kg/s]
      rrainm_evap, & ! Evaporation budget for rain water mixing ratio    [kg/kg/s]
      Nrm_auto,    & ! Change in Nrm due to autoconversion               [num/kg/s]
      Nrm_evap       ! Change in Nrm due to evaporation                  [num/kg/s]

  end subroutine microphys_sub
end interface
