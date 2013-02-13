! Interface specification for any microphysics scheme that we want to use latin
! hypercube sampling for.  New microphysics drivers will need to have the same
! number and type of arguments as this does.
interface
  subroutine microphys_sub( dt, nnzp, l_stats_sample, l_local_kk, &
                            l_latin_hypercube, thlm, wm, p_in_Pa, &
                            exner, rho, cloud_frac, pdf_params, w_std_dev, &
                            dzq, rcm, Ncm, s_mellor, rvm, Nc0, hydromet, &
                            hydromet_mc, hydromet_vel, &
                            rcm_mc, rvm_mc, thlm_mc, hydromet_vel_covar, &
                            wprtp_mc_tndcy, wpthlp_mc_tndcy, &
                            rtp2_mc_tndcy, thlp2_mc_tndcy, rtpthlp_mc_tndcy, &
                            rrainm_auto, rrainm_accr )

    use pdf_parameter_module, only: pdf_parameter

    use parameters_model, only: hydromet_dim

    use clubb_precision, only: core_rknd, time_precision

    real( kind = time_precision ), intent(in) :: dt ! Model timestep        [s]

    integer, intent(in) :: nnzp ! Points in the Vertical        [-]

    logical, intent(in) :: &
      l_stats_sample,   & ! Whether to accumulate statistics [T/F]
      l_local_kk,       & ! Whether we're using the local formulas
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
      Ncm,           & ! Cloud droplet number conc. in cloud [#/kg]
      s_mellor,      & ! The variable 's' from Mellor        [kg/kg]
      rvm              ! Vapor water mixing ratio            [kg/kg]

    ! Nc0 is a constant cloud droplet concentration NcV0 (#/m^3) / rho.
    ! This is normally used for cases that specify a constant cloud droplet
    ! concentration within cloud.  Set it to 0 when this is not required.
    ! Additionally, the flag l_const_Nc_in_cloud may have to be turned on to
    ! use a constant cloud droplet concentration within cloud. 
    real( kind = core_rknd ), dimension(nnzp), intent(in) :: &
      Nc0    ! Constant cloud droplet concentration within cloud  [#/kg] 

    real( kind = core_rknd ), dimension(nnzp,hydromet_dim), target, intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real( kind = core_rknd ), dimension(nnzp,hydromet_dim), target, intent(inout) :: &
      hydromet_mc, & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel   ! Hydrometeor sedimentation velocity [m/s]

    real( kind = core_rknd ), dimension(nnzp), intent(out) :: &
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

    real( kind = core_rknd ), dimension(nnzp,hydromet_dim), &
    target, intent(out) :: &
      hydromet_vel_covar    ! Covariance of V_xx and x_x (m-levs) [(m/s)(units)]

    real( kind = core_rknd ), dimension(nnzp), intent(out) :: &
      wprtp_mc_tndcy,   & ! Microphysics tendency for <w'rt'>   [m*(kg/kg)/s^2]
      wpthlp_mc_tndcy,  & ! Microphysics tendency for <w'thl'>  [m*K/s^2]
      rtp2_mc_tndcy,    & ! Microphysics tendency for <rt'^2>   [(kg/kg)^2/s]
      thlp2_mc_tndcy,   & ! Microphysics tendency for <thl'^2>  [K^2/s]
      rtpthlp_mc_tndcy    ! Microphysics tendency for <rt'thl'> [K*(kg/kg)/s]

    real( kind = core_rknd ), dimension(nnzp), intent(out) :: &
      rrainm_auto, & ! Autoconversion budget for rain water mixing ratio [kg/kg/s]
      rrainm_accr    ! Accretion budget for rain water mixing ratio      [kg/kg/s]

  end subroutine microphys_sub
end interface
