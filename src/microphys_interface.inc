! Interface specification for any microphysics scheme that we want to use latin
! hypercube sampling for.  New microphysics drivers will need to have the same
! number and type of arguments as this does.
interface
  subroutine microphys_sub( dt, nnzp, l_sample, l_latin_hypercube, thlm, p_in_Pa, exner, rho, pdf_params, &
                            wm, w_std_dev, dzq, rcm, s_mellor, rvm, hydromet, hydromet_mc, &
                            hydromet_vel, rcm_mc, rvm_mc, thlm_mc )

    use variables_prognostic_module, only: pdf_parameter

    use parameters_model, only: hydromet_dim

    real, intent(in) :: dt ! Model timestep        [s]

    integer, intent(in) :: nnzp ! Points in the Vertical        [-]

    logical, intent(in) :: &
      l_sample,   &     ! Whether to accumulate statistics [T/F]
      l_latin_hypercube ! Whether we're using latin hypercube sampling

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
      dzq          ! Difference in height       [m]

    real, dimension(nnzp), intent(in) :: &
      rcm,      & ! Liquid water mixing ratio        [kg/kg]
      s_mellor, & ! The variable 's' from Mellor     [kg/kg]
      rvm         ! Vapor water mixing ratio         [kg/kg]

    real, dimension(nnzp,hydromet_dim), target, intent(in) :: &
      hydromet ! Hydrometeor species    [units vary]

    real, dimension(nnzp,hydromet_dim), target, intent(inout) :: &
      hydromet_mc, & ! Hydrometeor time tendency          [(units vary)/s]
      hydromet_vel   ! Hydrometeor sedimentation velocity [m/s]

    real, dimension(nnzp), intent(out) :: &
      rcm_mc, & ! Time tendency of liquid water mixing ratio    [kg/kg/s]
      rvm_mc, & ! Time tendency of vapor water mixing ratio     [kg/kg/s]
      thlm_mc   ! Time tendency of liquid potential temperature [K/s]

  end subroutine microphys_sub
end interface
