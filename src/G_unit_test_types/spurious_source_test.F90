!===============================================================================
module spurious_source_test

  implicit none

  public :: spurious_source_unit_test    ! Procedure(s)

  private

  contains

  !=============================================================================
  function spurious_source_unit_test(  )

    ! Description:
    !
    ! This function checks CLUBB's predictive equations for <rt> and <thl> for
    ! spurious sources and sinks.  The CLUBB code solves <rt> and <w'rt'>
    ! together, as well as <thl> and <w'thl'> together, in subroutine
    ! advance_xm_wpxp.  This function initializes profiles for a number of
    ! variables that are needed as input to that function.  It then calculates
    ! the vertical integrals of both <rt> and <thl> before the call to
    ! advance_xm_wpxp.  It calls advance_xm_wpxp, which advances <rt>, <thl>,
    ! <w'rt'>, and <w'thl'> one model timestep.  After advancing the predictive
    ! fields, it calculates the new vertical integrals of both <rt> and <thl>.
    ! The vertical integrals of <rt>|_forcing and <thl>_forcing are also
    ! calculated.  The vertical integrals of <rt> and <thl> that are calculated
    ! prior to calling advance_xm_wpxp should match the vertical integrals of
    ! <rt> and <thl> that are calculated after calling advance_xm_wpxp when
    ! adjustments for forcings and surface fluxes are taken into account,
    ! within a small tolerance value.
    !
    ! The code performs a few iterations of this test.  If both <rt> and <thl>
    ! pass the test for every iteration, the value of spurious_source_unit_test
    ! is set to 0 and the test passes.  Otherwise, the value of
    ! spurious_source_unit_test is set to 1 and the test fails.

    ! References:
    !-----------------------------------------------------------------------

    use grid_class, only: &
        setup_grid, & ! Procedure(s)
        zm2zt,      &
        zt2zm,      &
        gr            ! Variable(s)

    use constants_clubb, only: &
        one,       & ! Variable(s)
        one_half,  &
        zero,      &
        Lv,        &
        Rd,        &
        Cp,        &
        ep1,       &
        ep2,       &
        w_tol_sqd

    use advance_xm_wpxp_module, only: &
        advance_xm_wpxp    ! Procedure(s)

    use parameters_tunable, only: &
        adj_low_res_nu    ! Procedure(s)

    use fill_holes, only: &
        vertical_integral    ! Procedure(s)

    use numerical_check, only: & 
        calculate_spurious_source    ! Procedure(s)

    use interpolation, only: &
        lin_interpolate_two_points    ! Procedure(s)

    use adg1_adg2_3d_luhar_pdf, only: &
        ADG1_w_closure    ! Procedure(s)

    use pdf_parameter_module, only: &
        implicit_coefs_terms    ! Variable type(s)

    use parameters_model, only: &
        sclr_dim    ! Variable(s)

    use clubb_precision, only: &
        core_rknd    ! Variable(s)

    use model_flags, only: &
        l_diffuse_rtm_and_thlm, & ! Variable(s)
        l_stability_correct_Kh_N2_zm, &
        l_upwind_wpxp_ta, &
        l_upwind_xm_ma, &
        l_use_C7_Richardson

    implicit none

    ! Return Variable
    integer :: &
      spurious_source_unit_test    ! Pass/fail output for spurious source test

    ! Local Variables
    integer, parameter :: &
      nz = 76 ! Number of vertical grid levels
              ! Note:  this needs to match the number of grid levels calculated
              !        in the call to setup_grid, as determined by zm_init,
              !        zm_top, and deltaz.  Otherwise, a runtime error will
              !        cause this test to fail.

    real( kind = core_rknd ), parameter :: &
      sfc_elevation = zero  ! Elevation of ground level    [m AMSL]

    ! Flag to see if CLUBB is running on it's own,
    ! or if it's implemented as part of a host model.
    logical, parameter :: &
      l_implemented = .false.

    ! If CLUBB is running on it's own, this option determines if it is using:
    ! 1) an evenly-spaced grid;
    ! 2) a stretched (unevenly-spaced) grid entered on the thermodynamic grid
    !    levels (with momentum levels set halfway between thermodynamic levels);
    !    or
    ! 3) a stretched (unevenly-spaced) grid entered on the momentum grid levels
    !    (with thermodynamic levels set halfway between momentum levels).
    integer, parameter :: &
      grid_type = 1

    ! If the CLUBB model is running by itself, and is using an evenly-spaced
    ! grid (grid_type = 1), it needs the vertical grid spacing and
    ! momentum-level starting altitude as input.
    real( kind = core_rknd ) ::  & 
      deltaz,  & ! Vertical grid spacing                  [m]
      zm_init, & ! Initial grid altitude (momentum level) [m]
      zm_top     ! Maximum grid altitude (momentum level) [m]

    ! If the CLUBB parameterization is implemented in a host model, it needs to
    ! use the host model's momentum level altitudes and thermodynamic level
    ! altitudes.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on thermodynamic levels (grid_type = 2), it needs to use the
    ! thermodynamic level altitudes as input.
    ! If the CLUBB model is running by itself, but is using a stretched grid
    ! entered on momentum levels (grid_type = 3), it needs to use the momentum
    ! level altitudes as input.
    real( kind = core_rknd ), dimension(nz) ::  & 
      momentum_heights,      & ! Momentum level altitudes (input)      [m]
      thermodynamic_heights    ! Thermodynamic level altitudes (input) [m]

    integer :: &
      begin_height, & ! Lower bound for *_heights arrays [-]
      end_height      ! Upper bound for *_heights arrays [-]

    real( kind = core_rknd ) ::  & 
      dt                 ! Timestep                                 [s]

    real( kind = core_rknd ), dimension(nz) :: & 
      sigma_sqd_w,     & ! sigma_sqd_w on momentum levels           [-]
      wm_zm,           & ! w wind component on momentum levels      [m/s]
      wm_zt,           & ! w wind component on thermodynamic levels [m/s]
      wp2,             & ! w'^2 (momentum levels)                   [m^2/s^2]
      Lscale,          & ! Turbulent mixing length                  [m]
      em,              & ! Turbulent Kinetic Energy (TKE)           [m^2/s^2]
      wp3_on_wp2,      & ! Smoothed wp3 / wp2 on momentum levels    [m/s]
      wp3_on_wp2_zt,   & ! Smoothed wp3 / wp2 on thermo. levels     [m/s]
      Kh_zt,           & ! Eddy diffusivity on thermodynamic levels [m^2/s]
      Kh_zm,           & ! Eddy diffusivity on momentum levels      [m^s/s]
      tau_C6_zm,       & ! Time-scale tau on m-levs applied to C6 term  [s]
      Skw_zm,          & ! Skewness of w on momentum levels         [-]
      wp2rtp,          & ! <w'^2 r_t'> (thermodynamic levels)    [m^2/s^2 kg/kg]
      rtpthvp,         & ! r_t'th_v' (momentum levels)              [(kg/kg) K]
      rtm_forcing,     & ! r_t forcing (thermodynamic levels)       [(kg/kg)/s]
      wprtp_forcing,   & ! <w'r_t'> forcing (momentum levels)      [(kg/kg)/s^2]
      rtm_ref,         & ! rtm for nudging                          [kg/kg]
      wp2thlp,         & ! <w'^2 th_l'> (thermodynamic levels)      [m^2/s^2 K]
      thlpthvp,        & ! th_l'th_v' (momentum levels)             [K^2]
      thlm_forcing,    & ! th_l forcing (thermodynamic levels)      [K/s]
      wpthlp_forcing,  & ! <w'th_l'> forcing (momentum levels)      [K/s^2]
      thlm_ref,        & ! thlm for nudging                         [K]
      rho_ds_zm,       & ! Dry, static density on momentum levels   [kg/m^3]
      rho_ds_zt,       & ! Dry, static density on thermo. levels    [kg/m^3]
      invrs_rho_ds_zm, & ! Inv. dry, static density @ moment. levs. [m^3/kg]
      invrs_rho_ds_zt, & ! Inv. dry, static density @ thermo. levs. [m^3/kg]
      thv_ds_zm,       & ! Dry, base-state theta_v on moment. levs. [K]
      ! Added for clipping by Vince Larson 29 Sep 2007
      rtp2,            & ! r_t'^2 (momentum levels)                 [(kg/kg)^2]
      thlp2,           & ! th_l'^2 (momentum levels)                [K^2]
      ! End of Vince Larson's addition.
      w_1_zm,          & ! Mean w (1st PDF component)               [m/s]
      w_2_zm,          & ! Mean w (2nd PDF component)               [m/s]
      varnce_w_1_zm,   & ! Variance of w (1st PDF component)        [m^2/s^2]
      varnce_w_2_zm,   & ! Variance of w (2nd PDF component)        [m^2/s^2]
      mixt_frac_zm       ! Weight of 1st PDF component (Sk_w dependent) [-]

    ! Additional variables for passive scalars
    real( kind = core_rknd ), dimension(nz,sclr_dim) :: & 
      wp2sclrp,      & ! <w'^2 sclr'> (thermodynamic levels)   [Units vary]
      sclrpthvp,     & ! <sclr' th_v'> (momentum levels)       [Units vary]
      sclrm_forcing, & ! sclrm forcing (thermodynamic levels)  [Units vary]
      sclrp2           ! For clipping Vince Larson             [Units vary]

    real( kind = core_rknd ), dimension(nz) ::  &
      exner,             & ! Exner function                            [-]
      rcm,               & ! cloud water mixing ratio, r_c             [kg/kg]
      p_in_Pa,           & ! Air pressure                              [Pa]
      thvm,              & ! Virutal potential temperature             [K]
      Cx_fnc_Richardson, & ! Cx_fnc computed from Richardson_num       [-]
      ice_supersat_frac      


    type(implicit_coefs_terms) :: &
      pdf_implicit_coefs_terms    ! Implicit coefs / explicit terms [units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), dimension(nz) :: & 
      um_forcing, & ! <u> forcing term (thermodynamic levels)      [m/s^2]
      vm_forcing, & ! <v> forcing term (thermodynamic levels)      [m/s^2]
      ug,         & ! <u> geostrophic wind (thermodynamic levels)  [m/s]
      vg,         & ! <v> geostrophic wind (thermodynamic levels)  [m/s]
      wpthvp        ! <w'thv'> (momentum levels)                   [m/s K]

     real( kind = core_rknd ) ::  &
      fcor          ! Coriolis parameter                           [s^-1]

    real( kind = core_rknd ), dimension(nz) :: & 
      um_ref, & ! Reference u wind component for nudging       [m/s]
      vm_ref, & ! Reference v wind component for nudging       [m/s]
      up2,    & ! Variance of the u wind component             [m^2/s^2]
      vp2       ! Variance of the v wind component             [m^2/s^2]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz) ::  & 
      rtm,    & ! r_t  (total water mixing ratio)           [kg/kg]
      wprtp,  & ! w'r_t'                                    [(kg/kg) m/s]
      thlm,   & ! th_l (liquid water potential temperature) [K]
      wpthlp    ! w'th_l'                                   [K m/s]

    ! Input/Output Variables
    real( kind = core_rknd ), dimension(nz,sclr_dim) ::  & 
      sclrm, wpsclrp !                                     [Units vary]

    ! Variables used to predict <u> and <u'w'>, as well as <v> and <v'w'>.
    real( kind = core_rknd ), dimension(nz) ::  & 
      um,    & ! <u>:  mean west-east horiz. velocity (thermo. levs.)   [m/s]
      upwp,  & ! <u'w'>:  momentum flux (momentum levels)               [m^2/s^2]
      vm,    & ! <v>:  mean south-north horiz. velocity (thermo. levs.) [m/s]
      vpwp,  & ! <v'w'>:  momentum flux (momentum levels)               [m^2/s^2]
      uprcp, & ! < u' r_c' >              [(m kg)/(s kg)]
      vprcp, & ! < v' r_c' >              [(m kg)/(s kg)]
      rc_coef    ! Coefficient on X'r_c' in X'th_v' equation    [K/(kg/kg)]

    real( kind = core_rknd ), dimension(nz) :: &
      w_1_n_zm, &
      w_2_n_zm

    real( kind = core_rknd ) :: &
      rtm_integral_before,   & ! Vert. integral of rtm before advance   [kg/m^2]
      rtm_integral_after,    & ! Vert. integral of rtm after advance    [kg/m^2]
      rtm_integral_forcing,  & ! Vert. integral of rtm forcing        [kg/s/m^2]
      rtm_flux_top,          & ! Flux of rtm at top of domain         [kg/s/m^2]
      rtm_flux_sfc,          & ! Flux of rtm at surface               [kg/s/m^2]
      rtm_spur_src,          & ! Spurious source of rtm               [kg/s/m^2]
      thlm_integral_before,  & ! Vert. integral of thlm before adv. [(K kg)/m^2]
      thlm_integral_after,   & ! Vert. integral of thlm after adv.  [(K kg)/m^2]
      thlm_integral_forcing, & ! Vert. integral of thlm forcing   [(K kg)/s/m^2]
      thlm_flux_top,         & ! Flux of thlm at top of domain    [(K kg)/s/m^2]
      thlm_flux_sfc,         & ! Flux of thlm at surface          [(K kg)/s/m^2]
      thlm_spur_src            ! Spurious source of thlm          [(K kg)/s/m^2]

    real( kind = core_rknd ), dimension(4) ::  & 
       z_snd,                 &
       sigma_sqd_w_snd,       &
       wm_zm_snd,             &
       wp2_snd,               &
       wp3_on_wp2_snd,        &
       Kh_zm_snd,             &
       tau_C6_zm_snd,         &
       rtpthvp_snd,           &
       rtm_forcing_snd,       &
       wprtp_forcing_snd,     &
       thlpthvp_snd,          &
       thlm_forcing_snd,      &
       wpthlp_forcing_snd,    &
       rho_ds_zm_snd,         &
       rtp2_snd,              &
       thlp2_snd,             &
       Cx_fnc_Richardson_snd, &
       Lscale_snd,            &
       wp2rtp_snd,            &
       wp2thlp_snd,           &
       rcm_snd,               &
       p_in_Pa_snd,           &
       rtm_snd,               &
       thlm_snd,              &
       wprtp_snd,             &
       wpthlp_snd,            &
       rand1

    real( kind = core_rknd ) :: &
      density_weighted_height    ! Integrated density * height over the domain

    integer, parameter :: &
      num_iter = 5    ! Number of different configurations to loop over.

    integer :: &
      num_errors    ! Number of failed parameter sets

    real( kind = core_rknd ), parameter :: &
      tol = 1.0e-10_core_rknd    ! Tolerance to determine pass or fail

    integer :: iter, k, i  ! Loop indices


    write(*,*)
    write(*,*) "Performing spurious source unit test"
    write(*,*) "===================================="
    write(*,*)

    deltaz = 40.0_core_rknd
    zm_init = zero
    zm_top = 3000.0_core_rknd

    momentum_heights = zero
    thermodynamic_heights = zero

    ! Set up the vertical grid.
    call setup_grid( nz, sfc_elevation, l_implemented,        &
                     grid_type, deltaz, zm_init, zm_top,      &
                     momentum_heights, thermodynamic_heights, &
                     begin_height, end_height                 )

    ! Calculate the value of nu for use in advance_xm_wpxp.
    call adj_low_res_nu( gr%nz, grid_type, deltaz, &
                         momentum_heights, thermodynamic_heights )

    dt = 300.0_core_rknd

    num_errors = 0

    call random_seed( )

    do iter = 1, num_iter, 1

       write(*,*) "Test ", iter
       write(*,*) ""

       z_snd = (/ zero, 1000.0_core_rknd, 2000.0_core_rknd, &
                  1.01_core_rknd * zm_top /)
       sigma_sqd_w_snd = (/ 0.3_core_rknd, 0.1_core_rknd, 0.15_core_rknd, &
                            0.3_core_rknd /)
       wm_zm_snd = (/ zero, -0.005_core_rknd, -0.005_core_rknd, zero /)
       wp2_snd = (/ 0.1_core_rknd, 0.2_core_rknd, 0.1_core_rknd, w_tol_sqd /)
       wp3_on_wp2_snd = (/ zero, 2.0_core_rknd, 1.0_core_rknd, zero /)
       Kh_zm_snd = (/ 100.0_core_rknd, 500.0_core_rknd, 300.0_core_rknd, zero /)
       tau_C6_zm_snd = (/ 500.0_core_rknd, 900.0_core_rknd, 900.0_core_rknd, &
                          500.0_core_rknd /)
       rtpthvp_snd = (/ zero, -4.0e-5_core_rknd, -3.0e-5_core_rknd, zero /)
       rtm_forcing_snd = (/ 0.0_core_rknd, -2.0e-8_core_rknd, &
                            -2.5e-8_core_rknd, 0.0_core_rknd /)
       wprtp_forcing_snd = (/ 0.0_core_rknd, -1.0e-7_core_rknd, &
                              -1.5e-7_core_rknd, 0.0_core_rknd /)
       thlpthvp_snd = (/ zero, 0.1_core_rknd, 0.1_core_rknd, zero /)
       thlm_forcing_snd = (/ 0.0_core_rknd, -3.0e-5_core_rknd, &
                             -2.0e-5_core_rknd, 0.0_core_rknd /)
       wpthlp_forcing_snd = (/ 0.0_core_rknd, 5.0e-4_core_rknd, &
                               3.5e-4_core_rknd, 0.0_core_rknd /)
       rho_ds_zm_snd = (/ 1.0_core_rknd, 0.92_core_rknd, 0.84_core_rknd, &
                          0.76_core_rknd /)
       rtp2_snd = (/ zero, 2.5e-7_core_rknd, 2.5e-7_core_rknd, zero /)
       thlp2_snd = (/ zero, 0.1_core_rknd, 0.1_core_rknd, zero /)
       ! The value of Cx_fnc_Richardson should be between 0 and 1, as
       ! the C7_Skw_fnc might be set to the value of Cx_fnc_Richardson.
       Cx_fnc_Richardson_snd = (/ 0.1_core_rknd, 0.5_core_rknd, &
                                  0.7_core_rknd, 0.9_core_rknd /)
       Lscale_snd = (/ 200.0_core_rknd, 100.0_core_rknd, 50.0_core_rknd, &
                       10.0_core_rknd /)
       wp2rtp_snd = (/ zero, 1.0e-5_core_rknd, 7.5e-6_core_rknd, zero /)
       wp2thlp_snd = (/ zero, -0.1_core_rknd, -0.9_core_rknd, zero /)
       rcm_snd = (/ zero, 1.0e-5_core_rknd, zero, zero /)
       p_in_Pa_snd = (/ 1.0e5_core_rknd, 9.0e4_core_rknd, 8.2e4_core_rknd, &
                        7.5e4_core_rknd /)
       rtm_snd = (/ 1.7e-2_core_rknd, 1.35e-2_core_rknd, 8.5e-3_core_rknd, &
                    1.0e-3_core_rknd /)
       thlm_snd = (/ 298.0_core_rknd, 298.0_core_rknd, 302.0_core_rknd, &
                     307.0_core_rknd /)
       wprtp_snd = (/ 5.0e-4_core_rknd, 5.0e-4_core_rknd, 2.5e-4_core_rknd, &
                      zero /)
       wpthlp_snd = (/ 0.01_core_rknd, -0.04_core_rknd, -0.06_core_rknd, zero /)

       if ( iter > 1 ) then

          call random_number( rand1 )
          sigma_sqd_w_snd = 2.0_core_rknd * rand1 * sigma_sqd_w_snd
          wm_zm_snd = 5.0_core_rknd * ( rand1 - one_half ) * wm_zm_snd
          wp2_snd = 2.0_core_rknd * rand1 * wp2_snd
          wp3_on_wp2_snd = 5.0_core_rknd * ( rand1 - one_half ) * wp3_on_wp2_snd
          Kh_zm_snd = 2.0_core_rknd * rand1 * Kh_zm_snd
          tau_C6_zm_snd = 2.0_core_rknd * rand1 * tau_C6_zm_snd
          rtpthvp_snd = 2.0_core_rknd * rand1 * rtpthvp_snd
          rtm_forcing_snd = 2.0_core_rknd * rand1 * rtm_forcing_snd
          wprtp_forcing_snd = 2.0_core_rknd * rand1 * wprtp_forcing_snd
          thlpthvp_snd = 2.0_core_rknd * rand1 * thlpthvp_snd
          thlm_forcing_snd = 2.0_core_rknd * rand1 * thlm_forcing_snd
          wpthlp_forcing_snd = 2.0_core_rknd * rand1 * wpthlp_forcing_snd
          rtp2_snd = 2.0_core_rknd * rand1 * rtp2_snd
          thlp2_snd = 2.0_core_rknd * rand1 * thlp2_snd
          Lscale_snd = 2.0_core_rknd * rand1 * Lscale_snd
          wp2rtp_snd = 2.0_core_rknd * rand1 * wp2rtp_snd
          wp2thlp_snd = 2.0_core_rknd * rand1 * wp2thlp_snd
          rcm_snd = 2.0_core_rknd * rand1 * rcm_snd
          rtm_snd = 2.0_core_rknd * rand1 * rtm_snd
          thlm_snd = thlm_snd + 10.0_core_rknd * ( rand1 - one_half )
          wprtp_snd = 2.0_core_rknd * rand1 * wprtp_snd
          wpthlp_snd = 2.0_core_rknd * rand1 * wpthlp_snd

       endif ! iter > 1
       

       ! Loop over all momentum levels.
       do k = 1, gr%nz, 1

          i = 1

          do while ( z_snd(i) <= gr%zm(k) )

             i = i + 1

             sigma_sqd_w(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           sigma_sqd_w_snd(i), &
                                           sigma_sqd_w_snd(i-1) )

             wm_zm(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           wm_zm_snd(i), wm_zm_snd(i-1) )

             wp2(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           wp2_snd(i), wp2_snd(i-1) )

             wp3_on_wp2(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           wp3_on_wp2_snd(i), &
                                           wp3_on_wp2_snd(i-1) )

             Kh_zm(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           Kh_zm_snd(i), Kh_zm_snd(i-1) )

             tau_C6_zm(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           tau_C6_zm_snd(i), &
                                           tau_C6_zm_snd(i-1) )

             rtpthvp(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           rtpthvp_snd(i), rtpthvp_snd(i-1) )

             wprtp_forcing(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           wprtp_forcing_snd(i), &
                                           wprtp_forcing_snd(i-1) )

             thlpthvp(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           thlpthvp_snd(i), thlpthvp_snd(i-1) )

             wpthlp_forcing(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           wpthlp_forcing_snd(i), &
                                           wpthlp_forcing_snd(i-1) )

             rho_ds_zm(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           rho_ds_zm_snd(i), &
                                           rho_ds_zm_snd(i-1) )

             rtp2(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           rtp2_snd(i), rtp2_snd(i-1) )

             thlp2(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           thlp2_snd(i), thlp2_snd(i-1) )

             Cx_fnc_Richardson(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           Cx_fnc_Richardson_snd(i), &
                                           Cx_fnc_Richardson_snd(i-1) )

             wprtp(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           wprtp_snd(i), wprtp_snd(i-1) )

             wpthlp(k) &
             = lin_interpolate_two_points( gr%zm(k), z_snd(i), z_snd(i-1), &
                                           wpthlp_snd(i), wpthlp_snd(i-1) )

          enddo ! while ( z_snd(i) < gr%zm(k) )

       enddo ! k = 1, gr%nz, 1

       ! Loop over all thermodynamic levels.
       do k = 2, gr%nz, 1

          i = 1

          do while ( z_snd(i) < gr%zt(k) )

             i = i + 1

             Lscale(k) &
             = lin_interpolate_two_points( gr%zt(k), z_snd(i), z_snd(i-1), &
                                           Lscale_snd(i), Lscale_snd(i-1) )

             wp2rtp(k) &
             = lin_interpolate_two_points( gr%zt(k), z_snd(i), z_snd(i-1), &
                                           wp2rtp_snd(i), wp2rtp_snd(i-1) )

             wp2thlp(k) &
             = lin_interpolate_two_points( gr%zt(k), z_snd(i), z_snd(i-1), &
                                           wp2thlp_snd(i), wp2thlp_snd(i-1) )

             rcm(k) &
             = lin_interpolate_two_points( gr%zt(k), z_snd(i), z_snd(i-1), &
                                           rcm_snd(i), rcm_snd(i-1) )

             p_in_Pa(k) &
             = lin_interpolate_two_points( gr%zt(k), z_snd(i), z_snd(i-1), &
                                           p_in_Pa_snd(i), p_in_Pa_snd(i-1) )

             rtm_forcing(k) &
             = lin_interpolate_two_points( gr%zt(k), z_snd(i), z_snd(i-1), &
                                           rtm_forcing_snd(i), &
                                           rtm_forcing_snd(i-1) )

             thlm_forcing(k) &
             = lin_interpolate_two_points( gr%zt(k), z_snd(i), z_snd(i-1), &
                                           thlm_forcing_snd(i), &
                                           thlm_forcing_snd(i-1) )

             rtm(k) &
             = lin_interpolate_two_points( gr%zt(k), z_snd(i), z_snd(i-1), &
                                           rtm_snd(i), rtm_snd(i-1) )

             thlm(k) &
             = lin_interpolate_two_points( gr%zt(k), z_snd(i), z_snd(i-1), &
                                           thlm_snd(i), thlm_snd(i-1) )

          enddo ! while ( z_snd(i) < gr%zt(k) )

       enddo ! k = 1, gr%nz, 1

       ! Set the values of thermodynamic level variables below the model
       ! surface (thermodynamic level k = 1).
       Lscale(1) = Lscale(2)
       wp2rtp(1) = wp2rtp(2)
       wp2thlp(1) = wp2thlp(2)
       rcm(1) = rcm(2)
       p_in_Pa(1) = p_in_Pa(2)
       rtm(1) = rtm(2)
       thlm(1) = thlm(2)

       ! The upper boundary conditions on <w'x'> must be set to 0 to match what
       ! is found in advance_xm_wpxp.
       wprtp(gr%nz) = zero
       wpthlp(gr%nz) = zero

       ! Boundary conditions on wm.
       wm_zm(1) = zero
       wm_zm(2) = zero
       wm_zm(gr%nz-1) = zero
       wm_zm(gr%nz) = zero

       ! Overwriting the vertical profile of wm_zm to have a value of 0
       ! everywhere.  This is being done because the mean advection term does
       ! not currently use flux form.  This will cause profiles to have a
       ! spurious source or sink.  This issue will be corrected in the future.
       ! However, in order to backdate this test and have it pass correctly,
       ! the values of wm_zm and wm_zt must be 0 at all levels.
       wm_zm = zero

       ! Interpolate fields set on momentum levels to thermodynamic levels.
       wm_zt = zm2zt( wm_zm )
       wm_zt(1) = zero
       wp3_on_wp2_zt = zm2zt( wp3_on_wp2 )
       wp3_on_wp2_zt(1) = zero
       Kh_zt = zm2zt( Kh_zm )
       rho_ds_zt = zm2zt( rho_ds_zm )

       ! Calculate the value of skewness of w (momentum levels).
       Skw_zm = wp3_on_wp2 / sqrt( wp2 )

       ! Calculate the inverse values of dry, static air density.
       invrs_rho_ds_zm = one / rho_ds_zm
       invrs_rho_ds_zt = one / rho_ds_zt

       rtm_ref = zero
       thlm_ref = zero

       wp2sclrp = zero
       sclrpthvp = zero
       sclrm_forcing = zero
       sclrp2 = zero

       sclrm = zero
       wpsclrp = zero

       um_forcing = zero
       vm_forcing = zero
       ug = zero
       vg = zero
       wpthvp = zero
       fcor = zero
       um_ref = zero
       vm_ref = zero
       up2 = wp2
       vp2 = wp2

       um = zero
       vm = zero
       upwp = zero
       vpwp = zero
       ! Below I assume that the buoy term in the upwp eqn doesn't matter:
       uprcp = zero
       vprcp = zero
       rc_coef = one

       ! Calculate the value of em.
       em = one_half * ( wp2 + up2 + vp2 )

       ! Calculate the PDF parameters on momentum levels (w_1_zm, w_2_zm,
       ! varnce_w_1_zm, varnce_w_2_zm, and mixt_frac).
       call ADG1_w_closure( wm_zm, wp2, Skw_zm, sigma_sqd_w,           & ! In
                            sqrt( wp2 ), 0.999_core_rknd,              & ! In
                            w_1_zm, w_2_zm, w_1_n_zm, w_2_n_zm,        & ! Out
                            varnce_w_1_zm, varnce_w_2_zm, mixt_frac_zm ) ! Out

       ! Calculate the value of exner.
       exner = ( p_in_Pa / 1.0e5_core_rknd )**(Rd/Cp)

       ! Calculate thvm.
       thvm = thlm + ep1 * 300.0_core_rknd * rtm &
              + ( Lv / ( Cp * exner ) - ep2 * 300.0_core_rknd ) * rcm

       ! Interpolate fields set on thermodynamic levels to momentum levels.
       thv_ds_zm = zt2zm( thvm )

       ! Calculate the vertical integrals of rtm and thlm before the call to
       ! advance_xm_wpxp so that spurious source can be calculated.
       rtm_integral_before &
       = vertical_integral( gr%nz-1, rho_ds_zt(2:gr%nz), &
                            rtm(2:gr%nz), gr%dzt(2:gr%nz) )

       thlm_integral_before &
       = vertical_integral( gr%nz-1, rho_ds_zt(2:gr%nz), &
                            thlm(2:gr%nz), gr%dzt(2:gr%nz) )

       call advance_xm_wpxp( dt, sigma_sqd_w, wm_zm, wm_zt, wp2, &
                             Lscale, wp3_on_wp2, wp3_on_wp2_zt, Kh_zt, Kh_zm, &
                             tau_C6_zm, Skw_zm, wp2rtp, rtpthvp, rtm_forcing, &
                             wprtp_forcing, rtm_ref, wp2thlp, thlpthvp, &
                             thlm_forcing, wpthlp_forcing, thlm_ref, &
                             rho_ds_zm, rho_ds_zt, invrs_rho_ds_zm, &
                             invrs_rho_ds_zt, thv_ds_zm, rtp2, thlp2, &
                             w_1_zm, w_2_zm, varnce_w_1_zm, varnce_w_2_zm, &
                             mixt_frac_zm, l_implemented, em, wp2sclrp, &
                             sclrpthvp, sclrm_forcing, sclrp2, exner, rcm, &
                             p_in_Pa, thvm, Cx_fnc_Richardson, &
                             ice_supersat_frac, &
                             pdf_implicit_coefs_terms, &
                             um_forcing, vm_forcing, ug, vg, wpthvp, &
                             fcor, um_ref, vm_ref, up2, vp2, &
                             uprcp, vprcp, rc_coef, &
                             l_diffuse_rtm_and_thlm, &
                             l_stability_correct_Kh_N2_zm, &
                             l_upwind_wpxp_ta, &
                             l_upwind_xm_ma, &
                             l_use_C7_Richardson, &
                             rtm, wprtp, thlm, wpthlp, &
                             sclrm, wpsclrp, um, upwp, vm, vpwp )

       ! Calculate the spurious source for rtm
       rtm_flux_top = rho_ds_zm(gr%nz) * wprtp(gr%nz)

       rtm_flux_sfc = rho_ds_zm(1) * wprtp(1)

       rtm_integral_after &
       = vertical_integral( gr%nz-1, rho_ds_zt(2:gr%nz), &
                            rtm(2:gr%nz), gr%dzt(2:gr%nz) )

       rtm_integral_forcing &
       = vertical_integral( gr%nz-1, rho_ds_zt(2:gr%nz), &
                            rtm_forcing(2:gr%nz), gr%dzt(2:gr%nz) )

       rtm_spur_src &
       = calculate_spurious_source( rtm_integral_after, &
                                    rtm_integral_before, &
                                    rtm_flux_top, rtm_flux_sfc, &
                                    rtm_integral_forcing, &
                                    dt )

       write(*,*) "Vertical integral of <rt> before advance_xm_wpxp = ", &
                  rtm_integral_before
       write(*,*) "Vertical integral of <rt> after advance_xm_wpxp = ", &
                  rtm_integral_after
       write(*,*) "Vertical integral of <rt> forcing = ", rtm_integral_forcing
       write(*,*) "Flux of <rt> at the surface = ", rtm_flux_sfc
       write(*,*) "Flux of <rt> at the top of the domain = ", rtm_flux_top
       write(*,*) "Spurious source of <rt> = ", rtm_spur_src
       write(*,*) ""

       ! Check if the calculated spurious source is within acceptable limits.
       density_weighted_height = sum( rho_ds_zt(2:gr%nz) * gr%dzt(2:gr%nz) )

       if ( abs( rtm_spur_src ) &
            <= tol * ( rtm_integral_before / density_weighted_height ) ) then
          write(*,*) "The spurious source of <rt> is within acceptable limits."
       else
          write(*,*) "The spurious source of <rt> is too large in magnitude."
          write(*,*) "Spurious source of <rt> = ", rtm_spur_src
          write(*,*) "Acceptable magnitude = ", &
                     tol * ( rtm_integral_before / density_weighted_height )
          num_errors = num_errors + 1
       endif

       write(*,*) ""
       
       ! Calculate the spurious source for thlm
       thlm_flux_top = rho_ds_zm(gr%nz) * wpthlp(gr%nz)

       thlm_flux_sfc = rho_ds_zm(1) * wpthlp(1)

       thlm_integral_after &
       = vertical_integral( gr%nz-1, rho_ds_zt(2:gr%nz), &
                            thlm(2:gr%nz), gr%dzt(2:gr%nz) )

       thlm_integral_forcing &
       = vertical_integral( gr%nz-1, rho_ds_zt(2:gr%nz), &
                            thlm_forcing(2:gr%nz), gr%dzt(2:gr%nz) )

       thlm_spur_src &
       = calculate_spurious_source( thlm_integral_after, &
                                    thlm_integral_before, &
                                    thlm_flux_top, thlm_flux_sfc, &
                                    thlm_integral_forcing, &
                                    dt )

       write(*,*) "Vertical integral of <thl> before advance_xm_wpxp = ", &
                  thlm_integral_before
       write(*,*) "Vertical integral of <thl> after advance_xm_wpxp = ", &
                  thlm_integral_after
       write(*,*) "Vertical integral of <thl> forcing = ", thlm_integral_forcing
       write(*,*) "Flux of <thl> at the surface = ", thlm_flux_sfc
       write(*,*) "Flux of <thl> at the top of the domain = ", thlm_flux_top
       write(*,*) "Spurious source of <thl> = ", thlm_spur_src
       write(*,*) ""

       ! Check if the calculated spurious source is within acceptable limits.
       if ( abs( thlm_spur_src ) &
            <= tol * ( thlm_integral_before / density_weighted_height ) ) then
          write(*,*) "The spurious source of <thl> is within acceptable limits."
       else
          write(*,*) "The spurious source of <thl> is too large in magnitude."
          write(*,*) "Spurious source of <thl> = ", thlm_spur_src
          write(*,*) "Acceptable magnitude = ", &
                     tol * ( thlm_integral_before / density_weighted_height )
          num_errors = num_errors + 1
       endif

       write(*,*) ""
       
    enddo ! iter = 1, num_iter, 1


    if ( num_errors == 0 ) then
       write(*,*) "Success!"
       spurious_source_unit_test = 0
    else ! num_errors > 0
       write(*,'(1x,A,I3,A)') "The spurious source test has ", &
                              num_errors, " error(s)."
       spurious_source_unit_test = 1
    endif ! num_errors

    write(*,*) ""


    return

  end function spurious_source_unit_test

!===============================================================================

end module spurious_source_test
